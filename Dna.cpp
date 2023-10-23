//
// Created by arrouan on 01/10/18.
//

#include "Dna.h"

#include <cassert>

Dna::Dna(int length, Threefry::Gen &&rng) : seq_(length) {
    // Generate a random genome
    for (int32_t i = 0; i < length; i++) {
        seq_[i] = '0' + rng.random(NB_BASE);
    }
}

int Dna::length() const {
    return seq_.size();
}

void Dna::save(gzFile backup_file) {
    int dna_length = length();
    gzwrite(backup_file, &dna_length, sizeof(dna_length));

    char tmp_seq[dna_length];
    for (size_t i = 0; i < dna_length; ++i) {
		tmp_seq[i] = seq_[i] + '0';
    }

    gzwrite(backup_file, tmp_seq, dna_length * sizeof(seq_[0]));
}

void Dna::load(gzFile backup_file) {
    int dna_length;
    gzread(backup_file, &dna_length, sizeof(dna_length));

    char tmp_seq[dna_length];
    gzread(backup_file, tmp_seq, dna_length * sizeof(tmp_seq[0]));

    seq_ = boost::dynamic_bitset<>(dna_length);
    for (size_t i = 0; i < dna_length; ++i) {
		seq_[i] = tmp_seq - '0';
    }
}

void Dna::set(int pos, char c) {
    seq_[pos] = c;
}

/**
 * Remove the DNA inbetween pos_1 and pos_2
 *
 * @param pos_1
 * @param pos_2
 */
void Dna::remove(int pos_1, int pos_2) {
    assert(pos_1 >= 0 && pos_2 >= pos_1 && pos_2 <= seq_.size());

	auto seq_copy = boost::dynamic_bitset<>(seq_);

	if (pos_2 < seq_.size()) {
		seq_.reset(pos_1, seq_.size() - pos_1);
		seq_copy.reset(0, pos_2);
		seq_copy >> pos_2 - pos_1;

		seq_ |= seq_copy;
	}

    seq_.resize(seq_.size() - (pos_2 - pos_1));
}

/**
 * Insert a sequence of a given length at a given position into the DNA of the Organism
 *
 * @param pos : where to insert the sequence
 * @param seq : the sequence itself
 * @param seq_length : the size of the sequence
 */
void Dna::insert(int pos, boost::dynamic_bitset<> &seq) {
// Insert sequence 'seq' at position 'pos'
    assert(pos >= 0 && pos < seq_.size());

    seq_.resize(seq_.size() + seq.size(), false);

	auto post_bitmask = boost::dynamic_bitset<>(seq_);
	post_bitmask.reset(0, pos);
	post_bitmask <<= seq.size();

	auto seq_bitmask = boost::dynamic_bitset<>(seq);
	seq_bitmask.resize(seq_.size(), false);
	seq_bitmask <<= pos;

	seq_.reset(pos, seq_.size() - pos);

	seq_ |= seq_bitmask;
	seq_ |= post_bitmask;
}

/**
 * Insert a sequence of a given length at a given position into the DNA of the Organism
 *
 * @param pos : where to insert the sequence
 * @param seq : the sequence itself
 * @param seq_length : the size of the sequence
 */
void Dna::insert(int pos, Dna *seq) {
// Insert sequence 'seq' at position 'pos'
    assert(pos >= 0 && pos < seq_.size());

    insert(pos, seq->seq_);
}

void Dna::do_switch(int pos) {
	seq_.flip(pos);
}

void Dna::do_duplication(int pos_1, int pos_2, int pos_3) {
    // Duplicate segment [pos_1; pos_2[ and insert the duplicate before pos_3
    char *duplicate_segment = NULL;

    int32_t seg_length;

    if (pos_1 < pos_2) {
        //
        //       pos_1         pos_2                   -> 0-
        //         |             |                   -       -
        // 0--------------------------------->      -         -
        //         ===============                  -         - pos_1
        //           tmp (copy)                      -       -
        //                                             -----      |
        //                                             pos_2    <-'
        //
        auto seq_dupl = boost::dynamic_bitset<>(seq_);
        printf("%lu\n", seq_.size());
        seq_dupl >>= pos_1;
        seq_dupl.resize(pos_2 - pos_1);

        insert(pos_3, seq_dupl);
    } else { // if (pos_1 >= pos_2)
        // The segment to duplicate includes the origin of replication.
        // The copying process will be done in two steps.
        //
        //                                            ,->
        //    pos_2                 pos_1            |      -> 0-
        //      |                     |                   -       - pos_2
        // 0--------------------------------->     pos_1 -         -
        // ======                     =======            -         -
        //  tmp2                        tmp1              -       -
        //                                                  -----
        //
        //
        printf("%lu\n", seq_.size());
        auto seq_begin = boost::dynamic_bitset<>(seq_);
		seq_begin.resize(pos_2);

        auto seq_end = boost::dynamic_bitset<>(seq_);
        seq_end >>= pos_1;
		seq_begin.resize(seq_end.size() - pos_1);

        insert(pos_3, seq_end);
        insert(pos_3, seq_begin);
    }
}

int Dna::promoter_at(int pos) {
    int prom_dist[PROM_SIZE];

	int dist_lead = 0;
    for (int motif_id = 0; motif_id < PROM_SIZE; motif_id++) {
        int search_pos = pos + motif_id;
        if (search_pos >= seq_.size())
            search_pos -= seq_.size();
        // Searching for the promoter
        dist_lead += (PROM_SEQ[motif_id] - '0') != seq_[search_pos];

    }

    return dist_lead;
}

// Given a, b, c, d boolean variable and X random boolean variable,
// a terminator look like : a b c d X X !d !c !b !a
int Dna::terminator_at(int pos) {
    int term_dist[TERM_STEM_SIZE];
    for (int motif_id = 0; motif_id < TERM_STEM_SIZE; motif_id++) {
        int right = pos + motif_id;
        int left = pos + (TERM_SIZE - 1) - motif_id;

        // loop back the dna inf needed
        if (right >= length()) right -= length();
        if (left >= length()) left -= length();

        // Search for the terminators
        term_dist[motif_id] = seq_[right] != seq_[left] ? 1 : 0;
    }
    int dist_term_lead = term_dist[0] +
                         term_dist[1] +
                         term_dist[2] +
                         term_dist[3];

    return dist_term_lead;
}

bool Dna::shine_dal_start(int pos) {
    bool start = false;
    int t_pos, k_t;

    for (int k = 0; k < SHINE_DAL_SIZE + CODON_SIZE; k++) {
        k_t = k >= SHINE_DAL_SIZE ? k + SD_START_SPACER : k;
        t_pos = pos + k_t;
        if (t_pos >= seq_.size())
            t_pos -= seq_.size();

        if (seq_[t_pos] == (SHINE_DAL_SEQ[k_t] - '0')) {
            start = true;
        } else {
            start = false;
            break;
        }
    }

    return start;
}

bool Dna::protein_stop(int pos) {
    bool is_protein;
    int t_k;

    for (int k = 0; k < CODON_SIZE; k++) {
        t_k = pos + k;
        if (t_k >= seq_.size())
            t_k -= seq_.size();

        if (seq_[t_k] == PROTEIN_END[k] - '0') {
            is_protein = true;
        } else {
            is_protein = false;
            break;
        }
    }

    return is_protein;
}

int Dna::codon_at(int pos) {
    int value = 0;

    int t_pos;

    for (int i = 0; i < CODON_SIZE; i++) {
        t_pos = pos + i;
        if (t_pos >= seq_.size())
            t_pos -= seq_.size();
        if (seq_[t_pos])
            value += 1 << (CODON_SIZE - i - 1);
    }

    return value;
}