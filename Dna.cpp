//
// Created by arrouan on 01/10/18.
//

#include "Dna.h"

#include <cassert>

Dna::Dna(int length, Threefry::Gen &&rng) : seq_(length) {
    // Generate a random genome
    //std::cout << "Dna : " << length << std::endl;
    for (int32_t i = 0; i < length; i++) {
        seq_[i] = rng.random(NB_BASE);
        //std::cout << "seq_["<<i<<"]" << seq_[i] <<  " ";
    }
    //std::cout << std::endl;
}

int Dna::length() const {
    return seq_.size();
}

void Dna::save(gzFile backup_file) {
    size_t dna_length = length();
    gzwrite(backup_file, &dna_length, sizeof(dna_length));
    char tmp_seq[dna_length];
    for (size_t i = 0; i < dna_length; ++i) {
		tmp_seq[i] = seq_[i];
    }

    gzwrite(backup_file, tmp_seq, dna_length * sizeof(tmp_seq[0]));
}

void Dna::load(gzFile backup_file) {
    int dna_length;
    gzread(backup_file, &dna_length, sizeof(dna_length));

    char tmp_seq[dna_length];
    gzread(backup_file, tmp_seq, dna_length * sizeof(tmp_seq[0]));

    seq_ = boost::dynamic_bitset<>(dna_length);
    for (size_t i = 0; i < dna_length; ++i) {
		seq_[i] = tmp_seq;
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
std::cout << "Insert" << std::endl;
// Insert sequence 'seq' at position 'pos'
    assert(pos >= 0 && pos < seq_.size());

    insert(pos, seq->seq_);
}

void Dna::do_switch(int pos) {
	seq_.flip(pos);
}

void Dna::do_duplication(int pos_1, int pos_2, int pos_3) {
    // Duplicate segment [pos_1; pos_2[ and insert the duplicate before pos_3
    std::cout << "Do_duplication" << std::endl;
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
    //std::cout << "promoter_at:" << pos << std::endl;
    int prom_dist[PROM_SIZE];

	int dist_lead = 0;
    for (int motif_id = 0; motif_id < PROM_SIZE; motif_id++) {
        int search_pos = pos + motif_id;
        // Circular search in a array
        if (search_pos >= seq_.size())
        {
            search_pos -= seq_.size();
        }
            
        // Searching for the promoter
        //std::cout << "bitset access: " << search_pos << std::endl;
        dist_lead += (PROM_SEQ[motif_id]) != seq_[search_pos];
    }

    //std::cout << "distance:" << dist_lead << std::endl;
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
    bool start = true;

	if (pos + SHINE_DAL_SIZE + CODON_SIZE + SD_START_SPACER > seq_.size()) {
		int t_pos, k_t;
		t_pos = pos;
		for (int k = 0; k < SHINE_DAL_SIZE; k++, t_pos++) {
			if (t_pos >= seq_.size())
				t_pos -= seq_.size();

			if (seq_[t_pos] != SHINE_DAL_SEQ[k]) {
				start = false;
				break;
			}
		}

		if (!start)
			return start;
		t_pos += SD_START_SPACER;

		for (int k = 0; k < CODON_SIZE; ++k, ++t_pos) {
			if (t_pos >= seq_.size())
				t_pos -= seq_.size();

			if (seq_[t_pos] != SHINE_DAL_SEQ[k + SHINE_DAL_SIZE + SD_START_SPACER]) {
				start = false;
				break;
			}
		}
	} else {
		auto shine_bitmask = boost::dynamic_bitset<>(SHINE_DAL_SEQ);
		shine_bitmask.resize(seq_.size());
		shine_bitmask >>= pos;

		shine_bitmask &= seq_;
		start = shine_bitmask.any();
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

        if (seq_[t_k] == PROTEIN_END[k]) {
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