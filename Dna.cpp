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
		tmp_seq[i] = seq_[i] + '0';
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
    auto val = seq_.get_subset(pos, PROM_SIZE);
    auto comparation = val ^ prom_seq;
    int dist_lead = std::popcount(comparation);

    return dist_lead;
}

// Given a, b, c, d boolean variable and X random boolean variable,
// a terminator look like : a b c d X X !d !c !b !a
int Dna::terminator_at(int pos) {
	auto right = seq_.get_subset(pos, TERM_STEM_SIZE);
	auto left = seq_.get_subset(pos + TERM_SIZE - TERM_STEM_SIZE, TERM_STEM_SIZE);

#if 0
	left = (left & 0xC) >> 2 | (left & 0x3) << 2;
	left = (left & 0xA) >> 1 | (left & 0x5) << 1;
	return std::popcount(left ^ right);

#elif 0
	static unsigned char lookup[16] = {
		0x0, 0x8, 0x4, 0xc, 0x2, 0xa, 0x6, 0xe,
		0x1, 0x9, 0x5, 0xd, 0x3, 0xb, 0x7, 0xf,
	};

	left = lookup[left];
	return std::popcount(left ^ right);
#else

	int dist_term_lead = 0;
	dist_term_lead += (left & 0b0001) != ((right & 0b1000) >> 3);
	dist_term_lead += (left & 0b0010) != ((right & 0b0100) >> 1);
	dist_term_lead += (left & 0b0100) != ((right & 0b0010) << 1);
	dist_term_lead += (left & 0b1000) != ((right & 0b0001) << 3);

    return dist_term_lead;
#endif
}

bool Dna::shine_dal_start(int pos) {
#if 1
	auto start = seq_.get_subset(pos, SHINE_DAL_SIZE);
	if (start != shine_dal_begin)
		return false;
	auto end = seq_.get_subset(pos + SHINE_DAL_SIZE + SD_START_SPACER, CODON_SIZE);
	if (end != shine_dal_end)
		return false;

	return true;
#else
	auto val = seq_.get_subset(pos, SD_TO_START);
	val &= ~seq_.bit_mask(SHINE_DAL_SIZE, SHINE_DAL_SIZE + SD_START_SPACER - 1);
	return val == shine_dal_seq;
#endif
}

bool Dna::protein_stop(int pos) {
	auto val = seq_.get_subset(pos, CODON_SIZE);
	return val == protein_end;
}

int Dna::codon_at(int pos) {
    return seq_.get_subset(pos, CODON_SIZE);
}