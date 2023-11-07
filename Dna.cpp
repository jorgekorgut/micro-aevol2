//
// Created by arrouan on 01/10/18.
//

#include "Dna.h"

#include <bit>
#include <bitset>
#include <cstdint>
#include <iostream>

#include <cassert>

Dna::Dna(int length, Threefry::Gen &&rng)
{
    // Generate a random genome
    char sequenceChar[length+1];
    for (int32_t i = 0; i < length; i++)
    {
        sequenceChar[i] = ('0' + rng.random(NB_BASE));
    }
    sequenceChar[length] = '\0';

    std::string randomSequence(sequenceChar);
    seq_ = Bitset(randomSequence, length);
}

Dna::Dna(const Dna &clone) : seq_(clone.seq_)
{

}

int Dna::length() const
{
    return seq_.bitsetSize();
}

void Dna::save(gzFile backup_file)
{
    size_t dna_length = length();
    gzwrite(backup_file, &dna_length, sizeof(dna_length));
    char tmp_seq[dna_length];
    for (size_t i = 0; i < dna_length; ++i)
    {
        tmp_seq[i] = seq_[i] + '0';
    }

    gzwrite(backup_file, tmp_seq, dna_length * sizeof(tmp_seq[0]));
}

void Dna::load(gzFile backup_file)
{
    int dna_length;
    gzread(backup_file, &dna_length, sizeof(dna_length));

    char tmp_seq[dna_length];
    gzread(backup_file, tmp_seq, dna_length * sizeof(tmp_seq[0]));

    std::string sequence = "";

    for (size_t i = 0; i < dna_length; ++i)
    {
        sequence += tmp_seq;
    }

    seq_ = Bitset(sequence, dna_length);
}

void Dna::set(int pos, char c)
{
    std::cout << "Dna::set Not implemented!" << std::endl;
    seq_.set(pos, (c == '1') ? 1 : 0);
}

/**
 * Remove the DNA inbetween pos_1 and pos_2
 *
 * @param pos_1
 * @param pos_2
 */
void Dna::remove(int pos_1, int pos_2)
{
    std::cout << "Dna::remove Not implemented!" << std::endl;
    /*
    assert(pos_1 >= 0 && pos_2 >= pos_1 && pos_2 <= seq_.size());

    auto seq_copy = boost::dynamic_bitset<>(seq_);

    if (pos_2 < seq_.size())
    {
        seq_.reset(pos_1, seq_.size() - pos_1);
        seq_copy.reset(0, pos_2);
        seq_copy >> pos_2 - pos_1;

        seq_ |= seq_copy;
    }

    seq_.resize(seq_.size() - (pos_2 - pos_1));
    */
}

/**
 * Insert a sequence of a given length at a given position into the DNA of the Organism
 *
 * @param pos : where to insert the sequence
 * @param seq : the sequence itself
 * @param seq_length : the size of the sequence
 */
void Dna::insert(int pos, Bitset &seq)
{
    std::cout << "Dna::insert Not implemented!" << std::endl;
    /*
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
    */
}

/**
 * Insert a sequence of a given length at a given position into the DNA of the Organism
 *
 * @param pos : where to insert the sequence
 * @param seq : the sequence itself
 * @param seq_length : the size of the sequence
 */
void Dna::insert(int pos, Dna *seq)
{
    // Insert sequence 'seq' at position 'pos'
    assert(pos >= 0 && pos < seq_.bitsetSize());

    insert(pos, seq->seq_);
}

void Dna::do_switch(int pos)
{
    seq_.flip(pos);
}

void Dna::do_duplication(int pos_1, int pos_2, int pos_3)
{
    // Duplicate segment [pos_1; pos_2[ and insert the duplicate before pos_3
    std::cout << "Dna::do_duplication Not implemented!" << std::endl;
    /*
    char *duplicate_segment = NULL;

    int32_t seg_length;
    if (pos_1 < pos_2)
    {
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
    }
    else
    { // if (pos_1 >= pos_2)
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
    */
}

int Dna::promoter_at(int pos)
{
    u_int64_t mask = seq_.getMask(pos, PROM_SIZE);
    u_int64_t comparation = mask ^ prom_seq;
    int dist_lead = std::popcount(comparation);

    return dist_lead;
}

// Given a, b, c, d boolean variable and X random boolean variable,
// a terminator look like : a b c d X X !d !c !b !a
int Dna::terminator_at(int pos)
{
    u_int64_t leftMask = seq_.getMask(pos, TERM_STEM_SIZE);
    u_int64_t rightMask = seq_.getMask(pos + TERM_SIZE - TERM_STEM_SIZE, TERM_STEM_SIZE);

    int dist_term_lead = 0;
    dist_term_lead += (leftMask & 0b0001) != ((rightMask & 0b1000) >> 3);
    dist_term_lead += (leftMask & 0b0010) != ((rightMask & 0b0100) >> 1);
    dist_term_lead += (leftMask & 0b0100) != ((rightMask & 0b0010) << 1);
    dist_term_lead += (leftMask & 0b1000) != ((rightMask & 0b0001) << 3);

    return dist_term_lead;
}

bool Dna::shine_dal_start(int pos)
{
    u_int64_t mask_start = seq_.getMask(pos, SHINE_DAL_SIZE);
    if (shine_dal_seq_start_ != mask_start)
    {
        return false;
    }

    u_int64_t mask_end = seq_.getMask(pos + SHINE_DAL_SIZE + SD_START_SPACER, CODON_SIZE);
    if (shine_dal_seq_end_ != mask_end)
    {
        return false;
    }

    return true;
}

bool Dna::protein_stop(int pos)
{
    u_int64_t mask = seq_.getMask(pos, CODON_SIZE);

    bool is_protein = (protein_end == mask);

    return is_protein;
}

int Dna::codon_at(int pos)
{
    int value = seq_.getMask(pos, CODON_SIZE);

    return value;
}