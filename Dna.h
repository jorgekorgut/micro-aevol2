//
// Created by arrouan on 01/10/18.
//

#pragma once

#include <vector>
#include <zlib.h>

#include "Threefry.h"

#include <cstdint>
constexpr int8_t NB_BASE = 2;
constexpr int8_t CODON_SIZE = 3;
// promoter
constexpr int8_t MAX_PROM_ERROR = 4;
constexpr int8_t PROM_SIZE = 22;
constexpr const char *PROM_SEQ = "0101011001110010010110";
// terminator
constexpr int8_t TERM_STEM_SIZE = 4;
constexpr int8_t TERM_LOOP_SIZE = 2;
constexpr int8_t TERMINATOR_SIZE = TERM_STEM_SIZE + TERM_LOOP_SIZE + TERM_STEM_SIZE;
// shine dalgardo
constexpr int8_t SHINE_DAL_SIZE = 6;
constexpr int8_t SD_START_SPACER = 4;
constexpr int8_t SD_TO_START = SHINE_DAL_SIZE + SD_START_SPACER + CODON_SIZE;
constexpr const char *SHINE_DAL_SEQ = "011011****000";
// stop
constexpr const char *PROTEIN_END = "001"; // CODON_STOP

// Codon
constexpr int8_t CODON_START = 0b000;
constexpr int8_t CODON_M0    = 0b100;
constexpr int8_t CODON_M1    = 0b101;
constexpr int8_t CODON_W0    = 0b010;
constexpr int8_t CODON_W1    = 0b011;
constexpr int8_t CODON_H0    = 0b110;
constexpr int8_t CODON_H1    = 0b111;

class Dna {

public:
    Dna() = default;

    Dna(const Dna &clone) = default;

    Dna(int length, Threefry::Gen &&rng);

    ~Dna() = default;

    int length() const;

    void save(gzFile backup_file);

    void load(gzFile backup_file);

    void set(int pos, char c);

    /// Remove the DNA inbetween pos_1 and pos_2
    void remove(int pos_1, int pos_2);

    /// Insert a sequence of a given length at a given position into the DNA of the Organism
    void insert(int pos, std::vector<char> seq);

    /// Insert a sequence of a given length at a given position into the DNA of the Organism
    void insert(int pos, Dna *seq);

    void do_switch(int pos);

    void do_duplication(int pos_1, int pos_2, int pos_3);

    int promoter_at(int pos);

    int terminator_at(int pos);

    bool shine_dal_start(int pos);

    bool protein_stop(int pos);

    int codon_at(int pos);

    std::vector<char> seq_;
};
