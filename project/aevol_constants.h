//
// Created by elturpin on 04/12/2020.
//

#pragma once

#include <cstdint>

#include <boost/dynamic_bitset.hpp>
#include <boost/align/aligned_allocator.hpp>

typedef u_int64_t block;
typedef boost::dynamic_bitset<block, boost::alignment::aligned_allocator<block, 16>> bitset;
#define blockSizeBites (sizeof(block) * 8)

constexpr int8_t NB_BASE = 2;
constexpr int8_t CODON_SIZE = 3;
// promoter
constexpr int8_t PROM_MAX_DIFF = 4;
constexpr int8_t PROM_SIZE = 22;
const std::string PROM_SEQ("0101011001110010010110");
// terminator
constexpr int8_t TERM_STEM_SIZE = 4;
constexpr int8_t TERM_LOOP_SIZE = 3;
constexpr int8_t TERM_SIZE = TERM_STEM_SIZE + TERM_LOOP_SIZE + TERM_STEM_SIZE;
// shine dalgardo
constexpr int8_t SHINE_DAL_SIZE = 6;
constexpr int8_t SD_START_SPACER = 4;
constexpr int8_t SD_TO_START = SHINE_DAL_SIZE + SD_START_SPACER + CODON_SIZE;
const std::string SHINE_DAL_SEQ("0110110000000");

const std::string SHINE_DAL_SEQ_START("011011");
const std::string SHINE_DAL_SEQ_END("000");
// stop
const std::string PROTEIN_END("001");

constexpr int32_t DO_TRANSLATION_LOOP = SHINE_DAL_SIZE + SD_START_SPACER + 3 * CODON_SIZE;

// Codon
constexpr block CODON_START = 0b000;
constexpr block CODON_STOP = 0b100; //0b001
constexpr block CODON_M0 = 0b001; //0b100; 
constexpr block CODON_M1 = 0b101;
constexpr block CODON_W0 = 0b010;
constexpr block CODON_W1 = 0b110; //0b011;
constexpr block CODON_H0 = 0b011;//0b110;
constexpr block CODON_H1 = 0b111;

// Protein / Fuzzy space
constexpr double X_MIN = 0.0;
constexpr double X_MAX = 1.0;
constexpr double Y_MIN = 0.0;
constexpr double Y_MAX = 1.0;
constexpr double H_MIN = -1.0;
constexpr double H_MAX = 1.0;
constexpr double W_MIN = 0.0;
constexpr double W_MAX = 0.1;

constexpr int FUZZY_SAMPLING = 300;
constexpr int SELECTION_PRESSURE = 1000;

// Selection
constexpr int8_t NEIGHBORHOOD_WIDTH = 3;
constexpr int8_t NEIGHBORHOOD_HEIGHT = 3;
constexpr int8_t NEIGHBORHOOD_SIZE   = NEIGHBORHOOD_HEIGHT * NEIGHBORHOOD_WIDTH;

const block prom_seq = 0b0110100100111001101010lu;
const block shine_dal_seq = 0b0000000110110;
const block shine_dal_begin = 0b110110;
const block shine_dal_end = 0b000;
const block protein_end = 0b100;
