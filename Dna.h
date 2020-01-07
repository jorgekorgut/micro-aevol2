//
// Created by arrouan on 01/10/18.
//

#pragma once

#include <cstdio>
#include <cstring>
#include <cassert>
#include <cstdint>
#include <vector>
#include <zlib.h>

#include "Threefry.h"

constexpr int8_t CODON_SIZE = 3;

constexpr const char* PROM_SEQ = "0101011001110010010110";
constexpr const char* SHINE_DAL_SEQ = "011011000";
constexpr const char* PROTEIN_END = "001"; // CODON_STOP

class ExpManager;

class Dna {

 public:
  Dna() = default;

  Dna(const Dna& clone);

  Dna(int length, Threefry::Gen& rng);

  Dna(char* genome, int length);

  Dna(int length);

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
  void insert(int pos, Dna* seq);

  void do_switch(int pos);

  void do_duplication(int pos_1, int pos_2, int pos_3);

  int promoter_at(int pos);

  int terminator_at(int pos);

  bool shine_dal_start(int pos);

  bool protein_stop(int pos);

  int codon_at(int pos);

  std::vector<char> seq_;
};
