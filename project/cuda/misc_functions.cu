//
// Created by elturpin on 16/11/2020.
//

#include "misc_functions.cuh"

#include "aevol_constants.h"

__device__ uint8_t is_promoter(const char* sequence) {
  uint8_t distance = 0;
  for (int offset = 0; offset < PROM_SIZE; ++offset) {
    if (sequence[offset] != PROM_SEQ[offset]){
      distance++;
      if (distance > PROM_MAX_DIFF)
        return 0b1111;
    }
  }
  return distance;
}

__device__ bool is_terminator(const char* sequence) {
  int left, right;
  for (left = 0, right = TERM_SIZE - 1; left < TERM_STEM_SIZE; ++left, --right) {
    if (sequence[left] == sequence[right])
      return false;
  }
  return true;
}

__device__ bool is_prot_start(const char* sequence) {
  for (int offset = 0; offset < SHINE_DAL_SIZE; ++offset) {
    if (sequence[offset] != SHINE_DAL_SEQ[offset])
      return false;
  }
  for (int offset = 0; offset < CODON_SIZE; ++offset) {
    if (sequence[offset + SHINE_DAL_SIZE + SD_START_SPACER] != '0')
      return false;
  }
  return true;
}

__device__ uint8_t translate_to_codon(const char* seq) {
  uint8_t codon = 0;

  for (uint8_t i = 0; i < CODON_SIZE; ++i) {
    codon += (seq[i] - '0') << (CODON_SIZE - 1 - i);
  }

  return codon;
}
