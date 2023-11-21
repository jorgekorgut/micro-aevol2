//
// Created by elturpin on 16/11/2020.
//

#include "misc_functions.cuh"

#include "aevol_constants.h"

__device__ uint8_t is_promoter(const block sequence) {
  uint32_t comparation = sequence ^ prom_seq;
  int dist_lead = __popc(comparation);
  // int dist_lead = std::popcount(comparation);

  if (dist_lead > PROM_MAX_DIFF)
    return 0b1111;
  return dist_lead;
}

__device__ bool is_terminator(const block sequence) {
  int left, right;
  for (left = 0, right = TERM_SIZE - 1; left < TERM_STEM_SIZE; ++left, --right) {
    if (sequence & (1 << left) == sequence & (1 << right))
      return false;
  }
  return true;
}

__device__ bool is_prot_start(const block sequence) {
  if (sequence & 0b111111 != shine_dal_begin)
    return false;
  if ((sequence >> SHINE_DAL_SIZE + SD_START_SPACER) & 0b111 != shine_dal_end)
    return false;

  return true;
}

__device__ uint8_t translate_to_codon(const block* seq) {
  uint8_t codon = 0;

  for (uint8_t i = 0; i < CODON_SIZE; ++i) {
    // codon += seq[i] << (CODON_SIZE - 1 - i);
    codon += seq[i] << i;
  }

  return codon;
}
