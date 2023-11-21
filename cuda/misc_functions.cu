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

__device__
uint
count_bitset(block* set, uint size)
{
  uint num = 0;
  uint i;
  uint mod = size & 63;

  // TODO: do not hardcode 64
  for (i = 0; i < size / 64; ++i)
    num += __popc(set[i]);
  if (mod)
    num += __popc(set[i] & ((1 << mod) - 1));

  return num;
}

__device__
uint
sparse_bitset(block* set, uint size, uint* idcs)
{
  uint idx = 0;

  for (uint i = 0; i < size; ++i) {
    auto val = set[i];

    block bitmask = 1;
    for (uint j = 0; j < 64; ++j, bitmask <<= 1) {
      if (val & bitmask)
        idcs[idx++] = i * 64 + j;
    }
  }

  return idx;
}

