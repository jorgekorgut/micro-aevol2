//
// Created by elturpin on 16/11/2020.
//

#include "misc_functions.cuh"

#include "aevol_constants.h"

__device__ uint8_t is_promoter(const block sequence) {
  uint32_t comparation = prom_seq ^ (sequence & ((1ul << PROM_SIZE) - 1));

  // int dist_lead = __popc(comparation);
  // int dist_lead = std::popcount(comparation);

  uint dist_lead = 0;
  for (; comparation; comparation >>= 1) {
    dist_lead += comparation & 1;

    if (dist_lead > PROM_MAX_DIFF)
      return 0b1111;
  }

  return dist_lead;
}

__device__ bool is_terminator(const block sequence) {
  int left, right;
  for (left = 0, right = TERM_SIZE - 1; left < TERM_STEM_SIZE; ++left, --right) {
    if (((sequence >> left) & 1) == ((sequence >> right) & 1))
      return false;
  }
  return true;
}

__device__ bool is_prot_start(const block sequence) {
  if ((sequence & 0b111111) != shine_dal_begin)
    return false;
  if (((sequence >> SHINE_DAL_SIZE + SD_START_SPACER) & 0b111) != shine_dal_end)
    return false;

  return true;
}

// TODO
__device__ uint8_t translate_to_codon(const block* seq) {
  uint8_t codon = 0;

  for (uint8_t i = 0; i < CODON_SIZE; ++i) {
    // codon += seq[i] << (CODON_SIZE - 1 - i);
    codon += seq[i] << i;
  }

  return codon;
}

__device__
inline int
fast_mod(const int input, const int ceil)
{
	// TODO: Use (& (ceil - 1)) for powers of 2
    // apply the modulo operator only when needed
    // (i.e. when the input is greater than the ceiling)
    return input >= ceil ? input % ceil : input;
    // NB: the assumption here is that the numbers are positive
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
    num += __popc(set[i] & ((1llu << mod) - 1));

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

// Inspired by
// https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#atomic-functions
__device__
uint64_t
atomicOr(uint64_t* address, uint64_t val)
{
	unsigned long long int* address_as_ull = (unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;

	do {
		assumed = old;
		old = atomicCAS(address_as_ull, assumed, val | assumed);

	// Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
	} while (assumed != old);

	return old;
}

__device__
uint64_t
atomicAnd(uint64_t* address, uint64_t val)
{
	unsigned long long int* address_as_ull = (unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;

	do {
		assumed = old;
		old = atomicCAS(address_as_ull, assumed, val & assumed);

	// Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
	} while (assumed != old);

	return old;
}

__device__
const bool
set_bit_to(block* bitset, uint pos, bool value)
{
    // TODO: use a shift and & - 1
    uint bidx = pos / blockSizeBites;
    uint idx = pos % blockSizeBites;

    if (value)
        atomicOr(bitset + bidx, 1llu << idx);
    else
        atomicAnd(bitset + bidx, ~(1llu << idx));

    return value;
}

__device__
const bool
set_bit_to_unsafe(block* bitset, uint pos, bool value)
{
    // TODO: use a shift and & - 1
    uint bidx = pos / blockSizeBites;
    uint idx = pos % blockSizeBites;

    if (value)
        bitset[bidx] |= (1ull << idx);
    else
        bitset[bidx] &= ~(1ull << idx);

    return value;
}

__device__
inline void
set_bit(block* bitset, uint pos)
{
    // TODO: use a shift and & - 1
    uint bidx = pos / blockSizeBites;
    uint idx = pos % blockSizeBites;

	atomicOr(bitset + bidx, 1llu << idx);
}

__device__
inline void
set_bit_unsafe(block* bitset, uint pos)
{
    // TODO: use a shift and & - 1
    uint bidx = pos / blockSizeBites;
    uint idx = pos % blockSizeBites;

    bitset[bidx] |= (1llu << idx);
}

__device__
void
flip_bit(block* set, uint pos)
{
	// TODO: do not hardcode
	// TODO: use a shift and & - 1
	uint bid = pos / 64;
	uint idx = pos % 64;

	uint mask = 1llu << idx;

	if (set[bid] & mask)
		atomicAnd(set + bid, ~mask);
	else
		atomicOr(set + bid, mask);
}

__device__
void
flip_bit_unsafe(block* set, uint pos)
{
	// TODO: do not hardcode
	// TODO: use a shift and & - 1
	uint bid = pos / 64;
	uint idx = pos % 64;

	uint mask = 1llu << idx;

	if (set[bid] & mask)
		set[bid] &= ~mask;
	else
		set[bid] |= mask;
}

__device__
const
block get_block_circ(block* genome, uint size, uint bit_index, uint length)
{

  if(bit_index >= size){
      bit_index -= size;
  }

  block value = 0;
  int blockIndex = bit_index / blockSizeBites;
  int bitIndex = fast_mod(bit_index, blockSizeBites);

  // If the mask is truncated by blocks or end of bitset
  if (bitIndex + length >= blockSizeBites || bit_index + length >= size)
  {
      int nextBlockIndex = blockIndex + 1;
      // uint blockCount = std::ceil(size / (float)blockSizeBites);
      // uint blockCount = size / blockSizeBites + !!(fast_mod(size, blockSizeBites));
      uint blockCount = size / blockSizeBites + !!(size & (blockSizeBites - 1));

      if (nextBlockIndex > blockCount - 1)
      {
          nextBlockIndex = 0;
      }

      block leftMask = 1;
      leftMask <<= length;
      leftMask -= 1;

      value = genome[blockIndex] >> bitIndex | genome[nextBlockIndex] << (blockSizeBites - bitIndex);

      if (bit_index + length >= size)
      {
          int lastBlockRealSize = fast_mod(size, blockSizeBites);
          value |= genome[nextBlockIndex] << (blockSizeBites - bitIndex) + lastBlockRealSize;
      }

      value = value & leftMask;
  }
  else
  {
      block rightMask = 1;
      rightMask <<= bitIndex + length;
      rightMask -= 1;

      value = (genome[blockIndex] & rightMask) >> bitIndex;
  }

  return value;
}

// No safety checks.
__device__
const block
get_block(block* genome, uint pos, uint len)
{
	// TODO: use a shift and & - 1
	uint bidx = pos / blockSizeBites;
	uint pos_idx = pos % blockSizeBites;

	block value = (genome[bidx] >> pos_idx) & ((1llu << len) - 1);

	if (pos_idx + len > blockSizeBites) {
		uint start_nbits = blockSizeBites - pos_idx;
		value |= (genome[bidx + 1] << start_nbits) & ((1llu << len) - 1);
	}

	return value;
}

__device__
void
convert_char_to_bitset(const char* arr, uint size, block* set)
{
    // TODO: do not hardcode 64
    block new_block;
    uint i, bid;

    for (i = bid = 0; i < (size & ~63); i += 64, ++bid) {
        new_block = 0;
        for (uint j = 0; j < 64; ++j)
            new_block |= (arr[i + j] - '0') << j;

        set[bid] = new_block;
    }

    new_block = 0;
    for (uint j = 0; i < size; ++i, ++j)
        new_block |= (arr[i] - '0') << j;
    set[bid] = new_block;
}

__device__
void
print_bitset(block* set, uint bit_size)
{
    uint size = std::ceil(bit_size / 64.0);
    for (uint idx = (bit_size % 64); idx; --idx)
        printf("%u", 1 & (set[size - 1] >> (idx - 1)));
    --size;

    for (; size; --size) {
        for (uint idx = 64; idx; --idx) {
            printf("%u", 1 & (set[size - 1] >> (idx - 1)));
        }
    }
    printf("\n");
}

__global__
void
print_indivs(uint nb_indivs, cuIndividual* indivs)
{
    for (int indiv_idx = 0; indiv_idx < 1 /* nb_indivs */; ++indiv_idx) {
        const auto& indiv = indivs[indiv_idx];

        printf("size: %u\n", indiv.size);
        printf("block_size: %u\n", indiv.block_size);
        printf("genom: ");
        print_bitset(indiv.genome, indiv.size);

        printf("promoters: ");
        for (uint i = indiv.size - 1; i; --i) {
                printf("%u|", indiv.promoters[i]);
        }
        printf("%u\n", indiv.promoters[0]);

        printf("terminators: ");
        if (!indiv.terminator_idxs[0]) {
            print_bitset(indiv.terminators, indiv.size);
        } else {
            for (uint i = indiv.nb_terminator - 1; i; --i) {
                printf("%u, ", indiv.terminator_idxs[i]);
            }
            printf("%u\n", indiv.terminator_idxs[0]);
        }
        printf("prot_start: ");
        if (!indiv.prot_start_idxs[0]) {
            print_bitset(indiv.prot_start, indiv.size);
        } else {
            for (uint i = indiv.nb_prot_start - 1; i; --i) {
                printf("%u, ", indiv.prot_start_idxs[i]);
            }
            printf("%u\n", indiv.prot_start_idxs[0]);
        }

        printf("nb_terminator: %u\n", indiv.nb_terminator);
        printf("nb_prot_start: %u\n", indiv.nb_prot_start);
        // terminator_idxs
        // prot_start_idxs

        printf("nb_rnas: %u\n", indiv.nb_rnas);
        indiv.print_rnas();

        printf("nb_gene: %u\n", indiv.nb_gene);
        indiv.print_gathered_genes();
        indiv.print_proteins();

        indiv.print_phenotype();
        printf("fitness: %1.10e\n", indiv.fitness);
    }
    printf("\n");
}
