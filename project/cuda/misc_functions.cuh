//
// Created by elturpin on 16/11/2020.
//

#pragma once

#include <cuda.h>
#include <cuda_runtime_api.h>

/// General Purpose
template <typename T>
__device__ int sparse(int size, T* sparse_collection){
  int insert_position = 0;

  for (int read_position = 0; read_position < size; ++read_position) {
    auto read_value = sparse_collection[read_position];
    if (read_value) {
      sparse_collection[insert_position] = read_position;
      insert_position++;
    }
  }
  if (insert_position < size) {
    sparse_collection[insert_position] = 0;
  }
  return insert_position;
}

template <typename T>
__device__ uint find_smallest_greater(T value, const T* array, uint size){
  if (value > array[size-1])
    return 0;
  uint min = 0;
  uint max = size - 1;
  while (min < max) {
    uint mid = (max + min) / 2;
    if (value > array[mid]) {
      min = mid + 1;
    } else {
      max = mid;
    }
  }
  return min;
}

template<typename T>
__device__ T clamp(T x, T a, T b)
{
    return max(a, min(b, x));
}


/// Specific to Aevol model

__device__ uint8_t is_promoter(const char* sequence);

__device__ bool is_terminator(const char* sequence);

__device__ bool is_prot_start(const char* sequence);

__device__ uint8_t translate_to_codon(const char* seq);

inline __device__ uint get_distance(uint a, uint b, uint size){
  if (a > b)
    return (b + size) - a;
  return b - a;
}

inline __device__ uint get_distance_ori(uint a, uint b, uint size){
  if (a >= b)
    return (b + size) - a;
  return b - a;
}