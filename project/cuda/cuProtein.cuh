//
// Created by elturpin on 26/11/2020.
//

#pragma once

#include <cstdint>

struct cuProtein {
  uint wmh[3]{};    // { Width, Mean, Height }
  uint wmh_nb[3]{};

  float concentration{};
  double width{};
  double mean{};
  double height{};

  __device__ void add_codon(uint8_t codon);
  __device__ void normalize();

  inline __device__ bool is_functional() const {
    return concentration != 0.0f;
  };
};

