//
// Created by elturpin on 26/11/2020.
//

#include "cuProtein.cuh"
#include "aevol_constants.h"

#include <cassert>

// CODON_START = 0b000; 0
// CODON_STOP  = 0b001; 1
// CODON_W0    = 0b010; 2
// CODON_W1    = 0b011; 3
// CODON_M0    = 0b100; 4
// CODON_M1    = 0b101; 5
// CODON_H0    = 0b110; 6
// CODON_H1    = 0b111; 7

__device__ void cuProtein::add_codon(uint8_t codon) {
    assert(codon != 1);
    auto get_codon_type = [](uint8_t codon) -> uint8_t {
        codon = codon >> 1;
        if (codon)
            return codon - 1;
        return 2;
    };

    uint8_t codon_type = get_codon_type(codon);
    auto &parameter = wmh[codon_type];
    wmh_nb[codon_type]++;

    // Convert Gray code to Binary code 1 bit at a time
    // cf : https://www.elprocus.com/code-converter-binary-to-gray-code-and-gray-code-to-binary-conversion/
    uint8_t gray_bit_value = codon & 1;
    uint8_t previous_value = parameter & 1;
    parameter <<= 1;
    parameter += previous_value ^ gray_bit_value;
}

__device__ void cuProtein::normalize() {
    // normalize to x_double = x_int / (2^n - 1)
    if (wmh_nb[0])
        width = (double) wmh[0] / ((1u << wmh_nb[0]) - 1);
    else
        width = 0.0;
    if (wmh_nb[1])
        mean = (double) wmh[1] / ((1u << wmh_nb[1]) - 1);
    else
        mean = 0.5;
    if (wmh_nb[2])
        height = (double) wmh[2] / ((1u << wmh_nb[2]) - 1);
    else
        height = 0.5;

    width = (W_MAX - W_MIN) * width + W_MIN;
    mean = (X_MAX - X_MIN) * mean + X_MIN;
    height = (H_MAX - H_MIN) * height + H_MIN;

    if (wmh_nb[0] == 0 || wmh_nb[1] == 0 || wmh_nb[2] == 0 || width == 0 || height == 0)
        concentration = 0.0f; // Protein is not functional
}