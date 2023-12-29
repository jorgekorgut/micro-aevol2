//
// Created by elturpin on 08/12/2020.
//

#pragma once

#include "Random123/threefry.h"

using RNG = r123::Threefry2x64;
using ctr_type = RNG::ctr_type;
using key_type = RNG::key_type;
using ctr_value_type = ctr_type::value_type;
using key_value_type = key_type::value_type;

enum Phase {
    SELECTION = 0, MUTATION = 1, NPHASES = 2
};

inline R123_CUDA_DEVICE double int64_to_double(int64_t number) {
    return (number&((1llu<<48)-1))/double(1llu<<48);
}

struct RandService {
    RNG generator{};

    uint phase_size{};

    key_type seed{};
    ctr_value_type* rng_counters{};

    inline R123_CUDA_DEVICE ctr_value_type gen_number(uint idx, Phase phase) {
        auto counter_idx = idx + phase_size * phase;
        auto counter = ++rng_counters[counter_idx];
        ctr_type ctr = {{counter_idx, counter}};
        return generator(ctr, seed)[0];
    }

    inline R123_CUDA_DEVICE double gen_double(uint idx, Phase phase) {
        return int64_to_double(gen_number(idx, phase));
    }

    inline R123_CUDA_DEVICE unsigned int gen_number_max(unsigned int max, uint idx, Phase phase) {
        return gen_double(idx, phase) * max;
    }

    R123_CUDA_DEVICE int32_t binomial_random(int32_t nb_drawings, double prob, uint idx, Phase phase);
    R123_CUDA_DEVICE int random_roulette(double* probability_array, int nb_elements, uint idx, Phase phase);
};
