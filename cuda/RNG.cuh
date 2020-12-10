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
    SELECTION = 0, MUTATION = 1, NPHASES
};

inline R123_CUDA_DEVICE double int64_to_double(int64_t number) {
    return (number&((1llu<<48)-1))/double(1llu<<48);
}

inline R123_CUDA_DEVICE int random_roulette(double rand_number ,double* probability_array, int nb_elements) {
    int selection = -1;
    do
    {
        assert(selection < nb_elements - 1);
        rand_number -= probability_array[++selection];
    } while (rand_number > 0);

    return selection;
}

struct RandService {
    RNG generator;

    uint nb_phase;
    uint phase_size;

    key_type seed;
    ctr_value_type* rng_counters;

    inline R123_CUDA_DEVICE ctr_value_type gen_number(uint idx, int phase) {
        auto counter_idx = idx + phase * nb_phase;
        auto counter = ++rng_counters[counter_idx];
        ctr_type ctr = {{counter_idx, counter}};
        return generator(ctr, seed)[0];
    }

    inline R123_CUDA_DEVICE double gen_double(uint idx, int phase) {
        return int64_to_double(gen_number(idx, phase));
    }
};
