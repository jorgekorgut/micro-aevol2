//
// Created by elturpin on 03/12/2020.
//

#pragma once

#include "Abstract_ExpManager.h"
#include "ExpManager.h"
#include "RNG.cuh"

class cuIndividual;

class cuExpManager: public Abstract_ExpManager {
public:
    explicit cuExpManager(const ExpManager* cpu_exp);

    ~cuExpManager() override;

    void run_evolution(int nb_gen) override;

    void save(int t) override;

    void load(int t) override;

private:
    void run_a_step();

    // Interface Host - Device
    void transfer_to_device();

    // Host Data
    int nb_indivs_;

    int dna_length_;
    char** host_organisms_;

    key_value_type seed_;
    size_t nb_counter_;
    ctr_value_type* counters_;

    double* target_;

    int grid_height_;
    int grid_width_;

    double mutation_rate_;

    int backup_step_;

    // Device data
    cuIndividual* device_organisms_;

    double* device_target_;

    key_type* device_seed_;
    ctr_value_type* device_rng_counters;
};


