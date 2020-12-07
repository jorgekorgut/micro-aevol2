//
// Created by elturpin on 03/12/2020.
//

#pragma once

#include "Abstract_ExpManager.h"
#include "ExpManager.h"

class cuIndividual;

class cuExpManager: public Abstract_ExpManager {
public:
    explicit cuExpManager(const ExpManager* cpu_exp);

    void save(int t) override;

    void run_evolution(int nb_gen) override;

    void load(int t) override;

private:
    // Host Data
    int nb_indivs_;

    int dna_length_;
    char **host_organisms_;

    unsigned int seed_;
    size_t nb_counter_;
    unsigned long long *counters_;

    double *target_;

    int grid_height_;
    int grid_width_;

    double mutation_rate_;

    int backup_step_;

    // Interface Host - Device
    void transfer_to_device();

    // Device data
    cuIndividual *device_organisms_;

};


