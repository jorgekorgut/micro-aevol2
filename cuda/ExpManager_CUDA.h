//
// Created by elturpin on 03/12/2020.
//

#pragma once

#include "Abstract_ExpManager.h"
#include "ExpManager.h"

class ExpManager_CUDA: public Abstract_ExpManager {
public:
    explicit ExpManager_CUDA(const ExpManager* cpu_exp);

    void save(int t) override;

    void run_evolution(int nb_gen) override;

    void load(int t) override;

private:

    void transfer_to_device();

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
};


