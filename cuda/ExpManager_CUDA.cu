//
// Created by elturpin on 03/12/2020.
//

#include "ExpManager_CUDA.h"

#include <cstdio>
#include <chrono>
#include "nvToolsExt.h"
#include <cuda_profiler_api.h>
#include <AeTime.h>

#include "ExpManager.h"

using namespace std::chrono;
using namespace std;

ExpManager_CUDA::ExpManager_CUDA(const ExpManager* cpu_exp) {
    grid_height_ = cpu_exp->grid_height_;
    grid_width_ = cpu_exp->grid_height_;

    mutation_rate_ = cpu_exp->mutation_rate_;

    backup_step_ = cpu_exp->grid_height_;

    nb_indivs_ = grid_height_ * grid_width_;

    dna_length_ = cpu_exp->internal_organisms_[0]->length();
    host_organisms_ = new char *[nb_indivs_];
    for (int i = 0; i < nb_indivs_; ++i) {
        const auto& org = cpu_exp->internal_organisms_[i];
        host_organisms_[i] = org->dna_->seq_.data();
    }

    target_ = cpu_exp->target;

    seed_ = cpu_exp->seed_;
    nb_counter_ = cpu_exp->rng_->counters().size();
    counters_ = new unsigned long long[nb_counter_];
    for (int i = 0; i < nb_counter_; ++i) {
        counters_[i] = cpu_exp->rng_->counters()[i];
    }
}

void ExpManager_CUDA::run_evolution(int nb_gen) {
    cudaProfilerStart();
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    cout << "Transfer" << endl;
    transfer_to_device();
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration_transfer_in = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    cout << "Transfer done in " << duration_transfer_in << endl;

    printf("Running evolution GPU from %d to %d\n",AeTime::time(),AeTime::time()+nb_gen);
    bool firstGen = true;
    for (int gen = 0; gen < nb_gen+1; gen++) {
        if(gen == 91) nvtxRangePushA("generation 91 to 100");
        AeTime::plusplus();

        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        run_a_step_on_GPU(nb_indivs_, w_max_, selection_pressure_, grid_width_, grid_height_,mutation_rate_);

        t2 = high_resolution_clock::now();
        auto duration_transfer_in = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

        std::cout<<"LOG,"<<duration_transfer_in<<std::endl;

        firstGen = false;
        if(gen == 100) nvtxRangePop();
        printf("Generation %d : \n",AeTime::time());
    }
    cudaProfilerStop();
}

void ExpManager_CUDA::save(int t) {
    printf("Oups, not supported !\n");
}

void ExpManager_CUDA::load(int t) {
    printf("Oups, not supported !\n");
}

void ExpManager_CUDA::transfer_to_device() {
    // Transfer data from individual to device

    // Evaluate all the individuals

    // Transfer data from prng
}
