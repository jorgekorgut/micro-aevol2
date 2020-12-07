//
// Created by elturpin on 03/12/2020.
//

#include "cuExpManager.h"

#include <cstdio>
#include <chrono>
#include "nvToolsExt.h"
#include <cuda_profiler_api.h>

#include "AeTime.h"
#include "ExpManager.h"

#include "cuIndividual.cuh"

using namespace std::chrono;
using namespace std;

cuExpManager::cuExpManager(const ExpManager* cpu_exp) {
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

void cuExpManager::run_evolution(int nb_gen) {
//    cudaProfilerStart();
//    high_resolution_clock::time_point t1 = high_resolution_clock::now();
//    cout << "Transfer" << endl;
//    transfer_to_device();
//    high_resolution_clock::time_point t2 = high_resolution_clock::now();
//    auto duration_transfer_in = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
//    cout << "Transfer done in " << duration_transfer_in << endl;
//
//    printf("Running evolution GPU from %d to %d\n",AeTime::time(),AeTime::time()+nb_gen);
//    bool firstGen = true;
//    for (int gen = 0; gen < nb_gen+1; gen++) {
//        if(gen == 91) nvtxRangePushA("generation 91 to 100");
//        AeTime::plusplus();
//
//        high_resolution_clock::time_point t1 = high_resolution_clock::now();
//        run_a_step_on_GPU(nb_indivs_, w_max_, selection_pressure_, grid_width_, grid_height_,mutation_rate_);
//
//        t2 = high_resolution_clock::now();
//        auto duration_transfer_in = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
//
//        std::cout<<"LOG,"<<duration_transfer_in<<std::endl;
//
//        firstGen = false;
//        if(gen == 100) nvtxRangePop();
//        printf("Generation %d : \n",AeTime::time());
//    }
//    cudaProfilerStop();
}

void cuExpManager::save(int t) {
    printf("Oups, not supported !\n");
}

void cuExpManager::load(int t) {
    printf("Oups, not supported !\n");
}

__global__
void init_device_population(int nb_indiv, int dna_length, cuIndividual* all_individuals, char* all_genomes,
                            uint8_t* all_promoters, uint* all_terminators, uint* all_prot_start) {
    auto idx = threadIdx.x + blockIdx.x * blockDim.x;
    auto rr_width = blockDim.x * gridDim.x;

    for (int i = idx; i < nb_indiv; i += rr_width) {
        auto& local_indiv = all_individuals[i];
        local_indiv.size = dna_length;
        auto offset = dna_length * i;
        local_indiv.genome = all_genomes + offset + i * PROM_SIZE;
        local_indiv.promoters = all_promoters + offset;
        local_indiv.terminators = all_terminators + offset;
        local_indiv.prot_start = all_prot_start + offset;
    }
}

void cuExpManager::transfer_to_device() {
    // Allocate memory in device world
    cudaMalloc(&(device_organisms_), nb_indivs_ * sizeof(cuIndividual));
    char* all_genomes;
    auto all_genomes_size = nb_indivs_ * dna_length_;
    // For each genome, we add a phantom space at the end.
    auto all_genomes_size_w_phantom = all_genomes_size + nb_indivs_ * PROM_SIZE;

    cudaMalloc(&(all_genomes), all_genomes_size_w_phantom * sizeof(char));

    uint8_t* all_promoters;
    uint* all_terminators;
    uint* all_prot_start;
    cudaMalloc(&(all_promoters), all_genomes_size * sizeof(uint8_t));
    cudaMalloc(&(all_terminators), all_genomes_size * sizeof(uint));
    cudaMalloc(&(all_prot_start), all_genomes_size * sizeof(uint));

    // Transfer data from individual to device
    for (int i = 0; i < nb_indivs_; ++i) {
        auto offset = dna_length_ + PROM_SIZE;
        auto indiv_genome_pointer = all_genomes + (i * offset);
        auto indiv_genome_phantom_pointer = indiv_genome_pointer + dna_length_;
        cudaMemcpy(indiv_genome_pointer, host_organisms_[i], dna_length_, cudaMemcpyHostToDevice);
        cudaMemcpy(indiv_genome_phantom_pointer, host_organisms_[i], PROM_SIZE, cudaMemcpyHostToDevice);
    }

    init_device_population<<<1, 1>>>(nb_indivs, dna_length_, all_genomes)

    // Evaluate all the individuals

    // Transfer data from prng
}
