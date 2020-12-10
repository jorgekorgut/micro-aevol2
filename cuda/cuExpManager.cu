//
// Created by elturpin on 03/12/2020.
//

#include "cuExpManager.h"

#include <cstdio>
#include <chrono>
#include "nvToolsExt.h"
#include <cuda_profiler_api.h>

#include "AeTime.h"

#include "cuIndividual.cuh"

using namespace std::chrono;
using namespace std;

#if !defined(_NDEBUG)
#define checkCuda(X) { \
    auto result = X; \
    if (result != cudaSuccess) { \
        fprintf(stderr, "CUDA Runtime Error: %s in file %s line %d\\n\n", \
                cudaGetErrorString(result), __FILE__, __LINE__); \
        assert(result == cudaSuccess); \
    } }
#else
#define checkCuda(X) X
#endif

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
    counters_ = new ctr_value_type[nb_counter_];
    for (int i = 0; i < nb_counter_; ++i) {
        counters_[i] = cpu_exp->rng_->counters()[i];
    }
}

cuExpManager::~cuExpManager() {

}

__global__
void evaluate_population(uint nb_indivs, cuIndividual* individuals, const double* target) {
    // one block per individual
    auto indiv_idx = blockIdx.x;
    if (indiv_idx < nb_indivs) {
        individuals[indiv_idx].evaluate(target);
    }
    if (indiv_idx == 0) {
        if (threadIdx.x == 0) {
            printf("Fitness: %0.10e\n", individuals[0].fitness);
        }
    }
}

void cuExpManager::run_a_step() {
    // Selection

    // Mutation

    // Evaluation
    evaluate_population<<<nb_indivs_, 32>>>(nb_indivs_, device_organisms_, device_target_);

}

void cuExpManager::run_evolution(int nb_gen) {
    cudaProfilerStart();
    cout << "Transfer" << endl;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    transfer_to_device();
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration_transfer_in = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    cout << "Transfer done in " << duration_transfer_in << " µs" << endl;

    evaluate_population<<<nb_indivs_, 32>>>(nb_indivs_, device_organisms_, device_target_);
    cudaDeviceSynchronize();
    checkCuda(cudaGetLastError());

//    printf("Running evolution GPU from %d to %d\n", AeTime::time(), AeTime::time() + nb_gen);
//    bool firstGen = true;
//    for (int gen = 0; gen < nb_gen+1; gen++) {
//        if(gen == 91) nvtxRangePushA("generation 91 to 100");
//        AeTime::plusplus();
//
//        high_resolution_clock::time_point t1 = high_resolution_clock::now();
//        run_a_step();
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
void init_device_population(int nb_indivs, int dna_length, cuIndividual* all_individuals, char* all_genomes,
                            uint8_t* all_promoters, uint* all_terminators, uint* all_prot_start, cuRNA* all_rnas) {
    auto idx = threadIdx.x + blockIdx.x * blockDim.x;
    auto rr_width = blockDim.x * gridDim.x;

    for (int i = idx; i < nb_indivs; i += rr_width) {
        auto& local_indiv = all_individuals[i];
        local_indiv.size = dna_length;
        auto offset = dna_length * i;
        local_indiv.genome = all_genomes + offset + i * PROM_SIZE;
        local_indiv.promoters = all_promoters + offset;
        local_indiv.terminators = all_terminators + offset;
        local_indiv.prot_start = all_prot_start + offset;
        local_indiv.list_rnas = all_rnas + offset;
    }
//    __syncthreads();
//    if (idx == 0) {
//        for (int i = 0; i < nb_indivs; ++i) {
//            auto& local_indiv = all_individuals[i];
//            for (int j = 0; j < 20; ++j) {
//                printf("%c", local_indiv.genome[j]);
//            }
//            printf("\n");
//        }
//    }
}

__global__
void check_rng(key_type* seed, ctr_value_type* counter, int nb_indivs) {
    auto idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx == 0) {
        printf("seed: %lu, counter[nb_indivs]: %lu\n", (*seed)[1], counter[nb_indivs]);
    }
}

__global__
void check_target(double* target) {
    for (int i = 0; i < FUZZY_SAMPLING; ++i) {
        if (target[i] == 0.0){
            printf("0|");
        } else {
            printf("%f|", target[i]);
        }
    }
    printf("\n");
}

void cuExpManager::transfer_to_device() {
    // Allocate memory in device world
    checkCuda(cudaMalloc(&(device_organisms_), nb_indivs_ * sizeof(cuIndividual)));
    char* all_genomes;
    auto all_genomes_size = nb_indivs_ * dna_length_;
    // For each genome, we add a phantom space at the end.
    auto all_genomes_size_w_phantom = all_genomes_size + nb_indivs_ * PROM_SIZE;

    checkCuda(cudaMalloc(&(all_genomes), all_genomes_size_w_phantom * sizeof(char)));

    uint8_t* all_promoters;
    uint* all_terminators;
    uint* all_prot_start;
    cuRNA* all_rnas;
    checkCuda(cudaMalloc(&(all_promoters), all_genomes_size * sizeof(uint8_t)));
    checkCuda(cudaMalloc(&(all_terminators), all_genomes_size * sizeof(uint)));
    checkCuda(cudaMalloc(&(all_prot_start), all_genomes_size * sizeof(uint)));
    checkCuda(cudaMalloc(&(all_rnas), all_genomes_size * sizeof(cuRNA)));

    // Transfer data from individual to device
    for (int i = 0; i < nb_indivs_; ++i) {
        auto offset = dna_length_ + PROM_SIZE;
        auto indiv_genome_pointer = all_genomes + (i * offset);
        auto indiv_genome_phantom_pointer = indiv_genome_pointer + dna_length_;
        checkCuda(cudaMemcpy(indiv_genome_pointer, host_organisms_[i], dna_length_, cudaMemcpyHostToDevice));
        checkCuda(cudaMemcpy(indiv_genome_phantom_pointer, host_organisms_[i], PROM_SIZE, cudaMemcpyHostToDevice));
    }

    init_device_population<<<1, 1>>>(nb_indivs_, dna_length_, device_organisms_, all_genomes,
                                     all_promoters, all_terminators, all_prot_start, all_rnas);
//    cudaDeviceSynchronize();
//    checkCuda(cudaGetLastError());

    // Transfer phenotypic target
    checkCuda(cudaMalloc(&(device_target_), FUZZY_SAMPLING * sizeof(double)));
    checkCuda(cudaMemcpy(device_target_, target_, FUZZY_SAMPLING * sizeof(double), cudaMemcpyHostToDevice));

    // Transfer data from prng
    checkCuda(cudaMalloc(&(device_rng_counters), nb_counter_ * sizeof(ctr_value_type)));
    checkCuda(cudaMemcpy(device_rng_counters, counters_, nb_counter_ * sizeof(ctr_value_type), cudaMemcpyHostToDevice));
    key_type tmp = {{0, seed_}};
    checkCuda(cudaMalloc(&(device_seed_), sizeof(key_type)));
    checkCuda(cudaMemcpy(device_seed_, &tmp, sizeof(key_type), cudaMemcpyHostToDevice));

//    check_rng<<<1, 1>>>(device_seed_, device_rng_counters, nb_indivs_);
//    check_target<<<1, 1>>>(device_target_);
    cudaDeviceSynchronize();
    checkCuda(cudaGetLastError());
}
