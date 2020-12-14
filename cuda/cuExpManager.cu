//
// Created by elturpin on 03/12/2020.
//

#include "cuExpManager.h"

#include <cstdio>
#include <cassert>
#include <chrono>
#include <cmath>
#include "nvToolsExt.h"
#include <cuda_profiler_api.h>
#include <zlib.h>

#include "AeTime.h"

#include "cuIndividual.cuh"
#include "cuda_kernels.cuh"

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

#define CHECK_KERNEL \
cudaDeviceSynchronize(); \
checkCuda(cudaGetLastError());
#else
#define checkCuda(X) X
#define CHECK_KERNEL
#endif

cuExpManager::cuExpManager(const ExpManager* cpu_exp) {
    grid_height_ = cpu_exp->grid_height_;
    grid_width_ = cpu_exp->grid_height_;

    mutation_rate_ = cpu_exp->mutation_rate_;

    backup_step_ = cpu_exp->backup_step_;

    nb_indivs_ = grid_height_ * grid_width_;

    dna_length_ = cpu_exp->internal_organisms_[0]->length();
    host_organisms_ = new char *[nb_indivs_];
    for (int i = 0; i < nb_indivs_; ++i) {
        host_organisms_[i] = new char[dna_length_];
        const auto& org = cpu_exp->internal_organisms_[i];
        memcpy(host_organisms_[i], org->dna_->seq_.data(), dna_length_ * sizeof(char));
    }

    target_ = new double[FUZZY_SAMPLING];
    memcpy(target_, cpu_exp->target, FUZZY_SAMPLING * sizeof(double));

    seed_ = cpu_exp->seed_;
    nb_counter_ = cpu_exp->rng_->counters().size();
    counters_ = new ctr_value_type[nb_counter_];
    for (int i = 0; i < nb_counter_; ++i) {
        counters_[i] = cpu_exp->rng_->counters()[i];
    }
}

cuExpManager::~cuExpManager() {

}

void cuExpManager::run_a_step() {
    auto threads_per_block = 64; // arbitrary : better if multiple of 32
    // Selection
    dim3 bloc_dim(threads_per_block / 2, threads_per_block / 2);
    auto grid_x = ceil((float) grid_width_ / (float) bloc_dim.x);
    auto grid_y = ceil((float) grid_height_ / (float) bloc_dim.y);
    dim3 grid_dim(grid_x, grid_y);
    selection<<<grid_dim, bloc_dim>>>(grid_height_,
                                      grid_width_,
                                      device_organisms_,
                                      rand_service_,
                                      reproducers_);

    // Reproduction
    reproduction<<<nb_indivs_, threads_per_block>>>(nb_indivs_, device_organisms_, reproducers_, all_parent_dna_);

    // Mutation
    auto grid_dim_1d = ceil((float)nb_indivs_ / (float)threads_per_block);
    do_mutation<<<grid_dim_1d, threads_per_block>>>(nb_indivs_, device_organisms_, mutation_rate_, rand_service_);

    // Evaluation
    evaluate_population<<<nb_indivs_, threads_per_block>>>(nb_indivs_, device_organisms_, device_target_);

    // Swap genome information
    swap_parent_child_genome<<<grid_dim_1d, threads_per_block>>>(nb_indivs_, device_organisms_, all_parent_dna_);
    swap(all_parent_dna_, all_child_dna_);
}

void cuExpManager::run_evolution(int nb_gen) {
    cudaProfilerStart();
    cout << "Transfer" << endl;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    transfer_to_device();
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration_transfer_in = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    cout << "Transfer done in " << duration_transfer_in << " µs" << endl;

    // Evaluation of population at generation 0
    auto threads_per_block = 64;
    auto grid_dim_1d = ceil((float) nb_indivs_ / (float) threads_per_block);
    evaluate_population<<<nb_indivs_, threads_per_block>>>(nb_indivs_, device_organisms_, device_target_);
    swap_parent_child_genome<<<grid_dim_1d, threads_per_block>>>(nb_indivs_, device_organisms_, all_parent_dna_);
    CHECK_KERNEL
    swap(all_parent_dna_, all_child_dna_);
    check_result<<<1,1>>>(nb_indivs_, device_organisms_);
    CHECK_KERNEL

    printf("Running evolution GPU from %d to %d\n", AeTime::time(), AeTime::time() + nb_gen);
    for (int gen = 0; gen < nb_gen; gen++) {
        AeTime::plusplus();
        printf("Generation %d : \n",AeTime::time());

        run_a_step();

        if (AeTime::time() % backup_step_ == 0) {
            save(AeTime::time());
            cout << "Backup for generation " << AeTime::time() << " done !" << endl;
        }
    }

    check_result<<<1,1>>>(nb_indivs_, device_organisms_);
    CHECK_KERNEL
    cudaProfilerStop();
}

void cuExpManager::save(int t) const {
    char exp_backup_file_name[255];
    sprintf(exp_backup_file_name, "backup/backup_%d.zae", t);

    // -------------------------------------------------------------------------
    // Open backup files
    // -------------------------------------------------------------------------
    gzFile exp_backup_file = gzopen(exp_backup_file_name, "w");

    // -------------------------------------------------------------------------
    // Check that files were correctly opened
    // -------------------------------------------------------------------------
    if (exp_backup_file == Z_NULL) {
        printf("Error: could not open backup file %s\n",
               exp_backup_file_name);
        exit(EXIT_FAILURE);
    }

    // -------------------------------------------------------------------------
    // Write the backup file
    // -------------------------------------------------------------------------
    gzwrite(exp_backup_file, &t, sizeof(t));

    gzwrite(exp_backup_file, &grid_height_, sizeof(grid_height_));
    gzwrite(exp_backup_file, &grid_width_, sizeof(grid_width_));
    gzwrite(exp_backup_file, &backup_step_, sizeof(backup_step_));
    gzwrite(exp_backup_file, &mutation_rate_, sizeof(mutation_rate_));


    for (int i = 0; i < FUZZY_SAMPLING; i++) {
        double tmp = target_[i];
        gzwrite(exp_backup_file, &tmp, sizeof(tmp));
    }

    transfer_to_host();

    for (int indiv_id = 0; indiv_id < nb_indivs_; indiv_id++) {
        gzwrite(exp_backup_file, &dna_length_, sizeof(dna_length_));
        gzwrite(exp_backup_file, host_organisms_[indiv_id], dna_length_ * sizeof(char));
    }

    gzwrite(exp_backup_file, &seed_, sizeof(seed_));
    gzwrite(exp_backup_file, counters_, nb_counter_ * sizeof(counters_[0]));

    if (gzclose(exp_backup_file) != Z_OK) {
        cerr << "Error while closing backup file" << endl;
    }
}

void cuExpManager::load(int t) {
    // Unused method because we load data from ExpManager
}


void cuExpManager::transfer_to_device() {
    // Allocate memory for individuals in device world
    checkCuda(cudaMalloc(&(device_organisms_), nb_indivs_ * sizeof(cuIndividual)));
    auto all_genomes_size = nb_indivs_ * dna_length_;
    // For each genome, we add a phantom space at the end.
    auto all_genomes_size_w_phantom = all_genomes_size + nb_indivs_ * PROM_SIZE;

    checkCuda(cudaMalloc(&(all_child_dna_), all_genomes_size_w_phantom * sizeof(char)));
    checkCuda(cudaMalloc(&(all_parent_dna_), all_genomes_size_w_phantom * sizeof(char)));

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
        auto indiv_genome_pointer = all_child_dna_ + (i * offset);
        auto indiv_genome_phantom_pointer = indiv_genome_pointer + dna_length_;
        checkCuda(cudaMemcpy(indiv_genome_pointer, host_organisms_[i], dna_length_, cudaMemcpyHostToDevice));
        checkCuda(cudaMemcpy(indiv_genome_phantom_pointer, host_organisms_[i], PROM_SIZE, cudaMemcpyHostToDevice));
    }

    init_device_population<<<1, 1>>>(nb_indivs_, dna_length_, device_organisms_, all_child_dna_,
                                     all_promoters, all_terminators, all_prot_start, all_rnas);
    CHECK_KERNEL

    // Transfer phenotypic target
    checkCuda(cudaMalloc(&(device_target_), FUZZY_SAMPLING * sizeof(double)));
    checkCuda(cudaMemcpy(device_target_, target_, FUZZY_SAMPLING * sizeof(double), cudaMemcpyHostToDevice));

    // Allocate memory for reproduction data
    checkCuda(cudaMalloc(&(reproducers_), nb_indivs_ * sizeof(int)));

    // Initiate Random Number generator
    RandService tmp;
    checkCuda(cudaMalloc(&(tmp.rng_counters), nb_counter_ * sizeof(ctr_value_type)));
    checkCuda(cudaMemcpy(tmp.rng_counters, counters_, nb_counter_ * sizeof(ctr_value_type), cudaMemcpyHostToDevice));
    tmp.seed = {{0, seed_}};
    tmp.phase_size = nb_indivs_;
    assert(nb_counter_ == tmp.phase_size * NPHASES);

    checkCuda(cudaMalloc(&(rand_service_), sizeof(RandService)));
    checkCuda(cudaMemcpy(rand_service_, &tmp, sizeof(RandService), cudaMemcpyHostToDevice));

//    check_rng<<<1, 1>>>(rand_service_);
//    check_target<<<1, 1>>>(device_target_);
//    CHECK_KERNEL
}

void cuExpManager::transfer_to_host() const {
    for (int i = 0; i < nb_indivs_; ++i) {
        auto offset = dna_length_ + PROM_SIZE;
        auto indiv_genome_pointer = all_parent_dna_ + (i * offset);
        checkCuda(cudaMemcpy(host_organisms_[i], indiv_genome_pointer, dna_length_, cudaMemcpyDeviceToHost));
    }

    RandService tmp;
    checkCuda(cudaMemcpy(&tmp, rand_service_, sizeof(RandService), cudaMemcpyDeviceToHost));
    checkCuda(cudaMemcpy(counters_, tmp.rng_counters, nb_counter_ * sizeof(ctr_value_type), cudaMemcpyDeviceToHost));
}
