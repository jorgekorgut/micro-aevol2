//
// Created by elturpin on 14/12/2020.
//

#include "cuda_kernels.cuh"

#include <cstdio>

// Evolution

__global__
void check_result(uint nb_indivs, cuIndividual* individuals) {
    for (int indiv_idx = 0; indiv_idx < nb_indivs; ++indiv_idx) {
        const auto& indiv = individuals[indiv_idx];
        printf("%d: %1.10e | ", indiv_idx, indiv.fitness);
    }
    printf("\n");
}

__global__
void evaluate_population(uint nb_indivs, cuIndividual* individuals, const double* target) {
    // One block per individual
    auto indiv_idx = blockIdx.x;
    if (indiv_idx < nb_indivs) {
        individuals[indiv_idx].evaluate(target);
    }
}

__global__
void do_mutation(uint nb_indivs, cuIndividual* individuals, double mutation_rate, RandService* rand_service) {
    // One thread per individual
    auto indiv_idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (indiv_idx >= nb_indivs)
        return;

    auto& indiv = individuals[indiv_idx];
    auto nb_switch = rand_service->binomial_random(indiv.size, mutation_rate, indiv_idx, MUTATION);

    for (int i = 0; i < nb_switch; ++i) {
        auto position = rand_service->gen_number_max(indiv.size, indiv_idx, MUTATION);
        if (indiv.genome[position] == '0') indiv.genome[position] = '1';
        else indiv.genome[position] = '0';

        if (position < PROM_SIZE) // Do not forget phantom space
            indiv.genome[position + indiv.size] = indiv.genome[position];

//        printf("indiv %u switch at %d\n", indiv_idx, position);
    }
}

__global__
void reproduction(uint nb_indivs, cuIndividual* individuals, const int* reproducers, const char* all_parent_genome) {
    // One block per individuals
    auto indiv_idx = blockIdx.x;
    auto idx = threadIdx.x;
    auto rr_width = blockDim.x;

    if (indiv_idx >= nb_indivs) {
        return;
    }

    __shared__ const char* parent_genome;
    __shared__ char* child_genome;
    __shared__ uint size;

    if (threadIdx.x == 0) {
        auto& child = individuals[indiv_idx];
        // size do not change, what a chance !
        size = child.size;
        auto offset = reproducers[indiv_idx] * (child.size + PROM_SIZE); // Do not forget phantom space
        parent_genome = all_parent_genome + offset;
        child_genome = child.genome;
    }
    __syncthreads();

    for (int position = idx; position < size + PROM_SIZE; position += rr_width) { // Do not forget phantom space
        child_genome[position] = parent_genome[position];
    }
}

__global__
void selection(uint grid_height, uint grid_width, const cuIndividual* individuals, RandService* rand_service,
          int* next_reproducers) {
    // One thread per grid cell
    int grid_x = threadIdx.x + blockIdx.x * blockDim.x;
    int grid_y = threadIdx.y + blockIdx.y * blockDim.y;

    if (grid_x >= grid_width || grid_y >= grid_height)
        return;

    double local_fit_array[NEIGHBORHOOD_SIZE];
    int count = 0;
    double sum_local_fit = 0.0;

    for (int8_t i = -1; i < NEIGHBORHOOD_WIDTH - 1; i++) {
        for (int8_t j = -1; j < NEIGHBORHOOD_HEIGHT - 1; j++) {
            // Toric topology
            int cur_x = (grid_x + i + grid_width) % grid_width;
            int cur_y = (grid_y + j + grid_height) % grid_height;

            local_fit_array[count] = individuals[cur_x * grid_width + cur_y].fitness;
            sum_local_fit += local_fit_array[count];

            count++;
        }
    }

    for(int8_t i = 0 ; i < NEIGHBORHOOD_SIZE ; i++) {
        local_fit_array[i] /= sum_local_fit;
    }

    uint grid_idx = grid_x * grid_width + grid_y;

    auto selected_cell = rand_service->random_roulette(local_fit_array, NEIGHBORHOOD_SIZE, grid_idx, SELECTION);

    int x_offset = (selected_cell / NEIGHBORHOOD_WIDTH) - 1;
    int y_offset = (selected_cell % NEIGHBORHOOD_HEIGHT) - 1;
    int selected_x = (grid_x + x_offset + grid_width) % grid_width;
    int selected_y = (grid_y + y_offset + grid_height) % grid_height;

    next_reproducers[grid_idx] = selected_x * grid_height + selected_y;
}

__global__
void swap_parent_child_genome(uint nb_indivs, cuIndividual* individuals, char* all_parent_genome) {
    // One thread per individual
    auto indiv_idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (indiv_idx >= nb_indivs)
        return;

    auto& indiv = individuals[indiv_idx];
    auto offset = indiv_idx * (indiv.size + PROM_SIZE); // Do not forget phantom space
    indiv.genome = all_parent_genome + offset;
}

// Interface Host | Device

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

__global__
void check_rng(RandService* rand_service) {
    auto idx = threadIdx.x + blockIdx.x * blockDim.x;
    auto number = rand_service->generator({{0, rand_service->rng_counters[rand_service->phase_size]}}, rand_service->seed);
    if (idx == 0) {
        printf("seed: %lu, counter[nb_indivs]: %lu, test: %lu\n",
               rand_service->seed[1],
               rand_service->rng_counters[rand_service->phase_size],
               number[0]);
    }
}

__global__
void init_device_population(int nb_indivs, int genome_length, cuIndividual* all_individuals, char* all_genomes,
                            uint8_t* all_promoters, uint* all_terminators, uint* all_prot_start, cuRNA* all_rnas) {
    auto idx = threadIdx.x + blockIdx.x * blockDim.x;
    auto rr_width = blockDim.x * gridDim.x;

    for (int i = idx; i < nb_indivs; i += rr_width) {
        auto& local_indiv = all_individuals[i];
        local_indiv.size = genome_length;
        auto offset = genome_length * i;
        local_indiv.genome = all_genomes + offset + i * PROM_SIZE;
        local_indiv.promoters = all_promoters + offset;
        local_indiv.terminators = all_terminators + offset;
        local_indiv.prot_start = all_prot_start + offset;
        local_indiv.list_rnas = all_rnas + offset;
        local_indiv.nb_terminator = 0;
        local_indiv.nb_prot_start = 0;
        local_indiv.nb_rnas = 0;
        local_indiv.nb_gene = 0;
        local_indiv.list_gene = nullptr;
        local_indiv.list_protein = nullptr;
        local_indiv.fitness = 0.0;
        for (int j = 0; j < FUZZY_SAMPLING; ++j) {
            local_indiv.phenotype[j] = 0.0;
        }
    }
}
