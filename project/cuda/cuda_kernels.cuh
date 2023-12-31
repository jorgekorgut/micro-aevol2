//
// Created by elturpin on 14/12/2020.
//

#pragma once

#include "cuIndividual.cuh"
#include "RandService.h"
#include "../aevol_constants.h"

// Evolution

__global__
void selection(uint grid_height, uint grid_width, const cuIndividual* individuals, RandService* rand_service,
               int* next_reproducers);

__global__
void reproduction(uint nb_indivs, cuIndividual* individuals, const int* reproducers, const block* all_parent_genome);

__global__
void do_mutation(uint nb_indivs, cuIndividual* individuals, double mutation_rate, RandService* rand_service);


__global__
void swap_parent_child_genome(uint nb_indivs, cuIndividual* individuals, uint64_t* all_parent_genome);

__global__
void check_result(uint nb_indivs, cuIndividual* individuals);

// Evaluation process
__global__
void clean_metadata(uint nb_indivs, cuIndividual* individual);

__global__
void search_patterns(uint nb_indivs, cuIndividual* individual);

__global__
void sparse_meta(uint nb_indivs, cuIndividual* individual);

__global__
void transcription(uint nb_indivs, cuIndividual* individual);

__global__
void find_gene_per_RNA(uint nb_indivs, cuIndividual* individual);

__global__
void gather_genes(uint nb_indivs, cuIndividual* individual);

__global__
void translation(uint nb_indivs, cuIndividual* individual);

__global__
void compute_phenotype(uint nb_indivs, cuIndividual* individual);

__global__
void compute_fitness(uint nb_indivs, cuIndividual* individual, const double* target);

// Interface Host | Device

__global__
void clean_population_metadata(uint nb_indivs, cuIndividual* individuals);

__global__
void init_device_population(int nb_indivs, uint block_length,
                            uint block_length_phantom, int genome_length,
                            cuIndividual* all_individuals, block* all_genomes,
                            uint8_t* all_promoters, block* all_terminators,
                            uint* all_terms_idx, block* all_prot_start,
                            uint* all_prots_idxs, cuRNA* all_rnas);

__global__
void check_rng(RandService* rand_service);

__global__
void check_target(double* target);


