//
// Created by elturpin on 22/10/2020.
//

#pragma once

#include <cstdint>

#include "aevol_constants.h"
#include "cuProtein.cuh"

struct cuGene {
    // value function of the error of the RNA hosting the gene
    uint8_t concentration{};
    // position at which translation will start, after the START
    uint start{};
    // position of the terminators of the RNA
    uint length_limit{};
};

struct cuRNA {
    uint8_t errors{};
    uint start_transcription{};
    uint transcription_length{};

    uint nb_gene{};
    cuGene *list_gene{};
};

struct cuIndividual {
    __device__ void search_patterns();
    __device__ void sparse_meta();
    __device__ void transcription();
    __device__ void find_gene_per_RNA();
    __device__ void translation();


    __device__ void prepare_rnas();

    __device__ void compute_rna(uint rna_idx) const;

    __device__ void prepare_gene(uint rna_idx) const;

    __device__ void gather_genes();

    __device__ void translate_gene(uint gene_idx) const;

    __device__ void compute_phenotype();

    __device__ void compute_fitness(const double* target);

    __device__ void clean_metadata();

    inline __device__ uint get_distance(uint a, uint b) const {
        if (a > b)
            return (b + size) - a;
        return b - a;
    }

    inline __device__ uint get_distance_ori(uint a, uint b) const {
        if (a >= b)
            return (b + size) - a;
        return b - a;
    }

    // Printing
    __device__ void print_metadata_summary() const;

    __device__ void print_rnas() const;

    __device__ void print_gathered_genes() const;

    __device__ void print_proteins() const;

    __device__ void print_phenotype() const;

    uint size{};
    char *genome{};
    uint8_t *promoters{};
    uint nb_terminator{};
    uint *terminators{};
    uint nb_prot_start{};
    uint *prot_start{};

    uint nb_rnas{};
    cuRNA *list_rnas{};

    uint nb_gene{};
    cuGene *list_gene{};
    cuProtein *list_protein{};

    double phenotype[FUZZY_SAMPLING]{};
    double fitness{};
};