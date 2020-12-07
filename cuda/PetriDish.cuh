//
// Created by elturpin on 03/12/2020.
//

#pragma once

constexpr int8_t SELECTION_SCOPE_X = 3;
constexpr int8_t SELECTION_SCOPE_Y = 3;
constexpr int8_t NEIGHBORHOOD_SIZE = SELECTION_SCOPE_X*SELECTION_SCOPE_Y;

// RNA structures
struct pRNA {
    int32_t begin;
    int32_t end;
    int8_t dist;
    bool transcribed;
    int32_t length;

    int32_t indiv_id;
};

// Protein structures
struct pProtein {
    int32_t protein_start;
    int32_t protein_end;
    int32_t protein_length;
    double m;
    double w;
    double h;
    double e;
    bool is_functional;

    int32_t indiv_id;
    int32_t stop_RNA;
    bool translated;
};

struct TypeMutation {
    int32_t type_;

    int32_t pos_1_,pos_2_,pos_3_;

    int32_t number_;

    char seq[6];

    size_t transient_size;
};

// DNAMutator
struct GPUDnaMutator {
    int nb_swi_;
    int nb_mut_;
    int cpt_mut_;

};

class PetriDish {
    double **host_phenotype;
    double **host_phenotype_activ;
    double **host_phenotype_inhib;

    int32_t* protein_offset;

    int32_t* protein_idx;

    pProtein* protein;

    int32_t* nb_proteins;
    int current_size_protein_list;

    int current_size_rna_list;
    pRNA* rna;

    int32_t* rna_offset;

    int32_t* rna_idx;

    GPUDnaMutator* dna_mutator_list;

    TypeMutation* tab_mutation;

    int* nb_mutations;
    int* mutations_idx;
    int* mutations_offset;

    int current_size_tab_mutation;


/**
 * Structure to transfer from host to device memory
 */
// All DNAs
    char* dna;
    char* next_gen_dna;

    size_t global_dna_size;
    size_t allocated_global_dna_size;

// All DNA size
    size_t* dna_size;
    size_t* next_gen_dna_size;


    size_t* dna_offset;
    size_t* next_gen_dna_offset;

// Current maximum DNA size
    size_t* max_dna_size;
    size_t host_max_dna_size;

    unsigned long long int* nb_mut_bp;

// Terminator
    int8_t* dna_term;

// Number of promoters
    int* nb_promoters;

// Start protein
    int8_t* start_protein;

// Fuzzy structures
    double** phenotype;
    double** phenotype_activ;
    double** phenotype_inhib;
    double** delta;

// Environment (static first)
    double* target;

    int* next_generation_reproducer;

// PRNG
    unsigned long long* gpu_counters;

//
/**
 * Structure to transfer from device to host
 */
    double* metaerror;
    double* fitness;

    void run_a_step_on_GPU(int nb_indiv, double w_max, double selection_pressure, int grid_width, int grid_height, double mutation_rate);
};


