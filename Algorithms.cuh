#pragma once

constexpr int32_t RNA_LIST_INCR_SIZE = 200;
constexpr int32_t RNA_LIST_PROTEIN_INCR_SIZE = 250;

constexpr int32_t PROTEIN_LIST_INCR_SIZE = 400;

constexpr int8_t SELECTION_SCOPE_X = 3;
constexpr int8_t SELECTION_SCOPE_Y = 3;
constexpr int8_t NEIGHBORHOOD_SIZE = SELECTION_SCOPE_X*SELECTION_SCOPE_Y;

char *host_dna;
char *host_next_gen_dna;

double **host_phenotype;
double **host_phenotype_activ;
double **host_phenotype_inhib;

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

int32_t* protein_offset;

int32_t* protein_idx;

pProtein* protein;

int32_t* nb_proteins;
int current_size_protein_list;

// RNA structures
struct pRNA {
    int32_t begin;
    int32_t end;
    int8_t dist;
    bool transcribed;
    int32_t length;

    int32_t indiv_id;
};

int current_size_rna_list;
pRNA* rna;

struct transcribedRNA {
    int32_t* start_prot;
    int32_t max_protein_elements;
    int32_t nb_protein;
    int32_t start_lenght;
};

int32_t* rna_offset;

int32_t* rna_idx;

// DNAMutator
struct GPUDnaMutator {
    int nb_swi_;
    int nb_mut_;
    int cpt_mut_;

};

GPUDnaMutator* dna_mutator_list;

struct TypeMutation {
    int32_t type_;

    int32_t pos_1_,pos_2_,pos_3_;

    int32_t number_;

    char seq[6];

    size_t transient_size;
};

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


/**
 * Kernels
 */
__global__
// Copy first indiv's dna into all the other indivs' dna
void clone_init_indiv(size_t* dna_size, char* dna);

__global__
void search_start_stop_RNA(size_t* dna_size, char* dna, size_t* dna_offset, int* nb_promoters,
                           int8_t* dna_term, int nb_indivs, unsigned long long* nb_mut_bp);

__global__
void compute_RNA_offset(int* nb_promoters, int* rna_offset);

__global__
void fill_RNA( int8_t* dna_term, size_t* dna_size, size_t* dna_offset, int* nb_promoters, int* rna_offset, pRNA* rnas, int* rna_idx, int nb_indiv);

__global__
void compute_RNA( int8_t* dna_term, size_t* dna_size, size_t* dna_offset, pRNA* rnas,  int global_nb_rna);

__global__
void compute_start_protein(int8_t* start_protein, size_t* dna_size, size_t* dna_offset, pRNA* rna, char* dna, int32_t* nb_proteins,
                           int32_t global_nb_rna, int nb_indiv);

__global__
void compute_protein_offset(int32_t* nb_proteins, int* protein_offset);

__global__
void fill_protein(int8_t* start_protein,size_t* dna_offset, int* protein_idx, int* protein_offset, pRNA* rna, pProtein* protein,
                  size_t* dna_size, int32_t global_nb_rna, int nb_indiv);

__global__
void compute_proteins( int8_t* start_protein, size_t* dna_size, size_t* dna_offset, pProtein* protein, char* dna,
                       int32_t global_nb_protein);

__global__
void translate_proteins( pProtein* protein, size_t* dna_size, char* dna,  size_t* dna_offset, int32_t global_nb_protein, double w_max);

__global__
void compute_phenotype( pProtein* protein, int32_t global_nb_protein, double** phenotype,
                        double** phenotype_activ, double** phenotype_inhib, int nb_indiv);

__global__
void compute_metaerror_fitness(double selection_pressure,double** phenotype,
                               double** phenotype_activ,double** phenotype_inhib,
                                          double* target,
                                          double* metaerror, double* fitness);
__global__
void selection(double* fitness, int* next_generation_reproducer, unsigned long long* gpu_counters,
                          int grid_width, int grid_height, int nb_indiv);

__global__
void generate_mutations(unsigned long long* gpu_counters, size_t* dna_size, int* nb_mutations,
                        GPUDnaMutator* dna_mutator_list,int* next_generation_reproducer,
                        int nb_indivs, double mutation_rate, int nb_indiv);
__global__
void do_memset(double** phenotype_activ, double** phenotype_inhib, int* nb_mutations, int32_t* rna_idx,
               int32_t* protein_idx, int* nb_proteins, int* nb_promoters,
               size_t* dna_size,
               int nb_indiv);
