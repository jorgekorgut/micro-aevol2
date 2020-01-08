#include "Algorithms.h"
#include "Algorithms.cuh"

#include "ExpManager.h"
#include "ThreefryGPU.h"
#include "GPUDna.cuh"

#include <cstdint>
#include <stdio.h>
#include <unistd.h>

#include <iostream>

#include<cuda.h>
#include<cuda_profiler_api.h>

using namespace std;

#define DEBUG 1
// Convenience function for checking CUDA runtime API results
// can be wrapped around any runtime API call. No-op in release builds.
inline
cudaError_t checkCuda(cudaError_t result)
{
#if defined(DEBUG) || defined(_DEBUG)
    if (result != cudaSuccess) {
        fprintf(stderr, "CUDA Runtime Error: %s\n",
                cudaGetErrorString(result));
        assert(result == cudaSuccess);
    }
#endif
    return result;
}


constexpr int32_t PROMOTER_ARRAY_SIZE = 10000;

void transfer_in(ExpManager* exp_m, bool first_gen) {
    exp_m->rng_->initDevice();

  std::vector<size_t> host_dna_size(exp_m->nb_indivs_);
  std::vector<size_t> host_dna_offset(exp_m->nb_indivs_);

  // Compute sizes:
  // * global_dna_size
  // * host_dna_offset[]
  // * host_max_dna_size
  // * host_dna_size[]
  global_dna_size = 0;
  for (int i = 0; i < exp_m->nb_indivs_; i++) {
    host_dna_offset[i] = global_dna_size;
    global_dna_size += exp_m->internal_organisms_[i]->dna_->seq_.size();
    host_max_dna_size =
        host_max_dna_size < exp_m->internal_organisms_[i]->dna_->seq_.size() ?
        exp_m->internal_organisms_[i]->dna_->seq_.size() : host_max_dna_size;
    host_dna_size[i] = exp_m->internal_organisms_[i]->dna_->seq_.size();
  }

  // Create shorthands
  auto seq0 = exp_m->internal_organisms_[0]->dna_->seq_.data();
  auto len0 = exp_m->internal_organisms_[0]->dna_->seq_.size();

    allocated_global_dna_size = global_dna_size*5;

  // Allocate mem for the meta dna

    checkCuda(cudaMalloc((void **) &next_gen_dna, allocated_global_dna_size * sizeof(char)));
  checkCuda(cudaMalloc((void**) &dna, allocated_global_dna_size * sizeof(char)));
  // Tranfer **the first** indiv's sequence
  checkCuda(cudaMemcpy(dna,
                       seq0,
                       len0 * sizeof(char),
                       cudaMemcpyHostToDevice));

  // Send dna_size array
  checkCuda(cudaMalloc((void**) &dna_size,
                       exp_m->nb_indivs_ * sizeof(size_t)));
  checkCuda(cudaMemcpy(dna_size,
                       host_dna_size.data(), exp_m->nb_indivs_ * sizeof(size_t),
                       cudaMemcpyHostToDevice));

        checkCuda(cudaMalloc((void **) &nb_mut_bp, 1 * sizeof(unsigned long long int)));
        checkCuda(cudaMemset(nb_mut_bp, 0, 1 * sizeof(unsigned long long int)));


  // Launch kernel to clone initial genome into the whole pop
  int x_dim_size = (len0 / 128)+1;
  int y_dim_size = exp_m->nb_indivs_;
  dim3 dimGrid(x_dim_size,y_dim_size);

  clone_init_indiv<<<dimGrid,128>>>(dna_size, dna);

  checkCuda(cudaMalloc((void**) &dna_term, allocated_global_dna_size * sizeof(int8_t*)));

  checkCuda(cudaMalloc((void**) &start_protein,
                       allocated_global_dna_size * sizeof(int8_t*)));


  checkCuda(cudaMalloc((void**) &dna_offset,
                       exp_m->nb_indivs_ * sizeof(size_t)));
  checkCuda(cudaMemcpy(dna_offset,
                       host_dna_offset.data(),
                       exp_m->nb_indivs_ * sizeof(size_t),
                       cudaMemcpyHostToDevice));

  checkCuda(cudaMalloc((void**) &next_gen_dna_offset,
                       exp_m->nb_indivs_ * sizeof(size_t)));


  checkCuda(cudaMalloc((void**) &next_gen_dna_size,
                       exp_m->nb_indivs_ * sizeof(size_t)));

  checkCuda(cudaMalloc((void**) &nb_mutations,
                       (exp_m->nb_indivs_ + 1) * sizeof(int)));
  checkCuda(cudaMemset(nb_mutations, 0, (exp_m->nb_indivs_ + 1) * sizeof(int)));

  checkCuda(cudaMalloc((void**) &mutations_offset,
                       exp_m->nb_indivs_ * sizeof(int)));
  checkCuda(cudaMemset(mutations_offset, 0, exp_m->nb_indivs_ * sizeof(int)));

  checkCuda(cudaMalloc((void**) &mutations_idx,
                       exp_m->nb_indivs_ * sizeof(int)));
  checkCuda(cudaMemset(mutations_idx, 0, exp_m->nb_indivs_ * sizeof(int)));

  checkCuda(cudaMalloc((void**) &dna_mutator_list,
                       exp_m->nb_indivs_ * sizeof(GPUDnaMutator)));

  current_size_tab_mutation = exp_m->nb_indivs_ * 100;
  checkCuda(cudaMalloc(&tab_mutation,
                       current_size_tab_mutation * sizeof(TypeMutation)));

  checkCuda(cudaMalloc((void**) &rna_idx,
                       (exp_m->nb_indivs_ + 1) * sizeof(int32_t)));
  checkCuda(cudaMemset(rna_idx, 0, (exp_m->nb_indivs_ + 1) * sizeof(int32_t)));

  checkCuda(cudaMalloc((void**) &rna_offset,
                       exp_m->nb_indivs_ * sizeof(int32_t)));
  checkCuda(cudaMemset(rna_offset, 0, exp_m->nb_indivs_ * sizeof(int32_t)));

  checkCuda(cudaMalloc((void**) &protein_idx,
                       (exp_m->nb_indivs_ + 1) * sizeof(int32_t)));
  checkCuda(
      cudaMemset(protein_idx, 0, (exp_m->nb_indivs_ + 1) * sizeof(int32_t)));

  checkCuda(cudaMalloc((void**) &protein_offset,
                       exp_m->nb_indivs_ * sizeof(int32_t)));
  checkCuda(cudaMemset(protein_offset, 0, exp_m->nb_indivs_ * sizeof(int32_t)));

  checkCuda(cudaMalloc((void**) &next_generation_reproducer,
                       exp_m->nb_indivs_ * sizeof(size_t)));

  checkCuda(cudaMalloc((void**) &nb_promoters,
                       (exp_m->nb_indivs_ + 1) * sizeof(int)));
  checkCuda(cudaMemset(nb_promoters, 0, (exp_m->nb_indivs_ + 1) * sizeof(int)));

  checkCuda(cudaMalloc((void**) &nb_proteins,
                       (exp_m->nb_indivs_ + 1) * sizeof(int)));
  checkCuda(cudaMemset(nb_proteins, 0, (exp_m->nb_indivs_ + 1) * sizeof(int)));

  host_phenotype = (double**) malloc(exp_m->nb_indivs_ * sizeof(double*));
  checkCuda(
      cudaMalloc((void***) &phenotype, exp_m->nb_indivs_ * sizeof(double*)));


  host_phenotype_activ = (double**) malloc(exp_m->nb_indivs_ * sizeof(double*));
  checkCuda(cudaMalloc((void***) &phenotype_activ,
                       exp_m->nb_indivs_ * sizeof(double*)));


  host_phenotype_inhib = (double**) malloc(exp_m->nb_indivs_ * sizeof(double*));
  checkCuda(cudaMalloc((void***) &phenotype_inhib,
                       exp_m->nb_indivs_ * sizeof(double*)));

  for (int indiv_id = 0; indiv_id < exp_m->nb_indivs_; indiv_id++) {
    checkCuda(
        cudaMalloc((void**) &host_phenotype[indiv_id], 300 * sizeof(double)));
    checkCuda(cudaMemset(host_phenotype[indiv_id], 0.0, 300 * sizeof(double)));

    checkCuda(cudaMalloc((void**) &host_phenotype_activ[indiv_id],
                         300 * sizeof(double)));
    checkCuda(
        cudaMemset(host_phenotype_activ[indiv_id], 0.0, 300 * sizeof(double)));

    checkCuda(cudaMalloc((void**) &host_phenotype_inhib[indiv_id],
                         300 * sizeof(double)));
    checkCuda(
        cudaMemset(host_phenotype_inhib[indiv_id], 0.0, 300 * sizeof(double)));
  }

  current_size_rna_list = exp_m->nb_indivs_ * 10000;
  checkCuda(cudaMalloc(&rna, current_size_rna_list * sizeof(pRNA)));

  current_size_protein_list = exp_m->nb_indivs_ * 1000;
  checkCuda(cudaMalloc(&protein, current_size_protein_list * sizeof(pProtein)));

  checkCuda(
      cudaMemcpy(phenotype, host_phenotype, exp_m->nb_indivs_ * sizeof(double*),
                 cudaMemcpyHostToDevice));

  checkCuda(cudaMemcpy(phenotype_activ, host_phenotype_activ,
                       exp_m->nb_indivs_ * sizeof(double*),
                       cudaMemcpyHostToDevice));

  checkCuda(cudaMemcpy(phenotype_inhib, host_phenotype_inhib,
                       exp_m->nb_indivs_ * sizeof(double*),
                       cudaMemcpyHostToDevice));

  checkCuda(cudaMalloc((void**) &target,
                       300 * sizeof(double)));

  double target_host[300];
  for (int i = 0; i < 300; i++) {
    target_host[i] = exp_m->target[i];
  }


  checkCuda(cudaMemcpy(target,
                       target_host,
                       300 * sizeof(double), cudaMemcpyHostToDevice));

  checkCuda(cudaMalloc((void**) &metaerror,
                       exp_m->nb_indivs_ * sizeof(double)));


  checkCuda(cudaMalloc((void**) &fitness,
                       exp_m->nb_indivs_ * sizeof(double)));


  //printf("GPU Counter %d\n",exp_m->rng_->counters().size());

  checkCuda(cudaMalloc((void**) &gpu_counters,
                       exp_m->rng_->counters().size() *
                       sizeof(unsigned long long)));

  checkCuda(cudaMemcpy(gpu_counters, exp_m->rng_->counters().data(),
                       exp_m->rng_->counters().size() *
                       sizeof(unsigned long long), cudaMemcpyHostToDevice));

}

__global__
// Copy first indiv's dna into all the other indivs' dna
void clone_init_indiv(size_t* dna_size, char* dna) {
  int dna_chunk_idx = blockIdx.x;
  int indiv_id = blockIdx.y;
  if(indiv_id == 0) return; // don't copy indiv 0 onto itself

  int pos = (dna_chunk_idx*128)+threadIdx.x;

  if (pos < dna_size[0]) {
    dna[indiv_id*dna_size[0] + pos] = dna[pos];
  }
}


__global__
void search_start_stop_RNA(size_t* dna_size, char* dna, size_t* dna_offset, int* nb_promoters,
                           int8_t* dna_term, int nb_indivs, int global_dna_size, unsigned long long* nb_mut_bp) {

    int dna_pos_block = blockIdx.x;
    int indiv_id = blockIdx.y;

    int dna_pos = (dna_pos_block*128)+threadIdx.x;

    __shared__ int nb_prom_block;
    if (threadIdx.x == 0) {
        nb_prom_block = 0;
    }
    __syncthreads();

    if (dna_pos < dna_size[indiv_id] && dna_size[indiv_id] >= PROM_SIZE) {
        dna_term[dna_offset[indiv_id]+dna_pos] = 22;
        //atomicAdd(nb_mut_bp,1);

        int prom_dist[22];
        int term_dist[4];

        for (int motif_id = 0; motif_id < 26; motif_id++) {
            if (motif_id < 22) {
                prom_dist[motif_id] =
                        PROM_SEQ[motif_id] ==
                        dna[dna_pos + motif_id >= dna_size[indiv_id] ? dna_offset[indiv_id]+ dna_pos + motif_id - dna_size[indiv_id]
                                                                     : dna_offset[indiv_id]+ dna_pos + motif_id]
                        ? 0
                        : 1;
            } else if (motif_id >= 22) {
                int t_motif_id = motif_id - 22;
                term_dist[t_motif_id] =
                        dna[dna_pos + t_motif_id >= dna_size[indiv_id] ?
                            dna_offset[indiv_id]+dna_pos + t_motif_id - dna_size[indiv_id] :
                            dna_offset[indiv_id]+ dna_pos + t_motif_id] !=
                        dna[dna_pos - t_motif_id + 10 >= dna_size[indiv_id] ?
                            dna_offset[indiv_id]+ dna_pos - t_motif_id + 10 - dna_size[indiv_id] :
                            dna_offset[indiv_id]+ dna_pos - t_motif_id + 10] ? 1 : 0;
            }
        }

        int8_t dist_prom = prom_dist[0] +
                        prom_dist[1] +
                        prom_dist[2] +
                        prom_dist[3] +
                        prom_dist[4] +
                        prom_dist[5] +
                        prom_dist[6] +
                        prom_dist[7] +
                        prom_dist[8] +
                        prom_dist[9] +
                        prom_dist[10] +
                        prom_dist[11] +
                        prom_dist[12] +
                        prom_dist[13] +
                        prom_dist[14] +
                        prom_dist[15] +
                        prom_dist[16] +
                        prom_dist[17] +
                        prom_dist[18] +
                        prom_dist[19] +
                        prom_dist[20] +
                        prom_dist[21];



        dna_term[dna_offset[indiv_id]+dna_pos] = dist_prom;


        if (dist_prom <= 4) {
            int rna_idx = atomicAdd(&nb_prom_block, 1);
        }


        int dist_term = term_dist[0] +
                        term_dist[1] +
                        term_dist[2] +
                        term_dist[3];
        dna_term[dna_offset[indiv_id]+dna_pos] |= dist_term == 4 ? 1<<7 : 0;

    }

    __syncthreads();

    if (threadIdx.x == 0) {
        atomicAdd(nb_promoters+indiv_id,nb_prom_block);
        atomicAdd(nb_promoters+nb_indivs,nb_prom_block);
    }
}



__global__
void compute_RNA_offset(int* nb_promoters, int* rna_offset) {

    const int indiv_id = blockIdx.x;
    __shared__ int grid_rna_offset;

    if (threadIdx.x == 0) {
        grid_rna_offset = 0;
    }
    __syncthreads();

    {
        int local_rna_offset = 0;
        for (int cpt = threadIdx.x; cpt < indiv_id; cpt += blockDim.x) {
            local_rna_offset += nb_promoters[cpt];
        }

        if (local_rna_offset > 0)
            atomicAdd(&grid_rna_offset, local_rna_offset);
    }
    __syncthreads();
    if (threadIdx.x == 0) {
        rna_offset[indiv_id] = grid_rna_offset;
    }
}



__global__
void fill_RNA( int8_t* dna_term, size_t* dna_size, size_t* dna_offset, int* nb_promoters, int* rna_offset, pRNA* rnas,
               int32_t* rna_idx, int nb_indiv) {
    int dna_pos_block = blockIdx.x;
    int indiv_id = blockIdx.y;

    int dna_pos = (dna_pos_block * 128) + threadIdx.x;

    if (dna_pos < dna_size[indiv_id] && dna_size[indiv_id] >= PROM_SIZE) {
        // Masque le bit de poid fort
        int8_t dist = dna_term[dna_offset[indiv_id]+dna_pos] & (0x7F);


        if (dist <= 4) {
            int local_rna_idx = atomicAdd(rna_idx + indiv_id, 1);
            atomicAdd(rna_idx + nb_indiv, 1);

            rnas[rna_offset[indiv_id] + local_rna_idx].begin = dna_pos;
            rnas[rna_offset[indiv_id] + local_rna_idx].dist = dist;
            rnas[rna_offset[indiv_id] + local_rna_idx].transcribed = false;
            rnas[rna_offset[indiv_id] + local_rna_idx].indiv_id = indiv_id;
        }
    }
}

__global__
void compute_RNA( int8_t* dna_term, size_t* dna_size, size_t* dna_offset, pRNA* rnas,  int global_nb_rna) {
	const int globalIdx = blockIdx.x*blockDim.x+threadIdx.x;

    if (globalIdx < global_nb_rna ) {
        int indiv_id = rnas[globalIdx].indiv_id;
        if (dna_size[indiv_id] >= PROM_SIZE) {
        int k = rnas[globalIdx].begin + 22;
        k = k >= dna_size[indiv_id] ? k - dna_size[indiv_id] : k;
        int k_end = k;
        bool found=false;

        do {

            //printf("%d -- %d %ld\n",indiv_id,k,dna_size[indiv_id]);

            if (dna_term[dna_offset[indiv_id]+k] & (1<<7)) {
                int32_t rna_end =
                        k + 10 >= dna_size[indiv_id] ? k + 10 - dna_size[indiv_id] :
                        k +
                        10;

                int32_t rna_length = 0;

                if (rnas[globalIdx].begin > rna_end)
                    rna_length = dna_size[indiv_id] - rnas[globalIdx].begin + rna_end;
                else
                    rna_length = rna_end - rnas[globalIdx].begin;

                if (rna_length < 19) {
                    rnas[globalIdx].begin = 0;
                    rnas[globalIdx].end = 0;
                    rnas[globalIdx].length = 0;
                    rnas[globalIdx].transcribed = false;
                    break;
                }




                rnas[globalIdx].end = rna_end;
                rnas[globalIdx].transcribed = true;
                rnas[globalIdx].length = rna_length;

                if (rnas[globalIdx].end>=dna_size[indiv_id]) {
                    printf("Termin %d %d S %d %ld\n",
                           rnas[globalIdx].begin,
                           rnas[globalIdx].end,indiv_id,dna_size[indiv_id]);
                    //assert(rnas[globalIdx].end<dna_size[indiv_id]);
                }

                found=true;
                break;
            }

            k++;
            k = k >= dna_size[indiv_id] ? k - dna_size[indiv_id] : k;
        } while (k != k_end);
        }
    } else {
        rnas[globalIdx].begin = 0;
        rnas[globalIdx].end = 0;
        rnas[globalIdx].length = 0;
        rnas[globalIdx].transcribed = false;
    }
}

__global__ void display_RNA( pRNA* rna, size_t* dna_size, int32_t global_nb_rna) {
    for(int i = 0; i < global_nb_rna; i++) {
        if (rna[i].transcribed)
            if (rna[i].end>=dna_size[rna[i].indiv_id]) {
                printf("UIIH %d %d S %d -- %ld\n",rna[i].begin,rna[i].end,rna[i].indiv_id,dna_size[rna[i].indiv_id]);
            }
    }

}

__global__
void compute_start_protein(int8_t* start_protein, size_t* dna_size, size_t* dna_offset, pRNA* rna, char* dna, int32_t* nb_proteins,
                           int32_t global_nb_rna, int nb_indiv) {

    const int globalIdx = blockIdx.x*blockDim.x+threadIdx.x;
    int nb_prot = 0;

    if (globalIdx < global_nb_rna) {
        if (rna[globalIdx].transcribed) {
            const int indiv_id = rna[globalIdx].indiv_id;
            int c_pos = rna[globalIdx].begin;
            if (rna[globalIdx].length > 22) {
                c_pos += 22;
                c_pos =
                            c_pos >= dna_size[indiv_id] ? c_pos - dna_size[indiv_id] : c_pos;
                int count_loop=0;

                if (rna[globalIdx].end>=dna_size[indiv_id]) {
                    printf("ator %d S %d\n",rna[globalIdx].end,indiv_id);
                    assert(rna[globalIdx].end<dna_size[indiv_id]);
                }

                while (c_pos != rna[globalIdx].end) {
                    //if (indiv_id==606) printf("%d -- %d %d\n",indiv_id,c_pos,rna[globalIdx].end);
                    bool start = false;
                    int t_pos, k_t;
                    for (int k = 0; k < 9; k++) {
                        k_t = k >= 6 ? k + 4 : k;
                        t_pos = c_pos + k_t >= dna_size[indiv_id] ? c_pos + k_t -
                                                                    dna_size[indiv_id] :
                                c_pos + k_t;
                        count_loop++;
                        if (count_loop>10000) {printf("%d %d %d %d %d %d %d %ld\n",indiv_id,globalIdx,k,
                                                     c_pos,t_pos,
                                                     rna[globalIdx].begin,rna[globalIdx].end,
                                                     dna_size[indiv_id]);assert(0);}
                        if (dna[dna_offset[indiv_id]+t_pos] == SHINE_DAL_SEQ[k]) {
                            start = true;
                        } else {
                            start = false;
                            break;
                        }
                    }


                    start_protein[dna_offset[indiv_id]+c_pos] = start;

                    if (start) {
                        nb_prot++;
                    }

                    c_pos++;
                    c_pos =
                            c_pos >= dna_size[indiv_id] ? c_pos - dna_size[indiv_id] : c_pos;
                }
            }

            atomicAdd(nb_proteins+indiv_id,nb_prot);
            atomicAdd(nb_proteins+nb_indiv,nb_prot);
        }
    }

}


__global__
void compute_protein_offset(int32_t* nb_proteins, int* protein_offset) {

    const int indiv_id = blockIdx.x;
    __shared__ int grid_protein_offset;

    if (threadIdx.x == 0) {
        grid_protein_offset = 0;
    }
    __syncthreads();

    {
        int local_protein_offset = 0;
        for (int cpt = threadIdx.x; cpt < indiv_id; cpt += blockDim.x) {
            local_protein_offset += nb_proteins[cpt];
        }

        if (local_protein_offset > 0)
            atomicAdd(&grid_protein_offset, local_protein_offset);
    }
    __syncthreads();
    if (threadIdx.x == 0) {
        protein_offset[indiv_id] = grid_protein_offset;
    }
}


__global__ void fill_protein(int8_t* start_protein, size_t* dna_offset, int* protein_idx,
                             int* protein_offset, pRNA* rna, pProtein* protein, size_t* dna_size,
                             int32_t global_nb_rna, int nb_indiv) {
	const int globalIdx = blockIdx.x*blockDim.x+threadIdx.x;

    if (globalIdx < global_nb_rna) {
        if (rna[globalIdx].transcribed) {
            int indiv_id = rna[globalIdx].indiv_id;
            int c_pos = rna[globalIdx].begin;
            if (rna[globalIdx].length > 22) {
                c_pos += 22;
                c_pos =
                        c_pos >= dna_size[indiv_id] ? c_pos - dna_size[indiv_id] : c_pos;

                while (c_pos != rna[globalIdx].end) {
                    if (start_protein[dna_offset[indiv_id]+c_pos] == 1) {

                        int local_protein_idx = atomicAdd(protein_idx + indiv_id, 1);
                        atomicAdd(protein_idx + nb_indiv, 1);

                        protein[protein_offset[indiv_id] + local_protein_idx].protein_start = c_pos;
                        protein[protein_offset[indiv_id] + local_protein_idx].indiv_id = rna[globalIdx].indiv_id;
                        protein[protein_offset[indiv_id] + local_protein_idx].stop_RNA = rna[globalIdx].end;
                        protein[protein_offset[indiv_id] + local_protein_idx].translated = false;
                        protein[protein_offset[indiv_id] + local_protein_idx].e = 1.0 -
                                                       fabs(((double) rna[globalIdx].dist)) /
                                                       5.0;
                    }

                    c_pos++;
                    c_pos =
                            c_pos >= dna_size[indiv_id] ? c_pos - dna_size[indiv_id] : c_pos;
                }
            }

        }
    }
}


__global__
void compute_proteins( int8_t* start_protein, size_t* dna_size, size_t* dna_offset, pProtein* protein, char* dna,
                       int32_t global_nb_protein) {
    __shared__ int next_protein_idx;
    if (threadIdx.x == 0) {
        next_protein_idx = 0;
    }
    __syncthreads();

    int local_protein_idx = atomicAdd(&next_protein_idx,1);

    while (local_protein_idx < global_nb_protein) {
        int indiv_id = protein[local_protein_idx].indiv_id;

        int start_protein_pos =  protein[local_protein_idx].protein_start + 13;
        int length = -1;
            start_protein_pos = start_protein_pos >= dna_size[indiv_id] ?
                                start_protein_pos - dna_size[indiv_id]
                                                                        : start_protein_pos;

            if (protein[local_protein_idx].protein_start < protein[local_protein_idx].stop_RNA) {
                length = protein[local_protein_idx].stop_RNA - protein[local_protein_idx].protein_start;
            } else {
                length = dna_size[indiv_id] - protein[local_protein_idx].protein_start + protein[local_protein_idx].stop_RNA + 1;
            }

            length -= 13;


        bool is_protein = false;
        length+=1;
        length = length - (length%3);

        for (int loop_i = 0; length - loop_i >= 2; loop_i+=3) {
            int t_k;

            start_protein_pos = start_protein_pos >= dna_size[indiv_id] ?
                                start_protein_pos - dna_size[indiv_id]
                                                                        : start_protein_pos;
            is_protein = false;

            for (int k = 0; k < 3; k++) {
                t_k = start_protein_pos + k >= dna_size[indiv_id] ?
                      start_protein_pos - dna_size[indiv_id] + k :
                      start_protein_pos + k;

                if (dna[dna_offset[indiv_id]+t_k] == PROTEIN_END[k]) {
                    is_protein = true;
                } else {
                    is_protein = false;
                    break;
                }
            }

            if (is_protein) {
                int prot_length = -1;
                if (protein[local_protein_idx].protein_start + 13 < t_k) {
                    prot_length = t_k - (protein[local_protein_idx].protein_start + 13);
                } else {
                    prot_length = dna_size[indiv_id] - (protein[local_protein_idx].protein_start + 13) + t_k;
                }

                if (prot_length >= 3) {
                    protein[local_protein_idx].protein_end = t_k;
                    protein[local_protein_idx].protein_length = prot_length;
                    protein[local_protein_idx].translated = true;
                }
                break;
            }


            start_protein_pos += 3;
            start_protein_pos = start_protein_pos >= dna_size[indiv_id] ?
                                start_protein_pos - dna_size[indiv_id]
                                                                        : start_protein_pos;

        }



        local_protein_idx = atomicAdd(&next_protein_idx,1);
    }
}


__global__
void translate_proteins( pProtein* protein, size_t* dna_size, char* dna,  size_t* dna_offset, int32_t global_nb_protein, double w_max) {
    __shared__ int next_protein_idx;
    if (threadIdx.x == 0) {
        next_protein_idx = 0;
    }
    __syncthreads();

    int local_protein_idx = atomicAdd(&next_protein_idx,1);

    while (local_protein_idx < global_nb_protein) {
        int indiv_id = protein[local_protein_idx].indiv_id;

        if (protein[local_protein_idx].translated) {

            int c_pos = protein[local_protein_idx].protein_start, t_pos;
            int end_pos = protein[local_protein_idx].protein_end;
            c_pos += 13;
            end_pos -= 3;

            c_pos = c_pos >= dna_size[indiv_id] ? c_pos - dna_size[indiv_id] : c_pos;
            end_pos = end_pos < 0 ? dna_size[indiv_id] + end_pos : end_pos;

            int8_t value = 0;
            int8_t codon_list[64] = {};
            int8_t codon_idx = 0;
            int32_t count_loop = 0;

            bool contin = true;


            while (count_loop < protein[local_protein_idx].protein_length / 3 && codon_idx < 64) {
                value = 0;
                for (int8_t i = 0; i < 3; i++) {
                    t_pos = c_pos + i >= dna_size[indiv_id] ? c_pos + i - dna_size[indiv_id] : c_pos + i;
                    if (dna[dna_offset[indiv_id]+t_pos] == '1') value += 1 << (CODON_SIZE - i - 1);
                }
                codon_list[codon_idx] = value;
                codon_idx++;

                count_loop++;
                c_pos += 3;
                c_pos = c_pos >= dna_size[indiv_id] ? c_pos - dna_size[indiv_id] : c_pos;
            }


            double M = 0.0;
            double W = 0.0;
            double H = 0.0;

            int32_t nb_m = 0;
            int32_t nb_w = 0;
            int32_t nb_h = 0;

            bool bin_m = false; // Initializing to false will yield a conservation of the high weight bit
            bool bin_w = false; // when applying the XOR operator for the Gray to standard conversion
            bool bin_h = false;


            for (int i = 0; i < codon_idx; i++) {
                switch (codon_list[i]) {
                    case CODON_M0 : {
                        // M codon found
                        nb_m++;

                        // Convert Gray code to "standard" binary code
                        bin_m ^= false; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

                        // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
                        //~ M <<= 1;
                        M *= 2;

                        // Add this nucleotide's contribution to M
                        if (bin_m) M += 1;

                        break;
                    }
                    case CODON_M1 : {
                        // M codon found
                        nb_m++;

                        // Convert Gray code to "standard" binary code
                        bin_m ^= true; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

                        // A lower-than-the-previous-lowest bit was found, make a left bitwise shift
                        //~ M <<= 1;
                        M *= 2;

                        // Add this nucleotide's contribution to M
                        if (bin_m) M += 1;

                        break;
                    }
                    case CODON_W0 : {
                        // W codon found
                        nb_w++;

                        // Convert Gray code to "standard" binary code
                        bin_w ^= false; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

                        // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
                        //~ W <<= 1;
                        W *= 2;

                        // Add this nucleotide's contribution to W
                        if (bin_w) W += 1;

                        break;
                    }
                    case CODON_W1 : {
                        // W codon found
                        nb_w++;

                        // Convert Gray code to "standard" binary code
                        bin_w ^= true; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

                        // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
                        //~ W <<= 1;
                        W *= 2;

                        // Add this nucleotide's contribution to W
                        if (bin_w) W += 1;

                        break;
                    }
                    case CODON_H0 :
                    case CODON_START : // Start codon codes for the same amino-acid as H0 codon
                    {
                        // H codon found
                        nb_h++;

                        // Convert Gray code to "standard" binary code
                        bin_h ^= false; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

                        // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
                        //~ H <<= 1;
                        H *= 2;

                        // Add this nucleotide's contribution to H
                        if (bin_h) H += 1;

                        break;
                    }
                    case CODON_H1 : {
                        // H codon found
                        nb_h++;

                        // Convert Gray code to "standard" binary code
                        bin_h ^= true; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

                        // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
                        //~ H <<= 1;
                        H *= 2;

                        // Add this nucleotide's contribution to H
                        if (bin_h) H += 1;

                        break;
                    }
                }
            }



            //  ----------------------------------------------------------------------------------
            //  2) Normalize M, W and H values in [0;1] according to number of codons of each kind
            //  ----------------------------------------------------------------------------------
            protein[local_protein_idx].m = nb_m != 0 ? M / (pow(2, nb_m) - 1) : 0.5;
            protein[local_protein_idx].w = nb_w != 0 ? W / (pow(2, nb_w) - 1) : 0.0;
            protein[local_protein_idx].h = nb_h != 0 ? H / (pow(2, nb_h) - 1) : 0.5;

            //  ------------------------------------------------------------------------------------
            //  3) Normalize M, W and H values according to the allowed ranges (defined in macros.h)
            //  ------------------------------------------------------------------------------------
            // x_min <= M <= x_max
            // w_min <= W <= w_max
            // h_min <= H <= h_max
            protein[local_protein_idx].m = (X_MAX - X_MIN) * protein[local_protein_idx].m + X_MIN;
            protein[local_protein_idx].w = (w_max - W_MIN) * protein[local_protein_idx].w + W_MIN;
            protein[local_protein_idx].h = (H_MAX - H_MIN) * protein[local_protein_idx].h + H_MIN;

            if (nb_m == 0 || nb_w == 0 || nb_h == 0 || protein[local_protein_idx].w == 0.0 ||
                protein[local_protein_idx].h == 0.0) {
                protein[local_protein_idx].is_functional = false;
            } else {
                protein[local_protein_idx].is_functional = true;
            }

        }

        local_protein_idx = atomicAdd(&next_protein_idx,1);
    }
}


__global__
void compute_phenotype( pProtein* protein, int32_t global_nb_protein, double** phenotype,
                        double** phenotype_activ, double** phenotype_inhib, int nb_indiv) {
    __shared__ int next_protein_idx;

    if (threadIdx.x == 0) {
        next_protein_idx = 0;
    }
    __syncthreads();

    int local_protein_idx = atomicAdd(&next_protein_idx,1);

    while (local_protein_idx < global_nb_protein) {

        int indiv_id = protein[local_protein_idx].indiv_id;

        if (protein[local_protein_idx].translated) {
            if (fabs(protein[local_protein_idx].w) < 1e-15 ||
                fabs(protein[local_protein_idx].h) < 1e-15) {

            } else {
                if (protein[local_protein_idx].is_functional) {

                    // Compute triangle points' coordinates
                    double x0 = protein[local_protein_idx].m -
                                protein[local_protein_idx].w;
                    double x1 = protein[local_protein_idx].m;
                    double x2 = protein[local_protein_idx].m +
                                protein[local_protein_idx].w;

                    int ix0 = (int) (x0 * 300);
                    int ix1 = (int) (x1 * 300);
                    int ix2 = (int) (x2 * 300);

                    if (ix0 < 0) ix0 = 0; else if (ix0 > (299)) ix0 = 299;
                    if (ix1 < 0) ix1 = 0; else if (ix1 > (299)) ix1 = 299;
                    if (ix2 < 0) ix2 = 0; else if (ix2 > (299)) ix2 = 299;

                    // Compute the first equation of the triangle
                    double incY = (protein[local_protein_idx].h *
                                   protein[local_protein_idx].e) / (ix1 - ix0);
                    int count = 1;

                    // Updating value between x0 and x1
                    for (int i = ix0 + 1; i < ix1; i++) {
                        if (protein[local_protein_idx].h > 0)
                            atomicAdd(&phenotype_activ[indiv_id][i], (incY * (count++)));
                        else
                            atomicAdd(&phenotype_inhib[indiv_id][i], (incY * (count++)));

                    }


                    if (protein[local_protein_idx].h > 0) {
                        atomicAdd(&phenotype_activ[indiv_id][ix1],
                                  (protein[local_protein_idx].h *
                                   protein[local_protein_idx].e));
                    } else
                        atomicAdd(&phenotype_inhib[indiv_id][ix1],
                                  (protein[local_protein_idx].h *
                                   protein[local_protein_idx].e));


                    // Compute the second equation of the triangle
                    incY =
                            (protein[local_protein_idx].h *
                             protein[local_protein_idx].e) /
                            (ix2 - ix1);
                    count = 1;

                    // Updating value between x1 and x2
                    for (int i = ix1 + 1; i < ix2; i++) {
                        if (protein[local_protein_idx].h > 0)
                            atomicAdd(&phenotype_activ[indiv_id][i],
                                      ((protein[local_protein_idx].h *
                                        protein[local_protein_idx].e) -
                                       (incY * (count++))));
                        else
                            atomicAdd(&phenotype_inhib[indiv_id][i],
                                      ((protein[local_protein_idx].h *
                                        protein[local_protein_idx].e) -
                                       (incY * (count++))));
                    }

                }
            }
        }

        local_protein_idx = atomicAdd(&next_protein_idx,1);
    }

    __syncthreads();

}


__global__ void compute_metaerror_fitness(double selection_pressure,double** phenotype,
                                          double** phenotype_activ,double** phenotype_inhib,
                                          double* target,
                                          double* metaerror, double* fitness) {
    int indiv_id = blockIdx.x;

    int fuzzy_idx = threadIdx.x;

        if (phenotype_activ[indiv_id][fuzzy_idx] > 1.0)
            phenotype_activ[indiv_id][fuzzy_idx] = 1.0;
        if (phenotype_inhib[indiv_id][fuzzy_idx] < -1.0)
            phenotype_inhib[indiv_id][fuzzy_idx] = -1.0;

        phenotype[indiv_id][fuzzy_idx] = phenotype_activ[indiv_id][fuzzy_idx] +
                                                   phenotype_inhib[indiv_id][fuzzy_idx];

    __shared__ double delta[300];

    if (phenotype[indiv_id][fuzzy_idx] > 1) phenotype[indiv_id][fuzzy_idx] = 1;
    if (phenotype[indiv_id][fuzzy_idx] < 0) phenotype[indiv_id][fuzzy_idx] = 0;

    delta[fuzzy_idx] = phenotype[indiv_id][fuzzy_idx] - target[fuzzy_idx];

    __syncthreads();

    if (threadIdx.x == 0) {
        metaerror[indiv_id] = 0;

        for (int i = 0; i < 299; i++) {
            metaerror[indiv_id] +=
                    ((fabs(delta[i]) +
                      fabs(delta[i + 1])) / (600.0));
        }

        fitness[indiv_id] = exp(
                -selection_pressure * ((double)metaerror[indiv_id]));
    }
}


__device__ int32_t Threefry::Device::roulette_random(double* probs, int32_t nb_elts)
{
    double pick_one = 0.0;

    while (pick_one == 0.0)
    {
        pick_one = randomDouble();
    }

    int32_t found_org = 0;

    pick_one -= probs[0];
    while (pick_one > 0)
    {
        assert(found_org<nb_elts-1);

        pick_one -= probs[++found_org];
    }
    return found_org;
}

__global__ void selection(double* fitness, int* next_generation_reproducer, unsigned long long* gpu_counters,
                          int grid_width, int grid_height, int nb_indiv) {
    int indiv_id = blockIdx.x;
    int neightbor = threadIdx.x;

    __shared__ double local_fit_array[NEIGHBORHOOD_SIZE];
    __shared__ double probs[NEIGHBORHOOD_SIZE];
    __shared__ int   count;
    __shared__ double    sum_local_fit;

    int32_t x = indiv_id / grid_height;
    int32_t y = indiv_id % grid_height;

    int cur_x,cur_y;

    if (threadIdx.x == 0) {
        count             = 0;
        sum_local_fit     = 0.0;
    }

    __syncthreads();

    if (threadIdx.x == 0) {

        for (int8_t i = -1; i < SELECTION_SCOPE_X - 1; i++) {
            for (int8_t j = -1; j < SELECTION_SCOPE_Y - 1; j++) {
                cur_x = (x + i + grid_width) % grid_width;
                cur_y = (y + j + grid_height) % grid_height;

                local_fit_array[count] = fitness[cur_x * grid_height + cur_y];
                atomicAdd(&sum_local_fit, local_fit_array[count]);

                count++;
            }
        }
    }

    __syncthreads();


    //for(int16_t i = 0 ; i < NEIGHBORHOOD_SIZE ; i++) {

        probs[neightbor] = local_fit_array[neightbor]/sum_local_fit;

    __syncthreads();

    if (threadIdx.x == 0) {
        Threefry::Device rng(gpu_counters,indiv_id,Threefry::Phase::REPROD,nb_indiv);
        int found_org = rng.roulette_random(probs, NEIGHBORHOOD_SIZE);

        int x_offset = (found_org / SELECTION_SCOPE_X) - 1;
        int y_offset = (found_org % SELECTION_SCOPE_Y) - 1;

        next_generation_reproducer[indiv_id] = ((x + x_offset + grid_width) % grid_width) * grid_height +
                                               ((y + y_offset + grid_height) % grid_height);
    }
}


__constant__ double cof[6] = {  76.18009172947146,
                                -86.50532032941677,
                                24.01409824083091,
                                -1.231739572450155,
                                0.1208650973866179e-2,
                                -0.5395239384953e-5 };



// Returns the value ln[gamma(X)] for X.
// The gamma function is defined by the integral  gamma(z) = int(0, +inf, t^(z-1).e^(-t)dt).
// When the argument z is an integer, the gamma function is just the familiar factorial
// function, but offset by one, n! = gamma(n + 1).
__device__ static double gammln(double X)
{
    double x, y, tmp, ser;

    y = x = X;
    tmp = x + 5.5;
    tmp -= (x+0.5) * log(tmp);
    ser = 1.000000000190015;

    for (int8_t j = 0 ; j <= 5 ; j++)
    {
        ser += cof[j] / ++y;
    }

    return -tmp + log(2.5066282746310005 * ser / x);
}


__device__ 
int32_t Threefry::Device::binomial_random(int32_t nb_drawings, double prob)
{
    int32_t nb_success;

    // The binomial distribution is invariant under changing
    // ProbSuccess to 1-ProbSuccess, if we also change the answer to
    // NbTrials minus itself; we ll remember to do this below.
    double p;
    if (prob <= 0.5) p = prob;
    else p = 1.0 - prob;

    // mean of the deviate to be produced
    double mean = nb_drawings * p;


    if (nb_drawings < 25)
        // Use the direct method while NbTrials is not too large.
        // This can require up to 25 calls to the uniform random.
    {
        nb_success = 0;
        for (int32_t j = 1 ; j <= nb_drawings ; j++)
        {
            if (randomDouble() < p) nb_success++;
        }
    }
    else if (mean < 1.0)
        // If fewer than one event is expected out of 25 or more trials,
        // then the distribution is quite accurately Poisson. Use direct Poisson method.
    {
        double g = exp(-mean);
        double t = 1.0;
        int32_t j;
        for (j = 0; j <= nb_drawings ; j++)
        {
            t = t * randomDouble();
            if (t < g) break;
        }

        if (j <= nb_drawings) nb_success = j;
        else nb_success = nb_drawings;
    }

    else
        // Use the rejection method.
    {
        double en     = nb_drawings;
        double oldg   = gammln(en + 1.0);
        double pc     = 1.0 - p;
        double plog   = log(p);
        double pclog  = log(pc);

        // rejection method with a Lorentzian comparison function.
        double sq = sqrt(2.0 * mean * pc);
        double angle, y, em, t;
        do
        {
            do
            {
                angle = M_PI * randomDouble();
                y = tan(angle);
                em = sq*y + mean;
            } while (em < 0.0 || em >= (en + 1.0)); // Reject.

            em = floor(em); // Trick for integer-valued distribution.
            t = 1.2 * sq * (1.0 + y*y)
                * exp(oldg - gammln(em + 1.0) - gammln(en - em + 1.0) + em * plog + (en - em) * pclog);

        } while (randomDouble() > t); // Reject. This happens about 1.5 times per deviate, on average.

        nb_success = (int32_t) rint(em);
    }


    // Undo the symmetry transformation.
    if (p != prob) nb_success = nb_drawings - nb_success;

    return nb_success;
}


__global__
void generate_mutations(unsigned long long* gpu_counters, size_t* dna_size, int* nb_mutations,
                        GPUDnaMutator* dna_mutator_list,int* next_generation_reproducer,
                        int nb_indivs, double mutation_rate) {
    int indiv_id = blockIdx.x;

    Threefry::Device rng(gpu_counters,indiv_id,Threefry::Phase::MUTATION,nb_indivs);

    double mutation_r = mutation_rate;
    int prev_gen_id = next_generation_reproducer[indiv_id];
    size_t prev_gen_size = dna_size[prev_gen_id];

    // Small mutations
    dna_mutator_list[indiv_id].nb_swi_ = rng.
            binomial_random(prev_gen_size, mutation_r);
    dna_mutator_list[indiv_id].nb_mut_ = dna_mutator_list[indiv_id].nb_swi_;
    dna_mutator_list[indiv_id].cpt_mut_ = dna_mutator_list[indiv_id].nb_mut_;

    nb_mutations[indiv_id] = dna_mutator_list[indiv_id].nb_mut_;
    atomicAdd(nb_mutations+nb_indivs,nb_mutations[indiv_id]);
}

__global__
void compute_tab_mutations_offset(int* nb_mutations, int* mutations_offset) {
    const int indiv_id = blockIdx.x;
    __shared__ int grid_mutation_offset;

    if (threadIdx.x == 0) {
        grid_mutation_offset = 0;
    }
    __syncthreads();

    {
        int local_mutation_offset = 0;
        for (int cpt = threadIdx.x; cpt < indiv_id; cpt += blockDim.x) {
            local_mutation_offset += nb_mutations[cpt];
        }

        if (local_mutation_offset > 0)
            atomicAdd(&grid_mutation_offset, local_mutation_offset);
    }
    __syncthreads();
    if (threadIdx.x == 0) {
        mutations_offset[indiv_id] = grid_mutation_offset;
    }
}

__device__ static int mod(int a, int b)
{

    assert(b > 0);

    while (a < 0)  a += b;
    while (a >= b) a -= b;

    return a;
}

__global__
void predict_size_v2(size_t* dna_size, size_t* next_gen_dna_size, GPUDnaMutator* dna_mutator_list,
                     TypeMutation* tab_mut,
                     int* nb_mutations, int* mutations_offset,
                     unsigned long long* gpu_counters,int* next_generation_reproducer,
                     int max_genome_length,
                     int min_genome_length, int nb_indiv) {
    const int indiv_id = blockIdx.x;

    int random_value;

    int transient_size = dna_size[next_generation_reproducer[indiv_id]];

    Threefry::Device rng(gpu_counters,indiv_id,Threefry::Phase::MUTATION,nb_indiv);

    for (int mut_idx = 0; mut_idx <  dna_mutator_list[indiv_id].nb_mut_; mut_idx++) {
            dna_mutator_list[indiv_id].cpt_mut_--;

                dna_mutator_list[indiv_id].nb_swi_--;

                int pos = rng.random(transient_size);

                tab_mut[mutations_offset[indiv_id]+mut_idx].type_ = MutationEventType::DO_SWITCH;
                tab_mut[mutations_offset[indiv_id]+mut_idx].pos_1_ = pos;

    }

    next_gen_dna_size[indiv_id] = transient_size;
}


__global__ void display_mut(TypeMutation* tab_mut,
                            int* nb_mutations, int* mutations_offset) {
    for (int indiv_id = 0; indiv_id < 25; indiv_id++) {
        printf("nb mut %d : %d %d\n",indiv_id,nb_mutations[indiv_id],mutations_offset[indiv_id]);
        for (int i = 0; i < nb_mutations[indiv_id]; i++) {
            printf("%d -- %d %d %d %d %d %c%c%c%c%c%c %d\n", i, tab_mut[mutations_offset[indiv_id] + i].type_,
                   tab_mut[mutations_offset[indiv_id] + i].pos_1_,
                   tab_mut[mutations_offset[indiv_id] + i].pos_2_,
                   tab_mut[mutations_offset[indiv_id] + i].pos_3_,
                   tab_mut[mutations_offset[indiv_id] + i].number_,
                   tab_mut[mutations_offset[indiv_id] + i].seq[0],
                   tab_mut[mutations_offset[indiv_id] + i].seq[1],
                   tab_mut[mutations_offset[indiv_id] + i].seq[2],
                   tab_mut[mutations_offset[indiv_id] + i].seq[3],
                   tab_mut[mutations_offset[indiv_id] + i].seq[4],
                   tab_mut[mutations_offset[indiv_id] + i].seq[5],
                   tab_mut[mutations_offset[indiv_id] + i].transient_size);
        }
    }

}


__global__
void compute_next_gen_dna_offset(size_t* next_gen_dna_size, size_t* next_gen_dna_offset) {

    const int indiv_id = blockIdx.x;
    __shared__ int grid_dna_offset;

    if (threadIdx.x == 0) {
        grid_dna_offset = 0;
    }
    __syncthreads();

    {
        int local_dna_offset = 0;
        for (int cpt = threadIdx.x; cpt < indiv_id; cpt += blockDim.x) {
            local_dna_offset += next_gen_dna_size[cpt];
        }

        if (local_dna_offset > 0)
            atomicAdd(&grid_dna_offset, local_dna_offset);
    }
    __syncthreads();
    if (threadIdx.x == 0) {
        next_gen_dna_offset[indiv_id] = grid_dna_offset;
    }
}

__global__ void do_mutation_v2(TypeMutation* tab_mut,
                               int* nb_mutations, size_t* dna_size, size_t* dna_offset, char* dna,
                               char* next_gen_dna, size_t* next_gen_dna_size, size_t* next_gen_dna_offset,
                               int* next_generation_reproducer,  int* mutations_offset, unsigned long long int* nb_mut_bp) {

    int dna_pos_block = blockIdx.x;
    int indiv_id = blockIdx.y;
    int32_t locus = (dna_pos_block*128)+threadIdx.x;
    int32_t next_locus = locus;

    if (locus < next_gen_dna_size[indiv_id]) {

        int8_t mutate = 0;
        int nb_events = nb_mutations[indiv_id];


        for (; nb_events > 0; nb_events--) {
            auto &mut = tab_mut[mutations_offset[indiv_id]+nb_events - 1];

            switch (mut.type_) {
                case DO_SWITCH:
                    if (locus == mut.pos_1_)
                        mutate = not mutate;
                    break;
            }
        }

        assert(locus >= 0);

        assert(locus < dna_size[next_generation_reproducer[indiv_id]]);

        auto base = dna[dna_offset[next_generation_reproducer[indiv_id]]+locus];
        if (mutate) base = (base == '0') ? '1' : '0';

        next_gen_dna[next_gen_dna_offset[indiv_id]+next_locus] = base;

    }
}


void run_a_step_on_GPU(int nb_indiv, double w_max, double selection_pressure, int grid_width, int grid_height, double mutation_rate) {
    int x_dim_size = (host_max_dna_size / 128)+1;

    int y_dim_size = nb_indiv;

    dim3 dimGrid(x_dim_size,y_dim_size);

    search_start_stop_RNA<<<dimGrid,128>>>(dna_size,dna,dna_offset,
            nb_promoters,dna_term,nb_indiv,global_dna_size,nb_mut_bp);

    int total_nb_promoters_host;
    checkCuda(cudaMemcpy(&total_nb_promoters_host,
                         nb_promoters+nb_indiv, sizeof(int), cudaMemcpyDeviceToHost));

    if (total_nb_promoters_host > current_size_rna_list) {
        checkCuda(cudaFree(rna));
        current_size_rna_list = total_nb_promoters_host * 1.1;
        checkCuda(cudaMalloc(&rna,current_size_rna_list* sizeof(pRNA)));
    }

    compute_RNA_offset<<<nb_indiv,128>>>(nb_promoters,rna_offset);

    fill_RNA<<<dimGrid,128>>>( dna_term, dna_size,dna_offset, nb_promoters, rna_offset, rna, rna_idx,nb_indiv);

    int global_nb_rna;
    checkCuda(cudaMemcpy(&global_nb_rna,
                         rna_idx+nb_indiv, sizeof(int), cudaMemcpyDeviceToHost));


    compute_RNA<<<global_nb_rna/128+1,128>>>( dna_term,dna_size, dna_offset, rna, global_nb_rna);

    cudaDeviceSynchronize();
    compute_start_protein<<<global_nb_rna,1>>>(start_protein, dna_size, dna_offset, rna, dna, nb_proteins,
            global_nb_rna, nb_indiv);
    cudaDeviceSynchronize();

        int total_nb_protein_host;
    checkCuda(cudaMemcpy(&total_nb_protein_host,
                         nb_proteins+nb_indiv, sizeof(int), cudaMemcpyDeviceToHost));

    if (total_nb_protein_host > current_size_protein_list) {
        checkCuda(cudaFree(protein));
        current_size_protein_list = total_nb_protein_host * 1.1;
        checkCuda(cudaMalloc(&protein,current_size_protein_list* sizeof(pProtein)));
    }

    compute_protein_offset<<<nb_indiv,128>>>(nb_proteins, protein_offset);

    fill_protein<<<global_nb_rna/128+1,128>>>(start_protein,dna_offset, protein_idx, protein_offset, rna, protein,
            dna_size, global_nb_rna, nb_indiv);

    int global_nb_protein;
    checkCuda(cudaMemcpy(&global_nb_protein,
                         protein_idx+nb_indiv, sizeof(int), cudaMemcpyDeviceToHost));

    compute_proteins<<<1,128>>>( start_protein, dna_size, dna_offset,protein, dna, global_nb_protein);

    translate_proteins<<<1,128>>>( protein, dna_size, dna, dna_offset, global_nb_protein, w_max);

    compute_phenotype<<<1,128>>>( protein,global_nb_protein, phenotype,
            phenotype_activ,phenotype_inhib, nb_indiv);

    compute_metaerror_fitness<<<nb_indiv,300>>>(selection_pressure,phenotype,
                                    phenotype_activ,phenotype_inhib,
                                   target,
                                   metaerror, fitness);

    // SELECTION
    selection<<<nb_indiv,NEIGHBORHOOD_SIZE>>>(fitness,next_generation_reproducer,gpu_counters,
            grid_width,grid_height,nb_indiv);

    // GENERATE MUTATION + PREDICT
    generate_mutations<<<nb_indiv,1>>>(gpu_counters,dna_size,nb_mutations,dna_mutator_list,
            next_generation_reproducer,
            nb_indiv,mutation_rate);


    compute_tab_mutations_offset<<<nb_indiv,1>>>(nb_mutations,mutations_offset);

    int total_nb_mutations_host;
    checkCuda(cudaMemcpy(&total_nb_mutations_host,
                         nb_mutations+nb_indiv, sizeof(int), cudaMemcpyDeviceToHost));

    if (total_nb_mutations_host > current_size_tab_mutation) {
        checkCuda(cudaFree(tab_mutation));
        current_size_tab_mutation = total_nb_mutations_host * 1.1;
        checkCuda(cudaMalloc(&tab_mutation,current_size_tab_mutation* sizeof(TypeMutation)));
    }

    int min_genome_length_  = 10;
    int max_genome_length_  = 10000000;

    predict_size_v2<<<nb_indiv,1>>>(dna_size, next_gen_dna_size, dna_mutator_list,
            tab_mutation,nb_mutations,mutations_offset,gpu_counters,next_generation_reproducer,
            max_genome_length_,min_genome_length_,nb_indiv);
    cudaDeviceSynchronize();
    // DO MUTATION

    std::vector <size_t> host_dna_size(
            nb_indiv);

    checkCuda(cudaMemcpy(host_dna_size.data(),
                         next_gen_dna_size, nb_indiv * sizeof(size_t), cudaMemcpyDeviceToHost));

    global_dna_size=0;
    for (int i = 0; i < nb_indiv; i++) {
        global_dna_size += host_dna_size[i];
        host_max_dna_size = host_max_dna_size < host_dna_size[i] ?
                            host_dna_size[i] : host_max_dna_size;
    }

    bool haveChange = false;
    if (global_dna_size >= allocated_global_dna_size) {
        haveChange = true;
        allocated_global_dna_size = global_dna_size*2;

        checkCuda(cudaMalloc((void **) &next_gen_dna, allocated_global_dna_size * sizeof(char)));
        checkCuda(cudaFree(dna_term));
        checkCuda(cudaMalloc((void **) &dna_term, allocated_global_dna_size * sizeof(int8_t * )));

        checkCuda(cudaFree(start_protein));
        checkCuda(cudaMalloc((void **) &start_protein, allocated_global_dna_size * sizeof(int8_t * )));
    }



    compute_next_gen_dna_offset<<<nb_indiv,128>>>(next_gen_dna_size, next_gen_dna_offset);

    x_dim_size = (host_max_dna_size / 128)+1;
    y_dim_size = nb_indiv;

    dim3 dimGrid2(x_dim_size,y_dim_size);


    do_mutation_v2<<<dimGrid2,128>>>(tab_mutation,
                                   nb_mutations, dna_size, dna_offset, dna,
                                   next_gen_dna, next_gen_dna_size, next_gen_dna_offset,next_generation_reproducer,
            mutations_offset,nb_mut_bp);

    //printf("DNA 1 %p\n",dna);
    //next_generation_dna_read<<<1,1>>>(next_gen_dna, next_gen_dna_offset,next_gen_dna_size, global_dna_size);

    // SWITCH STRUCTURE

    int block = ceil(nb_indiv/32);
    do_memset<<<block,32>>>(phenotype_activ,phenotype_inhib,nb_mutations,rna_idx,protein_idx,nb_proteins,
            nb_promoters,next_gen_dna_size,
            nb_indiv);

    //allocate_next_gen(nb_indiv);
    //printf("DNA 2 %p\n",dna);

    size_t* tmp_dna_size = dna_size;
    dna_size = next_gen_dna_size;
    next_gen_dna_size = tmp_dna_size;


    size_t* tmp_dna_offset = dna_offset;
    dna_offset = next_gen_dna_offset;
    next_gen_dna_offset = tmp_dna_offset;

    //global_dna_size = global_next_gen_dna_size;
    cudaDeviceSynchronize();

    assert(dna!=0);
    //printf("DNA 3 %p\n",dna);

    if (haveChange) {
        checkCuda(cudaFree(dna));
        checkCuda(cudaMalloc((void **) &dna, allocated_global_dna_size * sizeof(char)));
    }


    //printf("DNA 4 %p\n",dna);

    cudaDeviceSynchronize();


    char* dna_tmp = dna;
    dna = next_gen_dna;
    next_gen_dna = dna_tmp;

    //  clean(exp_m);
}

void allocate_next_gen(int nb_indiv) {
    for (int indiv_id = 0; indiv_id < nb_indiv; indiv_id++) {
        checkCuda(cudaMemset(host_phenotype[indiv_id], 0.0, 300 * sizeof(double)));
        checkCuda(cudaMemset(host_phenotype_activ[indiv_id], 0.0, 300 * sizeof(double)));
        checkCuda(cudaMemset(host_phenotype_inhib[indiv_id], 0.0, 300 * sizeof(double)));
    }

    checkCuda(cudaMemset(nb_mutations, 0, (nb_indiv+1) * sizeof(int)));
    checkCuda(cudaMemset(mutations_offset, 0, nb_indiv * sizeof(int)));
    checkCuda(cudaMemset(mutations_idx, 0, nb_indiv * sizeof(int)));

    checkCuda(cudaMemset(rna_idx, 0, (nb_indiv+1) * sizeof(int32_t)));
    checkCuda(cudaMemset(rna_offset, 0, nb_indiv * sizeof(int32_t)));

    checkCuda(cudaMemset(protein_idx, 0, (nb_indiv+1) * sizeof(int32_t)));
    checkCuda(cudaMemset(protein_offset, 0, nb_indiv * sizeof(int32_t)));
    checkCuda(cudaMemset(nb_proteins, 0, (nb_indiv+1) * sizeof(int)));

    checkCuda(cudaMemset(nb_promoters, 0, (nb_indiv+1) * sizeof(int)));
}

__global__
void do_memset(double** phenotype_activ, double** phenotype_inhib, int* nb_mutations, int32_t* rna_idx,
               int32_t* protein_idx, int* nb_proteins, int* nb_promoters,
               size_t* dna_size,
               int nb_indiv) {
    const int indiv_id = blockIdx.x * blockDim.x + threadIdx.x;

    if (indiv_id < nb_indiv) {
        for (int i = 0; i < 300; i++) {
            phenotype_inhib[indiv_id][i] = 0;
            phenotype_activ[indiv_id][i] = 0;
        }

        rna_idx[indiv_id] = 0;
        protein_idx[indiv_id] = 0;
        nb_proteins[indiv_id] = 0;
        nb_promoters[indiv_id] = 0;

        if (indiv_id == 0) {
            nb_mutations[nb_indiv] = 0;

            rna_idx[nb_indiv] = 0;
            protein_idx[nb_indiv] = 0;
            nb_proteins[nb_indiv] = 0;
            nb_promoters[nb_indiv] = 0;
        }
    }
}
