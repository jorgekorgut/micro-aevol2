//
// Created by elturpin on 03/12/2020.
//

#include <vector>
#include "PetriDish.cuh"
#include "Algorithms.cuh"

inline
void checkCuda(cudaError_t result)
{
#if !defined(_NDEBUG)
    if (result != cudaSuccess) {
        fprintf(stderr, "CUDA Runtime Error: %s\n",
                cudaGetErrorString(result));
        assert(result == cudaSuccess);
    }
#endif
}

void PetriDish::run_a_step_on_GPU(int nb_indiv, double w_max, double selection_pressure, int grid_width, int grid_height, double mutation_rate) {
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
