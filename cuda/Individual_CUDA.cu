//
// Created by elturpin on 22/10/2020.
//

#include "Individual_CUDA.cuh"
#include "misc_functions.cuh"

#include <cstdio>
#include <cuda.h>
#include <cuda_runtime_api.h>

__device__ void Individual_CUDA::evaluate() {
    uint idx = threadIdx.x;
    uint rr_width = blockDim.x;
    for (uint position = idx; position < size; position += rr_width) {
        const char *genome_at_pos_leading = genome + position;
        const char *genome_at_pos_lagging = genome + (size - 1 - position);

        promoters[position]   = is_promoter(genome_at_pos_leading);
        terminators[position] = is_terminator(genome_at_pos_leading);
        terminators[position] = is_terminator(genome_at_pos_lagging);
        prot_start[position]  = is_prot_start(genome_at_pos_leading);
        prot_start[position]  = is_prot_start(genome_at_pos_lagging);
    }
    __syncthreads();

    if (idx == 0) {
        prepare_rnas();
    }
    if (idx == 1) {
        nb_terminator = sparse(size, terminators);
    }
    if (idx == 2) {
        nb_prot_start = sparse(size, prot_start);
    }
    __syncthreads();
    if (nb_prot_start <= 0 || nb_terminator <= 0)
        return;

    for (uint rna_idx = idx; rna_idx < nb_rnas; rna_idx += rr_width) {
        compute_rna(rna_idx);
        prepare_gene(rna_idx);
    }

    __syncthreads();
    if (idx == 0) {
        gather_genes();
    }

    __syncthreads();
    for (uint gene_idx = idx; gene_idx < nb_gene; gene_idx += rr_width) {
        translate_gene(gene_idx);
    }
}

__device__ void Individual_CUDA::prepare_rnas() {
    int insert_position = 0;

    for (uint read_position = 0; read_position < size; ++read_position) {
        uint8_t read_value = promoters[read_position];
        if (read_value <= PROM_MAX_DIFF) {
            auto &rna = list_rnas[insert_position];
            rna.errors = read_value;
            rna.start_transcription = read_position + PROM_SIZE;
            if (rna.start_transcription >= size)
                rna.start_transcription -= size;
            insert_position++;
        }
    }

    nb_rnas = insert_position;
}

__device__ void Individual_CUDA::compute_rna(uint rna_idx) const {
    auto &rna = list_rnas[rna_idx];
    uint start_transcript = rna.start_transcription;
    if (not nb_terminator) {
        rna.errors = 0b1111u;
        return;
    }
    // get end of transcription
    // find the smallest element greater than start
    uint idx_term = find_smallest_greater(start_transcript, terminators, nb_terminator);
    uint term_position = terminators[idx_term];
    uint transcript_length = get_distance(start_transcript, term_position) + TERM_SIZE;

    if (transcript_length < DO_TRANSLATION_LOOP) {
        rna.errors = 0b1111u;
    } else {
        rna.transcription_length = transcript_length;
    }
}

__device__ void Individual_CUDA::prepare_gene(uint rna_idx) const {
    auto &rna = list_rnas[rna_idx];
    if (rna.errors > 5) {
        rna.nb_gene = 0;
        rna.list_gene = nullptr;
        return;
    }

    uint nb_ps = nb_prot_start;
    if (not nb_ps) {
        return;
    }
    uint *list_ps = prot_start;
    uint nb_gene = 0;

    // Correctly setup the research
    uint max_distance = rna.transcription_length - DO_TRANSLATION_LOOP;
    uint first_next_ps = find_smallest_greater(rna.start_transcription, list_ps, nb_ps);
    uint ps_idx = first_next_ps;

    uint distance = get_distance(rna.start_transcription, list_ps[ps_idx]);
    while (distance <= max_distance) {
        nb_gene++;
        uint prev_idx = ps_idx++;
        if (ps_idx == nb_ps)
            ps_idx = 0;
        distance += get_distance_ori(list_ps[prev_idx], list_ps[ps_idx]);
    }
    // all potential genes are counted
    // Let us put their position in a list

    rna.nb_gene = nb_gene;
    rna.list_gene = nb_gene ? new Gene[nb_gene] : nullptr;
    for (int i = 0; i < nb_gene; ++i) {
        uint start = list_ps[first_next_ps] + SD_TO_START;
        if (start >= size) {
            start -= size;
        }

        rna.list_gene[i].start = start;
        if (++first_next_ps >= nb_ps) {
            first_next_ps = 0;
        }
    }
}

__device__ void Individual_CUDA::gather_genes() {
    nb_gene = 0;
    for (int idx_rna = 0; idx_rna < nb_rnas; ++idx_rna) {
        nb_gene += list_rnas[idx_rna].nb_gene;
    }

    list_gene = new Gene[nb_gene];
    list_protein = new Protein[nb_gene];
    uint insert_idx = 0;

    for (int idx_rna = 0; idx_rna < nb_rnas; ++idx_rna) {
        auto &rna = list_rnas[idx_rna];
        for (int i = 0; i < rna.nb_gene; ++i) {
            list_gene[insert_idx] = rna.list_gene[i];
            list_gene[insert_idx].concentration = PROM_MAX_DIFF + 1 - rna.errors;
            list_gene[insert_idx].length_limit = rna.transcription_length;
            insert_idx++;
        }
    }
}

__device__ void Individual_CUDA::translate_gene(uint gene_idx) const {
    const auto &gene = list_gene[gene_idx];

    uint max_distance = gene.length_limit;

    auto next = [this](uint &it) -> void {
        it += 3;
        if (it >= this->size)
            it -= this->size;
    };

    uint it = gene.start;
    uint distance = 0;
    auto &new_protein = list_protein[gene_idx];

    while (true) {
        uint8_t codon = translate_to_codon(genome + it);
        if (codon == CODON_STOP)
            break;
        distance += CODON_SIZE;
        if (distance > max_distance) {
            distance = false;
            break;
        }

        new_protein.add_codon(codon);
        next(it);
    }

    if (distance) {
        // Gene has been translated into a protein
        new_protein.concentration = (float) gene.concentration / (float) (PROM_MAX_DIFF + 1);
        new_protein.normalize();
    } else
        new_protein.concentration = 0.0;
}

// ************ PRINTING
#define PRINT_HEADER(STRING)  printf("\n<%s>\n", STRING)
#define PRINT_CLOSING(STRING) printf("</%s>\n", STRING)

__device__ void Individual_CUDA::print_metadata_summary() const {
    __syncthreads();
    if (threadIdx.x == 0) {
        PRINT_HEADER("METADATA_SUMMARY");
        printf("proms: %d, terms: %d, prots: %d\n", nb_rnas, nb_terminator, nb_prot_start);
        PRINT_CLOSING("METADATA_SUMMARY");
    }
    __syncthreads();
}

__device__ void Individual_CUDA::print_rnas() const {
    __syncthreads();
    if (threadIdx.x == 0) {
        PRINT_HEADER("ARNS");
        uint nb_coding_rna = 0;
        for (int i = 0; i < nb_rnas; ++i) {
            const auto &rna = list_rnas[i];
            if (rna.errors <= PROM_MAX_DIFF) { // is LEADING
                nb_coding_rna++;
                uint start = rna.start_transcription;
                uint end = (start + rna.transcription_length) % size;
                printf("%d -> %d | %d\n", start, end, rna.errors);
            }
        }

        printf("\nnumber of terminated rna: %u\n", nb_coding_rna);
        PRINT_CLOSING("ARNS");
    }
    __syncthreads();
}

__device__ void Individual_CUDA::print_gathered_genes() const {
    __syncthreads();
    if (threadIdx.x == 0) {
        PRINT_HEADER("GENES");
        uint nb_gene = 0;
        for (int i = 0; i < nb_gene; ++i) {// is LEADING
            auto start = list_gene[i].start;
            nb_gene++;
            printf("\t%d: concentration: %d\n", start, list_gene[i].concentration);
        }

        printf("\nnumber of potential gene: %u <> %u|%u\n", nb_gene);
        PRINT_CLOSING("GENES");
    }
    __syncthreads();
}

__device__ void Individual_CUDA::print_proteins() const {
    __syncthreads();
    if (threadIdx.x == 0) {
        PRINT_HEADER("PROTEINS");
        uint nb_prot = 0;
        for (int i = 0; i < nb_gene; ++i) {
            const auto &prot = list_protein[i];
            nb_prot++;
            printf("%d: %d %f %f %f %f\n",
                   list_gene[i].start, prot.is_functional(),
                   prot.concentration, prot.width, prot.mean, prot.height);

        }

        printf("\nnumber of proteins: %u <> %u|%u\n", nb_prot);
        PRINT_CLOSING("PROTEINS");
    }
    __syncthreads();
}
