// ***************************************************************************************************************
//
//          Mini-Aevol is a reduced version of Aevol -- An in silico experimental evolution platform
//
// ***************************************************************************************************************
//
// Copyright: See the AUTHORS file provided with the package or <https://gitlab.inria.fr/rouzaudc/mini-aevol>
// Web: https://gitlab.inria.fr/rouzaudc/mini-aevol
// E-mail: See <jonathan.rouzaud-cornabas@inria.fr>
// Original Authors : Jonathan Rouzaud-Cornabas
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ***************************************************************************************************************



#include <cmath>
#include "Organism.h"

using namespace std;

/**
 * Constructor to generate a random organism (i.e. an organism with a random DNA)
 *
 * @param length : Length of the generated random DNA
 */
Organism::Organism(int length, std::shared_ptr<std::mt19937_64> rng) {
    rna_count_ = 0;

    dna_ = new Dna(length, std::move(rng));
}

/**
 * Constructor to create a clone of a given Organism
 *
 * @param clone : The organism to clone
 */
Organism::Organism(const std::shared_ptr<Organism> &clone) {
    rna_count_ = 0;
    dna_ = new Dna(*(clone->dna_));
    promoters_ = clone->promoters_;
}

/**
 * Create an Organism from a backup/checkpointing file
 *
 * @param backup_file : gzFile to read from
 */
Organism::Organism(gzFile backup_file) {
    rna_count_ = 0;

    load(backup_file);
}

/**
 * Destructor of an organism
 */
Organism::~Organism() {
    for (auto rna : rnas) {
        delete (rna);
    }
    rnas.clear();

    for (auto prot : proteins) {
        delete (prot);
    }
    proteins.clear();

    terminators.clear();

    delete dna_;
}

/**
 * Save the organism to backup/checkpointing file
 *
 * @param backup_file : where to the save the organism
 */
void Organism::save(gzFile backup_file) const {
    dna_->save(backup_file);
}

/**
 * Load the organism from backup/checkpointing file
 *
 * @param backup_file : from where restore the organism
 */
void Organism::load(gzFile backup_file) {
    dna_ = new Dna();
    dna_->load(backup_file);
}

/**
 * Reset the stats variable (used when an organism is a perfect clone of its parent, it means no mutation)
 */
void Organism::reset_mutation_stats() {
    nb_swi_ = 0;
    nb_mut_ = 0;
}

void Organism::compute_protein_stats() {
    nb_genes_activ = 0;
    nb_genes_inhib = 0;
    nb_func_genes = 0;
    nb_non_func_genes = 0;
    nb_coding_RNAs = 0;
    nb_non_coding_RNAs = 0;

    for (int i = 0; i < rna_count_; i++) {
        if (rnas[i] != nullptr) {
            if (rnas[i]->is_coding_)
                nb_coding_RNAs++;
            else
                nb_non_coding_RNAs++;
        }
    }

    for (int i = 0; i < protein_count_; i++) {
        if (rnas[i] != nullptr) {
            if (proteins[i]->is_functional) {
                nb_func_genes++;
            } else {
                nb_non_func_genes++;
            }
            if (proteins[i]->h > 0) {
                nb_genes_activ++;
            } else {
                nb_genes_inhib++;
            }
        }
    }
}

/** LOOK **/
void Organism::locate_promoters() {
    look_for_new_promoters_starting_between(0, length());
}

/**
 * Apply all the mutation events of the organism on its DNA
 */
void Organism::apply_mutations(const list<MutationEvent *> &mutation_list) {
    for (const auto *mutation: mutation_list) {
        switch (mutation->type()) {
            case DO_SWITCH:
                do_switch(mutation->pos_1());
                nb_swi_++;
                nb_mut_++;
                break;
        }
    }
}

void Organism::evaluate(const double *target) {
    compute_RNA();
    search_start_protein();
    compute_protein();
    translate_protein();
    compute_phenotype();
    compute_fitness(target);
}

void Organism::compute_RNA() {
    proteins.clear();
    rnas.clear();
    terminators.clear();

    rnas.resize(promoters_.size());

    for (const auto &prom_pair: promoters_) {
        int prom_pos = prom_pair.first;

        /* Search for terminators */
        int cur_pos = prom_pos + PROM_SIZE;
        loop_back(cur_pos);

        int start_pos = cur_pos;

        bool terminator_found = false;

        while (!terminator_found) {
            int term_dist_leading = dna_->terminator_at(cur_pos);

            if (term_dist_leading == TERM_STEM_SIZE)
                terminator_found = true;
            else {
                cur_pos++;
                loop_back(cur_pos);

                if (cur_pos == start_pos) {
                    break;
                }
            }
        }

        if (terminator_found) {
            int32_t rna_end = cur_pos + TERM_SIZE;
            loop_back(rna_end);

            int32_t rna_length;

            if (start_pos > rna_end)
                rna_length = (length() + rna_end) - start_pos;
            else
                rna_length = rna_end - start_pos;

            if (rna_length > 0) {
                int glob_rna_idx = rna_count_;
                rna_count_ = rna_count_ + 1;

                rnas[glob_rna_idx] = new RNA(
                        prom_pos,
                        rna_end,
                        1.0 - std::fabs(((float) prom_pair.second)) / 5.0,
                        rna_length);
            }
        }
    }
}

void Organism::search_start_protein() {
    for (int rna_idx = 0; rna_idx < rna_count_; rna_idx++) {
        const auto &rna = rnas[rna_idx];
        int c_pos = rna->begin;

        if (rna->length >= PROM_SIZE) {
            c_pos += PROM_SIZE;
            loop_back(c_pos);

            while (c_pos != rna->end) {
                if (dna_->shine_dal_start(c_pos)) {
                    rna->start_prot.push_back(c_pos);
                }

                c_pos++;
                loop_back(c_pos);
            }
        }
    }
}

void Organism::compute_protein() {
    int resize_to = 0;

    for (int rna_idx = 0; rna_idx < rna_count_; rna_idx++) {
        resize_to += rnas[rna_idx]->start_prot.size();
    }

    proteins.resize(resize_to);

    for (int rna_idx = 0; rna_idx < rna_count_; rna_idx++) {
        auto* rna = rnas[rna_idx];
        int transcribed_start = rna->begin + PROM_SIZE;
        loop_back(transcribed_start);
        for (int protein_idx = 0; protein_idx < rna->start_prot.size(); protein_idx++) {
            int protein_start = rna->start_prot[protein_idx];
            int current_position = protein_start + SD_TO_START;
            loop_back(current_position);

            int transcription_length;
            if (transcribed_start <= protein_start) {
                transcription_length = protein_start - transcribed_start;
            } else {
                transcription_length = length() - transcribed_start + protein_start;
            }
            transcription_length += SD_TO_START;


            while (rna->length - transcription_length >= CODON_SIZE) {
                if (dna_->protein_stop(current_position)) {
                    int prot_length;

                    int protein_end = current_position + CODON_SIZE;
                    loop_back(protein_end);

                    if (protein_start + SD_TO_START < protein_end) {
                        prot_length = protein_end - (protein_start + SD_TO_START);
                    } else {
                        prot_length = (length() + protein_end) - (protein_start + SD_TO_START);
                    }

                    if (prot_length > CODON_SIZE) { // it has at least 2 codons, among them a STOP
                        int glob_prot_idx = protein_count_;
                        protein_count_ += 1;

                        proteins[glob_prot_idx] =
                                new Protein(protein_start,
                                            protein_end,
                                            prot_length,
                                            rna->e);

                        rna->is_coding_ = true;
                    }
                    break;
                }

                current_position += CODON_SIZE;
                loop_back(current_position);
                transcription_length += CODON_SIZE;
            }
        }
    }
}

void Organism::translate_protein() {
    for (int protein_idx = 0; protein_idx < protein_count_; protein_idx++) {
        auto* protein = proteins[protein_idx];
        if (protein->is_init_) {
            int c_pos = protein->protein_start;
            c_pos += SD_TO_START;
            loop_back(c_pos);

            int nb_codon = (protein->protein_length / 3) - 1; // Do not count the STOP codon
            // Arbitrary limit to the number of codon in one gene
            // It has black magic reasons
            constexpr int FRACTION_SIZE = 52;
            nb_codon = min(nb_codon, FRACTION_SIZE);
            int codon_list[FRACTION_SIZE]{};

            for (int codon_idx = 0; codon_idx < nb_codon; ++codon_idx) {
                codon_list[codon_idx] = dna_->codon_at(c_pos);

                c_pos += 3;
                loop_back(c_pos);
            }

/** This part of the code translate a Gray code binary to standard
 * It looks like black magic (again) but the idea of the implementation can be found hear:
 * https://www.elprocus.com/code-converter-binary-to-gray-code-and-gray-code-to-binary-conversion/ **/
            double M = 0.0;
            double W = 0.0;
            double H = 0.0;

            int nb_m = 0;
            int nb_w = 0;
            int nb_h = 0;

            bool bin_m = false; // Initializing to false will yield a conservation of the high weight bit
            bool bin_w = false; // when applying the XOR operator for the Gray to standard conversion
            bool bin_h = false;


            for (int i = 0; i < nb_codon; i++) {
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
/// End of Black Magic

            protein->protein_length = nb_codon;


            //  ----------------------------------------------------------------------------------
            //  2) Normalize M, W and H values in [0;1] according to number of codons of each kind
            //  ----------------------------------------------------------------------------------
            protein->m = nb_m != 0 ? M / (pow(2, nb_m) - 1) : 0.5;
            protein->w = nb_w != 0 ? W / (pow(2, nb_w) - 1) : 0.0;
            protein->h = nb_h != 0 ? H / (pow(2, nb_h) - 1) : 0.5;

            //  ------------------------------------------------------------------------------------
            //  3) Normalize M, W and H values according to the allowed ranges (defined in macros.h)
            //  ------------------------------------------------------------------------------------
            // x_min <= M <= x_max
            // w_min <= W <= w_max
            // h_min <= H <= h_max
            protein->m = (X_MAX - X_MIN) * protein->m + X_MIN;
            protein->w = (W_MAX - W_MIN) * protein->w + W_MIN;
            protein->h = (H_MAX - H_MIN) * protein->h + H_MIN;

            if (nb_m == 0 || nb_w == 0 || nb_h == 0 ||
                protein->w == 0.0 ||
                protein->h == 0.0) {
                protein->is_functional = false;
            } else {
                protein->is_functional = true;
            }
        }
    }


    std::map<int, Protein *> lookup;

    for (int protein_idx = 0; protein_idx < protein_count_; protein_idx++) {
        auto* protein = proteins[protein_idx];
        if (protein->is_init_) {
            if (lookup.find(protein->protein_start) == lookup.end()) {
                lookup[protein->protein_start] = protein;
            } else {
                lookup[protein->protein_start]->e += protein->e;
                protein->is_init_ = false;
            }
        }
    }
}

void Organism::compute_phenotype() {
    double activ_phenotype[FUZZY_SAMPLING]{};
    double inhib_phenotype[FUZZY_SAMPLING]{};

    for (int protein_idx = 0; protein_idx < protein_count_; protein_idx++) {
        const auto* protein = proteins[protein_idx];
        if (protein->is_init_ && protein->is_functional) {
            // Compute triangle points' coordinates
            double x0 = protein->m - protein->w;
            double x1 = protein->m;
            double x2 = protein->m + protein->w;

            // Interface between continuous world (up) and discrete world (down)
            int ix0 = (int) (x0 * FUZZY_SAMPLING);
            int ix1 = (int) (x1 * FUZZY_SAMPLING);
            int ix2 = (int) (x2 * FUZZY_SAMPLING);

            // active contribution is positive and inhib is negative
            double height = protein->h * protein->e;
            auto* local_phenotype = protein->h > 0 ? activ_phenotype : inhib_phenotype;

            // Compute the first equation of the triangle
            double slope = height / (double)(ix1 - ix0);
            double y_intercept = -(double)ix0 * slope;

            // Updating value between x0 and x1
            for (int i = ix0; i < ix1; i++) {
                if(i >= 0) {
                    local_phenotype[i] += slope * (double)i + y_intercept;
                }
            }

            // Compute the second equation of the triangle
            slope = height / (double)(ix1 - ix2);
            y_intercept = -(double)ix2 * slope;

            // Updating value between x1 and x2
            for (int i = ix1; i < ix2; i++) {
                if(i < FUZZY_SAMPLING) {
                    local_phenotype[i] += slope * (double)i + y_intercept;
                }
            }
        }
    }


    for (int fuzzy_idx = 0; fuzzy_idx < FUZZY_SAMPLING; fuzzy_idx++) {
        if (activ_phenotype[fuzzy_idx] > 1)
            activ_phenotype[fuzzy_idx] = 1;
        if (inhib_phenotype[fuzzy_idx] < -1)
            inhib_phenotype[fuzzy_idx] = -1;
    }

    for (int fuzzy_idx = 0; fuzzy_idx < FUZZY_SAMPLING; fuzzy_idx++) {
        phenotype[fuzzy_idx] = activ_phenotype[fuzzy_idx] + inhib_phenotype[fuzzy_idx];
        if (phenotype[fuzzy_idx] < 0)
            phenotype[fuzzy_idx] = 0;
        if (phenotype[fuzzy_idx] > 1)
            phenotype[fuzzy_idx] = 1;
    }
}

void Organism::compute_fitness(const double *target) {
    metaerror = 0.0;

    for (int fuzzy_idx = 0; fuzzy_idx < FUZZY_SAMPLING; fuzzy_idx++) {
        delta[fuzzy_idx] = fabs(phenotype[fuzzy_idx] - target[fuzzy_idx]);
        delta[fuzzy_idx] /= (double) FUZZY_SAMPLING;
        metaerror += delta[fuzzy_idx];
    }

    fitness = exp(-SELECTION_PRESSURE * ((double) metaerror));
}

/**
 * Switch the DNA base-pair at a given position
 *
 * @param pos : the position where to switch the base-pair
 * @return
 */
bool Organism::do_switch(int pos) {
    dna_->do_switch(pos);

    // Remove promoters containing the switched base
    remove_promoters_around(pos, mod(pos + 1, length()));

    // Look for potential new promoters containing the switched base
    if (length() >= PROM_SIZE)
        look_for_new_promoters_around(pos, mod(pos + 1, length()));

    return true;
}

void Organism::remove_promoters_around(int32_t pos) {
    if (dna_->length() >= PROM_SIZE) {
        remove_promoters_starting_between(mod(pos - PROM_SIZE + 1, dna_->length()),
                                          pos);
    } else {
        remove_all_promoters();
    }
}

void Organism::remove_promoters_around(int32_t pos_1, int32_t pos_2) {
    if (mod(pos_1 - pos_2, dna_->length()) >= PROM_SIZE) {
        remove_promoters_starting_between(mod(pos_1 - PROM_SIZE + 1, dna_->length()),
                                          pos_2);
    } else {
        remove_all_promoters();
    }
}

void Organism::look_for_new_promoters_around(int32_t pos_1, int32_t pos_2) {
    if (dna_->length() >= PROM_SIZE) {
        look_for_new_promoters_starting_between(mod(pos_1 - PROM_SIZE + 1, dna_->length()),
                                                pos_2);
    }
}

void Organism::look_for_new_promoters_around(int32_t pos) {
    if (dna_->length() >= PROM_SIZE) {
        look_for_new_promoters_starting_between(mod(pos - PROM_SIZE + 1, dna_->length()),
                                                pos);
    }
}

void Organism::remove_all_promoters() {
    promoters_.clear();
}

/** REMOVE **/
void Organism::remove_promoters_starting_between(int32_t pos_1, int32_t pos_2) {
    if (pos_1 > pos_2) {
        remove_promoters_starting_after(pos_1);
        remove_promoters_starting_before(pos_2);
    } else {
        // suppression is in [pos1, pos_2[, pos_2 is excluded
        promoters_.erase(promoters_.lower_bound(pos_1), promoters_.upper_bound(pos_2 - 1));
    }
}

void Organism::remove_promoters_starting_after(int32_t pos) {
    promoters_.erase(promoters_.lower_bound(pos), promoters_.end());
}

void Organism::remove_promoters_starting_before(int32_t pos) {
    // suppression is in [0, pos[, pos is excluded
    promoters_.erase(promoters_.begin(), promoters_.upper_bound(pos - 1));
}

void Organism::add_new_promoter(int32_t position, int8_t error) {
    // TODO: Insertion should not always occur, especially if promoter become better or worse ?
    // Promoters are deleted anyway if victim of mutation. the IF stays unnecessary
    if (promoters_.find(position) == promoters_.end())
        promoters_[position] = error;
}

void Organism::look_for_new_promoters_starting_between(int32_t pos_1, int32_t pos_2) {
    // When pos_1 > pos_2, we will perform the search in 2 steps.
    // As positions  0 and dna_->length() are equivalent, it's preferable to
    // keep 0 for pos_1 and dna_->length() for pos_2.

    if (pos_1 >= pos_2) {
        look_for_new_promoters_starting_after(pos_1);
        look_for_new_promoters_starting_before(pos_2);
        return;
    }
    // Hamming distance of the sequence from the promoter consensus

    for (int32_t i = pos_1; i < pos_2; i++) {
        int8_t dist = dna_->promoter_at(i);

        if (dist <= PROM_MAX_DIFF) { // dist takes the hamming distance of the sequence from the consensus
            add_new_promoter(i, dist);
        }
    }
}

void Organism::look_for_new_promoters_starting_after(int32_t pos) {
    for (int32_t i = pos; i < dna_->length(); i++) {
        int dist = dna_->promoter_at(i);

        if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
            add_new_promoter(i, dist);
        }
    }
}

void Organism::look_for_new_promoters_starting_before(int32_t pos) {
    // Hamming distance of the sequence from the promoter consensus

    for (int32_t i = 0; i < pos; i++) {

        int dist = dna_->promoter_at(i);

        if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
            add_new_promoter(i, dist);
        }
    }
}

// Printings

void Organism::print_info() {
    printf("Fitness: %1.10e\n", fitness);
    for (int i = 0; i < protein_count_; i++) {
        const auto& prot = proteins[i];
        printf("%d: %d %f %f %f %f\n", prot->protein_start, prot->is_functional,
               prot->e, prot->w, prot->m, prot->h);
    }
    printf("\nnumber of proteins: %d\n", protein_count_);
    for (int i = 0; i < FUZZY_SAMPLING; ++i) {
        if (phenotype[i] == 0.0){
            printf("0|");
        } else {
            printf("%f|", phenotype[i]);
        }
    }
    printf("\n");
}

