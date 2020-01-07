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



#include <cstring>
#include "Organism.h"
#include "ExpManager.h"

#include <iostream>

using namespace std;

/**
 * Constructor to generate a random organism (i.e. an organism with a random DNA)
 *
 * @param exp_m : Related ExpManager object
 * @param length : Length of the generated random DNA
 * @param indiv_id : Unique Identification Number
 */
Organism::Organism(ExpManager *exp_m, int length, int indiv_id) {
    exp_m_ = exp_m;

    rna_count_ = 0;

    auto rng = exp_m->rng_->gen(indiv_id, Threefry::MUTATION);
    dna_ = new Dna(length, rng);
    parent_length_ = length;
    indiv_id_ = indiv_id;
}

/**
 * Create an organism with a given genome
 *
 * @param exp_m : Related ExpManager object
 * @param genome : Genome to assign to the organism
 * @param indiv_id : Unique Identification Number
 */
Organism::Organism(ExpManager *exp_m, char *genome, int indiv_id) {
    exp_m_ = exp_m;

    rna_count_ = 0;

    dna_ = new Dna(genome, strlen(genome));
    parent_length_ = strlen(genome);
    indiv_id_ = indiv_id;

}

/**
 * Constructor to create a clone of a given Organism
 *
 * @param exp_m : Related ExpManager object
 * @param clone : The organism to clone
 */
Organism::Organism(ExpManager *exp_m, const std::shared_ptr<Organism> &clone) {
    exp_m_ = exp_m;

    rna_count_ = 0;

    parent_length_ = clone->length();
    dna_ = new Dna(*(clone->dna_));

    promoters_ = clone->promoters_;
}

/**
 * Create an Organism from a backup/checkpointing file
 *
 * @param exp_m : Related ExpManager object
 * @param backup_file : gzFile to read from
 */
Organism::Organism(ExpManager *exp_m, gzFile backup_file) {
    exp_m_ = exp_m;

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
void Organism::save(gzFile backup_file) {
    gzwrite(backup_file, &indiv_id_, sizeof(indiv_id_));
    gzwrite(backup_file, &parent_id_, sizeof(parent_id_));
    gzwrite(backup_file, &global_id, sizeof(global_id));

    gzwrite(backup_file, &parent_length_, sizeof(parent_length_));

    dna_->save(backup_file);
}

/**
 * Load the organism from backup/checkpointing file
 *
 * @param backup_file : from where restore the organism
 */
void Organism::load(gzFile backup_file) {
    gzread(backup_file, &indiv_id_, sizeof(indiv_id_));
    gzread(backup_file, &parent_id_, sizeof(parent_id_));
    gzread(backup_file, &global_id, sizeof(global_id));

    gzread(backup_file, &parent_length_, sizeof(parent_length_));

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

/**
 * Apply all the mutation events of the organism on its DNA
 */
void Organism::apply_mutations() {
    auto mutation_list = exp_m_->dna_mutator_array_[indiv_id_]->mutation_list_;

    for (const auto mutation: mutation_list) {
        switch (mutation->type()) {
            case DO_SWITCH:
                do_switch(mutation->pos_1());
                nb_swi_++;
                nb_mut_++;
                break;
        }
    }

}

/**
Optimize promoters search
 **/


void Organism::remove_promoters_around(int32_t pos) {
    if (dna_->length() >= PROM_SIZE) {
        remove_promoters_starting_between(mod(pos - PROM_SIZE + 1,
                                              dna_->length()),
                                          pos);
    } else {
        remove_all_promoters();
    }
}

void Organism::remove_promoters_around(int32_t pos_1, int32_t pos_2) {
    if (mod(pos_1 - pos_2, dna_->length()) >= PROM_SIZE) {
        remove_promoters_starting_between(mod(pos_1 - PROM_SIZE + 1,
                                              dna_->length()),
                                          pos_2);
    } else {
        remove_all_promoters();
    }
}

void Organism::look_for_new_promoters_around(int32_t pos_1, int32_t pos_2) {
    if (dna_->length() >= PROM_SIZE) {
        look_for_new_promoters_starting_between(
                mod(pos_1 - PROM_SIZE + 1,
                    dna_->length()), pos_2);
    }
}

void Organism::look_for_new_promoters_around(int32_t pos) {
    if (dna_->length() >= PROM_SIZE) {
        look_for_new_promoters_starting_between(
                mod(pos - PROM_SIZE + 1, dna_->length()),
                pos);
    }
}

void Organism::remove_all_promoters() {
    promoters_.clear();
}

/** LEADING promoters **/
/** REMOVE **/
void Organism::remove_promoters_starting_between(int32_t pos_1, int32_t pos_2) {
    if (pos_1 > pos_2) {
        remove_promoters_starting_after(pos_1);
        remove_promoters_starting_before(pos_2);
    } else {
        // suppression is in [pos1, pos_2[, pos_2 is excluded
        promoters_.erase(promoters_.lower_bound(pos_1), promoters_.upper_bound(pos_2-1));
    }
}

void Organism::remove_promoters_starting_after(int32_t pos) {
    promoters_.erase(promoters_.lower_bound(pos), promoters_.end());
}

void Organism::remove_promoters_starting_before(int32_t pos) {
    // suppression is in [0, pos[, pos is excluded
    promoters_.erase(promoters_.begin(), promoters_.upper_bound(pos-1));
}


/** LOOK **/
void Organism::locate_promoters() {
    look_for_new_promoters_starting_between(0, dna_->length());
}

void Organism::add_new_promoter(int32_t position, int8_t error) {
    // TODO: Insertion should not always occur, especially if promoter become better or worse ?
    // Promoters are deleted anyway if victim of mutation. the IF stays unnecessary
    if(promoters_.find(position) == promoters_.end())
        promoters_[position] = Promoter(position, error);
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

        if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
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

