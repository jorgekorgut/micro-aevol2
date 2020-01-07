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


#ifndef PDC_MINI_AEVOL_ORGANISM_H
#define PDC_MINI_AEVOL_ORGANISM_H

#include <map>
#include <memory>
#include <set>
#include <zlib.h>
#include <list>

#include "Promoter.h"
#include "RNA.h"
#include "Protein.h"
#include "Dna.h"

class ExpManager;

/**
 * Class that implements an organism and its related DNAs, RNAs, Protein and Phenotype
 */
class Organism {

public:
    Organism(ExpManager *exp_m, int length, int indiv_id);

    /// Create an organism with a given genome
    Organism(ExpManager *exp_m, char *genome, int indiv_id);

    Organism(ExpManager *exp_m, const std::shared_ptr<Organism> &clone);

    Organism(ExpManager *exp_m, gzFile backup_file);

    ~Organism();

    void save(gzFile backup_file);

    void load(gzFile backup_file);

    int length() { return dna_->length(); };

    void apply_mutations();

    void reset_mutation_stats();

    void compute_protein_stats();


    // Map position (int) to Promoter
    std::map<int, Promoter> promoters_;

    std::set<int> terminators;
    std::vector<RNA *> rnas;
    std::vector<Protein *> proteins;

    double phenotype[300];
    double delta[300];

    double fitness;
    double metaerror;

    Dna *dna_;
    int parent_length_;

    int indiv_id_;
    int parent_id_;

    int protein_count_ = 0;
    int rna_count_ = 0;

    ExpManager *exp_m_;

    int global_id = -1;

    int usage_count_ = 1;

    // Stats
    int nb_genes_activ = 0;
    int nb_genes_inhib = 0;
    int nb_func_genes = 0;
    int nb_non_func_genes = 0;
    int nb_coding_RNAs = 0;
    int nb_non_coding_RNAs = 0;

    int nb_swi_ = 0;
    int nb_mut_ = 0;

//private:

    bool do_switch(int pos);

    void remove_all_promoters();

    void remove_promoters_around(int32_t pos);

    void remove_promoters_around(int32_t pos_1, int32_t pos_2);

    void remove_promoters_starting_between(int32_t pos_1, int32_t pos_2);

    void remove_promoters_starting_after(int32_t pos);

    void remove_promoters_starting_before(int32_t pos);

    void locate_promoters();

    void look_for_new_promoters_around(int32_t pos_1, int32_t pos_2);

    void look_for_new_promoters_around(int32_t pos);

    void look_for_new_promoters_starting_between(int32_t pos_1, int32_t pos_2);

    void look_for_new_promoters_starting_after(int32_t pos);

    void look_for_new_promoters_starting_before(int32_t pos);

    void add_new_promoter(int32_t position, int8_t error);

    inline int32_t mod(int32_t a, int32_t b) {

        assert(b > 0);

        while (a < 0) a += b;
        while (a >= b) a -= b;

        return a;
        //return m >= 0 ? m % n : ( n - abs ( m%n ) ) % n;
    }

    inline int64_t mod(int64_t a, int64_t b) {

        assert(b > 0);

        while (a < 0) a += b;
        while (a >= b) a -= b;

        return a;
        //return m >= 0 ? m % n : ( n - abs ( m%n ) ) % n;
    }
};


#endif //PDC_MINI_AEVOL_ORGANISM_H
