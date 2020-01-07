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


#ifndef PDC_MINI_AEVOL_STATS_H
#define PDC_MINI_AEVOL_STATS_H

#include <cstdint>
#include <fstream>
#include <limits>

#include "Organism.h"
#include "ExpManager.h"

/**
 * Class to manage and generate the Stats (and the related file) of a simulation
 */
class Stats {
public:
    Stats(ExpManager* exp_m, int generation, bool best_or_not);

    ~Stats() {
        if (is_indiv_) {
            statfile_best_.flush();
            statfile_best_.close();
        } else {
            statfile_mean_.flush();
            statfile_mean_.close();
        }
    }

    void compute_best();
    void compute_average();

    void write_best();
    void write_average();

    void reinit(int generation);

    bool is_indiv() { return is_indiv_; }


protected:
    ExpManager* exp_m_;

    int generation_;

    bool is_indiv_;

    int pop_size_;

    double fitness_ = 0;
    double mean_fitness_ = 0;
    double metabolic_error_ = 0;
    double mean_metabolic_error_ = 0;

    int amount_of_dna_ = 0;
    float mean_amount_of_dna_ = 0;
    int nb_coding_rnas_ = 0;
    float mean_nb_coding_rnas_ = 0;
    int nb_non_coding_rnas_ = 0;
    float mean_nb_non_coding_rnas_ = 0;

    int nb_functional_genes_ = 0;
    float mean_nb_functional_genes_ = 0;
    int nb_non_functional_genes_ = 0;
    float mean_nb_non_functional_genes_ = 0;

    int nb_mut_ = 0;
    float mean_nb_mut_ = 0;
    int nb_switch_ = 0;
    float mean_nb_switch_ = 0;

    bool is_computed_ = false;

    // Stats

    std::ofstream statfile_best_;
    std::ofstream statfile_mean_;
};


#endif //PDC_MINI_AEVOL_STATS_H
