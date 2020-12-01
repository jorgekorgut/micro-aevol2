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


#include "Stats.h"

/**
 * Create the data structure and open the file for the statistics of a simulation
 *
 * @param generation : Create statistics beginning from this generation (or resuming from this generation)
 * @param best_or_not : Statistics for the best organisms or mean of all the organisms
 */
Stats::Stats(int generation, bool best_or_not) {
    is_indiv_ = best_or_not;
    generation_ = generation;

    pop_size_ = 0;

    fitness_ = 0;
    metabolic_error_ = 0;

    amount_of_dna_ = 0;
    nb_coding_rnas_ = 0;
    nb_non_coding_rnas_ = 0;

    nb_functional_genes_ = 0;
    nb_non_functional_genes_ = 0;

    nb_mut_ = 0;
    nb_switch_ = 0;

    if (generation_ == 0) {
        if (is_indiv_)
            statfile_best_.open("stats/stats_simd_best.csv", std::ofstream::trunc);
        else
            statfile_mean_.open("stats/stats_simd_mean.csv", std::ofstream::trunc);

        if (is_indiv_) {
            statfile_best_ << "Generation" << "," << "fitness" << "," << "metabolic_error" << "," <<
                           "amount_of_dna" << "," << "nb_coding_rnas" << "," << "nb_non_coding_rnas" << "," <<
                           "nb_functional_genes" << "," << "nb_non_functional_genes" << "," << "nb_mut"
                           << "," << "nb_switch" << std::endl;
            statfile_best_.flush();
        } else {
            statfile_mean_ << "Generation" << "," << "fitness" << "," << "metabolic_error" << "," <<
                           "amount_of_dna" << "," << "nb_coding_rnas" << "," << "nb_non_coding_rnas" << "," <<
                           "nb_functional_genes" << "," << "nb_non_functional_genes" << "," << "nb_mut"
                           << "," << "nb_switch" << std::endl;
            statfile_mean_.flush();
        }
    } else {
        printf("Resume without rheader\n");
        std::ifstream tmp_mean;
        std::ifstream tmp_best;

        if (is_indiv_) {
            tmp_best.open("stats/stats_simd_best.csv", std::ifstream::in);
            statfile_best_.open("stats/stats_simd_best.csv.tmp", std::ofstream::trunc);
        } else {
            tmp_mean.open("stats/stats_simd_mean.csv", std::ifstream::in);
            statfile_mean_.open("stats/stats_simd_mean.csv.tmp", std::ofstream::trunc);
        }

        std::string str;
        for (int i = 0; i < generation_; i++) {
            if (is_indiv_) {
                std::getline(tmp_best, str);
                statfile_best_ << str << std::endl;
            } else {
                std::getline(tmp_mean, str);
                statfile_mean_ << str << std::endl;
            }
        }

        if (is_indiv_) {
            statfile_best_.flush();
            statfile_best_.close();
        } else {
            statfile_mean_.flush();
            statfile_mean_.close();
        }

        if (is_indiv_) {
            statfile_best_.open("stats/stats_simd_best.csv", std::ofstream::trunc);
            tmp_best.close();
            tmp_best.open("stats/stats_simd_best.csv.tmp", std::ifstream::in);
            tmp_best.seekg(0, std::ios::beg);
        } else {
            statfile_mean_.open("stats/stats_simd_mean.csv", std::ofstream::trunc);
            tmp_mean.close();
            tmp_mean.open("stats/stats_simd_mean.csv.tmp", std::ifstream::in);
            tmp_mean.seekg(0, std::ios::beg);
        }

        for (int i = 0; i < generation_; i++) {
            if (is_indiv_) {
                std::getline(tmp_best, str);
                statfile_best_ << str << std::endl;
            } else {
                std::getline(tmp_mean, str);
                statfile_mean_ << str << std::endl;
            }
        }

        tmp_best.close();
        tmp_mean.close();
        is_indiv_
        ? std::remove("stats/stats_simd_best.csv.tmp")
        : std::remove("stats/stats_simd_mean.csv.tmp");
    }
}

/**
 * Compute the statistics for the best organism
 */
void Stats::compute_best(const std::shared_ptr<Organism> &best) {
    is_indiv_ = true;

    fitness_ = best->fitness;
    metabolic_error_ = best->metaerror;

    amount_of_dna_ = best->length();

    nb_coding_rnas_ = best->nb_coding_RNAs;
    nb_non_coding_rnas_ = best->nb_non_coding_RNAs;

    nb_functional_genes_ = best->nb_func_genes;
    nb_non_functional_genes_ = best->nb_non_func_genes;


    nb_mut_ = best->nb_mut_;
    nb_switch_ = best->nb_swi_;

    is_computed_ = true;
}

/**
 * Compute the statistics of the mean of the whole population
 */
void Stats::compute_average(const std::shared_ptr<Organism> *population, int population_size) {
    is_indiv_ = false;
    pop_size_ = population_size;

    mean_fitness_ = 0;
    mean_metabolic_error_ = 0;
    mean_amount_of_dna_ = 0;
    mean_nb_coding_rnas_ = 0;
    mean_nb_non_coding_rnas_ = 0;
    mean_nb_functional_genes_ = 0;
    mean_nb_non_functional_genes_ = 0;

    mean_nb_mut_ = 0;
    mean_nb_switch_ = 0;

    for (int indiv_id = 0; indiv_id < pop_size_; indiv_id++) {
        const auto &organism = population[indiv_id];
        mean_fitness_ += organism->fitness;
        mean_metabolic_error_ += organism->metaerror;

        mean_amount_of_dna_ += organism->length();

        mean_nb_coding_rnas_ += organism->nb_coding_RNAs;
        mean_nb_non_coding_rnas_ += organism->nb_non_coding_RNAs;

        mean_nb_functional_genes_ += organism->nb_func_genes;
        mean_nb_non_functional_genes_ += organism->nb_non_func_genes;

        mean_nb_mut_ += organism->nb_mut_;
        mean_nb_switch_ += organism->nb_swi_;
    }


    //printf("pop_size : %d \n",pop_size_);

    mean_fitness_ /= pop_size_;
    mean_metabolic_error_ /= pop_size_;

    mean_amount_of_dna_ /= pop_size_;
    mean_nb_coding_rnas_ /= pop_size_;
    mean_nb_non_coding_rnas_ /= pop_size_;

    mean_nb_functional_genes_ /= pop_size_;
    mean_nb_non_functional_genes_ /= pop_size_;

    mean_nb_mut_ /= pop_size_;
    mean_nb_switch_ /= pop_size_;

    is_computed_ = true;
}

/**
 * Write the statistics of the best organism to the related statistics file
 */
void Stats::write_best(const std::shared_ptr<Organism> &best) {
    if (is_indiv_ && !is_computed_)
        compute_best(best);

    if (is_indiv_ && is_computed_) {
        // Write best stats
        statfile_best_ << generation_ << "," << fitness_ << "," << metabolic_error_ << "," <<
                       amount_of_dna_ << "," << nb_coding_rnas_ << "," << nb_non_coding_rnas_ << "," <<
                       nb_functional_genes_ << "," << nb_non_functional_genes_ << "," << nb_mut_
                       << "," << nb_switch_ << std::endl;
        statfile_best_.flush();
    }
}

/**
 * Write the statistics of the mean of the population to the related statistics file
 */
void Stats::write_average(const std::shared_ptr<Organism> *population, int population_size) {
    if (!is_indiv_ && !is_computed_)
        compute_average(population, population_size);

    if (!is_indiv_ && is_computed_) {
        // Write average stats
        statfile_mean_ << generation_ << "," << mean_fitness_ << "," << mean_metabolic_error_ << "," <<
                       mean_amount_of_dna_ << "," << mean_nb_coding_rnas_ << "," << mean_nb_non_coding_rnas_ << "," <<
                       mean_nb_functional_genes_ << "," << mean_nb_non_functional_genes_ << "," << mean_nb_mut_
                       << "," << mean_nb_switch_ << std::endl;
        statfile_mean_.flush();
    }
}

/**
 * Reinitilized the statistics variable before computing next generation
 */
void Stats::reinit(int generation) {
    generation_ = generation;

    pop_size_ = 0;

    fitness_ = 0;
    metabolic_error_ = 0;

    amount_of_dna_ = 0;
    nb_coding_rnas_ = 0;
    nb_non_coding_rnas_ = 0;

    nb_functional_genes_ = 0;
    nb_non_functional_genes_ = 0;

    nb_mut_ = 0;
    nb_switch_ = 0;

    is_computed_ = false;
}
