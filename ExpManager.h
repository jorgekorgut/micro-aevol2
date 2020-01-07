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


#ifndef PDC_MINI_AEVOL_EXPMANAGER_H
#define PDC_MINI_AEVOL_EXPMANAGER_H


#include "Threefry.h"
#include "DnaMutator.h"
#include "Stats.h"


#include <memory>

constexpr int8_t NB_BASE = 2;


constexpr int PROM_SIZE = 22;

constexpr int8_t CODON_START = 0b000;
constexpr int8_t CODON_M0    = 0b100;
constexpr int8_t CODON_M1    = 0b101;
constexpr int8_t CODON_W0    = 0b010;
constexpr int8_t CODON_W1    = 0b011;
constexpr int8_t CODON_H0    = 0b110;
constexpr int8_t CODON_H1    = 0b111;

constexpr double X_MIN = 0.0;
constexpr double X_MAX = 1.0;
constexpr double Y_MIN = 0.0;
constexpr double Y_MAX = 1.0;
constexpr double H_MIN = -1.0;
constexpr double H_MAX = 1.0;
constexpr double W_MIN = 0.0;

class Organism;
class Stats;


/**
 * Main class of the simulator.
 * ExpManager is in charge of running the simulation and maintaining all the data.
 * It is also that class that implements checkpointing and restore mechanisms.
 */
class ExpManager {


    public:
        ExpManager(int grid_height, int grid_width, int seed, double mutation_rate, int init_length_dna,
                   double w_max, int selection_pressure, int backup_step);
        ExpManager(int time);
        ~ExpManager();

        void create_directory();
        void save(int t);

        void load(int t);

        void run_evolution(int nb_gen);

#ifdef USE_CUDA
        void run_evolution_on_gpu(int nb_gen);
#endif

        void run_a_step(double w_max, double selection_pressure, bool first_gen);

        void prepare_mutation(int indiv_id);
        void selection(int indiv_id);

        inline void apply_mutation(int indiv_id) { internal_organisms_[indiv_id]->apply_mutations(); }

        void start_stop_RNA(int indiv_id);
        void compute_RNA(int indiv_id);

        void opt_prom_compute_RNA(int indiv_id);

        void start_protein(int indiv_id);
        void compute_protein(int indiv_id);

        void translate_protein(int indiv_id, double w_max);

        void compute_phenotype(int indiv_id);
        void compute_fitness(int indiv_id, double selection_pressure);

        std::shared_ptr<Organism>* internal_organisms_;
        std::shared_ptr<Organism>* prev_internal_organisms_;
        std::shared_ptr<Organism> best_indiv;

        int* next_generation_reproducer_;
        DnaMutator** dna_mutator_array_;

        int nb_indivs_;

        std::unique_ptr<Threefry> rng_;

        double geometric_area_;

        double* target;

        int selection_pressure_;
    //private:
        Stats* stats_best = nullptr;
        Stats* stats_mean = nullptr;


        int grid_height_;
        int grid_width_;

        double mutation_rate_;

        double w_max_;

        int backup_step_;
};


#endif //PDC_MINI_AEVOL_EXPMANAGER_H
