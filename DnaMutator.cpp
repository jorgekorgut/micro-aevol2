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

#include "DnaMutator.h"

/**
 * Constructor for DnaMutator class
 *
 * Generate mutations of the DNA of an Organism
 *
 * @param mut_prng : PRNG to simulate the mutation
 * @param length  : Size of the DNA at the initialization
 * @param mutation_rate : Mutation rate of the organisms
 */
DnaMutator::DnaMutator(std::shared_ptr<std::mt19937_64> mut_prng, int length, double mutation_rate) {
    mut_prng_ = mut_prng;
    length_ = length;
    mutation_rate_ = mutation_rate;
}

/**
 * Generate both type of the mutations (see below)
 */
void DnaMutator::generate_mutations() {
    hasMutate_ = false;
    std::binomial_distribution<> d(length_, mutation_rate_);
    nb_mut_ = d(*mut_prng_);

    if (nb_mut_ > 0) {
        for (int i = 0; i < nb_mut_; ++i) {
            generate_next_mutation(length_);
        }
        hasMutate_ = true;
    }
}

/**
 * Generate the next mutation event for an organism.
 *
 * @param length : Update size of the DNA of the Organism
 * @return The generated mutation event (or nullptr if none was created)
 */
MutationEvent *DnaMutator::generate_next_mutation(int length) {
    std::uniform_int_distribution<> distrib(0, length-1);
    int pos = distrib(*mut_prng_);

    auto mevent = new MutationEvent();
    mevent->switch_pos(pos);
    mutation_list_.push_back(mevent);
    return mevent;
};
