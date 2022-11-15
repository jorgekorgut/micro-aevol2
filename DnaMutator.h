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

#pragma once

#include <list>
#include<bits/stdc++.h>

#include "Organism.h"
#include "MutationEvent.h"


/**
 * Class that generates the mutation events for a given Organism
 */
class DnaMutator {
public:

    DnaMutator(std::shared_ptr<std::mt19937_64> rng, int length, double mutation_rate);

    ~DnaMutator() {
        for (auto repl : mutation_list_) {
            delete repl;
        }
        mutation_list_.clear();
        //delete mut_prng_;
    }

    void generate_mutations();

    MutationEvent *generate_next_mutation(int length);

    std::list<MutationEvent *> mutation_list_;

    bool hasMutate() const { return hasMutate_; }

    void setMutate(bool mutate) { hasMutate_ = mutate; }

    // private:
    std::shared_ptr<std::mt19937_64> mut_prng_;
    int length_;

    double mutation_rate_;

    //--------------------------- Mutation counters
    int nb_mut_{};
    bool hasMutate_ = false;
};

