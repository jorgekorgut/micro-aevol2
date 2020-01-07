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


#ifndef PDC_MINI_AEVOL_RNA_H
#define PDC_MINI_AEVOL_RNA_H

#include <vector>

/**
 * Class to store a RNA and its related variable
 */
class RNA {
public:
    RNA() {};
    RNA(int t_begin, int t_end, double t_e, int t_length) {
        begin = t_begin;
        end = t_end;
        e = t_e;
        length = t_length;
        is_coding_ = false;
        // is_init_ = true;
    }

    int begin;
    int end;

    double e;
    std::vector<int> start_prot;
    int length;
    bool is_coding_;

    // in the code is_init_ is always true and never turn to false
    // bool is_init_ = false;
};

#endif //PDC_MINI_AEVOL_RNA_H
