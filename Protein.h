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

/**
 * Class to store a protein and its related variable
 */
class Protein {
public:
    Protein() = default;
    Protein(int t_protein_start, int t_protein_end, int t_protein_length, double t_e) {
        protein_start = t_protein_start;
        protein_end = t_protein_end;
        protein_length = t_protein_length;
        e = t_e;
        is_init_ = true;
    }

    int protein_start = 0;
    int protein_end = 0;
    int protein_length = 0;

    double m = 0.0;
    double w = 0.0;
    double h = 0.0;
    double e = 0.0;
    bool is_functional = false;

    bool is_init_ = false;
};
