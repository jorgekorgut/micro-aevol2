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



#ifndef PDC_MINI_AEVOL_PROMOTER_H
#define PDC_MINI_AEVOL_PROMOTER_H

/**
 * Class to store a promoter located on the DNA of an Organism
 */
class Promoter {
public:

    //necessary to use them in std::map
    Promoter() = default;

    Promoter(int t_pos, int t_error) {
        pos = t_pos; error = t_error;
    }

    Promoter(const Promoter& clone) {
        pos=clone.pos;error=clone.error;
    }

    Promoter(Promoter* clone) {
        pos=clone->pos;error=clone->error;
    }

    int pos = -1;
    int error = -1;
};


#endif //PDC_MINI_AEVOL_PROMOTER_H
