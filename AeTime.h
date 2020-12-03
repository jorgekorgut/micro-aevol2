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

class AeTime {
public :
    // =================================================================
    //                             Constructors
    // =================================================================
    AeTime() = delete; //< Default ctor
    AeTime(const AeTime &) = delete; //< Copy ctor
    AeTime(AeTime &&) = delete; //< Move ctor

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~AeTime() = delete;

    // =================================================================
    //                        Accessors: getters
    // =================================================================
    static inline int time() { return time_; }

    // =================================================================
    //                        Accessors: setters
    // =================================================================
    static inline void set_time(int t) { time_ = t; }

    // =================================================================
    //                              Operators
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    static inline void plusplus() { time_++; }

    // =================================================================
    //                           Public Attributes
    // =================================================================

protected :
    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    static int time_;
};
