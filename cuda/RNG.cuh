//
// Created by elturpin on 08/12/2020.
//

#pragma once

#include "Random123/threefry.h"

using RNG = r123::Threefry2x64;
using ctr_type = RNG::ctr_type;
using key_type = RNG::key_type;
using ctr_value_type = ctr_type::value_type;
using key_value_type = key_type::value_type;
