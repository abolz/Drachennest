// Copyright 2019 Alexander Bolz
//
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <cstdint>

namespace dragon4 {

void Dragon4(uint64_t& digits, int& exponent, uint64_t f, int e, bool accept_bounds, bool lower_boundary_is_closer);

} // namespace dragon4
