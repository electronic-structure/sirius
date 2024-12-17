/* This file is part of SIRIUS electronic structure library.
 *
 * Copyright (c), ETH Zurich.  All rights reserved.
 *
 * Please, refer to the LICENSE file in the root directory.
 * SPDX-License-Identifier: BSD-3-Clause
 */

#include <sirius.hpp>
#include <testing.hpp>

using namespace sirius;

int
test_roundoff(cmd_args const& args)
{
    int err{0};
    if (std::abs(round(1.525252525, 4) - 1.5253) > 1e-20) {
        err++;
    }
    if (std::abs(round(2.12345678, 4) - 2.1235) > 1e-20) {
        err++;
    }
    if (std::abs(round(2.12344678, 4) - 2.1234) > 1e-20) {
        err++;
    }
    if (std::abs(round(1.999, 0) - 2) > 1e-20) {
        err++;
    }
    if (std::abs(round(1.999, 1) - 2) > 1e-20) {
        err++;
    }
    return err;
}

int
main(int argn, char** argv)
{
    cmd_args args;
    return call_test(argv[0], test_roundoff, args);
}
