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
test_pstdout()
{
    mpi::pstdout pout(mpi::Communicator::world());
    pout << "Hello from rank : " << mpi::Communicator::world().rank() << std::endl;

    /* this is a collective operation */
    std::cout << pout.flush(0);

    return 0;
}

int
main(int argn, char** argv)
{
    sirius::initialize(1);
    int result = call_test("test_pstdout", test_pstdout);
    sirius::finalize();
    return result;
}
