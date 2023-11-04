// Copyright (c) 2013-2023 Anton Kozhevnikov, Thomas Schulthess
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that
// the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the
//    following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
//    and the following disclaimer in the documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/** \file acc.cpp
 *
 *  \brief Definition of the functions for the acc:: namespace.
 *
 */
#include <atomic>
#include "acc.hpp"

namespace sirius {

namespace acc {

int num_devices()
{
#if defined(SIRIUS_CUDA) || defined(SIRIUS_ROCM)
    static std::atomic<int> count(-1);
    if (count.load(std::memory_order_relaxed) == -1) {
        int c;
        if (GPU_PREFIX(GetDeviceCount)(&c) != GPU_PREFIX(Success)) {
            count.store(0, std::memory_order_relaxed);
        } else {
            count.store(c, std::memory_order_relaxed);
        }
    }
    return count.load(std::memory_order_relaxed);
#else
    return 0;
#endif
}

std::vector<acc_stream_t>& streams()
{
    static std::vector<acc_stream_t> streams_;
    return streams_;
}

} // namespace acc

} // namespace sirius