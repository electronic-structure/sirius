// Copyright (c) 2013-2018 Anton Kozhevnikov, Thomas Schulthess
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

/** \file field4d.cpp
 *
 *  \brief Implementation of sirius::Field4D class.
 */

#include "field4d.hpp"
#include "periodic_function.hpp"
#include "symmetry/symmetrize.hpp"

namespace sirius {

void Field4D::symmetrize(Periodic_function<double>* f__, Periodic_function<double>* gz__,
                         Periodic_function<double>* gx__, Periodic_function<double>* gy__)
{
    PROFILE("sirius::Field4D::symmetrize");

    /* quick exit: the only symmetry operation is identity */
    if (ctx_.unit_cell().symmetry().size() == 1) {
        return;
    }

    auto& comm = ctx_.comm();

    if (ctx_.cfg().control().print_hash()) {
        auto h = f__->rg().hash_f_pw();
        if (ctx_.comm().rank() == 0) {
            utils::print_hash("f_unsymmetrized(G)", h);
        }
    }

    /* symmetrize PW components */
    switch (ctx_.num_mag_dims()) {
        case 0: {
            sirius::symmetrize(ctx_.unit_cell().symmetry(), ctx_.remap_gvec(), ctx_.sym_phase_factors(),
                &f__->rg().f_pw_local(0), nullptr, nullptr, nullptr);
            if (ctx_.cfg().control().print_hash()) {
                auto h = f__->rg().hash_f_pw();
                if (ctx_.comm().rank() == 0) {
                    utils::print_hash("f_symmetrized(G)", h);
                }
            }
            break;
        }
        case 1: {
            sirius::symmetrize(ctx_.unit_cell().symmetry(), ctx_.remap_gvec(), ctx_.sym_phase_factors(),
                &f__->rg().f_pw_local(0), nullptr, nullptr, &gz__->rg().f_pw_local(0));
            break;
        }
        case 3: {
            if (ctx_.cfg().control().print_hash()) {
                auto h1 = gx__->rg().hash_f_pw();
                auto h2 = gy__->rg().hash_f_pw();
                auto h3 = gz__->rg().hash_f_pw();
                if (ctx_.comm().rank() == 0) {
                    utils::print_hash("fx_unsymmetrized(G)", h1);
                    utils::print_hash("fy_unsymmetrized(G)", h2);
                    utils::print_hash("fz_unsymmetrized(G)", h3);
                }
            }

            sirius::symmetrize(ctx_.unit_cell().symmetry(), ctx_.remap_gvec(), ctx_.sym_phase_factors(),
                &f__->rg().f_pw_local(0), &gx__->rg().f_pw_local(0),
                &gy__->rg().f_pw_local(0), &gz__->rg().f_pw_local(0));

            if (ctx_.cfg().control().print_hash()) {
                auto h1 = gx__->rg().hash_f_pw();
                auto h2 = gy__->rg().hash_f_pw();
                auto h3 = gz__->rg().hash_f_pw();
                if (ctx_.comm().rank() == 0) {
                    utils::print_hash("fx_symmetrized(G)", h1);
                    utils::print_hash("fy_symmetrized(G)", h2);
                    utils::print_hash("fz_symmetrized(G)", h3);
                }
            }
            break;
        }
    }

    if (ctx_.full_potential()) {
        /* symmetrize MT components */
        symmetrize_function(ctx_.unit_cell().symmetry(), comm, f__->f_mt());
        switch (ctx_.num_mag_dims()) {
            case 1: {
                symmetrize_vector_function(ctx_.unit_cell().symmetry(), comm, gz__->f_mt());
                break;
            }
            case 3: {
                symmetrize_vector_function(ctx_.unit_cell().symmetry(), comm, gx__->f_mt(), gy__->f_mt(), gz__->f_mt());
                break;
            }
        }
    }
}

Field4D::Field4D(Simulation_context& ctx__, int lmmax__)
    : lmmax_(lmmax__)
    , ctx_(ctx__)
{
    for (int i = 0; i < ctx_.num_mag_dims() + 1; i++) {
        if (ctx_.full_potential()) {
            components_[i] = std::make_unique<Periodic_function<double>>(ctx_, lmmax__);
            /* allocate global MT array */
            components_[i]->allocate_mt(true);
        } else {
            components_[i] = std::make_unique<Periodic_function<double>>(ctx_, lmmax__);
        }
    }
}

Periodic_function<double>& Field4D::scalar()
{
    return *(components_[0]);
}

Periodic_function<double> const& Field4D::scalar() const
{
    return *(components_[0]);
}

void Field4D::zero()
{
    for (int i = 0; i < ctx_.num_mag_dims() + 1; i++) {
        component(i).zero();
    }
}

void Field4D::fft_transform(int direction__)
{
    for (int i = 0; i < ctx_.num_mag_dims() + 1; i++) {
        component(i).rg().fft_transform(direction__);
    }
}

} // namespace sirius
