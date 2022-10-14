// Copyright (c) 2013-2019 Anton Kozhevnikov, Mathieu Taillefumier, Thomas Schulthess
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

/** \file local_operator.cpp
 *
 *  \brief Implementation of sirius::Local_operator class.
 */

#include "local_operator.hpp"
#include "potential/potential.hpp"
#include "function3d/smooth_periodic_function.hpp"
#include "utils/profiler.hpp"

namespace sirius {

template <typename T>
Local_operator<T>::Local_operator(Simulation_context const& ctx__, spfft_transform_type<T>& fft_coarse__,
                                  std::shared_ptr<sddk::Gvec_fft> gvec_coarse_p__, Potential* potential__)
    : ctx_(ctx__)
    , fft_coarse_(fft_coarse__)
    , gvec_coarse_p_(gvec_coarse_p__)

{
    PROFILE("sirius::Local_operator");

    /* allocate functions */
    for (int j = 0; j < ctx_.num_mag_dims() + 1; j++) {
        veff_vec_[j] = std::unique_ptr<Smooth_periodic_function<T>>(
            new Smooth_periodic_function<T>(fft_coarse__, gvec_coarse_p__, &ctx_.mem_pool(sddk::memory_t::host)));
        #pragma omp parallel for schedule(static)
        for (int ir = 0; ir < fft_coarse__.local_slice_size(); ir++) {
            veff_vec_[j]->f_rg(ir) = 2.71828;
        }
    }
    /* map Theta(r) to the coarse mesh */
    if (ctx_.full_potential()) {
        auto& gvec_dense_p = ctx_.gvec_fft();
        veff_vec_[v_local_index_t::theta] = std::unique_ptr<Smooth_periodic_function<T>>(
            new Smooth_periodic_function<T>(fft_coarse__, gvec_coarse_p__, &ctx_.mem_pool(sddk::memory_t::host)));
        /* map unit-step function */
        #pragma omp parallel for schedule(static)
        for (int igloc = 0; igloc < gvec_coarse_p_->gvec().count(); igloc++) {
            /* map from fine to coarse set of G-vectors */
            veff_vec_[v_local_index_t::theta]->f_pw_local(igloc) =
                ctx_.theta_pw(gvec_dense_p.gvec().gvec_base_mapping(igloc) + gvec_dense_p.gvec().offset());
        }
        veff_vec_[v_local_index_t::theta]->fft_transform(1);
        if (fft_coarse_.processing_unit() == SPFFT_PU_GPU) {
            veff_vec_[v_local_index_t::theta]->f_rg().allocate(ctx_.mem_pool(sddk::memory_t::device)).copy_to(sddk::memory_t::device);
        }
        if (ctx_.print_checksum()) {
            auto cs1 = veff_vec_[v_local_index_t::theta]->checksum_pw();
            auto cs2 = veff_vec_[v_local_index_t::theta]->checksum_rg();
            if (ctx_.comm().rank() == 0) {
                utils::print_checksum("theta_pw", cs1);
                utils::print_checksum("theta_rg", cs2);
            }
        }
    }

    /* map potential */
    if (potential__) {

        if (ctx_.full_potential()) {

            auto& fft_dense    = ctx_.spfft<T>();
            auto& gvec_dense_p = ctx_.gvec_fft();

            Smooth_periodic_function<T> ftmp(const_cast<Simulation_context&>(ctx_).spfft<T>(), ctx_.gvec_fft_sptr(),
                                             &ctx_.mem_pool(sddk::memory_t::host));

            for (int j = 0; j < ctx_.num_mag_dims() + 1; j++) {
                /* multiply potential by step function theta(r) */
                for (int ir = 0; ir < fft_dense.local_slice_size(); ir++) {
                    ftmp.f_rg(ir) = potential__->component(j).f_rg(ir) * ctx_.theta(ir);
                }
                /* transform to plane-wave domain */
                ftmp.fft_transform(-1);
                if (j == 0) {
                    v0_[0] = ftmp.f_0().real();
                }
                /* loop over local set of coarse G-vectors */
                #pragma omp parallel for schedule(static)
                for (int igloc = 0; igloc < gvec_coarse_p_->gvec().count(); igloc++) {
                    /* map from fine to coarse set of G-vectors */
                    veff_vec_[j]->f_pw_local(igloc) = ftmp.f_pw_local(gvec_dense_p.gvec().gvec_base_mapping(igloc));
                }
                /* transform to real space */
                veff_vec_[j]->fft_transform(1);
            }
            if (ctx_.valence_relativity() == relativity_t::zora) {
                veff_vec_[v_local_index_t::rm_inv] = std::unique_ptr<Smooth_periodic_function<T>>(
                    new Smooth_periodic_function<T>(fft_coarse__, gvec_coarse_p__, &ctx_.mem_pool(sddk::memory_t::host)));
                /* loop over local set of coarse G-vectors */
                #pragma omp parallel for schedule(static)
                for (int igloc = 0; igloc < gvec_coarse_p_->gvec().count(); igloc++) {
                    /* map from fine to coarse set of G-vectors */
                    veff_vec_[v_local_index_t::rm_inv]->f_pw_local(igloc) =
                        potential__->rm_inv_pw(gvec_dense_p.gvec().offset() + gvec_dense_p.gvec().gvec_base_mapping(igloc));
                }
                /* transform to real space */
                veff_vec_[v_local_index_t::rm_inv]->fft_transform(1);
            }

        } else {

            for (int j = 0; j < ctx_.num_mag_dims() + 1; j++) {
                /* loop over local set of coarse G-vectors */
                #pragma omp parallel for schedule(static)
                for (int igloc = 0; igloc < gvec_coarse_p_->gvec().count(); igloc++) {
                    /* map from fine to coarse set of G-vectors */
                    veff_vec_[j]->f_pw_local(igloc) =
                        potential__->component(j).f_pw_local(potential__->component(j).gvec().gvec_base_mapping(igloc));
                }
                /* transform to real space */
                veff_vec_[j]->fft_transform(1);
            }

            /* change to canonical form */
            if (ctx_.num_mag_dims()) {
                #pragma omp parallel for schedule(static)
                for (int ir = 0; ir < fft_coarse_.local_slice_size(); ir++) {
                    T v0             = veff_vec_[v_local_index_t::v0]->f_rg(ir);
                    T v1             = veff_vec_[v_local_index_t::v1]->f_rg(ir);
                    veff_vec_[v_local_index_t::v0]->f_rg(ir) = v0 + v1; // v + Bz
                    veff_vec_[v_local_index_t::v1]->f_rg(ir) = v0 - v1; // v - Bz
                }
            }

            if (ctx_.num_mag_dims() == 0) {
                v0_[0] = potential__->component(0).f_0().real();
            } else {
                v0_[0] = potential__->component(0).f_0().real() + potential__->component(1).f_0().real();
                v0_[1] = potential__->component(0).f_0().real() - potential__->component(1).f_0().real();
            }
        }

        if (ctx_.print_checksum()) {
            for (int j = 0; j < ctx_.num_mag_dims() + 1; j++) {
                auto cs1 = veff_vec_[j]->checksum_pw();
                auto cs2 = veff_vec_[j]->checksum_rg();
                if (ctx_.comm().rank() == 0) {
                    utils::print_checksum("veff_pw", cs1);
                    utils::print_checksum("veff_rg", cs2);
                }
            }
        }
    }

    buf_rg_ = sddk::mdarray<std::complex<T>, 1>(fft_coarse_.local_slice_size(), ctx_.mem_pool(sddk::memory_t::host),
                                         "Local_operator::buf_rg_");
    /* move functions to GPU */
    if (fft_coarse_.processing_unit() == SPFFT_PU_GPU) {
        for (int j = 0; j < 6; j++) {
            if (veff_vec_[j]) {
                veff_vec_[j]->f_rg().allocate(ctx_.mem_pool(sddk::memory_t::device)).copy_to(sddk::memory_t::device);
            }
        }
        buf_rg_.allocate(ctx_.mem_pool(sddk::memory_t::device));
    }
}

template <typename T>
void Local_operator<T>::prepare_k(sddk::Gvec_fft const& gkvec_p__)
{
    PROFILE("sirius::Local_operator::prepare_k");

    int ngv_fft = gkvec_p__.gvec_count_fft();

    /* cache kinteic energy of plane-waves */
    if (static_cast<int>(pw_ekin_.size()) < ngv_fft) {
        pw_ekin_ = sddk::mdarray<T, 1>(ngv_fft, ctx_.mem_pool(sddk::memory_t::host), "Local_operator::pw_ekin");
    }
    #pragma omp parallel for schedule(static)
    for (int ig_loc = 0; ig_loc < ngv_fft; ig_loc++) {
        /* get G+k in Cartesian coordinates */
        auto gv          = gkvec_p__.gkvec_cart(ig_loc);
        pw_ekin_[ig_loc] = 0.5 * dot(gv, gv);
    }

    if (static_cast<int>(vphi_.size(0)) < ngv_fft) {
        vphi_ = sddk::mdarray<std::complex<T>, 1>(ngv_fft, ctx_.mem_pool(sddk::memory_t::host), "Local_operator::vphi");
    }

    if (fft_coarse_.processing_unit() == SPFFT_PU_GPU) {
        pw_ekin_.allocate(ctx_.mem_pool(sddk::memory_t::device)).copy_to(sddk::memory_t::device);
        vphi_.allocate(ctx_.mem_pool(sddk::memory_t::device));
    }
}

#ifdef SIRIUS_GPU
void mul_by_veff_real_real_gpu(int nr__, float* buf__, float* veff__)
{
    mul_by_veff_real_real_gpu_float(nr__, buf__, veff__);
}

void mul_by_veff_real_real_gpu(int nr__, double* buf__, double* veff__)
{
    mul_by_veff_real_real_gpu_double(nr__, buf__, veff__);
}

void mul_by_veff_complex_real_gpu(int nr__, std::complex<float>* buf__, float* veff__)
{
    mul_by_veff_complex_real_gpu_float(nr__, buf__, veff__);
}

void mul_by_veff_complex_real_gpu(int nr__, double_complex* buf__, double* veff__)
{
    mul_by_veff_complex_real_gpu_double(nr__, buf__, veff__);
}

void mul_by_veff_complex_complex_gpu(int nr__, std::complex<float>* buf__, float pref__, float* vx__, float* vy__)
{
    mul_by_veff_complex_complex_gpu_float(nr__, buf__, pref__, vx__, vy__);
}

void mul_by_veff_complex_complex_gpu(int nr__, double_complex* buf__, double pref__, double* vx__, double* vy__)
{
    mul_by_veff_complex_complex_gpu_double(nr__, buf__, pref__, vx__, vy__);
}
#endif

/// Multiply FFT buffer by the effective potential.
template <typename T>
static inline
std::enable_if_t<std::is_scalar<T>::value, void>
mul_by_veff(spfft_transform_type<T>& spfftk__, T* buff__,
    std::array<std::unique_ptr<Smooth_periodic_function<T>>, 6> const& veff_vec__, int idx_veff__)
{
    PROFILE("sirius::mul_by_veff");

    int nr = spfftk__.local_slice_size();

    switch (spfftk__.processing_unit()) {
        case SPFFT_PU_HOST: {
            if (idx_veff__ <= 1 || idx_veff__ >= 4) { /* up-up or dn-dn block or Theta(r) */
                switch (spfftk__.type()) {
                    case SPFFT_TRANS_R2C: {
                        #pragma omp parallel for schedule(static)
                        for (int ir = 0; ir < nr; ir++) {
                            /* multiply by V+Bz or V-Bz (in PP-PW case) or by V(r), B_z(r) or Theta(r) (in LAPW case) */
                            buff__[ir] *= veff_vec__[idx_veff__]->f_rg(ir);
                        }
                        break;
                    }
                    case SPFFT_TRANS_C2C: {
                        auto wf = reinterpret_cast<std::complex<T>*>(buff__);
                        #pragma omp parallel for schedule(static)
                        for (int ir = 0; ir < nr; ir++) {
                            /* multiply by V+Bz or V-Bz (in PP-PW case) or by V(r), B_z(r) or Theta(r) (in LAPW case) */
                            wf[ir] *= veff_vec__[idx_veff__]->f_rg(ir);
                        }
                        break;
                    }
                }
            } else { /* special case for idx_veff = 2 or idx_veff__ = 3 */
                T pref  = (idx_veff__ == 2) ? -1 : 1;
                auto wf = reinterpret_cast<std::complex<T>*>(buff__);
                #pragma omp parallel for schedule(static)
                for (int ir = 0; ir < nr; ir++) {
                    /* multiply by Bx +/- i*By */
                    wf[ir] *= std::complex<T>(veff_vec__[2]->f_rg(ir), pref * veff_vec__[3]->f_rg(ir));
                }
            }
            break;
        }
        case SPFFT_PU_GPU: {
#if defined(SIRIUS_GPU)
            if (idx_veff__ <= 1 || idx_veff__ >= 4) { /* up-up or dn-dn block or Theta(r) */
                switch (spfftk__.type()) {
                    case SPFFT_TRANS_R2C: {
                        /* multiply by V+Bz or V-Bz (in PP-PW case) or by V(r), B_z(r) or Theta(r) (in LAPW case) */
                        mul_by_veff_real_real_gpu(nr, buff__, veff_vec__[idx_veff__]->f_rg().at(sddk::memory_t::device));
                        break;
                    }
                    case SPFFT_TRANS_C2C: {
                        auto wf = reinterpret_cast<std::complex<T>*>(buff__);
                        /* multiply by V+Bz or V-Bz (in PP-PW case) or by V(r), B_z(r) or Theta(r) (in LAPW case) */
                        mul_by_veff_complex_real_gpu(nr, wf, veff_vec__[idx_veff__]->f_rg().at(sddk::memory_t::device));
                        break;
                    }
                }
            } else {
                /* multiply by Bx +/- i*By */
                T pref  = (idx_veff__ == 2) ? -1 : 1;
                auto wf = reinterpret_cast<std::complex<T>*>(buff__);
                mul_by_veff_complex_complex_gpu(nr, wf, pref, veff_vec__[2]->f_rg().at(sddk::memory_t::device),
                    veff_vec__[3]->f_rg().at(sddk::memory_t::device));
            }
            break;
#endif
        }
        break;
    }
}

/// Multiply wave-function by effective potential and store in FFT buffer.
template <typename T>
static inline
std::enable_if_t<std::is_scalar<T>::value, void>
mul_phi_by_veff(spfft_transform_type<T>& spfftk__, T const* phi_r__,
        std::array<std::unique_ptr<Smooth_periodic_function<T>>, 6> const& veff_vec__, int idx_veff__)
{
    /* assume the location of data on the current processing unit */
    auto spfft_pu = spfftk__.processing_unit();

    /* number of real-space points in the local part of FFT buffer */
    int nr = spfftk__.local_slice_size();

    /* pointer to memory where SpFFT stores real-space data */
    auto spfft_buf = spfftk__.space_domain_data(spfft_pu);

    switch (spfft_pu) {
        case SPFFT_PU_HOST: {
            if (idx_veff__ <= 1 || idx_veff__ >= 4) { /* up-up or dn-dn block or Theta(r) */
                switch (spfftk__.type()) {
                    case SPFFT_TRANS_R2C: { /* real-space buffer is real */
                        #pragma omp parallel for
                        for (int ir = 0; ir < nr; ir++) {
                            /* multiply by V+Bz or V-Bz (in PP-PW case) or by V(r), B_z(r) or Theta(r) (in LAPW case) */
                            spfft_buf[ir] = phi_r__[ir] * veff_vec__[idx_veff__]->f_rg(ir);
                        }
                        break;
                    }
                    case SPFFT_TRANS_C2C: {
                        auto spfft_buf_c = reinterpret_cast<std::complex<T>*>(spfft_buf);
                        auto phi_r_c = reinterpret_cast<std::complex<T> const*>(phi_r__);
                        #pragma omp parallel for
                        for (int ir = 0; ir < nr; ir++) {
                            /* multiply by V+Bz or V-Bz (in PP-PW case) or by V(r), B_z(r) or Theta(r) (in LAPW case) */
                            spfft_buf_c[ir] = phi_r_c[ir] * veff_vec__[idx_veff__]->f_rg(ir);
                        }
                        break;
                    }
                }
            } else { /* special case for idx_veff = 2 or idx_veff__ = 3 */
                T pref  = (idx_veff__ == 2) ? -1 : 1;
                auto spfft_buf_c = reinterpret_cast<std::complex<T>*>(spfft_buf);
                auto phi_r_c = reinterpret_cast<std::complex<T> const*>(phi_r__);
                #pragma omp parallel for
                for (int ir = 0; ir < nr; ir++) {
                    /* multiply by Bx +/- i*By */
                    spfft_buf_c[ir] = phi_r_c[ir] * std::complex<T>(veff_vec__[2]->f_rg(ir), pref * veff_vec__[3]->f_rg(ir));
                }
            }
            break;
        }
        case SPFFT_PU_GPU: {
            RTE_THROW("implement this");
//#if defined(SIRIUS_GPU)
//            if (idx_veff__ <= 1 || idx_veff__ >= 4) { /* up-up or dn-dn block or Theta(r) */
//                switch (spfftk__.type()) {
//                    case SPFFT_TRANS_R2C: {
//                        /* multiply by V+Bz or V-Bz (in PP-PW case) or by V(r), B_z(r) or Theta(r) (in LAPW case) */
//                        mul_by_veff_real_real_gpu(nr, buff__, veff_vec__[idx_veff__]->f_rg().at(sddk::memory_t::device));
//                        break;
//                    }
//                    case SPFFT_TRANS_C2C: {
//                        auto wf = reinterpret_cast<std::complex<T>*>(buff__);
//                        /* multiply by V+Bz or V-Bz (in PP-PW case) or by V(r), B_z(r) or Theta(r) (in LAPW case) */
//                        mul_by_veff_complex_real_gpu(nr, wf, veff_vec__[idx_veff__]->f_rg().at(sddk::memory_t::device));
//                        break;
//                    }
//                }
//            } else {
//                /* multiply by Bx +/- i*By */
//                T pref  = (idx_veff__ == 2) ? -1 : 1;
//                auto wf = reinterpret_cast<std::complex<T>*>(buff__);
//                mul_by_veff_complex_complex_gpu(nr, wf, pref, veff_vec__[2]->f_rg().at(sddk::memory_t::device),
//                    veff_vec__[3]->f_rg().at(sddk::memory_t::device));
//            }
//            break;
//#endif
        }
        break;
    }
}
#ifdef SIRIUS_GPU
void add_pw_ekin_gpu(int num_gvec__, float alpha__, float const* pw_ekin__, std::complex<float> const* phi__,
        std::complex<float> const* vphi__, std::complex<float>* hphi__)
{
    add_pw_ekin_gpu_float(num_gvec__, alpha__, pw_ekin__, phi__, vphi__, hphi__);
}

void add_pw_ekin_gpu(int num_gvec__, double alpha__, double const* pw_ekin__, std::complex<double> const* phi__,
        std::complex<double> const* vphi__, std::complex<double>* hphi__)
{
    add_pw_ekin_gpu_double(num_gvec__, alpha__, pw_ekin__, phi__, vphi__, hphi__);
}
#endif

template <typename T>
void
Local_operator<T>::apply_h(spfft_transform_type<T>& spfftk__, std::shared_ptr<sddk::Gvec_fft> gkvec_fft__,
    wf::spin_range spins__, wf::Wave_functions<T> const& phi__, wf::Wave_functions<T>& hphi__, wf::band_range br__)
{
    PROFILE("sirius::Local_operator::apply_h");

    if ((spfftk__.dim_x() != fft_coarse_.dim_x()) || (spfftk__.dim_y() != fft_coarse_.dim_y()) ||
        (spfftk__.dim_z() != fft_coarse_.dim_z())) {
        RTE_THROW("wrong FFT dimensions");
    }

    /* increment the counter by the number of wave-functions */
    ctx_.num_loc_op_applied(br__.size());

    /* this memory pool will be used to allocate extra storage in the host memory */
    //auto& mp = const_cast<Simulation_context&>(ctx_).mem_pool(ctx_.host_memory_t());
    /* this memory pool will be used to allocate extra storage in the device memory */
#if defined(SIRIUS_GPU)
    //sddk::memory_pool* mpd = &const_cast<Simulation_context&>(ctx_).mem_pool(sddk::memory_t::device);
#else
   // sddk::memory_pool* mpd{nullptr};
#endif

    /* local number of G-vectors for the FFT transformation */
    int ngv_fft = gkvec_fft__->gvec_count_fft();

    if (ngv_fft != spfftk__.num_local_elements()) {
        TERMINATE("wrong number of G-vectors");
    }

    std::array<wf::Wave_functions_fft_new<T>, 2> phi_fft;
    std::array<wf::Wave_functions_fft_new<T>, 2> hphi_fft;
    for (auto s = spins__.begin(); s != spins__.end(); s++) {
        phi_fft[s.get()] = wf::Wave_functions_fft_new<T>(gkvec_fft__, const_cast<wf::Wave_functions<T>&>(phi__), s,
                br__, wf::transform_layout::to);

        hphi_fft[s.get()] = wf::Wave_functions_fft_new<T>(gkvec_fft__, hphi__, s, br__, wf::transform_layout::from);
        auto hphi_mem = hphi_fft[s.get()].on_device() ? sddk::memory_t::device : sddk::memory_t::host;
        hphi_fft[s.get()].zero(hphi_mem, wf::spin_index(0), wf::band_range(0, hphi_fft[s.get()].num_wf_local()));
    }

    auto spl_num_wf = phi_fft[spins__.begin().get()].spl_num_wf();

    /* assume the location of data on the current processing unit */
    auto spfft_mem = spfftk__.processing_unit();

    /* number of real-space points in the local part of FFT buffer */
    int nr = spfftk__.local_slice_size();

    /* pointer to FFT buffer */
    auto spfft_buf = spfftk__.space_domain_data(spfft_mem);

    /* transform wave-function to real space; the result of the transformation is stored in the FFT buffer */
    auto phi_to_r = [&](wf::spin_index ispn,  wf::band_index i) {
        PROFILE("phi_to_r");
        auto phi_mem = phi_fft[ispn.get()].on_device() ? sddk::memory_t::device : sddk::memory_t::host;
        spfftk__.backward(phi_fft[ispn.get()].pw_coeffs_spfft(phi_mem, i), spfft_mem);
    };

    /* transform function to PW domain */
    auto vphi_to_G = [&]() {
        PROFILE("vphi_to_G");
        spfftk__.forward(spfft_mem, reinterpret_cast<T*>(vphi_.at(spfft_memory_t.at(spfft_mem))), SPFFT_FULL_SCALING);
    };

    /* store the resulting hphi
       spin block (ispn_block) is used as a bit mask:
        - first bit: spin component which is updated
        - second bit: add or not kinetic energy term */
    auto add_to_hphi = [&](int ispn_block, wf::band_index i) {
        PROFILE("add_to_hphi");
        /* index of spin component */
        int ispn = ispn_block & 1;
        /* add kinetic energy if this is a diagonal block */
        int ekin = (ispn_block & 2) ? 0 : 1;

        auto hphi_mem = hphi_fft[ispn].on_device() ? sddk::memory_t::device : sddk::memory_t::host;

        switch (hphi_mem) {
            case sddk::memory_t::host: {
                if (spfft_mem == SPFFT_PU_GPU) {
                    vphi_.copy_to(sddk::memory_t::host);
                }
                /* CPU case */
                if (ekin) {
                    #pragma omp parallel for
                    for (int ig = 0; ig < ngv_fft; ig++) {
                        hphi_fft[ispn].pw_coeffs(ig, i) += phi_fft[ispn].pw_coeffs(ig, i) * pw_ekin_[ig] + vphi_[ig];
                    }
                } else {
                    #pragma omp parallel for
                    for (int ig = 0; ig < ngv_fft; ig++) {
                        hphi_fft[ispn].pw_coeffs(ig, wf::band_index(i)) += vphi_[ig];
                    }
                }
                break;
            }
            case sddk::memory_t::device: {
#if defined(SIRIUS_GPU)
                T alpha = static_cast<T>(ekin);
                add_pw_ekin_gpu(ngv_fft, alpha, pw_ekin_.at(sddk::memory_t::device),
                                phi_fft[ispn].at(sddk::memory_t::device, 0, wf::spin_index(ispn), i),
                                vphi_.at(sddk::memory_t::device),
                                hphi_fft[ispn].at(sddk::memory_t::device, 0, wf::spin_index(ispn), i));
#endif
                break;
            }
            default: {
                break;
            }
        }
    };

    auto copy_phi = [&]()
    {
        switch (spfft_mem) {
            /* this is a non-collinear case, so the wave-functions and FFT buffer are complex and
               we can copy memory */
            case SPFFT_PU_HOST: {
                auto inp = reinterpret_cast<std::complex<T>*>(spfft_buf);
                std::copy(inp, inp + nr, buf_rg_.at(sddk::memory_t::host));
                break;
            }
            case SPFFT_PU_GPU: {
                acc::copy(buf_rg_.at(sddk::memory_t::device), reinterpret_cast<std::complex<T>*>(spfft_buf), nr);
                break;
            }
        }
    };

    PROFILE_START("sirius::Local_operator::apply_h|bands");
    for (int i = 0; i < spl_num_wf.local_size(); i++) {

        /* non-collinear case */
        /* 2x2 Hamiltonian in applied to spinor wave-functions
           .--------.--------.   .-----.   .------.
           |        |        |   |     |   |      |
           | H_{uu} | H_{ud} |   |phi_u|   |hphi_u|
           |        |        |   |     |   |      |
           .--------.--------. x .-----. = .------.
           |        |        |   |     |   |      |
           | H_{du} | H_{dd} |   |phi_d|   |hphi_d|
           |        |        |   |     |   |      |
           .--------.--------.   .-----.   .------.

           hphi_u = H_{uu} phi_u + H_{ud} phi_d
           hphi_d = H_{du} phi_u + H_{dd} phi_d

           The following indexing scheme will be used for spin-blocks
           .---.---.
           | 0 | 2 |
           .---.---.
           | 3 | 1 |
           .---.---.
        */
        if (spins__.size() == 2) {
            /* phi_u(G) -> phi_u(r) */
            phi_to_r(wf::spin_index(0), wf::band_index(i));
            /* save phi_u(r) in temporary buf_rg array */
            copy_phi();
            /* multiply phi_u(r) by effective potential */
            mul_by_veff<T>(spfftk__, spfft_buf, veff_vec_, v_local_index_t::v0);

            /* V_{uu}(r)phi_{u}(r) -> [V*phi]_{u}(G) */
            vphi_to_G();
            /* add kinetic energy */
            add_to_hphi(0, wf::band_index(i));
            /* multiply phi_{u} by V_{du} and copy to FFT buffer */
            switch (spfft_mem) {
                case SPFFT_PU_HOST: {
                    mul_by_veff<T>(spfftk__, reinterpret_cast<T*>(buf_rg_.at(sddk::memory_t::host)), veff_vec_, 3);
                    std::copy(buf_rg_.at(sddk::memory_t::host), buf_rg_.at(sddk::memory_t::host) + nr,
                              reinterpret_cast<std::complex<T>*>(spfft_buf));
                    break;
                }
                case SPFFT_PU_GPU: {
                    mul_by_veff<T>(spfftk__, reinterpret_cast<T*>(buf_rg_.at(sddk::memory_t::device)), veff_vec_, 3);
                    acc::copy(reinterpret_cast<std::complex<T>*>(spfft_buf), buf_rg_.at(sddk::memory_t::device), nr);
                    break;
                }
            }
            /* V_{du}(r)phi_{u}(r) -> [V*phi]_{d}(G) */
            vphi_to_G();
            /* add to hphi_{d} */
            add_to_hphi(3, wf::band_index(i));

            /* for the second spin component */

            /* phi_d(G) -> phi_d(r) */
            phi_to_r(wf::spin_index(1), wf::band_index(i));
            /* save phi_d(r) */
            copy_phi();
            /* multiply phi_d(r) by effective potential */
            mul_by_veff<T>(spfftk__, spfft_buf, veff_vec_,  v_local_index_t::v1);
           /* V_{dd}(r)phi_{d}(r) -> [V*phi]_{d}(G) */
            vphi_to_G();
            /* add kinetic energy */
            add_to_hphi(1, wf::band_index(i));
            /* multiply phi_{d} by V_{ud} and copy to FFT buffer */
            switch (spfft_mem) {
                case SPFFT_PU_HOST: {
                    mul_by_veff<T>(spfftk__, reinterpret_cast<T*>(buf_rg_.at(sddk::memory_t::host)), veff_vec_, 2);
                    std::copy(buf_rg_.at(sddk::memory_t::host), buf_rg_.at(sddk::memory_t::host) + nr,
                              reinterpret_cast<std::complex<T>*>(spfft_buf));
                    break;
                }
                case SPFFT_PU_GPU: {
                    mul_by_veff<T>(spfftk__, reinterpret_cast<T*>(buf_rg_.at(sddk::memory_t::device)), veff_vec_, 2);
                    acc::copy(reinterpret_cast<std::complex<T>*>(spfft_buf), buf_rg_.at(sddk::memory_t::device), nr);
                    break;
                }
            }
            /* V_{ud}(r)phi_{d}(r) -> [V*phi]_{u}(G) */
            vphi_to_G();
            /* add to hphi_{u} */
            add_to_hphi(2, wf::band_index(i));
        } else { /* spin-collinear or non-magnetic case */
            /* phi(G) -> phi(r) */
            phi_to_r(spins__.begin(), wf::band_index(i));
            /* multiply by effective potential */
            mul_by_veff<T>(spfftk__, spfft_buf, veff_vec_, spins__.begin().get());
            /* V(r)phi(r) -> [V*phi](G) */
            vphi_to_G();
            /* add kinetic energy */
            add_to_hphi(spins__.begin().get(), wf::band_index(i));
        }
    }
    PROFILE_STOP("sirius::Local_operator::apply_h|bands");
}

// This is full-potential case. Only C2C FFT transformation is considered here.
// TODO: document the data location on input/output
template <typename T>
void Local_operator<T>::apply_fplapw(spfft_transform_type<T>& spfftk__, std::shared_ptr<sddk::Gvec_fft> gkvec_fft__,
        wf::band_range b__, wf::Wave_functions<T>& phi__, wf::Wave_functions<T>* hphi__, wf::Wave_functions<T>* ophi__, 
        wf::Wave_functions<T>* bzphi__, wf::Wave_functions<T>* bxyphi__)
{
    PROFILE("sirius::Local_operator::apply_h_o");

    ctx_.num_loc_op_applied(b__.size());

    std::map<wf::Wave_functions<T>*, wf::Wave_functions_fft_new<T>> map_wf_fft;
    if (hphi__) {
        map_wf_fft[hphi__] = wf::Wave_functions_fft_new<T>(gkvec_fft__, *hphi__, wf::spin_index(0), b__,
                wf::transform_layout::from);
    }
    if (ophi__) {
        map_wf_fft[ophi__] = wf::Wave_functions_fft_new<T>(gkvec_fft__, *ophi__, wf::spin_index(0), b__,
                wf::transform_layout::from);
    }
    if (bzphi__) {
        map_wf_fft[bzphi__] = wf::Wave_functions_fft_new<T>(gkvec_fft__, *bzphi__, wf::spin_index(0), b__,
                wf::transform_layout::from);
    }
    if (bxyphi__) {
        map_wf_fft[bzphi__] = wf::Wave_functions_fft_new<T>(gkvec_fft__, *bxyphi__, wf::spin_index(0), b__,
                wf::transform_layout::from);
    }

    // TODO: need to pass temporaty functions or allocate from the pool
    wf::Wave_functions_fft_new<T> phi_fft(gkvec_fft__, phi__, wf::spin_index(0), b__, wf::transform_layout::to);

    //if (ctx_.processing_unit() == device_t::GPU) {
    //    phi__.pw_coeffs(0).copy_to(sddk::memory_t::host, N__, n__);
    //}
    // if (ctx_->control().print_checksum_) {
    //    auto cs = phi__.checksum_pw(N__, n__, ctx_->processing_unit());
    //    if (phi__.comm().rank() == 0) {
    //        DUMP("checksum(phi_pw): %18.10f %18.10f", cs.real(), cs.imag());
    //    }
    //}

    //auto& mp = const_cast<Simulation_context&>(ctx_).mem_pool(sddk::memory_t::host);

    auto spl_num_wf = phi_fft.spl_num_wf();

    //phi__.pw_coeffs(0).remap_forward(n__, N__, &mp);

    //if (hphi__ != nullptr) {
    //    hphi__->pw_coeffs(0).set_num_extra(n__, N__, &mp);
    //}

    //if (ophi__ != nullptr) {
    //    ophi__->pw_coeffs(0).set_num_extra(n__, N__, &mp);
    //}

    /* assume the location of data on the current processing unit */
    auto spfft_mem = spfftk__.processing_unit();

    /* number of real-space points in the local part of FFT buffer */
    int nr = spfftk__.local_slice_size();

    /* pointer to memory where SpFFT stores real-space data */
    auto spfft_buf = spfftk__.space_domain_data(spfft_mem);

    sddk::mdarray<std::complex<T>, 1> buf_pw(gkvec_fft__->gvec_count_fft(), ctx_.mem_pool(sddk::memory_t::host));

    auto phi_r = buf_rg_.at(sddk::memory_t::host);

    for (int j = 0; j < spl_num_wf.local_size(); j++) {
        /* phi(G) -> phi(r) */
        spfftk__.backward(phi_fft.pw_coeffs_spfft(sddk::memory_t::host, wf::band_index(j)), spfft_mem);

        /* we are going to apply parts of the Hamiltonian to the wave-function; save wave-function on the
         * real space grid first */

        /* save phi(r); real-space data is complex */
        auto inp = reinterpret_cast<std::complex<T>*>(spfft_buf);
        switch (spfft_mem) {
            case SPFFT_PU_HOST: {
                std::copy(inp, inp + nr, phi_r);
                break;
            }
            case SPFFT_PU_GPU: {
                acc::copy(phi_r, inp, nr);
                break;
            }
        }

        if (ophi__) {
            /* multiply phi(r) by step function */
            mul_phi_by_veff(spfftk__, reinterpret_cast<T const*>(phi_r), veff_vec_, v_local_index_t::theta);

            /* phi(r) * Theta(r) -> ophi(G) */
            spfftk__.forward(spfft_mem, map_wf_fft[ophi__].pw_coeffs_spfft(sddk::memory_t::host, wf::band_index(j)),
                    SPFFT_FULL_SCALING);
        }

        if (bzphi__) {
            mul_phi_by_veff(spfftk__, reinterpret_cast<T const*>(phi_r), veff_vec_, v_local_index_t::v1);
            /* phi(r) * Bz(r) -> bzphi(G) */
            spfftk__.forward(spfft_mem, map_wf_fft[bzphi__].pw_coeffs_spfft(sddk::memory_t::host, wf::band_index(j)),
                    SPFFT_FULL_SCALING);
        }

        if (bxyphi__) {
            mul_phi_by_veff(spfftk__, reinterpret_cast<T const*>(phi_r), veff_vec_, 2);
            /* phi(r) * (Bx(r) - iBy(r)) -> bxyphi(G) */
            spfftk__.forward(spfft_mem, map_wf_fft[bxyphi__].pw_coeffs_spfft(sddk::memory_t::host, wf::band_index(j)),
                    SPFFT_FULL_SCALING);
        }

        if (hphi__) {
            mul_phi_by_veff(spfftk__, reinterpret_cast<T const*>(phi_r), veff_vec_, v_local_index_t::v0);
            /* phi(r) * Theta(r) * V(r) -> hphi(G) */
            spfftk__.forward(spfft_mem, map_wf_fft[hphi__].pw_coeffs_spfft(sddk::memory_t::host, wf::band_index(j)),
                    SPFFT_FULL_SCALING);

            /* add kinetic energy */
            for (int x : {0, 1, 2}) {
                #pragma omp parallel for
                for (int igloc = 0; igloc < gkvec_fft__->gvec_count_fft(); igloc++) {
                    auto gvc = gkvec_fft__->gkvec_cart(igloc);
                    /* \hat P phi = phi(G+k) * (G+k), \hat P is momentum operator */
                    buf_pw[igloc] = phi_fft.pw_coeffs(igloc, wf::band_index(j)) * static_cast<T>(gvc[x]);
                }
                /* transform Cartesian component of wave-function gradient to real space */
                spfftk__.backward(reinterpret_cast<T const*>(&buf_pw[0]), spfft_mem);
                /* multiply by real-space function */
                switch (ctx_.valence_relativity()) {
                    case relativity_t::iora:
                    case relativity_t::zora: {
                        /* multiply be inverse relative mass */
                        mul_phi_by_veff(spfftk__, spfft_buf, veff_vec_, v_local_index_t::rm_inv);
                        break;
                    }
                    case relativity_t::none: {
                        /* multiply be step function */
                        mul_phi_by_veff(spfftk__, spfft_buf, veff_vec_, v_local_index_t::theta);
                        break;
                    }
                    default: {
                        break;
                    }
                }
                /* transform back to PW domain */
                spfftk__.forward(spfft_mem, reinterpret_cast<T*>(&buf_pw[0]), SPFFT_FULL_SCALING);
                #pragma omp parallel for
                for (int igloc = 0; igloc < gkvec_fft__->gvec_count_fft(); igloc++) {
                    auto gvc = gkvec_fft__->gkvec_cart(igloc);
                    map_wf_fft[hphi__].pw_coeffs(igloc, wf::band_index(j)) += buf_pw[igloc] * static_cast<T>(0.5 * gvc[x]);
                }
            } // x
        }
    }
    //if (hphi__) {
    //    wf::transform_from_fft_layout(map_wf_fft[hphi__], *hphi__, wf::spin_index(0), b__);
    //}
    //if (ophi__) {
    //    wf::transform_from_fft_layout(map_wf_fft[ophi__], *ophi__, wf::spin_index(0), b__);
    //}
    //if (bzphi__) {
    //    wf::transform_from_fft_layout(map_wf_fft[bzphi__], *bzphi__, wf::spin_index(0), b__);
    //}
    //if (bxyphi__) {
    //    wf::transform_from_fft_layout(map_wf_fft[bxyphi__], *bxyphi__, wf::spin_index(0), b__);
    //}

    //if (hphi__ != nullptr) {
    //    hphi__->pw_coeffs(0).remap_backward(n__, N__);
    //}
    //if (ophi__ != nullptr) {
    //    ophi__->pw_coeffs(0).remap_backward(n__, N__);
    //}

    //if (ctx_.processing_unit() == device_t::GPU) {
    //    if (hphi__ != nullptr) {
    //        hphi__->pw_coeffs(0).copy_to(sddk::memory_t::device, N__, n__);
    //    }
    //    if (ophi__ != nullptr) {
    //        ophi__->pw_coeffs(0).copy_to(sddk::memory_t::device, N__, n__);
    //    }
    //}
}

// instantiate for supported precision
template class Local_operator<double>;
#ifdef USE_FP32
template class Local_operator<float>;
#endif
} // namespace sirius
