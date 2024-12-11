/** \file generate_w90_coeffs.hpp
 *
 *  \brief Interface to W90 library.
 */
////////#ifdef SIRIUS_WANNIER90

#include <limits>
#include "dft/smearing.hpp"
#include "k_point/k_point.hpp"
#include "k_point/k_point_set.hpp"
#include "symmetry/get_irreducible_reciprocal_mesh.hpp"
#include "hamiltonian/non_local_operator.hpp"
#include "core/la/inverse_sqrt.hpp"
#include "generate_w90_coeffs.hpp"

namespace sirius {

int
find_line( const std::string& string__, const std::vector<std::string>& file_content__ )
{
    auto iterator = std::find_if(file_content__.begin(), file_content__.end(),
                         [&, string__]( const std::string& iter_file ) { return (string__ == iter_file); });
    return iterator - file_content__.begin();
}

void
read_nnkp(int& num_wann, int& nntot, mdarray<int, 2>& nnlist, mdarray<int32_t, 3>& nncell,
          mdarray<int32_t, 1>& exclude_bands, mdarray<double,2>& kp)
{
    /* read file as string set and split it in lines */
    auto fpath = "sirius.nnkp";
    if (!std::filesystem::exists(fpath)) {
        std::stringstream s;
        s << "Error in routine [read_nnkp]: Could not find file " << fpath;
        RTE_THROW(s);
    }
    std::ifstream readNNKP("sirius.nnkp");
    std::string line;
    std::vector<std::string> file_content;
    while (std::getline(readNNKP, line)) {
        file_content.push_back(line);
    }
    int iline;

    /* read kpoints contained in nnkp, to be matched with our K_point_set */
    iline = find_line("begin kpoints", file_content) + 1;
    auto nk_wannier = std::atoi( (*(file_content.begin() + iline) ).c_str() );

    kp = mdarray<double, 2>( { 3, nk_wannier } );
    for( int ik = 0; ik < nk_wannier; ++ik ) {
        ++iline;
        std::stringstream split_line;
        split_line << *( file_content.begin() + iline );
        split_line >> kp(0, ik) >> kp(1, ik) >> kp(2, ik);
    }

    /* read num_wann */
    iline = find_line("begin projections", file_content) + 1;
    num_wann = std::atoi((*(file_content.begin() + iline)).c_str());
    
    /* read nnlist and nncell */
    iline = find_line("begin nnkpts", file_content) + 1;
    nntot = std::atoi((*(file_content.begin() + iline)).c_str());

    iline++;
    int aux_int;
    for (int ik = 0; ik < nk_wannier; ik++) {
        for (int ib = 0; ib < nntot; ib++) {
            std::stringstream split_line;
            split_line << *(file_content.begin() + iline);
            split_line >> aux_int;
            assert(aux_int == ik + 1);
            split_line >> nnlist(ik, ib);
            split_line >> nncell(0, ik, ib) >> nncell(1, ik, ib) >> nncell(2, ik, ib);
            iline++;
        }
    }

    /* read exclude_bands */
    std::fill(exclude_bands.begin(), exclude_bands.end(), -1);
    iline = find_line("begin exclude_bands", file_content) + 1;
    int num_bands_excluded = std::atoi((*(file_content.begin() + iline)).c_str());
    for ( int ibnd = 0; ibnd < num_bands_excluded; ibnd++ ) {
        iline++;
        exclude_bands( ibnd ) = std::atoi((*(file_content.begin() + iline)).c_str());
    }
}

/*
 * This function creates a file with extension ".amn" that can eventually be read by wannier90
 * to set the matrix Amn (not needed if we want to use the library)
 */
void
write_Amn(mdarray<std::complex<double>, 3> const& Amn, int const& num_kpoints, int const& num_bands, int const& num_wann)
{
    std::ofstream writeAmn;
    writeAmn.open("sirius.amn");
    std::string line;
    writeAmn << "#produced in sirius" << std::endl;
    writeAmn << std::setw(12) << num_bands;
    writeAmn << std::setw(12) << num_kpoints;
    writeAmn << std::setw(12) << num_wann;
    writeAmn << std::endl;

    for (int ik = 0; ik < num_kpoints; ik++) {
        for (int n = 0; n < num_wann; n++) {
            for (int m = 0; m < num_bands; m++) {
                writeAmn << std::fixed << std::setw(10) << m + 1;
                writeAmn << std::fixed << std::setw(10) << n + 1;
                writeAmn << std::fixed << std::setw(10) << ik + 1;
                writeAmn << std::fixed << std::setprecision(12) << std::setw(18) << Amn(m, n, ik).real();
                writeAmn << std::fixed << std::setprecision(12) << std::setw(18) << Amn(m, n, ik).imag();
                // writeAmn << std::fixed << std::setprecision(12) << std::setw(18) << abs(Amn(m, n, ik));
                writeAmn << std::endl;
            }
        }
    }
}

/*
 * This function creates a file with extension ".mmn" that can eventually be read by wannier90
 * to set the matrix Mmn (not needed if we want to use the library)
 */
void
write_Mmn(mdarray<std::complex<double>, 4> const& M, mdarray<int, 2> const& nnlist, mdarray<int32_t, 3> const& nncell,
          int const& num_kpoints, int const& num_neighbors, int const& num_bands)
{
    std::ofstream writeMmn;
    writeMmn.open("sirius.mmn");
    writeMmn << "#produced in sirius" << std::endl;
    writeMmn << std::setw(12) << num_bands;
    writeMmn << std::setw(12) << num_kpoints;
    writeMmn << std::setw(12) << num_neighbors;
    writeMmn << std::endl;
    for (int ik = 0; ik < num_kpoints; ik++) {
        for (int ib = 0; ib < num_neighbors; ib++) {
            writeMmn << std::setw(10) << ik + 1;
            writeMmn << std::setw(10) << nnlist(ik, ib);
            writeMmn << std::setw(10) << nncell(0, ik, ib);
            writeMmn << std::setw(10) << nncell(1, ik, ib);
            writeMmn << std::setw(10) << nncell(2, ik, ib);
            writeMmn << std::endl;
            for (int n = 0; n < num_bands; n++) {
                for (int m = 0; m < num_bands; m++) {
                    writeMmn << std::fixed << std::setprecision(12) << std::setw(18) << M(m, n, ib, ik).real();
                    writeMmn << std::fixed << std::setprecision(12) << std::setw(18) << M(m, n, ib, ik).imag();
                    // writeMmn << std::fixed << std::setprecision(12) << std::setw(18) << abs(M(m, n, ib, ik));
                    writeMmn << std::endl;
                }
            }
        }
    }
    writeMmn.close();
}

/*
 * This function creates a file with extension ".eig" that can eventually be read by wannier90
 * to pass the energy eigenvalues if we need a window (not needed if we want to use the library)
 */
void
write_eig(mdarray<double, 2> const& eigval, int const& num_bands, int const& num_kpoints)
{
    std::ofstream writeEig;
    writeEig.open("sirius.eig");
    for (int ik = 0; ik < num_kpoints; ik++) {
        for (int iband = 0; iband < num_bands; iband++) {
            writeEig << std::setw(10) << iband + 1;
            writeEig << std::setw(10) << ik + 1;
            writeEig << std::fixed << std::setprecision(12) << std::setw(18) << eigval(iband, ik);
            writeEig << std::endl;
        }
    }
    writeEig.close();
}

/*
 * This function generates the Full Brillouin zone starting from the Irreducible wedge.
 * The equation to satisfy is:
 * \f[
 *     {\bf k}_{wann} + G = R.{\bf k}_{sirius}
 * \f]
 */
void
match_kpoints(K_point_set& kset_sirius, K_point_set& kset_wannier, std::vector<k_info>& k_temp)
{
    PROFILE_START("sirius::K_point_set::generate_w90_coeffs::match_kpoints");
    // Apply symmetry to all points of the IBZ. Save indices of ibz, fbz, sym
    k_temp.resize( kset_wannier.num_kpoints() );
    std::pair<r3::vector<double>, r3::vector<int>> Rk_reduced;

    for (int ikwann = 0; ikwann < kset_wannier.num_kpoints(); ikwann++) {
        auto kp_wann = r3::reduce_coordinates( kset_wannier.get<double>(ikwann)->vk() );
        bool found = false;
        for (int isym = 0; isym < kset_sirius.ctx().unit_cell().symmetry().size(); isym++) {
            for (int iksirius = 0; iksirius < kset_sirius.num_kpoints(); iksirius++) {
                auto& R = kset_sirius.ctx()
                                  .unit_cell()
                                  .symmetry()[isym]
                                  .spg_op.R; // point symmetry rotation in crystal coordinates

                auto Rk         = r3::dot(kset_sirius.get<double>(iksirius)->vk(), R);
                Rk_reduced      = r3::reduce_coordinates(Rk);
                found      = (( kp_wann.first - Rk_reduced.first ).length() < 1.e-08);
                if( found ) {
                    k_temp[ikwann].wannier    = kset_sirius.get<double>(ikwann)->vk();
                    k_temp[ikwann].isirius    = iksirius;
                    k_temp[ikwann].G          = kp_wann.second - Rk_reduced.second;
                    k_temp[ikwann].R          = kset_wannier.ctx().unit_cell().symmetry()[isym].spg_op.R;
                    k_temp[ikwann].invR       = kset_wannier.ctx().unit_cell().symmetry()[isym].spg_op.invR;
                    k_temp[ikwann].t          = kset_wannier.ctx().unit_cell().symmetry()[isym].spg_op.t;

                    std::cout << " new_kpt.ibz    =  " << k_temp[ikwann].wannier    << std::endl;
                    std::cout << " new_kpt.ik_ibz =  " << k_temp[ikwann].isirius << std::endl;
                    std::cout << " new_kpt.G      =  " << k_temp[ikwann].G      << std::endl;
                    std::cout << " new_kpt.R      =  " << k_temp[ikwann].R      << std::endl;
                    std::cout << " new_kpt.invR   =  " << k_temp[ikwann].invR   << std::endl;
                    std::cout << " new_kpt.t      =  " << k_temp[ikwann].t      << std::endl;
                    break;
                }
            } //end iksirius
            if ( found ) {
                break;
            }
        } // end isym
        if ( !found ) {
            std::stringstream s;
            s << "Error in match_kpoints. Could not find " << kset_wannier.get<double>(ikwann)->vk();
            s << " in sirius k-mesh\n";
            RTE_THROW(s);
        }
    } // end ikwann

    PROFILE_STOP("sirius::K_point_set::generate_w90_coeffs::unfold_fbz");
}

/*
 * This function generates the Full Brillouin zone starting from another set linked by symmetry.
 * The equation to satisfy is:
 * \f[
 *     \psi_{n, R \bf k}({\bf G}) = e^{-i {\bf \tau}\cdot ({R\bf k}+{\bf G})}\psi_{n, \bf k} (R^{-1} {\bf G})
 * \f]
 */
void
rotate_wavefunctions(K_point_set& kset_sirius, K_point_set& kset_wannier, std::vector<k_info> const& k_temp,
                     int const& num_bands, std::vector<int> const& band_index_tot)
{
    PROFILE_START("sirius::K_point_set::generate_w90_coeffs::unfold_wfs");
    int num_bands_tot            = kset_sirius.ctx().num_bands();
    std::complex<double> imtwopi = std::complex<double>(0., twopi);
    std::complex<double> exp1, exp2;
    //srand(time(NULL));
    for (int ik = 0; ik < kset_wannier.num_kpoints(); ik++) {
        int src_rank  = kset_sirius.spl_num_kpoints().location(typename kp_index_t::global(k_temp[ik].isirius)).ib;
        int dest_rank = kset_wannier.spl_num_kpoints().location(typename kp_index_t::global(ik)).ib;

        // send gvec
        auto gvec_IBZ = std::make_shared<fft::Gvec>(
                static_cast<fft::Gvec>(kset_sirius.get_gkvec(typename kp_index_t::global(k_temp[ik].isirius), dest_rank)));

        // send wf
        auto wf_IBZ = mdarray<std::complex<double>, 2>({gvec_IBZ->num_gvec(), num_bands_tot});
        int tag     = src_rank + kset_wannier.num_kpoints() * dest_rank;
        mpi::Request req;
        if (kset_wannier.ctx().comm_k().rank() == src_rank) {
            req = kset_wannier.ctx().comm_k().isend(
                    kset_sirius.get<double>(k_temp[ik].isirius)
                            ->spinor_wave_functions()
                            .at(memory_t::host, 0, wf::spin_index(0), wf::band_index(0)),
                    kset_sirius.get<double>(k_temp[ik].isirius)->gkvec().num_gvec() * num_bands_tot, dest_rank, tag);
        }
        if (kset_wannier.ctx().comm_k().rank() == dest_rank) {
            kset_wannier.ctx().comm_k().recv(&wf_IBZ(0, 0), gvec_IBZ->num_gvec() * num_bands_tot, src_rank, tag);
        }
        // rotate wf
        if (kset_wannier.ctx().comm_k().rank() == dest_rank) {
            // kset_wannier.get<double>(ik)->spinor_wave_functions_ = std::make_unique<wf::Wave_functions<double>>(
            //     kset_wannier.get<double>(ik)->gkvec_, wf::num_mag_dims(0), wf::num_bands(num_bands_tot),
            //     kset_wannier.ctx().host_memory_t());

            // kset_wannier.get<double>(ik)->spinor_wave_functions_->zero(memory_t::host);

            std::complex<double> exp1 = exp(-imtwopi * r3::dot(kset_wannier.get<double>(ik)->vk(), k_temp[ik].t));
            r3::vector<int> invRG;
            for (int ig = 0; ig < kset_wannier.get<double>(ik)->gkvec().num_gvec(); ig++) {
                auto G_shifted = k_temp[ik].G + kset_wannier.get<double>(ik)->gkvec().gvec(gvec_index_t::local(ig));
                invRG   = r3::dot(G_shifted, k_temp[ik].invR);
                exp2    = exp(-imtwopi *r3::dot(invRG, k_temp[ik].t));
                int ig_ = gvec_IBZ->index_by_gvec(invRG);
                assert((ig_ != -1));

                for (int iband = 0; iband < num_bands; iband++) {
                    kset_wannier.get<double>(ik)->spinor_wave_functions().pw_coeffs(ig, wf::spin_index(0),
                                                                                wf::band_index(iband)) =
                            exp1 * exp2 * wf_IBZ(ig_, band_index_tot[iband]); //+
                            //std::complex<double>(rand() % 1000, rand() % 1000) *
                            //        1.e-08; // needed to not get stuck on local
                            //                // minima. not working with 1.e-09
                }
            }
            for (int iband = 0; iband < num_bands; iband++) {
                kset_wannier.get<double>(ik)->band_energy(
                        iband, 0, kset_sirius.get<double>(k_temp[ik].isirius)->band_energy(band_index_tot[iband], 0));
            }
        }
        if (src_rank == kset_wannier.ctx().comm_k().rank()) {
            req.wait();
        }
    } // end ik loop
    PROFILE_STOP("sirius::K_point_set::generate_w90_coeffs::unfold_wfs");
}

/*
 * This function calculates the projection of the Bloch functions over an initial guess for the Wannier functions.
 * The matrix A has matrix elements:
 * \f[
 *                    A_{mn}({\bf k})   = \langle u_{m \bf k}|\hat{S}|w_{n \bf k}\rangle
 * \f]
 * where u is the periodic part of the Bloch function and w is the initial guess.
 * Here we set as initial guesses the atomic orbitals of the pseudopotential.
 */
void
calculate_Amn(K_point_set& kset_wannier, int const& num_bands, int const& num_wann, mdarray<std::complex<double>, 3>& A)
{

    A.zero();
    la::dmatrix<std::complex<double>> Ak(num_bands, num_wann); // matrix at the actual k point

    std::vector<int> atoms(kset_wannier.ctx().unit_cell().num_atoms());
    std::iota(atoms.begin(), atoms.end(), 0); // we need to understand which orbitals to pick up, I am using every here
    // int num_atomic_wf = kset_wannier.ctx().unit_cell().num_ps_atomic_wf().first;

    std::unique_ptr<wf::Wave_functions<double>> Swf_k;
    // mdarray<std::complex<double>, 3> psidotpsi(num_bands, num_bands, kp.size(1)); // sirius2wannier
    // mdarray<std::complex<double>, 3> atdotat(num_wann, num_wann, kp.size(1));     // sirius2wannier
    // psidotpsi.zero();
    // atdotat.zero();
    std::cout << "Calculating Amn...\n";
    auto mem = kset_wannier.ctx().processing_unit_memory_t();
    if( is_device_memory(mem) ) {
        Ak.allocate(mem);
    }

    PROFILE_START("sirius::K_point_set::generate_w90_coeffs::calculate_Amn");

    for (auto it : kset_wannier.spl_num_kpoints()) {
        int ik = it.i;

        // calculate atomic orbitals + orthogonalization
        auto q_op =
                (kset_wannier.ctx().unit_cell().augment()) ? std::make_unique<Q_operator<double>>(kset_wannier.ctx()) : nullptr;
        // kset_wannier.kpoints_[ik]->beta_projectors().prepare();
        //Swf_k = std::make_unique<wf::Wave_functions<double>>(kset_wannier.get<double>(ik)->gkvec_sptr(),
        //                                                     wf::num_mag_dims(0), wf::num_bands(num_bands),
        //                                                     kset_wannier.ctx().host_memory_t());

        auto bp_gen    = kset_wannier.get<double>(ik)->beta_projectors().make_generator();
        auto bp_coeffs = bp_gen.prepare();

        //apply_S_operator<double, std::complex<double>>(memory_t::host, wf::spin_range(0), wf::band_range(0, num_bands), bp_gen,
        //                                               bp_coeffs, (kset_wannier.get<double>(ik)->spinor_wave_functions()),
        //                                               q_op.get(), *Swf_k);

        kset_wannier.get<double>(ik)->generate_atomic_wave_functions(
                atoms, [&](int iat) { return &kset_wannier.ctx().unit_cell().atom_type(iat).indexb_wfs(); },
                *kset_wannier.ctx().ri().ps_atomic_wf_, kset_wannier.get<double>(ik)->atomic_wave_functions());

         /*  Pick up only needed atomic functions, with their proper linear combinations
        //define index in atomic_wave_functions for atom iat
        std::vector<int> offset(kset_wannier.ctx().unit_cell().num_atoms());
        offset[0]=0;
        for(int i=1; i<kset_wannier.ctx().unit_cell().num_atoms(); i++){
            offset[i] = offset[i-1] + kset_wannier.ctx_.unit_cell().atom_type(i-1).indexb_wfs()->size();
        }
       //reconstruct map i-th wann func -> atom, l, m
        std::vector<std::array<int,3>> atoms_info(num_wann);

        auto needed_atomic_wf = std::make_unique<wf::Wave_functions<double>>(
                               kset_wannier.kpoints_[ik]->gkvec_, wf::num_mag_dims(0), wf::num_bands(num_wann),
       ctx_.host_memory_t());

        for(int iw=0; iw<num_wann; iw++)
        {
            int iat__=-1;
            for(int iat=0; iat<kset_wannier.ctx().unit_cell().num_atoms(); iat++){
                //calculate norm of center_w - atomic_position to decide which atom is the correct one
                auto& frac = this->unit_cell().atom(iat).position();
                r3::vector<double> diff = {center_w(0,iw)-frac[0], center_w(1,iw)-frac[1], center_w(2,iw)-frac[2] }
                if(diff.length() < 1.e-08){
                    iat__ = iat;
                    break;
                }
            }
            if(iat__==-1){
                std::cout <<"\n\n\nWARNING!! Could not find center_w: " << center_w(0,iw) << "  " << center_w(1,iw);
                std::cout <<"  " << center_w(2,iw) << std::endl << std::endl;
            }

            atoms_info[iw][0] = offset[iat__];
            atoms_info[iw][1] = proj_l(iw);
            atoms_info[iw][2] = proj_m(iw);
        }//end definition of atoms_info
*/
        // TODO: what is going on here?
        // it this code is taken from generate_hubbard_orbitals() then it should be reused
        // no code repetition is allowed
        // also, why do we need orthogonalized atomic orbitals for the initial guess?

        // ORTHOGONALIZING -CHECK HUBBARD FUNCTION

      apply_S_operator<double, std::complex<double>>(memory_t::host, wf::spin_range(0), wf::band_range(0, num_wann), bp_gen,
                                                       bp_coeffs, kset_wannier.get<double>(ik)->atomic_wave_functions(),
                                                       q_op.get(), kset_wannier.get<double>(ik)->atomic_wave_functions_S());

        int BS = kset_wannier.ctx().cyclic_block_size();
        la::dmatrix<std::complex<double>> ovlp(num_wann, num_wann, kset_wannier.ctx().blacs_grid(), BS, BS);
        wf::inner(kset_wannier.ctx().spla_context(), memory_t::host, wf::spin_range(0),
                  kset_wannier.get<double>(ik)->atomic_wave_functions(), wf::band_range(0, num_wann),
                  kset_wannier.get<double>(ik)->atomic_wave_functions_S(), wf::band_range(0, num_wann), ovlp, 0, 0);

        auto B = std::get<0>(inverse_sqrt(ovlp, num_wann));
        wf::transform(kset_wannier.ctx().spla_context(), memory_t::host, *B, 0, 0, 1.0,
                      kset_wannier.get<double>(ik)->atomic_wave_functions(), wf::spin_index(0), wf::band_range(0, num_wann),
                      0.0, kset_wannier.get<double>(ik)->atomic_wave_functions_S(), wf::spin_index(0),
                      wf::band_range(0, num_wann));
        wf::copy(memory_t::host, kset_wannier.get<double>(ik)->atomic_wave_functions_S(), wf::spin_index(0),
                 wf::band_range(0, num_wann), kset_wannier.get<double>(ik)->atomic_wave_functions(), wf::spin_index(0),
                 wf::band_range(0, num_wann));
        // END of the orthogonalization.
        apply_S_operator<double, std::complex<double>>(memory_t::host, wf::spin_range(0), wf::band_range(0, num_wann), bp_gen,
                                                       bp_coeffs, kset_wannier.get<double>(ik)->atomic_wave_functions(),
                                                       q_op.get(), kset_wannier.get<double>(ik)->atomic_wave_functions_S());

        //copy to device
        if( is_device_memory(mem) ) {
            kset_wannier.get<double>(ik)->spinor_wave_functions().allocate(memory_t::device);
            kset_wannier.get<double>(ik)->spinor_wave_functions().copy_to(memory_t::device);
            kset_wannier.get<double>(ik)->atomic_wave_functions_S().allocate(memory_t::device);
            kset_wannier.get<double>(ik)->atomic_wave_functions_S().copy_to(memory_t::device);
        }

        wf::inner(kset_wannier.ctx().spla_context(), mem, wf::spin_range(0),
                  kset_wannier.get<double>(ik)->spinor_wave_functions(), wf::band_range(0, num_bands),
                  kset_wannier.get<double>(ik)->atomic_wave_functions_S(), wf::band_range(0, num_wann), Ak, 0, 0);
        // already in the correct way, we just copy in the bigger array. (alternative:: create dmatrix with an index
        // as multiindex to avoid copies) note!! we need +1 to copy the last element
        if( is_device_memory(mem) ) {
            kset_wannier.get<double>(ik)->spinor_wave_functions().deallocate(memory_t::device);
            kset_wannier.get<double>(ik)->atomic_wave_functions_S().deallocate(memory_t::device);
            Ak.copy_to(memory_t::host);
        }

        std::copy(Ak.begin(), Ak.end(), A.at(memory_t::host, 0, 0, ik));

        std::cout << "Calculated Amn in rank " << kset_wannier.ctx().comm().rank() << " ik: " << ik << std::endl;
    } // end ik loop for Amn

    for (int ik = 0; ik < kset_wannier.num_kpoints(); ik++) {
        int local_rank = kset_wannier.spl_num_kpoints().location(typename kp_index_t::global(ik)).ib;
        kset_wannier.ctx().comm_k().bcast(A.at(memory_t::host, 0, 0, ik), num_bands * num_wann, local_rank);
    }
    PROFILE_STOP("sirius::K_point_set::generate_w90_coeffs::calculate_Amn");
}

/*
 * This function uses MPI to send the wavefunction at k+b to the node that holds k.
 * All the wavefunctions and G vectors will be hold in vector structures, so that when calculating M each node is
 * independent, as the information has already been passed.
 */
void
send_receive_kpb(std::vector<std::shared_ptr<fft::Gvec>>& gvec_kpb,
                 std::vector<mdarray<std::complex<double>, 2>>& wf_kpb, K_point_set& kset_wannier,
                 std::vector<int>& ikpb_index, int const& nntot, mdarray<int, 2> const& nnlist, int const& num_bands)
{
    PROFILE_START("sirius::K_point_set::generate_w90_coeffs::send_k+b");
    int index = -1; // to keep track of the index to use
    bool found;

    mpi::Request req;
    for (int ik = 0; ik < kset_wannier.num_kpoints(); ik++) {
        for (int ib = 0; ib < nntot; ib++) {
            int ikpb      = nnlist(ik, ib) - 1;
            int src_rank  = kset_wannier.spl_num_kpoints().location(typename kp_index_t::global(ikpb)).ib;
            int dest_rank = kset_wannier.spl_num_kpoints().location(typename kp_index_t::global(ik)).ib;

            int tag = src_rank + kset_wannier.num_kpoints() * kset_wannier.num_kpoints() * dest_rank;
            if (kset_wannier.ctx().comm_k().rank() == dest_rank) {
                found = ikpb_index[ikpb] != -1; // std::find(ikpb2ik_.begin(), ikpb2ik_.end(), ikpb) != ikpb2ik_.end();
                                                // //false if ikpb is not in ikpb2ik_
                req = kset_wannier.ctx().comm_k().isend(&found, 1, src_rank, tag);
            }
            if (kset_wannier.ctx().comm_k().rank() == src_rank) {
                kset_wannier.ctx().comm_k().recv(&found, 1, dest_rank, tag);
            }
            if (kset_wannier.ctx().comm_k().rank() == dest_rank) {
                req.wait();
            }

            if (found) {
                continue;
            }

            tag = src_rank + kset_wannier.num_kpoints() * dest_rank;

            auto temp = std::make_shared<fft::Gvec>(
                    static_cast<fft::Gvec>(kset_wannier.get_gkvec(typename kp_index_t::global(ikpb), dest_rank)));

            if (kset_wannier.ctx().comm_k().rank() == src_rank) {
                req = kset_wannier.ctx().comm_k().isend(kset_wannier.get<double>(ikpb)->spinor_wave_functions().at(
                                                            memory_t::host, 0, wf::spin_index(0), wf::band_index(0)),
                                                    temp->num_gvec() * num_bands, dest_rank, tag);
            }
            if (kset_wannier.ctx().comm_k().rank() == dest_rank) {
                index++;
                gvec_kpb.push_back(temp);
                wf_kpb.push_back(mdarray<std::complex<double>, 2>({gvec_kpb[index]->num_gvec(), num_bands}));
                kset_wannier.ctx().comm_k().recv(&wf_kpb[index](0, 0), gvec_kpb[index]->num_gvec() * num_bands, src_rank,
                                             tag);
                ikpb_index[ikpb] = index;
            }
            if (kset_wannier.ctx().comm_k().rank() == src_rank) {
                req.wait();
            }
        } // end ib
    }     // end ik
    PROFILE_STOP("sirius::K_point_set::generate_w90_coeffs::send_k+b");
}

/*
 * This function calculates the projection of the periodic part of the Bloch functions at k over the periodic part of
 * the Bloch function at k+b. The matrix M has matrix elements: \f[ M_{mn}({\bf k},{\bf b})   = \langle u_{m, \bf
 * k}|\hat{S}|u_{n, \bf k+b}\rangle \f] where u is the periodic part of the Bloch function. The set of neighbors k+b for
 * each k is calculated with wannier_setup.
 */
void
calculate_Mmn(mdarray<std::complex<double>, 4>& M, K_point_set& kset_wannier, int const& num_bands,
              std::vector<std::shared_ptr<fft::Gvec>> const& gvec_kpb,
              std::vector<mdarray<std::complex<double>, 2>> const& wf_kpb, std::vector<int> const& ikpb_index,
              int const& nntot, mdarray<int, 2> const& nnlist, mdarray<int, 3> const& nncell)
{
    PROFILE("sirius::K_point_set::generate_w90_coeffs::calculate_Mmn");
    la::dmatrix<std::complex<double>> Mbk(num_bands, num_bands);
    Mbk.zero();
    auto mem = kset_wannier.ctx().processing_unit_memory_t();

    if( mem == memory_t::device ) {
        Mbk.allocate(mem);
    }

    for (auto it : kset_wannier.spl_num_kpoints()) {
        int ik = it.i;
        std::cout << "Calculating Mmn. ik = " << ik << std::endl;
        auto q_op      = (kset_wannier.unit_cell().augment())
                                 ? std::make_unique<Q_operator<double>>(kset_wannier.get<double>(ik)->ctx())
                                 : nullptr;
        auto bp_gen    = kset_wannier.get<double>(ik)->beta_projectors().make_generator();
        auto bp_coeffs = bp_gen.prepare();
        auto Swf_k     = std::make_unique<wf::Wave_functions<double>>(kset_wannier.get<double>(ik)->gkvec_sptr(),
                                                                  wf::num_mag_dims(0), wf::num_bands(num_bands),
                                                                  kset_wannier.ctx().host_memory_t());
        apply_S_operator<double, std::complex<double>>(memory_t::host, wf::spin_range(0), wf::band_range(0, num_bands), bp_gen,
                                                       bp_coeffs, (kset_wannier.get<double>(ik)->spinor_wave_functions()),
                                                       q_op.get(), *Swf_k);

        for (int ib = 0; ib < nntot; ib++) {
            int ikpb        = nnlist(ik, ib) - 1;
            auto index_ikpb = ikpb_index[ikpb];
            assert((index_ikpb != -1));

            std::unique_ptr<wf::Wave_functions<double>> aux_psi_kpb = std::make_unique<wf::Wave_functions<double>>(
                    kset_wannier.get<double>(ik)->gkvec_sptr(), wf::num_mag_dims(0), wf::num_bands(num_bands),
                    kset_wannier.ctx().host_memory_t());
            aux_psi_kpb->zero(memory_t::host);
            r3::vector<int> G;
            for (int ig = 0; ig < kset_wannier.get<double>(ik)->gkvec().num_gvec(); ig++) {
                // compute the total vector to use to get the index in kpb
                G = kset_wannier.get<double>(ik)->gkvec().gvec(gvec_index_t::local(ig));
                G += r3::vector<int>(nncell(0, ik, ib), nncell(1, ik, ib), nncell(2, ik, ib));
                int ig_ = gvec_kpb[index_ikpb]->index_by_gvec(G); // kpoints_[ikpb]->gkvec_->index_by_gvec(G);
                if (ig_ == -1) {
                    continue;
                }
                for (int iband = 0; iband < num_bands; iband++) {
                    aux_psi_kpb->pw_coeffs(ig, wf::spin_index(0), wf::band_index(iband)) =
                            wf_kpb[index_ikpb](ig_, iband);
                }
            } // end ig


            if( is_device_memory(mem) ) {
                Swf_k->allocate(memory_t::device);
                Swf_k->copy_to(memory_t::device);
                aux_psi_kpb->allocate(memory_t::device);
                aux_psi_kpb->copy_to(memory_t::device);
            }
            wf::inner(kset_wannier.ctx().spla_context(), mem, wf::spin_range(0), *aux_psi_kpb, wf::band_range(0, num_bands),
                      *Swf_k, wf::band_range(0, num_bands), Mbk, 0, 0);
                      
            if( is_device_memory(mem) ) {
                Swf_k->deallocate(memory_t::device);
                aux_psi_kpb->deallocate(memory_t::device);
                Mbk.copy_to(memory_t::host);
            }
            
            for (int n = 0; n < num_bands; n++) {
                for (int m = 0; m < num_bands; m++) {
                    M(m, n, ib, ik) = std::conj(Mbk(n, m));
                }
            }
        }
    } // end ik
    std::cout << "Mmn calculated.\n";
    std::cout << "starting broadcast...\n";
    for (int ik = 0; ik < kset_wannier.num_kpoints(); ik++) {
        int local_rank = kset_wannier.spl_num_kpoints().location(typename kp_index_t::global(ik)).ib;
        kset_wannier.ctx().comm_k().bcast(M.at(memory_t::host, 0, 0, 0, ik), num_bands * num_bands * nntot, local_rank);
    }
}

/// Generate the necessary data for the W90 input.
/** Wave-functions:
 * \f[
 *  \psi_{n{\bf k}} ({\bf r}) = \sum_{\bf G} e^{i({\bf G+k}){\bf r}} C_{n{\bf k}}({\bf G})
 * \f]
 *
 *  Matrix elements:
 *  \f{eqnarray*}{
 *  M_{nn'} &= \int e^{-i{\bf qr}}  \psi_{n{\bf k}}^{*} ({\bf r})  \psi_{n'{\bf k+q}} ({\bf r}) d{\bf r} =
 *    \sum_{\bf G} e^{-i({\bf G+k}){\bf r}} C_{n{\bf k}}^{*}({\bf G})
 *    \sum_{\bf G'} e^{i({\bf G'+k+q}){\bf r}} C_{n{\bf k+q}}({\bf G'}) e^{-i{\bf qr}} = \\
 *    &= \sum_{\bf GG'} \int e^{i({\bf G'-G}){\bf r}} d{\bf r}  C_{n{\bf k}}^{*}({\bf G}) C_{n{\bf k+q}}({\bf G'}) =
 *    \sum_{\bf G}  C_{n{\bf k}}^{*}({\bf G}) C_{n{\bf k+q}}({\bf G})
 *  \f}
 *
 *  Let's rewrite \f$ {\bf k + q} = {\bf \tilde G} + {\bf \tilde k} \f$. Now, through the property of plane-wave
 *  expansion coefficients \f$ C_{n{\bf k+q}}({\bf G}) = C_{n{\bf \tilde k}}({\bf G + \tilde G}) \f$ it follows that
 *  \f[
 *    M_{nn'} = \sum_{\bf G} C_{n{\bf k}}^{*}({\bf G}) C_{n{\bf \tilde k}}({\bf G + \tilde G})
 *  \f]
 */
void
K_point_set::generate_w90_coeffs() // sirius::K_point_set& k_set__)
{
    /* initializing variables and objects needed */
    PROFILE("sirius::K_point_set::generate_w90_coeffs");
    std::cout << "\n\n\nWannierization!!!!\n\n\n";
    K_point_set& kset_sirius = *this;
    K_point_set kset_wannier(this->ctx());
    std::vector<k_info> k_temp;

    int nntot, num_wann, num_bands;
    mdarray<int, 2> nnlist({kset_sirius.num_kpoints(), 12});        
    mdarray<int, 3> nncell({3, kset_sirius.num_kpoints(), 12}); 
    mdarray<int, 1> exclude_bands({ctx().num_bands()});   
    mdarray<double, 2> kp;  
    if( ctx().comm().rank() == 0 ) {
        read_nnkp( num_wann, nntot, nnlist, nncell, exclude_bands, kp );
    }
    std::cout << "read wannier90 .nnkp file" << std::endl;
    ctx().comm().bcast(&nntot, 1, 0);
    ctx().comm().bcast(nnlist.at(memory_t::host), nnlist.size(0) * nnlist.size(1), 0 );
    ctx().comm().bcast(nncell.at(memory_t::host), nncell.size(0) * nncell.size(1) * nncell.size(2) , 0 );
    ctx().comm().bcast(&num_wann, 1, 0 );
    ctx().comm().bcast(exclude_bands.at(memory_t::host), exclude_bands.size(0), 0 );
    ctx().comm().bcast( kp.at(memory_t::host), kp.size(0)*kp.size(1), 0 );

    for (int ik = 0; ik < kp.size(1); ik++) {
        kset_wannier.add_kpoint( &kp(0, ik), 1./kp.size(1) );
    }
    kset_wannier.initialize();

    /* initialize pw_coeffs in kset_wannier */
    match_kpoints( kset_sirius, kset_wannier, k_temp );

    std::vector<int> band_index_tot; // band_index_tot[iband] gives the index of iband in the full band vector
    for (int iband = 0; iband < exclude_bands.size(); iband++) {
        int band_fortran = iband + 1;
        bool is_excluded =
                (std::find_if(exclude_bands.at(memory_t::host), exclude_bands.at(memory_t::host) + exclude_bands.size(),
                              [&, band_fortran](int const& band_excluded) {
                                  return (band_excluded == band_fortran);
                              }) != exclude_bands.at(memory_t::host) + exclude_bands.size());

        if (!is_excluded) {
            band_index_tot.push_back(iband);
        }
        
        std::cout << "is excluded: " << is_excluded << " band_index_tot: " << band_index_tot << std::endl;
    }

    num_bands = band_index_tot.size();
    PROFILE_STOP("sirius::K_point_set::generate_w90_coeffs::wannier_setup");
    rotate_wavefunctions(*this, kset_wannier, k_temp, num_bands, band_index_tot);
    num_wann = ctx_.unit_cell().num_ps_atomic_wf().first;

    mdarray<std::complex<double>, 3> A({num_bands, num_wann, kset_wannier.num_kpoints()});
    A.zero();

    calculate_Amn(kset_wannier, num_bands, num_wann, A);

    if (ctx().comm().rank() == 0) {
        write_Amn(A, kp.size(1), num_bands, num_wann);
    }

    std::vector<std::shared_ptr<fft::Gvec>> gvec_kpb;
    std::vector<mdarray<std::complex<double>, 2>> wf_kpb;
    std::vector<int> ikpb_index(kset_wannier.num_kpoints(), -1);

    send_receive_kpb(gvec_kpb, wf_kpb, kset_wannier, ikpb_index, nntot, nnlist, num_bands);

    mdarray<std::complex<double>, 4> M({num_bands, num_bands, nntot, kset_wannier.num_kpoints()});
    M.zero();

    calculate_Mmn(M, kset_wannier, num_bands, gvec_kpb, wf_kpb, ikpb_index, nntot, nnlist, nncell);

    if (ctx().comm().rank() == 0) {
        write_Mmn(M, nnlist, nncell, kp.size(1), nntot, num_bands);
    }

    // Initialize eigval with the value of the energy dispersion

    mdarray<double, 2> eigval({num_bands, kp.size(1)}); // input

    for (int ik = 0; ik < kp.size(1); ik++) {
        int local_rank = kset_wannier.spl_num_kpoints().location(typename kp_index_t::global(ik)).ib;
        if (kset_wannier.ctx().comm_k().rank() == local_rank) {
            for (int iband = 0; iband < num_bands; iband++) {
                eigval(iband, ik) =
                        kset_wannier.get<double>(ik)->band_energies(0)[iband] * ha2ev; // sirius saves energy in
                                                                                 // Hartree, we need it in eV
            }
        }
        kset_wannier.ctx().comm_k().bcast(eigval.at(memory_t::host, 0, ik), num_bands, local_rank); // TODO: remove
    }

    if(ctx().comm().rank() == 0){
        write_eig(eigval, num_bands, kp.size(1));
    }
    /*
        if (kset_wannier.ctx().comm_k().rank() == 0) {
            std::cout << "Starting wannier_run..." << std::endl;

            // compute wannier orbitals
            // define additional arguments
            mdarray<std::complex<double>, 3> U_matrix(num_wann, num_wann, kp.size(1)); // output
            mdarray<std::complex<double>, 3> U_dis(num_bands, num_wann, kp.size(1));   // output
            mdarray<fortran_bool, 2> lwindow(num_bands, kp.size(1));                   // output
            mdarray<double, 2> wannier_centres(3, num_wann);                         // output
            mdarray<double, 1> wannier_spreads(num_wann);                            // output
            mdarray<double, 1> spread_loc(3);                                        // output-op

            write_eig(eigval, num_bands, kp.size(1));

            U_matrix.zero();
            U_dis.zero();
            lwindow.zero();
            wannier_centres.zero();
            wannier_spreads.zero();
            spread_loc.zero();

            PROFILE_START("sirius::K_point_set::generate_w90_coeffs::wannier_run");

            wannier_run_(seedname, this->ctx().cfg().parameters().ngridk().data(), &kp.size(1),
                         real_lattice.at(memory_t::host), recip_lattice.at(memory_t::host),
                         kpt_lattice.at(memory_t::host), &num_bands, &num_wann, &nntot, &num_atoms, atomic_symbol,
                         atoms_cart.at(memory_t::host), &gamma_only, M.at(memory_t::host),
                         A.at(memory_t::host), eigval.at(memory_t::host), U_matrix.at(memory_t::host),
                         U_dis.at(memory_t::host), lwindow.at(memory_t::host),
                         wannier_centres.at(memory_t::host), wannier_spreads.at(memory_t::host),
                         spread_loc.at(memory_t::host), length_seedname, length_atomic_symbol);
            std::cout << "Wannier_run succeeded. " << std::endl;
        }
    */
    PROFILE_STOP("sirius::K_point_set::generate_w90_coeffs::wannier_run");
}

} // namespace sirius
///////////#endif // SIRIUS_WANNIER90
