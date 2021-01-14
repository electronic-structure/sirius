// Copyright (c) 2013-2017 Anton Kozhevnikov, Thomas Schulthess
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

/** \file input.hpp
 *
 *  \brief Contains input parameters structures.
 *
 *  \todo Some of the parameters belong to SCF ground state mini-app. Mini-app should parse this values itself.
 *  \todo parse atomic coordinates and magnetic field separtely, not as 6D vector.
 */

#ifndef __INPUT_HPP__
#define __INPUT_HPP__

#include <list>
#include "constants.hpp"
#include "SDDK/geometry3d.hpp"
#include "utils/json.hpp"
#include <iostream>

using namespace geometry3d;
using namespace nlohmann;

namespace sirius {

/// Parse the parameters of iterative solver.
struct Iterative_solver_input
{
    /// Type of the iterative solver.
    std::string type_{""};

    /// Number of steps (iterations) of the solver.
    int num_steps_{20};

    /// Size of the variational subspace is this number times the number of bands.
    int subspace_size_{2};

    /// Lock eigenvectors of the smallest eigenvalues when they have converged at restart
    bool locking_{true};

    /// Restart early when the ratio unconverged vs lockable vectors drops below this threshold
    /** When there's just a few vectors left unconverged, it can be more efficient to lock the converged
        ones, such that the dense eigenproblem solved in each Davidson iteration has lower dimension.
        Restarting has some overhead in that it requires updating wave functions **/
    double early_restart_{0.5};

    /// Tolerance for the eigen-energy difference \f$ |\epsilon_i^{old} - \epsilon_i^{new} | \f$.
    /** This parameter is reduced during the SCF cycle to reach the high accuracy of the wave-functions. */
    double energy_tolerance_{1e-2};

    /// Tolerance for the residual L2 norm.
    double residual_tolerance_{1e-6};

    /// Relative tolerance for the residual L2 norm. (0 means this criterion is effectively not used)
    double relative_tolerance_{0};

    /// Additional tolerance for empty states.
    /** Setting this variable to 0 will treat empty states with the same tolerance as occupied states. */
    double empty_states_tolerance_{0};

    /// Defines the flavour of the iterative solver.
    /** If converge_by_energy is set to 0, then the residuals are estimated by their norm. If converge_by_energy
        is set to 1 then the residuals are estimated by the eigen-energy difference. This allows to estimate the
        unconverged residuals and then compute only the unconverged ones. */
    int converge_by_energy_{1}; // TODO: rename, this is meaningless

    /// Minimum number of residuals to continue iterative diagonalization process.
    int min_num_res_{0};

    /// Number of singular components for the LAPW Davidson solver.
    int num_singular_{-1};

    /// Initialize eigen-values with previous (old) values.
    bool init_eval_old_{true};

    /// Tell how to initialize the subspace.
    /** It can be either "lcao", i.e. start from the linear combination of atomic orbitals or "random" –- start from
        the randomized wave functions. */
    std::string init_subspace_{"lcao"};

    void read(json const& parser)
    {
        if (parser.count("iterative_solver")) {
            auto section            = parser["iterative_solver"];
            type_                   = section.value("type", type_);
            num_steps_              = section.value("num_steps", num_steps_);
            subspace_size_          = section.value("subspace_size", subspace_size_);
            locking_                = section.value("locking", locking_);
            energy_tolerance_       = section.value("energy_tolerance", energy_tolerance_);
            residual_tolerance_     = section.value("residual_tolerance", residual_tolerance_);
            relative_tolerance_     = section.value("relative_tolerance", relative_tolerance_);
            empty_states_tolerance_ = section.value("empty_states_tolerance", empty_states_tolerance_);
            converge_by_energy_     = section.value("converge_by_energy", converge_by_energy_);
            min_num_res_            = section.value("min_num_res", min_num_res_);
            num_singular_           = section.value("num_singular", num_singular_);
            init_eval_old_          = section.value("init_eval_old", init_eval_old_);
            init_subspace_          = section.value("init_subspace", init_subspace_);
            early_restart_          = section.value("early_restart", early_restart_);
            std::transform(init_subspace_.begin(), init_subspace_.end(), init_subspace_.begin(), ::tolower);
        }
    }
};

/// Parse control input section.
/** The following part of the input file is parsed:
    \code{.json}
    "control" : {
      "mpi_grid_dims" : (1- 2- or 3-dimensional vector<int>) MPI grid layout
      "cyclic_block_size" : (int) PBLAS / ScaLAPACK block size
      "reduce_gvec" : (bool) use reduced G-vector set (reduce_gvec = true) or full set (reduce_gvec = false)
      "std_evp_solver_type" : (string) type of eigen-solver for the standard eigen-problem
      "gen_evp_solver_type" : (string) type of eigen-solver for the generalized eigen-problem
      "processing_unit" : (string) primary processing unit
      "fft_mode" : (string) serial or parallel FFT
    }
    \endcode
    Parameters of the control input sections do not in general change the numerics, but instead control how the
    results are obtained. Changing parameters in control section should not change the significant digits in final
    results.
 */
struct Control_input
{
    /// Dimensions of the MPI grid (if used).
    std::vector<int> mpi_grid_dims_;

    /// Block size for ScaLAPACK and ELPA.
    int cyclic_block_size_{-1};

    /// Reduce G-vectors by inversion symmetry.
    /** For real-valued functions like density and potential it is sufficient to store only half of the G-vectors
     *  and use the relation f(G) = f^{*}(-G) to recover second half of the plane-wave expansion coefficients. */
    bool reduce_gvec_{true};

    /// Standard eigen-value solver to use.
    std::string std_evp_solver_name_{""};

    /// Generalized eigen-value solver to use.
    std::string gen_evp_solver_name_{""};

    /// Coarse grid FFT mode ("serial" or "parallel").
    std::string fft_mode_{"serial"};

    /// Main processing unit to run on.
    std::string processing_unit_{""};

    /// Maximum allowed muffin-tin radius in case of LAPW.
    double rmt_max_{2.2};

    /// Tolerance of the spglib in finding crystal symmetries.
    double spglib_tolerance_{1e-4};

    /// Level of verbosity.
    /** The following convention in proposed:
     *    - 0: silent mode (no output is printed) \n
     *    - 1: basic output (low level of output) \n
     *    - 2: extended output (medium level of output) \n
     *    - 3: extensive output (high level of output) */
    int verbosity_{0};

    /// Level of internal verification.
    int verification_{0};

    /// Number of eigen-values that are printed to the standard output.
    int num_bands_to_print_{10};

    /// If true then performance of some compute-intensive kernels will be printed to the standard output.
    bool print_performance_{false};

    /// If true then memory usage will be printed to the standard output.
    bool print_memory_usage_{false};

    /// If true then the checksums of some arrays will be printed (useful during debug).
    bool print_checksum_{false};

    /// If true then the hashsums of some arrays will be printed.
    bool print_hash_{false};

    /// If true then the stress tensor components are printed at the end of SCF run.
    bool print_stress_{false};

    /// If true then the atomic forces are printed at the end of SCF run.
    bool print_forces_{false};

    /// If true then the timer statistics is printed at the end of SCF run.
    bool print_timers_{true};

    /// If true then the list of nearest neighbours for each atom is printed to the standard output.
    bool print_neighbors_{false};

    /// True if second-variational diagonalization is used in LAPW method.
    bool use_second_variation_{true};

    /// Control the usage of the GPU memory.
    /** Possible values are: "low", "medium" and "high". */
    std::string memory_usage_{"high"};

    /// Number of atoms in the beta-projectors chunk.
    int beta_chunk_size_{256};

    void read(json const& parser)
    {
        if (parser.count("control")) {
            auto section         = parser["control"];
            mpi_grid_dims_       = section.value("mpi_grid_dims", mpi_grid_dims_);
            cyclic_block_size_   = section.value("cyclic_block_size", cyclic_block_size_);
            std_evp_solver_name_ = section.value("std_evp_solver_type", std_evp_solver_name_);
            gen_evp_solver_name_ = section.value("gen_evp_solver_type", gen_evp_solver_name_);
            processing_unit_     = section.value("processing_unit", processing_unit_);
            fft_mode_            = section.value("fft_mode", fft_mode_);
            reduce_gvec_         = section.value("reduce_gvec", reduce_gvec_);
            rmt_max_             = section.value("rmt_max", rmt_max_);
            spglib_tolerance_    = section.value("spglib_tolerance", spglib_tolerance_);
            verbosity_           = section.value("verbosity", verbosity_);
            verification_        = section.value("verification", verification_);
            num_bands_to_print_  = section.value("num_bands_to_print", num_bands_to_print_);
            print_performance_   = section.value("print_performance", print_performance_);
            print_memory_usage_  = section.value("print_memory_usage", print_memory_usage_);
            print_checksum_      = section.value("print_checksum", print_checksum_);
            print_hash_          = section.value("print_hash", print_hash_);
            print_stress_        = section.value("print_stress", print_stress_);
            print_forces_        = section.value("print_forces", print_forces_);
            print_timers_        = section.value("print_timers", print_timers_);
            print_neighbors_     = section.value("print_neighbors", print_neighbors_);
            memory_usage_        = section.value("memory_usage", memory_usage_);
            beta_chunk_size_     = section.value("beta_chunk_size", beta_chunk_size_);

            auto strings = {&std_evp_solver_name_, &gen_evp_solver_name_, &fft_mode_, &processing_unit_,
                            &memory_usage_};
            for (auto s : strings) {
                std::transform(s->begin(), s->end(), s->begin(), ::tolower);
            }

            std::list<std::string> kw;
            kw = {"low", "medium", "high"};
            if (std::find(kw.begin(), kw.end(), memory_usage_) == kw.end()) {
                throw std::runtime_error("wrong memory_usage input");
            }
        }
    }
};

/// Parse parameters input section.
/** Most of this parameters control the behavior of sirius::DFT_ground_state class. */
struct Parameters_input
{
    /// Electronic structure method.
    std::string electronic_structure_method_{"none"};

    /// List of XC functions (typically contains exchange term and correlation term).
    std::vector<std::string> xc_functionals_;

    /// Type of core-states relativity in full-potential LAPW case.
    std::string core_relativity_{"dirac"};

    /// Type of valence states relativity in full-potential LAPW case.
    std::string valence_relativity_{"zora"};

    /// Number of bands.
    /** In spin-collinear case this is the number of bands for each spin channel. */
    int num_bands_{-1};

    /// Number of first-variational states.
    int num_fv_states_{-1};

    /// Width of Gaussian smearing function in the units of [Ha].
    double smearing_width_{0.01};

    /// Type of occupancy smearing.
    std::string smearing_{"gaussian"};

    /// Cutoff for plane-waves (for density and potential expansion) in the units of [a.u.^-1].
    double pw_cutoff_{0.0};

    /// Cutoff for augmented-wave functions.
    double aw_cutoff_{0.0}; // this is R_{MT} * |G+k|_{max}

    /// Cutoff for |G+k| plane-waves in the units of [a.u.^-1].
    double gk_cutoff_{0.0};

    /// Maximum l for APW functions.
    int lmax_apw_{8};

    /// Maximum l for density.
    int lmax_rho_{8};

    /// Maximum l for potential
    int lmax_pot_{8};

    /// Number of dimensions of the magnetization and effective magnetic field (0, 1 or 3).
    int num_mag_dims_{0};

    /// Scale muffin-tin radii automatically.
    int auto_rmt_{1};

    /// Regular k-point grid for the SCF ground state.
    std::vector<int> ngridk_{1, 1, 1};

    /// Shift in the k-point grid.
    std::vector<int> shiftk_{0, 0, 0};

    /// optional k-point coordinates
    std::vector<vector3d<double>> vk_;

    /// Number of SCF iterations.
    int num_dft_iter_{100};

    /// Tolerance in total energy change (in units of [Ha]).
    double energy_tol_{1e-8};

    /// Tolerance for the density root mean square (in units of [a.u.^-3]).
    /** RMS is computed as Sqrt( (1/Omega) \int \delta rho(r) \delta rho(r) dr) */
    double density_tol_{1e-8};

    /// True if this is a molecule calculation.
    bool molecule_{false};

    /// True if gamma-point (real) version of the PW code is used.
    bool gamma_point_{false};

    /// True if spin-orbit correction is applied.
    bool so_correction_{false};

    /// True if Hubbard (or U) correction is applied.
    bool hubbard_correction_{false};

    /// True if symmetry is used.
    bool use_symmetry_{true};

    /// Radius of atom nearest-neighbour cluster.
    double nn_radius_{-1};

    /// Effective screening medium.
    bool enable_esm_{false};

    /// Type of periodic boundary conditions.
    std::string esm_bc_{"pbc"};

    /// Reduction of the auxiliary magnetic field at each SCF step.
    double reduce_aux_bf_{0.0};

    /// Introduce extra charge to the system. Positive charge means extra holes, negative charge - extra electrons.
    double extra_charge_{0.0};

    /// xc density threshold (debug purposes)
    double xc_dens_tre_{-1.0};

    void read(json const& parser)
    {
        if (parser.count("parameters")) {
            auto section = parser["parameters"];
            electronic_structure_method_ = section.value("electronic_structure_method", electronic_structure_method_);
            std::transform(electronic_structure_method_.begin(), electronic_structure_method_.end(),
                           electronic_structure_method_.begin(), ::tolower);
            xc_functionals_.clear();
            /* read list of XC functionals */
            if (section.count("xc_functionals")) {
                xc_functionals_.clear();
                for (auto& label : section["xc_functionals"]) {
                    xc_functionals_.push_back(std::string(label));
                }
            }

            if (section.count("vdw_functionals")) {
                xc_functionals_.push_back(section["vdw_functionals"].get<std::string>());
            }

            core_relativity_ = section.value("core_relativity", core_relativity_);
            std::transform(core_relativity_.begin(), core_relativity_.end(), core_relativity_.begin(), ::tolower);

            valence_relativity_ = section.value("valence_relativity", valence_relativity_);
            std::transform(valence_relativity_.begin(), valence_relativity_.end(), valence_relativity_.begin(),
                           ::tolower);

            num_fv_states_  = section.value("num_fv_states", num_fv_states_);
            smearing_width_ = section.value("smearing_width", smearing_width_);
            smearing_       = section.value("smearing", smearing_);
            pw_cutoff_      = section.value("pw_cutoff", pw_cutoff_);
            aw_cutoff_      = section.value("aw_cutoff", aw_cutoff_);
            gk_cutoff_      = section.value("gk_cutoff", gk_cutoff_);
            lmax_apw_       = section.value("lmax_apw", lmax_apw_);
            lmax_rho_       = section.value("lmax_rho", lmax_rho_);
            lmax_pot_       = section.value("lmax_pot", lmax_pot_);
            num_mag_dims_   = section.value("num_mag_dims", num_mag_dims_);
            auto_rmt_       = section.value("auto_rmt", auto_rmt_);
            use_symmetry_   = section.value("use_symmetry", use_symmetry_);
            gamma_point_    = section.value("gamma_point", gamma_point_);
            ngridk_         = section.value("ngridk", ngridk_);
            shiftk_         = section.value("shiftk", shiftk_);
            num_dft_iter_   = section.value("num_dft_iter", num_dft_iter_);
            auto vk         = section.value("vk", std::vector<std::vector<double>>{});
            for (auto& vki : vk) {
                if (vki.size() != 3) {
                    throw std::runtime_error("parameters.vk expected to be of size 3");
                }
                vk_.emplace_back(vector3d<double>(vki));
            }
            energy_tol_     = section.value("energy_tol", energy_tol_);
            /* potential_tol is obsolete */
            density_tol_    = section.value("potential_tol", density_tol_);
            density_tol_    = section.value("density_tol", density_tol_);
            molecule_       = section.value("molecule", molecule_);
            nn_radius_      = section.value("nn_radius", nn_radius_);
            reduce_aux_bf_  = section.value("reduce_aux_bf", reduce_aux_bf_);
            extra_charge_   = section.value("extra_charge", extra_charge_);
            xc_dens_tre_    = section.value("xc_density_threshold", xc_dens_tre_);

            so_correction_ = section.value("spin_orbit", so_correction_);

            /* spin-orbit correction requires non-collinear magnetism */
            if (so_correction_) {
                num_mag_dims_ = 3;
            }

            hubbard_correction_ = section.value("hubbard_correction", hubbard_correction_);
        }
    }
};

struct NLCG_input
{
    /// CG max iterations
    int maxiter_{300};
    /// CG restart
    int restart_{10};
    /// backtracking search, step parameter
    double tau_{0.1};
    /// temperature in Kelvin
    double T_{300};
    /// scalar preconditioning of pseudo Hamiltonian
    double kappa_{0.3};
    /// CG tolerance
    double tol_{1e-9};
    /// smearing
    std::string smearing_{"FD"};
    /// Main processing unit to run on.
    std::string processing_unit_{""};

    void read(json const& parser)
    {
        if (parser.count("nlcg")) {
            auto section     = parser["nlcg"];
            maxiter_         = section.value("maxiter", maxiter_);
            restart_         = section.value("restart", restart_);
            tau_             = section.value("tau", tau_);
            T_               = section.value("T", T_);
            kappa_           = section.value("kappa", kappa_);
            tol_             = section.value("tol", tol_);
            smearing_        = section.value("smearing", smearing_);
            processing_unit_ = section.value("processing_unit", processing_unit_);
        }
    }
};

struct Hubbard_input
{
    int number_of_species{0};
    bool hubbard_correction_{false};
    bool simplified_hubbard_correction_{false};
    bool orthogonalize_hubbard_orbitals_{false};
    bool normalize_hubbard_orbitals_{false};
    bool hubbard_U_plus_V_{false};

    /** by default we use the atomic orbitals given in the pseudo potentials */
    int projection_method_{0};

    struct hubbard_orbital_t
    {
        int l{-1};
        int n{-1};
        std::string level;
        std::array<double, 6> coeff{0, 0, 0, 0, 0, 0};
        double occupancy{0};
        std::vector<double> initial_occupancy;
    };

    std::string wave_function_file_;
    std::map<std::string, hubbard_orbital_t> species_with_U;

    bool hubbard_correction() const
    {
        return hubbard_correction_;
    }

    void read(json const& parser);
};

}; // namespace sirius

#endif // __INPUT_HPP__
