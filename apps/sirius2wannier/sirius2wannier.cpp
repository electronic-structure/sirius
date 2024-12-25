#include "apps.hpp"
//#include "k_point/generate_w90_coeffs.hpp"
//#include "nlcglib/apply_hamiltonian.hpp"
#include "hamiltonian/check_wave_functions.hpp"

std::vector<std::array<double, 3>>  
load_coordinates( const std::string& fname__ )
{
    std::vector<std::array<double, 3>> kp;

    //read k coordinates from hdf5
    HDF5_tree fin(fname__, hdf5_access_t::read_only);
    std::cout << "read num_kpoints" << std::endl;
    int num_kpoints;
    fin["K_point_set"].read("num_kpoints", &num_kpoints, 1);
    std::cout << "num_kpoints: " << num_kpoints << std::endl;
    kp.resize(num_kpoints);
    for( int ik = 0; ik < num_kpoints; ik++ ) {
        fin["K_point_set"][ik].read("vk", &kp[ik][0], 3);
        std::cout << "ik = " << ik << " kp = { " << kp[ik][0] << " " << kp[ik][1] << " " << kp[ik][2] << std::endl;
    }
    return kp;
}


int
main(int argn, char** argv)
{
    cmd_args args(argn, argv,
                  {{"input=", "{string} input file name"}});

    sirius::initialize(1);


    /* get the input file name */
    auto fpath = args.value<fs::path>("input", "state.h5");

    if (fs::is_directory(fpath)) {
        fpath /= "sirius.h5";
    }

    if (!fs::exists(fpath)) {
        if (mpi::Communicator::world().rank() == 0) {
            std::cout << "input file does not exist" << std::endl;
        }
        exit(1);
    }
    auto fname = fpath.string();
    
    /* create simulation context */
    auto ctx = create_sim_ctx(fname, args);
    ctx->initialize();

    /* read the wf */
    auto kp = load_coordinates(fname);
    K_point_set kset(*ctx, kp);
    std::cout << "kset initialized.\n";
    kset.load(fname);

    /* initialize the ground state */
    DFT_ground_state dft(kset);
    auto& potential = dft.potential();
    auto& density   = dft.density();
    density.load(fname);
    density.generate_paw_density();
    //potential.load(fname);
    potential.generate(density, ctx->use_symmetry(), true);
    Hamiltonian0<double> H0(potential, true);
    
    /* checksum over wavefunctions */
    //for (auto it : kset.spl_num_kpoints()) {
    //    int ik = it.i;
    //    auto Hk = H0(*kset.get<double>(ik));
    //    for (auto is=0; is< ctx->num_spins(); is++) {
    //      std::cout << "ik: " << ik << " ispn : "<< is << " " ; 
    //      std::cout << kset.get<double>(ik)->spinor_wave_functions().checksum(memory_t::host, wf::spin_index(is), wf::band_range(0, ctx->num_bands())) << std::endl;
    //    }
    //}
    /* check if the wfs diagonalize the hamiltonian and if the eigenvalues are correct */

    //la::dmatrix<std::complex<double>> psiHpsi(ctx->num_bands(), ctx->num_bands());

    std::cout << "num_spins " << kset.ctx().num_spins() << std::endl;

    for (auto it : kset.spl_num_kpoints()) {
        int ik = it.i;
        std::cout << "ik = " << ik << std::endl;
        auto kp = kset.get<double>(ik);
        auto Hk = H0(*kp);

        /* check wave-functions */
        if (true || ctx->cfg().control().verification() >= 2) {
            if (ctx->num_mag_dims() == 3) {
                auto eval = kp->band_energies(0);
                check_wave_functions<double, std::complex<double>>(Hk, kp->spinor_wave_functions(), wf::spin_range(0, 2),
                                           wf::band_range(0, ctx->num_bands()), eval.data());
            } else {
                for (int ispn = 0; ispn < ctx->num_spins(); ispn++) {
                    auto eval = kp->band_energies(ispn);
                    check_wave_functions<double, std::complex<double>>(Hk, kp->spinor_wave_functions(), wf::spin_range(ispn),
                                               wf::band_range(0, ctx->num_bands()), eval.data());
                }
            }
        }

        //bool nc_mag     = (kset.ctx().num_mag_dims() == 3);
        //int num_spinors = (kset.ctx().num_mag_dims() == 1) ? 2 : 1;
        //int num_sc      = nc_mag ? 2 : 1;

        //auto hpsi = std::make_unique<wf::Wave_functions<double>>(kp->gkvec_sptr(), 
        //                                                         wf::num_mag_dims(kset.ctx().num_mag_dims()),
        //                                                         wf::num_bands(ctx->num_bands()), 
        //                                                         memory_t::host);
        //auto spsi = std::make_unique<wf::Wave_functions<double>>(kp->gkvec_sptr(), 
        //                                                         wf::num_mag_dims(kset.ctx().num_mag_dims()),
        //                                                         wf::num_bands(ctx->num_bands()), 
        //                                                         memory_t::host);

        //for (int ispin_step = 0; ispin_step < num_spinors; ispin_step++) {
        //    auto sr = nc_mag ? wf::spin_range(0, 2) : wf::spin_range(ispin_step);
        //    std::cout << "ik= " << ik << " ispin_step = " << ispin_step; 
        //    /* get H|psi> */
        //    Hk.apply_h_s<std::complex<double>>(sr, wf::band_range(0, ctx->num_bands()), kp->spinor_wave_functions(), hpsi.get(), spsi.get());

        //    /* get <psi|H|psi> */
        //    wf::inner(kset.ctx().spla_context(), memory_t::host, sr, kp->spinor_wave_functions(), wf::band_range(0, ctx->num_bands()),
        //              *hpsi, wf::band_range(0, ctx->num_bands()), psiHpsi, 0, 0);

        //    /* check elements that are large compared to the threshold */
        //    std::vector<std::pair<int,int>> indices;
        //    std::vector<double> comp;
        //    for( int ibnd = 0; ibnd < ctx->num_bands(); ++ibnd ) {
        //        for( int jbnd = 0; jbnd < ctx->num_bands(); ++jbnd ) {
        //            auto comparison = (ibnd == jbnd) ? kset.get<double>(ik)->band_energies(ispin_step)[ibnd]  : 0.;

        //            if( std::abs ( psiHpsi(jbnd,ibnd) - comparison ) > 1.e-07  ) {
        //                indices.push_back(std::pair<int,int>(jbnd, ibnd));
        //                comp.push_back(comparison);
        //            }
        //        }
        //    }
        //    for( auto i = 0; i < indices.size(); ++i ) {
        //        auto& i_ = indices[i];
        //        std::cout << "      element = (";
        //        std::cout << i_.first << ", " << i_.second << ") = " << psiHpsi(i_.first, i_.second) << " comparison : " << comp[i] << std::endl;
        //    }
        //}//ispin_step
    }//kpoint

    exit(0);

    //kset.generate_w90_coeffs();
}
