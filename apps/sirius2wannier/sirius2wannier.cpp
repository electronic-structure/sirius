#include "apps.hpp"
#include "k_point/generate_w90_coeffs.hpp"
#include "nlcglib/apply_hamiltonian.hpp"

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
    DFT_ground_state dft(kset);
    
    auto& potential = dft.potential();
    auto& density   = dft.density();

    density.load(fname);
    density.generate_paw_density();
    //potential.load(fname);
    potential.generate(density, ctx->use_symmetry(), true);
    Hamiltonian0<double> H0(potential, true);

//    /* check if the wfs diagonalize the hamiltonian and if the eigenvalues are correct */
//    for (auto it : kset.spl_num_kpoints()) {
//        int ik = it.i;
//        auto Hk = H0(*kset.get<double>(ik));
//        wf::Wave_functions<double> hpsi(kset.get<double>(ik)->gkvec_sptr(), wf::num_mag_dims(ctx->num_mag_dims() == 3 ? 3 : 0),
//                               wf::num_bands(ctx->num_bands()), memory_t::host);
//        std::shared_ptr<wf::Wave_functions<double>> spsi = std::make_shared<wf::Wave_functions<double>>(kset.get<double>(ik)->gkvec_sptr(), wf::num_mag_dims(ctx->num_mag_dims() == 3 ? 3 : 0),
//                               wf::num_bands(ctx->num_bands()), memory_t::host);
//
//        apply_hamiltonian(H0, *kset.get<double>(ik), hpsi,
//                  kset.get<double>(ik)->spinor_wave_functions(), spsi);
//
//        la::dmatrix<std::complex<double>> psiHpsi(ctx->num_bands(), ctx->num_bands());
//        psiHpsi.zero();
//        wf::inner(kset.ctx().spla_context(), memory_t::host, wf::spin_range(0), kset.get<double>(ik)->spinor_wave_functions(), wf::band_range(0, ctx->num_bands()),
//                      hpsi, wf::band_range(0, ctx->num_bands()), psiHpsi, 0, 0);
//        std::cout << "ik = " << ik << std::endl;
//        for( int ibnd = 0; ibnd < ctx->num_bands(); ++ibnd ) {
//            for( int jbnd = 0; jbnd < ctx->num_bands(); ++jbnd ) {
//                if( ibnd != jbnd )  { 
//                    RTE_ASSERT( std::abs ( psiHpsi(jbnd,ibnd) ) < 1.e-12  )
//                }
//                else {
//                    std::cout << psiHpsi(jbnd,ibnd) << " " << kset.get <double>(ik)->band_energies(0)[ibnd] << std::endl;
//                    RTE_ASSERT( std::abs( psiHpsi(jbnd,ibnd) - kset.get <double>(ik)->band_energies(0)[ibnd] ) < 1.e-12 )
//                }
//            }
//        }
//    }


    kset.generate_w90_coeffs();
}
