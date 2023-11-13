import astropy.units as u

default_par = dict(
                 ##############
                 # COSMO params
                 ##############
                 cosmo_code = 'camb',
                 cosmo_input_camb=dict(f_NL=0,H0=67.36,cosmomc_theta=None,ombh2=0.02237, omch2=0.12, 
                               omk=0.0, neutrino_hierarchy='degenerate', 
                               num_massive_neutrinos=3, mnu=0.06, nnu=3.046, 
                               YHe=None, meffsterile=0.0, standard_neutrino_neff=3.046, 
                               TCMB=2.7255, tau=None, deltazrei=None, bbn_predictor=None, 
                               theta_H0_range=[10, 100],w=-1.0, wa=0., cs2=1.0, 
                               dark_energy_model='ppf',As=2.1e-09, ns=0.9649, nrun=0, 
                               nrunrun=0.0, r=0.0, nt=None, ntrun=0.0, 
                               pivot_scalar=0.05, pivot_tensor=0.05,
                               parameterization=2,halofit_version='mead2020'),
                 cosmo_input_class=dict(f_NL=0,H0=67.36,omega_b=0.02237, omega_cdm=0.12, 
                               A_s=2.1e-9,n_s=0.9649,
                               N_ncdm=3, m_ncdm='0.02,0.02,0.02', N_ur = 0.00641,
                               output='mPk,mTk'),
                 ###############
                 # ASTRO params
                 ###############
                 model_type='LF',
                 model_name='SchCut', 
                 model_par={'phistar':9.6e-11*u.Lsun**-1*u.Mpc**-3,
                 'Lstar':2.1e6*u.Lsun,'alpha':-1.87,'Lmin':5000*u.Lsun},
                 hmf_model='ST',
                 bias_model='ST99',
                 bias_par={}, #Otherwise, write a dict with the corresponding values
                 nu=115*u.GHz,
                 nuObs=30*u.GHz,
                 Mmin=1e9*u.Msun,
                 Mmax=1e15*u.Msun,
                 nM=500,
                 Lmin=10*u.Lsun,
                 Lmax=1e8*u.Lsun,
                 nL=5000,
                 v_of_M=None,
                 line_incli=True,
                 ###########
                 # Pk params
                 ###########
                 kmin = 1e-2*u.Mpc**-1,
                 kmax = 10.*u.Mpc**-1,
                 nk = 100,
                 k_kind = 'log',
                 sigma_scatter=0.,
                 fduty=1.,
                 do_onehalo=False,
                 do_Jysr=False,
                 do_RSD=True,
                 sigma_NL=7*u.Mpc,
                 nmu=1000,
                 FoG_damp='Lorentzian',
                 smooth=False,
                 do_conv_Wkmin = False,
                 nonlinear=False,
                 ############
                 #VID params
                 ############
                 smooth_VID = True,
                 Tmin_VID=1.0e-2*u.uK,
                 Tmax_VID=100.*u.uK,
                 Nbin_hist=100,
                 linear_VID_bin=False,
                 subtract_VID_mean=False,
                 #VID precision parameters
                 Lsmooth_tol=7,
                 T0_Nlogsigma=4,
                 fT0_min=1e-5*u.uK**-1,
                 fT0_max=1e4*u.uK**-1,
                 fT_min=1e-5*u.uK**-1,
                 fT_max=1e5*u.uK**-1,
                 nfT0=1000,
                 nT=2**18,
                 n_leggauss_nodes_FT='../nodes1e5.txt',
                 n_leggauss_nodes_IFT='../nodes1e4.txt',
                 sigma_PT_stable=0.0*u.uK,
                 #############
                 #OBS PARAMS #
                 #############
                 Tsys_NEFD=40*u.K,
                 Nfeeds=19,
                 beam_FWHM=4.1*u.arcmin,
                 Delta_nu=8*u.GHz,
                 dnu=15.6*u.MHz,
                 tobs=6000*u.hr, 
                 Omega_field=2.25*u.deg**2,
                 Nfield=1,
                 N_FG_par = 1,
                 N_FG_perp = 1,
                 do_FG_wedge = False,
                 a_FG = 0.*u.Mpc**-1,
                 b_FG = 0.,
                # PARAMETERS USED FOR LIMLAM SIMULATIONS
                catalogue_file = 
                    'limlam_mocker/catalogues/default_catalogue.npz',
                    # Location of peak-patch catalog file
                map_output_file = 'limlam_output.npz' # Output file location
              )
