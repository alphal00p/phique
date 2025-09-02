MODULE global
  IMPLICIT NONE
  ! beam configurations
  ! colliding particles
  ! 1: hadron-hadron; 2: e+p collisions; 3: e+e- collisions
  INTEGER,PUBLIC::colpar=1
  REAL(KIND(1d0))::energy_beam1=7000d0
  REAL(KIND(1d0))::energy_beam2=7000d0
  REAL(KIND(1d0)),DIMENSION(2)::ebeam0
  INTEGER::nuclearA1=0
  INTEGER::nuclearA2=0
  INTEGER::nuclearZ1=0
  INTEGER::nuclearZ2=0
  INTEGER::UPC_photon_flux=1 ! only colpar=1
  REAL(KIND(1d0)),DIMENSION(2)::NNLO_RA,NNLO_aA,NNLO_wA
  ! order
  INTEGER::order=2
  ! coulomb
  ! 0: no Coulomb resummation
  ! 1: Leading-power Coulomb resummation
  ! 2: Next-to-leading power Coulomb resummation
  INTEGER::coulomb=0
  
  INTEGER,PUBLIC::varnum
  LOGICAL,PUBLIC::lunwei2
  INTEGER,PUBLIC::lwmax
  INTEGER,PUBLIC::nprint
  INTEGER,PARAMETER::maxprint=8
  ! For phase space parametrisation
  ! assuming the flux is tau**(-ntau+1)
  ! or dL/dW is tau**(-ntau)
  ! where tau=x1*x2
  REAL(KIND(1d0)),PUBLIC::xp1,xp2
  REAL(KIND(1d0)),PUBLIC::ntau=2d0
  INTEGER,PUBLIC::n_NNLO_energies
  INTEGER::IE_NLO_350, IE_NLO_500, IE_NLO_1000, IE_NLO_5000, IE_NLO_10000
  INTEGER::IE_NNLO_350, IE_NNLO_500, IE_NNLO_1000
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE,PUBLIC::total_DNNLO_xs
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE,PUBLIC::total_DNNLO_xserr
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE,PUBLIC::NNLO_energies
  INTEGER,PUBLIC::n_NLO_energies
  REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE,PUBLIC::NLO_energies
  REAL(KIND(1d0)),DIMENSION(:,:),ALLOCATABLE,PUBLIC::total_NLO_xs
  ! masses from PDG
  REAL(KIND(1d0)),PARAMETER::emass_PDG=0.51099895d-3
  REAL(KIND(1d0)),PARAMETER::mumass_PDG=0.1056583755d0
  REAL(KIND(1d0)),PARAMETER::taumass_PDG=1.77686d0
  REAL(KIND(1d0)),PARAMETER::umass_PDG=41d-3
  REAL(KIND(1d0)),PARAMETER::dmass_PDG=41d-3
  REAL(KIND(1d0)),PARAMETER::smass_PDG=0.150d0
  REAL(KIND(1d0)),PARAMETER::cmass_PDG=1.5d0
  REAL(KIND(1d0)),PARAMETER::bmass_PDG=4.5d0
  REAL(KIND(1d0)),PARAMETER::tmass_PDG=172.69d0
  REAL(KIND(1d0)),PARAMETER::wmass_PDG=80.377d0
  REAL(KIND(1d0)),PARAMETER::zmass_PDG=91.1876d0
  REAL(KIND(1d0)),PARAMETER::zmass_default=91.188d0
  ! coupling constants
  REAL(KIND(1d0))::alphaemm1=137.036d0
  REAL(KIND(1d0))::alphaemm1_default=137d0
  REAL(KIND(1d0))::alphasMZ=0.118d0
  REAL(KIND(1d0))::alphasMZ_default=0.118d0
  INTEGER::alphas_nloop=4
  REAL(KIND(1d0))::Gfermi=1.1663787d-5
  REAL(KIND(1d0))::alphaMZm1=127.955d0
  REAL(KIND(1d0))::DalphahadMZ_PDG=0.02766d0 ! 5 flavour number scheme
                                             ! uncertainty is +-0.00007
  INTEGER::alpha_scheme=0
  INTEGER::scale=0
  REAL(KIND(1d0))::fscalevalue=173d0
  REAL(KIND(1d0))::muR_over_ref=1d0
  REAL(KIND(1d0))::rw_RScale_down=0.5d0
  REAL(KIND(1d0))::rw_RScale_up=2.0d0
  LOGICAL::reweight_Scale=.FALSE.
  REAL(KIND(1d0))::scale_default=91.188d0
  REAL(KIND(1d0))::Wgaga=1d0
  ! event generation
  INTEGER::nmc=100000
  ! histogram output
  LOGICAL::topdrawer_output=.FALSE.
  LOGICAL::gnuplot_output=.TRUE.
  LOGICAL::root_output=.FALSE.
  LOGICAL::hwu_output=.TRUE.
  LOGICAL::plot_output=.FALSE.
  INTEGER::ho_nscale=0 ! number of scales
  INTEGER::tot_nrwgt=0
  ! scale reweight stuff
  REAL(KIND(1d0)),DIMENSION(3)::rw_Rscales
  !REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::rwgt_values
  REAL(KIND(1d0)),DIMENSION(3)::LO_wgtxsecmu,NLO_wgtxsecmu,NNLO_wgtxsecmu
  ! input files
  CHARACTER(len=20),PARAMETER::Input_File="run.inp"
  CHARACTER(len=24),PARAMETER::input_dir="./input/"
  CHARACTER(len=24),PARAMETER::tmp_dir="./tmp/"
  CHARACTER(len=24),PARAMETER::output_dir="./output/"
  CHARACTER(len=24),PARAMETER::grid_dir="./grid/"
  ! q2max is for iww
  REAL(KIND(1d0))::q2max=1d0
  INTEGER::nunit1=6
  ! ,nunit2,nunit3,nunit30
  ! heavy quark related
  ! quark flavor (4: c; 5: b; 6: t)
  INTEGER::quark=6
  ! 3*eQ
  INTEGER::Q_charge=2
  REAL(KIND(1d0))::MQ=173d0
  REAL(KIND(1d0)),PARAMETER::MQ_default=173d0
  ! number of light quarks
  INTEGER::nq=5 ! NLO does not depend on nq
END MODULE global
