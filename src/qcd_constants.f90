MODULE qcd_constants
  IMPLICIT NONE
  PUBLIC
  
  ! default values
  INTEGER, PARAMETER:: nf_default = 5
  REAL(KIND(1d0)),PARAMETER:: CA_default = 3d0
  REAL(KIND(1d0)),PARAMETER:: CF_default = 4d0/3d0
  REAL(KIND(1d0)),PARAMETER:: TF_default = 0.5d0
  REAL(KIND(1d0)),PARAMETER:: nfTF_default = nf_default*TF_default

  ! the used values
  REAL(KIND(1d0)):: CA = CA_default
  REAL(KIND(1d0)):: CF = CF_default
  REAL(KIND(1d0)):: TF = TF_default
  REAL(KIND(1d0)):: nfTF = nfTF_default
  
  INTEGER:: nf = nf_default
  INTEGER:: nf_up = nf_default/2
  INTEGER:: nf_down = (nf_default+1)/2
  REAL(KIND(1d0)):: rnf = DBLE(nf_default)

  ! beta functions in terms of as/2/pi expansion
  ! should be evaluated by routine GET_QCD_BETA
  REAL(KIND(1d0)):: beta0, beta1, beta2, beta3, beta4
  
  ! gamma_m functions in terms of as/2/pi expansion
  ! should be evaluated by routine GET_QCD_GAMMA_MASS
  REAL(KIND(1d0)):: gammam0, gammam1, gammam2, gammam3
  
  ! inverse decoupling relation for alpha_s
  ! alpha_s(mu,nf+1)=alpha_s(mu,nf)*(1+\sum_{m,n}zeta_alphas_inv_X_mn*(alpha_s(mu,nf)/2pi)^m* log^n(mu^2/m_X^2)
  ! and the decoupling relation
  ! alpha_s(mu,nf-1)=alpha_s(mu,nf)*(1+\sum_{m,n}zeta_alphas_X_mn*(alpha_s(mu,nf)/2/pi)^m*log^n(mu^2/m_X^2)
  ! where X=MSbar,OS,SI 
  ! they will be used only when as_decoupling_on=.TRUE.
  LOGICAL::as_decoupling_on=.TRUE.
  REAL(KIND(1d0))::zeta_alphas_MSbar_11, zeta_alphas_MSbar_22, zeta_alphas_MSbar_21, zeta_alphas_MSbar_20
  REAL(KIND(1d0))::zeta_alphas_MSbar_33, zeta_alphas_MSbar_32, zeta_alphas_MSbar_31, zeta_alphas_MSbar_30
  REAL(KIND(1d0))::zeta_alphas_OS_11, zeta_alphas_OS_22, zeta_alphas_OS_21, zeta_alphas_OS_20
  REAL(KIND(1d0))::zeta_alphas_OS_33, zeta_alphas_OS_32, zeta_alphas_OS_31, zeta_alphas_OS_30
  REAL(KIND(1d0))::zeta_alphas_SI_11, zeta_alphas_SI_22, zeta_alphas_SI_21, zeta_alphas_SI_20
  REAL(KIND(1d0))::zeta_alphas_SI_33, zeta_alphas_SI_32, zeta_alphas_SI_31, zeta_alphas_SI_30
  REAL(KIND(1d0))::zeta_alphas_inv_MSbar_11, zeta_alphas_inv_MSbar_22, zeta_alphas_inv_MSbar_21, zeta_alphas_inv_MSbar_20
  REAL(KIND(1d0))::zeta_alphas_inv_MSbar_33, zeta_alphas_inv_MSbar_32, zeta_alphas_inv_MSbar_31, zeta_alphas_inv_MSbar_30
  REAL(KIND(1d0))::zeta_alphas_inv_OS_11, zeta_alphas_inv_OS_22, zeta_alphas_inv_OS_21, zeta_alphas_inv_OS_20
  REAL(KIND(1d0))::zeta_alphas_inv_OS_33, zeta_alphas_inv_OS_32, zeta_alphas_inv_OS_31, zeta_alphas_inv_OS_30
  REAL(KIND(1d0))::zeta_alphas_inv_SI_11, zeta_alphas_inv_SI_22, zeta_alphas_inv_SI_21, zeta_alphas_inv_SI_20
  REAL(KIND(1d0))::zeta_alphas_inv_SI_33, zeta_alphas_inv_SI_32, zeta_alphas_inv_SI_31, zeta_alphas_inv_SI_30
  ! Note: the 4-loop (inverse) decoupling relations are only known in QCD (SU(3) gauge theory)
  REAL(KIND(1d0))::zeta_alphas_MSbar_44,zeta_alphas_MSbar_43,zeta_alphas_MSbar_42,zeta_alphas_MSbar_41,zeta_alphas_MSbar_40
  REAL(KIND(1d0))::zeta_alphas_OS_44,zeta_alphas_OS_43,zeta_alphas_OS_42,zeta_alphas_OS_41,zeta_alphas_OS_40
  REAL(KIND(1d0))::zeta_alphas_SI_44,zeta_alphas_SI_43,zeta_alphas_SI_42,zeta_alphas_SI_41,zeta_alphas_SI_40
  REAL(KIND(1d0))::zeta_alphas_inv_MSbar_44,zeta_alphas_inv_MSbar_43,zeta_alphas_inv_MSbar_42,&
       zeta_alphas_inv_MSbar_41,zeta_alphas_inv_MSbar_40
  REAL(KIND(1d0))::zeta_alphas_inv_OS_44,zeta_alphas_inv_OS_43,zeta_alphas_inv_OS_42,&
       zeta_alphas_inv_OS_41,zeta_alphas_inv_OS_40
  REAL(KIND(1d0))::zeta_alphas_inv_SI_44,zeta_alphas_inv_SI_43,zeta_alphas_inv_SI_42,&
       zeta_alphas_inv_SI_41,zeta_alphas_inv_SI_40

  ! inverse decoupling relation for MSbar mass
  ! m(mu,nf+1)=m(mu,nf)*(1+\sum_{m,n}zeta_mass_inv_X_mn*(alpha_s(mu,nf)/2pi)^m* log^n(mu^2/m_X^2)
  ! and the decoupling relation
  ! m(mu,nf-1)=m(mu,nf)*(1+\sum_{m,n}zeta_mass_X_mn*(alpha_s(mu,nf)/2/pi)^m*log^n(mu^2/m_X^2)
  ! where X=MSbar,OS,SI
  ! they will be used only when mass_decoupling_on=.TRUE.
  LOGICAL::mass_decoupling_on=.TRUE.
  REAL(KIND(1d0))::zeta_mass_MSbar_20, zeta_mass_MSbar_21, zeta_mass_MSbar_22
  REAL(KIND(1d0))::zeta_mass_MSbar_33, zeta_mass_MSbar_32, zeta_mass_MSbar_31, zeta_mass_MSbar_30
  REAL(KIND(1d0))::zeta_mass_OS_20, zeta_mass_OS_21, zeta_mass_OS_22
  REAL(KIND(1d0))::zeta_mass_OS_33, zeta_mass_OS_32, zeta_mass_OS_31, zeta_mass_OS_30
  REAL(KIND(1d0))::zeta_mass_SI_20, zeta_mass_SI_21, zeta_mass_SI_22
  REAL(KIND(1d0))::zeta_mass_SI_33, zeta_mass_SI_32, zeta_mass_SI_31, zeta_mass_SI_30
  REAL(KIND(1d0))::zeta_mass_inv_MSbar_20, zeta_mass_inv_MSbar_21, zeta_mass_inv_MSbar_22
  REAL(KIND(1d0))::zeta_mass_inv_MSbar_33, zeta_mass_inv_MSbar_32, zeta_mass_inv_MSbar_31, zeta_mass_inv_MSbar_30
  REAL(KIND(1d0))::zeta_mass_inv_OS_20, zeta_mass_inv_OS_21, zeta_mass_inv_OS_22
  REAL(KIND(1d0))::zeta_mass_inv_OS_33, zeta_mass_inv_OS_32, zeta_mass_inv_OS_31, zeta_mass_inv_OS_30
  REAL(KIND(1d0))::zeta_mass_inv_SI_20, zeta_mass_inv_SI_21, zeta_mass_inv_SI_22
  REAL(KIND(1d0))::zeta_mass_inv_SI_33, zeta_mass_inv_SI_32, zeta_mass_inv_SI_31, zeta_mass_inv_SI_30

  ! quark mass threshold
  REAL(KIND(1d0)),DIMENSION(6),PARAMETER::quark_masses_default = (/1d-10,1d-10,0.15d0,&
       1.45d0, 4.5d0, 175.0d0/)
  ! type of quark masses
  CHARACTER(len=2),DIMENSION(3),PARAMETER::masses_type_string = (/"MS","OS","SI"/)
  ! default Lambda values
  ! one-loop
  REAL(KIND(1d0)),DIMENSION(6),PARAMETER::LAMBDA_1L_default = (/0.2470d0,0.2470d0,0.2470d0,&
       0.2150d0,0.1650d0,0.0847d0/)
  ! two-loop
  REAL(KIND(1d0)),DIMENSION(6),PARAMETER::LAMBDA_2L_default = (/0.3500d0,0.3500d0,0.3500d0,&
       0.3260d0,0.2262d0,0.1510d0/)
  ! three-loop
  REAL(KIND(1d0)),DIMENSION(6),PARAMETER::LAMBDA_3L_default = (/0.3500d0,0.3500d0,0.3500d0,&
       0.3260d0,0.2262d0,0.1510d0/)
  ! four-loop
  REAL(KIND(1d0)),DIMENSION(6),PARAMETER::LAMBDA_4L_default = (/0.3500d0,0.3500d0,0.3500d0,&
       0.3260d0,0.2262d0,0.1510d0/)
  ! five-loop
  REAL(KIND(1d0)),DIMENSION(6),PARAMETER::LAMBDA_5L_default = (/0.3500d0,0.3500d0,0.3500d0,&
       0.3260d0,0.2262d0,0.1510d0/)
  REAL(KIND(1d0)),DIMENSION(6)::LAMBDA_1L=LAMBDA_1L_default,LAMBDA_2L=LAMBDA_2L_default,&
       LAMBDA_3L=LAMBDA_3L_default,LAMBDA_4L=LAMBDA_4L_default,LAMBDA_5L=LAMBDA_5L_default

  ! the anomalous dimensions in terms of as/2/pi expansion
  ! cusp anomalous dimension
  ! three-loop see arXiv:1007.4005
  ! four-loop see arXiv:2002.04617
  ! gluon
  REAL(KIND(1d0))::QCD_gamma1g_cusp, QCD_gamma2g_cusp, QCD_gamma3g_cusp, QCD_gamma4g_cusp
  ! quark
  REAL(KIND(1d0))::QCD_gamma1q_cusp, QCD_gamma2q_cusp, QCD_gamma3q_cusp, QCD_gamma4q_cusp
  ! collinear/normal anomalous dimension
  ! see arXiv:0903.1126
  REAL(KIND(1d0))::QCD_gamma1g_collinear,QCD_gamma2g_collinear,QCD_gamma3g_collinear
  REAL(KIND(1d0))::QCD_gamma1q_collinear,QCD_gamma2q_collinear,QCD_gamma3q_collinear
END MODULE qcd_constants
