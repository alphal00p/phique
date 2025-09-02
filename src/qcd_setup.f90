MODULE qcd_setup
  IMPLICIT NONE
  PRIVATE
  PUBLIC::SET_QCD_NF, SET_QCD_GROUP, SET_QCD_BETA, SET_QCD_GAMMA_MASS, SET_QCD_as_decoupling, SET_QCD_mass_decoupling
  PUBLIC::SET_QCD_GAMMA_CUSP,SET_QCD_GAMMA_COLLINEAR
CONTAINS
  SUBROUTINE SET_QCD_NF(nf_in)
    USE qcd_constants
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nf_in
    INTEGER::init=0
    SAVE init
    INTEGER::n_loop
    nf = nf_in
    rnf = DBLE(nf_in)
    nf_up = nf_in/2
    nf_down = (nf_in+1)/2
    nfTF = rnf*TF
    CALL SET_QCD_BETA(4)
    IF(CA.NE.CA_default.OR.CF.NE.CF_default.OR.TF.NE.TF_default)THEN
       ! we only have 3-loop decoupling relations for SU(N) N=/=3
       n_loop=3
    ELSE
       ! we have 4-loop decoupling relations for QCD SU(3)
       n_loop=4
    ENDIF
    IF(init.EQ.0)THEN
       CALL SET_QCD_as_decoupling(n_loop,.FALSE.)
    ELSE
       CALL SET_QCD_as_decoupling(n_loop,.TRUE.)
    ENDIF
    init=1
    RETURN
  END SUBROUTINE SET_QCD_NF

  SUBROUTINE SET_QCD_GROUP(CA_in, CF_in, TF_in)
    USE qcd_constants
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::CA_in, CF_in, TF_in
    CA = CA_in
    CF = CF_in
    TF = TF_in
    nfTF = rnf*TF
    CALL SET_QCD_BETA(4) ! beta function up to 5-loop is known 
                         ! caveat: beta_3 is correct only for SU(N) not for others like SO(N)
    IF(CA.NE.CA_default.OR.CF.NE.CF_default.OR.TF.NE.TF_default)THEN
       CALL SET_QCD_as_decoupling(3,.FALSE.) ! alpha_s decoupling relation is only known up to 3-loop for SU(N)
       ! caveat: the relation works for SU(N) not for others like SO(N)
    ELSE
       CALL SET_QCD_as_decoupling(4,.FALSE.) ! alpha_s decoupling relation is only known up to 3-loop for QCD SU(3)
    ENDIF
    RETURN
  END SUBROUTINE SET_QCD_GROUP

  SUBROUTINE SET_QCD_BETA(n)
    ! Set QCD beta functions up to beta_n (n+1 loop) in MSbar scheme
    ! beta_i is the coefficeints of the beta function expansion in terms of as/2pi
    ! das(t)/dt = beta(as(t)) = -as(t)(beta_0*as(t)/2pi+beta_1*(as(t)/2pi)**2+...)
    ! where t=log(Q**2)
    ! IN:
    ! n - get beta_i, i=0,...,n (max n <= 4)
    ! OUT:
    ! beta0=beta_0, beta1=beta_1,...
    USE qcd_constants
    USE global_constants
    IMPLICIT NONE
    INTEGER,INTENT(IN)::n
    REAL(KIND(1d0))::dAdA,dFdA,dFdF
    beta0=(11d0*CA-4d0*nfTF)/6d0
    IF(n.LE.0)RETURN
    beta1=(17d0*CA**2-nfTF*(6d0*CF+10d0*CA))/6d0
    IF(n.LE.1)RETURN
    beta2=(2857d0*CA**3-nfTF*(2830d0*CA**2+1230d0*CA*CF-108d0*CF**2)+&
         nfTF**2*(316d0*CA+264d0*CF))/432d0
    IF(n.LE.2)RETURN
    ! Eq.(8) in hep-ph/9701390
    beta3=(CA**4*(173448d0*zeta3+143993d0)+2592d0*CA**3*CF*(5d0-132d0*zeta3))&
         +nfTF*(-6d0*CA**3*(9432d0*zeta3+37799d0)+46d0*CA**2*CF*(2160d0*zeta3+7d0)&
                +72d0*CA*CF**2*(264d0*zeta3-1051d0)+22356d0*CF**3)&
         +nfTF**2*(12d0*CA**2*(4752d0*zeta3+2249d0)-64d0*CA*CF*(2862d0*zeta3-2021d0)&
                   +144d0*CF**2*(1464d0*zeta3-623d0))+16d0*nfTF**3*(53d0*CA+154d0*CF)
    beta3=beta3/7776d0
    IF(n.LE.3)RETURN
    ! Eq.(2.14) in arXiv:1701.01404
    dAdA=-3D0*CA**3*CF+37D0/24D0*CA**4
    dFdA=(CA**3*7D0/24D0+CA**2*CF*(-1D0/2D0))/2D0
    dFdF=(CA**2*13D0/24D0+CA*CF*(-5D0/2D0)+CF**2*3D0)/4D0
    ! Eq.(3.5) in arXiv:1701.01404
    beta4=CA**5*(8296235D0/3888D0-1630D0/81D0*zeta3+121D0/6D0*zeta4-1045D0/9D0*zeta5)&
         +CA*dAdA*(-514D0/3D0+18716D0/3D0*zeta3-968D0*zeta4-15400D0/3D0*zeta5)&
         +CA**4*nfTF*(-5048959D0/972D0+10505D0/81D0*zeta3-583D0/3D0*zeta4+1230D0*zeta5)&
         +CA**3*CF*nfTF*(8141995D0/1944D0+146D0*zeta3+902D0/3D0*zeta4-8720D0/3D0*zeta5)&
         +CA**2*CF**2*nfTF*(-548732D0/81D0-50581D0/27D0*zeta3-484D0/3D0*zeta4+12820D0/3D0*zeta5)&
         +CA*CF**3*nfTF*(3717D0+5696D0/3D0*zeta3-7480D0/3D0*zeta5)&
         -CF**4*nfTF*(4157D0/6D0+128D0*zeta3)&
         +dAdA*nfTF*(904D0/9D0-20752D0/9D0*zeta3+352D0*zeta4+4000D0/9D0*zeta5)&
         +2D0*dFdA*CA*nfTF*(11312D0/9D0-127736D0/9D0*zeta3+2288D0*zeta4+67520D0/9D0*zeta5)&
         +2D0*dFdA*CF*nfTF*(-320D0+1280D0/3D0*zeta3+6400D0/3D0*zeta5)&
         +CA**3*nfTF**2*(843067D0/486D0+18446D0/27D0*zeta3-104D0/3D0*zeta4-2200D0/3D0*zeta5)&
         +CA**2*CF*nfTF**2*(5701D0/162D0+26452D0/27D0*zeta3-944D0/3D0*zeta4+1600D0/3D0*zeta5)&
         +CF**2*CA*nfTF**2*(31583D0/18D0-28628D0/27D0*zeta3+1144D0/3D0*zeta4-4400D0/3D0*zeta5)&
         +CF**3*nfTF**2*(-5018D0/9D0-2144D0/3D0*zeta3+4640D0/3D0*zeta5)&
         +2D0*dFdA*nfTF**2*(-3680D0/9D0+40160D0/9D0*zeta3-832D0*zeta4-1280D0/9D0*zeta5)&
         +4D0*dFdF*CA*nfTF**2*(-7184D0/3D0+40336D0/9D0*zeta3-704D0*zeta4+2240D0/9D0*zeta5)&
         +4D0*dFdF*CF*nfTF**2*(4160D0/3D0+5120D0/3D0*zeta3-12800D0/3D0*zeta5)&
         +CA**2*nfTF**3*(-2077D0/27D0-9736D0/81D0*zeta3+112D0/3D0*zeta4+320D0/9D0*zeta5)&
         +CA*CF*nfTF**3*(-736D0/81D0-5680D0/27D0*zeta3+224D0/3D0*zeta4)&
         +CF**2*nfTF**3*(-9922D0/81D0+7616D0/27D0*zeta3-352D0/3D0*zeta4)&
         +4D0*dFdF*nfTF**3*(3520D0/9D0-2624D0/3D0*zeta3+256D0*zeta4+1280D0/3D0*zeta5)&
         +CA*nfTF**4*(916D0/243D0-640D0/81D0*zeta3)&
         -CF*nfTF**4*(856D0/243D0+128D0/27D0*zeta3)
    beta4=beta4/32D0
    IF(n.LE.4)RETURN
    WRITE(*,*)"ERROR: beta function is only available with n<=4"
    STOP
    RETURN
  END SUBROUTINE SET_QCD_BETA

  SUBROUTINE SET_QCD_as_decoupling(n,only_nf)
    ! Set the decoupling and inverse decoupling relations for as with the flavour changing
    ! inverse decoupling relation:
    ! alpha_s(mu,nf+1)=alpha_s(mu,nf)*(1+\sum_{m,n}zeta_alphas_inv_X_nm*(alpha_s(mu,nf)/2pi)^n* log^m(mu^2/m_X^2)
    ! decoupling realtion:
    ! alpha_s(mu,nf-1)=alpha_s(mu,nf)*(1+\sum_{m,n}zeta_alphas_X_nm*(alpha_s(mu,nf)/2pi)^n*log^m(mu^2/m_X^2)
    ! where X=MSbar,OS,SI
    ! IN:
    ! n - get zeta_alphas_(inv_)X_1m,...,zeta_alphas_(inv_)X_nm (max n <= 4)
    ! OUT:
    ! zeta_alphas_(inv_)X_1m,..., zeta_alphas_(inv_)X_nm
    USE qcd_constants
    USE global_constants
    IMPLICIT NONE
    INTEGER,INTENT(IN)::n
    REAL(KIND(1d0))::common_term,common_term_inv
    ! try to be efficient
    LOGICAL,INTENT(IN),OPTIONAL::only_nf
    REAL(KIND(1d0))::zeta_alphas_MSbar_32_nf,zeta_alphas_inv_MSbar_32_nonf,zeta_alphas_MSbar_32_nonf
    SAVE zeta_alphas_MSbar_32_nf,zeta_alphas_MSbar_32_nonf,zeta_alphas_inv_MSbar_32_nonf
    REAL(KIND(1d0))::common_term_31_nf,common_term_31_nonf,common_term_inv_31_nonf
    SAVE common_term_31_nf,common_term_31_nonf,common_term_inv_31_nonf
    REAL(KIND(1d0))::zeta_alphas_MSbar_31_noncommon_nf,zeta_alphas_OS_31_noncommon_nf,zeta_alphas_SI_31_noncommon_nf
    SAVE zeta_alphas_MSbar_31_noncommon_nf,zeta_alphas_OS_31_noncommon_nf,zeta_alphas_SI_31_noncommon_nf
    REAL(KIND(1d0))::zeta_alphas_MSbar_31_noncommon_nonf,zeta_alphas_OS_31_noncommon_nonf,zeta_alphas_SI_31_noncommon_nonf
    SAVE zeta_alphas_MSbar_31_noncommon_nonf,zeta_alphas_OS_31_noncommon_nonf,zeta_alphas_SI_31_noncommon_nonf
    REAL(KIND(1d0))::zeta_alphas_inv_MSbar_31_noncommon_nonf,zeta_alphas_inv_OS_31_noncommon_nonf
    REAL(KIND(1d0))::zeta_alphas_inv_SI_31_noncommon_nonf
    SAVE zeta_alphas_inv_MSbar_31_noncommon_nonf,zeta_alphas_inv_OS_31_noncommon_nonf,zeta_alphas_inv_SI_31_noncommon_nonf
    REAL(KIND(1d0))::common_term_30_nf,common_term_30_nonf
    SAVE common_term_30_nf,common_term_30_nonf
    REAL(KIND(1d0))::zeta_alphas_MSbar_30_noncommon_nf,zeta_alphas_OS_30_noncommon_nf
    SAVE zeta_alphas_MSbar_30_noncommon_nf,zeta_alphas_OS_30_noncommon_nf
    REAL(KIND(1d0))::zeta_alphas_MSbar_30_noncommon_nonf,zeta_alphas_OS_30_noncommon_nonf
    SAVE zeta_alphas_MSbar_30_noncommon_nonf,zeta_alphas_OS_30_noncommon_nonf
    REAL(KIND(1d0))::zeta_alphas_MSbar_43_nf2,zeta_alphas_MSbar_43_nf1,zeta_alphas_MSbar_43_nf0
    SAVE zeta_alphas_MSbar_43_nf2,zeta_alphas_MSbar_43_nf1,zeta_alphas_MSbar_43_nf0
    REAL(KIND(1d0))::zeta_alphas_OS_43_nf1,zeta_alphas_OS_43_nf0,zeta_alphas_inv_OS_43_nf0
    SAVE zeta_alphas_OS_43_nf1,zeta_alphas_OS_43_nf0,zeta_alphas_inv_OS_43_nf0
    REAL(KIND(1d0))::zeta_alphas_inv_MSbar_43_nf0,zeta_alphas_inv_MSbar_43_nf1
    SAVE zeta_alphas_inv_MSbar_43_nf0,zeta_alphas_inv_MSbar_43_nf1
    REAL(KIND(1d0))::zeta_alphas_MSbar_42_nf2,zeta_alphas_MSbar_42_nf1,zeta_alphas_MSbar_42_nf0
    SAVE zeta_alphas_MSbar_42_nf2,zeta_alphas_MSbar_42_nf1,zeta_alphas_MSbar_42_nf0
    REAL(KIND(1d0))::zeta_alphas_OS_42_nf2,zeta_alphas_OS_42_nf1,zeta_alphas_OS_42_nf0
    SAVE zeta_alphas_OS_42_nf2,zeta_alphas_OS_42_nf1,zeta_alphas_OS_42_nf0
    REAL(KIND(1d0))::zeta_alphas_SI_42_nf2,zeta_alphas_SI_42_nf1,zeta_alphas_SI_42_nf0
    SAVE zeta_alphas_SI_42_nf2,zeta_alphas_SI_42_nf1,zeta_alphas_SI_42_nf0
    REAL(KIND(1d0))::zeta_alphas_inv_MSbar_42_nf1,zeta_alphas_inv_MSbar_42_nf0
    SAVE zeta_alphas_inv_MSbar_42_nf1,zeta_alphas_inv_MSbar_42_nf0
    REAL(KIND(1d0))::zeta_alphas_inv_OS_42_nf1,zeta_alphas_inv_OS_42_nf0
    SAVE zeta_alphas_inv_OS_42_nf1,zeta_alphas_inv_OS_42_nf0
    REAL(KIND(1d0))::zeta_alphas_inv_SI_42_nf1,zeta_alphas_inv_SI_42_nf0
    SAVE zeta_alphas_inv_SI_42_nf1,zeta_alphas_inv_SI_42_nf0
    REAL(KIND(1d0))::zeta_alphas_MSbar_41_nf2,zeta_alphas_MSbar_41_nf1,zeta_alphas_MSbar_41_nf0
    SAVE zeta_alphas_MSbar_41_nf2,zeta_alphas_MSbar_41_nf1,zeta_alphas_MSbar_41_nf0
    REAL(KIND(1d0))::zeta_alphas_OS_41_nf2,zeta_alphas_OS_41_nf1,zeta_alphas_OS_41_nf0
    SAVE zeta_alphas_OS_41_nf2,zeta_alphas_OS_41_nf1,zeta_alphas_OS_41_nf0
    REAL(KIND(1d0))::zeta_alphas_SI_41_nf2,zeta_alphas_SI_41_nf1,zeta_alphas_SI_41_nf0
    SAVE zeta_alphas_SI_41_nf2,zeta_alphas_SI_41_nf1,zeta_alphas_SI_41_nf0
    REAL(KIND(1d0))::zeta_alphas_inv_MSbar_41_nf1,zeta_alphas_inv_MSbar_41_nf0
    SAVE zeta_alphas_inv_MSbar_41_nf1,zeta_alphas_inv_MSbar_41_nf0
    REAL(KIND(1d0))::zeta_alphas_inv_OS_41_nf1,zeta_alphas_inv_OS_41_nf0
    SAVE zeta_alphas_inv_OS_41_nf1,zeta_alphas_inv_OS_41_nf0
    REAL(KIND(1d0))::zeta_alphas_inv_SI_41_nf1,zeta_alphas_inv_SI_41_nf0
    SAVE zeta_alphas_inv_SI_41_nf1,zeta_alphas_inv_SI_41_nf0
    REAL(KIND(1d0))::zeta_alphas_MSbar_40_nf2,zeta_alphas_MSbar_40_nf1,zeta_alphas_MSbar_40_nf0,zeta_alphas_MSbar_40_nf0_w0
    SAVE zeta_alphas_MSbar_40_nf2,zeta_alphas_MSbar_40_nf1,zeta_alphas_MSbar_40_nf0,zeta_alphas_MSbar_40_nf0_w0
    REAL(KIND(1d0))::zeta_alphas_OS_40_nf2,zeta_alphas_OS_40_nf1,zeta_alphas_OS_40_nf0,zeta_alphas_OS_40_nf0_w0
    SAVE zeta_alphas_OS_40_nf2,zeta_alphas_OS_40_nf1,zeta_alphas_OS_40_nf0,zeta_alphas_OS_40_nf0_w0
    REAL(KIND(1d0))::zeta_alphas_inv_MSbar_40_nf0_w0,zeta_alphas_inv_OS_40_nf0_w0
    SAVE zeta_alphas_inv_MSbar_40_nf0_w0,zeta_alphas_inv_OS_40_nf0_w0
    REAL(KIND(1d0)),PARAMETER::a4=0.5174790616738993863307581618988629456224D0 ! Li4(0.5)
    REAL(KIND(1d0)),PARAMETER::a5=0.5084005792422687074591088492585899413195D0 ! Li5(0.5)
    REAL(KIND(1d0)),PARAMETER::X0=1.808879546208334741426364595086952090D0 ! Eq.(3.2) of hep-ph/0512058
                                                                           ! the analytic X0 is derieved in Eq.(8) of hep-ph/0607202
    IF(PRESENT(only_nf).AND.n.GE.3.AND.only_nf)THEN
       ! update the nfTF terms only
       zeta_alphas_inv_MSbar_32=zeta_alphas_inv_MSbar_32_nonf+zeta_alphas_MSbar_32_nf*nfTF
       zeta_alphas_MSbar_32=zeta_alphas_MSbar_32_nonf-zeta_alphas_MSbar_32_nf*nfTF

       common_term_inv=common_term_31_nf*nfTF+common_term_inv_31_nonf
       common_term=-common_term_31_nf*nfTF+common_term_31_nonf
       zeta_alphas_inv_MSbar_31=common_term_inv+zeta_alphas_MSbar_31_noncommon_nf*nfTF+zeta_alphas_inv_MSbar_31_noncommon_nonf
       zeta_alphas_MSbar_31=common_term-zeta_alphas_MSbar_31_noncommon_nf*nfTF+zeta_alphas_MSbar_31_noncommon_nonf
       zeta_alphas_inv_OS_31=common_term_inv+zeta_alphas_OS_31_noncommon_nf*nfTF+zeta_alphas_inv_OS_31_noncommon_nonf
       zeta_alphas_OS_31=common_term-zeta_alphas_OS_31_noncommon_nf*nfTF+zeta_alphas_OS_31_noncommon_nonf
       zeta_alphas_inv_SI_31=common_term_inv+zeta_alphas_SI_31_noncommon_nf*nfTF+zeta_alphas_inv_SI_31_noncommon_nonf
       zeta_alphas_SI_31=common_term-zeta_alphas_SI_31_noncommon_nf*nfTF+zeta_alphas_SI_31_noncommon_nonf

       common_term_inv=common_term_30_nf*nfTF+common_term_30_nonf
       common_term=-common_term_inv
       zeta_alphas_inv_MSbar_30=common_term_inv+zeta_alphas_MSbar_30_noncommon_nf*nfTF+zeta_alphas_MSbar_30_noncommon_nonf
       zeta_alphas_MSbar_30=-zeta_alphas_inv_MSbar_30
       zeta_alphas_inv_OS_30=common_term_inv+zeta_alphas_OS_30_noncommon_nf*nfTF+zeta_alphas_OS_30_noncommon_nonf
       zeta_alphas_OS_30=-zeta_alphas_inv_OS_30
       zeta_alphas_inv_SI_30=zeta_alphas_inv_MSbar_30
       zeta_alphas_SI_30=-zeta_alphas_inv_SI_30
       IF(n.LE.3)RETURN

       IF(CA.NE.CA_default.OR.CF.NE.CF_default.OR.TF.NE.TF_default)THEN
          WRITE(*,*)"ERROR: 4-loop as decoupling relations in SU(N) N=/=3 are not available"
          STOP
       ENDIF
       ! known 4-loop relations for QCD SU(3)
       zeta_alphas_MSbar_43=zeta_alphas_MSbar_43_nf0+rnf*zeta_alphas_MSbar_43_nf1+rnf**2*zeta_alphas_MSbar_43_nf2
       zeta_alphas_OS_43=zeta_alphas_OS_43_nf0+rnf*zeta_alphas_OS_43_nf1
       zeta_alphas_SI_43=zeta_alphas_OS_43
       zeta_alphas_inv_MSbar_43=zeta_alphas_inv_MSbar_43_nf0&
            +rnf*zeta_alphas_inv_MSbar_43_nf1-rnf**2*zeta_alphas_MSbar_43_nf2
       zeta_alphas_inv_OS_43=zeta_alphas_inv_OS_43_nf0-rnf*zeta_alphas_OS_43_nf1
       zeta_alphas_inv_SI_43=zeta_alphas_inv_OS_43

       zeta_alphas_MSbar_42=zeta_alphas_MSbar_42_nf0+rnf*zeta_alphas_MSbar_42_nf1+rnf**2*zeta_alphas_MSbar_42_nf2
       zeta_alphas_OS_42=zeta_alphas_OS_42_nf0+rnf*zeta_alphas_OS_42_nf1+rnf**2*zeta_alphas_OS_42_nf2
       zeta_alphas_SI_42=zeta_alphas_SI_42_nf0+rnf*zeta_alphas_SI_42_nf1+rnf**2*zeta_alphas_SI_42_nf2
       zeta_alphas_inv_MSbar_42=zeta_alphas_inv_MSbar_42_nf0&
            +rnf*zeta_alphas_inv_MSbar_42_nf1-rnf**2*zeta_alphas_MSbar_42_nf2
       zeta_alphas_inv_OS_42=zeta_alphas_inv_OS_42_nf0+rnf*zeta_alphas_inv_OS_42_nf1-rnf**2*zeta_alphas_OS_42_nf2
       zeta_alphas_inv_SI_42=zeta_alphas_inv_SI_42_nf0+rnf*zeta_alphas_inv_SI_42_nf1-rnf**2*zeta_alphas_SI_42_nf2

       zeta_alphas_MSbar_41=zeta_alphas_MSbar_41_nf0+rnf*zeta_alphas_MSbar_41_nf1+rnf**2*zeta_alphas_MSbar_41_nf2
       zeta_alphas_OS_41=zeta_alphas_OS_41_nf0+rnf*zeta_alphas_OS_41_nf1+rnf**2*zeta_alphas_OS_41_nf2
       zeta_alphas_SI_41=zeta_alphas_SI_41_nf0+rnf*zeta_alphas_SI_41_nf1+rnf**2*zeta_alphas_SI_41_nf2
       zeta_alphas_inv_MSbar_41=zeta_alphas_inv_MSbar_41_nf0&
            +rnf*zeta_alphas_inv_MSbar_41_nf1-rnf**2*zeta_alphas_MSbar_41_nf2
       zeta_alphas_inv_OS_41=zeta_alphas_inv_OS_41_nf0+rnf*zeta_alphas_inv_OS_41_nf1-rnf**2*zeta_alphas_OS_41_nf2
       zeta_alphas_inv_SI_41=zeta_alphas_inv_SI_41_nf0+rnf*zeta_alphas_inv_SI_41_nf1-rnf**2*zeta_alphas_SI_41_nf2

       zeta_alphas_MSbar_40=zeta_alphas_MSbar_40_nf0_w0+zeta_alphas_MSbar_40_nf0&
            +rnf*zeta_alphas_MSbar_40_nf1+rnf**2*zeta_alphas_MSbar_40_nf2
       zeta_alphas_OS_40=zeta_alphas_OS_40_nf0_w0+zeta_alphas_OS_40_nf0+rnf*zeta_alphas_OS_40_nf1+rnf**2*zeta_alphas_OS_40_nf2
       zeta_alphas_SI_40=zeta_alphas_MSbar_40
       zeta_alphas_inv_MSbar_40=zeta_alphas_inv_MSbar_40_nf0_w0-zeta_alphas_MSbar_40_nf0&
            -rnf*zeta_alphas_MSbar_40_nf1-rnf**2*zeta_alphas_MSbar_40_nf2
       zeta_alphas_inv_OS_40=zeta_alphas_inv_OS_40_nf0_w0-zeta_alphas_OS_40_nf0&
            -rnf*zeta_alphas_OS_40_nf1-rnf**2*zeta_alphas_OS_40_nf2
       zeta_alphas_inv_SI_40=zeta_alphas_inv_MSbar_40
       IF(n.LE.4)RETURN

       WRITE(*,*)"ERROR: as decoupling relation is only available with n<=4"
       STOP
       RETURN
    ENDIF
    IF(n.LE.0)RETURN
    zeta_alphas_inv_MSbar_11=2d0/3d0*TF
    zeta_alphas_inv_OS_11=zeta_alphas_inv_MSbar_11
    zeta_alphas_inv_SI_11=zeta_alphas_inv_MSbar_11
    zeta_alphas_MSbar_11=-zeta_alphas_inv_MSbar_11
    zeta_alphas_OS_11=zeta_alphas_MSbar_11
    zeta_alphas_SI_11=zeta_alphas_MSbar_11
    IF(n.LE.1)RETURN
    zeta_alphas_MSbar_22=4d0/9d0*TF**2
    zeta_alphas_OS_22=zeta_alphas_MSbar_22
    zeta_alphas_SI_22=zeta_alphas_MSbar_22
    zeta_alphas_inv_MSbar_22=zeta_alphas_MSbar_22
    zeta_alphas_inv_OS_22=zeta_alphas_MSbar_22
    zeta_alphas_inv_SI_22=zeta_alphas_MSbar_22

    zeta_alphas_inv_MSbar_21=TF*(5d0/3d0*CA-CF)
    zeta_alphas_inv_OS_21=TF*(5d0/3d0*CA+CF)
    zeta_alphas_inv_SI_21=zeta_alphas_inv_OS_21
    zeta_alphas_MSbar_21=-zeta_alphas_inv_MSbar_21
    zeta_alphas_OS_21=-zeta_alphas_inv_OS_21
    zeta_alphas_SI_21=-zeta_alphas_inv_SI_21

    zeta_alphas_inv_MSbar_20=TF*(13d0/12d0*CF-8d0/9d0*CA)
    zeta_alphas_inv_OS_20=TF*(15d0/4d0*CF-8d0/9d0*CA)
    zeta_alphas_inv_SI_20=zeta_alphas_inv_MSbar_20
    zeta_alphas_MSbar_20=-zeta_alphas_inv_MSbar_20
    zeta_alphas_OS_20=-zeta_alphas_inv_OS_20
    zeta_alphas_SI_20=-zeta_alphas_inv_SI_20
    IF(n.LE.2)RETURN
    zeta_alphas_inv_MSbar_33=TF**3*8d0/27d0
    zeta_alphas_inv_OS_33=zeta_alphas_inv_MSbar_33
    zeta_alphas_inv_SI_33=zeta_alphas_inv_MSbar_33
    zeta_alphas_MSbar_33=-zeta_alphas_inv_MSbar_33
    zeta_alphas_OS_33=-zeta_alphas_inv_OS_33
    zeta_alphas_SI_33=-zeta_alphas_inv_SI_33

    zeta_alphas_MSbar_32_nf=TF*2d0/3d0*CF
    zeta_alphas_inv_MSbar_32_nonf=TF*(7d0/12d0*CA**2-11d0/12d0*CA*CF+25d0/9d0*CA*TF-5d0/3d0*CF*TF)
    zeta_alphas_MSbar_32_nonf=TF*(-7d0/12d0*CA**2+11d0/12d0*CA*CF+25d0/9d0*CA*TF-5d0/3d0*CF*TF)
    zeta_alphas_inv_MSbar_32=zeta_alphas_inv_MSbar_32_nonf+zeta_alphas_MSbar_32_nf*nfTF
    zeta_alphas_inv_OS_32=TF*(7d0/12d0*CA**2+11d0/12d0*CA*CF+25d0/9d0*CA*TF+5d0/3d0*CF*TF)
    zeta_alphas_inv_SI_32=zeta_alphas_inv_OS_32
    zeta_alphas_MSbar_32=zeta_alphas_MSbar_32_nonf-zeta_alphas_MSbar_32_nf*nfTF
    zeta_alphas_OS_32=TF*(-7d0/12d0*CA**2-11d0/12d0*CA*CF+25d0/9d0*CA*TF+5d0/3d0*CF*TF)
    zeta_alphas_SI_32=zeta_alphas_OS_32

    common_term_31_nf=TF*CA*(-47d0/54d0)
    common_term_inv_31_nonf=TF*CA*(1063d0/216d0*CA-23d0/12d0*TF)
    common_term_31_nonf=TF*CA*(-1063d0/216d0*CA-113d0/108d0*TF)
    common_term_inv=common_term_31_nf*nfTF+common_term_inv_31_nonf
    common_term=-common_term_31_nf*nfTF+common_term_31_nonf
    zeta_alphas_MSbar_31_noncommon_nf=TF*(-5d0/6d0*CF)
    zeta_alphas_inv_MSbar_31_noncommon_nonf=TF*(-50d0/9d0*CA*CF+9d0/4d0*CF**2+35d0/18d0*CF*TF)
    zeta_alphas_MSbar_31_noncommon_nonf=TF*(50d0/9d0*CA*CF-9d0/4d0*CF**2+5d0/3d0*CF*TF)
    zeta_alphas_inv_MSbar_31=common_term_inv+zeta_alphas_MSbar_31_noncommon_nf*nfTF+zeta_alphas_inv_MSbar_31_noncommon_nonf
    zeta_alphas_MSbar_31=common_term-zeta_alphas_MSbar_31_noncommon_nf*nfTF+zeta_alphas_MSbar_31_noncommon_nonf
    zeta_alphas_OS_31_noncommon_nf=TF*(-67d0/18d0*CF)
    zeta_alphas_inv_OS_31_noncommon_nonf=TF*(175d0/18d0*CA*CF-CF**2/4d0+79d0/18d0*CF*TF)
    zeta_alphas_OS_31_noncommon_nonf=TF*(-175d0/18d0*CA*CF+CF**2/4d0+73d0/9d0*CF*TF)
    zeta_alphas_inv_OS_31=common_term_inv+zeta_alphas_OS_31_noncommon_nf*nfTF+zeta_alphas_inv_OS_31_noncommon_nonf
    zeta_alphas_OS_31=common_term-zeta_alphas_OS_31_noncommon_nf*nfTF+zeta_alphas_OS_31_noncommon_nonf
    zeta_alphas_SI_31_noncommon_nf=TF*(-35d0/18d0*CF)
    zeta_alphas_inv_SI_31_noncommon_nonf=TF*(29d0/6d0*CA*CF-CF**2/4d0+5d0/6d0*CF*TF)
    zeta_alphas_SI_31_noncommon_nonf=TF*(-29d0/6d0*CA*CF+CF**2/4d0+25d0/9d0*CF*TF)
    zeta_alphas_inv_SI_31=common_term_inv+zeta_alphas_SI_31_noncommon_nf*nfTF+zeta_alphas_inv_SI_31_noncommon_nonf
    zeta_alphas_SI_31=common_term-zeta_alphas_SI_31_noncommon_nf*nfTF+zeta_alphas_SI_31_noncommon_nonf

    common_term_30_nf=TF*CA/324d0
    common_term_30_nonf=TF*CA*((5d0/192d0*zeta3-11347d0/2592d0)*CA+(7d0/16d0*zeta3-245d0/648d0)*TF)
    common_term_inv=common_term_30_nf*nfTF+common_term_30_nonf
    common_term=-common_term_inv
    zeta_alphas_MSbar_30_noncommon_nf=TF*164d0/81d0*CF
    zeta_alphas_MSbar_30_noncommon_nonf=TF*((1273d0/96d0*zeta3-2999d0/324d0)*CA*CF+(-95d0/24d0*zeta3+97d0/36d0)*CF**2&
         +(7d0/8d0*zeta3-103d0/162d0)*CF*TF)
    zeta_alphas_inv_MSbar_30=common_term_inv+zeta_alphas_MSbar_30_noncommon_nf*nfTF+zeta_alphas_MSbar_30_noncommon_nonf
    zeta_alphas_MSbar_30=-zeta_alphas_inv_MSbar_30
    zeta_alphas_OS_30_noncommon_nf=TF*(-(8d0/3d0*zeta2+311d0/162d0)*CF)
    zeta_alphas_OS_30_noncommon_nonf=TF*((8d0*LOGTWO*zeta2+1081d0/96d0*zeta3-8d0/3d0*zeta2+8321d0/648d0)*CA*CF&
         +(-16d0*LOGTWO*zeta2+zeta3/24d0+10d0*zeta2+77d0/72d0)*CF**2&
         +(7d0/8d0*zeta3+16d0/3d0*zeta2-695d0/81d0)*CF*TF)
    zeta_alphas_inv_OS_30=common_term_inv+zeta_alphas_OS_30_noncommon_nf*nfTF+zeta_alphas_OS_30_noncommon_nonf
    zeta_alphas_OS_30=-zeta_alphas_inv_OS_30
    zeta_alphas_inv_SI_30=zeta_alphas_inv_MSbar_30
    zeta_alphas_SI_30=-zeta_alphas_inv_SI_30
    IF(n.LE.3)RETURN

    IF(CA.NE.CA_default.OR.CF.NE.CF_default.OR.TF.NE.TF_default)THEN
       WRITE(*,*)"ERROR: 4-loop as decoupling relations in SU(N) N=/=3 are not available"
       STOP
    ENDIF

    ! QCD only, i.e. SU(N) N==3
    zeta_alphas_MSbar_44=1D0/81D0
    zeta_alphas_OS_44=zeta_alphas_MSbar_44
    zeta_alphas_SI_44=zeta_alphas_MSbar_44
    zeta_alphas_inv_MSbar_44=zeta_alphas_MSbar_44
    zeta_alphas_inv_OS_44=zeta_alphas_MSbar_44
    zeta_alphas_inv_SI_44=zeta_alphas_MSbar_44

    zeta_alphas_MSbar_43_nf2=4D0/81D0
    zeta_alphas_MSbar_43_nf1=-127D0/324D0
    zeta_alphas_MSbar_43_nf0=-1883D0/648D0
    zeta_alphas_MSbar_43=zeta_alphas_MSbar_43_nf0+rnf*zeta_alphas_MSbar_43_nf1+rnf**2*zeta_alphas_MSbar_43_nf2
    zeta_alphas_OS_43_nf1=107D0/108D0
    zeta_alphas_OS_43_nf0=-8371D0/648D0
    zeta_alphas_OS_43=zeta_alphas_OS_43_nf0+rnf*zeta_alphas_OS_43_nf1
    zeta_alphas_SI_43=zeta_alphas_OS_43
    zeta_alphas_inv_MSbar_43_nf0=2909D0/648D0
    zeta_alphas_inv_MSbar_43_nf1=271D0/324D0
    zeta_alphas_inv_MSbar_43=zeta_alphas_inv_MSbar_43_nf0&
         +rnf*zeta_alphas_inv_MSbar_43_nf1-rnf**2*zeta_alphas_MSbar_43_nf2
    zeta_alphas_inv_OS_43_nf0=14149D0/648D0
    zeta_alphas_inv_OS_43=zeta_alphas_inv_OS_43_nf0-rnf*zeta_alphas_OS_43_nf1
    zeta_alphas_inv_SI_43=zeta_alphas_inv_OS_43

    zeta_alphas_MSbar_42_nf2=-77D0/1296D0
    zeta_alphas_MSbar_42_nf1=-1483D0/648D0
    zeta_alphas_MSbar_42_nf0=2177D0/216D0
    zeta_alphas_MSbar_42=zeta_alphas_MSbar_42_nf0+rnf*zeta_alphas_MSbar_42_nf1+rnf**2*zeta_alphas_MSbar_42_nf2
    zeta_alphas_OS_42_nf2=-493D0/1296D0
    zeta_alphas_OS_42_nf1=6661D0/648D0
    zeta_alphas_OS_42_nf0=-7693D0/72D0
    zeta_alphas_OS_42=zeta_alphas_OS_42_nf0+rnf*zeta_alphas_OS_42_nf1+rnf**2*zeta_alphas_OS_42_nf2
    zeta_alphas_SI_42_nf2=-79D0/432D0
    zeta_alphas_SI_42_nf1=983D0/216D0
    zeta_alphas_SI_42_nf0=-14023D0/216D0
    zeta_alphas_SI_42=zeta_alphas_SI_42_nf0+rnf*zeta_alphas_SI_42_nf1+rnf**2*zeta_alphas_SI_42_nf2
    zeta_alphas_inv_MSbar_42_nf1=277D0/648D0
    zeta_alphas_inv_MSbar_42_nf0=1837D0/72D0
    zeta_alphas_inv_MSbar_42=zeta_alphas_inv_MSbar_42_nf0&
         +rnf*zeta_alphas_inv_MSbar_42_nf1-rnf**2*zeta_alphas_MSbar_42_nf2
    zeta_alphas_inv_OS_42_nf1=-9115D0/648D0
    zeta_alphas_inv_OS_42_nf0=47039D0/216D0
    zeta_alphas_inv_OS_42=zeta_alphas_inv_OS_42_nf0+rnf*zeta_alphas_inv_OS_42_nf1-rnf**2*zeta_alphas_OS_42_nf2
    zeta_alphas_inv_SI_42_nf1=-515D0/72D0
    zeta_alphas_inv_SI_42_nf0=33887D0/216D0
    zeta_alphas_inv_SI_42=zeta_alphas_inv_SI_42_nf0+rnf*zeta_alphas_inv_SI_42_nf1-rnf**2*zeta_alphas_SI_42_nf2

    zeta_alphas_MSbar_41_nf2=6865D0/11664D0
    zeta_alphas_MSbar_41_nf1=-110341D0/23328D0+110779D0/5184D0*zeta3
    zeta_alphas_MSbar_41_nf0=7391699D0/46656D0-2529743D0/10368D0*zeta3
    zeta_alphas_MSbar_41=zeta_alphas_MSbar_41_nf0+rnf*zeta_alphas_MSbar_41_nf1+rnf**2*zeta_alphas_MSbar_41_nf2
    zeta_alphas_OS_41_nf2=-1679D0/11664D0-16D0/27D0*zeta2
    zeta_alphas_OS_41_nf1=1110443D0/23328D0+328D0/27D0*zeta2+32D0/27D0*zeta2*LOGTWO+132283D0/5184D0*zeta3
    zeta_alphas_OS_41_nf0=-19696909D0/46656D0-464D0/9D0*zeta2-464D0/27D0*zeta2*LOGTWO-2439119D0/10368D0*zeta3
    zeta_alphas_OS_41=zeta_alphas_OS_41_nf0+rnf*zeta_alphas_OS_41_nf1+rnf**2*zeta_alphas_OS_41_nf2
    zeta_alphas_SI_41_nf2=8545D0/11664D0
    zeta_alphas_SI_41_nf1=190283D0/23328D0+133819D0/5184D0*zeta3
    zeta_alphas_SI_41_nf0=-2398621D0/46656D0-2483663D0/10368D0*zeta3
    zeta_alphas_SI_41=zeta_alphas_SI_41_nf0+rnf*zeta_alphas_SI_41_nf1+rnf**2*zeta_alphas_SI_41_nf2
    zeta_alphas_inv_MSbar_41_nf1=141937D0/23328D0-110779D0/5184D0*zeta3
    zeta_alphas_inv_MSbar_41_nf0=-11093717D0/46656D0+3022001D0/10368D0*zeta3
    zeta_alphas_inv_MSbar_41=zeta_alphas_inv_MSbar_41_nf0&
         +rnf*zeta_alphas_inv_MSbar_41_nf1-rnf**2*zeta_alphas_MSbar_41_nf2
    zeta_alphas_inv_OS_41_nf1=-1140191D0/23328D0-376D0/27D0*zeta2-32D0/27D0*zeta2*LOGTWO-132283D0/5184D0*zeta3
    zeta_alphas_inv_OS_41_nf0=21084715D0/46656D0+560D0/9D0*zeta2+560D0/27D0*zeta2*LOGTWO+2922161D0/10368D0*zeta3
    zeta_alphas_inv_OS_41=zeta_alphas_inv_OS_41_nf0+rnf*zeta_alphas_inv_OS_41_nf1-rnf**2*zeta_alphas_OS_41_nf2
    zeta_alphas_inv_SI_41_nf1=-158687D0/23328D0-133819D0/5184D0*zeta3
    zeta_alphas_inv_SI_41_nf0=-1531493D0/46656D0+2975921D0/10368D0*zeta3
    zeta_alphas_inv_SI_41=zeta_alphas_inv_SI_41_nf0+rnf*zeta_alphas_inv_SI_41_nf1-rnf**2*zeta_alphas_SI_41_nf2

    zeta_alphas_MSbar_40_nf2=-271883D0/279936D0+167D0/324D0*zeta3
    zeta_alphas_MSbar_40_nf1=-4770941D0/139968D0-541549D0/10368D0*zeta4+3645913D0/62208D0*zeta3&
         +115D0/36D0*zeta5-685D0/1296D0*LOGTWO**2*zeta2+685D0/7776D0*LOGTWO**4+685D0/324D0*a4
    zeta_alphas_MSbar_40_nf0_w0=291716893D0/382725D0
    zeta_alphas_MSbar_40_nf0=3031309D0/81648D0*LOGTWO**4-121D0/270D0*LOGTWO**5+121D0/27D0*zeta2*LOGTWO**3&
         -3031309D0/13608D0*zeta2*LOGTWO**2+2057D0/36D0*zeta4*LOGTWO+1389D0/16D0*zeta5-76940219D0/136080D0*zeta4&
         -2362581983D0/5443200D0*zeta3+3031309D0/3402D0*a4+484D0/9D0*a5-151369D0/136080D0*X0
    zeta_alphas_MSbar_40=zeta_alphas_MSbar_40_nf0_w0+zeta_alphas_MSbar_40_nf0&
         +rnf*zeta_alphas_MSbar_40_nf1+rnf**2*zeta_alphas_MSbar_40_nf2
    zeta_alphas_OS_40_nf2=-140825D0/93312D0-104D0/81D0*zeta2-19D0/108D0*zeta3
    zeta_alphas_OS_40_nf1=1773073D0/46656D0+173D0/7776D0*LOGTWO**4-1709D0/1296D0*LOGTWO**2*zeta2&
         +352D0/81D0*LOGTWO*zeta2+4756441D0/62208D0*zeta3+4456D0/81D0*zeta2-697709D0/10368D0*zeta4&
         +115D0/36D0*zeta5+173D0/324D0*a4
    zeta_alphas_OS_40_nf0_w0=-141841753D0/1530900D0
    zeta_alphas_OS_40_nf0=-121D0/270D0*LOGTWO**5+3179149D0/81648D0*LOGTWO**4+121D0/27D0*LOGTWO**3*zeta2&
         -2913037D0/13608D0*zeta2*LOGTWO**2+8216D0/81D0*zeta2*LOGTWO+2057D0/36D0*LOGTWO*zeta4+2878D0/27D0*zeta2*zeta3&
         +49309D0/1296D0*zeta5-2408412383D0/5443200D0*zeta3-697121D0/1215D0*zeta2-71102219D0/136080D0*zeta4&
         +3179149D0/3402D0*a4+484D0/9D0*a5-151369D0/136080D0*X0
    zeta_alphas_OS_40=zeta_alphas_OS_40_nf0_w0+zeta_alphas_OS_40_nf0+rnf*zeta_alphas_OS_40_nf1+rnf**2*zeta_alphas_OS_40_nf2
    zeta_alphas_SI_40=zeta_alphas_MSbar_40
    zeta_alphas_inv_MSbar_40_nf0_w0=-1165152397D0/1530900D0
    zeta_alphas_inv_MSbar_40=zeta_alphas_inv_MSbar_40_nf0_w0-zeta_alphas_MSbar_40_nf0&
         -rnf*zeta_alphas_MSbar_40_nf1-rnf**2*zeta_alphas_MSbar_40_nf2
    zeta_alphas_inv_OS_40_nf0_w0=37023232D0/382725D0
    zeta_alphas_inv_OS_40=zeta_alphas_inv_OS_40_nf0_w0-zeta_alphas_OS_40_nf0&
         -rnf*zeta_alphas_OS_40_nf1-rnf**2*zeta_alphas_OS_40_nf2
    zeta_alphas_inv_SI_40=zeta_alphas_inv_MSbar_40

    IF(n.LE.4)RETURN
    WRITE(*,*)"ERROR: as decoupling relation is only available with n<=4"
    STOP
    RETURN
  END SUBROUTINE SET_QCD_as_decoupling

  SUBROUTINE SET_QCD_GAMMA_MASS(n)
    ! Set QCD gamma_m functions up to gamma_{m,n} (n+1 loop) in the MSbar scheme
    ! where gamma_{m,n} appears in 
    ! dm/dlog(mu)=-2m*(gamma_{m,0} as/2/pi+ gamma_{m,1} (as/2/pi)**2 + ...)
    ! IN:
    ! n - get gamma_{m,i}, i=0,...,n (max n<= 3, see Eq.(12) in hep-ph/9703284)
    ! OUT:
    ! gammam0=gamma_{m,0}, gammam1=gamma_{m,1},...
    USE qcd_constants
    USE global_constants
    IMPLICIT NONE
    INTEGER,INTENT(IN)::n
    gammam0=1.5d0*CF
    IF(n.LE.0)RETURN
    gammam1=0.25d0*(1.5d0*CF**2+97d0/6d0*CF*CA-10d0/3d0*nfTF*CF)
    IF(n.LE.1)RETURN
    gammam2=0.125d0*(129d0/2d0*CF**3-129d0/4d0*CF**2*CA+11413d0/108d0*CF*CA**2&
                    +nfTF*(CF**2*(-46d0+48d0*zeta3)+CF*CA*(-556d0/27d0-48d0*zeta3))&
                    -140d0/27d0*nfTF**2*CF)
    IF(n.LE.2)RETURN
    gammam3=1d0/16d0*(CF**4*(-1261d0/8d0-336d0*zeta3)+CF**3*CA*(15349d0/12d0+316d0*zeta3)&
                      +CF**2*CA**2*(-33469d0/36d0-272d0*zeta3+440d0*zeta5)&
                      +CF*CA**3*(69383d0/72d0+2048d0/9d0*zeta3-440d0*zeta5)&
                      +2d0*nfTF*(CF**3*(148d0/3d0-444d0*zeta3-240d0*zeta5)&
                      +CF**2*CA*(-13139d0/54d0+784d0*zeta3-132d0*zeta4+40d0*zeta5)&
                      +CF*CA**2*(-59843d0/324d0-1732d0/3d0*zeta3+132d0*zeta4+200d0*zeta5))&
                      +4d0*nfTF**2*(CF**2*(76d0/27d0-40d0*zeta3+24d0*zeta4)&
                      +CF*CA*(671d0/162d0+40d0*zeta3-24d0*zeta4))&
                      +nfTF**3*CF*(-664d0/81d0+128d0/9d0*zeta3))
    IF(n.LE.3)RETURN
    WRITE(*,*)"ERROR: gamma_m function is only available with n<=3"
    STOP
    RETURN
  END SUBROUTINE SET_QCD_GAMMA_MASS

  SUBROUTINE SET_QCD_mass_decoupling(n,only_nf)
    ! Set the decoupling and inverse decoupling relations for MSbar mass with the flavour changing
    ! inverse decoupling relation:
    ! m(mu,nf+1)=m(mu,nf)*(1+\sum_{m,n}zeta_mass_inv_X_nm*(alpha_s(mu,nf)/2pi)^n* log^m(mu^2/m_X^2)
    ! decoupling realtion:
    ! m(mu,nf-1)=m(mu,nf)*(1+\sum_{m,n}zeta_mass_X_nm*(alpha_s(mu,nf)/2pi)^n*log^m(mu^2/m_X^2)
    ! where X=MSbar,OS,SI
    ! IN:
    ! n - get zeta_mass_(inv_)X_2m,...,zeta_mass_(inv_)X_nm (max n <= 3)
    ! OUT:
    ! zeta_mass_(inv_)X_2m,..., zeta_mass_(inv_)X_nm
    USE qcd_constants
    USE global_constants
    IMPLICIT NONE
    INTEGER,INTENT(IN)::n
    REAL(KIND(1d0))::common_term,common_term_inv
    ! try to be efficient
    LOGICAL,INTENT(IN),OPTIONAL::only_nf
    REAL(KIND(1d0))::zeta_mass_MSbar_33_nf,zeta_mass_inv_MSbar_33_nonf,zeta_mass_MSbar_33_nonf
    SAVE zeta_mass_MSbar_33_nf,zeta_mass_MSbar_33_nonf,zeta_mass_inv_MSbar_33_nonf
    REAL(KIND(1d0))::common_term_31_nf,common_term_31_nonf,common_term_inv_31_nonf
    SAVE common_term_31_nf,common_term_31_nonf,common_term_inv_31_nonf
    REAL(KIND(1d0))::zeta_mass_inv_MSbar_31_noncommon_nonf,zeta_mass_inv_OS_31_noncommon_nonf
    REAL(KIND(1d0))::zeta_mass_inv_SI_31_noncommon_nonf
    SAVE zeta_mass_inv_MSbar_31_noncommon_nonf,zeta_mass_inv_OS_31_noncommon_nonf,zeta_mass_inv_SI_31_noncommon_nonf
    REAL(KIND(1d0))::common_term_30_nf,common_term_30_nonf
    SAVE common_term_30_nf,common_term_30_nonf
    REAL(KIND(1d0))::zeta_mass_inv_MSbar_30_noncommon_nonf,zeta_mass_inv_OS_30_noncommon_nonf
    SAVE zeta_mass_inv_MSbar_30_noncommon_nonf,zeta_mass_inv_OS_30_noncommon_nonf
    IF(PRESENT(only_nf).AND.n.GE.3.AND.only_nf)THEN
       ! update the nfTF terms only
       zeta_mass_inv_MSbar_33=zeta_mass_inv_MSbar_33_nonf-nfTF*zeta_mass_MSbar_33_nf
       zeta_mass_inv_OS_33=zeta_mass_inv_MSbar_33
       zeta_mass_inv_SI_33=zeta_mass_inv_MSbar_33
       zeta_mass_MSbar_33=zeta_mass_MSbar_33_nonf+nfTF*zeta_mass_MSbar_33_nf
       zeta_mass_OS_33=zeta_mass_MSbar_33
       zeta_mass_SI_33=zeta_mass_MSbar_33


       common_term_inv=common_term_31_nf*nfTF+common_term_inv_31_nonf
       common_term=-common_term_31_nf*nfTF+common_term_31_nonf
       zeta_mass_inv_MSbar_31=common_term_inv+zeta_mass_inv_MSbar_31_noncommon_nonf
       zeta_mass_MSbar_31=common_term-zeta_mass_inv_MSbar_31_noncommon_nonf
       zeta_mass_inv_OS_31=common_term_inv+zeta_mass_inv_OS_31_noncommon_nonf
       zeta_mass_OS_31=common_term-zeta_mass_inv_OS_31_noncommon_nonf
       zeta_mass_inv_SI_31=common_term_inv+zeta_mass_inv_SI_31_noncommon_nonf
       zeta_mass_SI_31=common_term-zeta_mass_inv_SI_31_noncommon_nonf

       common_term_inv=common_term_30_nf*nfTF+common_term_30_nonf
       zeta_mass_inv_MSbar_30=common_term_inv+zeta_mass_inv_MSbar_30_noncommon_nonf
       zeta_mass_MSbar_30=-zeta_mass_inv_MSbar_30
       zeta_mass_inv_OS_30=common_term_inv+zeta_mass_inv_OS_30_noncommon_nonf
       zeta_mass_OS_30=-zeta_mass_inv_OS_30
       zeta_mass_inv_SI_30=zeta_mass_inv_MSbar_30
       zeta_mass_SI_30=zeta_mass_MSbar_30
       IF(n.LE.3)RETURN

       WRITE(*,*)"ERROR: mass decoupling relation is only available with n<=3"
       STOP
       RETURN
    ENDIF
    IF(n.LE.1)RETURN
    zeta_mass_MSbar_22=1d0/2d0*TF*CF
    zeta_mass_OS_22=zeta_mass_MSbar_22
    zeta_mass_SI_22=zeta_mass_MSbar_22
    zeta_mass_inv_MSbar_22=-zeta_mass_MSbar_22
    zeta_mass_inv_OS_22=zeta_mass_inv_MSbar_22
    zeta_mass_inv_SI_22=zeta_mass_inv_MSbar_22

    zeta_mass_MSbar_21=-5d0/6d0*CF*TF
    zeta_mass_OS_21=zeta_mass_MSbar_21
    zeta_mass_SI_21=zeta_mass_MSbar_21
    zeta_mass_inv_MSbar_21=-zeta_mass_MSbar_21
    zeta_mass_inv_OS_21=zeta_mass_inv_MSbar_21
    zeta_mass_inv_SI_21=zeta_mass_inv_MSbar_21

    zeta_mass_inv_MSbar_20=-89d0/72d0*CF*TF
    zeta_mass_inv_OS_20=zeta_mass_inv_MSbar_20
    zeta_mass_inv_SI_20=zeta_mass_inv_MSbar_20
    zeta_mass_MSbar_20=-zeta_mass_inv_MSbar_20
    zeta_mass_OS_20=zeta_mass_MSbar_20
    zeta_mass_SI_20=zeta_mass_MSbar_20
    IF(n.LE.2)RETURN
    zeta_mass_MSbar_33_nf=-2d0/9d0*CF*TF
    zeta_mass_inv_MSbar_33_nonf=TF*(-11d0/18d0*CA*CF-2d0/9d0*CF*TF)
    zeta_mass_MSbar_33_nonf=TF*(11d0/18d0*CA*CF-4d0/9d0*CF*TF)
    zeta_mass_inv_MSbar_33=zeta_mass_inv_MSbar_33_nonf-nfTF*zeta_mass_MSbar_33_nf
    zeta_mass_inv_OS_33=zeta_mass_inv_MSbar_33
    zeta_mass_inv_SI_33=zeta_mass_inv_MSbar_33
    zeta_mass_MSbar_33=zeta_mass_MSbar_33_nonf+nfTF*zeta_mass_MSbar_33_nf
    zeta_mass_OS_33=zeta_mass_MSbar_33
    zeta_mass_SI_33=zeta_mass_MSbar_33

    zeta_mass_inv_MSbar_32=TF*CF*(-29d0/12d0*CA+2d0*CF+5d0/9d0*TF)
    zeta_mass_inv_OS_32=TF*CF*(-29d0/12d0*CA-CF+5d0/9d0*TF)
    zeta_mass_inv_SI_32=zeta_mass_inv_OS_32
    zeta_mass_MSbar_32=TF*CF*(29d0/12d0*CA-2d0*CF+5d0/9d0*TF)
    zeta_mass_OS_32=TF*CF*(29d0/12d0*CA+CF+5d0/9d0*TF)
    zeta_mass_SI_32=zeta_mass_OS_32

    common_term_31_nf=TF*CF*(53d0/18d0)
    common_term_inv_31_nonf=-TF*CF*(CA*(5d0/8d0-6d0*zeta3)+CF*6d0*zeta3-35d0/54d0*Tf)
    common_term_31_nonf=TF*CF*(CA*(5d0/8d0-6d0*zeta3)+CF*6d0*zeta3-62d0/27d0*Tf)
    common_term_inv=common_term_31_nf*nfTF+common_term_inv_31_nonf
    common_term=-common_term_31_nf*nfTF+common_term_31_nonf
    zeta_mass_inv_MSbar_31_noncommon_nonf=13d0/8d0*TF*CF**2
    zeta_mass_inv_MSbar_31=common_term_inv+zeta_mass_inv_MSbar_31_noncommon_nonf
    zeta_mass_MSbar_31=common_term-zeta_mass_inv_MSbar_31_noncommon_nonf
    zeta_mass_inv_OS_31_noncommon_nonf=1d0/8d0*TF*CF**2
    zeta_mass_inv_OS_31=common_term_inv+zeta_mass_inv_OS_31_noncommon_nonf
    zeta_mass_OS_31=common_term-zeta_mass_inv_OS_31_noncommon_nonf
    zeta_mass_inv_SI_31_noncommon_nonf=33d0/8d0*TF*CF**2
    zeta_mass_inv_SI_31=common_term_inv+zeta_mass_inv_SI_31_noncommon_nonf
    zeta_mass_SI_31=common_term-zeta_mass_inv_SI_31_noncommon_nonf
    
    common_term_30_nf=TF*CF*(16d0/9d0*zeta3-1327d0/486d0)
    common_term_30_nonf=TF*CF*((-16627d0/1944d0-9d0*zeta4+629d0/72d0*zeta3+B4_CONST)*CA&
         +(-2d0*B4_CONST-57d0/4d0*zeta3+9d0*zeta4)*CF+(1685d0/972d0-28d0/9d0*zeta3)*TF)
    common_term_inv=common_term_30_nf*nfTF+common_term_30_nonf
    zeta_mass_inv_MSbar_30_noncommon_nonf=TF*683d0/72d0*CF**2
    zeta_mass_inv_MSbar_30=common_term_inv+zeta_mass_inv_MSbar_30_noncommon_nonf
    zeta_mass_MSbar_30=-zeta_mass_inv_MSbar_30
    zeta_mass_inv_OS_30_noncommon_nonf=TF*923d0/72d0*CF**2
    zeta_mass_inv_OS_30=common_term_inv+zeta_mass_inv_OS_30_noncommon_nonf
    zeta_mass_OS_30=-zeta_mass_inv_OS_30
    zeta_mass_inv_SI_30=zeta_mass_inv_MSbar_30
    zeta_mass_SI_30=zeta_mass_MSbar_30
    IF(n.LE.3)RETURN
    WRITE(*,*)"ERROR: mass decoupling relation is only available with n<=3"
    STOP
    RETURN
  END SUBROUTINE SET_QCD_mass_decoupling

  SUBROUTINE SET_QCD_GAMMA_CUSP(n)
    ! Set QCD cusp anomalous dimensions up to gamma_{cusp,n} (n loop) in the MSbar scheme
    ! where gamma_{cusp,n} appears in
    ! gamma_{cusp}=gamma_{cusp,1}*as/2/pi+ gamma_{cusp,2} (as/2/pi)**2 + ...
    ! IN:
    ! n - get gamma_{cusp,i}, i=1,...,n (max n<= 4)
    ! OUT:
    ! QCD_gamma1g_cusp=gamma_{cusp,g,1}, QCD_gamma2g_cusp=gamma_{cusp,g,2},...
    ! QCD_gamma1q_cusp=gamma_{cusp,q,1}, QCD_gamma2q_cusp=gamma_{cusp,q,2},...
    USE qcd_constants
    USE global_constants
    IMPLICIT NONE
    INTEGER,INTENT(IN)::n
    REAL(KIND(1d0))::d44FAoNF,d44FFoNF,d44FAoNA,d44AAoNA
    QCD_gamma1q_cusp=2d0*CF
    QCD_gamma1g_cusp=QCD_gamma1q_cusp*CA/CF
    IF(n.LE.1)RETURN
    QCD_gamma2q_cusp=cf*(ca*(67d0/9d0-2d0*zeta2)-20d0/9d0*nfTF)
    QCD_gamma2g_cusp=QCD_gamma2q_cusp*CA/CF
    IF(n.LE.2)RETURN
    QCD_gamma3q_cusp=0.5d0*cf*(-16d0/27d0*nfTF**2+ca*nfTF*(-418d0/27d0+80d0/9d0*zeta2-56d0/3d0*zeta3)&
         +ca**2*(245d0/6d0-268d0/9d0*zeta2+44d0/5d0*zeta2**2+22d0/3d0*zeta3)+cf*nfTF*(-55d0/3d0+16d0*zeta3))
    QCD_gamma3g_cusp=QCD_gamma3q_cusp*CA/CF
    IF(n.LE.3)RETURN
    QCD_gamma4q_cusp=cf/16d0*(nfTF**3*8d0*(64d0/27d0*zeta3-32d0/81d0)&
         +nfTF**2*4D0*ca*(-224d0/15d0*zeta2**2+2240d0/27d0*zeta3-608d0/81d0*zeta2+923d0/81d0)&
         +nfTF**2*4d0*cf*(64d0/5d0*zeta2**2-640d0/9d0*zeta3+2392d0/81d0)&
         +nfTF*2d0*ca**2*(2096d0/9d0*zeta5+448d0/3d0*zeta3*zeta2-352d0/15d0*zeta2**2-23104d0/27d0*zeta3&
         +20320d0/81d0*zeta2-24137d0/81d0)+nfTF*2d0*ca*cf*(160d0*zeta5-128d0*zeta3*zeta2-352d0/5d0*zeta2**2&
         +3712d0/9d0*zeta3+440d0/3d0*zeta2-34066d0/81d0)+nfTF*2d0*cf**2*(-320d0*zeta5+592d0/3d0*zeta3+572d0/9d0)&
         +ca**3*(-16d0*zeta3**2-20032d0/105d0*zeta2**3-3608d0/9d0*zeta5-352d0/3d0*zeta3*zeta2+3608d0/5d0*zeta2**2&
         +20944d0/27d0*zeta3-88400d0/81d0*zeta2+84278d0/81d0))
    QCD_gamma4g_cusp=QCD_gamma4q_cusp*CA/CF
    ! generalised Casimir scaling
    ! d44FF/NF/16
    d44FFoNF=1d0/32d0*cf*(13d0/24d0*ca**2-5d0/2d0*ca*cf+3d0*cf**2)
    ! d44FA/NF/16
    d44FAoNF=1d0/16d0*cf*(7d0/24d0*ca**3-ca**2*cf/2d0)
    ! d44FA/NA/16
    d44FAoNA=1d0/32d0*(7d0/24d0*ca**3-ca**2*cf/2d0)
    ! d44AA/NA/16
    d44AAoNA=1d0/16d0*(37d0/24d0*ca**4-3d0*ca**3*cf)
    QCD_gamma4q_cusp=QCD_gamma4q_cusp+nfTF*2d0*d44FFoNF*(-1280d0/3d0*zeta5-256d0/3d0*zeta3+256d0*zeta2)&
         +d44FAoNF*(-384d0*zeta3**2-7936d0/35d0*zeta2**3+3520d0/3d0*zeta5+128d0/3d0*zeta3-128d0*zeta2)
    QCD_gamma4g_cusp=QCD_gamma4g_cusp+nfTF*2d0*d44FAoNA*(-1280d0/3d0*zeta5-256d0/3d0*zeta3+256d0*zeta2)&
         +d44AAoNA*(-384d0*zeta3**2-7936d0/35d0*zeta2**3+3520d0/3d0*zeta5+128d0/3d0*zeta3-128d0*zeta2)
    IF(n.LE.4)RETURN
    WRITE(*,*)"ERROR: cusp anomalous dimensions are only available with n<=4"
    STOP
    RETURN
  END SUBROUTINE SET_QCD_GAMMA_CUSP
  
  SUBROUTINE SET_QCD_GAMMA_COLLINEAR(n)
    ! Set QCD collinear/normal anomalous dimensions up to gamma_{n} (n loop) in the MSbar scheme
    ! where gamma_{n} appears in
    ! gamma=gamma_{1}*as/2/pi+ gamma_{2} (as/2/pi)**2 + ...
    ! IN:                                                                                                                      
    ! n - get gamma_{i}, i=1,...,n (max n<= 3)
    ! OUT:
    ! QCD_gamma1g_collinear=gamma_{g,1}, QCD_gamma2g_collinear=gamma_{g,2},...
    ! QCD_gamma1q_collinear=gamma_{q,1}, QCD_gamma2q_collinear=gamma_{q,2},...
    USE qcd_constants
    USE global_constants
    IMPLICIT NONE
    INTEGER,INTENT(IN)::n
    QCD_gamma1q_collinear=-3d0/2d0*cf
    QCD_gamma1g_collinear=-11d0/6d0*ca+2d0/3d0*nfTF
    IF(n.LE.1)RETURN
    QCD_gamma2q_collinear=1d0/4d0*(cf**2*(-3d0/2d0+12d0*zeta2-24d0*zeta3)&
         +cf*ca*(-961d0/54d0-11d0*zeta2+26d0*zeta3)+cf*nfTF*(130d0/27d0+4d0*zeta2))
    QCD_gamma2g_collinear=1d0/4d0*(ca**2*(-692d0/27d0+11d0/3d0*zeta2+2d0*zeta3)&
         +ca*nfTF*(256d0/27d0-4d0/3d0*zeta2)+4d0*cf*nfTF)
    IF(n.LE.2)RETURN
    QCD_gamma3q_collinear=1d0/8d0*(cf**3*(-29d0/2d0-18d0*zeta2-288d0/5d0*zeta2**2-68d0*zeta3&
         +32d0*zeta2*zeta3+240d0*zeta5)+cf**2*ca*(-151d0/4d0+410d0/3d0*zeta2+988d0/15d0*zeta2**2&
         -844d0/3d0*zeta3-16d0*zeta2*zeta3-120d0*zeta5)+cf*ca**2*(-139345d0/2916d0-7163d0/81d0*zeta2&
         -166d0/5d0*zeta2**2+3526d0/9d0*zeta3-88d0/3d0*zeta2*zeta3-136d0*zeta5)&
         +cf**2*nfTF*(2953d0/27d0-52d0/3d0*zeta2-112d0/3d0*zeta2**2+512d0/9d0*zeta3)&
         +cf*ca*nfTF*(-17318d0/729d0+5188d0/81d0*zeta2+88d0/5d0*zeta2**2-1928d0/27d0*zeta3)&
         +cf*nfTF**2*(9668d0/729d0-80d0/9d0*zeta2-32d0/27d0*zeta3))
    QCD_gamma3g_collinear=1d0/8d0*(ca**3*(-97186d0/729d0+6109d0/81d0*zeta2-638d0/15d0*zeta2**2&
         +122d0/3d0*zeta3-40d0/3d0*zeta2*zeta3-16d0*zeta5)+ca**2*nfTF*(30715d0/729d0-2396d0/81d0*zeta2&
         +328d0/15d0*zeta2**2+712d0/27d0*zeta3)+ca*cf*nfTF*(2434d0/27d0-4d0*zeta2-32d0/5d0*zeta2**2-304d0/9d0*zeta3)&
         -2d0*cf**2*nfTF+ca*nfTF**2*(-538d0/729d0+80d0/27d0*zeta2-224d0/27d0*zeta3)-44d0/9d0*cf*nfTF**2)
    IF(n.LE.3)RETURN
    WRITE(*,*)"ERROR: collinear/normal anomalous dimensions are only available with n<=3"
    STOP
    RETURN
  END SUBROUTINE SET_QCD_GAMMA_COLLINEAR
END MODULE qcd_setup
