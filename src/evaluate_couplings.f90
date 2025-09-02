MODULE evaluate_couplings
  IMPLICIT NONE
CONTAINS
    FUNCTION QMASS_RG_RUNNING(nloop,as0,asmu,mass0) RESULT(res)
    ! the RG running of MSbar quark mass
    ! it has the relation wrt the MSbar Yukawa coupling via 
    ! y(mu)=sqrt(2)*m(mu)/vev, vev=246.2 GeV at scale mu
    ! It will not vary the number of quark flavours in RG running
    USE qcd_constants
    USE global_constants
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nloop ! from 0 to 3
    REAL(KIND(1d0)),INTENT(IN)::as0,asmu
    REAL(KIND(1d0)),INTENT(IN)::mass0 ! scale invariant mass, i.e. 
                                      ! m^{MSar}(m^{SI})=m^{SI} can be read from PDG
    REAL(KIND(1d0))::res
    REAL(KIND(1d0))::c0,c1,c2,c3
    SAVE c0,c1,c2,c3
    INTEGER::init=0
    SAVE init
    REAL(KIND(1d0))::casmu,cas0
    REAL(KIND(1d0))::asmuo2pi,as0o2pi,fac
    IF(nloop.LT.0)THEN
       res=mass0
       RETURN
    ENDIF
    IF(init.EQ.0)THEN
       IF(beta0.EQ.0d0.OR.gammam0.EQ.0d0)THEN
          WRITE(*,*)"Error: Please update beta or gamma_m functions first #0"
          STOP
       ENDIF
       IF(nloop.GE.1.AND.(beta1.EQ.0.OR.gammam1.EQ.0d0))THEN
          WRITE(*,*)"Error: Please update beta or gamma_m functions first #1"
          STOP
       ENDIF
       IF(nloop.GE.2.AND.(beta2.EQ.0.OR.gammam2.EQ.0d0))THEN
          WRITE(*,*)"Error: Please update beta or gamma_m functions first #2"
          STOP
       ENDIF
       IF(nloop.GE.3.AND.(beta3.EQ.0.OR.gammam3.EQ.0d0))THEN
          WRITE(*,*)"Error: Please update beta or gamma_m functions first #3"
          STOP
       ENDIF
       c0=gammam0/beta0
       c1=-beta1*gammam0/beta0**2+gammam1/beta0
       c2=gammam0/beta0**2*(beta1**2/beta0-beta2)-beta1*gammam1/beta0**2+gammam2/beta0
       c3=gammam0/beta0**2*(TWO*beta1*beta2/beta0-beta1**3/beta0**2-beta3)&
            +gammam1/beta0**2*(beta1**2/beta0-beta2)-beta1*gammam2/beta0**2+gammam3/beta0
       init=1
    ENDIF
    as0o2pi=as0/TWOPI
    asmuo2pi=asmu/TWOPI
    cas0=one
    casmu=one
    IF(nloop.LE.0)THEN
       cas0=cas0*as0o2pi**c0
       casmu=casmu*asmuo2pi**c0
       res=casmu/cas0*mass0
       RETURN
    ENDIF
    cas0=cas0+c1*as0o2pi
    casmu=casmu+c1*asmuo2pi
    IF(nloop.LE.1)THEN
       cas0=cas0*as0o2pi**c0
       casmu=casmu*asmuo2pi**c0
       res=casmu/cas0*mass0
       RETURN
    ENDIF
    fac=HALF*(c1**2+c2)
    cas0=cas0+fac*as0o2pi**2
    casmu=casmu+fac*asmuo2pi**2
    IF(nloop.LE.2)THEN
       cas0=cas0*as0o2pi**c0
       casmu=casmu*asmuo2pi**c0
       res=casmu/cas0*mass0
       RETURN
    ENDIF
    fac=THIRD*(HALF*c1**3+THREE*HALF*c1*c2+c3)
    cas0=cas0+fac*as0o2pi**3
    casmu=casmu+fac*asmuo2pi**3
    IF(nloop.LE.3)THEN
       cas0=cas0*as0o2pi**c0
       casmu=casmu*asmuo2pi**c0
       res=casmu/cas0*mass0
       RETURN
    ENDIF
    WRITE(*,*)"ERROR: cannot run quark mass with nloop = ", nloop
    STOP
  END FUNCTION QMASS_RG_RUNNING

  FUNCTION ALPHAS(Q)
    USE qcd_coupling
    USE coupling_global
    IMPLICIT NONE
    REAL(KIND(1d0))::ALPHAS
    REAL(KIND(1d0)),INTENT(IN)::Q
    ALPHAS=as_value(as_box,Q)
  END FUNCTION ALPHAS

  FUNCTION ALPHAEW(Q)
    USE global
    IMPLICIT NONE
    REAL(KIND(1d0))::ALPHAEW
    REAL(KIND(1d0)),INTENT(IN)::Q
    REAL(KIND(1d0)),SAVE::ALPHAEW_SAVE
    REAL(KIND(1d0)),SAVE::DalphalepMZ=0d0,DalphahadMZ=0d0
    REAL(KIND(1d0))::DalphalepQ,DDalphahad
    INTEGER,SAVE::init=0
    REAL(KIND(1d0)),PARAMETER::pipi=3.14159265358979323846264338328d0
    IF(init.EQ.0)THEN
       IF(alpha_scheme.EQ.0)THEN
          ! alpha(0) scheme
          ALPHAEW_SAVE=1d0/alphaemm1
       ELSEIF(alpha_scheme.EQ.1)THEN
          ! Gmu scheme
          ALPHAEW_SAVE=DSQRT(2d0)*Gfermi*wmass_PDG**2&
               *(zmass_PDG**2-wmass_PDG**2)/zmass_PDG**2/pipi
       ELSE
          ALPHAEW_SAVE=0d0
          DalphalepMZ=DeltaAlphaLep(zmass_PDG)
          DalphahadMZ=DeltaAlphaHad(zmass_PDG)
       ENDIF
       init=1
    ENDIF
    IF(alpha_scheme.LE.1)THEN
       ! alpha(0) or Gmu scheme
       ! we use the saved value
       ALPHAEW=ALPHAEW_SAVE
    ELSE
       ! alpha(Q) scheme
       ! alpha(Q)=alpha(0)/(1-Delta alpha(Q))
       ! Delta alpha(Q)=Delta alpha_{had}(Q)+Delta alpha_{lep}(Q)
       ! Delta alpha(Q)=Sigma_{AA}(0)-Re(Sigma_{AA}(Q))
       ! alpha(Q)=alpha(MZ)*(1-Delta alpha_{had}(MZ)-Delta alpha_{lep}(MZ))/
       !    (1-Delta alpha_{had}(MZ)-Delta alpha_{lep}(Q)+(Delta alpha_{had}(MZ)-Delta alpha_{had}(Q))),
       ! where Delta alpha_{had}(MZ)-Delta alpha_{had}(Q) we use the perturbative one-loop results [insensitive to light quark mass when Q > 1 GeV]
       DalphalepQ=DeltaAlphaLep(Q)
       DDalphahad=DalphahadMZ-DeltaAlphaHad(Q)
       ALPHAEW=1d0/alphaMZm1*(1d0-DalphahadMZ_PDG-DalphalepMZ)&
            /(1d0-DalphahadMZ_PDG-DalphalepQ+DDalphahad)
    ENDIF
  END FUNCTION ALPHAEW

  FUNCTION ReSIGMAAA(Q,M)
    ! one-loop perturbation result of
    ! ReSigma_AA(Q)-Sigma_AA(0) for a given fermion
    ! It has been dropped alpha/(3*Pi)*Ncf*Qf**2
    IMPLICIT NONE
    REAL(KIND(1d0))::ReSIGMAAA
    REAL(KIND(1d0)),INTENT(IN)::Q,M
    COMPLEX(KIND(1d0))::SIGMAAA
    REAL(KIND(1d0))::MQ2Q2
    COMPLEX(KIND(1d0))::sqrtt
    MQ2Q2=M**2/Q**2
    sqrtt=SQRT(DCMPLX(1d0-4d0*MQ2Q2,0d0))
    SIGMAAA=5d0/3d0+4d0*MQ2Q2+sqrtt*(1d0+2d0*MQ2Q2)*LOG((sqrtt-1d0)/(sqrtt+1d0))
    ReSIGMAAA=DREAL(SIGMAAA)
    RETURN
  END FUNCTION ReSIGMAAA

  FUNCTION DeltaAlphaLep(Q)
    USE global
    ! one-loop
    ! Delta alpha_{lep}(Q)
    IMPLICIT NONE
    REAL(KIND(1d0))::DeltaAlphaLep
    REAL(KIND(1d0)),INTENT(IN)::Q
    REAL(KIND(1d0)),PARAMETER::pipi=3.14159265358979323846264338328d0
    DeltaAlphaLep=ReSIGMAAA(Q,emass_PDG)
    DeltaAlphaLep=DeltaAlphaLep+ReSIGMAAA(Q,mumass_PDG)
    DeltaAlphaLep=DeltaAlphaLep+ReSIGMAAA(Q,taumass_PDG)
    DeltaAlphaLep=DeltaAlphaLep*(-1d0/(3d0*pipi*alphaemm1))
    RETURN
  END FUNCTION DeltaAlphaLep

  FUNCTION DeltaAlphaHad(Q)
    use global
    ! one-loop
    ! Delta alpha_{had}(Q) [only 5 quark flavours]
    IMPLICIT NONE
    REAL(KIND(1d0))::DeltaAlphaHad
    REAL(KIND(1d0)),INTENT(IN)::Q
    REAL(KIND(1d0)),PARAMETER::pipi=3.14159265358979323846264338328d0
    DeltaAlphaHad=ReSIGMAAA(Q,umass_PDG)*(2d0/3d0)**2
    DeltaAlphaHad=DeltaAlphaHad+ReSIGMAAA(Q,dmass_PDG)*(-1d0/3d0)**2
    DeltaAlphaHad=DeltaAlphaHad+ReSIGMAAA(Q,smass_PDG)*(-1d0/3d0)**2
    DeltaAlphaHad=DeltaAlphaHad+ReSIGMAAA(Q,cmass_PDG)*(2d0/3d0)**2
    DeltaAlphaHad=DeltaAlphaHad+ReSIGMAAA(Q,bmass_PDG)*(-1d0/3d0)**2
    DeltaAlphaHad=DeltaAlphaHad*(-1d0/(pipi*alphaemm1))
    RETURN
  END FUNCTION DeltaAlphaHad
END MODULE evaluate_couplings
