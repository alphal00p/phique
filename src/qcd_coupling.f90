MODULE qcd_coupling
  USE assertions
  IMPLICIT NONE
  PRIVATE

  TYPE as_segment
     ! it respresents the segmentation for a given number of flavours
     PRIVATE
     REAL(KIND(1d0))::t_low,t_high,dt
     REAL(KIND(1d0)),DIMENSION(0:4):: beta
     REAL(KIND(1d0))::zeta_as_11, zeta_as_22, zeta_as_21, zeta_as_20
     REAL(KIND(1d0))::zeta_as_33, zeta_as_32, zeta_as_31, zeta_as_30
     REAL(KIND(1d0))::zeta_as_44, zeta_as_43, zeta_as_42, zeta_as_41, zeta_as_40
     REAL(KIND(1d0))::zeta_as_inv_11, zeta_as_inv_22, zeta_as_inv_21, zeta_as_inv_20
     REAL(KIND(1d0))::zeta_as_inv_33, zeta_as_inv_32, zeta_as_inv_31, zeta_as_inv_30
     REAL(KIND(1d0))::zeta_as_inv_44, zeta_as_inv_43, zeta_as_inv_42, zeta_as_inv_41, zeta_as_inv_40
     REAL(KIND(1d0)),DIMENSION(:),POINTER:: as_m1 ! as^(-1)
     INTEGER:: iflip ! above iflip -> upward evolution, below iflip -> downward evolution
  END TYPE as_segment

  TYPE as_container
     PRIVATE
     TYPE(as_segment), DIMENSION(:), POINTER :: seg
     ! n_low: the lowest number of nf
     ! n_high: the highest number of nf
     INTEGER:: n_low, n_high
     ! n_loop: number of loop for beta function (i.e. beta_(n_loop-1))
     ! fix_nf: fix the number of flavour, otherwise it is simply no_fix_nf below
     INTEGER:: n_loop, fix_nf
     ! alpha_s: initial value of alpha_s at initial scale Q
     ! Q: initial scale Q
     REAL(KIND(1d0))::alpha_s, Q
     ! quark_masses: quark masses (to determin nf)
     REAL(KIND(1d0)),DIMENSION(6)::quark_masses
     ! muMatch_mQuark: the flavour threshold is muMatch_mQuark*quark_masses(i)
     REAL(KIND(1d0))::muMatch_mQuark
     INTEGER:: masses_type ! 1: MSbar mass; 2: OS mass; 3: SI mass
  END TYPE as_container

  TYPE as_nf_interval
     PRIVATE
     REAL(KIND(1d0))::Q_low,Q_high
     REAL(KIND(1d0))::zeta_as_11, zeta_as_22, zeta_as_21, zeta_as_20
     REAL(KIND(1d0))::zeta_as_33, zeta_as_32, zeta_as_31, zeta_as_30
     REAL(KIND(1d0))::zeta_as_44, zeta_as_43, zeta_as_42, zeta_as_41, zeta_as_40
     REAL(KIND(1d0))::zeta_as_inv_11, zeta_as_inv_22, zeta_as_inv_21, zeta_as_inv_20
     REAL(KIND(1d0))::zeta_as_inv_33, zeta_as_inv_32, zeta_as_inv_31, zeta_as_inv_30
     REAL(KIND(1d0))::zeta_as_inv_44, zeta_as_inv_43, zeta_as_inv_42, zeta_as_inv_41, zeta_as_inv_40
     REAL(KIND(1d0)),DIMENSION(5)::as_at_Q_low,as_at_Q_high ! 1-5 loops of as at Q_low and Q_high
  END TYPE as_nf_interval

  TYPE as_framework
     PRIVATE
     TYPE(as_nf_interval), DIMENSION(:), POINTER :: seg
     ! n_low: the lowest number of nf
     ! n_high: the highest number of nf
     INTEGER:: n_low, n_high
     ! alpha_s: initial value of alpha_s at initial scale Q
     ! Q: initial scale Q
     REAL(KIND(1d0))::alpha_s,Q
     ! quark_masses: quark masses (to determin nf)
     REAL(KIND(1d0)),DIMENSION(6)::quark_masses
     ! muMatch_mQuark: the flavour threshold is muMatch_mQuark*quark_masses(i)
     REAL(KIND(1d0))::muMatch_mQuark
     INTEGER:: masses_type ! 1: MSbar mass; 2: OS mass; 3: SI mass
     REAL(KIND(1d0)),DIMENSION(6),PUBLIC::Q_THR ! threshold for heavy quarks
  END TYPE as_framework

  ! t=Log(Q^2)
  REAL(KIND(1d0))::dt_used = 0.2d0
  REAL(KIND(1d0)),PARAMETER::t_low=-1.3862943611198906d0 ! Q_low=0.5 GeV
  REAL(KIND(1d0)),PARAMETER::t_high=92.1034037197618273607196581873745d0 ! Q_high=10^20 GeV
  INTEGER,PARAMETER::no_fix_nf=-1000000045 ! fix_nf == no_fix_nf means we used the varied flavour number scheme

  TYPE(as_segment),POINTER:: seg ! put it here necessary to use the runge kutta method as it will be used in dasm1_dt (the scalar convolution function)

  PUBLIC::no_fix_nf
  PUBLIC::as_container
  PUBLIC::as_run_init,as_value
  PUBLIC::SET_AS_DT,GET_AS_DT
  PUBLIC::as_clean,Get_as_QuarkMass,GET_as_NumberOfLoop
  PUBLIC::Get_as_QuarkMassType,Get_as_nf_range
  PUBLIC::Get_as_nfAtQ,Get_as_QRangeAtnf
  ! a few special NLO ALPHAS running routines
  ! works only up to NLO (2-loop) running
  PUBLIC::as_framework
  PUBLIC::ALPHAS_FROM_LAMBDA,LAMBDA_FROM_ALPHAS
  PUBLIC::ALPHAS_LAMBDA_INIT

CONTAINS

!=============================================================================!
!=                                                                           =!
!=  BEGINING OF GENERAL alpha_s running with runge kutta method              =!
!=                                                                           =!
!=============================================================================!
   
  SUBROUTINE as_run_init(as_box,alpha_s,Q,n_loop,fix_nf,quark_masses,masses_type,muMatch_mQuark)
    USE global_constants
    USE qcd_constants
    USE qcd_setup
    IMPLICIT NONE
    TYPE(as_container),INTENT(OUT),TARGET::as_box
    REAL(KIND(1d0)),INTENT(IN),OPTIONAL::alpha_s,Q
    INTEGER,INTENT(IN),OPTIONAL::n_loop,fix_nf
    REAL(KIND(1d0)),INTENT(IN),DIMENSION(4:6),OPTIONAL::quark_masses
    REAL(KIND(1d0)),INTENT(IN),OPTIONAL::muMatch_mQuark
    INTEGER,INTENT(IN),OPTIONAL::masses_type
    !  local variables
    INTEGER::i,j,nbin_total,nbin,nseg,istart
    REAL(KIND(1d0))::tstart,t,as_m1,as_start,aso2pi_start,log_match
    TYPE(as_segment),POINTER:: seg2
    INTEGER::nf_orig

    nf_orig = nf

    as_box%alpha_s = default_or_opt(0.118d0, alpha_s) ! initial value of alpha_s
    as_box%Q = default_or_opt(91.2d0, Q) ! initial value of scale
    as_box%n_loop = default_or_opt(2, n_loop) ! running via loop of beta function
    as_box%fix_nf = default_or_opt(no_fix_nf,fix_nf) ! whether need to fix the number of flavour
    as_box%muMatch_mQuark = default_or_opt(ONE, muMatch_mQuark) ! the threshold is muMatch_mQuark*quark_mass
    as_box%masses_type = default_or_opt(2, masses_type) ! 1: MSbar mass; 2: OS; 3: SI
    as_box%quark_masses(1:3) = quark_masses_default(1:3) ! quark masses
    as_box%quark_masses(4:6) = default_or_opt(quark_masses_default(4:6),quark_masses)
    
    DO i=LBOUND(as_box%quark_masses, dim=1), UBOUND(as_box%quark_masses, dim=1)
       IF(t_low.GT.tOfQ(as_box%muMatch_mQuark*as_box%quark_masses(i))) as_box%n_low = i
       IF(t_high.GT.tOfQ(as_box%muMatch_mQuark*as_box%quark_masses(i))) as_box%n_high = i
    ENDDO

    ! now set up the parameters of the segments
    ALLOCATE(as_box%seg(as_box%n_low:as_box%n_high))
    nbin_total = 0
    DO i = as_box%n_low, as_box%n_high
       seg => as_box%seg(i)
       seg%t_low = MAX(t_low, tOfQ(as_box%muMatch_mQuark*as_box%quark_masses(i)))
       IF(i.LT.as_box%n_high ) THEN
          seg%t_high = tOfQ(as_box%muMatch_mQuark*as_box%quark_masses(i+1))
       ELSE
          seg%t_high = t_high
       ENDIF

       ! now check the number of flavours
       IF(as_box%fix_nf.EQ.no_fix_nf)THEN
          CALL SET_QCD_NF(i)
       ELSE
          CALL SET_QCD_NF(as_box%fix_nf)
       ENDIF
       
       seg%beta(0)=beta0
       SELECT CASE(as_box%n_loop)
       CASE(1)
          seg%beta(1)=0d0
          seg%beta(2)=0d0
          seg%beta(3)=0d0
          seg%beta(4)=0d0
       CASE(2)
          seg%beta(1)=beta1
          seg%beta(2)=0d0
          seg%beta(3)=0d0
          seg%beta(4)=0d0
       CASE(3)
          seg%beta(1)=beta1
          seg%beta(2)=beta2
          seg%beta(3)=0d0
          seg%beta(4)=0d0
       CASE(4)
          seg%beta(1)=beta1
          seg%beta(2)=beta2
          seg%beta(3)=beta3
          seg%beta(4)=0d0
       CASE(5)
          seg%beta(1)=beta1
          seg%beta(2)=beta2
          seg%beta(3)=beta3
          seg%beta(4)=beta4
       CASE default
          PRINT *, "ERROR: only 5-loop beta function is supported in as_run_init"
          STOP
       END SELECT

       ! the decoupling relation when changing nf
       SELECT CASE(as_box%masses_type)
       CASE(1)
          ! quark masses are MSbar masses
          seg%zeta_as_11=zeta_alphas_MSbar_11
          seg%zeta_as_22=zeta_alphas_MSbar_22
          seg%zeta_as_21=zeta_alphas_MSbar_21
          seg%zeta_as_20=zeta_alphas_MSbar_20
          seg%zeta_as_33=zeta_alphas_MSbar_33
          seg%zeta_as_32=zeta_alphas_MSbar_32
          seg%zeta_as_31=zeta_alphas_MSbar_31
          seg%zeta_as_30=zeta_alphas_MSbar_30
          seg%zeta_as_44=zeta_alphas_MSbar_44
          seg%zeta_as_43=zeta_alphas_MSbar_43
          seg%zeta_as_42=zeta_alphas_MSbar_42
          seg%zeta_as_41=zeta_alphas_MSbar_41
          seg%zeta_as_40=zeta_alphas_MSbar_40

          seg%zeta_as_inv_11=zeta_alphas_inv_MSbar_11
          seg%zeta_as_inv_22=zeta_alphas_inv_MSbar_22
          seg%zeta_as_inv_21=zeta_alphas_inv_MSbar_21
          seg%zeta_as_inv_20=zeta_alphas_inv_MSbar_20
          seg%zeta_as_inv_33=zeta_alphas_inv_MSbar_33
          seg%zeta_as_inv_32=zeta_alphas_inv_MSbar_32
          seg%zeta_as_inv_31=zeta_alphas_inv_MSbar_31
          seg%zeta_as_inv_30=zeta_alphas_inv_MSbar_30
          seg%zeta_as_inv_44=zeta_alphas_inv_MSbar_44
          seg%zeta_as_inv_43=zeta_alphas_inv_MSbar_43
          seg%zeta_as_inv_42=zeta_alphas_inv_MSbar_42
          seg%zeta_as_inv_41=zeta_alphas_inv_MSbar_41
          seg%zeta_as_inv_40=zeta_alphas_inv_MSbar_40
       CASE(2)
          ! quark masses are OS masses
          seg%zeta_as_11=zeta_alphas_OS_11
          seg%zeta_as_22=zeta_alphas_OS_22
          seg%zeta_as_21=zeta_alphas_OS_21
          seg%zeta_as_20=zeta_alphas_OS_20
          seg%zeta_as_33=zeta_alphas_OS_33
          seg%zeta_as_32=zeta_alphas_OS_32
          seg%zeta_as_31=zeta_alphas_OS_31
          seg%zeta_as_30=zeta_alphas_OS_30
          seg%zeta_as_44=zeta_alphas_OS_44
          seg%zeta_as_43=zeta_alphas_OS_43
          seg%zeta_as_42=zeta_alphas_OS_42
          seg%zeta_as_41=zeta_alphas_OS_41
          seg%zeta_as_40=zeta_alphas_OS_40

          seg%zeta_as_inv_11=zeta_alphas_inv_OS_11
          seg%zeta_as_inv_22=zeta_alphas_inv_OS_22
          seg%zeta_as_inv_21=zeta_alphas_inv_OS_21
          seg%zeta_as_inv_20=zeta_alphas_inv_OS_20
          seg%zeta_as_inv_33=zeta_alphas_inv_OS_33
          seg%zeta_as_inv_32=zeta_alphas_inv_OS_32
          seg%zeta_as_inv_31=zeta_alphas_inv_OS_31
          seg%zeta_as_inv_30=zeta_alphas_inv_OS_30
          seg%zeta_as_inv_44=zeta_alphas_inv_OS_44
          seg%zeta_as_inv_43=zeta_alphas_inv_OS_43
          seg%zeta_as_inv_42=zeta_alphas_inv_OS_42
          seg%zeta_as_inv_41=zeta_alphas_inv_OS_41
          seg%zeta_as_inv_40=zeta_alphas_inv_OS_40
       CASE(3)
          ! quark masses are SI masses
          seg%zeta_as_11=zeta_alphas_SI_11
          seg%zeta_as_22=zeta_alphas_SI_22
          seg%zeta_as_21=zeta_alphas_SI_21
          seg%zeta_as_20=zeta_alphas_SI_20
          seg%zeta_as_33=zeta_alphas_SI_33
          seg%zeta_as_32=zeta_alphas_SI_32
          seg%zeta_as_31=zeta_alphas_SI_31
          seg%zeta_as_30=zeta_alphas_SI_30
          seg%zeta_as_44=zeta_alphas_SI_44
          seg%zeta_as_43=zeta_alphas_SI_43
          seg%zeta_as_42=zeta_alphas_SI_42
          seg%zeta_as_41=zeta_alphas_SI_41
          seg%zeta_as_40=zeta_alphas_SI_40

          seg%zeta_as_inv_11=zeta_alphas_inv_SI_11
          seg%zeta_as_inv_22=zeta_alphas_inv_SI_22
          seg%zeta_as_inv_21=zeta_alphas_inv_SI_21
          seg%zeta_as_inv_20=zeta_alphas_inv_SI_20
          seg%zeta_as_inv_33=zeta_alphas_inv_SI_33
          seg%zeta_as_inv_32=zeta_alphas_inv_SI_32
          seg%zeta_as_inv_31=zeta_alphas_inv_SI_31
          seg%zeta_as_inv_30=zeta_alphas_inv_SI_30
          seg%zeta_as_inv_44=zeta_alphas_inv_SI_44
          seg%zeta_as_inv_43=zeta_alphas_inv_SI_43
          seg%zeta_as_inv_42=zeta_alphas_inv_SI_42
          seg%zeta_as_inv_41=zeta_alphas_inv_SI_41
          seg%zeta_as_inv_40=zeta_alphas_inv_SI_40
       CASE DEFAULT
          PRINT *, "ERROR: masses_type can only be MSbar(1)/OS(2)/SI(3) masses"
          STOP
       END SELECT
    
       seg%dt = dt_used*0.5d0**(6-i)
       nbin = CEILING((seg%t_high-seg%t_low)/seg%dt)
       nbin_total = nbin_total+nbin
       seg%dt = (seg%t_high-seg%t_low)/DBLE(nbin)
       ALLOCATE(seg%as_m1(0:nbin))
    ENDDO
    
    tstart = tOfQ(as_box%Q)
    DO i=as_box%n_low, as_box%n_high
       IF(tstart.LE.as_box%seg(i)%t_high.AND.tstart.GE.as_box%seg(i)%t_low)EXIT
    ENDDO
    IF(i.GT.as_box%n_high)THEN
       PRINT *, "ERROR: the initial scale Q is not in the supported range"
       STOP
    ENDIF
    nseg=i
    
    ! we start from i-th segmentation as the initial scale Q lies in
    as_m1 = 1d0/as_box%alpha_s
    seg => as_box%seg(nseg)
    istart=NINT((tstart-seg%t_low)/seg%dt)
    t = tstart
    seg%as_m1(istart) = as_run_onestep(as_m1, seg%t_low+istart*seg%dt-tstart)
    DO i=istart+1, UBOUND(seg%as_m1, dim=1)
       seg%as_m1(i) = as_run_onestep(seg%as_m1(i-1), seg%dt)
    ENDDO
    DO i=istart-1, LBOUND(seg%as_m1, dim=1), -1
       seg%as_m1(i) = as_run_onestep(seg%as_m1(i+1), -seg%dt)
    ENDDO
    seg%iflip = istart

    ! for other segments
    log_match=2d0*DLOG(as_box%muMatch_mQuark) ! the ratio of log(mu^2/quark_mass^2) with the matching mu=muMatch_mQuark*quark_mass
    ! the above segments first
    DO j=nseg+1, as_box%n_high
       seg => as_box%seg(j)
       seg2 => as_box%seg(j-1)
       ! nf -> nf+1 (inverse decoupling relation)
       ! the as^(nf) at threshold muMatch_mQuark*m^(nf+1)
       as_start = 1d0/as_box%seg(j-1)%as_m1(UBOUND(as_box%seg(j-1)%as_m1,dim=1))
       aso2pi_start = as_start/TWOPI
       IF(as_decoupling_on.AND.as_box%fix_nf.EQ.no_fix_nf)THEN
          SELECT CASE(as_box%n_loop)
          CASE(2)
             as_start=as_start*(1d0+seg2%zeta_as_inv_11*log_match*aso2pi_start)
          CASE(3)
             as_start=as_start*(1d0+seg2%zeta_as_inv_11*log_match*aso2pi_start&
                  +(seg2%zeta_as_inv_22*log_match**2+seg2%zeta_as_inv_21*log_match&
                  +seg2%zeta_as_inv_20)*aso2pi_start**2)
          CASE(4)
             as_start=as_start*(1d0+seg2%zeta_as_inv_11*log_match*aso2pi_start&
                  +(seg2%zeta_as_inv_22*log_match**2+seg2%zeta_as_inv_21*log_match&
                  +seg2%zeta_as_inv_20)*aso2pi_start**2+(seg2%zeta_as_inv_33*log_match**3&
                  +seg2%zeta_as_inv_32*log_match**2+seg2%zeta_as_inv_31*log_match&
                  +seg2%zeta_as_inv_30)*aso2pi_start**3)
          CASE(5)
             as_start=as_start*(1d0+seg2%zeta_as_inv_11*log_match*aso2pi_start&
                  +(seg2%zeta_as_inv_22*log_match**2+seg2%zeta_as_inv_21*log_match&
                  +seg2%zeta_as_inv_20)*aso2pi_start**2+(seg2%zeta_as_inv_33*log_match**3&
                  +seg2%zeta_as_inv_32*log_match**2+seg2%zeta_as_inv_31*log_match&
                  +seg2%zeta_as_inv_30)*aso2pi_start**3+(seg2%zeta_as_inv_44*log_match**4&
                  +seg2%zeta_as_inv_43*log_match**3+seg2%zeta_as_inv_42*log_match**2&
                  +seg2%zeta_as_inv_41*log_match+seg2%zeta_as_inv_40)*aso2pi_start**4)
          CASE(1)
             CONTINUE ! do nothing
          CASE DEFAULT
             PRINT *, "ERROR: only 5-loop beta expression (=> 4 loop inverse decoupling relation) is allowed"
             STOP
          END SELECT
       ENDIF
       seg%as_m1(0)=1d0/as_start
       DO i=1, UBOUND(seg%as_m1,dim=1)
          seg%as_m1(i)=as_run_onestep(seg%as_m1(i-1),seg%dt)
       ENDDO
       seg%iflip=-1 ! all i (0-UBOUND(seg%as_m1,dim=1)) is above iflip, upward evolution !
    ENDDO
    ! the below segments then
    DO j=nseg-1, as_box%n_low, -1
       seg => as_box%seg(j)
       ! nf -> nf-1 (decoupling relation)
       ! the as^(nf) at threshold muMatch_mQuark*m^(nf) 
       as_start = 1d0/as_box%seg(j+1)%as_m1(0)
       aso2pi_start = as_start/TWOPI
       IF(as_decoupling_on.AND.as_box%fix_nf.EQ.no_fix_nf)THEN
          SELECT CASE(as_box%n_loop)
          CASE(2)
             as_start=as_start*(1d0+seg%zeta_as_11*log_match*aso2pi_start)
          CASE(3)
             as_start=as_start*(1d0+seg%zeta_as_11*log_match*aso2pi_start&
                  +(seg%zeta_as_22*log_match**2+seg%zeta_as_21*log_match&
                  +seg%zeta_as_20)*aso2pi_start**2)
          CASE(4)
             as_start=as_start*(1d0+seg%zeta_as_11*log_match*aso2pi_start&
                  +(seg%zeta_as_22*log_match**2+seg%zeta_as_21*log_match&
                  +seg%zeta_as_20)*aso2pi_start**2+(seg%zeta_as_33*log_match**3&
                  +seg%zeta_as_32*log_match**2+seg%zeta_as_31*log_match&
                  +seg%zeta_as_30)*aso2pi_start**3)
          CASE(5)
             as_start=as_start*(1d0+seg%zeta_as_11*log_match*aso2pi_start&
                  +(seg%zeta_as_22*log_match**2+seg%zeta_as_21*log_match&
                  +seg%zeta_as_20)*aso2pi_start**2+(seg%zeta_as_33*log_match**3&
                  +seg%zeta_as_32*log_match**2+seg%zeta_as_31*log_match&
                  +seg%zeta_as_30)*aso2pi_start**3+(seg%zeta_as_44*log_match**4&
                  +seg%zeta_as_43*log_match**3+seg%zeta_as_42*log_match**2&
                  +seg%zeta_as_41*log_match+seg%zeta_as_40)*aso2pi_start**4)
          CASE(1)
             CONTINUE ! do nothing
          CASE DEFAULT
             PRINT *, "ERROR: only 5-loop beta expression (=> 4 loop decoupling relation) is allowed"
             STOP
          END SELECT
       ENDIF
       seg%as_m1(UBOUND(seg%as_m1,dim=1))=1d0/as_start
       DO i=UBOUND(seg%as_m1,dim=1)-1, 0, -1
          seg%as_m1(i)=as_run_onestep(seg%as_m1(i+1), -seg%dt)
       ENDDO
       seg%iflip=UBOUND(seg%as_m1,dim=1)+1 ! all i (0-UBOUND(seg%as_m1,dim=1)) is below iflip, downward evolution ! 
    ENDDO

    ! restore the original nf setup
    CALL SET_QCD_NF(nf_orig)

    RETURN
  END SUBROUTINE as_run_init

  ! get the as value at Q with a given as_box generated with as_run_init
  FUNCTION as_value(as_box, Q, fix_nf) RESULT(res)
    IMPLICIT NONE
    TYPE(as_container),INTENT(IN), TARGET::as_box
    REAL(KIND(1d0)),INTENT(IN)::Q
    INTEGER,INTENT(IN),OPTIONAL::fix_nf
    REAL(KIND(1d0))::res
    REAL(KIND(1d0))::t
    INTEGER::nseg,i,n
    REAL(KIND(1d0))::t_rescaled,delta_t,as_m1
    INTEGER::n_warn=0
    SAVE n_warn
    INTEGER,PARAMETER::max_warn=2
    
    t=tOfQ(Q)
    IF(PRESENT(fix_nf).AND.as_box%fix_nf.EQ.no_fix_nf)THEN
       IF(fix_nf.LT.as_box%n_low.OR.fix_nf.GT.as_box%n_high)THEN
          PRINT *, "ERROR: the fix_nf requested is outside the range of as_box"
          STOP
       ENDIF
       IF(t.LT.t_low.OR.t.GT.t_high)THEN
          PRINT *, "ERROR: the Q value is outside the range of as_box"
          STOP
       ENDIF
       nseg=fix_nf
    ELSE
       IF(PRESENT(fix_nf).AND.as_box%fix_nf.NE.no_fix_nf)THEN
          ! as_box is generated with a fixed flavour of nf
          IF(fix_nf.NE.as_box%fix_nf)THEN
             PRINT *, "ERROR: the fix_nf requested does not match with that in as_box"
             STOP
          ENDIF
       ENDIF
       DO nseg=as_box%n_low,as_box%n_high
          IF(t.LE.as_box%seg(nseg)%t_high.AND.t.GE.as_box%seg(nseg)%t_low)EXIT
       ENDDO
       IF(nseg.GT.as_box%n_high)THEN
          PRINT *, "ERROR: A given Q is not supported by as_box"
          STOP
       ENDIF
    ENDIF

    seg => as_box%seg(nseg)

    t_rescaled=(t-seg%t_low)/seg%dt
    IF(t_rescaled.GT.seg%iflip)THEN
       i=FLOOR(t_rescaled)
    ELSE
       i=CEILING(t_rescaled)
    ENDIF
    i=MAX(0,MIN(UBOUND(seg%as_m1,dim=1),i))

    delta_t=t-(seg%t_low+i*seg%dt)
    IF(DABS(delta_t).LE.1.3d0*seg%dt)THEN
       res=1d0/as_run_onestep(seg%as_m1(i), delta_t)
    ELSE
       IF(n_warn.LT.max_warn)THEN
          PRINT *, "WARNING: it will evolve fixed-nf alpha_s beyond calculated range. To update may take a while"
          n_warn=n_warn+1
       ENDIF
       n=CEILING(DABS(delta_t)/seg%dt)
       delta_t=delta_t/n
       as_m1=seg%as_m1(i)
       DO i=1,n
          as_m1=as_run_onestep(as_m1,delta_t)
       ENDDO
       res=1d0/as_m1
    ENDIF
    RETURN
  END FUNCTION as_value

  FUNCTION as_run_onestep(as_m1, dt) result(res)
    USE runge_kutta
    IMPLICIT NONE
    ! as_m1=as^(-1)
    REAL(KIND(1d0)),INTENT(IN):: as_m1,dt
    REAL(KIND(1d0))::res,t
    t=0d0
    res=as_m1
    CALL rkstp(dt,t,res,dasm1_dt)
    RETURN
  END FUNCTION as_run_onestep

  SUBROUTINE dasm1_dt(t, as_m1, dasm1)
    ! the right-handed side of das^(-1)/dt
    USE global_constants
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::t,as_m1
    REAL(KIND(1d0)),INTENT(OUT)::dasm1
    dasm1 = seg%beta(0)/TWOPI+seg%beta(1)/as_m1/TWOPI2+seg%beta(2)/as_m1**2/TWOPI3&
         +seg%beta(3)/as_m1**3/TWOPI4+seg%beta(4)/as_m1**4/TWOPI5
    RETURN
  END SUBROUTINE dasm1_dt

  FUNCTION tOfQ(Q)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::Q
    REAL(KIND(1d0))::tofQ
    tOfQ=2d0*DLOG(Q)
    RETURN
  END FUNCTION tOfQ

  FUNCTION QOft(t)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::t
    REAL(KIND(1d0))::QOft
    QOft=DEXP(0.5d0*t)
    RETURN
  END FUNCTION QOft

  ! HELPER ROUTINES AND FUNCTIONS TO GET AND SET CONDITIONS

  ! set the value of dt to decide the spacing of the alpha_s evolution
  SUBROUTINE SET_AS_DT(new_dt)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::new_dt
    dt_used = new_dt
    RETURN
  END SUBROUTINE SET_AS_DT

  ! get the value of dt to decide the spacing of the alpha_s evolution
  FUNCTION GET_AS_DT()
    IMPLICIT NONE
    REAL(KIND(1d0))::GET_AS_DT
    GET_AS_DT=dt_used
    RETURN
  END FUNCTION GET_AS_DT

  ! get the number of loops in running of this coupling
  FUNCTION GET_as_NumberOfLoop(as_box)
    IMPLICIT NONE
    INTEGER::GET_as_NumberOfLoop
    TYPE(as_container),INTENT(IN)::as_box
    GET_as_NumberOfLoop=as_box%n_loop
    RETURN
  END FUNCTION GET_as_NumberOfLoop

  ! Cleaning up
  SUBROUTINE as_clean(as_box)
    IMPLICIT NONE
    TYPE(as_container),INTENT(INOUT)::as_box
    INTEGER::i, status
    DO i=as_box%n_low, as_box%n_high
       DEALLOCATE(as_box%seg(i)%as_m1,stat=status)
    END DO
    DEALLOCATE(as_box%seg, stat=status)
    RETURN
  END SUBROUTINE as_clean

  ! get the quark mass of iflv-th flavour
  FUNCTION Get_as_QuarkMass(as_box, iflv) RESULT(mass)
    IMPLICIT NONE
    REAL(KIND(1d0))::mass
    TYPE(as_container),INTENT(IN)::as_box
    INTEGER,INTENT(IN)::iflv
    IF(iflv.GT.6.OR.iflv.LT.1)THEN
       PRINT *, "ERROR: iflv should be within 1-6"
       STOP
    ENDIF
    mass=as_box%quark_masses(iflv)
    RETURN
  END FUNCTION Get_as_QuarkMass

  ! get the type of quark mass
  FUNCTION Get_as_QuarkMassType(as_box) RESULT(type)
    IMPLICIT NONE
    INTEGER::type
    TYPE(as_container),INTENT(IN)::as_box
    type=as_box%masses_type
    RETURN
  END FUNCTION Get_as_QuarkMassType
  
  ! get the range of nf
  SUBROUTINE Get_as_nf_range(as_box,nf_min,nf_max)
    IMPLICIT NONE
    TYPE(as_container),INTENT(IN)::as_box
    INTEGER,INTENT(OUT)::nf_min,nf_max
    IF(as_box%fix_nf.EQ.no_fix_nf)THEN
       nf_min=as_box%n_low
       nf_max=as_box%n_high
    ELSE
       nf_min=as_box%fix_nf
       nf_max=nf_min
    ENDIF
    RETURN
  END SUBROUTINE Get_as_nf_range

  ! get nf at a given scale Q
  ! it will also return the range of Q_min, Q_max in that segment
  SUBROUTINE Get_as_nfAtQ(as_box, Q, nfAtQ, Q_min, Q_max, muM_mQ)
    IMPLICIT NONE
    TYPE(as_container),INTENT(IN)::as_box
    REAL(KIND(1d0)),INTENT(IN)::Q
    INTEGER,INTENT(OUT)::nfAtQ
    REAL(KIND(1d0)),INTENT(OUT),OPTIONAL::Q_min,Q_max
    REAL(KIND(1d0)),INTENT(IN),OPTIONAL::muM_mQ
    REAL(KIND(1d0))::deltat_match,t
    INTEGER::nseg
    REAL(KIND(1d0)),DIMENSION(as_box%n_low:as_box%n_high)::mt_low,mt_high
    
    deltat_match=tofQ(default_or_opt(as_box%muMatch_mQuark,muM_mQ)/as_box%muMatch_mQuark)
    
    ! redefine limits for the new matching threshold muM_mQ
    mt_low(as_box%n_low+1:as_box%n_high)=as_box%seg(as_box%n_low+1:as_box%n_high)%t_low+deltat_match
    mt_high(as_box%n_low:as_box%n_high-1)=as_box%seg(as_box%n_low:as_box%n_high-1)%t_high+deltat_match
    mt_low(as_box%n_low) = as_box%seg(as_box%n_low)%t_low ! do not change the outter limit
    mt_high(as_box%n_high) = as_box%seg(as_box%n_high)%t_high ! do not change the outter limit
    
    t=tOfQ(Q)
    IF(as_box%fix_nf.EQ.no_fix_nf)THEN
       ! in a varied flavour number scheme
       ! get the segments of the t first with the new matching threshold
       DO nseg=as_box%n_low, as_box%n_high
          IF(t.LE.mt_high(nseg).AND.t.GE.mt_low(nseg))EXIT
       ENDDO
       IF(nseg.GT.as_box%n_high)THEN
          PRINT *, "ERROR: the Q is outside the supported range of as_box #1"
          STOP
       ENDIF
       nfAtQ=nseg
       IF(PRESENT(Q_min).AND.PRESENT(Q_max))THEN
          Q_min=QOft(mt_low(nseg))
          Q_max=QOft(mt_high(nseg))
       ENDIF
    ELSE
       IF(t.GT.t_high.OR.t.LT.t_low)THEN
          PRINT *, "ERROR: the Q is outside the supported range of as_box #2"
          STOP
       ENDIF
       nfAtQ=as_box%fix_nf
       IF(PRESENT(Q_min).AND.PRESENT(Q_max))THEN
          Q_min=QOft(t_low)
          Q_max=QOft(t_high)
       ENDIF
    ENDIF
    RETURN
  END SUBROUTINE Get_as_nfAtQ

  ! get the Q range for a given value of nf=iflv
  SUBROUTINE Get_as_QRangeAtnf(as_box, iflv, Q_min, Q_max, muM_mQ)
    IMPLICIT NONE
    TYPE(as_container), INTENT(IN)::as_box
    INTEGER, INTENT(IN)::iflv
    REAL(KIND(1d0)),INTENT(OUT)::Q_min, Q_max
    REAL(KIND(1d0)),INTENT(IN),OPTIONAL::muM_mQ
    REAL(KIND(1d0))::deltat_match
    
    deltat_match = tofQ(default_or_opt(as_box%muMatch_mQuark, muM_mQ)/as_box%muMatch_mQuark)
    
    IF(as_box%fix_nf.EQ.no_fix_nf)THEN
       ! varied flavour number scheme
       IF(iflv.LT.as_box%n_low.OR.iflv.GT.as_box%n_high)THEN
          PRINT *, "ERROR: nf value iflv is not in range of as_box"
          STOP
       ENDIF
       ! special treatment of the outter limits
       IF(iflv.EQ.as_box%n_low)THEN
          Q_min=QOft(as_box%seg(iflv)%t_low)
       ELSE
          Q_min=QOft(as_box%seg(iflv)%t_low+deltat_match)
       ENDIF
       IF(iflv.EQ.as_box%n_high)THEN
          Q_max=QOft(as_box%seg(iflv)%t_high)
       ELSE
          Q_max=QOft(as_box%seg(iflv)%t_high+deltat_match)
       ENDIF
    ELSE
       IF(iflv.NE.as_box%fix_nf)THEN
          PRINT *, "ERROR: nf vlaue iflv does not match to as_box"
          STOP
       ENDIF
       Q_min=QOft(t_low)
       Q_max=QOft(t_high)
    ENDIF
    RETURN
  END SUBROUTINE Get_as_QRangeAtnf

!=============================================================================!
!=                                                                           =!
!=  ENDING OF GENERAL alpha_s running with runge kutta method                =!
!=                                                                           =!
!=============================================================================!

!=============================================================================!
!=                                                                           =!
!=  BEGINING OF special or approximated alpha_s running routines             =!
!=                                                                           =!
!=============================================================================!

  FUNCTION ALPHAS_FROM_LAMBDA(Q,n_loop,Q_THR)
    USE global_constants
    USE qcd_constants
    USE qcd_setup
    USE newton_method
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::Q
    INTEGER,INTENT(IN)::n_loop
    REAL(KIND(1d0)),DIMENSION(6),OPTIONAL::Q_THR
    REAL(KIND(1d0))::ALPHAS_FROM_LAMBDA
    INTEGER::nq, k,i
    REAL(KIND(1d0))::Q2
    REAL(KIND(1d0)),DIMENSION(6)::QTHR
    ! Heavy Quark Lambda Values
    !REAL(KIND(1d0)),DIMENSION(3:6)::LAMBDAL=(/0.2470d0,0.2150d0,0.1650d0,0.0847d0 /),&
    !     LAMBDAN=(/ 0.3500d0, 0.3260d0, 0.2262d0, 0.1510d0 /)
    ! Heavy Quark Mass Thresholds mc=1.4,mb=4.5,mt=180.0
    !REAL(KIND(1d0)),DIMENSION(3)::Q2THR=(/1.690d0,20.25d0,32400d0 /)
    REAL(KIND(1d0)),DIMENSION(6),PARAMETER::Q_THR_default=(/1d-10,1d-10,0.15d0,&
         1.4d0,4.5d0,180d0 /)
    REAL(KIND(1d0))::ALP,ALP0,LAM2,XLP,XL,XLM,Y,Y1
    REAL(KIND(1d0))::B1,B2,B3,B4,beta0LQ2
    REAL(KIND(1d0))::LQ2,LLQ2
    INTEGER::NLOOP
    COMMON/ALPHAS_LAMBDA/LQ2,NLOOP
    INTEGER::nf_orig
    REAL(KIND(1d0))::FPALP
    REAL(KIND(1d0)),DIMENSION(2)::FPALP2

    nf_orig=nf
    
    Q2=Q**2
    NLOOP=n_loop
    QTHR(1:6)=default_or_opt(Q_THR_default(1:6),Q_THR)
    !...DETERMINATION OF THE APPROPRIATE NUMBER OF FLAVOURS :
    nq=0
    DO K=1,6
       IF(Q.GT.QTHR(K))THEN
          nq=nq+1
       ELSE
          EXIT
       ENDIF
    ENDDO

    CALL SET_QCD_NF(nq)

    SELECT CASE(n_loop)
    CASE(1)
       LAM2=LAMBDA_1L(nf)*LAMBDA_1L(nf)
    CASE(2)
       LAM2=LAMBDA_2L(nf)*LAMBDA_2L(nf)
    CASE(3)
       LAM2=LAMBDA_3L(nf)*LAMBDA_3L(nf)
    CASE(4)
       LAM2=LAMBDA_4L(nf)*LAMBDA_4L(nf)
    CASE(5)
       LAM2=LAMBDA_5L(nf)*LAMBDA_5L(nf)
    CASE DEFAULT
       PRINT *, "ERROR: only <= 5-loop as running supported by ALPHAS_FROM_LAMBDA #1"
       STOP
    END SELECT
    
    ! check Eq.(5) in hep-ph/0004189
    ! ...LO ALPHA_S
    IF(n_loop.EQ.1)THEN
       LQ2=DLOG(Q2/LAM2)
       ! it is alpha_s/2/pi
       ALP=1.d0/(beta0*LQ2)
    ! ...BEYOND LO ALPHA_S (up to N^4LO=5-loop)
    ELSEIF(n_loop.GE.2)THEN
       B1=beta1/beta0
       B2=beta2/beta0
       B3=beta3/beta0
       B4=beta4/beta0
       LQ2=DLOG(Q2/LAM2)
       beta0LQ2=beta0*LQ2
       LLQ2=DLOG(LQ2)
       ! Eq(5) of hep-ph/0004189
       ! note that the following approximated is obtained via expansion in 1/LQ2
       SELECT CASE(n_loop)
       CASE(2)
          ALP=1D0/(beta0LQ2)*(1D0-B1/beta0LQ2*LLQ2)
       CASE(3)
          ALP=1D0/(beta0LQ2)*(1D0-B1/beta0LQ2*LLQ2&
               +1D0/beta0LQ2**2*(B1**2*(LLQ2**2-LLQ2-1D0)+B2))
       CASE(4)
          ! the exact expression does not matter as it will 
          ALP=1D0/(beta0LQ2)*(1D0-B1/beta0LQ2*LLQ2&
               +1D0/beta0LQ2**2*(B1**2*(LLQ2**2-LLQ2-1D0)+B2)&
               +1D0/beta0LQ2**3*(B1**3*(-LLQ2**3+2.5D0*LLQ2**2+2D0*LLQ2-0.5D0)&
               -3D0*B1*B2*LLQ2+0.5D0*B3))
       CASE(5)
          ! the exact expression does not matter as it will
          ALP=1D0/(beta0LQ2)*(1D0-B1/beta0LQ2*LLQ2&
               +1D0/beta0LQ2**2*(B1**2*(LLQ2**2-LLQ2-1D0)+B2)&
               +1D0/beta0LQ2**3*(B1**3*(-LLQ2**3+2.5D0*LLQ2**2+2D0*LLQ2-0.5D0)&
               -3D0*B1*B2*LLQ2+0.5D0*B3)&
               +1D0/beta0LQ2**4*(B1**4*(LLQ2**4-13D0/3D0*LLQ2**3-1.5D0*LLQ2**2+4D0*LLQ2+7D0/6D0)&
               +3D0*B1**2*B2*(LLQ2-1D0)*(2D0*LLQ2+1D0)&
               -B1*B3*(2D0*LLQ2+1D0/6D0)+5D0/3D0*B2**2+1D0/3D0*B4))
       CASE DEFAULT
          PRINT *, "ERROR: only <= 5-loop as running supported by ALPHAS_FROM_LAMBDA #2"
          STOP
       END SELECT
       ALP0=ALP
       ! ... EXACT NLO VALUE, FOUND VIA NEWTON PROCEDURE (TO RECOVER MORE LQ2 TERMS)
       ! f(x)=LHS(Eq.(4))-RHS(Eq.(4)), where Eq.(4) is the equation in hep-ph/0004189
       ! Find root of f(x)=0, use iteration of x_n=x_(n-1)-f(x_(n-1))/f'(x_(n-1))
       ! Here, f(ALP)=Y. Note, ALP=alpha_s/2/pi
       !IF(n_loop.EQ.2)THEN
       !   ! Note: it is necessary to use quadratic Netown's method than the linear Newton's method
       !   ! Especially for Q close to 0.5 GeV at 2-loop level calculation. Otherwise, the linear method
       !   ! has too bad convergence !!!
       !   CALL newtonsolver(ALP,FPALP2,LAMBDA_IN_ALP_EQ_QUA,1D-14,10)
       !ELSE
       !   CALL newtonsolver(ALP,FPALP, LAMBDA_IN_ALP_EQ, 1D-14,-1)
       !ENDIF
       CALL newtonsolver(ALP,FPALP, LAMBDA_IN_ALP_EQ, 1D-14,-1)
       IF(FPALP.EQ.-1D99.OR.ISNAN(FPALP).OR.1D0/FPALP.EQ.0d0)THEN
         PRINT *, "INFO: STARTING TO USE Quadratic Newton procedure"
          ! Note: it is necessary to use quadratic Netown's method than the linear Newton's method
          ! Especially for Q close to 0.5 GeV beyond 1-loop level calculation. Otherwise, the linear method
          ! has too bad convergence !!!
          ALP=ALP0
          CALL newtonsolver(ALP,FPALP2,LAMBDA_IN_ALP_EQ_QUA,1D-14,30)
       ENDIF
    ELSE
       PRINT *, "ERROR: < 1-loop as run does not support in ALPHAS_FROM_LAMBDA"
       STOP
    ENDIF

    CALL SET_QCD_NF(nf_orig)

    ALPHAS_FROM_LAMBDA=ALP*2D0*PI
    RETURN

  END FUNCTION ALPHAS_FROM_LAMBDA


  SUBROUTINE LAMBDA_IN_ALP_EQ(as,f,fp)
    ! f(x)=LHS(Eq.(4))-RHS(Eq.(4)), where Eq.(4) is the equation in hep-ph/0004189
    USE qcd_constants
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::as ! as is alpha_s/2/pi
    REAL(KIND(1d0)),INTENT(OUT)::f,fp
    REAL(KIND(1d0))::LQ2
    INTEGER::NLOOP
    COMMON/ALPHAS_LAMBDA/LQ2,NLOOP
    REAL(KIND(1d0))::B1,B2,B3,B4,beta0as,BB1,BB2
    REAL(KIND(1d0))::BB3,BB4,BB5,BB6
    IF(as.LT.0d0.OR.ISNAN(as).OR.1D0/as.EQ.0d0)THEN
       PRINT *, "WARNING: THE CONVERGENCE OF LINEAR NEWTWON METHOD IS TOO BAD !"
       PRINT *, "as=",as
       f=0d0
       fp=-1d99
       RETURN
    ENDIF
    SELECT CASE(NLOOP)
    CASE(1)
       f=LQ2-1D0/(beta0*as)
       fp=1D0/(beta0*as**2)
    CASE(2)
       B1=beta1/beta0
       beta0as=beta0*as
       ! it can be analytically solved with betai=0 (i>=2)
       ! We do not use the approximation anymore !!!
       !f=LQ2-1D0/(beta0as)-B1/beta0*DLOG(beta0as)
       !fp=1D0/(beta0*as**2)-B1/(beta0as)
       f=LQ2-1D0/(beta0as)+B1/beta0*DLOG(1D0/(beta0as)+B1/beta0)
       fp=1D0/(beta0*as**2+beta1*as**3)
    CASE(3)
       B1=beta1/beta0
       B2=beta2/beta0
       beta0as=beta0*as
       ! it can be analytically solved with betai=0 (i>=3)
       ! We do not use the approximation anymore !!!
       !BB1=B2-B1**2
       !f=LQ2-1D0/(beta0as)-B1/beta0*DLOG(beta0as)-BB1/beta0*as
       !fp=1D0/(beta0*as**2)-B1/(beta0as)-BB1/beta0
       BB1=4D0*B2-B1**2
       BB2=2D0*B2-B1**2
       IF(BB1.GE.0d0)THEN
          f=LQ2-1D0/(beta0as)+B1/beta0/2D0*DLOG(1D0/(beta0as**2)+B1/beta0/beta0as+B2/beta0**2)&
               -BB2/DSQRT(BB1)/beta0*DATAN(DSQRT(BB1)/(B1+2D0/as))
       ELSE
          f=LQ2-1D0/(beta0as)+B1/beta0/2D0*DLOG(1D0/(beta0as**2)+B1/beta0/beta0as+B2/beta0**2)&
               -BB2/DSQRT(-BB1)/beta0*DATANH(DSQRT(-BB1)/(B1+2D0/as))
       ENDIF
       fp=1D0/(beta0*as**2+beta1*as**3+beta2*as**4)
    CASE(4,5)
       B1=beta1/beta0
       B2=beta2/beta0
       B3=beta3/beta0
       IF(NLOOP.EQ.4)THEN
          B4=0D0
       ELSE
          B4=beta4/beta0
       ENDIF
       beta0as=beta0*as
       BB1=B2-B1**2
       BB2=B3/2D0-B1*B2+B1**3/2D0
       BB3=-1D0/3D0*(-B4+2d0*B1*B3+B2**2-3D0*B1**2*B2+B1**4)
       BB4=-1D0/4D0*(2D0*B1*B4+2D0*B2*B3-3D0*B1**2*B3-3D0*B1*B2**2+4D0*B1**3*B2-B1**5)
       BB5=-1D0/5D0*(2D0*B2*B4-3D0*B1**2*B4&
            +B3**2-6D0*B1*B2*B3+4D0*B1**3*B3-B2**3+6D0*B1**2*B2**2-5D0*B1**4*B2+B1**6)
       BB6=-1D0/6D0*(2D0*B3*B4-6D0*B1*B2*B4+4D0*B1**3*B4&
            -3D0*B1*B3**2-3D0*B2**2*B3+12D0*B1**2*B2*B3-5D0*B1**4*B3+4D0*B1*B2**3&
            -10D0*B1**3*B2**2+6D0*B1**5*B2-B1**7)
       ! it is hard to get the exact solution
       ! Try to keep more terms in order to improve the convergence
       f=LQ2-1D0/(beta0as)-B1/beta0*DLOG(beta0as)-BB1/beta0*as-BB2/beta0*as**2-BB3/beta0*as**3&
            -BB4/beta0*as**4-BB5/beta0*as**5-BB6/beta0*as**6
       fp=1D0/(beta0*as**2)-B1/(beta0as)-BB1/beta0-BB2/beta0*2D0*as-BB3/beta0*3D0*as**2&
            -BB4/beta0*4D0*as**3-BB5/beta0*5D0*as**4-BB6/beta0*6D0*as**5
    CASE DEFAULT
       PRINT *, "ERROR: only maximal 5-loop LAMBDA_IN_ALP_EQ relation supports"
       STOP
    END SELECT
    RETURN
  END SUBROUTINE LAMBDA_IN_ALP_EQ

  SUBROUTINE LAMBDA_IN_ALP_EQ_QUA(as,f,fp,fpp)
    ! f(x)=LHS(Eq.(4))-RHS(Eq.(4)), where Eq.(4) is the equation in hep-ph/0004189
    USE qcd_constants
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::as ! as is alpha_s/2/pi
    REAL(KIND(1d0)),INTENT(OUT)::f,fp,fpp
    REAL(KIND(1d0))::LQ2
    INTEGER::NLOOP
    COMMON/ALPHAS_LAMBDA/LQ2,NLOOP
    REAL(KIND(1d0))::B1,B2,B3,B4,beta0as,BB1,BB2
    REAL(KIND(1d0))::BB3,BB4,BB5,BB6
    IF(as.LT.0d0.OR.ISNAN(as).OR.1D0/as.EQ.0d0)THEN
       PRINT *, "ERROR: Quadratic Newton procedure cannot be applied here! Check it first !"
       PRINT *, "NLOOP=",NLOOP
       PRINT *, "LQ2=",LQ2
       PRINT *, "as=",as
       STOP
    ENDIF
    SELECT CASE(NLOOP)
    CASE(1)
       f=LQ2-1D0/(beta0*as)
       fp=1D0/(beta0*as**2)
       fpp=-fp**2*2D0*beta0*as
    CASE(2)
       B1=beta1/beta0
       beta0as=beta0*as
       ! it can be analytically solved with betai=0 (i>=2)
       ! We do not use the approximation anymore !!!
       !f=LQ2-1D0/(beta0as)-B1/beta0*DLOG(beta0as)
       !fp=1D0/(beta0*as**2)-B1/(beta0as)
       f=LQ2-1D0/(beta0as)+B1/beta0*DLOG(1D0/(beta0as)+B1/beta0)
       fp=1D0/(beta0*as**2+beta1*as**3)
       fpp=-fp**2*(2D0*beta0*as+3D0*beta1*as**2)
    CASE(3)
       B1=beta1/beta0
       B2=beta2/beta0
       beta0as=beta0*as
       ! it can be analytically solved with betai=0 (i>=3)
       ! We do not use the approximation anymore !!!
       !BB1=B2-B1**2
       !f=LQ2-1D0/(beta0as)-B1/beta0*DLOG(beta0as)-BB1/beta0*as
       !fp=1D0/(beta0*as**2)-B1/(beta0as)-BB1/beta0
       BB1=4D0*B2-B1**2
       BB2=2D0*B2-B1**2
       IF(BB1.GE.0d0)THEN
          f=LQ2-1D0/(beta0as)+B1/beta0/2D0*DLOG(1D0/(beta0as**2)+B1/beta0/beta0as+B2/beta0**2)&
               -BB2/DSQRT(BB1)/beta0*DATAN(DSQRT(BB1)/(B1+2D0/as))
       ELSE
          f=LQ2-1D0/(beta0as)+B1/beta0/2D0*DLOG(1D0/(beta0as**2)+B1/beta0/beta0as+B2/beta0**2)&
               -BB2/DSQRT(-BB1)/beta0*DATANH(DSQRT(-BB1)/(B1+2D0/as))
       ENDIF
       fp=1D0/(beta0*as**2+beta1*as**3+beta2*as**4)
       fpp=-fp**2*(2D0*beta0*as+3D0*beta1*as**2+4D0*beta2*as**3)
    CASE(4,5)
       B1=beta1/beta0
       B2=beta2/beta0
       B3=beta3/beta0
       IF(NLOOP.EQ.4)THEN
          B4=0D0
       ELSE
          B4=beta4/beta0
       ENDIF
       beta0as=beta0*as
       BB1=B2-B1**2
       BB2=B3/2D0-B1*B2+B1**3/2D0
       BB3=-1D0/3D0*(-B4+2d0*B1*B3+B2**2-3D0*B1**2*B2+B1**4)
       BB4=-1D0/4D0*(2D0*B1*B4+2D0*B2*B3-3D0*B1**2*B3-3D0*B1*B2**2+4D0*B1**3*B2-B1**5)
       BB5=-1D0/5D0*(2D0*B2*B4-3D0*B1**2*B4&
            +B3**2-6D0*B1*B2*B3+4D0*B1**3*B3-B2**3+6D0*B1**2*B2**2-5D0*B1**4*B2+B1**6)
       BB6=-1D0/6D0*(2D0*B3*B4-6D0*B1*B2*B4+4D0*B1**3*B4&
            -3D0*B1*B3**2-3D0*B2**2*B3+12D0*B1**2*B2*B3-5D0*B1**4*B3+4D0*B1*B2**3&
            -10D0*B1**3*B2**2+6D0*B1**5*B2-B1**7)
       ! it is hard to get the exact solution
       ! Try to keep more terms in order to improve the convergence
       f=LQ2-1D0/(beta0as)-B1/beta0*DLOG(beta0as)-BB1/beta0*as-BB2/beta0*as**2-BB3/beta0*as**3&
            -BB4/beta0*as**4-BB5/beta0*as**5-BB6/beta0*as**6
       fp=1D0/(beta0*as**2)-B1/(beta0as)-BB1/beta0-BB2/beta0*2D0*as-BB3/beta0*3D0*as**2&
            -BB4/beta0*4D0*as**3-BB5/beta0*5D0*as**4-BB6/beta0*6D0*as**5
       fpp=-2D0/(beta0*as**3)+B1/(beta0*as**2)-BB2/beta0*2D0-BB3/beta0*6D0*as&
            -BB4/beta0*12D0*as**2-BB5/beta0*20D0*as**3-BB6/beta0*30D0*as**4
    CASE DEFAULT
       PRINT *, "ERROR: only maximal 5-loop LAMBDA_IN_ALP_EQ_QUA relation supports"
       STOP
    END SELECT
    RETURN
  END SUBROUTINE LAMBDA_IN_ALP_EQ_QUA

  FUNCTION LAMBDA_FROM_ALPHAS(alpha_s, Q, n_loop, nf_in)
    ! Get Lambda up to 4-loop via Eq.(4) in hep-ph/0004189
    USE global_constants
    USE qcd_setup
    USE qcd_constants
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::alpha_s,Q
    INTEGER,INTENT(IN)::n_loop
    INTEGER,INTENT(IN),OPTIONAL::nf_in
    REAL(KIND(1d0))::LAMBDA_FROM_ALPHAS
    INTEGER::nf_orig,nf_used
    REAL(KIND(1d0))::b1, b2, b3, b4,as, RHS, BB1
    INTEGER::init=0
    SAVE init
    nf_orig=nf
    nf_used=default_or_opt(nf,nf_in)
    IF(nf_used.NE.nf_orig.OR.init.EQ.0)CALL SET_QCD_NF(nf_used)
    as=alpha_s/2D0/PI
    SELECT CASE(n_loop)
    CASE(1)
       b1=0d0
       b2=0d0
       b3=0d0
       b4=0d0
       ! EXACT INTEGRATION
       RHS=1d0/beta0*(1d0/as)+b1/beta0*DLOG(beta0)
    CASE(2)
       b1=beta1/beta0
       b2=0d0
       b3=0d0
       b4=0d0
       ! EXACT INTEGRATION
       RHS=1d0/beta0*(1d0/as-b1*DLOG(1D0/(beta0*as)+b1/beta0))
    CASE(3)
       b1=beta1/beta0
       b2=beta2/beta0
       b3=0d0
       b4=0d0
       BB1=4D0*b2-b1**2
       ! EXACT INTEGRATION
       IF(BB1.GE.0d0)THEN
          RHS=1d0/beta0*(1d0/as-0.5D0*b1*DLOG(1D0/beta0**2/as**2+b1/beta0**2/as+b2/beta0**2)&
               +(2D0*b2-b1**2)/DSQRT(BB1)*DATAN(DSQRT(BB1)/(b1+2D0/as)))
       ELSE
          RHS=1d0/beta0*(1d0/as-0.5D0*b1*DLOG(1D0/beta0**2/as**2+b1/beta0**2/as+b2/beta0**2)&
               +(2D0*b2-b1**2)/DSQRT(-BB1)*DATANH(DSQRT(-BB1)/(b1+2D0/as)))
       ENDIF
    CASE(4)
       b1=beta1/beta0
       b2=beta2/beta0
       b3=beta3/beta0
       b4=0d0
       ! NOTE THE FOLLOWING IS AN APPROXIMATION
       RHS=1d0/beta0*(1d0/as+b1*DLOG(as)+(b2-b1**2)*as+(b3/2d0-b1*b2+b1**3/2d0)*as**2)+b1/beta0*DLOG(beta0)
    CASE(5)
       b1=beta1/beta0
       b2=beta2/beta0
       b3=beta3/beta0
       b4=beta4/beta0
       ! NOTE THE FOLLOWING IS AN APPROXIMATION
       RHS=1d0/beta0*(1d0/as+b1*DLOG(as)+(b2-b1**2)*as+(b3/2d0-b1*b2+b1**3/2d0)*as**2&
            +1D0/3D0*(b4-2D0*b1*b3-b2**2+3D0*b1**2*b2-b1**4)*as**3)+b1/beta0*DLOG(beta0)
    CASE DEFAULT
       PRINT *, "ERROR: do not support > 5-loop beta function"
       STOP
    END SELECT
    LAMBDA_FROM_ALPHAS=Q*DEXP(-RHS/2d0)
    IF(nf_used.NE.nf_orig.OR.init.EQ.0)CALL SET_QCD_NF(nf_orig)
    init=1
    RETURN
  END FUNCTION LAMBDA_FROM_ALPHAS

  SUBROUTINE ALPHAS_LAMBDA_INIT(as_frame,alpha_s, Q, quark_masses, masses_type, muMatch_mQuark)
    USE qcd_constants
    USE global_constants
    USE qcd_setup
    IMPLICIT NONE
    TYPE(as_framework),INTENT(OUT),TARGET::as_frame
    REAL(KIND(1d0)),INTENT(IN),OPTIONAL::alpha_s,Q
    REAL(KIND(1d0)),DIMENSION(4:6),INTENT(IN),OPTIONAL::quark_masses
    INTEGER,INTENT(IN),OPTIONAL::masses_type
    REAL(KIND(1d0)),INTENT(IN),OPTIONAL::muMatch_mQuark
    ! local variables
    INTEGER::nf_orig,i,j,nseg,n_loop
    INTEGER,PARAMETER::n_loop_max=5 ! maximal is 5 loop
    TYPE(as_nf_interval),POINTER::segi,segi2
    REAL(KIND(1d0))::log_match
    REAL(KIND(1d0)),DIMENSION(6),PARAMETER::Q_THR_default=(/1d-10,1d-10,0.15d0,&
         1.4d0,4.5d0,180d0 /)
    REAL(KIND(1d0)),DIMENSION(6)::Q_THR
    REAL(KIND(1d0))::as_start,aso2pi_start
    
    nf_orig=nf

    as_frame%alpha_s = default_or_opt(0.118d0, alpha_s) ! initial value of alpha_s
    as_frame%Q = default_or_opt(91.2d0, Q) ! initial value of scale
    as_frame%muMatch_mQuark = default_or_opt(ONE, muMatch_mQuark) ! the threshold is muMatch_mQuark*quark_mass
    as_frame%masses_type = default_or_opt(2, masses_type) ! 1: MSbar mass; 2: OS; 3: SI
    as_frame%quark_masses(1:3) = quark_masses_default(1:3) ! quark masses
    as_frame%quark_masses(4:6) = default_or_opt(quark_masses_default(4:6),quark_masses)
    as_frame%Q_THR(1:6)=Q_THR_default(1:6)


    DO i=LBOUND(as_frame%quark_masses, dim=1), UBOUND(as_frame%quark_masses, dim=1)
       IF(t_low.GT.tOfQ(as_frame%muMatch_mQuark*as_frame%quark_masses(i))) as_frame%n_low = i
       IF(t_high.GT.tOfQ(as_frame%muMatch_mQuark*as_frame%quark_masses(i))) as_frame%n_high = i
    ENDDO

    ALLOCATE(as_frame%seg(as_frame%n_low:as_frame%n_high))
    DO i=as_frame%n_low,as_frame%n_high
       as_frame%Q_THR(i)=as_frame%muMatch_mQuark*as_frame%quark_masses(i)
       segi => as_frame%seg(i)
       segi%Q_low=MAX(QOft(t_low), as_frame%Q_THR(i))
       IF(i.LT.as_frame%n_high)THEN
          segi%Q_high=as_frame%muMatch_mQuark*as_frame%quark_masses(i+1)
       ELSE
          segi%Q_high=QOft(t_high)
       ENDIF

       CALL SET_QCD_NF(i)

       ! the decoupling relation when changing nf
       SELECT CASE(as_frame%masses_type)
       CASE(1)
          ! quark masses are MSbar masses
          segi%zeta_as_11=zeta_alphas_MSbar_11
          segi%zeta_as_22=zeta_alphas_MSbar_22
          segi%zeta_as_21=zeta_alphas_MSbar_21
          segi%zeta_as_20=zeta_alphas_MSbar_20
          segi%zeta_as_33=zeta_alphas_MSbar_33
          segi%zeta_as_32=zeta_alphas_MSbar_32
          segi%zeta_as_31=zeta_alphas_MSbar_31
          segi%zeta_as_30=zeta_alphas_MSbar_30
          segi%zeta_as_44=zeta_alphas_MSbar_44
          segi%zeta_as_43=zeta_alphas_MSbar_43
          segi%zeta_as_42=zeta_alphas_MSbar_42
          segi%zeta_as_41=zeta_alphas_MSbar_41
          segi%zeta_as_40=zeta_alphas_MSbar_40

          segi%zeta_as_inv_11=zeta_alphas_inv_MSbar_11
          segi%zeta_as_inv_22=zeta_alphas_inv_MSbar_22
          segi%zeta_as_inv_21=zeta_alphas_inv_MSbar_21
          segi%zeta_as_inv_20=zeta_alphas_inv_MSbar_20
          segi%zeta_as_inv_33=zeta_alphas_inv_MSbar_33
          segi%zeta_as_inv_32=zeta_alphas_inv_MSbar_32
          segi%zeta_as_inv_31=zeta_alphas_inv_MSbar_31
          segi%zeta_as_inv_30=zeta_alphas_inv_MSbar_30
          segi%zeta_as_inv_44=zeta_alphas_inv_MSbar_44
          segi%zeta_as_inv_43=zeta_alphas_inv_MSbar_43
          segi%zeta_as_inv_42=zeta_alphas_inv_MSbar_42
          segi%zeta_as_inv_41=zeta_alphas_inv_MSbar_41
          segi%zeta_as_inv_40=zeta_alphas_inv_MSbar_40
       CASE(2)
          ! quark masses are OS masses
          segi%zeta_as_11=zeta_alphas_OS_11
          segi%zeta_as_22=zeta_alphas_OS_22
          segi%zeta_as_21=zeta_alphas_OS_21
          segi%zeta_as_20=zeta_alphas_OS_20
          segi%zeta_as_33=zeta_alphas_OS_33
          segi%zeta_as_32=zeta_alphas_OS_32
          segi%zeta_as_31=zeta_alphas_OS_31
          segi%zeta_as_30=zeta_alphas_OS_30
          segi%zeta_as_44=zeta_alphas_OS_44
          segi%zeta_as_43=zeta_alphas_OS_43
          segi%zeta_as_42=zeta_alphas_OS_42
          segi%zeta_as_41=zeta_alphas_OS_41
          segi%zeta_as_40=zeta_alphas_OS_40
       
          segi%zeta_as_inv_11=zeta_alphas_inv_OS_11
          segi%zeta_as_inv_22=zeta_alphas_inv_OS_22
          segi%zeta_as_inv_21=zeta_alphas_inv_OS_21
          segi%zeta_as_inv_20=zeta_alphas_inv_OS_20
          segi%zeta_as_inv_33=zeta_alphas_inv_OS_33
          segi%zeta_as_inv_32=zeta_alphas_inv_OS_32
          segi%zeta_as_inv_31=zeta_alphas_inv_OS_31
          segi%zeta_as_inv_30=zeta_alphas_inv_OS_30
          segi%zeta_as_inv_44=zeta_alphas_inv_OS_44
          segi%zeta_as_inv_43=zeta_alphas_inv_OS_43
          segi%zeta_as_inv_42=zeta_alphas_inv_OS_42
          segi%zeta_as_inv_41=zeta_alphas_inv_OS_41
          segi%zeta_as_inv_40=zeta_alphas_inv_OS_40
       CASE(3)
          ! quark masses are SI masses
          segi%zeta_as_11=zeta_alphas_SI_11
          segi%zeta_as_22=zeta_alphas_SI_22
          segi%zeta_as_21=zeta_alphas_SI_21
          segi%zeta_as_20=zeta_alphas_SI_20
          segi%zeta_as_33=zeta_alphas_SI_33
          segi%zeta_as_32=zeta_alphas_SI_32
          segi%zeta_as_31=zeta_alphas_SI_31
          segi%zeta_as_30=zeta_alphas_SI_30
          segi%zeta_as_44=zeta_alphas_SI_44
          segi%zeta_as_43=zeta_alphas_SI_43
          segi%zeta_as_42=zeta_alphas_SI_42
          segi%zeta_as_41=zeta_alphas_SI_41
          segi%zeta_as_40=zeta_alphas_SI_40

          segi%zeta_as_inv_11=zeta_alphas_inv_SI_11
          segi%zeta_as_inv_22=zeta_alphas_inv_SI_22
          segi%zeta_as_inv_21=zeta_alphas_inv_SI_21
          segi%zeta_as_inv_20=zeta_alphas_inv_SI_20
          segi%zeta_as_inv_33=zeta_alphas_inv_SI_33
          segi%zeta_as_inv_32=zeta_alphas_inv_SI_32
          segi%zeta_as_inv_31=zeta_alphas_inv_SI_31
          segi%zeta_as_inv_30=zeta_alphas_inv_SI_30
          segi%zeta_as_inv_44=zeta_alphas_inv_SI_44
          segi%zeta_as_inv_43=zeta_alphas_inv_SI_43
          segi%zeta_as_inv_42=zeta_alphas_inv_SI_42
          segi%zeta_as_inv_41=zeta_alphas_inv_SI_41
          segi%zeta_as_inv_40=zeta_alphas_inv_SI_40
       CASE DEFAULT
          PRINT *, "ERROR: masses_type can only be MSbar(1)/OS(2)/SI(3) masses in ALPHAS_LAMBDA_INIT"
          STOP
       END SELECT

    ENDDO

    ! determin which segmentation we are in for the initial scale
    DO i=as_frame%n_low,as_frame%n_high
       IF(Q.LE.as_frame%seg(i)%Q_high.AND.Q.GE.as_frame%seg(i)%Q_low)EXIT
    ENDDO
    IF(i.GT.as_frame%n_high)THEN
       PRINT *, "ERROR: the initial scale Q is not in the supported range in ALPHAS_LAMBDA_INIT"
       STOP
    ENDIF
    nseg=i

    Q_THR(1:6)=as_frame%Q_THR(1:6)
    Q_THR(nseg)=Q_THR(nseg)*0.99D0 ! in order to make sure Q_low and Q_high is still with nseg-th flavour
    IF(nseg.LT.6)Q_THR(nseg+1)=Q_THR(nseg+1)*1.01D0
    
    segi => as_frame%seg(nseg)
    ! enters into the different loop orders
    DO n_loop=1,n_loop_max
       SELECT CASE(n_loop)
       CASE(1)
          LAMBDA_1L(nseg)=LAMBDA_FROM_ALPHAS(alpha_s, Q, n_loop, nseg)
       CASE(2)
          LAMBDA_2L(nseg)=LAMBDA_FROM_ALPHAS(alpha_s, Q, n_loop, nseg)
       CASE(3)
          LAMBDA_3L(nseg)=LAMBDA_FROM_ALPHAS(alpha_s, Q, n_loop, nseg)
       CASE(4)
          LAMBDA_4L(nseg)=LAMBDA_FROM_ALPHAS(alpha_s, Q, n_loop, nseg)
       CASE(5)
          LAMBDA_5L(nseg)=LAMBDA_FROM_ALPHAS(alpha_s, Q, n_loop, nseg)
       CASE DEFAULT
          PRINT *, "ERROR: > 5-loop as running is not allowed"
          STOP
       END SELECT
       segi%as_at_Q_low(n_loop)=ALPHAS_FROM_LAMBDA(segi%Q_low,n_loop,Q_THR)
       segi%as_at_Q_high(n_loop)=ALPHAS_FROM_LAMBDA(segi%Q_high,n_loop,Q_THR)
    ENDDO

    ! for other segments
    log_match=2d0*DLOG(as_frame%muMatch_mQuark) ! the ratio of log(mu^2/quark_mass^2) with the matching mu=muMatch_mQuark*quark_mass
    ! the above segments first
    DO j=nseg+1, as_frame%n_high
       Q_THR(1:6)=as_frame%Q_THR(1:6)
       Q_THR(j)=Q_THR(j)*0.99D0 ! in order to make sure Q_low and Q_high is still with nseg-th flavour
       IF(j.LT.6)Q_THR(j+1)=Q_THR(j+1)*1.01D0
       
       segi => as_frame%seg(j)
       segi2 => as_frame%seg(j-1)
       ! nf -> nf+1 (inverse decoupling relation)
       ! the as^(nf) at threshold muMatch_mQuark*m^(nf+1)
       DO n_loop=1,n_loop_max
          as_start=segi2%as_at_Q_high(n_loop)
          aso2pi_start=as_start/TWOPI
          IF(as_decoupling_on)THEN
             SELECT CASE(n_loop)
             CASE(2)
                as_start=as_start*(1d0+segi2%zeta_as_inv_11*log_match*aso2pi_start)
             CASE(3)
                as_start=as_start*(1d0+segi2%zeta_as_inv_11*log_match*aso2pi_start&
                  +(segi2%zeta_as_inv_22*log_match**2+segi2%zeta_as_inv_21*log_match&
                  +segi2%zeta_as_inv_20)*aso2pi_start**2)
             CASE(4)
                as_start=as_start*(1d0+segi2%zeta_as_inv_11*log_match*aso2pi_start&
                  +(segi2%zeta_as_inv_22*log_match**2+segi2%zeta_as_inv_21*log_match&
                  +segi2%zeta_as_inv_20)*aso2pi_start**2+(segi2%zeta_as_inv_33*log_match**3&
                  +segi2%zeta_as_inv_32*log_match**2+segi2%zeta_as_inv_31*log_match&
                  +segi2%zeta_as_inv_30)*aso2pi_start**3)
             CASE(5)
                as_start=as_start*(1d0+segi2%zeta_as_inv_11*log_match*aso2pi_start&
                  +(segi2%zeta_as_inv_22*log_match**2+segi2%zeta_as_inv_21*log_match&
                  +segi2%zeta_as_inv_20)*aso2pi_start**2+(segi2%zeta_as_inv_33*log_match**3&
                  +segi2%zeta_as_inv_32*log_match**2+segi2%zeta_as_inv_31*log_match&
                  +segi2%zeta_as_inv_30)*aso2pi_start**3+(segi2%zeta_as_inv_44*log_match**4&
                  +segi2%zeta_as_inv_43*log_match**3+segi2%zeta_as_inv_42*log_match**2&
                  +segi2%zeta_as_inv_41*log_match+segi2%zeta_as_inv_40)*aso2pi_start**4)
             CASE(1)
                CONTINUE
             CASE DEFAULT
                PRINT *, "ERROR: only 5-loop beta expression (=> 4 loop inverse decoupling relation) is allowed"
                STOP
             END SELECT
          ENDIF
          segi%as_at_Q_low(n_loop)=as_start
          SELECT CASE(n_loop)
          CASE(1)
             LAMBDA_1L(j)=LAMBDA_FROM_ALPHAS(as_start,segi%Q_low,n_loop,j)
          CASE(2)
             LAMBDA_2L(j)=LAMBDA_FROM_ALPHAS(as_start,segi%Q_low,n_loop,j)
          CASE(3)
             LAMBDA_3L(j)=LAMBDA_FROM_ALPHAS(as_start,segi%Q_low,n_loop,j)
          CASE(4)
             LAMBDA_4L(j)=LAMBDA_FROM_ALPHAS(as_start,segi%Q_low,n_loop,j)
          CASE(5)
             LAMBDA_5L(j)=LAMBDA_FROM_ALPHAS(as_start,segi%Q_low,n_loop,j)
          CASE DEFAULT
             PRINT *, "ERROR: > 5-loop as running is not allowed"
             STOP
          END SELECT
          ! I can exclude the alphas(Q=1E20), which will not be used
          !IF(j.LT.as_frame%n_high)segi%as_at_Q_high(n_loop)=ALPHAS_FROM_LAMBDA(segi%Q_high,n_loop,Q_THR)
          segi%as_at_Q_high(n_loop)=ALPHAS_FROM_LAMBDA(segi%Q_high,n_loop,Q_THR)
       ENDDO
    ENDDO

    ! the below segments then
    DO j=nseg-1,as_frame%n_low,-1
       Q_THR(1:6)=as_frame%Q_THR(1:6)
       Q_THR(j)=Q_THR(j)*0.99D0 ! in order to make sure Q_low and Q_high is still with nseg-th flavour
       IF(j.LT.6)Q_THR(j+1)=Q_THR(j+1)*1.01D0

       segi => as_frame%seg(j)
       ! nf -> nf-1 (decoupling relation)
       ! the as^(nf) at threshold muMatch_mQuark*m^(nf)
       DO n_loop=1,n_loop_max
          as_start=as_frame%seg(j+1)%as_at_Q_low(n_loop)
          aso2pi_start=as_start/TWOPI
          IF(as_decoupling_on)THEN
             SELECT CASE(n_loop)
             CASE(2)
                as_start=as_start*(1d0+segi%zeta_as_11*log_match*aso2pi_start)
             CASE(3)
                as_start=as_start*(1d0+segi%zeta_as_11*log_match*aso2pi_start&
                     +(segi%zeta_as_22*log_match**2+segi%zeta_as_21*log_match&
                     +segi%zeta_as_20)*aso2pi_start**2)
             CASE(4)
                as_start=as_start*(1d0+segi%zeta_as_11*log_match*aso2pi_start&
                     +(segi%zeta_as_22*log_match**2+segi%zeta_as_21*log_match&
                     +segi%zeta_as_20)*aso2pi_start**2+(segi%zeta_as_33*log_match**3&
                     +segi%zeta_as_32*log_match**2+segi%zeta_as_31*log_match&
                     +segi%zeta_as_30)*aso2pi_start**3)
             CASE(5)
                as_start=as_start*(1d0+segi%zeta_as_11*log_match*aso2pi_start&
                     +(segi%zeta_as_22*log_match**2+segi%zeta_as_21*log_match&
                     +segi%zeta_as_20)*aso2pi_start**2+(segi%zeta_as_33*log_match**3&
                     +segi%zeta_as_32*log_match**2+segi%zeta_as_31*log_match&
                     +segi%zeta_as_30)*aso2pi_start**3+(segi%zeta_as_44*log_match**4&
                     +segi%zeta_as_43*log_match**3+segi%zeta_as_42*log_match**2&
                     +segi%zeta_as_41*log_match+segi%zeta_as_40)*aso2pi_start**4)
             CASE(1)
                CONTINUE
             CASE DEFAULT
                PRINT *, "ERROR: only 5-loop beta expression (=> 4 loop decoupling relation) is allowed"
                STOP
             END SELECT
          ENDIF
          segi%as_at_Q_high(n_loop)=as_start
          SELECT CASE(n_loop)
          CASE(1)
             LAMBDA_1L(j)=LAMBDA_FROM_ALPHAS(as_start,segi%Q_high,n_loop,j)
          CASE(2)
             LAMBDA_2L(j)=LAMBDA_FROM_ALPHAS(as_start,segi%Q_high,n_loop,j)
          CASE(3)
             LAMBDA_3L(j)=LAMBDA_FROM_ALPHAS(as_start,segi%Q_high,n_loop,j)
          CASE(4)
             LAMBDA_4L(j)=LAMBDA_FROM_ALPHAS(as_start,segi%Q_high,n_loop,j)
          CASE(5)
             LAMBDA_5L(j)=LAMBDA_FROM_ALPHAS(as_start,segi%Q_high,n_loop,j)
          CASE DEFAULT
             PRINT *, "ERROR: > 5-loop as running is not allowed"
             STOP
          END SELECT
          ! I have to exclude the alphas(Q=0.5), which will not be used
          ! because both linear and quadratic Newton's procedure fails for 3- and 4- loop running with Q ~ 0.5 GeV
          IF(j.GT.as_frame%n_low)segi%as_at_Q_low(n_loop)=ALPHAS_FROM_LAMBDA(segi%Q_low,n_loop,Q_THR)
          !segi%as_at_Q_low(n_loop)=ALPHAS_FROM_LAMBDA(segi%Q_low,n_loop,Q_THR)
       ENDDO
    ENDDO

    ! restore the original nf setup
    CALL SET_QCD_NF(nf_orig)
    
    RETURN
  END SUBROUTINE ALPHAS_LAMBDA_INIT

!=============================================================================!
!=                                                                           =!
!=  ENDING OF special or approximated alpha_s running routines               =!
!=                                                                           =!
!=============================================================================!

END MODULE qcd_coupling
