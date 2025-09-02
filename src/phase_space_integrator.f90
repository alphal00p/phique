MODULE phase_space_integrator
  USE global
  USE plot_NNLO
  USE Constants
  IMPLICIT NONE
CONTAINS

  SUBROUTINE calc_NNLO
    USE qcd_constants
    USE qcd_setup
    USE qcd_coupling
    USE coupling_global
    USE ElasticPhotonPhotonFlux
    USE evaluate_couplings
    IMPLICIT NONE
    INTEGER::i,i1,j1,k1,s1x
    INTEGER::nh,innum,outnum
    INTEGER::ioerror,icase=0
    REAL(KIND(1d0)),DIMENSION(0:2)::wme
    REAL(KIND(1d0)),DIMENSION(3)::rslt
    INTEGER::itmxn,ncalln
    LOGICAL::lexist
    INTEGER::lfound
    REAL(KIND(1d0))::temp,muR,aS
    CHARACTER(len=7)::A1name,A2name
    REAL(KIND(1d0))::Aval,Zval

    include 'banner.inc'
    CALL ReadElem_integer('colpar',colpar)
    if(colpar.LT.1.OR.colpar.GT.4)then
       WRITE(*,*)"ERROR: colpar = ",colpar
       STOP
    ENDIF

    IF(colpar.EQ.1)THEN
       ! heavy ion is only possible when colpar=1
       CALL ReadElem_integer('UPC_photon_flux',UPC_photon_flux)
       IF(UPC_photon_flux.LE.0.OR.UPC_photon_flux.GT.3)THEN
          WRITE(*,*)"ERROR: UPC_Photon_flux =",UPC_photon_flux
          STOP
       ENDIF

       CALL ReadElem_integer('nuclearA_beam1',nuclearA1)
       CALL ReadElem_integer('nuclearA_beam2',nuclearA2)
       CALL ReadElem_integer('nuclearZ_beam1',nuclearZ1)
       CALL ReadElem_integer('nuclearZ_beam2',nuclearZ2)
       if(nuclearA1.EQ.1.and.nuclearZ1.EQ.1)THEN
          ! for proton
          nuclearA1=0
          nuclearZ1=0
       endif
       if(nuclearA2.EQ.1.and.nuclearZ2.EQ.1)THEN
          ! for proton
          nuclearA2=0
          nuclearZ2=0
       endif

       IF(UPC_photon_flux.EQ.3.AND.(nuclearA1.GT.1.OR.&
            nuclearA2.GT.1))THEN
          WRITE(*,*)"ERROR: UPC_photon_flux=3 (iww) only applies to proton-proton case"
          STOP
       ENDIF
       IF(.NOT.(nuclearA1.EQ.0.and.nuclearZ1.EQ.0))THEN
          A1name=GetASymbol(nuclearA1,nuclearZ1)
          CALL GetNuclearInfo(A1name,Aval,Zval,NNLO_RA(1),NNLO_aA(1),NNLO_wA(1))
       ENDIF
       IF(.NOT.(nuclearA2.EQ.0.and.nuclearZ2.EQ.0))THEN
          A2name=GetASymbol(nuclearA2,nuclearZ2)
          CALL GetNuclearInfo(A2name,Aval,Zval,NNLO_RA(2),NNLO_aA(2),NNLO_wA(2))
       ENDIF
    ENDIF

    CALL ReadElem_real('energy_beam1',energy_beam1)
    CALL ReadElem_real('energy_beam2',energy_beam2)
    IF(energy_beam1.LE.0d0.OR.energy_beam2.LE.0d0)THEN
       WRITE(*,*)"ERROR: the energies of the beams cannot be non positive"
       STOP
    ENDIF

    CALL ReadElem_integer('quark',quark)
    IF(quark.LT.4.OR.quark.GT.6)THEN
       WRITE(*,*)"ERROR: quark must be between 4 to 6. quark=",quark
       STOP
    ENDIF
    ! number of massless quark flavors
    nq=quark-1
    IF(quark.EQ.4.OR.quark.EQ.6)THEN
       Q_charge=2
    ELSE
       Q_charge=-1
    ENDIF

    CALL ReadElem_integer('order',order)
    IF(order.LT.0.OR.order.GT.2)then
       WRITE(*,*)"ERROR: order must be between 0 to 2. order=",order
       STOP
    ENDIF

    CALL ReadElem_integer('coulomb',coulomb)
    IF(coulomb.LT.0.OR.coulomb.GT.2)THEN
       WRITE(*,*)"ERROR: coulomb must be between 0 to 2. coulomb=",coulomb
       STOP
    ENDIF

    CALL ReadElem_logic('topdrawer_output',topdrawer_output)
    CALL ReadElem_logic('gnuplot_output',gnuplot_output)
    CALL ReadElem_logic('root_output',root_output)
    CALL ReadElem_logic('hwu_output',hwu_output)
    plot_output=topdrawer_output.OR.gnuplot_output.OR.root_output.OR.hwu_output

    CALL ReadConst
    IF(alphasMZ.LT.0d0)THEN
       alphasMZ=0.118d0
    ENDIF
    IF(abs(order).GE.1)THEN
       ! Need QCD corrections
       ! set QCD colour factors
       CALL SET_QCD_GROUP(QCD_CA,QCD_CF,QCD_TF)
       ! initialisation for alpha_s running
       CALL ReadElem_integer('alphas_nloop',alphas_nloop)
       IF(alphas_nloop.LT.1.OR.alphas_nloop.GT.5)THEN
          WRITE(*,*)"ERROR: the alphas RG running can only be from 1 to 5"
          STOP
       ENDIF
       ! Q_THR_default is the OS mass (2)
       ! the threshold is 1d0*quark_mass (1d0)
       CALL as_run_init(as_box,alphasMZ,zmass_PDG,alphas_nloop,&
            no_fix_nf,Q_THR_default(4:6),2,1D0)
    ENDIF   

    IF(nunit1.NE.6)THEN
       OPEN(UNIT=nunit1,FILE=TRIM(output_dir)//"RESULT_NNLO.out")
    ENDIF

    ! Scale scheme
    CALL ReadElem_integer('Scale',scale)
    IF(scale.NE.0.AND.scale.NE.1)THEN
       WRITE(*,*)"ERROR: only 0 or 1 is allowed for Scale. Scale = ", scale
       STOP
    ENDIF
    CALL ReadElem_real('FScaleValue',fscalevalue)
    CALL ReadElem_real('muR_over_ref',muR_over_ref)

    ! reweighting stuff
    CALL ReadElem_logic('reweight_Scale',reweight_scale)
    IF(order.EQ.0)THEN
       ! no need to reweight the scale for LO
       reweight_scale=.FALSE.
    ENDIF
    IF(reweight_scale)THEN
       CALL ReadElem_real('rw_RScale_down',rw_Rscale_down)
       CALL ReadElem_real('rw_RScale_up',rw_Rscale_up)
       ho_nscale=2
       IF(rw_Rscale_down.LE.0d0.OR.rw_Rscale_up.LE.0d0)THEN
          reweight_scale=.FALSE.
          ho_nscale=0
          WRITE(*,*)"WARNING:Scale reweighting is off."
          WRITE(*,*)"WARNING:Please make sure rw_Rscale_up, rw_Rscale_down are larger than 0."
       ENDIF
       IF(reweight_scale.AND.rw_Rscale_down.GT.rw_Rscale_up)THEN
          temp=rw_Rscale_down
          rw_Rscale_down=rw_Rscale_up
          rw_Rscale_up=temp
       ENDIF
       rw_Rscales(1)=1d0
       rw_Rscales(2)=rw_Rscale_up
       rw_Rscales(3)=rw_Rscale_down
    ENDIF
    IF(reweight_scale)THEN
       tot_nrwgt=tot_nrwgt+ho_nscale+1
    ENDIF

    CALL ReadElem_integer("nmc",nmc)

    Wgaga=2d0*DSQRT(energy_beam1*energy_beam2)
    IF(Wgaga.LE.2d0*MQ)THEN
       rslt(1:3)=0d0
    ELSE
       ! we start to read the NLOxs.grid
       OPEN(23331,FILE=TRIM(grid_dir)//"NLOxs.grid",ACTION="read")
       READ(23331,*)n_NLO_energies
       ALLOCATE(NLO_energies(n_NLO_energies))
       ! In each row, there are LO, NLO, DNLO and NLO/LO
       ! cross sections are in unit of pb
       ALLOCATE(total_NLO_xs(n_NLO_energies,4))
       IE_NLO_350=0
       IE_NLO_500=0
       IE_NLO_1000=0
       IE_NLO_5000=0
       IE_NLO_10000=0
       DO i=1,n_NLO_energies
          READ(23331,*)NLO_energies(i),total_NLO_xs(i,1),total_NLO_xs(i,2)
          IF(IE_NLO_350.EQ.0.AND.NLO_energies(i).GT.350d0)THEN
             IE_NLO_350=i
          ENDIF
          IF(IE_NLO_500.EQ.0.AND.NLO_energies(i).GT.500d0)THEN
             IE_NLO_500=i
          ENDIF
          IF(IE_NLO_1000.EQ.0.AND.NLO_energies(i).GT.1000d0)THEN
             IE_NLO_1000=i
          ENDIF
          IF(IE_NLO_5000.EQ.0.AND.NLO_energies(i).GT.5000d0)THEN
             IE_NLO_5000=i
          ENDIF
          IF(IE_NLO_10000.EQ.0.AND.NLO_energies(i).GT.10000d0)THEN
             IE_NLO_10000=i
          ENDIF
          total_NLO_xs(i,3)=total_NLO_xs(i,2)-total_NLO_xs(i,1)
          total_NLO_xs(i,4)=total_NLO_xs(i,2)/total_NLO_xs(i,1)
       ENDDO
       IE_NLO_350=MAX(1,IE_NLO_350)
       IE_NLO_500=MAX(1,IE_NLO_500)
       IE_NLO_1000=MAX(1,IE_NLO_1000)
       IE_NLO_5000=MAX(1,IE_NLO_5000)
       IE_NLO_10000=MAX(1,IE_NLO_10000)
       CLOSE(23331)
       ! we start to read the DNNLOxs.grid
       ! with as=0.118 (alphasMZ_default) and
       ! alpha=1/137 (1/alphaemm1_default), MQ=173 GeV (MQ_default)
       ! NLO has no explict mu and nq dependence
       OPEN(23332,FILE=TRIM(grid_dir)//"DNNLOxs.grid",ACTION="read")
       READ(23332,*)n_NNLO_energies
       ALLOCATE(NNLO_energies(n_NNLO_energies))
       ! In each row, there are DNNLO_nonnf/LO, DNNLO_nf/LO in unit of pb
       ! DNNLO = LO*(DNNLO_nonnf+nq*DNNLO_nf)
       ! at scale 91.188 GeV (scale_default)
       ! and as = 0.118 (alphasMZ_default),
       ! alpha=1/137 (1/alphaemm1_default), MQ=173 GeV (MQ_default)
       ALLOCATE(total_DNNLO_xs(n_NNLO_energies,2))
       ALLOCATE(total_DNNLO_xserr(n_NNLO_energies,2))
       IE_NNLO_350=0
       IE_NNLO_500=0
       IE_NNLO_1000=0
       DO i=1,n_NNLO_energies
          READ(23332,*)NNLO_energies(i),total_DNNLO_xs(i,1),&
               total_DNNLO_xserr(i,1),total_DNNLO_xs(i,2),total_DNNLO_xserr(i,2)
          IF(IE_NNLO_350.EQ.0.AND.NNLO_energies(i).GT.350d0)THEN
             IE_NNLO_350=i
          ENDIF
          IF(IE_NNLO_500.EQ.0.AND.NNLO_energies(i).GT.500d0)THEN
             IE_NNLO_500=i
          ENDIF
          IF(IE_NNLO_1000.EQ.0.AND.NNLO_energies(i).GT.1000d0)THEN
             IE_NNLO_1000=i
          ENDIF
       ENDDO
       IE_NNLO_350=MAX(1,IE_NNLO_350)
       IE_NNLO_500=MAX(1,IE_NNLO_500)
       IE_NNLO_1000=MAX(1,IE_NNLO_1000)
       CLOSE(23332)
       
       ! colpar=4 should be treated separately since there is no integration
       IF(colpar.LE.3)THEN
          CALL eval_xs_VEGAS(rslt,nmc)
          WRITE(nunit1,*)"sigma (pb)                   sd (pb)"
          WRITE(nunit1,*)rslt(1),rslt(2)
       ELSE
          CALL NNLO_scale(muR)
          aS=ALPHAS(muR)
          CALL GetME2(aS,Wgaga,muR,wme)
          WRITE(nunit1,*)"LO sigma (pb) : ", wme(0)
          IF(coulomb.gt.0.and.reweight_scale)then
             WRITE(nunit1,*)" LO scale var (%) : ", (LO_wgtxsecmu(3)-1d0)*100,&
                  (LO_wgtxsecmu(2)-1d0)*100
          endif
          IF(order.GE.1)THEN
             WRITE(nunit1,*)"NLO sigma (pb) : ", wme(1)
             IF(reweight_scale)THEN
                WRITE(nunit1,*)"NLO scale var (%) : ", (NLO_wgtxsecmu(3)-1d0)*100,&
                     (NLO_wgtxsecmu(2)-1d0)*100
             ENDIF
             IF(order.GE.2)THEN
                WRITE(nunit1,*)"NNLO sigma (pb) : ", wme(2)
                IF(reweight_scale)THEN
                   WRITE(nunit1,*)"NNLO scale var (%) : ", (NNLO_wgtxsecmu(3)-1d0)*100,&
                        (NNLO_wgtxsecmu(2)-1d0)*100
                ENDIF
             ENDIF
          ENDIF
       ENDIF
    ENDIF

    RETURN
    
  END SUBROUTINE calc_NNLO

  SUBROUTINE eval_xs_VEGAS(rslt,ncalln)
    USE global
    USE MC_VEGAS
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(3),INTENT(OUT)::rslt
    INTEGER,INTENT(IN),OPTIONAL::ncalln
    REAL(KIND(1d0))::vfes,sd,chi2a
    INTEGER::ncalm,nc,ii
    CHARACTER(len=4),DIMENSION(20)::chchar
    lwmax=0
    NPRN=-1
    varnum=2
    DO ii=1,varnum
       XL(ii)=0.0d0
       XU(ii)=1.0d0
    ENDDO
    ITMX=1
    IF(PRESENT(ncalln))THEN
       ncalm=ncalln
    ELSE
       ncalm=5000
    ENDIF
    chchar(1)="40K"
    chchar(2)="80K"
    chchar(3)="160K"
    chchar(4)="320K"
    chchar(5)="640K"
    chchar(6)="1M"
    chchar(7)="2M"
    chchar(8)="4M"
    chchar(9)="8M"
    chchar(10)="16M"
    chchar(11)='32M'
    chchar(12)='64M'
    chchar(13)='120M'
    chchar(14)='240M'
    chchar(15)='480M'
    chchar(16)='960M'
    chchar(17)='2G'
    chchar(18)='4G'
    chchar(19)='8G'
    chchar(20)='16G'
    nprint=10000
    NCALL=20000
    IF(plot_output)CALL initplot_NNLO
    CALL VEGAS(varnum,eval_xs_fxn,vfes,sd,chi2a)
    IF(plot_output)CALL plotout_NNLO
    WRITE(*,*)' '
    WRITE(*,*)"====================NCALL=20K==========================="
    WRITE(*,*)" "
    ii=1
    WRITE(*,*)"ITERATION ",ii,":"
    WRITE(*,*)vfes,"+\-",sd
    WRITE(*,*)"precision:",sd/vfes
    DO ii=2,10
       IF(plot_output)CALL initplot_NNLO
       CALL VEGAS(varnum,eval_xs_fxn,vfes,sd,chi2a,1)
       IF(plot_output)CALL plotout_NNLO
       WRITE(*,*)"ITERATION ",ii,":"
       WRITE(*,*)vfes,"+\-",sd
       WRITE(*,*)"precision:",sd/vfes
    ENDDO
    ii=1
    DO
       nc=2*NCALL
       IF(nc.GT.ncalm)EXIT
       !IF(2*nc.GT.ncalm.AND.lunwei2.AND.lwmax.EQ.0)THEN
       !   lwmax=1
       !ENDIF
       NCALL=nc
       IF(NCALL/maxprint.GT.nprint)nprint=NCALL/maxprint
       WRITE(*,*)"====================NCALL="//chchar(ii)//"==========================="
       WRITE(*,*)" "
       IF(plot_output)CALL initplot_NNLO
       CALL VEGAS(varnum,eval_xs_fxn,vfes,sd,chi2a,1)
       IF(plot_output)CALL plotout_NNLO
       WRITE(*,*)vfes,"+\-",sd
       WRITE(*,*)"precision:",sd/vfes
       ii=ii+1
    ENDDO
    !IDBMUP1=2212
    !IDBMUP2=2212
    !NPRUP=NPRUP+1
    !WRITE(200,5100) IDBMUP1,IDBMUP2,EBMUP1,EBMUP2,&
    !     iPDFGUP1,iPDFGUP2,iPDFSUP1,iPDFSUP2,IDWTUP,NPRUP
    !WRITE(200,5200) vfes,sd,1d0, 92
    rslt(1)=vfes
    rslt(2)=sd
    rslt(3)=chi2a
    RETURN
5100 FORMAT(1P,2I8,2E14.6,6I8)
5200 FORMAT(1P,3E20.10,I6)
  END SUBROUTINE eval_xs_VEGAS

  FUNCTION eval_xs_fxn(x,wgt)
    USE global
    USE OpticalGlauber_Geometry
    USE ElasticPhotonPhotonFlux
    USE evaluate_couplings
    IMPLICIT NONE
    include '../vendor/gammaUPC/run90.inc'
    REAL(KIND(1d0))::eval_xs_fxn
    REAL(KIND(1d0)),DIMENSION(varnum),INTENT(IN)::x
    REAL(KIND(1d0)),INTENT(IN)::wgt
    INTEGER::init=0,nwarn,nwmax,nwri,nnn=0,nnntot=0
    SAVE init,nwarn,nwmax,nwri,nnn,nnntot
    INTEGER::i,j,ii,kk
    LOGICAL::icuts,lflag
    REAL(KIND(1d0))::sqs,sq,ycollcm,maa,jac,ycms
    SAVE sqs,sq,ycollcm
    REAL(KIND(1d0))::taumin,taumax,ymax0,ymin0,jac0,tau10,tau11,tau
    SAVE taumin,taumax,jac0,tau10,tau11
    !REAL(KIND(1d0)),PARAMETER::pipi=3.14159265358979323846264338328d0
    !REAL(KIND(1d0)),DIMENSION(0:3)::PBOO
    REAL(KIND(1d0))::wsf,wthis,w1
    REAL(KIND(1d0)),DIMENSION(0:2)::wtotal,wme
    REAL(KIND(1d0))::recmax=0,recmin=0
    REAL(KIND(1d0))::Q1,Q2,Q12,Q22,energy
    LOGICAL::labeqcoll
    SAVE labeqcoll
    REAL(KIND(1d0))::energy_min,energy_max
    SAVE energy_min,energy_max
    REAL(KIND(1d0))::eQ
    SAVE eQ
    REAL(KIND(1d0))::muR,aS
    IF(init.EQ.0)THEN
       nwmax=0
       nwarn=0
       nwri=0
       ebeam0(1)=ABS(energy_beam1)
       ebeam0(2)=ABS(energy_beam2)
       IF(ABS(ebeam0(1)-ebeam0(2))/MAX(ebeam0(1)+ebeam0(2),1d-17).LT.1d-8)THEN
          labeqcoll=.TRUE.
       ELSE
          labeqcoll=.FALSE.
       ENDIF
       sqs=2d0*DSQRT(ebeam0(1)*ebeam0(2)) ! we always neglect the mass of initial states
       ycollcm=DLOG(ABS(ebeam0(1))/ABS(ebeam0(2)))/2d0
       
       xp1=1d0
       xp2=1d0
       sq=sqs*sqs

       energy_min=2d0*MQ
       energy_max=sqs
       
       taumin=energy_min**2/sq
       taumax=MIN(1d0,energy_max**2/sq)
       CALL GET_ntau(taumin,ntau)
       if(ntau.LE.1d0)then
          write(*,*)"ERROR: ntau <= 1"
          stop
       endif
       tau10=taumin**(1d0-ntau)
       tau11=tau10-taumax**(1d0-ntau)
       jac0=tau11/(ntau-1d0)

       init=1
    ENDIF
    
    nnntot=nnntot+1
    wtotal(0:2)=0d0
    eval_xs_fxn=0d0
    jac=jac0
    ! tau=x1*x2, y=1/2*Log(x1/x2)
    ! x1=Sqrt(tau)*Exp(y), x2=Sqrt(tau)*Exp(-y)
    ! dx1 from taumin to taumax
    ! dx2 from taumin/x1 to taumax/x1
    ! is same as
    ! dtau from taumin to taumax
    ! dy from 1/2*Log(tau) to -1/2*Log(tau)
    ! with Jacobi unity
    ! with PDF
    tau=(tau10-tau11*x(1))**(1d0/(1d0-ntau))
    if(tau.GT.1d0)tau=1d0
    if(tau.LE.0d0)then
       write(*,*)"Error: tau <= 0"
       stop
    endif
    ycms=(1d0-2d0*x(2))/2d0*DLOG(tau) ! ycms is defined in the center-of-mass frame of the two intial beams
    xp1=DSQRT(tau)*DEXP(ycms)
    xp2=DSQRT(tau)*DEXP(-ycms)
    jac=jac*tau**(ntau)*(-DLOG(tau))
    
    ! partonic energy
    energy=DSQRT(sq*xp1*xp2)
    Wgaga=energy

    CALL NNLO_scale(muR)
    aS=ALPHAS(muR)

    CALL GetME2(aS,energy,muR,wme(0:2))
    IF(wme(0).EQ.0d0)THEN
       eval_xs_fxn=0d0
       return
    endif

    ! for rapidity observables one should change the bin by +ycms+ycollcm
    CALL strf_pdf(wsf)
    IF(wsf.LE.0d0)THEN
       eval_xs_fxn=0d0
       return
    endif
    wtotal(0)=wme(0)*wsf*jac
    wtotal(1)=wme(1)*wsf*jac
    wtotal(2)=wme(2)*wsf*jac
    
    if(plot_output.and.wtotal(0).ne.0d0)then
       w1=wtotal(0)*wgt
       call outfun_NNLO(0,w1)
       IF(order.GE.1)THEN
          w1=wtotal(1)*wgt
          call outfun_NNLO(1,w1)
          IF(order.GE.2)THEN
             w1=wtotal(2)*wgt
             call outfun_NNLO(2,w1)
          ENDIF
       ENDIF
    endif
    
    eval_xs_fxn=wtotal(order)
    nnn=nnn+1
    IF(recmax.LT.eval_xs_fxn)recmax=eval_xs_fxn
    IF(recmin.GT.eval_xs_fxn)recmin=eval_xs_fxn
    IF(MOD(nnn,nprint).EQ.0)THEN
       PRINT *,"max=",recmax
       PRINT *,"min=",recmin
       PRINT *, "      n_pass,     n_total"
       PRINT *,nnn,nnntot
    ENDIF
    RETURN
  END FUNCTION eval_xs_fxn

  SUBROUTINE GetME2(aS,energy,muR,wme)
    USE interpolation
    USE evaluate_couplings
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::aS,muR,energy
    REAL(KIND(1d0)),DIMENSION(0:2),INTENT(OUT)::wme
    INTEGER::init=0
    SAVE init
    REAL(KIND(1d0))::eQ
    SAVE eQ
    REAL(KIND(1d0))::energy0,LOxs,DNLOxs,NLOxs,DNLOxsCoul,DNLOnonCoulxs,NLOKfactor
    REAL(KIND(1d0))::DNNLOxsnf,DNNLOxsnonnf,DNNLOxs,NNLOxs
    REAL(KIND(1d0))::NLO_energy_min,NLO_energy_max
    SAVE NLO_energy_min, NLO_energy_max
    REAL(KIND(1d0))::NNLO_energy_Coul,NNLO_energy_HE
    SAVE NNLO_energy_Coul,NNLO_energy_HE
    REAL(KIND(1d0)),PARAMETER::twopi=6.28318530717958647692528676656d0
    REAL(KIND(1d0))::LOxsup,NLOxsup,NNLOxsup
    REAL(KIND(1d0))::LOxsdown,NLOxsdown,NNLOxsdown
    INTEGER,PARAMETER::n_interp=10
    REAL(KIND(1d0)),DIMENSION(n_interp)::energies_interp
    REAL(KIND(1d0)),DIMENSION(n_interp)::total_xs1_interp,total_xs2_interp
    INTEGER::i,j,k
    REAL(KIND(1d0))::DResxs,DResxsup,DResxsdown
    REAL(KIND(1d0))::aSCoul,aSCoulup,aSCouldown,aSMEup,aSMEdown
    REAL(KIND(1d0))::scaleCoul,scaleCoulup,scaleCouldown
    REAL(KIND(1d0))::muRup,muRdown
    REAL(KIND(1d0))::DNLOxsNonCoul,DNLOxsNonCoulup,DNLOxsNonCouldown
    REAL(KIND(1d0)),PARAMETER::scale_min=1d0 ! cannot below 1 GeV
    IF(init.EQ.0)THEN
       ! heavy quark related
       eQ=DBLE(Q_charge)/3d0
       NLO_energy_min=NLO_energies(1)
       !NLO_energy_max=NLO_energies(n_NLO_energies)
       !NLO_energy_min=346.05d0
       !NLO_energy_min=346.002d0
       !NLO_energy_min=346.3d0
       NLO_energy_max=NLO_energies(n_NLO_energies)
       !NNLO_energy_Coul=346.05d0
       NNLO_energy_Coul=346.3d0
       NNLO_energy_HE=1150d0
       init=1
    ENDIF
    IF(energy.LE.2d0*MQ)THEN
       wme(0:2)=0d0
       RETURN
    ENDIF

    LO_wgtxsecmu(1:3)=1d0
    NLO_wgtxsecmu(1:3)=1d0
    NNLO_wgtxsecmu(1:3)=1d0
    IF(reweight_scale)THEN
       muRup=MAX(muR*rw_Rscales(2),scale_min)
       muRdown=MAX(muR*rw_Rscales(3),scale_min)
    ENDIF

    ! can also use eQ=2/3 to obtain LO xs (in unit of pb)
    CALL Get_LOxs(energy,MQ,eQ,LOxs)
    wme(0)=LOxs
    IF(coulomb.gt.0)THEN
       CALL Get_CoulRes_Scale(energy,MQ,muR,scaleCoul)
       aSCoul=ALPHAS(scaleCoul)
       IF(coulomb.eq.1)THEN
          ! LP
          CALL Get_LPxs_CoulRes(0,energy,MQ,aSCoul,aS,LOxs,DResxs)
       ELSEIF(coulomb.eq.2)THEN
          ! NLP
          CALL Get_NLPxs_CoulRes(0,energy,MQ,scaleCoul,muR,&
               aSCoul,aS,LOxs,0d0,DResxs)
       ELSE
          WRITE(*,*)"ERROR: Cannot reach here #1"
          STOP
       ENDIF
       wme(0)=LOxs+DResxs
       IF(reweight_scale)THEN
          !CALL Get_CoulRes_Scale(energy,MQ,muRup,scaleCoulup)
          !CALL Get_CoulRes_Scale(energy,MQ,muRdown,scaleCouldown)
          scaleCoulup=MAX(scaleCoul*rw_Rscales(2),scale_min)
          scaleCouldown=MAX(scaleCoul*rw_Rscales(3),scale_min)
          aSCoulup=ALPHAS(scaleCoulup)
          aSCouldown=ALPHAS(scaleCouldown)
          aSMEup=ALPHAS(muRup)
          aSMEdown=ALPHAS(muRdown)
          IF(wme(0).NE.0d0)THEN
             IF(coulomb.eq.1)THEN
                ! LP
                CALL Get_LPxs_CoulRes(0,energy,MQ,aSCoulup,aSMEup,LOxs,DResxsup)
                CALL Get_LPxs_CoulRes(0,energy,MQ,aSCouldown,aSMEdown,LOxs,DResxsdown)
             ELSEIF(coulomb.eq.2)THEN
                ! NLP
                CALL Get_NLPxs_CoulRes(0,energy,MQ,scaleCoulup,muRup,&
                     aSCoulup,aSMEup,LOxs,0d0,DResxsup)
                CALL Get_NLPxs_CoulRes(0,energy,MQ,scaleCouldown,muRdown,&
                     aSCouldown,aSMEdown,LOxs,0d0,DResxsdown)
             ELSE
                WRITE(*,*)"ERROR: cannot reach here #2"
                STOP
             ENDIF
             LO_wgtxsecmu(2)=(LOxs+DResxsup)/wme(0)
             LO_wgtxsecmu(3)=(LOxs+DResxsdown)/wme(0)
          ENDIF
       ENDIF
    ENDIF
    IF(wme(0).EQ.0d0)RETURN
    
    IF(order.EQ.0)RETURN

    ! rescaled energy (to meet the grid with MQ=173 GeV)
    energy0=energy*MQ_default/MQ
    
    
    IF(energy0.LT.NLO_energy_min)THEN
       ! we should use the Coulomb approx
       CALL Get_DNLOxs_Coul(energy,MQ,aS,LOxs,DNLOxs)
    ELSEIF(energy0.GT.NLO_energy_max)THEN
       ! we should use HE fit extrapolation
       CALL Get_DNLOxs_HE(energy,MQ,aS,LOxs,DNLOxs)
    ELSE
       ! The spline approximation reference:
       !\bibitem{dierckx1995curve}
       !P.~Dierckx,{\em Curve and Surface Fitting with Splines}.\newblock Monographs on numerical analysis.Clarendon Press,1995.
       ! we just use interpolation
       !CALL SPLINE_INTERPOLATE(NLO_energies(1:n_NLO_energies),&
       !     total_NLO_xs(1:n_NLO_energies,4),n_NLO_energies,&
       !     energy0,NLOKfactor)
       j=1
       IF(energy0.GT.350d0.AND.energy0.LE.500d0)THEN
          j=IE_NLO_350
       ELSEIF(energy0.GT.500d0.AND.energy0.LE.1000d0)THEN
          j=IE_NLO_500
       ELSEIF(energy0.GT.1000d0.AND.energy0.LE.5000d0)THEN
          j=IE_NLO_1000
       ELSEIF(energy0.GT.5000d0.AND.energy0.LE.10000d0)THEN
          j=IE_NLO_5000
       ELSEIF(energy0.GT.10000d0)THEN
          j=IE_NLO_10000
       ENDIF
       k=0       
       DO i=j,n_NLO_energies
          IF(NLO_energies(i).GE.energy0)THEN
             k=i
             EXIT
          ENDIF
          IF(i.EQ.n_NLO_energies)THEN
             k=i+1
          ENDIF
       ENDDO
       k=k-n_interp/2-1
       IF(k.LT.0)THEN
          k=0
       ELSEIF(k+n_interp.GT.n_NLO_energies)THEN
          k=n_NLO_energies-n_interp
       ENDIF
       DO i=1,n_interp
          energies_interp(i)=NLO_energies(k+i)
          total_xs1_interp(i)=total_NLO_xs(k+i,4)
       ENDDO
       CALL SPLINE_INTERPOLATE(energies_interp,&
            total_xs1_interp,n_interp,&
            energy0,NLOKfactor)
       DNLOxs=LOxs*(NLOKfactor-1d0)*aS/alphasMZ_default
    ENDIF
    NLOxs=DNLOxs+LOxs
    wme(1)=NLOxs
    IF(coulomb.gt.0)THEN
       IF(coulomb.eq.1)THEN
          ! LP
          CALL Get_LPxs_CoulRes(1,energy,MQ,aSCoul,aS,LOxs,DResxs)
       ELSEIF(coulomb.eq.2)THEN
          ! NLP
          CALL Get_DNLOxs_Coul(energy,MQ,aS,LOxs,DNLOxsCoul)
          DNLOxsNonCoul=DNLOxs-DNLOxsCoul
          CALL Get_NLPxs_CoulRes(1,energy,MQ,scaleCoul,muR,&
               aSCoul,aS,LOxs,DNLOxsNonCoul,DResxs)
       ELSE
          WRITE(*,*)"ERROR: cannot reach here #3"
          STOP
       ENDIF
       wme(1)=NLOxs+DResxs
    ENDIF

    IF(order.EQ.1)THEN
       IF(reweight_scale)THEN
          LOxsup=LOxs
          NLOxsup=NLOxs
          CALL UPDATE_NLO_SCALE_ALPHAS(LOxsup,NLOxsup,muR,aS,muRup)
          LOxsdown=LOxs
          NLOxsdown=NLOxs
          CALL UPDATE_NLO_SCALE_ALPHAS(LOxsup,NLOxsup,muR,aS,muRdown)
          IF(coulomb.eq.0)THEN
             IF(NLOxs.NE.0d0)THEN
                NLO_wgtxsecmu(2)=NLOxsup/NLOxs
                NLO_wgtxsecmu(3)=NLOxsdown/NLOxs
             ENDIF
          ELSE
             IF(wme(1).NE.0d0)THEN
                IF(coulomb.eq.1)THEN
                   ! LP
                   CALL Get_LPxs_CoulRes(1,energy,MQ,aSCoulup,aSMEup,LOxs,DResxsup)
                   CALL Get_LPxs_CoulRes(1,energy,MQ,aSCouldown,aSMEdown,LOxs,DResxsdown)
                ELSEIF(coulomb.eq.2)THEN
                   ! NLP
                   DNLOxsNonCoulup=DNLOxsNonCoul*aSMEup/aS
                   DNLOxsNonCouldown=DNLOxsNonCoul*aSMEdown/aS
                   CALL Get_NLPxs_CoulRes(1,energy,MQ,scaleCoulup,muRup,&
                        aSCoulup,aSMEup,LOxs,DNLOxsNonCoulup,DResxsup)
                   CALL Get_NLPxs_CoulRes(1,energy,MQ,scaleCouldown,muRdown,&
                        aSCouldown,aSMEdown,LOxs,DNLOxsNonCouldown,DResxsdown)
                ELSE
                   WRITE(*,*)"ERROR: cannot reach here #3"
                   STOP
                ENDIF
                NLO_wgtxsecmu(2)=(NLOxsup+DResxsup)/wme(1)
                NLO_wgtxsecmu(3)=(NLOxsdown+DResxsdown)/wme(1)
             ENDIF
          ENDIF
       ENDIF
       RETURN
    ENDIF

    IF(energy0.LT.NNLO_energy_Coul)THEN
       ! we should use the Coulomb approx
       CALL Get_DNLOxs_Coul(energy,MQ,aS,LOxs,DNLOxsCoul)
       DNLOnonCoulxs=DNLOxs-DNLOxsCoul
       CALL Get_DNNLOxs_Coul(nq,energy,MQ,muR,aS,LOxs,DNLOnonCoulxs,DNNLOxs)
    ELSEIF(energy0.GT.NNLO_energy_HE)THEN
       ! we should use HE fit extrapolation
       CALL Get_DNNLOnqxs_HE(energy,MQ,muR,aS,LOxs,DNLOxs,DNNLOxsnf)
       CALL Get_DNNLOnonnqxs_HE(energy,MQ,muR,aS,LOxs,DNLOxs,DNNLOxsnonnf)
       DNNLOxs=DNNLOxsnonnf+nq*DNNLOxsnf
    ELSE
       ! we just use interpolation
       !CALL SPLINE_INTERPOLATE(NNLO_energies(1:n_NNLO_energies),&
       !     total_DNNLO_xs(1:n_NNLO_energies,1),n_NNLO_energies,&
       !     energy0,DNNLOxsnonnf)
       !CALL SPLINE_INTERPOLATE(NNLO_energies(1:n_NNLO_energies),&
       !     total_DNNLO_xs(1:n_NNLO_energies,2),n_NNLO_energies,&
       !     energy0,DNNLOxsnf)
       j=1
       IF(energy0.GT.350d0.AND.energy0.LE.500d0)THEN
          j=IE_NNLO_350
       ELSEIF(energy0.GT.500d0.AND.energy0.LE.1000d0)THEN
          j=IE_NNLO_500
       ELSEIF(energy0.GT.1000d0)THEN
          j=IE_NNLO_1000
       ENDIF
       k=0
       DO i=j,n_NNLO_energies
          IF(NNLO_energies(i).GE.energy0)THEN
             k=i
             EXIT
          ENDIF
          IF(i.EQ.n_NNLO_energies)THEN
             k=i+1
          ENDIF
       ENDDO
       k=k-n_interp/2-1
       IF(k.LT.0)THEN
          k=0
       ELSEIF(k+n_interp.GT.n_NNLO_energies)THEN
          k=n_NNLO_energies-n_interp
       ENDIF
       DO i=1,n_interp
          energies_interp(i)=NNLO_energies(k+i)
          total_xs1_interp(i)=total_DNNLO_xs(k+i,1)
          total_xs2_interp(i)=total_DNNLO_xs(k+i,2)
       ENDDO
       CALL SPLINE_INTERPOLATE(energies_interp,&
            total_xs1_interp,n_interp,&
            energy0,DNNLOxsnonnf)
       CALL SPLINE_INTERPOLATE(energies_interp,&
            total_xs2_interp,n_interp,&
            energy0,DNNLOxsnf)
       DNNLOxsnonnf=LOxs*DNNLOxsnonnf*(aS/alphasMZ_default)**2&
            +11d0/2d0*DNLOxs*(aS/twopi)*DLOG(muR**2/scale_default**2&
            *MQ_default**2/MQ**2)
       DNNLOxsnf=LOxs*DNNLOxsnf*(aS/alphasMZ_default)**2&
            -1d0/3d0*DNLOxs*(aS/twopi)*DLOG(muR**2/scale_default**2&
            *MQ_default**2/MQ**2)
       DNNLOxs=DNNLOxsnonnf+nq*DNNLOxsnf
    ENDIF
    NNLOxs=DNNLOxs+NLOxs
    wme(2)=NNLOxs
    IF(coulomb.gt.0)THEN
       IF(coulomb.eq.1)THEN
          ! LP
          CALL Get_LPxs_CoulRes(2,energy,MQ,aSCoul,aS,LOxs,DResxs)
       ELSEIF(coulomb.eq.2)THEN
          ! NLP
          CALL Get_NLPxs_CoulRes(2,energy,MQ,scaleCoul,muR,&
               aSCoul,aS,LOxs,DNLOxsNonCoul,DResxs)
       ELSE
          WRITE(*,*)"ERROR: cannot reach here #4"
          STOP
       ENDIF
       wme(2)=NNLOxs+DResxs
    ENDIF

    IF(reweight_scale)THEN
       LOxsup=LOxs
       NLOxsup=NLOxs
       NNLOxsup=NNLOxs
       CALL UPDATE_NNLO_SCALE_ALPHAS(LOxsup,NLOxsup,NNLOxsup,muR,aS,muRup)
       LOxsdown=LOxs
       NLOxsdown=NLOxs
       NNLOxsdown=NNLOxs
       CALL UPDATE_NNLO_SCALE_ALPHAS(LOxsdown,NLOxsdown,NNLOxsdown,muR,aS,muRdown)
       IF(coulomb.eq.0)THEN
          IF(NLOxs.NE.0d0)THEN
             NLO_wgtxsecmu(2)=NLOxsup/NLOxs
             NLO_wgtxsecmu(3)=NLOxsdown/NLOxs
          ENDIF
          IF(NNLOxs.NE.0d0)THEN
             NNLO_wgtxsecmu(2)=NNLOxsup/NNLOxs
             NNLO_wgtxsecmu(3)=NNLOxsdown/NNLOxs
          ENDIF
       ELSE
          IF(coulomb.EQ.2.AND.(wme(1).NE.0d0.OR.wme(2).NE.0d0))THEN
             DNLOxsNonCoulup=DNLOxsNonCoul*aSMEup/aS
             DNLOxsNonCouldown=DNLOxsNonCoul*aSMEdown/aS
          ENDIF
          IF(wme(1).NE.0d0)THEN
             IF(coulomb.EQ.1)THEN
                ! LP
                CALL Get_LPxs_CoulRes(1,energy,MQ,aSCoulup,aSMEup,LOxs,DResxsup)
                CALL Get_LPxs_CoulRes(1,energy,MQ,aSCouldown,aSMEdown,LOxs,DResxsdown)
             ELSEIF(coulomb.eq.2)THEN
                ! NLP
                CALL Get_NLPxs_CoulRes(1,energy,MQ,scaleCoulup,muRup,&
                     aSCoulup,aSMEup,LOxs,DNLOxsNonCoulup,DResxsup)
                CALL Get_NLPxs_CoulRes(1,energy,MQ,scaleCouldown,muRdown,&
                     aSCouldown,aSMEdown,LOxs,DNLOxsNonCouldown,DResxsdown)
             ELSE
                WRITE(*,*)"ERROR: cannot reach here #5"
                STOP
             ENDIF
             NLO_wgtxsecmu(2)=(NLOxsup+DResxsup)/wme(1)
             NLO_wgtxsecmu(3)=(NLOxsdown+DResxsdown)/wme(1)
          ENDIF
          IF(wme(2).NE.0d0)THEN
             IF(coulomb.EQ.1)THEN
                ! LP
                CALL Get_LPxs_CoulRes(2,energy,MQ,aSCoulup,aSMEup,LOxs,DResxsup)
                CALL Get_LPxs_CoulRes(2,energy,MQ,aSCouldown,aSMEdown,LOxs,DResxsdown)
             ELSEIF(coulomb.EQ.2)THEN
                ! NLP
                CALL Get_NLPxs_CoulRes(2,energy,MQ,scaleCoulup,muRup,&
                     aSCoulup,aSMEup,LOxs,DNLOxsNonCoulup,DResxsup)
		CALL Get_NLPxs_CoulRes(2,energy,MQ,scaleCouldown,muRdown,&
                     aSCouldown,aSMEdown,LOxs,DNLOxsNonCouldown,DResxsdown)
             ELSE
                WRITE(*,*)"ERROR: cannot reach here #6"
                STOP
             ENDIF
             NNLO_wgtxsecmu(2)=(NNLOxsup+DResxsup)/wme(2)
             NNLO_wgtxsecmu(3)=(NNLOxsdown+DResxsdown)/wme(2)
          ENDIF
       ENDIF
    ENDIF
    
    RETURN
  END SUBROUTINE GetME2

  SUBROUTINE strf_pdf(wsf)
    USE global
    use Photon_PDFs
    USE ElasticPhotonPhotonFlux
    IMPLICIT NONE
    include '../vendor/gammaUPC/run90.inc'
    REAL(KIND(1d0)),INTENT(OUT)::wsf
    INTEGER::init=0
    SAVE init
    wsf=1d0
    IF(init.eq.0)then
       nuclearA_beam1=nuclearA1
       nuclearZ_beam1=nuclearZ1
       nuclearA_beam2=nuclearA2
       nuclearZ_beam2=nuclearZ2
       ebeam(1)=ebeam0(1)
       ebeam(2)=ebeam0(2)
       alphaem_elasticphoton=1d0/alphaemm1
       ! whether or not to use MC-Glauber TAA for the survival probability
       ! when colpar=1 and UPC_photon_flux=1 or 2.
       CALL ReadElem_logic('use_MC_Glauber',use_MC_Glauber)
       init=1
    endif
    if(colpar.eq.4)then
     ! gamma-gamma
       wsf=1d0
       WRITE(*,*)"ERROR: cannot reach here !"
    elseif(colpar.eq.3)then
       ! e-e+
       wsf=epa_electron(xp1,q2max)
       wsf=wsf*epa_electron(xp2,q2max)
    elseif(colpar.eq.2)then
       ! e-p
       wsf=epa_proton(xp1,q2max)
       wsf=wsf*epa_electron(xp2,q2max)
    else
       ! hadron-hadron UPC
       if(UPC_photon_flux.EQ.3)then
          ! iww (no survival probability) for pp
          wsf=epa_proton(xp1,q2max)
          wsf=wsf*epa_proton(xp2,q2max)
       else
          if(UPC_photon_flux.eq.1)then
             ! ChFF (see gamma-UPC 2207.03012)
             USE_CHARGEFORMFACTOR4PHOTON=.TRUE.
          else
             ! EDFF (see gamma-UPC 2207.03012)
             USE_CHARGEFORMFACTOR4PHOTON=.FALSE.
          endif
          IF(nuclearA_beam1.EQ.0.AND.nuclearA_beam2.EQ.0)THEN
             ! pp
             wsf=PhotonPhotonFlux_pp(xp1,xp2)
          ELSEIF((nuclearA_beam1.NE.0.AND.nuclearA_beam2.EQ.0).OR.&
               (nuclearA_beam1.EQ.0.AND.nuclearA_beam2.NE.0))THEN
             ! pA
             wsf=PhotonPhotonFlux_pA_WoodsSaxon(xp1,xp2)
          ELSE
             ! AB
             wsf=PhotonPhotonFlux_AB_WoodsSaxon(xp1,xp2)
          ENDIF
       endif
    endif
    
    IF(init.EQ.0)init=1
    RETURN
  END SUBROUTINE strf_pdf
  
  SUBROUTINE GET_ntau(tau,ntau0)
    ! get ntau for a given tau
    use global
    IMPLICIT NONE
    include '../vendor/gammaUPC/run90.inc'
    REAL(KIND(1d0)),INTENT(IN)::tau
    REAL(KIND(1d0)),INTENT(OUT)::ntau0
    REAL(KIND(1d0))::x
    x=DLOG10(tau)
    if(x.GT.0d0)then
       write(*,*)"ERROR: tau > 1"
       stop
    endif
    IF(colpar.eq.1)then
       if(nuclearZ_beam1.EQ.82.and.nuclearZ_beam2.EQ.82)then
        ! Pb+Pb
          if(x.LE.-7.3d0)then
             ntau0=-1.67d0
          else
             ntau0=-94.3051532633065d0-73.44040922736761d0*x-24.20040079794743d0*x**2-&
                  4.078561946253477d0*x**3-0.3478111311655984d0*x**4-0.011926801483889686d0*x**5
          endif
          ntau0=-ntau0
       elseif(nuclearZ_beam1.GE.54.and.nuclearZ_beam2.GE.54)then
          ! Xe+Xe or other two heavier ions
          if(x.LE.-7.3d0)then
             ntau0=-1.65d0
          else
             ntau0=-80.40283638414304d0-61.87697086384057d0*x-20.190330720309685d0*x**2-&
                  3.3649524772624226d0*x**3-0.2833770322998644d0*x**4-0.009585167381630894d0*x**5
          endif
          ntau0=-ntau0
       elseif(nuclearZ_beam1.EQ.20.and.nuclearZ_beam2.EQ.20)then
          ! Ca+Ca
          if(x.LE.-7.6d0)then
             ntau0=-1.5586d0
          else
             ntau0=-54.557802450791044d0-41.00008257855877d0*x-13.247561234064998d0*x**2-&
                  2.1906729816860953d0*x**3-0.18310288951474804d0*x**4-0.006144522758036254d0*x**5
          endif
          ntau0=-ntau0
       elseif(nuclearZ_beam1.EQ.18.and.nuclearZ_beam2.EQ.18)then
          ! Ar+Ar
          if(x.LE.-7.6d0)then
             ntau0=-1.54506d0
          else
             ntau0=-56.73289386251643d0-43.29067661237275d0*x-14.192594861897161d0*x**2-&
                  2.3816880463718246d0*x**3-0.20203012870909018d0*x**4-0.0068804496003979365d0*x**5
          endif
          ntau0=-ntau0
       elseif(nuclearZ_beam1.GE.8.and.nuclearZ_beam2.GE.8)then
          ! O+O or other two heavier ions
          if(x.LE.-7.6d0)then
             ntau0=-1.53397d0
          else
             ntau0=-39.30707034802285d0-28.904825168893716d0*x-9.308769766515642d0*x**2-&
                  1.5405970204691584d0*x**3-0.12911565048937673d0*x**4-0.004347965155734601d0*x**5
          endif
          ntau0=-ntau0
       elseif((nuclearZ_beam1.EQ.0.and.nuclearZ_beam2.GE.54)&
            .or.(nuclearZ_beam2.EQ.0.and.nuclearZ_beam1.GE.54))then
          ! p+Pb or Pb+p or a heavy ion+p
          if(x.LE.-7.8d0)then
             ntau0=-1.49d0
          else
             ntau0=-20.85512449173276d0-13.434567465290675d0*x-4.010519989747026d0*x**2-&
                  0.6244230265646026d0*x**3-0.04965283409587662d0*x**4-0.0015947318297498464d0*x**5
          endif
          ntau0=-ntau0
       elseif(nuclearA_beam1.EQ.0.and.nuclearA_beam2.EQ.0.and.&
            nuclearZ_beam1.EQ.0.and.nuclearZ_beam2.EQ.0)then
          ! pp
          if(x.LE.-8.2d0)then
             ntau0=-1.39d0
          else
             ntau0=-8.050931677173754d0-3.963308342683446d0*x-1.0494503104856845d0*x**2-&
                  0.1481666911671101d0*x**3-0.010796574134484576d0*x**4-0.0003194433638012126d0*x**5
          endif
          ntau0=-ntau0
       else
          WRITE(*,*)"WARNING: do not implement ntau for (Z1,Z2)=",nuclearZ_beam1,nuclearZ_beam2
          WRITE(*,*)"WARNING: will simply use ntau=2"
          ntau0=2d0
       endif
       return
    elseif(colpar.eq.2)then
       ! e+p collisions
       ntau0=2d0
       if(x.LE.-5.6d0)then
          ntau0=-1.45662d0
       else
          ntau0=-32.264482248394486d0-102.18795453870347d0*x-174.36487147984568d0*x**2-&
               180.61250551307256d0*x**3-121.34418842063245d0*x**4-54.551887771742415d0*x**5-&
               16.542376180975708d0*x**6-3.336791516816648d0*x**7-0.4285875710923549d0*x**8-&
               0.03168449402584781d0*x**9-0.001025143571009115d0*x**10
       endif
       ntau0=-ntau0
       ntau0=MAX(ntau0,1.01d0)
    elseif(colpar.eq.3)then
       ! e+e-
       if(x.LE.-3.9d0)then
          ntau0=-1.46482d0
       else
          ntau0=-24930.214875089572d0-28641.899212638815d0*x-5274.155833095974d0*x**2-&
               2283.3990257961455d0*x**3-998.9639559499716d0*x**4-351.0699145262098d0*x**5-&
               89.33462851169318d0*x**6-15.136722123477405d0*x**7-1.5133829561232466d0*x**8-&
               0.06718609263939923d0*x**9-22260.339669715704d0*DLog(-x)-&
               9161.961548713729d0*DLog(-x)**2-2202.9674489393824d0*DLog(-x)**3-&
               321.82005540583435d0*DLog(-x)**4-26.551577150222602d0*DLog(-x)**5-&
               0.9468751293556315d0*DLog(-x)**6
       endif
       ntau0=-ntau0
    elseif(colpar.eq.4)then
       ntau0=2d0
    else
       WRITE(*,*)"WARNING: do not implement ntau for colpar=",colpar
       WRITE(*,*)"WARNING: will simply use ntau=2"
       ntau0=2d0
    endif
    return
  END SUBROUTINE GET_ntau

  SUBROUTINE UPDATE_NLO_SCALE_ALPHAS(LOxs,NLOxs,Q,asatQ,muR)
    USE global_constants
    USE qcd_constants
    USE evaluate_couplings
    IMPLICIT NONE
    ! the input LOxs and NLOxs are evaluated at scale Q
    ! with alpha_s=asatQ
    ! Now, we want to calculate LOxs and NLOxs at muR with our alpha_s run
    REAL(KIND(1d0)),INTENT(INOUT)::LOxs,NLOxs
    REAL(KIND(1d0)),INTENT(IN)::Q,asatQ,muR
    REAL(KIND(1d0))::LRQ,asotwopi
    INTEGER::b
    REAL(KIND(1d0))::xs0hatQ,xs1hatQ
    REAL(KIND(1d0))::xs0hatmuR,xs1hatmuR
    LRQ=DLOG(muR**2/Q**2)
    b=0 ! the power of alpha_s at LO
    asotwopi=asatQ/TWOPI
    xs0hatQ=LOxs/(asotwopi)**b
    xs1hatQ=(NLOxs-LOxs)/(asotwopi)**(b+1)
    asotwopi=ALPHAS(muR)
    asotwopi=asotwopi/TWOPI
    xs0hatmuR=xs0hatQ
    xs1hatmuR=xs1hatQ+b*beta0*xs0hatQ*LRQ
    LOxs=xs0hatmuR*asotwopi**b
    NLOxs=xs1hatmuR*asotwopi**(b+1)
    NLOxs=LOxs+NLOxs
    RETURN
  END SUBROUTINE UPDATE_NLO_SCALE_ALPHAS

  SUBROUTINE UPDATE_NNLO_SCALE_ALPHAS(LOxs,NLOxs,NNLOxs,Q,asatQ,muR)
    USE global_constants
    USE qcd_constants
    USE evaluate_couplings
    IMPLICIT NONE
    ! the input LOxs, NLOxs, NNLOxs are evaluated at scale Q
    ! with alpha_s=asatQ
    ! Now, we want to calculate LOxs, NLOxs, NNLOxs at muR with our alpha_s run
    REAL(KIND(1d0)),INTENT(INOUT)::LOxs,NLOxs,NNLOxs
    REAL(KIND(1d0)),INTENT(IN)::Q,asatQ,muR
    REAL(KIND(1d0))::LRQ,asotwopi
    INTEGER::b
    REAL(KIND(1d0))::xs0hatQ,xs1hatQ,xs2hatQ
    REAL(KIND(1d0))::xs0hatmuR,xs1hatmuR,xs2hatmuR
    LRQ=DLOG(muR**2/Q**2)
    b=0 ! the power of alpha_s at LO
    asotwopi=asatQ/TWOPI
    xs0hatQ=LOxs/(asotwopi)**b
    xs1hatQ=(NLOxs-LOxs)/(asotwopi)**(b+1)
    xs2hatQ=(NNLOxs-NLOxs)/(asotwopi)**(b+2)
    asotwopi=ALPHAS(muR)
    asotwopi=asotwopi/TWOPI
    xs0hatmuR=xs0hatQ
    xs1hatmuR=xs1hatQ+b*beta0*xs0hatQ*LRQ
    xs2hatmuR=xs2hatQ+((b+1d0)*beta0*xs1hatQ+b*beta1*xs0hatQ)*LRQ&
         +0.5d0*b*(b+1d0)*beta0**2*xs0hatQ*LRQ**2
    LOxs=xs0hatmuR*asotwopi**b
    NLOxs=xs1hatmuR*asotwopi**(b+1)
    NLOxs=LOxs+NLOxs
    NNLOxs=xs2hatmuR*asotwopi**(b+2)
    NNLOxs=NLOxs+NNLOxs
    RETURN
  END SUBROUTINE UPDATE_NNLO_SCALE_ALPHAS

  SUBROUTINE Get_LOxs(energy,MQ,eQ,LOxs)
    IMPLICIT NONE
    ! eQ is the electric charge of heavy quark Q in unit of the positron charge
    ! energy and MQ are in unit of GeV
    REAL(KIND(1d0)),INTENT(IN)::energy,MQ,eQ
    REAL(KIND(1d0)),INTENT(OUT)::LOxs
    REAL(KIND(1d0))::betaQ
    REAL(KIND(1d0)),PARAMETER::GeVm2Topb=389379660d0 ! GeV-2 to pb
    REAL(KIND(1d0)),PARAMETER::pipi=3.14159265358979323846264338328d0
!    REAL(KIND(1d0)),PARAMETER::aem=0.0072992700729927005d0 ! 1/137.0
    REAL(KIND(1d0))::aem
    IF(energy.LE.2d0*MQ)THEN
       LOxs=0d0
       RETURN
    ENDIF
    betaQ=DSQRT(MAX(0d0,1d0-4d0*MQ**2/energy**2))
    IF(betaQ.EQ.0d0.OR.betaQ.GE.1d0)THEN
       LOxs=0d0
       RETURN
    ENDIF
    aem=1d0/alphaemm1
    LOxs=3d0*pipi*eQ**4*aem**2/MQ**2*(1d0-betaQ**2)
    LOxs=LOxs*((3d0-betaQ**4)/2d0*DLOG((1d0+betaQ)/(1d0-betaQ))&
         -betaQ*(2d0-betaQ**2))
    LOxs=LOxs*GeVm2Topb
    RETURN
  END SUBROUTINE Get_LOxs

  ! evaluate the Coulomb approximation of NLO QCD correction term
  SUBROUTINE Get_DNLOxs_Coul(energy,MQ,aS,LOxs,DNLOxs)
    IMPLICIT NONE
    ! energy and MQ are in unit of GeV
    REAL(KIND(1d0)),INTENT(IN)::energy,MQ
    REAL(KIND(1d0)),INTENT(IN)::aS,LOxs
    REAL(KIND(1d0)),INTENT(OUT)::DNLOxs
    REAL(KIND(1d0))::CF
    REAL(KIND(1d0))::betaQ
    REAL(KIND(1d0)),PARAMETER::pipi=3.14159265358979323846264338328d0
    IF(energy.LE.2d0*MQ.OR.aS.EQ.0d0.OR.LOxs.EQ.0d0)THEN
       DNLOxs=0d0
       RETURN
    ENDIF
    betaQ=DSQRT(MAX(0d0,1d0-4d0*MQ**2/energy**2))
    IF(betaQ.EQ.0d0.OR.betaQ.GE.1d0)THEN
       DNLOxs=0d0
       RETURN
    ENDIF
    CF=4d0/3d0
    DNLOxs=CF*pipi*aS*(1d0-betaQ**2)/(2d0*betaQ)*LOxs
    RETURN
  END SUBROUTINE Get_DNLOxs_Coul

  ! evaluate the Coulomb approximation of NNLO QCD correction term
  SUBROUTINE Get_DNNLOxs_Coul(nq,energy,MQ,muR,aS,LOxs,DNLOnonCoulxs,DNNLOxs)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nq
    ! energy and MQ are in unit of GeV
    REAL(KIND(1d0)),INTENT(IN)::energy,MQ,muR
    REAL(KIND(1d0)),INTENT(IN)::aS,LOxs,DNLOnonCoulxs
    REAL(KIND(1d0)),INTENT(OUT)::DNNLOxs
    REAL(KIND(1d0))::CF,CA,beta00,a11
    REAL(KIND(1d0))::betaQ,sqrtobetaQ,alogv
    REAL(KIND(1d0)),PARAMETER::pipi=3.14159265358979323846264338328d0
    IF(energy.LE.2d0*MQ.OR.aS.EQ.0d0)THEN
       DNNLOxs=0d0
       RETURN
    ENDIF
    betaQ=DSQRT(MAX(0d0,1d0-4d0*MQ**2/energy**2))
    IF(betaQ.EQ.0d0.OR.betaQ.GE.1d0)THEN
       DNNLOxs=0d0
       RETURN
    ENDIF
    CF=4d0/3d0
    CA=3d0
    beta00=11d0/6d0*CA-1d0/3d0*nq
    a11=31d0/9d0*CA-10d0/9d0*nq
    sqrtobetaQ=DSQRT(1d0-betaQ**2)
    alogv=DLOG(8d0*MQ**2/muR**2*(1d0-sqrtobetaQ)/(1d0-betaQ**2))
    DNNLOxs=aS**2*(CF**2*pipi**2*sqrtobetaQ**3/(12d0*betaQ**2)&
         +CF*(a11-2d0*beta00*alogv)/(8d0*betaQ))*LOxs&
         +CF*pipi*aS*(1d0-betaQ**2)/(2d0*betaQ)*DNLOnonCoulxs
    RETURN
  END SUBROUTINE Get_DNNLOxs_Coul

  ! evalute the high-energy (Regge) limit of NLO correction term
  ! it is a fit instead of a prediction
  SUBROUTINE Get_DNLOxs_HE(energy,MQ,aS,LOxs,DNLOxs)
    IMPLICIT NONE
    ! energy and MQ are in unit of GeV
    REAL(KIND(1d0)),INTENT(IN)::energy,MQ
    REAL(KIND(1d0)),INTENT(IN)::aS,LOxs
    REAL(KIND(1d0)),INTENT(OUT)::DNLOxs
    REAL(KIND(1d0)),PARAMETER::logtwo=0.693147180559945309417232121458d0
    REAL(KIND(1d0)),PARAMETER::pipi=3.14159265358979323846264338328d0
    REAL(KIND(1d0))::logsomQ2
    IF(energy.LE.2d0*MQ.OR.aS.EQ.0d0.OR.LOxs.EQ.0d0)THEN
       DNLOxs=0d0
       RETURN
    ENDIF
    logsomQ2=DLOG(energy**2/MQ**2)
    DNLOxs=LOxs*aS/pipi**2*logtwo*(logsomQ2**2-0.5d0*logsomQ2)
    RETURN
  END SUBROUTINE Get_DNLOxs_HE

  ! evalute the high-energy (Regge) limit of NNLO correction nq term (taking nq=1)
  ! it is a fit instead of a prediction
  SUBROUTINE Get_DNNLOnqxs_HE(energy,MQ,muR,aS,LOxs,DNLOxs,DNNLOxs)
    IMPLICIT NONE
    ! energy and MQ are in unit of GeV
    REAL(KIND(1d0)),INTENT(IN)::energy,MQ,muR
    REAL(KIND(1d0)),INTENT(IN)::aS,LOxs,DNLOxs
    REAL(KIND(1d0)),INTENT(OUT)::DNNLOxs
    REAL(KIND(1d0)),PARAMETER::logtwo=0.693147180559945309417232121458d0
    REAL(KIND(1d0)),PARAMETER::pipi=3.14159265358979323846264338328d0
    REAL(KIND(1d0))::logsomQ2
    IF(energy.LE.2d0*MQ.OR.aS.EQ.0d0)THEN
       DNNLOxs=0d0
       RETURN
    ENDIF
    logsomQ2=DLOG(energy**2/MQ**2)
    DNNLOxs=LOxs*aS**2/(2d0*pipi)**2*(logtwo**2*(8d0/3d0*logsomQ2**2-39d0/2d0*logsomQ2)+14d0)&
         -1d0/3d0*DNLOxs*aS/(2d0*pipi)*DLOG(muR**2/MQ**2)
    RETURN
  END SUBROUTINE Get_DNNLOnqxs_HE

  ! evalute the high-energy (Regge) limit of NNLO correction non-nq term
  ! it is a fit instead of a prediction
  SUBROUTINE Get_DNNLOnonnqxs_HE(energy,MQ,muR,aS,LOxs,DNLOxs,DNNLOxs)
    IMPLICIT NONE
    ! energy and MQ are in unit of GeV
    REAL(KIND(1d0)),INTENT(IN)::energy,MQ,muR
    REAL(KIND(1d0)),INTENT(IN)::aS,LOxs,DNLOxs
    REAL(KIND(1d0)),INTENT(OUT)::DNNLOxs
    REAL(KIND(1d0)),PARAMETER::logtwo=0.693147180559945309417232121458d0
    REAL(KIND(1d0)),PARAMETER::pipi=3.14159265358979323846264338328d0
    REAL(KIND(1d0))::CA,logsomQ2
    IF(energy.LE.2d0*MQ.OR.aS.EQ.0d0)THEN
       DNNLOxs=0d0
       RETURN
    ENDIF
    CA=3d0
    logsomQ2=DLOG(energy**2/MQ**2)
    DNNLOxs=LOxs*aS**2/(2d0*pipi)**2*(logtwo**2*(1d0/3d0*logsomQ2**4&
         -8d0/3d0*logsomQ2**3-167d0/4d0*logsomQ2**2+775d0/2d0*logsomQ2)-328d0)&
         +11d0/6d0*CA*DNLOxs*aS/(2d0*pipi)*DLOG(muR**2/MQ**2)
    RETURN
  END SUBROUTINE Get_DNNLOnonnqxs_HE

  SUBROUTINE NNLO_scale(scale0)
    USE global
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(OUT)::scale0
    REAL(KIND(1d0))::px,py,ptw2,mas
    REAL(KIND(1d0)),DIMENSION(0:3)::pmom
    INTEGER::ij
    REAL(KIND(1d0)),PARAMETER::scale_min=1d0 ! cannot below 1 GeV
    SELECT CASE(scale)
    CASE(0)
       scale0=fscalevalue
    CASE(1)
       ! invariant mass of the photon
       scale0=Wgaga
    CASE DEFAULT
       WRITE(*,*)"ERROR: Do not understand the scale scheme=",scale
       STOP
    END SELECT
    scale0=MAX(scale0*muR_over_ref,scale_min)
    RETURN
  END SUBROUTINE NNLO_scale

  ! evaluate the Coulomb resummation at LP
  ! with width = 0 and energy >= 2*MQ
  ! according to order, we have subtracted the double counting part wrt fixed order
  SUBROUTINE Get_LPxs_CoulRes(iorder,energy,MQ,aSCoul,aSME,LOxs,DResxs)
    IMPLICIT NONE
    INTEGER::iorder
    ! energy and MQ are in unit of GeV
    REAL(KIND(1d0)),INTENT(IN)::energy,MQ
    ! aSCoul is aS value in the Coulomb resummation
    ! while aSME is aS value in the ME
    REAL(KIND(1d0)),INTENT(IN)::aSCoul,aSME
    REAL(KIND(1d0)),INTENT(IN)::LOxs
    REAL(KIND(1d0)),INTENT(OUT)::DResxs
    REAL(KIND(1d0))::CF
    REAL(KIND(1d0))::betaQ,sqrtonembetaQ2
    REAL(KIND(1d0)),PARAMETER::pipi=3.14159265358979323846264338328d0
    IF(energy.LE.2d0*MQ.OR.LOxs.EQ.0d0)THEN
       DResxs=0d0
       RETURN
    ENDIF
    betaQ=DSQRT(MAX(0d0,1d0-4d0*MQ**2/energy**2))
    IF(betaQ.EQ.0d0.OR.betaQ.GE.1d0)THEN
       DResxs=0d0
       RETURN
    ENDIF
    CF=4d0/3d0
    sqrtonembetaQ2=DSQRT(1d0-betaQ**2)
    ! Sommerfeld enhancement factor
    DResxs=CF*aSCoul*pipi*sqrtonembetaQ2**2/betaQ
    DResxs=DResxs/(1d0-DEXP(-CF*aSCoul*pipi*sqrtonembetaQ2/betaQ))
    DResxs=(DResxs-sqrtonembetaQ2)*LOxs
    IF(iorder.EQ.0)RETURN
    DResxs=DResxs-CF*aSME*pipi/2d0*sqrtonembetaQ2**2/betaQ*LOxs
    IF(iorder.EQ.1)RETURN
    DResxs=DResxs-CF**2*aSME**2*pipi**2*sqrtonembetaQ2**3/(12d0*betaQ**2)*LOxs
    IF(iorder.EQ.2)RETURN
    WRITE(*,*)"ERROR: cannot reach here"
    STOP
    RETURN
  END SUBROUTINE Get_LPxs_CoulRes

  ! evaluate the Coulomb resummation at NLP
  ! with width = 0 and energy >= 2*MQ
  ! according to order, we have subtracted the double counting part wrt fixed order
  SUBROUTINE Get_NLPxs_CoulRes(iorder,energy,MQ,muC,muR,&
       aSCoul,aSME,LOxs,DNLOxsNonCoul,DResxs)
    USE potentialfunction
    IMPLICIT NONE
    INTEGER::iorder
    ! energy, MQ, and muC are in unit of GeV
    REAL(KIND(1d0)),INTENT(IN)::energy,MQ,muC,muR
    ! aSCoul is aS value in the Coulomb resummation
    ! while aSME is aS value in the ME
    REAL(KIND(1d0)),INTENT(IN)::aSCoul,aSME
    REAL(KIND(1d0)),INTENT(IN)::LOxs,DNLOxsNonCoul
    REAL(KIND(1d0)),INTENT(OUT)::DResxs
    REAL(KIND(1d0))::CF,D1,J1NLP,J1NLPaS2,Eb
    REAL(KIND(1d0))::betaQ,sqrtonembetaQ2
    REAL(KIND(1d0)),PARAMETER::pipi=3.14159265358979323846264338328d0
    REAL(KIND(1d0))::J1LP,J1LPaS0,J1LPaS1,J1LPaS2,J1LPexp0,J1LPexp1,J1PSfact
    IF(energy.LE.2d0*MQ.OR.LOxs.EQ.0d0)THEN
       DResxs=0d0
       RETURN
    ENDIF
    betaQ=DSQRT(MAX(0d0,1d0-4d0*MQ**2/energy**2))
    IF(betaQ.EQ.0d0.OR.betaQ.GE.1d0)THEN
       DResxs=0d0
       RETURN
    ENDIF
    CF=4d0/3d0
    D1=-CF
    ! we make sure the only the potential function uses muC
    ! while the hard function and uses muR
    ! NLP potential function
    Eb=energy-2d0*MQ
    CALL Get_potentialJ_NLP(nq,aSCoul,muC,D1,Eb,MQ,0d0,J1NLP)
    J1NLPaS2=0d0
    IF(iorder.GE.2)THEN
       CALL Get_potentialJ_NLP_aS2(nq,aSME,muR,D1,Eb,MQ,0d0,J1NLPaS2)
    ENDIF
    sqrtonembetaQ2=DSQRT(1d0-betaQ**2)

    ! LP is simply Sommerfeld enhancement factor
    J1LP=CF*aSCoul*MQ**2/2d0/(1d0-DEXP(-CF*aSCoul*pipi*sqrtonembetaQ2/betaQ))
    ! this is LP J1 at O(aS^0)
    J1LPaS0=MQ**2*betaQ/(2d0*pipi*sqrtonembetaQ2)
    J1LPexp0=J1LPaS0
    J1LPexp1=0d0
    IF(iorder.GE.1)THEN
       ! this is LP J1 at O(aS^1)
       J1LPaS1=0.25d0*CF*MQ**2*aSME
       J1LPexp0=J1LPexp0+J1LPaS1
       J1LPexp1=J1LPexp1+J1LPaS0
    ENDIF
    IF(iorder.GE.2)THEN
       ! this is LP J1 at O(aS^2)
       J1LPaS2=aSME**2*CF**2*MQ**2*pipi*sqrtonembetaQ2/(24d0*betaQ)
       J1LPexp0=J1LPexp0+J1LPaS2
       J1LPexp1=J1LPexp1+J1LPaS1
    ENDIF
    ! this is the pure phase space factor in eq.(2.7) of arXiv:2506.23791
    ! it is amount to J1LPaS0 dividing by sqrtonembetaQ2
    J1PSfact=J1LPaS0/sqrtonembetaQ2
    
    DResxs=(J1LP-J1LPexp0+J1NLP-J1NLPaS2)/J1PSfact*LOxs
    DResxs=DResxs+(J1LP-J1LPexp1)/J1PSfact*DNLOxsNonCoul
    
    RETURN
  END SUBROUTINE Get_NLPxs_CoulRes

  SUBROUTINE Get_CoulRes_Scale(energy,MQ,muR,scale)
    USE newton_method
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::energy,MQ,muR
    REAL(KIND(1d0)),INTENT(OUT)::scale
    REAL(KIND(1d0))::betaQ
    INTEGER,SAVE::init=0
    REAL(KIND(1d0))::muC,FPMUC
    REAL(KIND(1d0)),DIMENSION(2)::FPMUC2
    REAL(KIND(1d0)),DIMENSION(4:6),PARAMETER::scale_min=(/1d0,1.8d0,32d0/)
    REAL(KIND(1d0)),SAVE::muC_min
    IF(init.EQ.0)THEN
       muC=scale_min(quark)
       CALL newtonsolver(muC,FPMUC,SCALE_IN_MUC_EQ,1D-14,-1)
       IF(FPMUC.EQ.-1D99.OR.ISNAN(FPMUC).OR.1D0/FPMUC.EQ.0d0)THEN
          PRINT *, "WARNING: Linear Newton procedure does not work !"
         !PRINT *, "INFO: STARTING TO USE Quadratic Newton procedure"
         ! Note: it is necessary to use quadratic Netown's method than the linear Newton's method
         ! Especially for Q close to 0.5 GeV beyond 1-loop level calculation. Otherwise, the linear method
         ! has too bad convergence !!!
         muC=scale_min(quark)
         !CALL newtonsolver(muC,FPMUC2,SCALE_IN_MUC_EQ_QUA,1D-14,30)
      ENDIF
      IF(muC.LT.0.5d0.OR.ISNAN(muC))THEN
         muC=scale_min(quark)
      ENDIF
      !muC_min=MAX(muC,1d0) ! forbid to below 1 GeV
      muC_min=muC
      init=1
    ENDIF
    IF(energy.LE.2d0*MQ)THEN
       scale=MIN(muR,muC_min)
       RETURN
    ENDIF
    betaQ=DSQRT(MAX(0d0,1d0-4d0*MQ**2/energy**2))
    IF(betaQ.EQ.0d0)THEN
       scale=MIN(muR,muC_min)
       RETURN
    ENDIF
    IF(betaQ.GE.1d0)THEN
       scale=muR
       RETURN
    ENDIF
    scale=MAX(muC_min,muR*2d0*betaQ)
    scale=MIN(muR,scale)
    RETURN
  END SUBROUTINE Get_CoulRes_Scale

  SUBROUTINE SCALE_IN_MUC_EQ(muC,f,fp)
    ! f(x)=CF*mQ*as(muC)-muC
    USE qcd_constants
    USE evaluate_couplings
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::muC
    REAL(KIND(1d0)),INTENT(OUT)::f,fp
    REAL(KIND(1d0)),PARAMETER::twopi=6.28318530717958647692528676656d0
    REAL(KIND(1d0))::asmuC,asmuCo2pi,betaasmuC
    IF(muC.LT.0.5d0.OR.ISNAN(muC))THEN
       PRINT *, "WARNING: THE CONVERGENCE OF LINEAR NEWTWON METHOD IS TOO BAD !"
       PRINT *, "muC=",muC
       f=0d0
       fp=-1d99
       RETURN
    ENDIF
    asmuC=ALPHAS(muC)
    asmuCo2pi=asmuC/twopi
    f=CF*MQ*asmuC-muC
    betaasmuC=beta0*asmuCo2pi
    IF(alphas_nloop.GE.2)THEN
       betaasmuC=betaasmuC+beta1*asmuCo2pi**2
    ENDIF
    IF(alphas_nloop.GE.3)THEN
       betaasmuC=betaasmuC+beta2*asmuCo2pi**3
    ENDIF
    IF(alphas_nloop.GE.4)THEN
       betaasmuC=betaasmuC+beta3*asmuCo2pi**4
    ENDIF
    IF(alphas_nloop.GE.5)THEN
       betaasmuC=betaasmuC+beta4*asmuCo2pi**5
    ENDIF
    fp=CF*MQ*2d0/muC*(-asmuC)*betaasmuC-1d0
    RETURN
  END SUBROUTINE SCALE_IN_MUC_EQ

END MODULE phase_space_integrator
