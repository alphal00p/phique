MODULE MC_VEGAS
!
!  This module is a modification f95 version of VEGA_ALPHA.for
!  by G.P. LEPAGE SEPT 1976/(REV)AUG 1979.
!  
!   author: Hua-Sheng Shao
!   Physics school of Peking University
!
IMPLICIT NONE
SAVE
INTEGER::vegas_grid_stage=0 ! 0: do nothing for the grid (do not use the outputed grid and do not write out the grid);
                            ! 1: calculate the grid and output it
                            ! 2: use the grid from the outputed one, i.e. grid_vegas
CHARACTER(len=24),PARAMETER::grid_vegas="./input/vegas.grid"
INTEGER,PARAMETER::grid_vegas_unit=38568
INTEGER,PARAMETER::MAX_SIZE=20            ! The max dimensions of the integrals
! The lower and upper integration limit
INTEGER,PRIVATE::i_vegas
REAL(KIND(1d0)),DIMENSION(MAX_SIZE),PUBLIC::XL=(/(0d0,i_vegas=1,MAX_SIZE)/),&
                                      XU=(/(1d0,i_vegas=1,MAX_SIZE)/)
INTEGER,PUBLIC::NCALL=5000,&              ! The number of integrand evaluations per iteration
                ITMX=5,&                  ! The maximum number of iterations
				NPRN=5,&                  ! printed or not
				NDEV=6,&                  ! device number for output
				IT=0,&                    ! number of iterations completed
				NDO=1,&                   ! number of subdivisions on an axis
				NDMX=50,&       ! determines the maximum number of increments along each axis
				MDS=1                     ! =0 use importance sampling only
				                          ! =\0 use importance sampling and stratified sampling
										  ! increments are concentrated either wehre the
										  ! integrand is largest in magnitude (MDS=1), or
										  ! where the contribution to the error is largest(MDS=-1)
INTEGER,PUBLIC::IINIP
REAL(KIND(1d0)),PUBLIC::ACC=-1d0           ! Algorithm stops when the relative accuracy,
                                          ! |SD/AVGI|, is less than ACC; accuracy is not 
										  ! cheched when ACC<0
REAL(KIND(1d0)),PUBLIC::MC_SI=0d0,&        ! sum(AVGI_i/SD_i^2,i=1,IT)
                SWGT=0d0,&                ! sum(1/SD_i^2,i=1,IT)
				SCHI=0d0,&                ! sum(AVGI_i^2/SD_i^2,i=1,IT)
				ALPH=1.5d0                ! controls the rate which the grid is modified from
				                          ! iteration to iteration; decreasing ALPH slows
										  ! modification of the grid
										  ! (ALPH=0 implies no modification)
REAL(KIND(1d0)),PUBLIC::DSEED=1234567d0    ! seed of 
! location of the I-th division on the J-th axis, normalized to lie between 0 and 1.
REAL(KIND(1d0)),DIMENSION(50,MAX_SIZE),PUBLIC::XI=1d0
REAL(KIND(1d0)),PUBLIC::CALLS,TI,TSI

CONTAINS

! Generate random numbers
SUBROUTINE RANDA(NR,R)
!SPECIFICATIONS FOR ARGUMENTS
INTEGER,INTENT(IN)::NR
REAL(KIND(1d0)),DIMENSION(NR),INTENT(OUT)::R
!SPECIFICATIONS FOR LOCAL VARIABLES
INTEGER::I
! D2P31M=(2**31) - 1
! D2P31 =(2**31)(OR AN ADJUSTED VALUE)
REAL(KIND(1d0))::D2P31M=2147483647.d0,D2P31=2147483711.d0
!FIRST EXECUTABLE STATEMENT
DO I=1,NR
   DSEED = DMOD(16807.d0*DSEED,D2P31M)
   R(I) = DSEED / D2P31
ENDDO
END SUBROUTINE RANDA

SUBROUTINE VEGAS(NDIM,FXN,AVGI,SD,CHI2A,INIT)
!
!     SUBROUTINE PERFORMS NDIM-DIMENSIONAL MONTE CARLO INTEG'N
!     - BY G.P. LEPAGE    SEPT 1976/(REV)AUG 1979
!     - ALGORITHM DESCRIBED IN J COMP PHYS 27,192(1978)
!
! Without INIT or INIT=0, CALL VEGAS
! INIT=1  CALL VEGAS1
! INIT=2  CALL VEGAS2
! INIT=3  CALL VEGAS3
INTEGER,INTENT(IN)::NDIM
REAL(KIND(1d0)),EXTERNAL::FXN
INTEGER,INTENT(IN),OPTIONAL::INIT
REAL(KIND(1d0)),INTENT(INOUT)::AVGI,SD,CHI2A
REAL(KIND(1d0)),DIMENSION(50,MAX_SIZE)::D,DI
REAL(KIND(1d0)),DIMENSION(50)::XIN,R
REAL(KIND(1d0)),DIMENSION(MAX_SIZE)::DX,X,DT,RAND
INTEGER,DIMENSION(MAX_SIZE)::IA,KG
INTEGER::initflag
REAL(KIND(1d0)),PARAMETER::ONE=1.d0
INTEGER::I,J,K,NPG,NG,ND,NDM,LABEL=0
REAL(KIND(1d0))::DXG,DV2G,XND,XJAC,RC,XN,DR,XO,TI2,WGT,FB,F2B,F,F2
LOGICAl::lexist
INTEGER::NDIM_READ,ND_READ
!SAVE AVGI,SD,CHI2A
!SQRT(A)=DSQRT(A)
!ALOG(A)=DLOG(A)
!ABS(A)=DABS(A)
IF(PRESENT(INIT))THEN
   initflag=INIT
ELSE
   initflag=0
ENDIF
! INIT=0  - INITIALIZES CUMULATIVE VARIABLES AND GRID
ini0:IF(initflag.LT.1) THEN
   NDO=1
   DO  J=1,NDIM
       XI(1,J)=ONE
   ENDDO
ENDIF ini0
IF(vegas_grid_stage.EQ.2)THEN
   INQUIRE(FILE=TRIM(grid_vegas),EXIST=lexist)
   IF(.NOT.lexist)THEN
      WRITE(*,*)"WARNING: the vegas.grid does not exist !"
      vegas_grid_stage=1
   ELSE
      OPEN(UNIT=grid_vegas_unit,FILE=TRIM(grid_vegas))
      READ(grid_vegas_unit,*)ND_READ,NDIM_READ
      IF(NDIM_READ.NE.NDIM)THEN
         WRITE(*,*)"WARNING: the dimension of vegas.grid does not match !"
         vegas_grid_stage=1
      ELSE
         DO J=1,NDIM_READ
            READ(grid_vegas_unit,*)(XI(I,J),I=1,ND_READ)
         ENDDO
         ND=ND_READ
         NDO=ND_READ
      ENDIF
      CLOSE(UNIT=grid_vegas_unit)
   ENDIF
ENDIF
!  INIT=1    - INITIALIZES CUMULATIVE VARIABLES, BUT NOT GRID     
ini1:IF(initflag.LT.2)THEN
  IT=0
  MC_SI=0.d0
  SWGT=MC_SI
  SCHI=MC_SI
ENDIF ini1
!  INIT=2   - NO INITIALIZATION
ini2:IF(initflag.LE.2)THEN
   ND=NDMX
   NG=1
   IF(MDS.NE.0) THEN
     NG=(NCALL/2.d0)**(1.d0/NDIM)
     MDS=1
     IF((2*NG-NDMX).GE.0) THEN
       MDS=-1
       NPG=NG/NDMX+1
       ND=NG/NPG
       NG=NPG*ND
    ENDIF
   ENDIF
   K=NG**NDIM                      ! K sub volumes
   NPG=NCALL/K                     ! The number of random numbers in per sub volumes Ms
   IF(NPG.LT.2)NPG=2
   CALLS=DBLE(NPG*K)               ! The total number of random numbers M
   DXG=ONE/NG
   DV2G=(CALLS*DXG**NDIM)**2/NPG/NPG/(NPG-ONE)  ! 1/(Ms-1)
   XND=ND                          ! ~NDMX! 
                                   ! determines the number of increments along each axis
   NDM=ND-1                        ! ~NDMX-1
   DXG=DXG*XND                     ! determines the number of increments along each axis per sub-v
   XJAC=ONE/CALLS
   DO J=1,NDIM
      DX(J)=XU(J)-XL(J)
      XJAC=XJAC*DX(J)              ! XJAC=Volume/M
   ENDDO
!     REBIN, PRESERVING BIN DENSITY
   IF(ND.NE.NDO) THEN
      RC=NDO/XND                   ! XND=ND
      outer:DO J=1, NDIM           ! Set the new division
          K=0
          XN=0.d0
          DR=XN
          I=K
          LABEL=0
          inner5:DO
             IF(LABEL.EQ.0) THEN
                inner4:DO
                   K=K+1
                   DR=DR+ONE
                   XO=XN
                   XN=XI(K,J)
                   IF(RC.LE.DR) EXIT
                ENDDO inner4
             ENDIF
             I=I+1
             DR=DR-RC
             XIN(I)=XN-(XN-XO)*DR
             IF(I.GE.NDM) THEN
                EXIT
             ELSEIF(RC.LE.DR) THEN
                LABEL=1
             ELSE
                LABEL=0
             ENDIF
          ENDDO inner5
          inner:DO I=1,NDM
             XI(I,J)=XIN(I)
          ENDDO inner
          XI(ND,J)=ONE
      ENDDO outer
      NDO=ND
   ENDIF

   IF(NPRN.GE.0) WRITE(NDEV,200) NDIM,CALLS,IT,ITMX,ACC,NPRN,&
                         ALPH,MDS,ND,(XL(J),XU(J),J=1,NDIM)
ENDIF ini2
!      ENTRY VEGAS3(NDIM,FXN,AVGI,SD,CHI2A)
!  INIT=3       - MAIN INTEGRATION LOOP
mainloop:DO
    IT=IT+1
    TI=0.d0
    TSI=TI
    DO J=1,NDIM
       KG(J)=1
       DO I=1,ND
          D(I,J)=TI
          DI(I,J)=TI
       ENDDO
    ENDDO
    
    LABEL=0
    level1:DO
       level2:DO
          ifla:IF(LABEL.EQ.0)THEN
             FB=0.d0
             F2B=FB
             level3:DO K=1,NPG
                CALL RANDA(NDIM,RAND)
                WGT=XJAC
                DO J=1,NDIM
                   XN=(KG(J)-RAND(J))*DXG+ONE
                   IA(J)=XN
                   IF(IA(J).LE.1) THEN
                      XO=XI(IA(J),J)
                      RC=(XN-IA(J))*XO
                   ELSE
                      XO=XI(IA(J),J)-XI(IA(J)-1,J)
                      RC=XI(IA(J)-1,J)+(XN-IA(J))*XO
                   ENDIF
                   X(J)=XL(J)+RC*DX(J)
                   WGT=WGT*XO*XND
                ENDDO
                
                F=WGT
                F=F*FXN(X,WGT)
                F2=F*F
                FB=FB+F
                F2B=F2B+F2
                DO J=1,NDIM
                   DI(IA(J),J)=DI(IA(J),J)+F
                   IF(MDS.GE.0) D(IA(J),J)=D(IA(J),J)+F2
                ENDDO
             ENDDO level3
             !          K=K-1                    !K=NPG
             
             F2B=DSQRT(F2B*DBLE(NPG))
             F2B=(F2B-FB)*(F2B+FB)
             TI=TI+FB
             TSI=TSI+F2B
             IF(MDS.LT.0) THEN
                DO J=1,NDIM
                   D(IA(J),J)=D(IA(J),J)+F2B
                ENDDO
             ENDIF
             K=NDIM
          ENDIF ifla
          KG(K)=MOD(KG(K),NG)+1
          IF(KG(K).EQ.1) THEN
             EXIT
          ELSE
             LABEL=0
          ENDIF
       ENDDO level2
       K=K-1
       IF(K.GT.0) THEN
          LABEL=1
       ELSE
          EXIT
       ENDIF
    ENDDO level1

!    COMPUTE FINAL RESULTS FOR THIS ITERATION
    TSI=TSI*DV2G
    TI2=TI*TI
    WGT=ONE/TSI
    MC_SI=MC_SI+TI*WGT
    SWGT=SWGT+WGT
    SCHI=SCHI+TI2*WGT
    AVGI=MC_SI/SWGT
    CHI2A=(SCHI-MC_SI*AVGI)/(IT-0.9999d0)
    SD=DSQRT(ONE/SWGT)
    IF(NPRN.GE.0) THEN
       TSI=DSQRT(TSI)
       WRITE(NDEV,201) IT,TI,TSI,AVGI,SD,CHI2A
    ENDIF
    IF(NPRN.GT.0) THEN
       DO J=1,NDIM
          WRITE(NDEV,202) J,(XI(I,J),DI(I,J),I=1+NPRN/2,ND,NPRN)
       ENDDO
    ENDIF

!   REFINE GRID
!   XI(k,j)=XI(k,j)-(XI(k,j)-XI(k-1,j))*(sum(R(i),i=1,k)-s*sum(R(i),i=1,ND)/M)/R(k)
!   divides the original k-th interval into s parts
    outer2:DO J=1,NDIM
       XO=D(1,J)
       XN=D(2,J)
       D(1,J)=(XO+XN)/2.d0
       DT(J)=D(1,J)
       inner2:DO I=2,NDM
          D(I,J)=XO+XN
          XO=XN
          XN=D(I+1,J)
          D(I,J)=(D(I,J)+XN)/3.d0
          DT(J)=DT(J)+D(I,J)
       ENDDO inner2
       D(ND,J)=(XN+XO)/2.d0
       DT(J)=DT(J)+D(ND,J)
    ENDDO outer2

    le1:DO J=1,NDIM
       RC=0.d0
       DO I=1,ND
          R(I)=0.d0
          IF(D(I,J).GT.0.) THEN
             XO=DT(J)/D(I,J)
             R(I)=((XO-ONE)/XO/DLOG(XO))**ALPH
          ENDIF
          RC=RC+R(I)
       ENDDO
       RC=RC/XND
       K=0
       XN=0.d0
       DR=XN
       I=K
       LABEL=0
       le2:DO
          le3:DO
             IF(LABEL.EQ.0)THEN
                K=K+1
                DR=DR+R(K)
                XO=XN
                XN=XI(K,J)
             ENDIF
             IF(RC.LE.DR) THEN
                EXIT
             ELSE
                LABEL=0
             ENDIF
          ENDDO le3
          I=I+1
          DR=DR-RC
          XIN(I)=XN-(XN-XO)*DR/R(K)
          IF(I.GE.NDM) THEN
             EXIT
          ELSE
             LABEL=1
          ENDIF
       ENDDO le2
       DO I=1,NDM
          XI(I,J)=XIN(I)
       ENDDO
       XI(ND,J)=ONE
    ENDDO le1

    IF(IT.GE.ITMX.OR.ACC*ABS(AVGI).GE.SD) EXIT
 ENDDO mainloop

! write out the grid if needed
IF(vegas_grid_stage.GE.1)THEN
   OPEN(UNIT=grid_vegas_unit,FILE=TRIM(grid_vegas))
   WRITE(grid_vegas_unit,*)ND,NDIM
   DO J=1,NDIM
      WRITE(grid_vegas_unit,*)(XI(I,J),I=1,ND)
   ENDDO
   CLOSE(UNIT=grid_vegas_unit)
ENDIF

200   FORMAT(/," INPUT PARAMETERS FOR MC_VEGAS: ",/," NDIM=",I3,"    NCALL=",F8.0,&
     "     IT=",I3,/," ITMX=",I3,"    ACC=   ",G9.3,&
     "   NPRN=",I3,/," ALPH=",F5.2,"    MDS=",I3,"          ND=",I4,/,&
     "(XL,XU)=",(T10,"(" G12.6,",",G12.6 ")"))
201   FORMAT(/," INTEGRATION BY MC_VEGAS ", " ITERATION NO. ",I3, /,&
     " INTEGRAL = ",G14.8, /," SQURE DEV  = ",G10.4,/,&
     " ACCUMULATED RESULTS:   INTEGRAL = ",G14.8,/,&
     " DEV  = ",G10.4, /," CHI**2 PER IT'N = ",G10.4)
! X is the division of the coordinate
! DELTA I is the sum of F in this interval 
202   FORMAT(/,"DATA FOR AXIS ",I2,/,"    X       DELTA I       ", &
     24H   X       DELTA I      ,18H   X       DELTA I, &
      /(1H ,F7.6,1X,G11.4,5X,F7.6,1X,G11.4,5X,F7.6,1X,G11.4))
END SUBROUTINE VEGAS
END MODULE MC_VEGAS
