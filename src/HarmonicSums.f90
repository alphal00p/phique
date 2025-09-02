MODULE HarmonicSums
  IMPLICIT NONE
CONTAINS
  ! single harmonic sums S_n(z)
  FUNCTION Harmonic_Sum_S(n,z)
    USE Euler_gamma
    IMPLICIT NONE
    INTEGER,INTENT(IN)::n
    COMPLEX(KIND(1d0)),INTENT(IN)::z
    COMPLEX(KIND(1d0))::Harmonic_Sum_S
    ! EulerGamma
    REAL(KIND(1d0)),PARAMETER::gammaE=0.577215664901532860606512090082d0
    COMPLEX(KIND(1d0)),EXTERNAL::polygamma
    COMPLEX(KIND(1d0)),PARAMETER::ipi=DCMPLX(0d0,3.14159265358979323846264338328d0)
    REAL(KIND(1d0)),PARAMETER::logtwo=0.693147180559945309417232121458d0
    INTEGER::an
    REAL(KIND(1d0))::rz
    rz=DREAL(z)
    IF(DIMAG(z).EQ.0d0.AND.rz.EQ.INT(rz).AND.rz.LT.0d0)THEN
       WRITE(*,*)"ERROR: Harmonic_Sum_S hits the pole !"
       WRITE(*,*)"(n,z)=",n,z
       STOP
    ENDIF
    IF(n.EQ.1)THEN
       Harmonic_Sum_S=polygamma(0,z+1d0)+gammaE
    ELSEIF(n.GT.1)THEN
       IF(MOD(n-1,2).EQ.0)THEN
          Harmonic_Sum_S=1
       ELSE
          Harmonic_Sum_S=-1
       ENDIF
       Harmonic_Sum_S=Harmonic_Sum_S/CGAMMA(DCMPLX(n,0d0))*polygamma(n-1,z+1d0)&
            +RIEMANNZETA(DBLE(n))
    ELSEIF(n.EQ.-1)THEN
       Harmonic_Sum_S=EXP(ipi*z)/2d0*(polygamma(0,z/2d0+1d0)-polygamma(0,z/2d0+0.5d0))-logtwo
    ELSEIF(n.LE.-2)THEN
       an=-n
       IF(MOD(an-1,2).EQ.0)THEN
          Harmonic_Sum_S=1
       ELSE
          Harmonic_Sum_S=-1
       ENDIF
       Harmonic_Sum_S=Harmonic_Sum_S/CGAMMA(DCMPLX(an,0d0))*EXP(ipi*z)&
            /2d0**an*(polygamma(an-1,z/2d0+1d0)-polygamma(an-1,z/2d0+0.5d0))&
            -(1d0-2d0**(1-an))*RIEMANNZETA(DBLE(an))
    ENDIF
    RETURN
  END FUNCTION Harmonic_Sum_S

  ! Nested Harmonic sum
  ! S_{2,1}(z)
  RECURSIVE FUNCTION Harmonic_Sum_S21(z) RESULT(res)
    IMPLICIT NONE
    COMPLEX(KIND(1d0)),INTENT(IN)::z
    COMPLEX(KIND(1d0))::res
    REAL(KIND(1d0))::rz
    INTEGER::ir,n1
    REAL(KIND(1d0)),DIMENSION(9)::B2k=(/1d0/6d0,-1d0/30d0,1d0/42d0,&
         -1d0/30d0,5d0/66d0,-691d0/2730d0,7d0/6d0,-3617d0/510d0,43867d0/798d0/)
    ! EulerGamma
    REAL(KIND(1d0)),PARAMETER::gammaE=0.577215664901532860606512090082d0
    ! Zeta[3]
    REAL(KIND(1d0)),PARAMETER::zeta3=1.20205690315959428539973816151d0
    COMPLEX(KIND(1d0))::NNm1,logNbar
    rz=DREAL(z)
    IF(DIMAG(z).EQ.0d0.AND.rz.EQ.INT(rz).AND.rz.LT.0d0)THEN
       WRITE(*,*)"ERROR: Harmonic_Sum_S21 hits the pole !"
       WRITE(*,*)"z=",z
       STOP
    ENDIF
    IF(rz.LT.5d0)THEN
       ir=5-FLOOR(rz)
       res=Harmonic_Sum_S21(z+ir)
       DO n1=1,ir
          res=res-1d0/(n1+z)**2*Harmonic_Sum_S(1,n1+z)
       ENDDO
       RETURN
    ENDIF
    logNbar=LOG(z)+gammaE
    NNm1=z**(-1)
    res=(1+logNbar)*NNm1
    res=res+(0.25d0-0.5d0*logNbar)*NNm1**2
    res=res+(-0.25d0-2d0/3d0*B2k(1)+logNbar*B2k(1))*NNm1**3
    res=res+B2k(1)*NNm1**4
    res=res+(-B2k(1)**2-17d0/15d0*B2k(2)+B2k(2)*logNbar)*NNm1**5
    res=res+11d0/8d0*B2k(2)*NNm1**6
    res=res+(-13d0/4d0*B2k(1)*B2k(2)-619d0/420d0*B2k(3)+logNbar*B2k(3))*NNm1**7
    res=res+11d0/6d0*B2k(3)*NNm1**8
    res=res+(-3.5d0*B2k(2)**2-16d0/3d0*B2k(1)*B2k(3)&
         -1091d0/630d0*B2k(4)+logNbar*B2k(4))*NNm1**9
    res=res+37d0/16d0*B2k(4)*NNm1**10
    res=res+(-31d0/2d0*B2k(2)*B2k(3)-65d0/8d0*B2k(1)*B2k(4)&
         -53723d0/27720d0*B2k(5)+logNbar*B2k(5))*NNm1**11
    res=res+14d0/5d0*B2k(5)*NNm1**12
    res=res+(-22d0*B2k(3)**2-253d0/8d0*B2k(2)*B2k(4)&
         -58d0/5d0*B2k(1)*B2k(5)-760223d0/360360d0*B2k(6)&
         +logNbar*B2k(6))*NNm1**13
    res=res+79d0/24d0*B2k(6)*NNm1**14
    res=res+(-2717d0/24d0*B2k(3)*B2k(4)-1183d0/20d0*B2k(2)*B2k(5)&
         -63d0/4d0*B2k(1)*B2k(6)-813089d0/360360d0*B2k(7)+logNbar*B2k(7))*NNm1**15
    res=res+53d0/14d0*B2k(7)*NNm1**16
    res=res+(-715d0/4d0*B2k(4)**2-3952d0/15d0*B2k(3)*B2k(5)&
         -308d0/3d0*B2k(2)*B2k(6)-144d0/7d0*B2k(1)*B2k(7)&
         -7303577d0/3063060d0*B2k(8)+logNbar*B2k(8))*NNm1**17
    res=res+137d0/32d0*B2k(8)*NNm1**18
    res=res+(-20111d0/20d0*B2k(4)*B2k(5)-561d0*B2k(3)*B2k(6)&
         -1173d0/7d0*B2k(2)*B2k(7)-417d0/16d0*B2k(1)*B2k(8)&
         -581523277d0/232792560d0*B2k(9)+logNbar*B2k(9))*NNm1**19
    res=res+43d0/9d0*B2k(9)*NNm1**20
    res=2d0*zeta3-res
    RETURN
  END FUNCTION Harmonic_Sum_S21

  !***********************************************************************************************************************************
  !  RIEMANNZETA
  !
  !  Riemann zeta function.
  !
  !  Algorithm from "Atlas for Computing Mathematical Functions" by W.J. Thompson, Wiley, 1997.
  !***********************************************************************************************************************************
  
  FUNCTION RIEMANNZETA(S,EPS0)
    ! Note that this is zeta(s) for s > 1
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::S
    REAL(KIND(1d0)),INTENT(IN),OPTIONAL::EPS0
    REAL(KIND(1d0))::EPS
    REAL(KIND(1d0))::RIEMANNZETA
    REAL(KIND(1d0))::NSTERM, SUM, FN, NEGS
    INTEGER::N,K
    IF(S.LE.1d0)THEN
       WRITE(*,*)"ERROR: S<=1 in RIEMANNZETA"
       STOP
    ENDIF
    IF(PRESENT(EPS0))THEN
       EPS=EPS0
    ELSE
       EPS=1D-14
    ENDIF
    
    !     Estimate N for accuracy  eps
    
    NSTERM = S*(S+1.0D00)*(S+2.0D00)* &
         (S+3.0D00)*(S+4.0D00)/30240.0D00
    N =INT((NSTERM*(2.0D00**S)/EPS)**(1d0/(S+5.0D00)))
    IF(N.LT.10)THEN
       N = 10
    END IF

    FN = N
    NEGS = -S
    !     Direct sum
    SUM = 0.0D00
    DO K =2, N-1
       SUM = SUM+K**NEGS
    END DO
    
    !     Add Euler-Maclaurin correction terms
    SUM = SUM+(FN**NEGS)*(0.5D00+FN/(S-1.0D00) &
         +S*(1D0-(S+1D0)*(S+2.0D00)/ &
         (60.0D00*FN*FN)) &
         /(12.0D00*FN))+NSTERM/(FN**(S+5.0D00))
    RIEMANNZETA=SUM+1d0
    RETURN
  END FUNCTION RIEMANNZETA
END MODULE HarmonicSums
