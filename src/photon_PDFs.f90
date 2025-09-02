MODULE photon_PDFs
  USE global
  IMPLICIT NONE
CONTAINS
  !/* ********************************************************* */
  !/*  Equivalent photon approximation structure function.   * */
  !/*     Improved Weizsaecker-Williams formula              * */
  !/*   V.M.Budnev et al., Phys.Rep. 15C (1975) 181          * */
  !/* ********************************************************* */ 
  !   provided by Tomasz Pierzchala - UCL
  
  FUNCTION epa_electron(x,q2max)
    ! photon PDF from electron
    IMPLICIT NONE
    REAL(KIND(1d0))::epa_electron
    REAL(KIND(1d0)),INTENT(IN)::x,q2max
    REAL(KIND(1d0))::xin=0.511d-3 ! electron mass in GeV
    REAL(KIND(1d0))::alpha
    REAL(KIND(1d0))::f,q2min
    REAL(KIND(1d0)),PARAMETER::PI=3.14159265358979323846d0

    alpha=1d0/alphaemm1

    !     // x = omega/E = (E-E')/E
    IF(x.LT.1)THEN
       q2min= xin*xin*x*x/(1-x)
       IF(q2min.LT.q2max)THEN
          f = alpha/2d0/PI*&
               (2d0*xin*xin*x*(-1/q2min+1/q2max)+&
               (2-2d0*x+x*x)/x*dlog(q2max/q2min))
       ELSE
          f = 0d0
       ENDIF
    ELSE
       f= 0d0
    ENDIF
    !      write (*,*) x,dsqrt(q2min),dsqrt(q2max),f
    IF(f.LT.0d0)f = 0d0
    epa_electron= f
    
  END FUNCTION epa_electron

  ! similar as Eq.(3) in 0909.3047
  ! also see Eq.(6) in nucl-ex/0502005
  FUNCTION epa_proton(x,q2max)
    ! photon PDF from elastic proton
    IMPLICIT NONE
    REAL(KIND(1d0))::epa_proton
    REAL(KIND(1d0)),INTENT(IN)::x,q2max
    REAL(KIND(1d0))::xin=0.938d0 ! proton mass in GeV
    REAL(KIND(1d0))::alpha,qz
    REAL(KIND(1d0))::f,qmi
    REAL(KIND(1d0)),PARAMETER::PI=3.14159265358979323846d0
    
    alpha = 1d0/alphaemm1

    qz = 0.71d0

    !     // x = omega/E = (E-E')/E
    IF(x.LT.1)THEN
       qmi= xin*xin*x*x/(1-x)
       IF(qmi.LT.q2max)THEN
          f = alpha/PI*(phi_f(x,q2max/qz)-phi_f(x,qmi/qz))*(1-x)/x
       ELSE
          f=0d0
       ENDIF
    ELSE
       f= 0d0
    ENDIF
    IF(f.LT.0d0) f = 0d0
    epa_proton= f
  END FUNCTION epa_proton

  FUNCTION phi_f(x,qq)
    IMPLICIT NONE
    REAL(KIND(1d0))::phi_f
    REAL(KIND(1d0)),INTENT(IN)::x,qq
    REAL(KIND(1d0))::y,qq1,f,a,b,c

    a = 7.16d0
    b = -3.96d0
    c = 0.028d0
    
    qq1=1+qq
    y= x*x/(1-x)
    f=(1+a*y)*(-DLOG(qq1/qq)+1/qq1+1/(2*qq1*qq1)+1/(3*qq1*qq1*qq1))
    f=f + (1-b)*y/(4*qq*qq1*qq1*qq1);
    f=f+ c*(1+y/4)*(DLOG((qq1-b)/qq1)+b/qq1+b*b/(2*qq1*qq1)+&
         b*b*b/(3*qq1*qq1*qq1))
    phi_f= f
  END FUNCTION phi_f

END MODULE photon_PDFs
