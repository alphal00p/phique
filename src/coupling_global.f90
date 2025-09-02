MODULE coupling_global
  USE qcd_coupling
  IMPLICIT NONE
  ! QCD colour factors
  REAL(KIND(1d0)),PARAMETER::QCD_CA=3d0,QCD_CF=1.333333333333333333333333333333d0,QCD_TF=0.5d0
  ! running of alpha_s
  ! the threshold for changing the flavours
  REAL(KIND(1d0)),DIMENSION(6)::Q_THR_default=(/1d-10,1d-10,0.15d0,&
       1.5d0,4.7d0,173.3d0 /)
  TYPE(as_container)::as_box
END MODULE coupling_global
