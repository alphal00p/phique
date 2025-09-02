! $Id: newton_method.f90,v 1.0 2018/06/21 18:34:17 hshao Exp $
module newton_method
  implicit none
  private

  interface newtonsolver
     module procedure newtonsolver_linear, newtonsolver_quadratic
  end interface
  public :: newtonsolver


contains
  !----------------------------------------------------------------------
  ! write the routines for the Newton's method
  ! This is a linear version of Newton's method
  subroutine newtonsolver_linear(x0,fpx0,f_series,eps,nmax)
    real(kind(1d0)), intent(in)    :: eps
    real(kind(1d0)), intent(out) :: fpx0 ! f'(x0)
    integer, intent(in) :: nmax
    real(kind(1d0)), intent(inout) :: x0
    interface
       subroutine f_series(x,f,fp)
         real(kind(1d0)), intent(in)  :: x
         real(kind(1d0)), intent(out) :: f,fp
       end subroutine f_series
    end interface
    !------------------------------------------------------------
    integer :: i
    logical :: done
    real(kind(1d0)) :: f, fp, x1
    if(nmax.LE.0.and.eps.Le.0d0)then
       write(*,*)"ERROR: cannot do Newton's method with nmax<=0 and eps <= 0"
       stop
    endif
    done=.FALSE.
    i=0
    do while(.NOT.done)
       call f_series(x0,f,fp)
       x1=x0-f/fp
       i=i+1
       if(nmax.GT.0.and.i.GT.nmax)then
          x0=x1
          done=.TRUE.
       elseif(eps.GT.0d0.and.dabs(x1-x0).LT.eps)then
          x0=x1
          done=.TRUE.
       elseif(isnan(x0).or.1d0/x0.EQ.0d0)then
          done=.TRUE.
       else
          x0=x1
       endif
    enddo
    ! since x0 is very close to x1
    ! take the one has been calculated (improve the efficiency)
    ! If one needs to be more precise, uncomment the below line.
    ! call f_series(x0, f,fp)
    fpx0=fp
    return
  end subroutine newtonsolver_linear

  !----------------------------------------------------------------------
  ! This is a quadratic version of Newton's method
  subroutine newtonsolver_quadratic(x0,fpx0,f_series,eps,nmax)
    real(kind(1d0)), intent(in)    :: eps
    real(kind(1d0)), dimension(2), intent(out) :: fpx0
    integer, intent(in) :: nmax
    real(kind(1d0)), intent(inout) :: x0
    interface
       subroutine f_series(x,f,fp,fpp)
         real(kind(1d0)), intent(in)  :: x
         real(kind(1d0)), intent(out) :: f,fp,fpp
       end subroutine f_series
    end interface
    !------------------------------------------------------------
    integer :: i
    logical :: done
    real(kind(1d0)) :: f, fp, fpp, x1
    if(nmax.LE.0.and.eps.Le.0d0)then
       write(*,*)"ERROR: cannot do Newton's method with nmax<=0 and eps <= 0"
       stop
    endif
    done=.FALSE.
    i=0
    do while(.NOT.done)
       call f_series(x0,f,fp,fpp)
       x1=x0-2d0*f*fp/(2d0*fp**2-f*fpp)
       i=i+1
       if(nmax.GT.0.and.i.GT.nmax)then
          x0=x1
          done=.TRUE.
       elseif(eps.GT.0d0.and.dabs(x1-x0).LT.eps)then
          x0=x1
          done=.TRUE.
       elseif(isnan(x0).or.1d0/x0.EQ.0d0)then
          done=.TRUE.
       else
          x0=x1
       endif
    enddo
    ! since x0 is very close to x1
    ! take the ones been calculated (improve the efficiency)
    ! If one needs to be more precise, uncomment the below line.
    ! call f_series(x0,f,fp,fpp)
    fpx0(1)=fp
    fpx0(2)=fpp
    return
  end subroutine newtonsolver_quadratic

end module newton_method
