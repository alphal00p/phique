MODULE plot_NNLO
  USE global
  IMPLICIT NONE
  INTEGER::max_weight=1
  INTEGER::nwgt_analysis
  INTEGER::num_plots=15 ! num_plots=inum_plots*3
  INTEGER::inum_plots=5
  INTEGER,PARAMETER::n_width=1,n_height=1
CONTAINS
  SUBROUTINE initplot_NNLO
    IMPLICIT NONE
    INTEGER::nwgt
    CHARACTER(len=15),DIMENSION(:),ALLOCATABLE::weights_info
    INTEGER::i
    max_weight=1+ho_nscale
    IF(ALLOCATED(weights_info))THEN
       DEALLOCATE(weights_info)
    ENDIF
    ALLOCATE(weights_info(max_weight))
    nwgt=1
    weights_info(nwgt)="central value  "
    IF(reweight_scale)THEN
       nwgt=nwgt+ho_nscale
       IF(ho_nscale.NE.2)THEN
          WRITE(*,*) 'ERROR #1 in initplot_NNLO:',ho_nscale
          STOP
       ENDIF
       WRITE(weights_info(nwgt-1),'(a4,f3.1,a2)')&
            "muR=",rw_Rscale_up
       WRITE(weights_info(nwgt  ),'(a4,f3.1,a2)')&
            "muR=",rw_Rscale_down
    END IF
    ! output plot files
    CALL plot_begin_NNLO(nwgt,max_weight,weights_info)
    RETURN
  END SUBROUTINE initplot_NNLO

  SUBROUTINE plotout_NNLO
    USE MC_VEGAS
    IMPLICIT NONE
    LOGICAL::usexinteg=.FALSE.,mint=.FALSE. ! I have not implemented MINT 
    LOGICAL::useitmax=.TRUE.
    REAL(KIND(1d0))::xnorm
    xnorm=1.d0
    xnorm=xnorm/float(itmx) ! use vegas with itmax
    ! output plot files
    CALL plot_end_NNLO(xnorm)
    RETURN
  END SUBROUTINE plotout_NNLO

  SUBROUTINE outfun_NNLO(iord,www)
    ! iord: 0, LO; 1, NLO; 2, NNLO
    IMPLICIT NONE
    INTEGER,INTENT(IN)::iord
    REAL(KIND(1d0)),INTENT(IN)::www
    INTEGER::i,j
    INTEGER::nwgt
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::wgts !,wgtden
    INTEGER::init=0
    SAVE init,wgts
    IF(www.EQ.0d0)RETURN
    IF(init.EQ.0)THEN
       max_weight=1+ho_nscale
       IF(ALLOCATED(wgts))THEN
          DEALLOCATE(wgts)
       ENDIF
       ALLOCATE(wgts(max_weight))
       init=1
    ENDIF
    nwgt=1
    wgts(1)=www
    IF(reweight_scale)THEN
       DO i=2,3
          nwgt=nwgt+1
          IF(iord.EQ.0)THEN
             wgts(nwgt)=www*LO_wgtxsecmu(i)
          ELSEIF(iord.EQ.1)THEN
             wgts(nwgt)=www*NLO_wgtxsecmu(i)
          ELSE
             wgts(nwgt)=www*NNLO_wgtxsecmu(i)
          ENDIF
       ENDDO
    ENDIF
    ! output plot file 
    CALL plot_fill_NNLO(iord,wgts)
    RETURN
  END SUBROUTINE outfun_NNLO

  SUBROUTINE plot_begin_NNLO(nwgt,max_weight0,weights_info)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nwgt,max_weight0
    CHARACTER(len=*),DIMENSION(max_weight0),INTENT(IN)::weights_info
    INTEGER::i,j,kk,l
    REAL(KIND(1d0)),PARAMETER::PI=3.14159265358979312D0
    character(len=5)::ord(3)=(/'   LO','  NLO',' NNLO'/)
    INCLUDE '../hbook/dbookf90.inc'
    INCLUDE '../hbook/histogram_common90.inc'
    top_yes=topdrawer_output
    gnu_yes=gnuplot_output
    root_yes=root_output
    hwu_yes=hwu_output
    CALL HISTOGRAM_INIT(nwgt,weights_info)
    nwgt_analysis=nwgt
    IF(nwgt_analysis*16.GT.nplots/4)THEN
       WRITE(*,*) 'error in plot_begin_NNLO: ',&
            nwgt_analysis*num_plots*4 ! 4 is internal in dbook
       STOP
    ENDIF
    DO i=1,3
       l=(i-1)*inum_plots
       CALL histogram_book(l+ 1,'total rate '//ord(i),&
            2,0d0,4d0,.FALSE.)
       CALL histogram_book(l+ 2,'W(aa) '//ord(i),&
            100,0d0,10000d0,.TRUE.)
       CALL histogram_book(l+ 3,'W(aa) zoom in '//ord(i),&
            100,0d0,1000d0,.TRUE.)
       CALL histogram_book(l+ 4,'W(aa) zoom in 2 '//ord(i),&
            100,0d0,100d0,.TRUE.)
       CALL histogram_book(l+ 5,'W(aa) zoom in 3 '//ord(i),&
            100,0d0,10d0,.TRUE.)
    ENDDO

    CALL HISTOGRAM_BOOK_FINISH
    RETURN
  END SUBROUTINE plot_begin_NNLO

  SUBROUTINE plot_end_NNLO(xnorm)
    IMPLICIT NONE
    CHARACTER(len=14)::ytit
    REAL(KIND(1d0)),INTENT(IN)::xnorm
    INTEGER::i
    INTEGER::kk,l
    INCLUDE '../hbook/dbookf90.inc'
    CHARACTER(len=400)::psfile,rootfile
    IF(topdrawer_output.OR.gnuplot_output.OR.root_output)THEN
       CALL mclear
       CALL mclear_3d
       DO i=1,NPLOTS
          CALL mopera(i,'+',i,i,xnorm,0.d0)
          CALL mfinal(i)
          CALL mopera_3d(i,'+',i,i,xnorm,0.d0)
          CALL mfinal_3d(i)
       ENDDO
    ENDIF
    ytit='sigma per bin '
    IF(topdrawer_output)THEN
       CALL open_topdrawer_file_NNLO
       CALL HISTOGRAM_END(1,ytit,' ',' ')
       CALL close_topdrawer_file_NNLO
    ENDIF
    IF(gnuplot_output)THEN
       CALL open_gnuplot_file_NNLO
       psfile="aa2QQbar_NNLO.ps"
       CALL HISTOGRAM_END(2,ytit,psfile,' ')
       CALL close_gnuplot_file_NNLO
    ENDIF
    IF(root_output)THEN
       CALL open_root_file_NNLO
       rootfile="aa2QQbar_NNLO.root"
       CALL HISTOGRAM_END(3,ytit,' ',rootfile)
       ! wirte the end of the .C file for root
       WRITE (96, 100)
100    FORMAT(/1X,&
            ' hohisto -> cd();',/1X,&
            ' if (histos -> GetEntries() > 0 ) then {',/1X,&
            '  histos->Write();',/1X,&
            '  hohisto -> Close();',/1X,&
            ' }',/1X,'}')
       CALL close_root_file_NNLO
    ENDIF
    IF(hwu_output)THEN
       CALL HISTOGRAM_END(4,ytit,' ',' ')
       CALL HwU_write_file_NNLO
    ENDIF
    RETURN
  END SUBROUTINE plot_end_NNLO

  SUBROUTINE plot_fill_NNLO(iord,wgts)
    IMPLICIT NONE
    ! iord: 0, LO; 1, NLO; 2, NNLO
    INTEGER,INTENT(IN)::iord
    REAL(KIND(1d0)),DIMENSION(:),INTENT(IN)::wgts
    REAL(KIND(1d0))::var,Waa,www
    INTEGER::i,l
    IF(wgts(1).EQ.0d0)RETURN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! DEFINE OBSERVABLES
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    var=1.0d0 ! total cross section
    Waa=Wgaga
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! START TO FILL HISTOGRAM
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO i=1,3
       IF(iord+1.NE.i)CYCLE
       l=(i-1)*inum_plots
       CALL HISTOGRAM_fill(l+1,var,wgts)
       CALL HISTOGRAM_fill(l+2,Waa,wgts)
       CALL HISTOGRAM_fill(l+3,Waa,wgts)
       CALL HISTOGRAM_fill(l+4,Waa,wgts)
       CALL HISTOGRAM_fill(l+5,Waa,wgts)
    ENDDO
    IF(hwu_output)THEN
       call HwU_add_points
    ENDIF
    ! boost back to c.m. frame
    RETURN
  END SUBROUTINE plot_fill_NNLO
  
  SUBROUTINE GetPTOrdered(npt,pt_ordered)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::npt
    REAL(KIND(1d0)),DIMENSION(npt),INTENT(INOUT)::pt_ordered
    INTEGER::i,j
    REAL(KIND(1d0)),DIMENSION(npt)::pt_orig
    LOGICAL,DIMENSION(npt)::used
    REAL(KIND(1d0))::ptmax
    INTEGER::indexmax
    pt_orig(1:npt)=pt_ordered(1:npt)
    used(1:npt)=.FALSE.
    DO i=1,npt
       ptmax=-1d0
       indexmax=0
       DO j=1,npt
          IF(used(j))CYCLE
          IF(pt_orig(j).GT.ptmax)THEN
             ptmax=pt_orig(j)
             indexmax=j
          ENDIF
       ENDDO
       pt_ordered(i)=ptmax
       used(indexmax)=.TRUE.
    ENDDO
    RETURN
  END SUBROUTINE GetPTOrdered

  SUBROUTINE HwU_write_file_NNLO
    IMPLICIT NONE
    REAL(KIND(1d0))::xnorm
    OPEN(unit=599,file=TRIM(output_dir)//'aa2QQbar_NNLO.HwU',&
         status='unknown')
    xnorm=1d0
    CALL HwU_output(599,xnorm)
    CLOSE(599)
    RETURN
  END SUBROUTINE HwU_write_file_NNLO
  
  SUBROUTINE open_topdrawer_file_NNLO
    IMPLICIT NONE
    LOGICAL::useitmax
    OPEN(unit=99,FILE=TRIM(output_dir)//'aa2QQbar_NNLO.top',&
         status='unknown')
    useitmax=.FALSE.
  END SUBROUTINE open_topdrawer_file_NNLO

  SUBROUTINE close_topdrawer_file_NNLO
    IMPLICIT NONE
    CLOSE(99)
    RETURN
  END SUBROUTINE close_topdrawer_file_NNLO

  SUBROUTINE open_gnuplot_file_NNLO
    IMPLICIT NONE
    LOGICAL::useitmax
    OPEN(unit=97,FILE=TRIM(output_dir)//'aa2QQbar_NNLO.gnu',&
         status='unknown')
    useitmax=.FALSE.
  END SUBROUTINE open_gnuplot_file_NNLO

  SUBROUTINE close_gnuplot_file_NNLO
    IMPLICIT NONE
    CLOSE(97)
    RETURN
  END SUBROUTINE close_gnuplot_file_NNLO

  SUBROUTINE open_root_file_NNLO
    USE global
    IMPLICIT NONE
    LOGICAL::useitmax
    OPEN(unit=96,FILE=TRIM(output_dir)//'aa2QQbar_NNLO.C',&
         status='unknown')
    useitmax=.FALSE.
  END SUBROUTINE open_root_file_NNLO

  SUBROUTINE close_root_file_NNLO
    IMPLICIT NONE
    CLOSE(96)
    RETURN
  END SUBROUTINE close_root_file_NNLO
END MODULE plot_NNLO
