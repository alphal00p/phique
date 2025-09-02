      SUBROUTINE HISTOGRAM_INIT(NWGT,WEIGHTS_INFO)
      implicit none
      include "histogram_common.inc"
      integer nwgt
      character*(*) weights_info(*)
      integer i
      IF(top_yes.or.gnu_yes.or.root_yes)THEN
         CALL INIHIST
         CALL INIHIST_3D
         IF(nwgt.GT.max_histogram_nwgt)THEN
            WRITE(*,*)"ERROR:please enlarge the number of "//
     $           "max_histogram_nwgt in histograms_common.inc"//
     $           " and histograms_common90.inc"
            STOP
         ENDIF
         hist_nwgt=nwgt
         do i=1,nwgt
            hist_info(i)=WEIGHTS_INFO(i)
         enddo
      ENDIF
      IF(hwu_yes)THEN
         call HwU_inithist(nwgt,weights_info)
         ! Set method for error estimation to '0', i.e., use Poisson statistics
         ! for the uncertainty estimate
         call set_error_estimation(0)
      ENDIF
      ihist_booked=0
      do i=1,max_histograms
         hist_maxid=0
         hist_book_id(i)=0
         hist_book_nbin1(i)=0
         hist_book_nbin2(i)=0
         hist_book_type(i)=0
         hist_book_title(i)=''
         hist_book_dx1(i)=0d0
         hist_book_dx2(i)=0d0
         hist_book_x1min(i)=0d0
         hist_book_x1max(i)=0d0
         hist_book_x2min(i)=0d0
         hist_book_x2max(i)=0d0
         hist_book_scale(i)='LIN'
      enddo
      END

      SUBROUTINE HISTOGRAM_BOOK(id,title,nbin,xmin,xmax,logscale)
      implicit none
      integer id,nbin
      character*(*) title
      double precision xmin,xmax
      logical logscale
      include "histogram_common.inc"
      ihist_booked=ihist_booked+1
      if(ihist_booked.GT.max_histograms)then
         write(*,*)"ERROR:please enlarge the number of "//
     $        "max_histograms in histograms_common.inc"//
     $        " and histograms_common90.inc"
         stop
      endif
      hist_maxid=max(hist_maxid,id)
      hist_book_dx1(ihist_booked)=(xmax-xmin)/DBLE(nbin)
      hist_book_id(ihist_booked)=id
      if(logscale)then
         hist_book_scale(ihist_booked)='LOG'
      else
         hist_book_scale(ihist_booked)='LIN'
      endif
      hist_book_title(ihist_booked)=title
      hist_book_nbin1(ihist_booked)=nbin
      hist_book_x1min(ihist_booked)=xmin
      hist_book_x1max(ihist_booked)=xmax
      hist_book_type(ihist_booked)=1 ! 1d plot
      if(hwu_yes)then
         call HwU_book(id,title,nbin,xmin,xmax)
      endif
      if(top_yes.or.gnu_yes.or.root_yes)then
         call bookup(id,title,hist_book_dx1(ihist_booked),xmin,xmax)
      endif
      END

      SUBROUTINE HISTOGRAM_BOOK3D(id,title,nbin1,x1min,x1max,
     $     nbin2,x2min,x2max,logscale)
      implicit none
      integer id,nbin1,nbin2
      character*(*) title
      double precision x1min,x1max,x2min,x2max
      integer init
      data init/0/
      save init
      logical logscale
      include "histogram_common.inc"
      ihist_booked=ihist_booked+1
      if(ihist_booked.GT.max_histograms)then
          write(*,*)"ERROR:please enlarge the number of "//
     $        "max_histograms in histograms_common.inc"
           stop
      endif
      hist_book_dx1(ihist_booked)=(x1max-x1min)/DBLE(nbin1)
      hist_book_dx2(ihist_booked)=(x2max-x2min)/DBLE(nbin2)
      hist_book_id(ihist_booked)=id
      if(logscale)then
         hist_book_scale(ihist_booked)='LOG'
      else
         hist_book_scale(ihist_booked)='LIN'
      endif
      hist_book_title(ihist_booked)=title
      hist_book_nbin1(ihist_booked)=nbin1
      hist_book_nbin2(ihist_booked)=nbin2
      hist_book_x1min(ihist_booked)=x1min
      hist_book_x1max(ihist_booked)=x1max
      hist_book_x2min(ihist_booked)=x2min
      hist_book_x2max(ihist_booked)=x2max
      hist_book_type(ihist_booked)=2 ! 2d plot
      if(gnu_yes.or.root_yes)then
         hist_maxid=max(hist_maxid,id)
         call bookup_3d(id,title,hist_book_dx1(ihist_booked),
     $        x1min,x1max,hist_book_dx2(ihist_booked),x2min,x2max)
       else
          if(init.eq.0)then
             WRITE(*,*)"WARNING:2D plot can be kept when gnu_yes=T"
             WRITE(*,*)"WARNING:or root_yes=T,we will skip this 2D plot"
             init=1
          endif
      endif
      END

      SUBROUTINE HISTOGRAM_BOOK_FINISH
      implicit none
      integer i,l,j,newid
      character*50 newtitle
      include "histogram_common.inc"
      if(top_yes.or.gnu_yes.or.root_yes)then
         if(hist_nwgt.GT.1.and.ihist_booked.ge.1)then
            do i=2,hist_nwgt
               l=(i-1)*hist_maxid
               do j=1,ihist_booked
                  newid=hist_book_id(j)+l
                  newtitle=TRIM(hist_book_title(j))//" "//hist_info(i)
                  if(hist_book_type(j).EQ.2)THEN
                     call bookup_3d(newid,newtitle,hist_book_dx1(j),
     $                    hist_book_x1min(j),hist_book_x1max(j),
     $                    hist_book_dx2(j),hist_book_x2min(j),
     $                    hist_book_x2max(j))
                  else
                     call bookup(newid,newtitle,hist_book_dx1(j),
     $                    hist_book_x1min(j),hist_book_x1max(j))
                  endif
               enddo
            enddo
         endif
      endif
      END

      SUBROUTINE HISTOGRAM_END(itype,ytit,gnu_file,root_file)
      implicit none
      ! itype: 1, topdrawer, 2, gnuplot, 3, root, 4, hwu
      integer itype
      character*(*) ytit,gnu_file,root_file
      integer i,j,l,newid
      character*50 newtitle
      include "histogram_common.inc"
      if(top_yes.and.itype.eq.1)then
         do i=1,ihist_booked
            if(hist_book_type(i).ne.1)cycle
            call multitop(hist_book_id(i),1,1,
     $           hist_book_title(i),ytit,hist_book_scale(i))
         enddo
         if(hist_nwgt.GT.1.and.ihist_booked.ge.1)then
            do i=2,hist_nwgt
               l=(i-1)*hist_maxid
               do j=1,ihist_booked
                  if(hist_book_type(j).ne.1)cycle
                  newid=hist_book_id(j)+l
                  newtitle=TRIM(hist_book_title(j))//" "//hist_info(i)
                  call multitop(newid,1,1,newtitle,ytit,
     $                 hist_book_scale(i))
               enddo
            enddo
         endif
      endif
      if(gnu_yes.and.itype.eq.2)then
         do i=1,ihist_booked
            if(hist_book_type(i).eq.2)then
               call MGNUPLOT_3D(hist_book_id(i),hist_book_title(i),
     $              hist_book_title(i),ytit,hist_book_scale(i),gnu_file)
            else
               call MGNUPLOT(hist_book_id(i),hist_book_title(i),
     $              ytit,hist_book_scale(i),gnu_file)
            endif
         enddo
         if(hist_nwgt.GT.1.and.ihist_booked.ge.1)then
            do i=2,hist_nwgt
               l=(i-1)*hist_maxid
               do j=1,ihist_booked
                  newid=hist_book_id(j)+l
                  newtitle=TRIM(hist_book_title(j))//" "//hist_info(i)
                  if(hist_book_type(j).eq.2)then
                     call MGNUPLOT_3D(newid,newtitle,newtitle,
     $                    ytit,hist_book_scale(j),gnu_file)
                  else
                     call MGNUPLOT(newid,newtitle,ytit,
     $                    hist_book_scale(i),gnu_file)
                  endif
               enddo
            enddo
         endif
      endif
      if(root_yes.and.itype.eq.3)then
         do i=1,ihist_booked
            if(hist_book_type(i).eq.2)then
               call MROOTPLOT_3D(hist_book_id(i),hist_book_title(i),
     $              hist_book_title(i),ytit,
     $              hist_book_scale(i),root_file)
            else
               call MROOTPLOT(hist_book_id(i),hist_book_title(i),
     $              ytit,hist_book_scale(i),root_file)
            endif
         enddo
         if(hist_nwgt.GT.1.and.ihist_booked.ge.1)then
            do i=2,hist_nwgt
               l=(i-1)*hist_maxid
               do j=1,ihist_booked
                  newid=hist_book_id(j)+l
                  newtitle=TRIM(hist_book_title(j))//" "//hist_info(i)
                  if(hist_book_type(j).eq.2)then
                     call MROOTPLOT_3D(newid,newtitle,newtitle,
     $                    ytit,hist_book_scale(j),root_file)
                  else
                     call MROOTPLOT(newid,newtitle,ytit,
     $                    hist_book_scale(i),root_file)
                  endif
               enddo
            enddo
         endif
      endif
      if(hwu_yes.and.itype.eq.4)then
         call finalize_histograms(1)
      endif
      return
      END

      SUBROUTINE HISTOGRAM_FILL(id,var,wgts)
      implicit none
      integer id
      double precision var
      double precision wgts(*)
      double precision www
      integer i,l,j,newid
      include "histogram_common.inc"
      if(top_yes.or.gnu_yes.or.root_yes)then
         www=wgts(1)
         CALL mfill(id,var,www)
         if(hist_nwgt.GT.1)then
            do i=2,hist_nwgt
               l=(i-1)*hist_maxid
               newid=id+l
               www=wgts(i)
               CALL mfill(newid,var,www)
            enddo
         endif
      endif
      if(hwu_yes)then
         CALL HwU_fill(id,var,wgts)
      endif
      return
      END

      SUBROUTINE HISTOGRAM_FILL3D(id,var1,var2,wgts)
      implicit none
      integer id
      double precision var1,var2
      double precision wgts(*)
      double precision www
      integer i,l,j,newid
      include "histogram_common.inc"
      if(gnu_yes.or.root_yes)then
         www=wgts(1)
         CALL mfill_3d(id,var1,var2,www)
         if(hist_nwgt.GT.1)then
            do i=2,hist_nwgt
               l=(i-1)*hist_maxid
               newid=id+l
               www=wgts(i)
               CALL mfill_3d(newid,var1,var2,www)
            enddo
         endif
      endif
      RETURN
      END
