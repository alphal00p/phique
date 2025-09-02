*     -*-fortran-*-

      integer max_histograms
      parameter (max_histograms=1000)
      integer max_histogram_nwgt
      parameter (max_histogram_nwgt=500)
      logical hwu_yes,top_yes,gnu_yes,root_yes
      integer ihist_booked,hist_nwgt,hist_maxid
      character*50 hist_info(max_histogram_nwgt)
      integer hist_book_id(max_histograms)
      integer hist_book_nbin1(max_histograms)
      integer hist_book_nbin2(max_histograms)
      integer hist_book_type(max_histograms)
      character*3 hist_book_scale(max_histograms)
      character*50 hist_book_title(max_histograms)
      double precision hist_book_dx1(max_histograms)
      double precision hist_book_dx2(max_histograms)
      double precision hist_book_x1min(max_histograms)
      double precision hist_book_x2min(max_histograms)
      double precision hist_book_x1max(max_histograms)
      double precision hist_book_x2max(max_histograms)

      common/histogram_common1/top_yes,gnu_yes,root_yes,hwu_yes
      common/histogram_common2/ihist_booked,hist_book_id,hist_book_nbin1
     $     ,hist_book_nbin2,hist_book_type,hist_book_title
      common/histogram_common3/hist_book_dx1,hist_book_dx2,
     $     hist_book_x1min,hist_book_x2min,
     $     hist_book_x1max,hist_book_x2max
      common/histogram_common4/hist_book_scale,
     $     hist_nwgt,hist_info,hist_maxid
