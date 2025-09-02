C************************************************************************
C
C                           3D HISTOGRAMMING PACKAGE 
C                              Hua-Sheng Shao
C
C************************************************************************
C
C**********************************************************************
C    SIMPLE HISTOGRAMMING PACKAGE --  SIMPLIFIED VERSION OF HBOOK
C    BY Hua-Sheng Shao    June 2014
C**********************************************************************
C
C Fills up to 500 histograms with up to 100*100 bins. 
C Gives a data file (to be specified in the calling program by assigning 
C a file name to unit 98) and a topdrawer file (to be specified in the 
C calling program by assigning a file name to unit 99-topdrawer,97-gnuplot,96-root).
C
C INITIALIZATION:
C Call once INIHIST_3D; this just resets a few counters and logicals
C Call MBOOK_3D(N,'TITLE',XDEL,XMIN,XMAX,YDEL,YMIN,YMAX) for each histogram to be booked.
C N (an integer) is the label of the histogram;
C 'TITLE' is the name of the histogram (no more then 100 characters);
C XDEL,YDEL (real*8) is the bin size;
C XMIN,YMIN (real*8) is the lower limit of the first bin;
C XMAX,YMAX (real*8)is the upper limit of the last  bin
C Example:
C      call mbook_3d(2,'pt1-pt2 distribution',1.,10,70,1.,10,70)
C This call initializes histogram number 2, called 'pt1-pt2 distribution';
C The bin size will be 1.*1. (possibly GeV, if that's what you want), the
C first bin being  10<x<11 and 10<y<11. and the last one being 69.<x<70 and 69<y<70
C
C FILLING:
C When it's time, call MFILL4_3D(N,X,Y,Z); this will add Z (real*8) to the bin 
C in which X*Y (real*8) happens to be, within histogram N. 
C
C PLAYING AROUND:
C At the end of the day you may want to sum, divide, cancel, etc.etc.
C various histograms (bin by bin). Then you call MOPERA_3D(I,'O',J,K,X,Y). 
C The 1-character string O can take the following values:
C +  : sums       X*(hist I) with Y*(hist J) and puts the result in hist K;
C -  : subtracts  X*(hist I) with Y*(hist J) and puts the result in hist K;
C *  : multiplies X*(hist I) with Y*(hist J) and puts the result in hist K;
C /  : divides    X*(hist I) with Y*(hist J) and puts the result in hist K;
C F  : multiplies hist I by the factor X, and puts the result in hist K;
C R  : takes the square root of  hist  I, and puts the result in hist K;if
C      the value at a given bin is less than or equal to 0, puts 0 in K
C S  : takes the square      of  hist  I, and puts the result in hist K;
C L  : takes the log_10 of  hist  I, and puts the result in hist K; if the
C      value at a given bin is less than or equal to 0, puts 0 in K
C M  : statistical analysis; if I contains the weights (let's say WGT),
C      J contains variable times weight (F*WGT) and K contains the
C      variable squared times the weight (F**2*WGT), then, after using 'M',
C      J will contain the average value of the variable <F> and K will 
C      contain the sigma of the average: sigma=sqrt(<F**2>-<F>**2).
C      If WGT=1. for all the entries, then it is enough to put I=J, and
C      it is not necessary to book a hist with the weights.
C V  : estimates errors for vegas evaluation of differential distributions.
C      Fill I with the values of
C      the functions do integrate times the Vegas weight (fun*wgt); fill
C      J with fun**2*wgt; then K will contain an estimate of the error
C      of the integration. Putting X=1/(#of iterations) performs the 
C      avegare over the iterations, and gives the right normalization to 
C      the differential distribution, I, and to the errors, K. J stays the same.
C
C FINAL ACCOUNTING:
C Now we can finalize our histograms; MFINAL_3D(N) will calculate the integral
C of the histogram N, the mean value of the X variable and its RMS.
C If we now want to renormalize the hist's, we can call MNORM_3D(N,X), which
C will normalize the integral to X  -- CAUTION: do not call MNORM_3D before
C MFINAL_3D, it will blow up.
C
C OUTPUT:
C To get a .dat file containing the values of the histograms, together with
C some information (like integral, mean values, etc.etc.) call MPRINT_3D(N),
C for each hist N that you want in the .dat file. Before the call to MPRINT_3D
C you want to open unit 98 and give it a name:                       
C     OPEN(UNIT=98,NAME='NAME.DAT',STATUS='NEW')
C If you want a topdrawer file with a plot of the hist values, call 
C MTOP(N,M,'X','Y','SCALE'). The points of the plot will be taken from histogram
C N, the error bars from histogram M. 'SCALE', character*(*), determines
C the scale for y, logarithmic or linear (SCALE=LOG,LIN). 
C If you do not want error bars, keep
C a histogram of zeros, or just call a hist that had not been booked.
C X will appear as a 'bottom title', and Y will appear as a 'left title'.
C The top title is by default the name of the histogram itself.
C A little box below the plot will contain some information on the plot
C itself. Before calling MTOP,
C     OPEN(UNIT=99,NAME='NAME.TOP',STATUS='NEW')
C--------------------------------------------------------------------------
C
C  COMMON/HISTO/  Histogram N
C                           
C   BOOK(N),      Three-letter character-string: 'NO' if histogram was not 
C		  Booked, 'YES' otherwise
C   TITLE(N),     Title of the histogram
C
C   HMIN(N),      Min value of x range
C   HMAX(N),      Max value of x range
C   HDEL(N),      Bin width
C   NBIN(N),      Total number of bins
C   USCORE(N),    Total integral of underscores with x < HMIN(N)
C   OSCORE(N),    Total integral of onderscores with x > HMAX(N)
C   IUSCORE(N),   Number of entries with x < HMIN(N)
C   IOSCORE(N),   Number of entries with x > HMAX(N)
C   IENT(N),      Total number of entries within x range HMIN(N)<x<HMAX(N)
C   HINT(N),      Integral of the histogram within HMIN(N)<x<HMAX(N)
C   HAVG(N),      Average value of x, weighted over the x range of the histo
C   HSIG(N),      Quadratic dispersion of x around the average
C   HIST(N,L),    Value of bin L-th
C   XHIS(N,L),    Central x value of bin L-th
C   IHIS(N,L),    Number of entries within bin L-th
C   NHIST         Total number of booked histograms
C


      SUBROUTINE INIHIST_3D
C*************************************************************
c     initialization
      IMPLICIT NONE
      include 'dbook-3d.inc'
c
c     LOCAL
c
      INTEGER I
C
C     START
C
      NHIST=0
      DO 1, I=1,NPLOTS             
   1  BOOK(I)=' NO'
      END  
 

    
      SUBROUTINE MBOOK_3D(N,TIT,XDEL,XMIN,XMAX,YDEL,YMIN,YMAX)
C*************************************************************
      IMPLICIT NONE
C
C     ARGUMENTS
C
      INTEGER N
      CHARACTER*(*) TIT
      REAL*8 XDEL,XMIN,XMAX
      REAL*8 YDEL,YMIN,YMAX
C
C     GLOBAL
C
      include 'dbook-3d.inc'
C
C     LOCAL
C
      INTEGER I,J
Cq
C     START
C
      NHIST = MAX(N,NHIST)
      IF(BOOK(N)(1:1).EQ.'Y') THEN
	 CALL MWARN_3D('MBOOK_3D')
         WRITE(*,*) 'Histogram',N,TITLE(N),' already in use. '
         WRITE(*,*) 'superseded by ',TIT
      ENDIF
      BOOK(N) = 'YES'
      TITLE(N) = ' '//TIT
1     HXDEL(N) = XDEL
      HYDEL(N) = YDEL
      NXBIN(N) = INT((XMAX-XMIN)/(XDEL*0.999999d0))
      NYBIN(N) = INT((YMAX-YMIN)/(YDEL*0.999999d0))
      IF(NXBIN(N).GT.100) THEN
	WRITE(*,*) 'TOO MANY X BINS (',NXBIN(N),') REQUIRED IN HIST ',N
	WRITE(*,*) 'RE-ENTER BIN SIZE DELTA (OLD BIN = ',XDEL,' ):'
	READ(*,*) XDEL
	GO TO 1
      ENDIF
      IF(NYBIN(N).GT.100)THEN
         WRITE(*,*) 'TOO MANY X BINS (',NYBIN(N),') REQUIRED IN HIST ',N
         WRITE(*,*) 'RE-ENTER BIN SIZE DELTA (OLD BIN = ',YDEL,' ):'
         READ(*,*) YDEL
         GO TO 1
      ENDIF
      HXMIN(N) = XMIN
      HYMIN(N) = YMIN
      HXMAX(N) = NXBIN(N)*XDEL+XMIN
      HYMAX(N) = NYBIN(N)*YDEL+YMIN
      IF(abs(HXMAX(N)-XMAX).gt.0.001d0*XDEL) THEN
	 CALL MWARN_3D('MBOOK_3D')
         WRITE(*,*)
     #'Histogram ', TIT, ' Change of upper limit:',xmax,'-->',HXMAX(N)
      ENDIF
      IF(abs(HYMAX(N)-YMAX).gt.0.001d0*YDEL) THEN
         CALL MWARN_3D('MBOOK_3D')
         WRITE(*,*)
     #'Histogram ', TIT, ' Change of upper limit:',ymax,'-->',HYMAX(N)
      ENDIF
      IENT(N) = 0
      IUUSCORE(N) = 0
      IUOSCORE(N) = 0
      IOUSCORE(N) = 0
      IOOSCORE(N) = 0
      IUESCORE(N) = 0
      IOESCORE(N) = 0
      IEUSCORE(N) = 0
      IEOSCORE(N) = 0
      UUSCORE(N) = 0
      UOSCORE(N) = 0
      OUSCORE(N) = 0
      OOSCORE(N) = 0
      UESCORE(N) = 0
      OESCORE(N) = 0
      EUSCORE(N) = 0
      EOSCORE(N) = 0
      HAVG(N) = 0
      HINT(N) = 0
      HSIG(N) = 0
      DO I=1,NXBIN(N)
         XHIS(N,I)=HXMIN(N)+HXDEL(N)*(DFLOAT(I)-0.5d0)
      ENDDO
      DO I=1,NYBIN(N)
         YHIS(N,I)=HYMIN(N)+HYDEL(N)*(DFLOAT(I)-0.5d0)
      ENDDO
      DO I=1,NXBIN(N)
         DO J=1,NYBIN(N)
            IHIS(N,I,J)=0
            HIST(N,I,J)=0
         ENDDO
      ENDDO
      END

      SUBROUTINE MFILL4_3D(N,X,Y,Z)
C*************************************************************
      IMPLICIT NONE
C
C     ARGUMENTS
C
      INTEGER N
      REAL*8 X,Y,Z
C
C     LOCAL
C
      REAL*8 XI,YI
      INTEGER I,J
C
C     GLOBAL
C
      INCLUDE 'dbook-3d.inc'
c
c     START
c

      IF(X.LT.HXMIN(N).AND.Y.LT.HYMIN(N)) THEN
         UUSCORE(N) = UUSCORE(N) + Z
         IUUSCORE(N) = IUUSCORE(N) + 1
      ELSEIF(X.GT.HXMAX(N).AND.Y.LT.HYMIN(N)) THEN
         OUSCORE(N) = OUSCORE(N) + Z
         IOUSCORE(N) = IOUSCORE(N) + 1
      ELSEIF(X.LT.HXMIN(N).AND.Y.GT.HYMAX(N))THEN
         UOSCORE(N) = UOSCORE(N) + Z
         IUOSCORE(N) = IUOSCORE(N) + 1
      ELSEIF(X.GT.HXMAX(N).AND.Y.GT.HYMAX(N))THEN
         OOSCORE(N) = OOSCORE(N) + Z
         IOOSCORE(N) = IOOSCORE(N) + 1
      ELSEIF(X.LT.HXMIN(N).AND.Y.GE.HYMIN(N).AND.Y.LE.HYMAX(N))THEN
         UESCORE(N) = UESCORE(N) + Z
         IUESCORE(N) = IUESCORE(N) + 1
      ELSEIF(X.GT.HXMAX(N).AND.Y.GE.HYMIN(N).AND.Y.LE.HYMAX(N))THEN
         OESCORE(N) = OESCORE(N) + Z
         IOESCORE(N) = IOESCORE(N) + 1
      ELSEIF(X.GE.HXMIN(N).AND.X.LE.HXMAX(N).AND.Y.LT.HYMIN(N))THEN
         EUSCORE(N) = EUSCORE(N) + Z
         IEUSCORE(N) = IEUSCORE(N) + 1
      ELSEIF(X.GE.HXMIN(N).AND.X.LE.HXMAX(N).AND.Y.GT.HYMAX(N))THEN
         EOSCORE(N) = EOSCORE(N) + Z
         IEOSCORE(N) = IEOSCORE(N) + 1
      ELSE
         XI=((X-HXMIN(N))/HXDEL(N))+1
         YI=((Y-HYMIN(N))/HYDEL(N))+1
         I=INT(XI)
         J=INT(YI)
         IENT(N)=IENT(N)+1
         IHIS(N,I,J)=IHIS(N,I,J)+1
         HIST(N,I,J)=HIST(N,I,J)+Z
      ENDIF
      END


C      SUBROUTINE MINTEG_3D(NIN,NOUT,IDIR,IPOW)
C If IPOW=1 performs the integral of the distribution contained in histogram
C NIN up to the value specified by the abscissa (if IDIR=1) or from this
C value on (if IDIR=-1). The resulting integral distribution is put into 
C NOUT, which is automatically booked if NOUT.ne.NIN .  Choosing IPOW=2
C the routine will return the square root of the integral of the squares,
C as is required, for example, for the propagation of the mean quadratic error
C of a given distribution. Overscores and underscores are included.
C      IMPLICIT NONE
C
C     ARGUMENTS
C     
C      INTEGER NIN,NOUT,IDIR,IPOW
C
C     GLOBAL
C
C      INCLUDE 'dbook-3d.inc'
c
c     LOCAL
c
C      INTEGER M,I,L
C      CHARACTER*14  C
C      DIMENSION C(2) 
C      DATA C/' INTEG BELOW X',' INTEG ABOVE X'/
C
C     EXTERNAL
C
C      REAL*8 SUMPOW
C
C     START
C
C      IF(BOOK(NIN)(1:1).NE.'Y') RETURN
C      M = NXBIN(NIN)                                           
C      I = (IDIR + 3)/2
C      IF(NOUT.NE.NIN.AND.BOOK(NOUT)(1:1).NE.'Y') THEN
C      	CALL MBOOK_3D(NOUT,TITLE(NIN)//C(I), 
C     &                HXDEL(NIN),HXMIN(NIN),HXMAX(NIN),
C     &                HYDEL(NIN),HYMIN(NIN),HYMAX(NIN))
C      ENDIF
C      IF(IDIR.EQ.1) THEN
C         HIST(NOUT,1) = SUMPOW_3D(HIST(NIN,1,1),UUSCORE(NIN),IPOW)
C         IHIS(NOUT,1) = IHIS(NIN,1) + IUSCORE(NIN)
C         XHIS(NOUT,1) = XHIS(NIN,1) + HDEL(NIN)/2
C         DO L=2,M                      
C            HIST(NOUT,L) = SUMPOW(HIST(NIN,L),HIST(NOUT,L-1),IPOW)
C            IHIS(NOUT,L) = IHIS(NIN,L) + IHIS(NOUT,L-1) 
C            XHIS(NOUT,L) = XHIS(NIN,L) + HDEL(NIN)/2
C         ENDDO
C         OSCORE(NOUT) = SUMPOW(OSCORE(NIN),HIST(NIN,M),IPOW)
C         IOSCORE(NOUT) = IOSCORE(NIN) + IHIS(NIN,M)
C      ELSEIF(IDIR.EQ.-1) THEN
C         HIST(NOUT,M) = SUMPOW(HIST(NIN,M),OSCORE(NIN),IPOW)
C         IHIS(NOUT,M) = IHIS(NIN,M) + IOSCORE(NIN)
C         XHIS(NOUT,M) = XHIS(NIN,M) - HDEL(NIN)/2
C         DO L=M-1,1,-1                        
C            HIST(NOUT,L) = SUMPOW(HIST(NIN,L),HIST(NOUT,L+1),IPOW)
C            IHIS(NOUT,L) = IHIS(NIN,L) + IHIS(NOUT,L+1)
C            XHIS(NOUT,L) = XHIS(NIN,L) - HDEL(NIN)/2
C         ENDDO
C         USCORE(NOUT) = SUMPOW(USCORE(NIN),HIST(NIN,1),IPOW)
C         IUSCORE(NOUT) = IUSCORE(NIN)+IHIS(NIN,1)
C      ELSE                                 
C         CALL MWARN_3D('MINTEG')
C         WRITE(*,*) 'Wrong idir in minteg: OPERATION NOT PERFORMED'
C         RETURN
C      ENDIF
C      END

      REAL*8 FUNCTION SUMPOW_3D(X,Y,IPOW)
C*************************************************************
      IMPLICIT NONE
C
C     ARGUMENTS
C
      REAL*8 X,Y
      INTEGER IPOW
C
C     START
C
      IF(IPOW.EQ.1) THEN
         SUMPOW_3D = X + Y
      ELSEIF(IPOW.EQ.2) THEN
         SUMPOW_3D = DSQRT(X**2+Y**2)
      ELSEIF(IPOW.EQ.0) THEN
         CALL MWARN_3D('SUMPOW_3D')
         WRITE(*,*)'Error: IPOW=0 not allowed in SUMPOW_3D'
      ELSE
         SUMPOW_3D = (X**IPOW+Y**IPOW)**(1./IPOW)
      ENDIF
      END

      SUBROUTINE MOPERA_3D(I,OPER,J,K,X,Y)
C*************************************************************
      IMPLICIT NONE
C
C     ARGUMENTS
C
      INTEGER I,J,K
      CHARACTER OPER*1
      REAL*8 X,Y
C
C     LOCAL
C
      REAL*8 XXX,XSUM,XSUMSQ,XNORM,XAVG,XSQAVG
      INTEGER L,L2
C
C     GLOBAL
C
      INCLUDE 'dbook-3d.inc'
C
C     START
C
      IF( (BOOK(I)(1:1).NE.'Y'.AND.BOOK(I).NE.'NST') .OR.
     &    (BOOK(J)(1:1).NE.'Y'.AND.BOOK(J).NE.'NST') ) RETURN
      IF((NXBIN(I).NE.NXBIN(J).OR.NYBIN(I).NE.NYBIN(J))
     &     .AND.(OPER.EQ.'+'.OR.OPER.EQ.'-'.OR.OPER.EQ.
     &    '*'.OR.OPER.EQ.'/'.OR.OPER.EQ.'M'.OR.OPER.EQ.'A')) THEN
	  CALL MWARN_3D('MOPERA_3D')
          WRITE(*,*) I,J                               
  20      FORMAT(' ****** INCOMPATIBLE OPERATION HIST ',I2,' &',I2,
     &    '*******'/)
          RETURN
      ENDIF
      IF(BOOK(K)(1:1).NE.'Y') 
     &  CALL MBOOK_3D(K,TITLE(I),HXDEL(I),HXMIN(I),HXMAX(I),
     &     HYDEL(I),HYMIN(I),HYMAX(I))
      IF(OPER.EQ.'E') THEN
c If I contains the accumulated weights, J the accumulated squares of the
c weights and IHIS(J,1,1) the number of accumulated entries, 'E' will add
c the average value of I to K and will put in J the quadratic dispersion.
         IF(IHIS(J,1,1).NE.0) THEN
            XXX = 1./IHIS(J,1,1)
         ELSE
            XXX = 0
         ENDIF
         DO L=1,NXBIN(I)
            DO L2=1,NYBIN(I)
               XSUM   = HIST(I,L,L2)
               XSUMSQ = HIST(J,L,L2)
               HIST(K,L,L2)=HIST(K,L,L2) + XXX*HIST(I,L,L2)
               HIST(J,L,L2)=XXX*DSQRT(ABS(XSUMSQ-XSUM**2*XXX))
            ENDDO
         ENDDO
         IENT(K)=IENT(K)+IENT(I)
         XSUM = UUSCORE(I)
         XSUMSQ = UUSCORE(J)
         UUSCORE(K) = UUSCORE(K)+XXX*XSUM
         UUSCORE(J) = XXX*DSQRT(ABS(XSUMSQ-XSUM**2*XXX))
         XSUM = UOSCORE(I)
         XSUMSQ = UOSCORE(J)
         UOSCORE(K) = UOSCORE(K)+XXX*XSUM
         UOSCORE(J) = XXX*DSQRT(ABS(XSUMSQ-XSUM**2*XXX))
         XSUM = OUSCORE(I)
         XSUMSQ = OUSCORE(J)
         OUSCORE(K) = OUSCORE(K)+XXX*XSUM
         OUSCORE(J) = XXX*DSQRT(ABS(XSUMSQ-XSUM**2*XXX))
         XSUM = OOSCORE(I)
         XSUMSQ = OOSCORE(J)
         OOSCORE(K) = OOSCORE(K)+XXX*XSUM
         OOSCORE(J) = XXX*DSQRT(ABS(XSUMSQ-XSUM**2*XXX))
         XSUM = UESCORE(I)
         XSUMSQ = UESCORE(J)
         UESCORE(K) = UESCORE(K)+XXX*XSUM
         UESCORE(J) = XXX*DSQRT(ABS(XSUMSQ-XSUM**2*XXX))
         XSUM = OESCORE(I)
         XSUMSQ = OESCORE(J)
         OESCORE(K) = OESCORE(K)+XXX*XSUM
         OESCORE(J) = XXX*DSQRT(ABS(XSUMSQ-XSUM**2*XXX))
         XSUM = EUSCORE(I)
         XSUMSQ = EUSCORE(J)
         EUSCORE(K) = EUSCORE(K)+XXX*XSUM
         EUSCORE(J) = XXX*DSQRT(ABS(XSUMSQ-XSUM**2*XXX))
         XSUM = EOSCORE(I)
         XSUMSQ = EOSCORE(J)
         EOSCORE(K) = EOSCORE(K)+XXX*XSUM
         EOSCORE(J) = XXX*DSQRT(ABS(XSUMSQ-XSUM**2*XXX))
      ELSEIF(OPER.EQ.'Q') THEN
         DO L=1,NXBIN(I)
            DO L2=1,NYBIN(I)
               HIST(K,L,L2) = SQRT(HIST(J,L,L2)**2+HIST(I,L,L2)**2)
            ENDDO
         ENDDO
         UUSCORE(K) = SQRT(UUSCORE(J)**2+UUSCORE(I)**2)
         OUSCORE(K) = SQRT(OUSCORE(J)**2+OUSCORE(I)**2)
         UOSCORE(K) = SQRT(UOSCORE(J)**2+UOSCORE(I)**2)
         OOSCORE(K) = SQRT(OOSCORE(J)**2+OOSCORE(I)**2)
         UESCORE(K) = SQRT(UESCORE(J)**2+UESCORE(I)**2)
         EUSCORE(K) = SQRT(EUSCORE(J)**2+EUSCORE(I)**2)
         OESCORE(K) = SQRT(OESCORE(J)**2+OESCORE(I)**2)
         EOSCORE(K) = SQRT(EOSCORE(J)**2+EOSCORE(I)**2)
      ELSEIF(OPER.EQ.'A') THEN
         DO L=1,NXBIN(I)
            DO L2=1,NYBIN(I)
               HIST(J,L,L2) = HIST(J,L,L2) + HIST(I,L,L2)
               IHIS(J,L,L2) = IHIS(J,L,L2) + IHIS(I,L,L2)
               HIST(K,L,L2) = HIST(K,L,L2) + HIST(I,L,L2)**2
               IHIS(K,L,L2) = IHIS(K,L,L2) + 1
               HIST(I,L,L2) = 0
               IHIS(I,L,L2) = 0
            ENDDO
         ENDDO
         IENT(J) = IENT(J)+IENT(I)
         IUUSCORE(J) = IUUSCORE(J) + IUUSCORE(I)
         UUSCORE(J) = UUSCORE(J) + UUSCORE(I)
         IUOSCORE(J) = IUOSCORE(J) + IUOSCORE(I)
         UOSCORE(J) = UOSCORE(J) + UOSCORE(I)
         IOUSCORE(J) = IOUSCORE(J) + IOUSCORE(I)
         OUSCORE(J) = OUSCORE(J) + OUSCORE(I)
         IOOSCORE(J) = IOOSCORE(J) + IOOSCORE(I)
         OOSCORE(J) = OOSCORE(J) + OOSCORE(I)
         IUESCORE(J) = IUESCORE(J) + IUESCORE(I)
         UESCORE(J) = UESCORE(J) + UESCORE(I)
         IEUSCORE(J) = IEUSCORE(J) + IEUSCORE(I)
         EUSCORE(J) = EUSCORE(J) + EUSCORE(I)
         IOESCORE(J) = IOESCORE(J) + IOESCORE(I)
         OESCORE(J) = OESCORE(J) + OESCORE(I)
         IEOSCORE(J) = IEOSCORE(J) + IEOSCORE(I)
         EOSCORE(J) = EOSCORE(J) + EOSCORE(I)
         IENT(K) = IENT(K)+1
         IUUSCORE(K) = IUUSCORE(K) + 1
         UUSCORE(K) = UUSCORE(K) + UUSCORE(I)**2
         IUOSCORE(K) = IUOSCORE(K) + 1
         UOSCORE(K) = UOSCORE(K) + UOSCORE(I)**2
         IOUSCORE(K) = IOUSCORE(K) + 1
         OUSCORE(K) = OUSCORE(K) + OUSCORE(I)**2
         IOOSCORE(K) = IOOSCORE(K) + 1
         OOSCORE(K) = OOSCORE(K) + OOSCORE(I)**2
         IUESCORE(K) = IUESCORE(K) + 1
         UESCORE(K) = UESCORE(K) + UESCORE(I)**2
         IOESCORE(K) = IOESCORE(K) + 1
         OESCORE(K) = OESCORE(K) + OESCORE(I)**2
         IEUSCORE(K) = IEUSCORE(K) + 1
         EUSCORE(K) = EUSCORE(K) + EUSCORE(I)**2
         IEOSCORE(K) = IEOSCORE(K) + 1
         EOSCORE(K) = EOSCORE(K) + EOSCORE(I)**2
         IENT(I) = 0
         IUUSCORE(I) = 0
         IUOSCORE(I) = 0
         IOUSCORE(I) = 0
         IOOSCORE(I) = 0
         IUESCORE(I) = 0
         IOESCORE(I) = 0
         IEUSCORE(I) = 0
         IEOSCORE(I) = 0
         UUSCORE(I) = 0
         UOSCORE(I) = 0
         OUSCORE(I) = 0
         OOSCORE(I) = 0
         UESCORE(I) = 0
         OESCORE(I) = 0
         EUSCORE(I) = 0
         EOSCORE(I) = 0
      ELSEIF(OPER.EQ.'X') THEN
         DO L=1,NXBIN(I)
            DO L2=1, NYBIN(I)
               HIST(I,L,L2) = 0
               IHIS(I,L,L2) = 0
            ENDDO
         ENDDO
         IENT(I) = 0
         IUUSCORE(I) = 0
         IUOSCORE(I) = 0
         IOUSCORE(I) = 0
         IOOSCORE(I) = 0
         IUESCORE(I) = 0
         IOESCORE(I) = 0
         IEUSCORE(I) = 0
         IEOSCORE(I) = 0
         UUSCORE(I) = 0
         UOSCORE(I) = 0
         OUSCORE(I) = 0
         OOSCORE(I) = 0
         UESCORE(I) = 0
         OESCORE(I) = 0
         EUSCORE(I) = 0
         EOSCORE(I) = 0
      ELSE
        DO L=1,NXBIN(I)
        DO L2=1,NYBIN(I)
      	IF(OPER.EQ.'+') THEN
       	  HIST(K,L,L2)=X*HIST(I,L,L2) + Y*HIST(J,L,L2)
      	ELSEIF(OPER.EQ.'-') THEN
      	  HIST(K,L,L2)=X*HIST(I,L,L2) - Y*HIST(J,L,L2)
      	ELSEIF(OPER.EQ.'*') THEN
      	  HIST(K,L,L2)=X*HIST(I,L,L2) * Y*HIST(J,L,L2)
      	ELSEIF(OPER.EQ.'/') THEN
          IF(Y.EQ.0..OR.HIST(J,L,L2).EQ.0.) THEN
            HIST(K,L,L2)=0.
          ELSE
            HIST(K,L,L2)=X*HIST(I,L,L2) / (Y*HIST(J,L,L2))
          ENDIF
       	ELSEIF(OPER.EQ.'F') THEN
      	  HIST(K,L,L2)=X*HIST(I,L,L2)
      	ELSEIF(OPER.EQ.'R') THEN
          IF(HIST(I,L,L2).GT.0.) THEN
            HIST(K,L,L2)=X*DSQRT(HIST(I,L,L2))
          ELSE                           
            HIST(K,L,L2)=0.
          ENDIF
      	ELSEIF(OPER.EQ.'S') THEN
          HIST(K,L,L2)=X*HIST(I,L,L2)**2
      	ELSEIF(OPER.EQ.'L') THEN  
          IF(HIST(I,L,L2).EQ.0..OR.Y.EQ.0.) THEN
             HIST(K,L,L2)=0.
           ELSE
             HIST(K,L,L2)=X*LOG10(Y*HIST(I,L,L2))
           ENDIF
      	ELSEIF(OPER.EQ.'M') THEN
           IF(I.NE.J) XNORM=HIST(I,L,L2)
           IF(I.EQ.J) XNORM=DFLOAT(IHIS(J,L,L2))
           IF(XNORM.NE.0.) THEN
             XAVG=HIST(J,L,L2)/XNORM
             HIST(K,L,L2)=
     &       DSQRT(ABS(-XAVG**2+HIST(K,L,L2)/XNORM)/DFLOAT(IHIS(I,L,L2)))
             HIST(J,L,L2)=XAVG 
           ELSE 
             HIST(K,L,L2)=0.
             HIST(J,L,L2)=0.
           ENDIF
      	ELSEIF(OPER.EQ.'V') THEN                 
           XAVG=HIST(I,L,L2)*X
           XSQAVG=HIST(J,L,L2)*X
           XNORM=DFLOAT(IHIS(I,L,L2))*X
           IF(XNORM.NE.0.) THEN
              HIST(K,L,L2)=DSQRT(ABS(XSQAVG-XAVG**2)/XNORM)
              HIST(I,L,L2)=XAVG
           ELSE  
              HIST(K,L,L2)=0.
           ENDIF 
      	ELSE 
	 CALL MWARN_3D('MOPERA_3D')
         WRITE(*,*) OPER
   5     FORMAT(' ****** OPERATION ="',A1,'" UNKNOWN ********'/)
         RETURN
        ENDIF
        ENDDO
        ENDDO
      	IF(OPER.EQ.'+') THEN
       	  UUSCORE(K)=X*UUSCORE(I) + Y*UUSCORE(J)
          UOSCORE(K)=X*UOSCORE(I) + Y*UOSCORE(J)
          OUSCORE(K)=X*OUSCORE(I) + Y*OUSCORE(J)
          OOSCORE(K)=X*OOSCORE(I) + Y*OOSCORE(J)
          UESCORE(K)=X*UESCORE(I) + Y*UESCORE(J)
          OESCORE(K)=X*OESCORE(I) + Y*OESCORE(J)
          EUSCORE(K)=X*EUSCORE(I) + Y*EUSCORE(J)
          EOSCORE(K)=X*EOSCORE(I) + Y*EOSCORE(J)
      	ELSEIF(OPER.EQ.'-') THEN     
      	  UUSCORE(K)=X*UUSCORE(I) - Y*UUSCORE(J)
          UOSCORE(K)=X*UOSCORE(I) - Y*UOSCORE(J)
          OUSCORE(K)=X*OUSCORE(I) - Y*OUSCORE(J)
          OOSCORE(K)=X*OOSCORE(I) - Y*OOSCORE(J)
          UESCORE(K)=X*UESCORE(I) - Y*UESCORE(J)
          OESCORE(K)=X*OESCORE(I) - Y*OESCORE(J)
          EUSCORE(K)=X*EUSCORE(I) - Y*EUSCORE(J)
          EOSCORE(K)=X*EOSCORE(I) - Y*EOSCORE(J)
      	ELSEIF(OPER.EQ.'*') THEN     
      	  UUSCORE(K)=X*UUSCORE(I) * Y*UUSCORE(J)
          UOSCORE(K)=X*UOSCORE(I) * Y*UOSCORE(J)
          OUSCORE(K)=X*OUSCORE(I) * Y*OUSCORE(J)
          OOSCORE(K)=X*OOSCORE(I) * Y*OOSCORE(J)
          UESCORE(K)=X*UESCORE(I) * Y*UESCORE(J)
          OESCORE(K)=X*OESCORE(I) * Y*OESCORE(J)
          EUSCORE(K)=X*EUSCORE(I) * Y*EUSCORE(J)
      	  EOSCORE(K)=X*EOSCORE(I) * Y*EOSCORE(J)
      	ELSEIF(OPER.EQ.'/') THEN     
          IF(Y.EQ.0..OR.UUSCORE(J).EQ.0.) THEN
            UUSCORE(K)=0.
          ELSE
            UUSCORE(K)=X*UUSCORE(I) / (Y*UUSCORE(J))
          ENDIF
          IF(Y.EQ.0..OR.UOSCORE(J).EQ.0.) THEN
            UOSCORE(K)=0.
          ELSE
            UOSCORE(K)=X*UOSCORE(I) / (Y*UOSCORE(J))
          ENDIF
          IF(Y.EQ.0..OR.OUSCORE(J).EQ.0.) THEN
            OUSCORE(K)=0.
          ELSE
            OUSCORE(K)=X*OUSCORE(I) / (Y*OUSCORE(J))
          ENDIF
          IF(Y.EQ.0..OR.OOSCORE(J).EQ.0.) THEN
            OOSCORE(K)=0.
          ELSE
            OOSCORE(K)=X*OOSCORE(I) / (Y*OOSCORE(J))
          ENDIF
          IF(Y.EQ.0..OR.UESCORE(J).EQ.0.) THEN
            UESCORE(K)=0.
          ELSE
            UESCORE(K)=X*UESCORE(I) / (Y*UESCORE(J))
          ENDIF
          IF(Y.EQ.0..OR.OESCORE(J).EQ.0.) THEN
            OESCORE(K)=0.
          ELSE
            OESCORE(K)=X*OESCORE(I) / (Y*OESCORE(J))
          ENDIF
          IF(Y.EQ.0..OR.EUSCORE(J).EQ.0.) THEN
            EUSCORE(K)=0.
          ELSE
            EUSCORE(K)=X*EUSCORE(I) / (Y*EUSCORE(J))
          ENDIF
          IF(Y.EQ.0..OR.EOSCORE(J).EQ.0.) THEN
            EOSCORE(K)=0.
          ELSE
            EOSCORE(K)=X*EOSCORE(I) / (Y*EOSCORE(J))
          ENDIF
       	ELSEIF(OPER.EQ.'F') THEN
      	  UUSCORE(K)=X*UUSCORE(I)
      	  UOSCORE(K)=X*UOSCORE(I)
          OUSCORE(K)=X*OUSCORE(I)
          OOSCORE(K)=X*OOSCORE(I)
          UESCORE(K)=X*UESCORE(I)
          OESCORE(K)=X*OESCORE(I)
          EUSCORE(K)=X*EUSCORE(I)
          EOSCORE(K)=X*EOSCORE(I)
      	ELSEIF(OPER.EQ.'R') THEN
          IF(UUSCORE(I).GT.0.) THEN
            UUSCORE(K)=X*DSQRT(UUSCORE(I))
          ELSE                           
            UUSCORE(K)=0.
          ENDIF     
          IF(UOSCORE(I).GT.0.) THEN
            UOSCORE(K)=X*DSQRT(UOSCORE(I))
          ELSE                           
            UOSCORE(K)=0.
          ENDIF
          IF(OUSCORE(I).GT.0.) THEN
            OUSCORE(K)=X*DSQRT(OUSCORE(I))
          ELSE
            OUSCORE(K)=0.
          ENDIF
          IF(OOSCORE(I).GT.0.) THEN
            OOSCORE(K)=X*DSQRT(OOSCORE(I))
          ELSE
            OOSCORE(K)=0.
          ENDIF
          IF(UESCORE(I).GT.0.) THEN
            UESCORE(K)=X*DSQRT(UESCORE(I))
          ELSE
            UESCORE(K)=0.
          ENDIF
          IF(OESCORE(I).GT.0.) THEN
            OESCORE(K)=X*DSQRT(OESCORE(I))
          ELSE
            OESCORE(K)=0.
          ENDIF
          IF(EOSCORE(I).GT.0.) THEN
            EOSCORE(K)=X*DSQRT(EOSCORE(I))
          ELSE
            EOSCORE(K)=0.
          ENDIF
          IF(EUSCORE(I).GT.0.) THEN
            EUSCORE(K)=X*DSQRT(EUSCORE(I))
          ELSE
            EUSCORE(K)=0.
          ENDIF
      	ELSEIF(OPER.EQ.'S') THEN
          UUSCORE(K)=X*UUSCORE(I)**2
          UOSCORE(K)=X*UOSCORE(I)**2
          OUSCORE(K)=X*OUSCORE(I)**2
          OOSCORE(K)=X*OOSCORE(I)**2
          UESCORE(K)=X*UESCORE(I)**2
          OESCORE(K)=X*OESCORE(I)**2
          EUSCORE(K)=X*EUSCORE(I)**2
          EOSCORE(K)=X*EOSCORE(I)**2
      	ELSEIF(OPER.EQ.'L') THEN  
          IF(UUSCORE(I).EQ.0..OR.Y.EQ.0.) THEN
             UUSCORE(K)=0.
           ELSE
             UUSCORE(K)=X*LOG10(Y*UUSCORE(I))
           ENDIF                         
          IF(UOSCORE(I).EQ.0..OR.Y.EQ.0.) THEN
             UOSCORE(K)=0.
           ELSE
             UOSCORE(K)=X*LOG10(Y*UOSCORE(I))
           ENDIF
          IF(OUSCORE(I).EQ.0..OR.Y.EQ.0.) THEN
             OUSCORE(K)=0.
           ELSE
             OUSCORE(K)=X*LOG10(Y*OUSCORE(I))
           ENDIF
          IF(OOSCORE(I).EQ.0..OR.Y.EQ.0.) THEN
             OOSCORE(K)=0.
           ELSE
             OOSCORE(K)=X*LOG10(Y*OOSCORE(I))
           ENDIF
          IF(UESCORE(I).EQ.0..OR.Y.EQ.0.) THEN
             UESCORE(K)=0.
           ELSE
             UESCORE(K)=X*LOG10(Y*UESCORE(I))
           ENDIF
          IF(OESCORE(I).EQ.0..OR.Y.EQ.0.) THEN
             OESCORE(K)=0.
           ELSE
             OESCORE(K)=X*LOG10(Y*OESCORE(I))
           ENDIF
          IF(EUSCORE(I).EQ.0..OR.Y.EQ.0.) THEN
             EUSCORE(K)=0.
           ELSE
             EUSCORE(K)=X*LOG10(Y*EUSCORE(I))
           ENDIF
          IF(EOSCORE(I).EQ.0..OR.Y.EQ.0.) THEN
             EOSCORE(K)=0.
           ELSE
             EOSCORE(K)=X*LOG10(Y*EOSCORE(I))
           ENDIF                         
        ENDIF
      ENDIF
      RETURN
      END

      SUBROUTINE MZERO_3D(N)
C*************************************************************
      IMPLICIT NONE
C
C     ARGUMENTS
C
      INTEGER N 
C
C     LOCAL
C
      INTEGER I,J
C     
C     GLOBAL
C
      INCLUDE 'dbook-3d.inc'
C
C     START
C
      BOOK(N)='RES'
      IENT(N)=0
      IUUSCORE(N)=0
      IUOSCORE(N)=0
      IOUSCORE(N)=0
      IOOSCORE(N)=0
      IUESCORE(N)=0
      IOESCORE(N)=0
      IEUSCORE(N)=0
      IEOSCORE(N)=0
      HAVG(N)=0.
      HINT(N)=0.
      DO I=1,NXBIN(N)
         DO J=1,NYBIN(N)
            HIST(N,I,J)=0.
         ENDDO
      ENDDO
      END

      SUBROUTINE MRESET_3D(N)
C*************************************************************
      IMPLICIT NONE
C
C     ARGUMENTS
C
      INTEGER N
C
C     GLOBAL
C
      INCLUDE 'dbook-3d.inc'
C
C     START
C
      BOOK(N)='RES'
      END

      SUBROUTINE PUTTAG_3D(J,NAME)
C*************************************************************
c Per marcare un istogramma
      IMPLICIT NONE
C
C     ARGUMENTS
C
      CHARACTER*(*) NAME
      INTEGER J
C
C     GLOBAL
C
      INCLUDE 'dbook-3d.inc'
C
C     LOCAL
C
      CHARACTER*(*) TAG
C
C     START
C
      BOOK(J) = NAME
      RETURN
      ENTRY GETTAG_3D(J,TAG)
      TAG = BOOK(J)
      END

      SUBROUTINE MFINAL_3D(N)
C*************************************************************
      IMPLICIT NONE
C
C     ARGUMENTS
C
      INTEGER N
C
C     LOCAL
C
      INTEGER J,J2,IF
      REAL*8 AVG,XIN,SIG,X
C
C     GLOBAL
C
      INCLUDE 'dbook-3d.inc'
C
C     START
C
      IF(BOOK(N)(1:1).NE.'Y') RETURN
      AVG=0
      XIN=0                                  
      SIG=0
      IF=0
      DO J=1,NXBIN(N)
         DO J2=1,NYBIN(N)
            X=HIST(N,J,J2)
            AVG=AVG+X*XHIS(N,J)*YHIS(N,J2)
            XIN=XIN+X
            IF(X.NE.0) IF=1
         ENDDO
      ENDDO             
      IF(XIN.EQ.0) GO TO 10
      AVG = AVG/XIN
      DO J=1,NXBIN(N)
         DO J2=1,NYBIN(N)
            SIG=HIST(N,J,J2)*(XHIS(N,J)*YHIS(N,J2)-AVG)**2+SIG
         ENDDO
      ENDDO
      SIG=DSQRT(ABS(SIG/XIN))
 10   CONTINUE
      HINT(N) = XIN
      HAVG(N) = AVG
      HSIG(N) = SIG
      END               

      SUBROUTINE MNORM_3D(N,X)
C*************************************************************
      IMPLICIT NONE
C
C     ARGUMENTS
C
      INTEGER N
      REAL*8  X
C
C     GLOBAL
C
      INCLUDE 'dbook-3d.inc'
C
C     LOCAL
C
      INTEGER I,J
      REAL*8 Y
C
C     START
C
      IF(BOOK(N)(:1).NE.'Y')RETURN
      IF(HINT(N).EQ.0.) THEN
	CALL MWARN_3D('MNORM_3D')
	WRITE(*,*)' INTEGRAL HIST ',N,' IS ZERO: CANNOT RENORMALIZE'
	RETURN               
      ELSE
	Y=X/HINT(N)
      ENDIF
      DO I=1,NXBIN(N)
         DO J=1,NYBIN(N)
            HIST(N,I,J)=HIST(N,I,J)*Y
         ENDDO
      ENDDO
      HINT(N)=X   
      UUSCORE(N)=UUSCORE(N)*Y
      UOSCORE(N)=UOSCORE(N)*Y
      OOSCORE(N)=OOSCORE(N)*Y
      OUSCORE(N)=OUSCORE(N)*Y
      UESCORE(N)=UESCORE(N)*Y
      OESCORE(N)=OESCORE(N)*Y
      EUSCORE(N)=EUSCORE(N)*Y
      EOSCORE(N)=EOSCORE(N)*Y
      END                  

      SUBROUTINE MPRINT_3D(N)
C*************************************************************
      IMPLICIT NONE
C
C     ARGUMENTS
C
      INTEGER N
C
C     LOCAL
C
      INTEGER INI,J,J2
C
C     GLOBAL
C
      INCLUDE 'dbook-3d.inc'
C
C     START
C
      DATA INI/0/
      IF(INI.EQ.0) THEN
c     CALL IDATE(IMON,IDAY,IYEAR)
c     CALL TIME(CTIME)
      INI=1
      ENDIF
      IF(BOOK(N)(:1).NE.'Y') RETURN
C      WRITE(98,7) N,IYEAR,IMON,IDAY,CTIME(1:5)
      WRITE(98,*) TITLE(N)
      DO J=1,NXBIN(N)
      DO J2=1,NYBIN(N)
      IF(HIST(N,J,J2).EQ.0.)CYCLE
      WRITE(98,'(3X,F10.4,2X,F10.4,2X,E15.4)')  
     &                            XHIS(N,J),YHIS(N,J2),HIST(N,J,J2)
      ENDDO
      ENDDO
      WRITE(98,15) HAVG(N),HSIG(N),HINT(N)
      WRITE(98,20) IENT(N),IUUSCORE(N),IOOSCORE(N)
C    7 FORMAT(4X,'HIST = ',I3,'   19',I2,'-',I2,'-',I2,1X,A5/)
   10 FORMAT(4X,2E10.3)
   15 FORMAT(/' AVG =',E10.3,4X,' RMS =',E10.3,' INTEGRAL =',E10.3,/)
   20 FORMAT(' ENTRIES=',I5,4X,'UNDERSCORE=',I5,4x,'OVERSCORE=',I5,//)
      END


      SUBROUTINE MTOP4_3D(N,M,BTIT,LTIT,SCALE)
      
C*************************************************************
      IMPLICIT NONE
C
C     ARGUMENTS
C
      INTEGER N,M
      CHARACTER*(*) LTIT,BTIT
      CHARACTER*(*) SCALE
C
C     GLOBAL
C
      INCLUDE 'dbook-3d.inc'
      logical plot_2d3d,plot_2d,plot_3d
      common/plot_top_2d3d/plot_2d3d,plot_2d,plot_3d
C
C     LOCAL
C
      INTEGER INI,J,J2
      real*8 dx,dy
C
C     START
C
      DATA INI/0/
      IF(INI.EQ.0) THEN
C      CALL IDATE(IMON,IDAY,IYEAR)
C      CALL TIME(CTIME)
      INI=1
      ENDIF
      IF(BOOK(N)(:1).NE.'Y') RETURN
cRF
      IF (N.eq.1.OR..NOT.plot_2d3d)THEN
         WRITE(99,'(A)')' SET DEVICE POSTSCRIPT ORIENT=3'
         plot_2d3d=.TRUE.
         plot_2d=.FALSE.
         plot_3d=.TRUE.
       ELSEIF(.NOT.plot_3d)THEN
          plot_3d=.TRUE.
       ENDIF

      WRITE(99,100) TITLE(N),BTIT,LTIT,SCALE,HXMIN(N),HXMAX(N)
     &     ,HYMIN(N),HYMAX(N)
  100 FORMAT( /1x,                               
     &' SET INTENSITY 4'/,1X,
     &' SET WINDOW Y 2.5 TO 9.'/,1X,
     &' SET WINDOW X 2.5 TO 12.'/,1X,
     &' SET FONT DUPLEX '/1X, 
     &' SET SYMBOL 5O SIZE 1.8'/,1X,
     &' TITLE TOP ','"',A20,'"',/1X,
     &' TITLE BOTTOM ','"',A20,'"',/1X,
     &' TITLE LEFT ','"',A20,'"',/1X,
     &' TITLE RIGHT ','"','HELAC-ONIA','"',/1X,
     &' SET SCALE Z ',A,/1X,
     &' (SET TICKS TOP OFF)   '/1x,     
     &' SET LIMITS X ',F10.5,' ',F10.5,/1X,
     &' SET LIMITS Y ',F10.5,' ',F10.5,/1X,
     &' READ MESH')
      dx=HXDEL(N)*0.02d0
      dy=HYDEL(N)*0.02d0
      WRITE(99,*)" FOR Y=",(YHIS(N,J2),J2=1,NYBIN(N))
      DO J=1,NXBIN(N)
         WRITE(99,*)
     &        "    X=",XHIS(N,J)," Z=",(HIST(N,J,J2),J2=1,NYBIN(N))
      ENDDO
      WRITE(99,200)dx,dy
  200 FORMAT( /1x,
     & ' HISTOGRAM HIDE xy dx=',G16.3,' dy=',G16.3,' blue',/1x
     & ' PLOT AXIS HIDE')
      WRITE(99,300)HINT(N),HAVG(N),HSIG(N),IENT(N),IENT(N)
     & ,IUUSCORE(N),IOOSCORE(N)
  300 FORMAT( /1x,
     &' (INFO-BOX'/,1X,                               
     &' BOX 7.3 1.3 SIZE 9.5 0.7'/,1X,
     &' SET WINDOW Y 0. TO 2.'/,1X,
     &' SET TITLE SIZE -1.5'/1X,
     &' SET FONT DUPLEX '/1X,
     &' TITLE 3.8 1.4 "X-sect =',1PE10.3,'(nb)   AVG =',1PE10.3,
     &             '   RMS =',1PE10.3,'"',/1X,
     &' TITLE 3.8 1.1 "Entries =',I8,2x,'Entries =',I8,2X
     &                ,'Undersc =',I6,2X,'Oversc =',I6,'"',/1X,
     &' SET TITLE SIZE -2')                            
      WRITE(99,400)
  400 FORMAT('   NEW PLOT')
      END

      SUBROUTINE MULTITOP4_3D(NH,NE,N,M,BTIT,LTIT,SCA)
C*************************************************************
      IMPLICIT NONE
C
C     ARGUMENTS
C
      INTEGER NH,NE,N,M ! NE stores the error histrogram
      CHARACTER*(*) LTIT,BTIT
      CHARACTER*(*) SCA
C
C     GLOBAL
C
      INCLUDE 'dbook-3d.inc'
C
C     LOCAL
C
      REAL*8 YTIT,XTIT,FEXP,FMAX,FMIN,X,XMX,FMX,FMN
      REAL*8 YU,XU,YL,XL,YTIT0,YD,XTIT0,TITX,TITS,TICS,SRED,XD
      INTEGER IXBIN,IYBIN,IPNS,J,I,NOLD,IFRAME,IFRMAX,IP,NS,MOLD
      INTEGER INI
C
C     START
C
      CHARACTER SCALE*3
      CHARACTER*7  PLOT(4)
      DATA PLOT/'SOLID','DASHES','DOTS','DOTDASH'/
C  PLOT SIZE, CORNERS
      REAL*8 WIDTH,HEIGHT,XCORN,YCORN
      DATA WIDTH,HEIGHT/11.5,8.5/,XCORN,YCORN/1.5,1./
C  PLOT VERSUS TEXT FRACTION                  
      REAL*8 XPFRAC,YPFRAC,XTFRAC,YTFRAC
      DATA XPFRAC,YPFRAC/0.75,0.75/,XTFRAC,YTFRAC/0.25,0.25/
C  DEFAULT SIZES                                           
      REAL*8 TIC0,LAB0,TIT0,LABS
      DATA TIT0,LAB0,TIC0/3d0,3d0,0.06d0/
      DATA INI/0/                                          
      IF(INI.EQ.0) THEN
C      CALL IDATE(IMON,IDAY,IYEAR)
C      CALL TIME(CTIME)
      IFRAME=0        
C      WRITE(99,71) IYEAR,IMON,IDAY,CTIME(1:5)
C   71 FORMAT(4X,' (   19',I2,' -',I2,' -',I2,1X,A5/)
      INI=1         
      ENDIF
      IF(SCA.EQ.'REF') THEN
	IFRAME=0
	RETURN
      ENDIF
      IF(BOOK(NH)(:1).NE.'Y') RETURN
      IFRMAX=N*M         
      IFRAME=IFRAME+1
      IF(IFRAME.GT.IFRMAX.OR.N.NE.NOLD.OR.M.NE.MOLD) THEN
      	IFRAME=1
        WRITE(99,202)   
C        WRITE(99,1) IMON,IDAY,CTIME(1:5)
C  1     FORMAT(' SET FONT DUPLEX',/,'  SET TITLE SIZE 2',/,
C     +      ' TITLE 12.8 9 ANGLE -90 ','" MLM   ',I2,'-',I2,1X,A5,'"')
      ENDIF                                
      IF(IFRAME.EQ.1) THEN
    	I=1
	J=1
      ELSEIF(IFRAME.LE.IFRMAX) THEN
	IF(I.LE.N) I=I+1
        IF(I.GT.N) THEN
		I=1
		J=J+1
	ENDIF
      ENDIF
      IF(N.EQ.NOLD) GO TO 10
      NS=N-1
      XD=WIDTH/DFLOAT(N)
      SRED=DSQRT(DFLOAT(N*M))
      TITS=TIT0/SRED          
      LABS=LAB0/SRED
      TICS=TIC0/SRED
      XTIT0=0.55*XPFRAC*XD
      NOLD=N            
10    IF(M.EQ.MOLD) GO TO 20
      YD=HEIGHT/DFLOAT(M)
      YTIT0=0.06*YD
      MOLD=M        
20    CONTINUE
      XL=(I-1)*XD + XCORN
      YL=(M-J)*YD + YCORN
      XU=XL+XD*XPFRAC
      YU=YL+YD*YPFRAC        
      IP=0
      XMX=0.
      DO IXBIN=1,NXBIN(NH)
      DO IYBIN=1,NYBIN(NH)
         X=HIST(NH,IXBIN,IYBIN)
         IF(X.NE.0.) THEN
            IF(XMX.EQ.0.) THEN
               FMX = X + HIST(NE,IXBIN,IYBIN) ! NE is the standard deviation
               FMN = X - HIST(NE,IXBIN,IYBIN)
           ELSE
               FMX=MAX(FMX,X + HIST(NE,IXBIN,IYBIN))
               FMN=MIN(FMN,X - HIST(NE,IXBIN,IYBIN))
           ENDIF
        ENDIF
        XMX=MAX(XMX,ABS(X)+ HIST(NE,IXBIN,IYBIN))
      ENDDO
      ENDDO
      SCALE=SCA
50    IF(SCALE.EQ.'LIN') THEN
        IF(FMN.GE.0.)   FMIN=0.
        IF(FMN.LT.0.)   FMIN=FMN*1.3
        IF(FMX.GT.0.)   FMAX=FMX*1.3
        IF(FMX.LT.0.)   FMAX=0.
      ELSEIF(SCALE.EQ.'LOG') THEN
cRF
c$$$        IF(FMN.LE.0.) THEN
c$$$                SCALE='LIN'
c$$$                GO TO 50
c$$$        ENDIF
        FMAX=10.**( AINT(LOG10(ABS(FMX))+1000001) - 1000000 )
        FMIN=10.**( AINT(LOG10(ABS(FMN))+1000000) - 1000000 )
      ENDIF                         
      WRITE(99,100) TITS,LABS,TICS,XL,XU,YL,YU
100   FORMAT(2X,'( SET FONT DUPLEX',/,
     *       2X,'SET TITLE SIZE ',F8.4,/,
     *       2X,'SET LABEL SIZE ',F8.4,/,
     *       2X,'SET TICKS TOP OFF SIZE ',F8.4,/,
     *       2X,'SET WINDOW X ',F8.4,' TO ',F8.4,/,
     *       2X,'SET WINDOW Y ',F8.4,' TO ',F8.4)
      XTIT=XL+XTIT0
      YTIT=YU+YTIT0
      WRITE(99,101) XL,YTIT,TITLE(NH)(1:20)
101   FORMAT('  TITLE ',2(F8.4,1X),'"',A,'"')                  
      YTIT=YTIT-2.*YTIT0
      WRITE(99,102) XTIT,YTIT,HINT(NH)
102   FORMAT('  TITLE ',2(F8.4,1X),'" INT=',1PE10.3,'"')                  
      YTIT=YTIT-YTIT0
      WRITE(99,103) XTIT,YTIT,IENT(NH)
103   FORMAT('  TITLE ',2(F8.4,1X),'" ENT=',I9,'"')                  
      YTIT=YTIT-YTIT0                         
      IF(IUUSCORE(NH).NE.0.) THEN
        WRITE(99,104) XTIT,YTIT,UUSCORE(NH)
104     FORMAT('  TITLE ',2(F8.4,1X),'" UFL=',1PE10.3,'"')                  
        YTIT=YTIT-YTIT0                      
      ENDIF
      IF(IOOSCORE(NH).NE.0.) THEN
        WRITE(99,105) XTIT,YTIT,OOSCORE(NH)
105     FORMAT('  TITLE ',2(F8.4,1X),'" OFL=',1PE10.3,'"')                  
        YTIT=YTIT-YTIT0                      
      ENDIF      
      WRITE(99,106) XTIT,YTIT,XU,YTIT,XTIT,YTIT,XTIT,YU
106   FORMAT(2X,'SET ORD X Y ',/,2(F8.4,1X),/,2(F8.4,1X),/,
     *       2X,'JOIN TEXT',/,
     *       2X,2(F8.4,1X),/,2(F8.4,1X),/,
     *       2X,'JOIN TEXT')                                    
      WRITE(99,108) TITS*1.5
108   FORMAT(2X,'SET TITLE SIZE ',F8.4)
      WRITE(99,107) BTIT,XL-0.75*XD*XTFRAC,YL+(YU-YL)/3.,LTIT,SCALE,
     * HXMIN(NH),HXMAX(NH),HYMIN(NH),HYMAX(NH),FMIN,FMAX
107   FORMAT(                                           
     &' TITLE BOTTOM ','"',A,'"',/1X,
     &' TITLE RIGHT ','"','HELAC-ONIA','"',/1X,
     &' TITLE ',f10.5,f10.5,' ANGLE 90 ','"',A,'"',/1X,
     &' SET SCALE Y ',A,/1X,
     &' SET TICKS TOP OFF   '/1x,     
     &' SET LIMITS X ',F10.5,' ',F10.5,/1X,
     &' SET LIMITS Y ',F10.5,' ',F10.5,/1X,
     &' SET LIMITS Z ',1PE10.3,' ',1PE10.3,/1X,
     &' READ MESH')               
C                       
C  END HEADER , FILL TOPDRAWER WITH DATA
C
      ENTRY MTFILL_3D(NH,NE,N,M,BTIT,LTIT,SCA)
      IP=IP+1                             
      IF(IP.GT.4) IP=1
      WRITE(99,110) TITLE(NH),HINT(NH),IENT(NH)
110   FORMAT(' ( ',A,/,' ( INT=',1PE10.3,'  ENTRIES=',I12)
C      DO 200 IXBIN=1,NXBIN(NH)           
C      WRITE(99,'(3X,F10.4,2(2X,E15.4))')  
C     &          XHIS(NH,IBIN),HIST(NH,IBIN),HIST(NE,IBIN)
200   CONTINUE                                           
      WRITE(99,201)  PLOT(IP)
      IF(BOOK(NE).NE.'NO')   WRITE(99,*)  '  PLOT'
201   FORMAT(2X,'HISTOGRAM HIDE xy  ',A)           
202   FORMAT('   NEW PLOT',/,/)
203   RETURN     
      END                     


      SUBROUTINE MWARN_3D(ROUT)
C*************************************************************
      CHARACTER*(*) ROUT
      WRITE(*,*) '***********************************************'
      WRITE(*,*) '***** WARNING CALLED FROM ROUTINE ',ROUT,':'
      END


C*******************************************************************
C     END OF THE HISTOGRAMMING PACKAGE
C*******************************************************************


c--F  Add gnuplot output
      
      SUBROUTINE MGNUPLOT4_3D(N,M,BXTIT,BYTIT,LTIT,SCALE,psfile)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*(*) LTIT,BXTIT,BYTIT,SCALE,psfile
C      double precision xpos,ypos,tmp
      include "dbook-3d.inc"
c--- added these variables to scale plots at intermediate steps
      logical scaleplots                  
      double precision scalefac
      common/scaleplots/scalefac,scaleplots
      logical plot_2d,plot_3d,plot_2d3d
      common/plot_gnu_2d3d/plot_2d3d,plot_2d,plot_3d
      logical plus_x
      double precision x_val,y_val
c--- set scaleplots to be false here
      scaleplots=.false.
      scalefac=1d0

      istring=LEN_TRIM(TITLE(N))

      IF(BOOK(N)(:1).NE.'Y') RETURN
      IF (N.eq.1.OR..NOT.plot_2d3d)THEN
         WRITE(97,102)TRIM(psfile)
 102     FORMAT(/1X,
     &' set terminal postscript col enhanced',/1X,
     &' set output ','"',A,'"',/1X,
     &' set style data points',/1X,
     &' set key off')
         plot_2d3d=.TRUE.
         plot_2d=.FALSE.
         plot_3d=.TRUE.
      ELSE
         WRITE(97,*)' unset label'
         IF(.NOT.plot_3d)plot_3d=.TRUE.
      ENDIF
      WRITE(97,101) TITLE(N)(1:istring), TRIM(BXTIT), TRIM(BYTIT),
     & TRIM(LTIT), HXMIN(N),HXMAX(N),HYMIN(N),HYMAX(N)
 101   FORMAT( /1x,
     &' set view map',/1X,
     &' set title ','"',A,' distribution" font "Helvetica, 20"',/1X,
     &' set xlabel ','"',A,'" font "Helvetica, 20"',/1X,
     &' set ylabel ','"',A,'" font "Helvetica, 20"',/1X,
     &' set zlabel ','"',A,
     &' [nb]" font "Helvetica, 20"',/1X,
C     &' set label ','"','HELAC-ONIA','" font "Helvetica, 15"',
C     &' rotate by 90 at 'F10.5,',',F10.5,/1X,
     &' set xrange [ ',F10.5,':',F10.5,']',/1X,
     &' set yrange [ ',F10.5,':',F10.5,']',/1X,
     &' set nologscale y')
      if(SCALE .eq. 'LOG') then
         write(97,*) ' set logscale z'
      endif
      write(97,*) ' splot "-" with pm3d'

      plus_x=.FALSE.
      DO 1 J=1,NXBIN(N)
 3      CONTINUE
         IF(plus_x)THEN
            x_val=XHIS(N,J)+HXDEL(N)/2d0
         ELSE
            x_val=XHIS(N,J)-HXDEL(N)/2d0
         ENDIF
         IF(J.NE.1.OR.plus_x)WRITE(97,*)'    '
         DO 2 J2=1,NYBIN(N)
c     comment by HSS
c      IF(HIST(N,J,J2).EQ.0.) GO TO 1
      if (scaleplots) then
      WRITE(97,'(4(2X,G13.6))')  
     & XHIS(N,J),YHIS(N,J2),scalefac*HIST(N,J,J2),HIST(M,J,J2)
      else
         y_val=YHIS(N,J2)-HYDEL(N)/2d0
      WRITE(97,'(3X,G13.6,2(2X,G13.6))')  
     &                            x_val,y_val,HIST(N,J,J2)
        y_val=YHIS(N,J2)+HYDEL(N)/2d0
      WRITE(97,'(3X,G13.6,2(2X,G13.6))')
     &                            x_val,y_val,HIST(N,J,J2)
      endif
 2     CONTINUE
      IF(plus_x)THEN
          plus_x=.NOT.plus_x
          GOTO 1
       ELSE
          plus_x=.NOT.plus_x
          GOTO 3
       ENDIF
 1     CONTINUE
      WRITE(97,200)
 200   FORMAT(' e',/1X,
     & ' set yrange [*:*]')
      END

      subroutine MGNUPLOT_3D(n,string1,string2,string3,string4,psfile)
      implicit none
      integer n,n_by4,m_by4
      character*(*) string1,string2,string3,string4,psfile
      n_by4=4*(n-1)+1
      m_by4=n_by4+3
c write the 'n' plots with the 'n+3' error bars
      call mgnuplot4_3D(n_by4,m_by4,string1,string2,string3,string4,
     &  psfile)
      return
      end
c--F  Add root output
      
      SUBROUTINE MROOTPLOT4_3D(N,M,BXTIT,BYTIT,LTIT,rootfile)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      INTEGER::J2
      CHARACTER*(*) LTIT,BXTIT,BYTIT,rootfile
      CHARACTER*7 histoid
      integer nn
      include 'dbook-3d.inc'
c--- added these variables to scale plots at intermediate steps
      logical scaleplots                  
      double precision scalefac
      common/scaleplots/scalefac,scaleplots
      logical plot_2d3d,plot_2d,plot_3d
      common/plot2d3d_root_init/plot_2d3d,plot_2d,plot_3d
c--- set scaleplots to be false here
      scaleplots=.false.
      scalefac=1d0

      istring=LEN_TRIM(TITLE(N))

      nn=N/4+1
      if (nn.lt.10) then
         write(histoid, '(A4,I1,A2)') 'id3d', nn, '  '
      elseif ((nn.ge.10).and.(nn.lt.100)) then
         write(histoid, '(A4,I2,A1)') 'id3d', nn, ' '
      elseif (nn.ge.100) then
         write(histoid, '(A4,I3)') 'id3d', nn
      endif
      idlength=LEN_TRIM(histoid)

      IF(BOOK(N)(:1).NE.'Y') RETURN
      IF(N.EQ.1.OR..NOT.plot_2d3d)THEN
         WRITE(96,130)TRIM(rootfile)
 130     FORMAT (/1X,
     & ' {',/1X,
     & ' Int_t colors[50];',/1X,
     & ' Int_t number = 3;',/1X,
     & ' Double_t red[number] = { 1.00, 0.00, 0.00 };',/1X,
     & ' Double_t green[number] = { 0.00, 1.00, 0.00 };',/1X,
     & ' Double_t blue[number] = { 1.00, 0.00, 1.00 };',/1X,
     & ' Double_t length[number] = { 0.00, 0.50, 1.00 };',/1X,
     & ' Int_t nb=50;',/1X,
     & ' Int_t fi = TColor::CreateGradientColorTable(number,'
     & //'length,red,green,blue,nb);',/1X,
     & ' for (int i=0;i<50;i++) colors[i]=fi+i;',/1X,
     & ' gStyle -> SetPalette(50,colors);',/1X,
     & ' hohisto = new TFile(','"',A,'","recreate");',/1X,
     & ' hohisto -> cd();',/1X,
     & ' histos = new TObjArray(0);')
         plot_2d3d=.TRUE.
         plot_2d=.FALSE.
         plot_3d=.TRUE.
      ELSEIF(.NOT.plot_3d)THEN
         WRITE(96,1300)
 1300    FORMAT (/1X,
     & ' Int_t colors[50];',/1X,
     & ' Int_t number = 3;',/1X,
     & ' Double_t red[number] = { 1.00, 0.00, 0.00 };',/1X,
     & ' Double_t green[number] = { 0.00, 1.00, 0.00 };',/1X,
     & ' Double_t blue[number] = { 1.00, 0.00, 1.00 };',/1X,
     & ' Double_t length[number] = { 0.00, 0.50, 1.00 };',/1X,
     & ' Int_t nb=50;',/1X,
     & ' Int_t fi = TColor::CreateGradientColorTable(number,'
     & //'length,red,green,blue,nb);',/1X,
     & ' for (int i=0;i<50;i++) colors[i]=fi+i;',/1X,
     & ' gStyle -> SetPalette(50,colors);')
         plot_3d=.TRUE.
      ENDIF
      WRITE(96,131) histoid(1:idlength), TITLE(N)(1:istring), NXBIN(N), 
     & HXMIN(N), HXMAX(N), NYBIN(N), HYMIN(N), HYMAX(N)
 131    FORMAT ( /1X,
     & ' hohisto -> cd();', /1X,
     & ' TH2F *hist_3d = new TH2F( "', A, '", "', A, '", ',
     & I0, ', ', F10.5, ', ', F10.5, ', ', I0, ', ', 
     & F10.5, ', ', F10.5, ');')

      WRITE(96, 132) histoid(1:idlength), TRIM(BXTIT),
     & histoid(1:idlength), TRIM(BYTIT),
     & histoid(1:idlength), TITLE(N)(1:istring)
 132    FORMAT ( /1X, 
     & ' ', A, ' -> GetXaxis() -> SetTitle("', A, '");', /1X,
     & ' ', A, ' -> GetYaxis() -> SetTitle("', A, '");', /1X,
     & ' ', A, ' -> GetZaxis() -> SetTitle(" ', A, 
     & ' [nb]");', /1X)
      WRITE (96,*) ' ', histoid(1:idlength), ' -> GetZaxis() -> ',
     & 'SetTitleOffset(1.2);'

      WRITE(96,*) ' ', histoid(1:idlength), ' -> SetStats(false);'

      DO 1 J=1,NXBIN(N)
         DO J2=1,NYBIN(N)
      IF(HIST(N,J,J2).EQ.0.) GO TO 2
      if (scaleplots) then
      WRITE(96,'(4(2X,G13.6))')  
     & XHIS(N,J),YHIS(N,J2),scalefac*HIST(N,J,J2),HIST(M,J,J2)
      else
C         write(96,*) ' ', histoid(1:idlength), ' -> Fill(', 
C     &        XHIS(N,J), ', ', HIST(N,J), ');'
         write(96,*) '  int xybin = ', histoid(1:idlength),'->FindBin(',
     & XHIS(n,j),', ',YHIS(n,j2),');' 
         write(96,*) ' ', histoid(1:idlength), ' -> SetBinContent(', 
     &       ' xybin', ', ', HIST(N,J,J2), ');'
         write(96,*) ' ', histoid(1:idlength), ' -> SetBinError(', 
     &       ' xybin', ', ', HIST(M,J,J2), ');'
      endif
 2    CONTINUE
      ENDDO
 1     CONTINUE

      WRITE (96, *) ' hist_3d -> SetOption("colz"); '
      WRITE (96, *) ' histos -> Add(hist_3d); '
      WRITE (96, *) ''
      WRITE (96, *) ''

      END

      subroutine MROOTPLOT_3D(n,string1,string2,string3,
     &     string4,rootfile)
      implicit none
      integer n,n_by4,m_by4
      character*(*) string1,string2,string3,string4,rootfile
      n_by4=4*(n-1)+1
      m_by4=n_by4+3
c write the 'n' plots with the 'n+3' error bars
      call mrootplot4_3d(n_by4,m_by4,string1,string2,string3,rootfile)
      return
      end

      subroutine mfill_3d(n,xvar,yvar,www)
      implicit none
      integer n,n_by4
      double precision xvar,yvar,www
      n_by4=4*(n-1)+1
      call mfill4_3d(n_by4,xvar,yvar,www)
      return
      end

      subroutine bookup_3d(n,string,xdel,xl,xu,ydel,yl,yu)
      implicit none
      integer n,n_by4
      character*(*) string
      double precision xdel,xl,xu,ydel,yl,yu
      n_by4=4*(n-1)+1
      call bookup4_3d(n_by4,string,xdel,xl,xu,ydel,yl,yu)
      return
      end

      subroutine MTOP_3D(n,string1,string2,string3)
      implicit none
      integer n,n_by4,m_by4
      character*(*) string1,string2,string3
      n_by4=4*(n-1)+1
      m_by4=n_by4+3
c write the 'n' plots with the 'n+3' error bars
      call mtop4_3d(n_by4,m_by4,string1,string2,string3)
      return
      end

      subroutine multitop_3d(n,lr,lh,string1,string2,string3)
      implicit none
      integer n,n_by4,lr,lh,m_by4
      character*(*) string1,string2,string3
      n_by4=4*(n-1)+1
      m_by4=n_by4+3
c write the 'n' plots with the 'n+3' error bars
      if(lr*lh.eq.1)then
         call mtop4_3d(n_by4,m_by4,string1,string2,string3)
      else
         call multitop4_3d(n_by4,m_by4,lr,lh,string1,string2,string3)
      endif
      return
      end



      subroutine bookup4_3d(n,string,xdel,xl,xu,ydel,yl,yu)
      implicit none
      integer n
      character*(*) string
      double precision xdel,xl,xu,ydel,yl,yu
c
c Per ogni istogramma da fare, ne sono richiesti quattro
c In n si accumulano i valori in outfun.
c A ogni iterazione di vegas l'istogramma n viene sommato a n+1,
c e il quadrato del suo valore viene sommato a n+2,
c L'istogramma n viene anche usato alla fine per combinare
c i totali dei vari contributi sig0,sig2, etc. mentre
c l'istogramma n+3 viene usato per combinare gli
c errori dei vari contributi sig0,sig2, etc.
c Non si vuole che n e n+3 vengano salvati o riesumati.
c Cambiando il tag in N,N+3, questo non avviene (si guardi
c in mbook e save/restart e anche mclear in questo programma.
c
      call mbook_3d(n,  string,xdel,xl,xu,ydel,yl,yu)
      call mbook_3d(n+1,'tmp ',xdel,xl,xu,ydel,yl,yu)
      call mbook_3d(n+2,'tmp square',xdel,xl,xu,ydel,yl,yu)
      call mbook_3d(n+3,'error ',xdel,xl,xu,ydel,yl,yu)
      call puttag_3d(n,'YST')
      call puttag_3d(n+3,'NST')
      return
      end

      subroutine accum_3d(inclde)
      implicit real * 8 (a-h,o-z)
      include 'dbook-3d.inc'
      PARAMETER (NMB=NPLOTS)
      character * 3 tag
      logical inclde
c
c     Accumula i valori e i valori al quadrato per l'analisi statistica,
c     e svuota l'istogramma di accumulo.
c
      do j=1,nmb-3
         call gettag_3d(j,tag)
         if(tag.eq.'YST') then
             if (inclde) then
c Sum the results in histos 'j' onto 'j+1' and the squares to 'j+2'
c This also empties histograms 'j'
                call mopera_3d(j,'A',j+1,j+2,dum,dum)
             else
c Only empty the histograms 'j'.
                call mopera_3d(j,'X',j,j,dum,dum)
             endif
         endif
      enddo
      return
      end

      subroutine mclear_3d
      implicit real * 8 (a-h,o-z)
      character * 70 files(20), filn
      data nfil/0/
c
c sum up files
      call sumfil_3d(files,nfil)
      nfil=0
      return
      entry addfil_3d(filn)
c adds file filn to the list of save files.
      nfil = nfil+1
      if(nfil.gt.20) then
         write(*,*) 'mclear_3d: too many files.'
         stop
      endif
      files(nfil) = filn
      end

      subroutine sumfil_3d(files,nfil)
      implicit real * 8 (a-h,o-z)
      include 'dbook-3d.inc'
      PARAMETER (NMB=NPLOTS)
      external dummyinit
      character * 3 tag
      character * 70 files(*)
c
      do j=1,nfil
c comment by HSS
c         call resume(files(j),dummyinit,'NO')
c
c completa l'analisi statistica
c
         do k=1,nmb-3
            call gettag_3d(k,tag)
            if(tag.eq.'YST') then
               call addup_3d(k)
            endif
         enddo
      enddo
      end

      subroutine dummyinit_3d()
      return
      end

      subroutine addup_3d(j)
      implicit none
      integer j
      real*8 dum
c
c accumula j+1 riscalato in j e pone in j+2 la stima dell'errore
c
      call mopera_3d(j+1,'E',j+2,j,dum,dum)
c
c accumula l'errore in quadratura
c
      call mopera_3d(j+2,'Q',j+3,j+3,dum,dum)
      end



