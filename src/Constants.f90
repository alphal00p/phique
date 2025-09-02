MODULE Constants
  USE global
  IMPLICIT NONE
  INTEGER,PARAMETER::uinput=45139
  LOGICAL::NFindQ
CONTAINS

  SUBROUTINE ReadConst
    IMPLICIT NONE
    CHARACTER(len=24)::file
    REAL(KIND(1d0))::pi
    INTEGER::iounit
    LOGICAL::lexist
    ! open user's input file
    INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),EXIST=lexist)
    IF(.NOT.lexist)THEN
       PRINT *,"Warning: the file "//TRIM(Input_File)//" does not exist ! STOP !"
       STOP
    ENDIF
    INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),OPENED=lexist)
    IF(lexist)THEN
       INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),NUMBER=iounit)
       IF(iounit.NE.uinput)THEN
          PRINT *,"WARNING: the "//TRIM(Input_File)//" has been linked with another unit ! Close and reopen !"
          CLOSE(UNIT=iounit)
          OPEN(UNIT=uinput,FILE=TRIM(input_dir)//TRIM(Input_File))
       ENDIF
    ELSE
       OPEN(UNIT=uinput,FILE=TRIM(input_dir)//TRIM(Input_File))
    ENDIF
    
    CALL ReadElem_r('qmass',MQ)
    CALL ReadElem_r('alphaemm1',alphaemm1)
    CALL ReadElem_r('alphasMZ',alphasMZ)
    
  END SUBROUTINE ReadConst

  SUBROUTINE ReadElem_r(keyword,aakey)
    IMPLICIT NONE
    CHARACTER(*),INTENT(IN)::keyword
    REAL(KIND(1d0)),INTENT(INOUT)::aakey
    REAL(KIND(1d0))::rrrgg
    NFindQ=.FALSE.
    rrrgg=readvalue_r(keyword,2)
    IF(.NOT.NFindQ)aakey=rrrgg
  END SUBROUTINE ReadElem_r

  FUNCTION readvalue_r(keyword,flag)
    IMPLICIT NONE
    LOGICAL::lexist
    CHARACTER(*),INTENT(IN)::keyword
    REAL(KIND(1d0))::readvalue_r
    INTEGER,INTENT(IN)::flag
    CHARACTER(100)::ctmp
    CHARACTER(LEN(TRIM(keyword)))::keyw
    INTEGER::LengthC
    INTEGER::iostate
    iostate=0
    keyw=TRIM(keyword)
    LengthC=LEN(keyw)
    INQUIRE(FILE=TRIM(tmp_dir)//"temp.inp",EXIST=lexist)

    IF(.NOT.lexist)THEN
       OPEN(UNIT=45141,FILE=TRIM(tmp_dir)//"temp.inp",STATUS="NEW")
    ELSE
       OPEN(UNIT=45141,FILE=TRIM(tmp_dir)//"temp.inp",STATUS="REPLACE")
    ENDIF
    ! search in user's input file
    IF(flag.EQ.0.OR.flag.EQ.2)THEN
       REWIND(uinput)
       DO WHILE(iostate.NE.-1.AND.iostate.LE.0)
          READ(uinput,'(A)',IOSTAT=iostate) ctmp
          IF((ctmp(1:LengthC)==keyw).AND.&
               (ctmp(LengthC+1:LengthC+1)==' '))THEN
             WRITE(45141,'(A)')ctmp(LengthC+1:100)
             CLOSE(UNIT=45141)
             OPEN(FILE=TRIM(tmp_dir)//'temp.inp',UNIT=45141,STATUS='OLD')
             READ(45141,*)readvalue_r
             CLOSE(UNIT=45141)
             RETURN
          ENDIF
       END DO
    ENDIF

    IF(flag.EQ.1.OR.flag.EQ.0)THEN
       PRINT *,"WARNING: keyword "//keyw//" is invalid ! STOP !"
       STOP
    ELSEIF(flag.EQ.2)THEN
       NFindQ=.TRUE.
       readvalue_r=0
       CLOSE(UNIT=45141)
       RETURN
    ENDIF

    PRINT *,"WARNING: keyword "//keyw//" is invalid ! STOP !"
    STOP
  END FUNCTION readvalue_r

  FUNCTION readvalue_l(keyword,flag)
    IMPLICIT NONE
    LOGICAL::lexist
    CHARACTER(*),INTENT(IN)::keyword
    LOGICAL::readvalue_l
    INTEGER,INTENT(IN)::flag
    CHARACTER(100)::ctmp
    CHARACTER(LEN(TRIM(keyword)))::keyw
    INTEGER::LengthC
    INTEGER::iostate

    iostate=0
    keyw=TRIM(keyword)
    LengthC=LEN(keyw)
    INQUIRE(FILE=TRIM(tmp_dir)//"temp.inp",EXIST=lexist)

    IF(.NOT.lexist)THEN
       OPEN(UNIT=45141,FILE=TRIM(tmp_dir)//"temp.inp",STATUS="NEW")
    ELSE
       OPEN(UNIT=45141,FILE=TRIM(tmp_dir)//"temp.inp",STATUS="REPLACE")
    ENDIF

    ! search in user's input file
    IF(flag.EQ.0.OR.flag.EQ.2)THEN
       REWIND(uinput)
       DO WHILE(iostate.NE.-1.AND.iostate.LE.0)
          READ(uinput,'(A)',IOSTAT=iostate) ctmp
          IF((ctmp(1:LengthC)==keyw).AND.&
               (ctmp(LengthC+1:LengthC+1)==' '))THEN
             WRITE(45141,'(A)')ctmp(LengthC+1:100)
             CLOSE(UNIT=45141)
             OPEN(FILE=TRIM(tmp_dir)//'temp.inp',UNIT=45141,STATUS='OLD')
             READ(45141,*)readvalue_l
             CLOSE(UNIT=45141)
             RETURN
          ENDIF
       END DO
    ENDIF
    
    IF(flag.EQ.1.OR.flag.EQ.0)THEN
       PRINT *,"WARNING: keyword "//keyw//" is invalid ! STOP !"
       STOP
    ELSEIF(flag.EQ.2)THEN
       NFindQ=.TRUE.
       readvalue_l=.FALSE.
       CLOSE(UNIT=45141)
       RETURN
    ENDIF
    
    PRINT *,"WARNING: keyword "//keyw//" is invalid ! STOP !"
    STOP
  END FUNCTION readvalue_l

  SUBROUTINE ReadElem_real(keyword,aakey)
    IMPLICIT NONE
    CHARACTER(*),INTENT(IN)::keyword
    REAL(KIND(1d0)),INTENT(INOUT)::aakey
    INTEGER::flag,iounit
    LOGICAL::lexist

    flag=0
    ! open user's input file
    INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),EXIST=lexist)
    IF(.NOT.lexist)THEN
       PRINT *,"Warning: the file "//TRIM(Input_File)//" does not exist ! STOP !"
       STOP
    ENDIF
    INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),OPENED=lexist)
    IF(lexist)THEN
       INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),NUMBER=iounit)
       IF(iounit.NE.uinput)THEN
          PRINT *,"WARNING: the "//TRIM(Input_File)//" has been linked with another unit ! Close and reopen !"
          CLOSE(UNIT=iounit)
          OPEN(UNIT=uinput,FILE=TRIM(input_dir)//TRIM(Input_File))
       ENDIF
    ELSE
       OPEN(UNIT=uinput,FILE=TRIM(input_dir)//TRIM(Input_File))
    ENDIF

    aakey=readvalue_r(keyword,flag)

  END SUBROUTINE ReadElem_real

  SUBROUTINE ReadElem_integer(keyword,aakey)
    IMPLICIT NONE
    CHARACTER(*),INTENT(IN)::keyword
    INTEGER,INTENT(INOUT)::aakey
    INTEGER::flag,iounit
    LOGICAL::lexist

    flag=0
    INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),EXIST=lexist)
    IF(.NOT.lexist)THEN
       PRINT *,"Warning: the file "//TRIM(Input_File)//" does not exist ! STOP !"
       STOP
    ENDIF
    INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),OPENED=lexist)
    IF(lexist)THEN
       INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),NUMBER=iounit)
       IF(iounit.NE.uinput)THEN
          PRINT *,"WARNING: the "//TRIM(Input_File)//" has been linked with another unit ! Close and reopen !"
          CLOSE(UNIT=iounit)
          OPEN(UNIT=uinput,FILE=TRIM(input_dir)//TRIM(Input_File))
       ENDIF
    ELSE
       OPEN(UNIT=uinput,FILE=TRIM(input_dir)//TRIM(Input_File))
    ENDIF

    aakey=readvalue_r(keyword,flag)
    
  END SUBROUTINE ReadElem_integer

  SUBROUTINE ReadElem_logic(keyword,aakey)
    IMPLICIT NONE
    CHARACTER(*),INTENT(IN)::keyword
    LOGICAL,INTENT(INOUT)::aakey
    INTEGER::flag,iounit
    LOGICAL::lexist

    flag=0
    INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),EXIST=lexist)
    IF(.NOT.lexist)THEN
       PRINT *,"Warning: the file "//TRIM(Input_File)//" does not exist ! STOP !"
       STOP
    ENDIF
    INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),OPENED=lexist)
    IF(lexist)THEN
       INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),NUMBER=iounit)
       IF(iounit.NE.uinput)THEN
          PRINT *,"WARNING: the "//TRIM(Input_File)//" has been linked with another unit ! Close and reopen !"
          CLOSE(UNIT=iounit)
          OPEN(UNIT=uinput,FILE=TRIM(input_dir)//TRIM(Input_File))
       ENDIF
    ELSE
       OPEN(UNIT=uinput,FILE=TRIM(input_dir)//TRIM(Input_File))
    ENDIF

    aakey=readvalue_l(keyword,flag)
    
  END SUBROUTINE ReadElem_logic
END MODULE Constants
