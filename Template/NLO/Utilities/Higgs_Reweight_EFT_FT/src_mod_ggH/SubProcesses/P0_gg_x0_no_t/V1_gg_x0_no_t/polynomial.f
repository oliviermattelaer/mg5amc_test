C     THE SUBROUTINE TO CREATE THE COEFFICIENTS FROM LAST LOOP WF AND 
C     MULTIPLY BY THE BORN

      SUBROUTINE CREATE_LOOP_COEFS(LOOP_WF,RANK,LCUT_SIZE,LOOP_GROUP_NU
     $ MBER,SYMFACT,COLOR_ID,HELCONFIG)
      IMPLICIT NONE
C     
C     CONSTANTS 
C     
      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=4)
      REAL*8 ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(ZERO,ONE))
      COMPLEX*16 CMPLX_ZERO
      PARAMETER (CMPLX_ZERO=(ZERO,ZERO))
      INTEGER MAXLWFSIZE
      PARAMETER (MAXLWFSIZE=4)
      INTEGER LOOPMAXCOEFS
      PARAMETER (LOOPMAXCOEFS=126)
      INTEGER    NCOLORROWS
      PARAMETER (NCOLORROWS=224)
      INTEGER    NLOOPGROUPS
      PARAMETER (NLOOPGROUPS=16)
      INTEGER    NCOMB
      PARAMETER (NCOMB=8)
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=1, NSQUAREDSO=1, NAMPSO=2)
C     
C     ARGUMENTS 
C     
      COMPLEX*16 LOOP_WF(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER RANK, COLOR_ID, SYMFACT, LCUT_SIZE, HELCONFIG, LOOP_GROUP
     $ _NUMBER
C     
C     LOCAL VARIABLES 
C     
      COMPLEX*16 CFTOT
      COMPLEX*16 CONST(NAMPSO)
      INTEGER I,J
C     
C     FUNCTIONS
C     
      INTEGER ML5SOINDEX_FOR_BORN_AMP, ML5SOINDEX_FOR_LOOP_AMP
     $ , ML5SQSOINDEX
C     
C     GLOBAL VARIABLES
C     
      INTEGER CF_D(NCOLORROWS,NBORNAMPS)
      INTEGER CF_N(NCOLORROWS,NBORNAMPS)
      COMMON/CF/CF_D,CF_N

      LOGICAL CHECKPHASE
      LOGICAL HELDOUBLECHECKED
      COMMON/INIT/CHECKPHASE, HELDOUBLECHECKED

      INTEGER HELOFFSET
      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/FILTERS/GOODAMP,GOODHEL,HELOFFSET

      COMPLEX*16 LOOPCOEFS(0:LOOPMAXCOEFS-1,NSQUAREDSO,NLOOPGROUPS)
      COMMON/LCOEFS/LOOPCOEFS

      INTEGER HELPICKED
      COMMON/HELCHOICE/HELPICKED

      COMPLEX*16 AMP(NBORNAMPS)
      COMMON/AMPS/AMP

      DO I=1,NAMPSO
        CONST(I)=CMPLX_ZERO
      ENDDO

      DO I=1,NBORNAMPS
        CFTOT=CMPLX(CF_N(COLOR_ID,I)/(ONE*ABS(CF_D(COLOR_ID,I))),ZERO
     $   ,KIND=8)
        IF(CF_D(COLOR_ID,I).LT.0) CFTOT=CFTOT*IMAG1
        CONST(ML5SOINDEX_FOR_BORN_AMP(I))=CONST(ML5SOINDEX_FOR_BORN_AMP
     $   (I))+CFTOT*CONJG(AMP(I))
      ENDDO

      DO I=1,NAMPSO
        IF (CONST(I).NE.CMPLX_ZERO) THEN
          CONST(I)=CONST(I)/SYMFACT
          IF (.NOT.CHECKPHASE.AND.HELDOUBLECHECKED.AND.HELPICKED.EQ.
     $     -1) THEN
            CONST(I)=CONST(I)*GOODHEL(HELCONFIG)
          ENDIF
          CALL MERGE_WL(LOOP_WF,RANK,LCUT_SIZE,CONST(I),LOOPCOEFS(0
     $     ,ML5SQSOINDEX(I,ML5SOINDEX_FOR_LOOP_AMP(COLOR_ID))
     $     ,LOOP_GROUP_NUMBER))
        ENDIF
      ENDDO

      END


C     THE SUBROUTINE TO CREATE THE COEFFICIENTS FROM LAST LOOP WF AND 
C     MULTIPLY BY THE BORN

      SUBROUTINE MP_CREATE_LOOP_COEFS(LOOP_WF,RANK,LCUT_SIZE,LOOP_GROUP
     $ _NUMBER,SYMFACT,COLOR_ID,HELCONFIG)
      IMPLICIT NONE
C     
C     CONSTANTS 
C     
      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=4)
      REAL*16 ZERO,ONE
      PARAMETER (ZERO=0.0E0_16,ONE=1.0E0_16)
      COMPLEX*32 IMAG1
      PARAMETER (IMAG1=(ZERO,ONE))
      COMPLEX*32 CMPLX_ZERO
      PARAMETER (CMPLX_ZERO=(ZERO,ZERO))
      INTEGER MAXLWFSIZE
      PARAMETER (MAXLWFSIZE=4)
      INTEGER LOOPMAXCOEFS
      PARAMETER (LOOPMAXCOEFS=126)
      INTEGER    NCOLORROWS
      PARAMETER (NCOLORROWS=224)
      INTEGER    NLOOPGROUPS
      PARAMETER (NLOOPGROUPS=16)
      INTEGER    NCOMB
      PARAMETER (NCOMB=8)
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=1, NSQUAREDSO=1, NAMPSO=2)
C     
C     ARGUMENTS 
C     
      COMPLEX*32 LOOP_WF(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER RANK, COLOR_ID, SYMFACT, LCUT_SIZE, HELCONFIG, LOOP_GROUP
     $ _NUMBER
C     
C     LOCAL VARIABLES 
C     
      COMPLEX*32 CFTOT
      COMPLEX*32 CONST(NAMPSO)
      INTEGER I,J
C     
C     FUNCTIONS
C     
      INTEGER ML5SOINDEX_FOR_BORN_AMP, ML5SOINDEX_FOR_LOOP_AMP
     $ , ML5SQSOINDEX
C     
C     GLOBAL VARIABLES
C     
      INTEGER CF_D(NCOLORROWS,NBORNAMPS)
      INTEGER CF_N(NCOLORROWS,NBORNAMPS)
      COMMON/CF/CF_D,CF_N

      LOGICAL CHECKPHASE
      LOGICAL HELDOUBLECHECKED
      COMMON/INIT/CHECKPHASE, HELDOUBLECHECKED

      INTEGER HELOFFSET
      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/FILTERS/GOODAMP,GOODHEL,HELOFFSET

      COMPLEX*32 LOOPCOEFS(0:LOOPMAXCOEFS-1,NSQUAREDSO,NLOOPGROUPS)
      COMMON/MP_LCOEFS/LOOPCOEFS

      INTEGER HELPICKED
      COMMON/HELCHOICE/HELPICKED

      COMPLEX*32 AMP(NBORNAMPS)
      COMMON/MP_AMPS/AMP

      DO I=1,NAMPSO
        CONST(I)=CMPLX_ZERO
      ENDDO

      DO I=1,NBORNAMPS
        CFTOT=CMPLX(CF_N(COLOR_ID,I)/(ONE*ABS(CF_D(COLOR_ID,I))),ZERO
     $   ,KIND=16)
        IF(CF_D(COLOR_ID,I).LT.0) CFTOT=CFTOT*IMAG1
        CONST(ML5SOINDEX_FOR_BORN_AMP(I))=CONST(ML5SOINDEX_FOR_BORN_AMP
     $   (I))+CFTOT*CONJG(AMP(I))
      ENDDO

      DO I=1,NAMPSO
        IF (CONST(I).NE.CMPLX_ZERO) THEN
          CONST(I)=CONST(I)/SYMFACT
          IF (.NOT.CHECKPHASE.AND.HELDOUBLECHECKED.AND.HELPICKED.EQ.
     $     -1) THEN
            CONST(I)=CONST(I)*GOODHEL(HELCONFIG)
          ENDIF
          CALL MP_MERGE_WL(LOOP_WF,RANK,LCUT_SIZE,CONST(I),LOOPCOEFS(0
     $     ,ML5SQSOINDEX(I,ML5SOINDEX_FOR_LOOP_AMP(COLOR_ID))
     $     ,LOOP_GROUP_NUMBER))
        ENDIF
      ENDDO

      END


      SUBROUTINE EVAL_POLY(C,R,Q,OUT)
      INCLUDE 'coef_specs.inc'
      COMPLEX*16 C(0:LOOP_MAXCOEFS-1)
      INTEGER R
      COMPLEX*16 Q(0:3)
      COMPLEX*16 OUT

      OUT=C(0)
      IF (R.GE.1) THEN
        OUT=OUT+C(1)*Q(0)+C(2)*Q(1)+C(3)*Q(2)+C(4)*Q(3)
      ENDIF
      IF (R.GE.2) THEN
        OUT=OUT+C(5)*Q(0)*Q(0)+C(6)*Q(0)*Q(1)+C(7)*Q(1)*Q(1)+C(8)*Q(0)
     $   *Q(2)+C(9)*Q(1)*Q(2)+C(10)*Q(2)*Q(2)+C(11)*Q(0)*Q(3)
     $   +C(12)*Q(1)*Q(3)+C(13)*Q(2)*Q(3)+C(14)*Q(3)*Q(3)
      ENDIF
      IF (R.GE.3) THEN
        OUT=OUT+C(15)*Q(0)*Q(0)*Q(0)+C(16)*Q(0)*Q(0)*Q(1)+C(17)*Q(0)
     $   *Q(1)*Q(1)+C(18)*Q(1)*Q(1)*Q(1)+C(19)*Q(0)*Q(0)*Q(2)
     $   +C(20)*Q(0)*Q(1)*Q(2)+C(21)*Q(1)*Q(1)*Q(2)+C(22)*Q(0)*Q(2)
     $   *Q(2)+C(23)*Q(1)*Q(2)*Q(2)+C(24)*Q(2)*Q(2)*Q(2)+C(25)*Q(0)
     $   *Q(0)*Q(3)+C(26)*Q(0)*Q(1)*Q(3)+C(27)*Q(1)*Q(1)*Q(3)
     $   +C(28)*Q(0)*Q(2)*Q(3)+C(29)*Q(1)*Q(2)*Q(3)+C(30)*Q(2)*Q(2)
     $   *Q(3)+C(31)*Q(0)*Q(3)*Q(3)+C(32)*Q(1)*Q(3)*Q(3)+C(33)*Q(2)
     $   *Q(3)*Q(3)+C(34)*Q(3)*Q(3)*Q(3)
      ENDIF
      IF (R.GE.4) THEN
        OUT=OUT+C(35)*Q(0)*Q(0)*Q(0)*Q(0)+C(36)*Q(0)*Q(0)*Q(0)*Q(1)
     $   +C(37)*Q(0)*Q(0)*Q(1)*Q(1)+C(38)*Q(0)*Q(1)*Q(1)*Q(1)
     $   +C(39)*Q(1)*Q(1)*Q(1)*Q(1)+C(40)*Q(0)*Q(0)*Q(0)*Q(2)
     $   +C(41)*Q(0)*Q(0)*Q(1)*Q(2)+C(42)*Q(0)*Q(1)*Q(1)*Q(2)
     $   +C(43)*Q(1)*Q(1)*Q(1)*Q(2)+C(44)*Q(0)*Q(0)*Q(2)*Q(2)
     $   +C(45)*Q(0)*Q(1)*Q(2)*Q(2)+C(46)*Q(1)*Q(1)*Q(2)*Q(2)
     $   +C(47)*Q(0)*Q(2)*Q(2)*Q(2)+C(48)*Q(1)*Q(2)*Q(2)*Q(2)
     $   +C(49)*Q(2)*Q(2)*Q(2)*Q(2)+C(50)*Q(0)*Q(0)*Q(0)*Q(3)
     $   +C(51)*Q(0)*Q(0)*Q(1)*Q(3)+C(52)*Q(0)*Q(1)*Q(1)*Q(3)
     $   +C(53)*Q(1)*Q(1)*Q(1)*Q(3)+C(54)*Q(0)*Q(0)*Q(2)*Q(3)
     $   +C(55)*Q(0)*Q(1)*Q(2)*Q(3)+C(56)*Q(1)*Q(1)*Q(2)*Q(3)
     $   +C(57)*Q(0)*Q(2)*Q(2)*Q(3)+C(58)*Q(1)*Q(2)*Q(2)*Q(3)
     $   +C(59)*Q(2)*Q(2)*Q(2)*Q(3)+C(60)*Q(0)*Q(0)*Q(3)*Q(3)
     $   +C(61)*Q(0)*Q(1)*Q(3)*Q(3)+C(62)*Q(1)*Q(1)*Q(3)*Q(3)
     $   +C(63)*Q(0)*Q(2)*Q(3)*Q(3)+C(64)*Q(1)*Q(2)*Q(3)*Q(3)
        OUT=OUT+C(65)*Q(2)*Q(2)*Q(3)*Q(3)+C(66)*Q(0)*Q(3)*Q(3)*Q(3)
     $   +C(67)*Q(1)*Q(3)*Q(3)*Q(3)+C(68)*Q(2)*Q(3)*Q(3)*Q(3)
     $   +C(69)*Q(3)*Q(3)*Q(3)*Q(3)
      ENDIF
      IF (R.GE.5) THEN
        OUT=OUT+C(70)*Q(0)*Q(0)*Q(0)*Q(0)*Q(0)+C(71)*Q(0)*Q(0)*Q(0)
     $   *Q(0)*Q(1)+C(72)*Q(0)*Q(0)*Q(0)*Q(1)*Q(1)+C(73)*Q(0)*Q(0)
     $   *Q(1)*Q(1)*Q(1)+C(74)*Q(0)*Q(1)*Q(1)*Q(1)*Q(1)+C(75)*Q(1)
     $   *Q(1)*Q(1)*Q(1)*Q(1)+C(76)*Q(0)*Q(0)*Q(0)*Q(0)*Q(2)+C(77)
     $   *Q(0)*Q(0)*Q(0)*Q(1)*Q(2)+C(78)*Q(0)*Q(0)*Q(1)*Q(1)*Q(2)
     $   +C(79)*Q(0)*Q(1)*Q(1)*Q(1)*Q(2)+C(80)*Q(1)*Q(1)*Q(1)*Q(1)
     $   *Q(2)+C(81)*Q(0)*Q(0)*Q(0)*Q(2)*Q(2)+C(82)*Q(0)*Q(0)*Q(1)
     $   *Q(2)*Q(2)+C(83)*Q(0)*Q(1)*Q(1)*Q(2)*Q(2)+C(84)*Q(1)*Q(1)
     $   *Q(1)*Q(2)*Q(2)+C(85)*Q(0)*Q(0)*Q(2)*Q(2)*Q(2)+C(86)*Q(0)
     $   *Q(1)*Q(2)*Q(2)*Q(2)+C(87)*Q(1)*Q(1)*Q(2)*Q(2)*Q(2)+C(88)
     $   *Q(0)*Q(2)*Q(2)*Q(2)*Q(2)+C(89)*Q(1)*Q(2)*Q(2)*Q(2)*Q(2)
     $   +C(90)*Q(2)*Q(2)*Q(2)*Q(2)*Q(2)+C(91)*Q(0)*Q(0)*Q(0)*Q(0)
     $   *Q(3)+C(92)*Q(0)*Q(0)*Q(0)*Q(1)*Q(3)+C(93)*Q(0)*Q(0)*Q(1)
     $   *Q(1)*Q(3)+C(94)*Q(0)*Q(1)*Q(1)*Q(1)*Q(3)+C(95)*Q(1)*Q(1)
     $   *Q(1)*Q(1)*Q(3)+C(96)*Q(0)*Q(0)*Q(0)*Q(2)*Q(3)+C(97)*Q(0)
     $   *Q(0)*Q(1)*Q(2)*Q(3)+C(98)*Q(0)*Q(1)*Q(1)*Q(2)*Q(3)+C(99)
     $   *Q(1)*Q(1)*Q(1)*Q(2)*Q(3)
        OUT=OUT+C(100)*Q(0)*Q(0)*Q(2)*Q(2)*Q(3)+C(101)*Q(0)*Q(1)*Q(2)
     $   *Q(2)*Q(3)+C(102)*Q(1)*Q(1)*Q(2)*Q(2)*Q(3)+C(103)*Q(0)*Q(2)
     $   *Q(2)*Q(2)*Q(3)+C(104)*Q(1)*Q(2)*Q(2)*Q(2)*Q(3)+C(105)*Q(2)
     $   *Q(2)*Q(2)*Q(2)*Q(3)+C(106)*Q(0)*Q(0)*Q(0)*Q(3)*Q(3)
     $   +C(107)*Q(0)*Q(0)*Q(1)*Q(3)*Q(3)+C(108)*Q(0)*Q(1)*Q(1)*Q(3)
     $   *Q(3)+C(109)*Q(1)*Q(1)*Q(1)*Q(3)*Q(3)+C(110)*Q(0)*Q(0)*Q(2)
     $   *Q(3)*Q(3)+C(111)*Q(0)*Q(1)*Q(2)*Q(3)*Q(3)+C(112)*Q(1)*Q(1)
     $   *Q(2)*Q(3)*Q(3)+C(113)*Q(0)*Q(2)*Q(2)*Q(3)*Q(3)+C(114)*Q(1)
     $   *Q(2)*Q(2)*Q(3)*Q(3)+C(115)*Q(2)*Q(2)*Q(2)*Q(3)*Q(3)
     $   +C(116)*Q(0)*Q(0)*Q(3)*Q(3)*Q(3)+C(117)*Q(0)*Q(1)*Q(3)*Q(3)
     $   *Q(3)+C(118)*Q(1)*Q(1)*Q(3)*Q(3)*Q(3)+C(119)*Q(0)*Q(2)*Q(3)
     $   *Q(3)*Q(3)+C(120)*Q(1)*Q(2)*Q(3)*Q(3)*Q(3)+C(121)*Q(2)*Q(2)
     $   *Q(3)*Q(3)*Q(3)+C(122)*Q(0)*Q(3)*Q(3)*Q(3)*Q(3)+C(123)*Q(1)
     $   *Q(3)*Q(3)*Q(3)*Q(3)+C(124)*Q(2)*Q(3)*Q(3)*Q(3)*Q(3)
     $   +C(125)*Q(3)*Q(3)*Q(3)*Q(3)*Q(3)
      ENDIF
      END

      SUBROUTINE MP_EVAL_POLY(C,R,Q,OUT)
      INCLUDE 'coef_specs.inc'
      COMPLEX*32 C(0:LOOP_MAXCOEFS-1)
      INTEGER R
      COMPLEX*32 Q(0:3)
      COMPLEX*32 OUT

      OUT=C(0)
      IF (R.GE.1) THEN
        OUT=OUT+C(1)*Q(0)+C(2)*Q(1)+C(3)*Q(2)+C(4)*Q(3)
      ENDIF
      IF (R.GE.2) THEN
        OUT=OUT+C(5)*Q(0)*Q(0)+C(6)*Q(0)*Q(1)+C(7)*Q(1)*Q(1)+C(8)*Q(0)
     $   *Q(2)+C(9)*Q(1)*Q(2)+C(10)*Q(2)*Q(2)+C(11)*Q(0)*Q(3)
     $   +C(12)*Q(1)*Q(3)+C(13)*Q(2)*Q(3)+C(14)*Q(3)*Q(3)
      ENDIF
      IF (R.GE.3) THEN
        OUT=OUT+C(15)*Q(0)*Q(0)*Q(0)+C(16)*Q(0)*Q(0)*Q(1)+C(17)*Q(0)
     $   *Q(1)*Q(1)+C(18)*Q(1)*Q(1)*Q(1)+C(19)*Q(0)*Q(0)*Q(2)
     $   +C(20)*Q(0)*Q(1)*Q(2)+C(21)*Q(1)*Q(1)*Q(2)+C(22)*Q(0)*Q(2)
     $   *Q(2)+C(23)*Q(1)*Q(2)*Q(2)+C(24)*Q(2)*Q(2)*Q(2)+C(25)*Q(0)
     $   *Q(0)*Q(3)+C(26)*Q(0)*Q(1)*Q(3)+C(27)*Q(1)*Q(1)*Q(3)
     $   +C(28)*Q(0)*Q(2)*Q(3)+C(29)*Q(1)*Q(2)*Q(3)+C(30)*Q(2)*Q(2)
     $   *Q(3)+C(31)*Q(0)*Q(3)*Q(3)+C(32)*Q(1)*Q(3)*Q(3)+C(33)*Q(2)
     $   *Q(3)*Q(3)+C(34)*Q(3)*Q(3)*Q(3)
      ENDIF
      IF (R.GE.4) THEN
        OUT=OUT+C(35)*Q(0)*Q(0)*Q(0)*Q(0)+C(36)*Q(0)*Q(0)*Q(0)*Q(1)
     $   +C(37)*Q(0)*Q(0)*Q(1)*Q(1)+C(38)*Q(0)*Q(1)*Q(1)*Q(1)
     $   +C(39)*Q(1)*Q(1)*Q(1)*Q(1)+C(40)*Q(0)*Q(0)*Q(0)*Q(2)
     $   +C(41)*Q(0)*Q(0)*Q(1)*Q(2)+C(42)*Q(0)*Q(1)*Q(1)*Q(2)
     $   +C(43)*Q(1)*Q(1)*Q(1)*Q(2)+C(44)*Q(0)*Q(0)*Q(2)*Q(2)
     $   +C(45)*Q(0)*Q(1)*Q(2)*Q(2)+C(46)*Q(1)*Q(1)*Q(2)*Q(2)
     $   +C(47)*Q(0)*Q(2)*Q(2)*Q(2)+C(48)*Q(1)*Q(2)*Q(2)*Q(2)
     $   +C(49)*Q(2)*Q(2)*Q(2)*Q(2)+C(50)*Q(0)*Q(0)*Q(0)*Q(3)
     $   +C(51)*Q(0)*Q(0)*Q(1)*Q(3)+C(52)*Q(0)*Q(1)*Q(1)*Q(3)
     $   +C(53)*Q(1)*Q(1)*Q(1)*Q(3)+C(54)*Q(0)*Q(0)*Q(2)*Q(3)
     $   +C(55)*Q(0)*Q(1)*Q(2)*Q(3)+C(56)*Q(1)*Q(1)*Q(2)*Q(3)
     $   +C(57)*Q(0)*Q(2)*Q(2)*Q(3)+C(58)*Q(1)*Q(2)*Q(2)*Q(3)
     $   +C(59)*Q(2)*Q(2)*Q(2)*Q(3)+C(60)*Q(0)*Q(0)*Q(3)*Q(3)
     $   +C(61)*Q(0)*Q(1)*Q(3)*Q(3)+C(62)*Q(1)*Q(1)*Q(3)*Q(3)
     $   +C(63)*Q(0)*Q(2)*Q(3)*Q(3)+C(64)*Q(1)*Q(2)*Q(3)*Q(3)
        OUT=OUT+C(65)*Q(2)*Q(2)*Q(3)*Q(3)+C(66)*Q(0)*Q(3)*Q(3)*Q(3)
     $   +C(67)*Q(1)*Q(3)*Q(3)*Q(3)+C(68)*Q(2)*Q(3)*Q(3)*Q(3)
     $   +C(69)*Q(3)*Q(3)*Q(3)*Q(3)
      ENDIF
      IF (R.GE.5) THEN
        OUT=OUT+C(70)*Q(0)*Q(0)*Q(0)*Q(0)*Q(0)+C(71)*Q(0)*Q(0)*Q(0)
     $   *Q(0)*Q(1)+C(72)*Q(0)*Q(0)*Q(0)*Q(1)*Q(1)+C(73)*Q(0)*Q(0)
     $   *Q(1)*Q(1)*Q(1)+C(74)*Q(0)*Q(1)*Q(1)*Q(1)*Q(1)+C(75)*Q(1)
     $   *Q(1)*Q(1)*Q(1)*Q(1)+C(76)*Q(0)*Q(0)*Q(0)*Q(0)*Q(2)+C(77)
     $   *Q(0)*Q(0)*Q(0)*Q(1)*Q(2)+C(78)*Q(0)*Q(0)*Q(1)*Q(1)*Q(2)
     $   +C(79)*Q(0)*Q(1)*Q(1)*Q(1)*Q(2)+C(80)*Q(1)*Q(1)*Q(1)*Q(1)
     $   *Q(2)+C(81)*Q(0)*Q(0)*Q(0)*Q(2)*Q(2)+C(82)*Q(0)*Q(0)*Q(1)
     $   *Q(2)*Q(2)+C(83)*Q(0)*Q(1)*Q(1)*Q(2)*Q(2)+C(84)*Q(1)*Q(1)
     $   *Q(1)*Q(2)*Q(2)+C(85)*Q(0)*Q(0)*Q(2)*Q(2)*Q(2)+C(86)*Q(0)
     $   *Q(1)*Q(2)*Q(2)*Q(2)+C(87)*Q(1)*Q(1)*Q(2)*Q(2)*Q(2)+C(88)
     $   *Q(0)*Q(2)*Q(2)*Q(2)*Q(2)+C(89)*Q(1)*Q(2)*Q(2)*Q(2)*Q(2)
     $   +C(90)*Q(2)*Q(2)*Q(2)*Q(2)*Q(2)+C(91)*Q(0)*Q(0)*Q(0)*Q(0)
     $   *Q(3)+C(92)*Q(0)*Q(0)*Q(0)*Q(1)*Q(3)+C(93)*Q(0)*Q(0)*Q(1)
     $   *Q(1)*Q(3)+C(94)*Q(0)*Q(1)*Q(1)*Q(1)*Q(3)+C(95)*Q(1)*Q(1)
     $   *Q(1)*Q(1)*Q(3)+C(96)*Q(0)*Q(0)*Q(0)*Q(2)*Q(3)+C(97)*Q(0)
     $   *Q(0)*Q(1)*Q(2)*Q(3)+C(98)*Q(0)*Q(1)*Q(1)*Q(2)*Q(3)+C(99)
     $   *Q(1)*Q(1)*Q(1)*Q(2)*Q(3)
        OUT=OUT+C(100)*Q(0)*Q(0)*Q(2)*Q(2)*Q(3)+C(101)*Q(0)*Q(1)*Q(2)
     $   *Q(2)*Q(3)+C(102)*Q(1)*Q(1)*Q(2)*Q(2)*Q(3)+C(103)*Q(0)*Q(2)
     $   *Q(2)*Q(2)*Q(3)+C(104)*Q(1)*Q(2)*Q(2)*Q(2)*Q(3)+C(105)*Q(2)
     $   *Q(2)*Q(2)*Q(2)*Q(3)+C(106)*Q(0)*Q(0)*Q(0)*Q(3)*Q(3)
     $   +C(107)*Q(0)*Q(0)*Q(1)*Q(3)*Q(3)+C(108)*Q(0)*Q(1)*Q(1)*Q(3)
     $   *Q(3)+C(109)*Q(1)*Q(1)*Q(1)*Q(3)*Q(3)+C(110)*Q(0)*Q(0)*Q(2)
     $   *Q(3)*Q(3)+C(111)*Q(0)*Q(1)*Q(2)*Q(3)*Q(3)+C(112)*Q(1)*Q(1)
     $   *Q(2)*Q(3)*Q(3)+C(113)*Q(0)*Q(2)*Q(2)*Q(3)*Q(3)+C(114)*Q(1)
     $   *Q(2)*Q(2)*Q(3)*Q(3)+C(115)*Q(2)*Q(2)*Q(2)*Q(3)*Q(3)
     $   +C(116)*Q(0)*Q(0)*Q(3)*Q(3)*Q(3)+C(117)*Q(0)*Q(1)*Q(3)*Q(3)
     $   *Q(3)+C(118)*Q(1)*Q(1)*Q(3)*Q(3)*Q(3)+C(119)*Q(0)*Q(2)*Q(3)
     $   *Q(3)*Q(3)+C(120)*Q(1)*Q(2)*Q(3)*Q(3)*Q(3)+C(121)*Q(2)*Q(2)
     $   *Q(3)*Q(3)*Q(3)+C(122)*Q(0)*Q(3)*Q(3)*Q(3)*Q(3)+C(123)*Q(1)
     $   *Q(3)*Q(3)*Q(3)*Q(3)+C(124)*Q(2)*Q(3)*Q(3)*Q(3)*Q(3)
     $   +C(125)*Q(3)*Q(3)*Q(3)*Q(3)*Q(3)
      ENDIF
      END

      SUBROUTINE ADD_COEFS(A,RA,B,RB)
      INCLUDE 'coef_specs.inc'
      INTEGER I
      COMPLEX*16 A(0:LOOP_MAXCOEFS-1),B(0:LOOP_MAXCOEFS-1)
      INTEGER RA,RB

      INTEGER NCOEF_R(0:5)
      DATA NCOEF_R/1,5,15,35,70,126/

      DO I=0,NCOEF_R(RB)-1
        A(I)=A(I)+B(I)
      ENDDO
      END

      SUBROUTINE MP_ADD_COEFS(A,RA,B,RB)
      INCLUDE 'coef_specs.inc'
      INTEGER I
      COMPLEX*32 A(0:LOOP_MAXCOEFS-1),B(0:LOOP_MAXCOEFS-1)
      INTEGER RA,RB

      INTEGER NCOEF_R(0:5)
      DATA NCOEF_R/1,5,15,35,70,126/

      DO I=0,NCOEF_R(RB)-1
        A(I)=A(I)+B(I)
      ENDDO
      END

      SUBROUTINE MERGE_WL(WL,R,LCUT_SIZE,CONST,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J
      COMPLEX*16 WL(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER R,LCUT_SIZE
      COMPLEX*16 CONST
      COMPLEX*16 OUT(0:LOOP_MAXCOEFS-1)

      INTEGER NCOEF_R(0:5)
      DATA NCOEF_R/1,5,15,35,70,126/

      DO I=1,LCUT_SIZE
        DO J=0,NCOEF_R(R)-1
          OUT(J)=OUT(J)+WL(I,J,I)*CONST
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_MERGE_WL(WL,R,LCUT_SIZE,CONST,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J
      COMPLEX*32 WL(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER R,LCUT_SIZE
      COMPLEX*32 CONST
      COMPLEX*32 OUT(0:LOOP_MAXCOEFS-1)

      INTEGER NCOEF_R(0:5)
      DATA NCOEF_R/1,5,15,35,70,126/

      DO I=1,LCUT_SIZE
        DO J=0,NCOEF_R(R)-1
          OUT(J)=OUT(J)+WL(I,J,I)*CONST
        ENDDO
      ENDDO
      END

      SUBROUTINE UPDATE_WL_0_1(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*16 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,4
            OUT(J,K,I)=(0.0D0,0.0D0)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_UPDATE_WL_0_1(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*32 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,4
            OUT(J,K,I)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE UPDATE_WL_1_2(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*16 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,34
            OUT(J,K,I)=(0.0D0,0.0D0)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,0,I)*B(J,5,K)+A(K,1,I)*B(J,1,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,0,I)*B(J,6,K)+A(K,1,I)*B(J,2,K)
     $       +A(K,2,I)*B(J,1,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,0,I)*B(J,7,K)+A(K,2,I)*B(J,2,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,0,I)*B(J,8,K)+A(K,1,I)*B(J,3,K)
     $       +A(K,3,I)*B(J,1,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,0,I)*B(J,9,K)+A(K,2,I)*B(J,3,K)
     $       +A(K,3,I)*B(J,2,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,0,I)*B(J,10,K)+A(K,3,I)
     $       *B(J,3,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,0,I)*B(J,11,K)+A(K,1,I)
     $       *B(J,4,K)+A(K,4,I)*B(J,1,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,0,I)*B(J,12,K)+A(K,2,I)
     $       *B(J,4,K)+A(K,4,I)*B(J,2,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,0,I)*B(J,13,K)+A(K,3,I)
     $       *B(J,4,K)+A(K,4,I)*B(J,3,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,0,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,4,K)
            OUT(J,15,I)=OUT(J,15,I)+A(K,1,I)*B(J,5,K)
            OUT(J,16,I)=OUT(J,16,I)+A(K,1,I)*B(J,6,K)+A(K,2,I)*B(J,5,K)
            OUT(J,17,I)=OUT(J,17,I)+A(K,1,I)*B(J,7,K)+A(K,2,I)*B(J,6,K)
            OUT(J,18,I)=OUT(J,18,I)+A(K,2,I)*B(J,7,K)
            OUT(J,19,I)=OUT(J,19,I)+A(K,1,I)*B(J,8,K)+A(K,3,I)*B(J,5,K)
            OUT(J,20,I)=OUT(J,20,I)+A(K,1,I)*B(J,9,K)+A(K,2,I)
     $       *B(J,8,K)+A(K,3,I)*B(J,6,K)
            OUT(J,21,I)=OUT(J,21,I)+A(K,2,I)*B(J,9,K)+A(K,3,I)*B(J,7,K)
            OUT(J,22,I)=OUT(J,22,I)+A(K,1,I)*B(J,10,K)+A(K,3,I)
     $       *B(J,8,K)
            OUT(J,23,I)=OUT(J,23,I)+A(K,2,I)*B(J,10,K)+A(K,3,I)
     $       *B(J,9,K)
            OUT(J,24,I)=OUT(J,24,I)+A(K,3,I)*B(J,10,K)
            OUT(J,25,I)=OUT(J,25,I)+A(K,1,I)*B(J,11,K)+A(K,4,I)
     $       *B(J,5,K)
            OUT(J,26,I)=OUT(J,26,I)+A(K,1,I)*B(J,12,K)+A(K,2,I)
     $       *B(J,11,K)+A(K,4,I)*B(J,6,K)
            OUT(J,27,I)=OUT(J,27,I)+A(K,2,I)*B(J,12,K)+A(K,4,I)
     $       *B(J,7,K)
            OUT(J,28,I)=OUT(J,28,I)+A(K,1,I)*B(J,13,K)+A(K,3,I)
     $       *B(J,11,K)+A(K,4,I)*B(J,8,K)
            OUT(J,29,I)=OUT(J,29,I)+A(K,2,I)*B(J,13,K)+A(K,3,I)
     $       *B(J,12,K)+A(K,4,I)*B(J,9,K)
            OUT(J,30,I)=OUT(J,30,I)+A(K,3,I)*B(J,13,K)+A(K,4,I)
     $       *B(J,10,K)
            OUT(J,31,I)=OUT(J,31,I)+A(K,1,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,11,K)
            OUT(J,32,I)=OUT(J,32,I)+A(K,2,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,12,K)
            OUT(J,33,I)=OUT(J,33,I)+A(K,3,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,13,K)
            OUT(J,34,I)=OUT(J,34,I)+A(K,4,I)*B(J,14,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_UPDATE_WL_1_2(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*32 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,34
            OUT(J,K,I)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,0,I)*B(J,5,K)+A(K,1,I)*B(J,1,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,0,I)*B(J,6,K)+A(K,1,I)*B(J,2,K)
     $       +A(K,2,I)*B(J,1,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,0,I)*B(J,7,K)+A(K,2,I)*B(J,2,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,0,I)*B(J,8,K)+A(K,1,I)*B(J,3,K)
     $       +A(K,3,I)*B(J,1,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,0,I)*B(J,9,K)+A(K,2,I)*B(J,3,K)
     $       +A(K,3,I)*B(J,2,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,0,I)*B(J,10,K)+A(K,3,I)
     $       *B(J,3,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,0,I)*B(J,11,K)+A(K,1,I)
     $       *B(J,4,K)+A(K,4,I)*B(J,1,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,0,I)*B(J,12,K)+A(K,2,I)
     $       *B(J,4,K)+A(K,4,I)*B(J,2,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,0,I)*B(J,13,K)+A(K,3,I)
     $       *B(J,4,K)+A(K,4,I)*B(J,3,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,0,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,4,K)
            OUT(J,15,I)=OUT(J,15,I)+A(K,1,I)*B(J,5,K)
            OUT(J,16,I)=OUT(J,16,I)+A(K,1,I)*B(J,6,K)+A(K,2,I)*B(J,5,K)
            OUT(J,17,I)=OUT(J,17,I)+A(K,1,I)*B(J,7,K)+A(K,2,I)*B(J,6,K)
            OUT(J,18,I)=OUT(J,18,I)+A(K,2,I)*B(J,7,K)
            OUT(J,19,I)=OUT(J,19,I)+A(K,1,I)*B(J,8,K)+A(K,3,I)*B(J,5,K)
            OUT(J,20,I)=OUT(J,20,I)+A(K,1,I)*B(J,9,K)+A(K,2,I)
     $       *B(J,8,K)+A(K,3,I)*B(J,6,K)
            OUT(J,21,I)=OUT(J,21,I)+A(K,2,I)*B(J,9,K)+A(K,3,I)*B(J,7,K)
            OUT(J,22,I)=OUT(J,22,I)+A(K,1,I)*B(J,10,K)+A(K,3,I)
     $       *B(J,8,K)
            OUT(J,23,I)=OUT(J,23,I)+A(K,2,I)*B(J,10,K)+A(K,3,I)
     $       *B(J,9,K)
            OUT(J,24,I)=OUT(J,24,I)+A(K,3,I)*B(J,10,K)
            OUT(J,25,I)=OUT(J,25,I)+A(K,1,I)*B(J,11,K)+A(K,4,I)
     $       *B(J,5,K)
            OUT(J,26,I)=OUT(J,26,I)+A(K,1,I)*B(J,12,K)+A(K,2,I)
     $       *B(J,11,K)+A(K,4,I)*B(J,6,K)
            OUT(J,27,I)=OUT(J,27,I)+A(K,2,I)*B(J,12,K)+A(K,4,I)
     $       *B(J,7,K)
            OUT(J,28,I)=OUT(J,28,I)+A(K,1,I)*B(J,13,K)+A(K,3,I)
     $       *B(J,11,K)+A(K,4,I)*B(J,8,K)
            OUT(J,29,I)=OUT(J,29,I)+A(K,2,I)*B(J,13,K)+A(K,3,I)
     $       *B(J,12,K)+A(K,4,I)*B(J,9,K)
            OUT(J,30,I)=OUT(J,30,I)+A(K,3,I)*B(J,13,K)+A(K,4,I)
     $       *B(J,10,K)
            OUT(J,31,I)=OUT(J,31,I)+A(K,1,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,11,K)
            OUT(J,32,I)=OUT(J,32,I)+A(K,2,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,12,K)
            OUT(J,33,I)=OUT(J,33,I)+A(K,3,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,13,K)
            OUT(J,34,I)=OUT(J,34,I)+A(K,4,I)*B(J,14,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE UPDATE_WL_3_2(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*16 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,125
            OUT(J,K,I)=(0.0D0,0.0D0)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,0,I)*B(J,5,K)+A(K,1,I)*B(J,1,K)
     $       +A(K,5,I)*B(J,0,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,0,I)*B(J,6,K)+A(K,1,I)*B(J,2,K)
     $       +A(K,2,I)*B(J,1,K)+A(K,6,I)*B(J,0,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,0,I)*B(J,7,K)+A(K,2,I)*B(J,2,K)
     $       +A(K,7,I)*B(J,0,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,0,I)*B(J,8,K)+A(K,1,I)*B(J,3,K)
     $       +A(K,3,I)*B(J,1,K)+A(K,8,I)*B(J,0,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,0,I)*B(J,9,K)+A(K,2,I)*B(J,3,K)
     $       +A(K,3,I)*B(J,2,K)+A(K,9,I)*B(J,0,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,0,I)*B(J,10,K)+A(K,3,I)
     $       *B(J,3,K)+A(K,10,I)*B(J,0,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,0,I)*B(J,11,K)+A(K,1,I)
     $       *B(J,4,K)+A(K,4,I)*B(J,1,K)+A(K,11,I)*B(J,0,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,0,I)*B(J,12,K)+A(K,2,I)
     $       *B(J,4,K)+A(K,4,I)*B(J,2,K)+A(K,12,I)*B(J,0,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,0,I)*B(J,13,K)+A(K,3,I)
     $       *B(J,4,K)+A(K,4,I)*B(J,3,K)+A(K,13,I)*B(J,0,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,0,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,4,K)+A(K,14,I)*B(J,0,K)
            OUT(J,15,I)=OUT(J,15,I)+A(K,1,I)*B(J,5,K)+A(K,5,I)
     $       *B(J,1,K)+A(K,15,I)*B(J,0,K)
            OUT(J,16,I)=OUT(J,16,I)+A(K,1,I)*B(J,6,K)+A(K,2,I)
     $       *B(J,5,K)+A(K,5,I)*B(J,2,K)+A(K,6,I)*B(J,1,K)+A(K,16,I)
     $       *B(J,0,K)
            OUT(J,17,I)=OUT(J,17,I)+A(K,1,I)*B(J,7,K)+A(K,2,I)
     $       *B(J,6,K)+A(K,6,I)*B(J,2,K)+A(K,7,I)*B(J,1,K)+A(K,17,I)
     $       *B(J,0,K)
            OUT(J,18,I)=OUT(J,18,I)+A(K,2,I)*B(J,7,K)+A(K,7,I)
     $       *B(J,2,K)+A(K,18,I)*B(J,0,K)
            OUT(J,19,I)=OUT(J,19,I)+A(K,1,I)*B(J,8,K)+A(K,3,I)
     $       *B(J,5,K)+A(K,5,I)*B(J,3,K)+A(K,8,I)*B(J,1,K)+A(K,19,I)
     $       *B(J,0,K)
            OUT(J,20,I)=OUT(J,20,I)+A(K,1,I)*B(J,9,K)+A(K,2,I)
     $       *B(J,8,K)+A(K,3,I)*B(J,6,K)+A(K,6,I)*B(J,3,K)+A(K,8,I)
     $       *B(J,2,K)+A(K,9,I)*B(J,1,K)+A(K,20,I)*B(J,0,K)
            OUT(J,21,I)=OUT(J,21,I)+A(K,2,I)*B(J,9,K)+A(K,3,I)
     $       *B(J,7,K)+A(K,7,I)*B(J,3,K)+A(K,9,I)*B(J,2,K)+A(K,21,I)
     $       *B(J,0,K)
            OUT(J,22,I)=OUT(J,22,I)+A(K,1,I)*B(J,10,K)+A(K,3,I)
     $       *B(J,8,K)+A(K,8,I)*B(J,3,K)+A(K,10,I)*B(J,1,K)+A(K,22,I)
     $       *B(J,0,K)
            OUT(J,23,I)=OUT(J,23,I)+A(K,2,I)*B(J,10,K)+A(K,3,I)
     $       *B(J,9,K)+A(K,9,I)*B(J,3,K)+A(K,10,I)*B(J,2,K)+A(K,23,I)
     $       *B(J,0,K)
            OUT(J,24,I)=OUT(J,24,I)+A(K,3,I)*B(J,10,K)+A(K,10,I)
     $       *B(J,3,K)+A(K,24,I)*B(J,0,K)
            OUT(J,25,I)=OUT(J,25,I)+A(K,1,I)*B(J,11,K)+A(K,4,I)
     $       *B(J,5,K)+A(K,5,I)*B(J,4,K)+A(K,11,I)*B(J,1,K)+A(K,25,I)
     $       *B(J,0,K)
            OUT(J,26,I)=OUT(J,26,I)+A(K,1,I)*B(J,12,K)+A(K,2,I)
     $       *B(J,11,K)+A(K,4,I)*B(J,6,K)+A(K,6,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,2,K)+A(K,12,I)*B(J,1,K)+A(K,26,I)*B(J,0,K)
            OUT(J,27,I)=OUT(J,27,I)+A(K,2,I)*B(J,12,K)+A(K,4,I)
     $       *B(J,7,K)+A(K,7,I)*B(J,4,K)+A(K,12,I)*B(J,2,K)+A(K,27,I)
     $       *B(J,0,K)
            OUT(J,28,I)=OUT(J,28,I)+A(K,1,I)*B(J,13,K)+A(K,3,I)
     $       *B(J,11,K)+A(K,4,I)*B(J,8,K)+A(K,8,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,1,K)+A(K,28,I)*B(J,0,K)
            OUT(J,29,I)=OUT(J,29,I)+A(K,2,I)*B(J,13,K)+A(K,3,I)
     $       *B(J,12,K)+A(K,4,I)*B(J,9,K)+A(K,9,I)*B(J,4,K)+A(K,12,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,2,K)+A(K,29,I)*B(J,0,K)
            OUT(J,30,I)=OUT(J,30,I)+A(K,3,I)*B(J,13,K)+A(K,4,I)
     $       *B(J,10,K)+A(K,10,I)*B(J,4,K)+A(K,13,I)*B(J,3,K)
     $       +A(K,30,I)*B(J,0,K)
            OUT(J,31,I)=OUT(J,31,I)+A(K,1,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,11,K)+A(K,11,I)*B(J,4,K)+A(K,14,I)*B(J,1,K)
     $       +A(K,31,I)*B(J,0,K)
            OUT(J,32,I)=OUT(J,32,I)+A(K,2,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,12,K)+A(K,12,I)*B(J,4,K)+A(K,14,I)*B(J,2,K)
     $       +A(K,32,I)*B(J,0,K)
            OUT(J,33,I)=OUT(J,33,I)+A(K,3,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,13,K)+A(K,13,I)*B(J,4,K)+A(K,14,I)*B(J,3,K)
     $       +A(K,33,I)*B(J,0,K)
            OUT(J,34,I)=OUT(J,34,I)+A(K,4,I)*B(J,14,K)+A(K,14,I)
     $       *B(J,4,K)+A(K,34,I)*B(J,0,K)
            OUT(J,35,I)=OUT(J,35,I)+A(K,5,I)*B(J,5,K)+A(K,15,I)
     $       *B(J,1,K)
            OUT(J,36,I)=OUT(J,36,I)+A(K,5,I)*B(J,6,K)+A(K,6,I)
     $       *B(J,5,K)+A(K,15,I)*B(J,2,K)+A(K,16,I)*B(J,1,K)
            OUT(J,37,I)=OUT(J,37,I)+A(K,5,I)*B(J,7,K)+A(K,6,I)
     $       *B(J,6,K)+A(K,7,I)*B(J,5,K)+A(K,16,I)*B(J,2,K)+A(K,17,I)
     $       *B(J,1,K)
            OUT(J,38,I)=OUT(J,38,I)+A(K,6,I)*B(J,7,K)+A(K,7,I)
     $       *B(J,6,K)+A(K,17,I)*B(J,2,K)+A(K,18,I)*B(J,1,K)
            OUT(J,39,I)=OUT(J,39,I)+A(K,7,I)*B(J,7,K)+A(K,18,I)
     $       *B(J,2,K)
            OUT(J,40,I)=OUT(J,40,I)+A(K,5,I)*B(J,8,K)+A(K,8,I)
     $       *B(J,5,K)+A(K,15,I)*B(J,3,K)+A(K,19,I)*B(J,1,K)
            OUT(J,41,I)=OUT(J,41,I)+A(K,5,I)*B(J,9,K)+A(K,6,I)
     $       *B(J,8,K)+A(K,8,I)*B(J,6,K)+A(K,9,I)*B(J,5,K)+A(K,16,I)
     $       *B(J,3,K)+A(K,19,I)*B(J,2,K)+A(K,20,I)*B(J,1,K)
            OUT(J,42,I)=OUT(J,42,I)+A(K,6,I)*B(J,9,K)+A(K,7,I)
     $       *B(J,8,K)+A(K,8,I)*B(J,7,K)+A(K,9,I)*B(J,6,K)+A(K,17,I)
     $       *B(J,3,K)+A(K,20,I)*B(J,2,K)+A(K,21,I)*B(J,1,K)
            OUT(J,43,I)=OUT(J,43,I)+A(K,7,I)*B(J,9,K)+A(K,9,I)
     $       *B(J,7,K)+A(K,18,I)*B(J,3,K)+A(K,21,I)*B(J,2,K)
            OUT(J,44,I)=OUT(J,44,I)+A(K,5,I)*B(J,10,K)+A(K,8,I)
     $       *B(J,8,K)+A(K,10,I)*B(J,5,K)+A(K,19,I)*B(J,3,K)+A(K,22,I)
     $       *B(J,1,K)
            OUT(J,45,I)=OUT(J,45,I)+A(K,6,I)*B(J,10,K)+A(K,8,I)
     $       *B(J,9,K)+A(K,9,I)*B(J,8,K)+A(K,10,I)*B(J,6,K)+A(K,20,I)
     $       *B(J,3,K)+A(K,22,I)*B(J,2,K)+A(K,23,I)*B(J,1,K)
            OUT(J,46,I)=OUT(J,46,I)+A(K,7,I)*B(J,10,K)+A(K,9,I)
     $       *B(J,9,K)+A(K,10,I)*B(J,7,K)+A(K,21,I)*B(J,3,K)+A(K,23,I)
     $       *B(J,2,K)
            OUT(J,47,I)=OUT(J,47,I)+A(K,8,I)*B(J,10,K)+A(K,10,I)
     $       *B(J,8,K)+A(K,22,I)*B(J,3,K)+A(K,24,I)*B(J,1,K)
            OUT(J,48,I)=OUT(J,48,I)+A(K,9,I)*B(J,10,K)+A(K,10,I)
     $       *B(J,9,K)+A(K,23,I)*B(J,3,K)+A(K,24,I)*B(J,2,K)
            OUT(J,49,I)=OUT(J,49,I)+A(K,10,I)*B(J,10,K)+A(K,24,I)
     $       *B(J,3,K)
            OUT(J,50,I)=OUT(J,50,I)+A(K,5,I)*B(J,11,K)+A(K,11,I)
     $       *B(J,5,K)+A(K,15,I)*B(J,4,K)+A(K,25,I)*B(J,1,K)
            OUT(J,51,I)=OUT(J,51,I)+A(K,5,I)*B(J,12,K)+A(K,6,I)
     $       *B(J,11,K)+A(K,11,I)*B(J,6,K)+A(K,12,I)*B(J,5,K)
     $       +A(K,16,I)*B(J,4,K)+A(K,25,I)*B(J,2,K)+A(K,26,I)*B(J,1,K)
            OUT(J,52,I)=OUT(J,52,I)+A(K,6,I)*B(J,12,K)+A(K,7,I)
     $       *B(J,11,K)+A(K,11,I)*B(J,7,K)+A(K,12,I)*B(J,6,K)
     $       +A(K,17,I)*B(J,4,K)+A(K,26,I)*B(J,2,K)+A(K,27,I)*B(J,1,K)
            OUT(J,53,I)=OUT(J,53,I)+A(K,7,I)*B(J,12,K)+A(K,12,I)
     $       *B(J,7,K)+A(K,18,I)*B(J,4,K)+A(K,27,I)*B(J,2,K)
            OUT(J,54,I)=OUT(J,54,I)+A(K,5,I)*B(J,13,K)+A(K,8,I)
     $       *B(J,11,K)+A(K,11,I)*B(J,8,K)+A(K,13,I)*B(J,5,K)
     $       +A(K,19,I)*B(J,4,K)+A(K,25,I)*B(J,3,K)+A(K,28,I)*B(J,1,K)
            OUT(J,55,I)=OUT(J,55,I)+A(K,6,I)*B(J,13,K)+A(K,8,I)
     $       *B(J,12,K)+A(K,9,I)*B(J,11,K)+A(K,11,I)*B(J,9,K)
     $       +A(K,12,I)*B(J,8,K)+A(K,13,I)*B(J,6,K)+A(K,20,I)*B(J,4,K)
     $       +A(K,26,I)*B(J,3,K)+A(K,28,I)*B(J,2,K)+A(K,29,I)*B(J,1,K)
            OUT(J,56,I)=OUT(J,56,I)+A(K,7,I)*B(J,13,K)+A(K,9,I)
     $       *B(J,12,K)+A(K,12,I)*B(J,9,K)+A(K,13,I)*B(J,7,K)
     $       +A(K,21,I)*B(J,4,K)+A(K,27,I)*B(J,3,K)+A(K,29,I)*B(J,2,K)
            OUT(J,57,I)=OUT(J,57,I)+A(K,8,I)*B(J,13,K)+A(K,10,I)
     $       *B(J,11,K)+A(K,11,I)*B(J,10,K)+A(K,13,I)*B(J,8,K)
     $       +A(K,22,I)*B(J,4,K)+A(K,28,I)*B(J,3,K)+A(K,30,I)*B(J,1,K)
            OUT(J,58,I)=OUT(J,58,I)+A(K,9,I)*B(J,13,K)+A(K,10,I)
     $       *B(J,12,K)+A(K,12,I)*B(J,10,K)+A(K,13,I)*B(J,9,K)
     $       +A(K,23,I)*B(J,4,K)+A(K,29,I)*B(J,3,K)+A(K,30,I)*B(J,2,K)
            OUT(J,59,I)=OUT(J,59,I)+A(K,10,I)*B(J,13,K)+A(K,13,I)
     $       *B(J,10,K)+A(K,24,I)*B(J,4,K)+A(K,30,I)*B(J,3,K)
            OUT(J,60,I)=OUT(J,60,I)+A(K,5,I)*B(J,14,K)+A(K,11,I)
     $       *B(J,11,K)+A(K,14,I)*B(J,5,K)+A(K,25,I)*B(J,4,K)
     $       +A(K,31,I)*B(J,1,K)
            OUT(J,61,I)=OUT(J,61,I)+A(K,6,I)*B(J,14,K)+A(K,11,I)
     $       *B(J,12,K)+A(K,12,I)*B(J,11,K)+A(K,14,I)*B(J,6,K)
     $       +A(K,26,I)*B(J,4,K)+A(K,31,I)*B(J,2,K)+A(K,32,I)*B(J,1,K)
            OUT(J,62,I)=OUT(J,62,I)+A(K,7,I)*B(J,14,K)+A(K,12,I)
     $       *B(J,12,K)+A(K,14,I)*B(J,7,K)+A(K,27,I)*B(J,4,K)
     $       +A(K,32,I)*B(J,2,K)
            OUT(J,63,I)=OUT(J,63,I)+A(K,8,I)*B(J,14,K)+A(K,11,I)
     $       *B(J,13,K)+A(K,13,I)*B(J,11,K)+A(K,14,I)*B(J,8,K)
     $       +A(K,28,I)*B(J,4,K)+A(K,31,I)*B(J,3,K)+A(K,33,I)*B(J,1,K)
            OUT(J,64,I)=OUT(J,64,I)+A(K,9,I)*B(J,14,K)+A(K,12,I)
     $       *B(J,13,K)+A(K,13,I)*B(J,12,K)+A(K,14,I)*B(J,9,K)
     $       +A(K,29,I)*B(J,4,K)+A(K,32,I)*B(J,3,K)+A(K,33,I)*B(J,2,K)
            OUT(J,65,I)=OUT(J,65,I)+A(K,10,I)*B(J,14,K)+A(K,13,I)
     $       *B(J,13,K)+A(K,14,I)*B(J,10,K)+A(K,30,I)*B(J,4,K)
     $       +A(K,33,I)*B(J,3,K)
            OUT(J,66,I)=OUT(J,66,I)+A(K,11,I)*B(J,14,K)+A(K,14,I)
     $       *B(J,11,K)+A(K,31,I)*B(J,4,K)+A(K,34,I)*B(J,1,K)
            OUT(J,67,I)=OUT(J,67,I)+A(K,12,I)*B(J,14,K)+A(K,14,I)
     $       *B(J,12,K)+A(K,32,I)*B(J,4,K)+A(K,34,I)*B(J,2,K)
            OUT(J,68,I)=OUT(J,68,I)+A(K,13,I)*B(J,14,K)+A(K,14,I)
     $       *B(J,13,K)+A(K,33,I)*B(J,4,K)+A(K,34,I)*B(J,3,K)
            OUT(J,69,I)=OUT(J,69,I)+A(K,14,I)*B(J,14,K)+A(K,34,I)
     $       *B(J,4,K)
            OUT(J,70,I)=OUT(J,70,I)+A(K,15,I)*B(J,5,K)
            OUT(J,71,I)=OUT(J,71,I)+A(K,15,I)*B(J,6,K)+A(K,16,I)
     $       *B(J,5,K)
            OUT(J,72,I)=OUT(J,72,I)+A(K,15,I)*B(J,7,K)+A(K,16,I)
     $       *B(J,6,K)+A(K,17,I)*B(J,5,K)
            OUT(J,73,I)=OUT(J,73,I)+A(K,16,I)*B(J,7,K)+A(K,17,I)
     $       *B(J,6,K)+A(K,18,I)*B(J,5,K)
            OUT(J,74,I)=OUT(J,74,I)+A(K,17,I)*B(J,7,K)+A(K,18,I)
     $       *B(J,6,K)
            OUT(J,75,I)=OUT(J,75,I)+A(K,18,I)*B(J,7,K)
            OUT(J,76,I)=OUT(J,76,I)+A(K,15,I)*B(J,8,K)+A(K,19,I)
     $       *B(J,5,K)
            OUT(J,77,I)=OUT(J,77,I)+A(K,15,I)*B(J,9,K)+A(K,16,I)
     $       *B(J,8,K)+A(K,19,I)*B(J,6,K)+A(K,20,I)*B(J,5,K)
            OUT(J,78,I)=OUT(J,78,I)+A(K,16,I)*B(J,9,K)+A(K,17,I)
     $       *B(J,8,K)+A(K,19,I)*B(J,7,K)+A(K,20,I)*B(J,6,K)+A(K,21,I)
     $       *B(J,5,K)
            OUT(J,79,I)=OUT(J,79,I)+A(K,17,I)*B(J,9,K)+A(K,18,I)
     $       *B(J,8,K)+A(K,20,I)*B(J,7,K)+A(K,21,I)*B(J,6,K)
            OUT(J,80,I)=OUT(J,80,I)+A(K,18,I)*B(J,9,K)+A(K,21,I)
     $       *B(J,7,K)
            OUT(J,81,I)=OUT(J,81,I)+A(K,15,I)*B(J,10,K)+A(K,19,I)
     $       *B(J,8,K)+A(K,22,I)*B(J,5,K)
            OUT(J,82,I)=OUT(J,82,I)+A(K,16,I)*B(J,10,K)+A(K,19,I)
     $       *B(J,9,K)+A(K,20,I)*B(J,8,K)+A(K,22,I)*B(J,6,K)+A(K,23,I)
     $       *B(J,5,K)
            OUT(J,83,I)=OUT(J,83,I)+A(K,17,I)*B(J,10,K)+A(K,20,I)
     $       *B(J,9,K)+A(K,21,I)*B(J,8,K)+A(K,22,I)*B(J,7,K)+A(K,23,I)
     $       *B(J,6,K)
            OUT(J,84,I)=OUT(J,84,I)+A(K,18,I)*B(J,10,K)+A(K,21,I)
     $       *B(J,9,K)+A(K,23,I)*B(J,7,K)
            OUT(J,85,I)=OUT(J,85,I)+A(K,19,I)*B(J,10,K)+A(K,22,I)
     $       *B(J,8,K)+A(K,24,I)*B(J,5,K)
            OUT(J,86,I)=OUT(J,86,I)+A(K,20,I)*B(J,10,K)+A(K,22,I)
     $       *B(J,9,K)+A(K,23,I)*B(J,8,K)+A(K,24,I)*B(J,6,K)
            OUT(J,87,I)=OUT(J,87,I)+A(K,21,I)*B(J,10,K)+A(K,23,I)
     $       *B(J,9,K)+A(K,24,I)*B(J,7,K)
            OUT(J,88,I)=OUT(J,88,I)+A(K,22,I)*B(J,10,K)+A(K,24,I)
     $       *B(J,8,K)
            OUT(J,89,I)=OUT(J,89,I)+A(K,23,I)*B(J,10,K)+A(K,24,I)
     $       *B(J,9,K)
            OUT(J,90,I)=OUT(J,90,I)+A(K,24,I)*B(J,10,K)
            OUT(J,91,I)=OUT(J,91,I)+A(K,15,I)*B(J,11,K)+A(K,25,I)
     $       *B(J,5,K)
            OUT(J,92,I)=OUT(J,92,I)+A(K,15,I)*B(J,12,K)+A(K,16,I)
     $       *B(J,11,K)+A(K,25,I)*B(J,6,K)+A(K,26,I)*B(J,5,K)
            OUT(J,93,I)=OUT(J,93,I)+A(K,16,I)*B(J,12,K)+A(K,17,I)
     $       *B(J,11,K)+A(K,25,I)*B(J,7,K)+A(K,26,I)*B(J,6,K)
     $       +A(K,27,I)*B(J,5,K)
            OUT(J,94,I)=OUT(J,94,I)+A(K,17,I)*B(J,12,K)+A(K,18,I)
     $       *B(J,11,K)+A(K,26,I)*B(J,7,K)+A(K,27,I)*B(J,6,K)
            OUT(J,95,I)=OUT(J,95,I)+A(K,18,I)*B(J,12,K)+A(K,27,I)
     $       *B(J,7,K)
            OUT(J,96,I)=OUT(J,96,I)+A(K,15,I)*B(J,13,K)+A(K,19,I)
     $       *B(J,11,K)+A(K,25,I)*B(J,8,K)+A(K,28,I)*B(J,5,K)
            OUT(J,97,I)=OUT(J,97,I)+A(K,16,I)*B(J,13,K)+A(K,19,I)
     $       *B(J,12,K)+A(K,20,I)*B(J,11,K)+A(K,25,I)*B(J,9,K)
     $       +A(K,26,I)*B(J,8,K)+A(K,28,I)*B(J,6,K)+A(K,29,I)*B(J,5,K)
            OUT(J,98,I)=OUT(J,98,I)+A(K,17,I)*B(J,13,K)+A(K,20,I)
     $       *B(J,12,K)+A(K,21,I)*B(J,11,K)+A(K,26,I)*B(J,9,K)
     $       +A(K,27,I)*B(J,8,K)+A(K,28,I)*B(J,7,K)+A(K,29,I)*B(J,6,K)
            OUT(J,99,I)=OUT(J,99,I)+A(K,18,I)*B(J,13,K)+A(K,21,I)
     $       *B(J,12,K)+A(K,27,I)*B(J,9,K)+A(K,29,I)*B(J,7,K)
            OUT(J,100,I)=OUT(J,100,I)+A(K,19,I)*B(J,13,K)+A(K,22,I)
     $       *B(J,11,K)+A(K,25,I)*B(J,10,K)+A(K,28,I)*B(J,8,K)
     $       +A(K,30,I)*B(J,5,K)
            OUT(J,101,I)=OUT(J,101,I)+A(K,20,I)*B(J,13,K)+A(K,22,I)
     $       *B(J,12,K)+A(K,23,I)*B(J,11,K)+A(K,26,I)*B(J,10,K)
     $       +A(K,28,I)*B(J,9,K)+A(K,29,I)*B(J,8,K)+A(K,30,I)*B(J,6,K)
            OUT(J,102,I)=OUT(J,102,I)+A(K,21,I)*B(J,13,K)+A(K,23,I)
     $       *B(J,12,K)+A(K,27,I)*B(J,10,K)+A(K,29,I)*B(J,9,K)
     $       +A(K,30,I)*B(J,7,K)
            OUT(J,103,I)=OUT(J,103,I)+A(K,22,I)*B(J,13,K)+A(K,24,I)
     $       *B(J,11,K)+A(K,28,I)*B(J,10,K)+A(K,30,I)*B(J,8,K)
            OUT(J,104,I)=OUT(J,104,I)+A(K,23,I)*B(J,13,K)+A(K,24,I)
     $       *B(J,12,K)+A(K,29,I)*B(J,10,K)+A(K,30,I)*B(J,9,K)
            OUT(J,105,I)=OUT(J,105,I)+A(K,24,I)*B(J,13,K)+A(K,30,I)
     $       *B(J,10,K)
            OUT(J,106,I)=OUT(J,106,I)+A(K,15,I)*B(J,14,K)+A(K,25,I)
     $       *B(J,11,K)+A(K,31,I)*B(J,5,K)
            OUT(J,107,I)=OUT(J,107,I)+A(K,16,I)*B(J,14,K)+A(K,25,I)
     $       *B(J,12,K)+A(K,26,I)*B(J,11,K)+A(K,31,I)*B(J,6,K)
     $       +A(K,32,I)*B(J,5,K)
            OUT(J,108,I)=OUT(J,108,I)+A(K,17,I)*B(J,14,K)+A(K,26,I)
     $       *B(J,12,K)+A(K,27,I)*B(J,11,K)+A(K,31,I)*B(J,7,K)
     $       +A(K,32,I)*B(J,6,K)
            OUT(J,109,I)=OUT(J,109,I)+A(K,18,I)*B(J,14,K)+A(K,27,I)
     $       *B(J,12,K)+A(K,32,I)*B(J,7,K)
            OUT(J,110,I)=OUT(J,110,I)+A(K,19,I)*B(J,14,K)+A(K,25,I)
     $       *B(J,13,K)+A(K,28,I)*B(J,11,K)+A(K,31,I)*B(J,8,K)
     $       +A(K,33,I)*B(J,5,K)
            OUT(J,111,I)=OUT(J,111,I)+A(K,20,I)*B(J,14,K)+A(K,26,I)
     $       *B(J,13,K)+A(K,28,I)*B(J,12,K)+A(K,29,I)*B(J,11,K)
     $       +A(K,31,I)*B(J,9,K)+A(K,32,I)*B(J,8,K)+A(K,33,I)*B(J,6,K)
            OUT(J,112,I)=OUT(J,112,I)+A(K,21,I)*B(J,14,K)+A(K,27,I)
     $       *B(J,13,K)+A(K,29,I)*B(J,12,K)+A(K,32,I)*B(J,9,K)
     $       +A(K,33,I)*B(J,7,K)
            OUT(J,113,I)=OUT(J,113,I)+A(K,22,I)*B(J,14,K)+A(K,28,I)
     $       *B(J,13,K)+A(K,30,I)*B(J,11,K)+A(K,31,I)*B(J,10,K)
     $       +A(K,33,I)*B(J,8,K)
            OUT(J,114,I)=OUT(J,114,I)+A(K,23,I)*B(J,14,K)+A(K,29,I)
     $       *B(J,13,K)+A(K,30,I)*B(J,12,K)+A(K,32,I)*B(J,10,K)
     $       +A(K,33,I)*B(J,9,K)
            OUT(J,115,I)=OUT(J,115,I)+A(K,24,I)*B(J,14,K)+A(K,30,I)
     $       *B(J,13,K)+A(K,33,I)*B(J,10,K)
            OUT(J,116,I)=OUT(J,116,I)+A(K,25,I)*B(J,14,K)+A(K,31,I)
     $       *B(J,11,K)+A(K,34,I)*B(J,5,K)
            OUT(J,117,I)=OUT(J,117,I)+A(K,26,I)*B(J,14,K)+A(K,31,I)
     $       *B(J,12,K)+A(K,32,I)*B(J,11,K)+A(K,34,I)*B(J,6,K)
            OUT(J,118,I)=OUT(J,118,I)+A(K,27,I)*B(J,14,K)+A(K,32,I)
     $       *B(J,12,K)+A(K,34,I)*B(J,7,K)
            OUT(J,119,I)=OUT(J,119,I)+A(K,28,I)*B(J,14,K)+A(K,31,I)
     $       *B(J,13,K)+A(K,33,I)*B(J,11,K)+A(K,34,I)*B(J,8,K)
            OUT(J,120,I)=OUT(J,120,I)+A(K,29,I)*B(J,14,K)+A(K,32,I)
     $       *B(J,13,K)+A(K,33,I)*B(J,12,K)+A(K,34,I)*B(J,9,K)
            OUT(J,121,I)=OUT(J,121,I)+A(K,30,I)*B(J,14,K)+A(K,33,I)
     $       *B(J,13,K)+A(K,34,I)*B(J,10,K)
            OUT(J,122,I)=OUT(J,122,I)+A(K,31,I)*B(J,14,K)+A(K,34,I)
     $       *B(J,11,K)
            OUT(J,123,I)=OUT(J,123,I)+A(K,32,I)*B(J,14,K)+A(K,34,I)
     $       *B(J,12,K)
            OUT(J,124,I)=OUT(J,124,I)+A(K,33,I)*B(J,14,K)+A(K,34,I)
     $       *B(J,13,K)
            OUT(J,125,I)=OUT(J,125,I)+A(K,34,I)*B(J,14,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_UPDATE_WL_3_2(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*32 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,125
            OUT(J,K,I)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,0,I)*B(J,5,K)+A(K,1,I)*B(J,1,K)
     $       +A(K,5,I)*B(J,0,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,0,I)*B(J,6,K)+A(K,1,I)*B(J,2,K)
     $       +A(K,2,I)*B(J,1,K)+A(K,6,I)*B(J,0,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,0,I)*B(J,7,K)+A(K,2,I)*B(J,2,K)
     $       +A(K,7,I)*B(J,0,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,0,I)*B(J,8,K)+A(K,1,I)*B(J,3,K)
     $       +A(K,3,I)*B(J,1,K)+A(K,8,I)*B(J,0,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,0,I)*B(J,9,K)+A(K,2,I)*B(J,3,K)
     $       +A(K,3,I)*B(J,2,K)+A(K,9,I)*B(J,0,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,0,I)*B(J,10,K)+A(K,3,I)
     $       *B(J,3,K)+A(K,10,I)*B(J,0,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,0,I)*B(J,11,K)+A(K,1,I)
     $       *B(J,4,K)+A(K,4,I)*B(J,1,K)+A(K,11,I)*B(J,0,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,0,I)*B(J,12,K)+A(K,2,I)
     $       *B(J,4,K)+A(K,4,I)*B(J,2,K)+A(K,12,I)*B(J,0,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,0,I)*B(J,13,K)+A(K,3,I)
     $       *B(J,4,K)+A(K,4,I)*B(J,3,K)+A(K,13,I)*B(J,0,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,0,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,4,K)+A(K,14,I)*B(J,0,K)
            OUT(J,15,I)=OUT(J,15,I)+A(K,1,I)*B(J,5,K)+A(K,5,I)
     $       *B(J,1,K)+A(K,15,I)*B(J,0,K)
            OUT(J,16,I)=OUT(J,16,I)+A(K,1,I)*B(J,6,K)+A(K,2,I)
     $       *B(J,5,K)+A(K,5,I)*B(J,2,K)+A(K,6,I)*B(J,1,K)+A(K,16,I)
     $       *B(J,0,K)
            OUT(J,17,I)=OUT(J,17,I)+A(K,1,I)*B(J,7,K)+A(K,2,I)
     $       *B(J,6,K)+A(K,6,I)*B(J,2,K)+A(K,7,I)*B(J,1,K)+A(K,17,I)
     $       *B(J,0,K)
            OUT(J,18,I)=OUT(J,18,I)+A(K,2,I)*B(J,7,K)+A(K,7,I)
     $       *B(J,2,K)+A(K,18,I)*B(J,0,K)
            OUT(J,19,I)=OUT(J,19,I)+A(K,1,I)*B(J,8,K)+A(K,3,I)
     $       *B(J,5,K)+A(K,5,I)*B(J,3,K)+A(K,8,I)*B(J,1,K)+A(K,19,I)
     $       *B(J,0,K)
            OUT(J,20,I)=OUT(J,20,I)+A(K,1,I)*B(J,9,K)+A(K,2,I)
     $       *B(J,8,K)+A(K,3,I)*B(J,6,K)+A(K,6,I)*B(J,3,K)+A(K,8,I)
     $       *B(J,2,K)+A(K,9,I)*B(J,1,K)+A(K,20,I)*B(J,0,K)
            OUT(J,21,I)=OUT(J,21,I)+A(K,2,I)*B(J,9,K)+A(K,3,I)
     $       *B(J,7,K)+A(K,7,I)*B(J,3,K)+A(K,9,I)*B(J,2,K)+A(K,21,I)
     $       *B(J,0,K)
            OUT(J,22,I)=OUT(J,22,I)+A(K,1,I)*B(J,10,K)+A(K,3,I)
     $       *B(J,8,K)+A(K,8,I)*B(J,3,K)+A(K,10,I)*B(J,1,K)+A(K,22,I)
     $       *B(J,0,K)
            OUT(J,23,I)=OUT(J,23,I)+A(K,2,I)*B(J,10,K)+A(K,3,I)
     $       *B(J,9,K)+A(K,9,I)*B(J,3,K)+A(K,10,I)*B(J,2,K)+A(K,23,I)
     $       *B(J,0,K)
            OUT(J,24,I)=OUT(J,24,I)+A(K,3,I)*B(J,10,K)+A(K,10,I)
     $       *B(J,3,K)+A(K,24,I)*B(J,0,K)
            OUT(J,25,I)=OUT(J,25,I)+A(K,1,I)*B(J,11,K)+A(K,4,I)
     $       *B(J,5,K)+A(K,5,I)*B(J,4,K)+A(K,11,I)*B(J,1,K)+A(K,25,I)
     $       *B(J,0,K)
            OUT(J,26,I)=OUT(J,26,I)+A(K,1,I)*B(J,12,K)+A(K,2,I)
     $       *B(J,11,K)+A(K,4,I)*B(J,6,K)+A(K,6,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,2,K)+A(K,12,I)*B(J,1,K)+A(K,26,I)*B(J,0,K)
            OUT(J,27,I)=OUT(J,27,I)+A(K,2,I)*B(J,12,K)+A(K,4,I)
     $       *B(J,7,K)+A(K,7,I)*B(J,4,K)+A(K,12,I)*B(J,2,K)+A(K,27,I)
     $       *B(J,0,K)
            OUT(J,28,I)=OUT(J,28,I)+A(K,1,I)*B(J,13,K)+A(K,3,I)
     $       *B(J,11,K)+A(K,4,I)*B(J,8,K)+A(K,8,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,1,K)+A(K,28,I)*B(J,0,K)
            OUT(J,29,I)=OUT(J,29,I)+A(K,2,I)*B(J,13,K)+A(K,3,I)
     $       *B(J,12,K)+A(K,4,I)*B(J,9,K)+A(K,9,I)*B(J,4,K)+A(K,12,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,2,K)+A(K,29,I)*B(J,0,K)
            OUT(J,30,I)=OUT(J,30,I)+A(K,3,I)*B(J,13,K)+A(K,4,I)
     $       *B(J,10,K)+A(K,10,I)*B(J,4,K)+A(K,13,I)*B(J,3,K)
     $       +A(K,30,I)*B(J,0,K)
            OUT(J,31,I)=OUT(J,31,I)+A(K,1,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,11,K)+A(K,11,I)*B(J,4,K)+A(K,14,I)*B(J,1,K)
     $       +A(K,31,I)*B(J,0,K)
            OUT(J,32,I)=OUT(J,32,I)+A(K,2,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,12,K)+A(K,12,I)*B(J,4,K)+A(K,14,I)*B(J,2,K)
     $       +A(K,32,I)*B(J,0,K)
            OUT(J,33,I)=OUT(J,33,I)+A(K,3,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,13,K)+A(K,13,I)*B(J,4,K)+A(K,14,I)*B(J,3,K)
     $       +A(K,33,I)*B(J,0,K)
            OUT(J,34,I)=OUT(J,34,I)+A(K,4,I)*B(J,14,K)+A(K,14,I)
     $       *B(J,4,K)+A(K,34,I)*B(J,0,K)
            OUT(J,35,I)=OUT(J,35,I)+A(K,5,I)*B(J,5,K)+A(K,15,I)
     $       *B(J,1,K)
            OUT(J,36,I)=OUT(J,36,I)+A(K,5,I)*B(J,6,K)+A(K,6,I)
     $       *B(J,5,K)+A(K,15,I)*B(J,2,K)+A(K,16,I)*B(J,1,K)
            OUT(J,37,I)=OUT(J,37,I)+A(K,5,I)*B(J,7,K)+A(K,6,I)
     $       *B(J,6,K)+A(K,7,I)*B(J,5,K)+A(K,16,I)*B(J,2,K)+A(K,17,I)
     $       *B(J,1,K)
            OUT(J,38,I)=OUT(J,38,I)+A(K,6,I)*B(J,7,K)+A(K,7,I)
     $       *B(J,6,K)+A(K,17,I)*B(J,2,K)+A(K,18,I)*B(J,1,K)
            OUT(J,39,I)=OUT(J,39,I)+A(K,7,I)*B(J,7,K)+A(K,18,I)
     $       *B(J,2,K)
            OUT(J,40,I)=OUT(J,40,I)+A(K,5,I)*B(J,8,K)+A(K,8,I)
     $       *B(J,5,K)+A(K,15,I)*B(J,3,K)+A(K,19,I)*B(J,1,K)
            OUT(J,41,I)=OUT(J,41,I)+A(K,5,I)*B(J,9,K)+A(K,6,I)
     $       *B(J,8,K)+A(K,8,I)*B(J,6,K)+A(K,9,I)*B(J,5,K)+A(K,16,I)
     $       *B(J,3,K)+A(K,19,I)*B(J,2,K)+A(K,20,I)*B(J,1,K)
            OUT(J,42,I)=OUT(J,42,I)+A(K,6,I)*B(J,9,K)+A(K,7,I)
     $       *B(J,8,K)+A(K,8,I)*B(J,7,K)+A(K,9,I)*B(J,6,K)+A(K,17,I)
     $       *B(J,3,K)+A(K,20,I)*B(J,2,K)+A(K,21,I)*B(J,1,K)
            OUT(J,43,I)=OUT(J,43,I)+A(K,7,I)*B(J,9,K)+A(K,9,I)
     $       *B(J,7,K)+A(K,18,I)*B(J,3,K)+A(K,21,I)*B(J,2,K)
            OUT(J,44,I)=OUT(J,44,I)+A(K,5,I)*B(J,10,K)+A(K,8,I)
     $       *B(J,8,K)+A(K,10,I)*B(J,5,K)+A(K,19,I)*B(J,3,K)+A(K,22,I)
     $       *B(J,1,K)
            OUT(J,45,I)=OUT(J,45,I)+A(K,6,I)*B(J,10,K)+A(K,8,I)
     $       *B(J,9,K)+A(K,9,I)*B(J,8,K)+A(K,10,I)*B(J,6,K)+A(K,20,I)
     $       *B(J,3,K)+A(K,22,I)*B(J,2,K)+A(K,23,I)*B(J,1,K)
            OUT(J,46,I)=OUT(J,46,I)+A(K,7,I)*B(J,10,K)+A(K,9,I)
     $       *B(J,9,K)+A(K,10,I)*B(J,7,K)+A(K,21,I)*B(J,3,K)+A(K,23,I)
     $       *B(J,2,K)
            OUT(J,47,I)=OUT(J,47,I)+A(K,8,I)*B(J,10,K)+A(K,10,I)
     $       *B(J,8,K)+A(K,22,I)*B(J,3,K)+A(K,24,I)*B(J,1,K)
            OUT(J,48,I)=OUT(J,48,I)+A(K,9,I)*B(J,10,K)+A(K,10,I)
     $       *B(J,9,K)+A(K,23,I)*B(J,3,K)+A(K,24,I)*B(J,2,K)
            OUT(J,49,I)=OUT(J,49,I)+A(K,10,I)*B(J,10,K)+A(K,24,I)
     $       *B(J,3,K)
            OUT(J,50,I)=OUT(J,50,I)+A(K,5,I)*B(J,11,K)+A(K,11,I)
     $       *B(J,5,K)+A(K,15,I)*B(J,4,K)+A(K,25,I)*B(J,1,K)
            OUT(J,51,I)=OUT(J,51,I)+A(K,5,I)*B(J,12,K)+A(K,6,I)
     $       *B(J,11,K)+A(K,11,I)*B(J,6,K)+A(K,12,I)*B(J,5,K)
     $       +A(K,16,I)*B(J,4,K)+A(K,25,I)*B(J,2,K)+A(K,26,I)*B(J,1,K)
            OUT(J,52,I)=OUT(J,52,I)+A(K,6,I)*B(J,12,K)+A(K,7,I)
     $       *B(J,11,K)+A(K,11,I)*B(J,7,K)+A(K,12,I)*B(J,6,K)
     $       +A(K,17,I)*B(J,4,K)+A(K,26,I)*B(J,2,K)+A(K,27,I)*B(J,1,K)
            OUT(J,53,I)=OUT(J,53,I)+A(K,7,I)*B(J,12,K)+A(K,12,I)
     $       *B(J,7,K)+A(K,18,I)*B(J,4,K)+A(K,27,I)*B(J,2,K)
            OUT(J,54,I)=OUT(J,54,I)+A(K,5,I)*B(J,13,K)+A(K,8,I)
     $       *B(J,11,K)+A(K,11,I)*B(J,8,K)+A(K,13,I)*B(J,5,K)
     $       +A(K,19,I)*B(J,4,K)+A(K,25,I)*B(J,3,K)+A(K,28,I)*B(J,1,K)
            OUT(J,55,I)=OUT(J,55,I)+A(K,6,I)*B(J,13,K)+A(K,8,I)
     $       *B(J,12,K)+A(K,9,I)*B(J,11,K)+A(K,11,I)*B(J,9,K)
     $       +A(K,12,I)*B(J,8,K)+A(K,13,I)*B(J,6,K)+A(K,20,I)*B(J,4,K)
     $       +A(K,26,I)*B(J,3,K)+A(K,28,I)*B(J,2,K)+A(K,29,I)*B(J,1,K)
            OUT(J,56,I)=OUT(J,56,I)+A(K,7,I)*B(J,13,K)+A(K,9,I)
     $       *B(J,12,K)+A(K,12,I)*B(J,9,K)+A(K,13,I)*B(J,7,K)
     $       +A(K,21,I)*B(J,4,K)+A(K,27,I)*B(J,3,K)+A(K,29,I)*B(J,2,K)
            OUT(J,57,I)=OUT(J,57,I)+A(K,8,I)*B(J,13,K)+A(K,10,I)
     $       *B(J,11,K)+A(K,11,I)*B(J,10,K)+A(K,13,I)*B(J,8,K)
     $       +A(K,22,I)*B(J,4,K)+A(K,28,I)*B(J,3,K)+A(K,30,I)*B(J,1,K)
            OUT(J,58,I)=OUT(J,58,I)+A(K,9,I)*B(J,13,K)+A(K,10,I)
     $       *B(J,12,K)+A(K,12,I)*B(J,10,K)+A(K,13,I)*B(J,9,K)
     $       +A(K,23,I)*B(J,4,K)+A(K,29,I)*B(J,3,K)+A(K,30,I)*B(J,2,K)
            OUT(J,59,I)=OUT(J,59,I)+A(K,10,I)*B(J,13,K)+A(K,13,I)
     $       *B(J,10,K)+A(K,24,I)*B(J,4,K)+A(K,30,I)*B(J,3,K)
            OUT(J,60,I)=OUT(J,60,I)+A(K,5,I)*B(J,14,K)+A(K,11,I)
     $       *B(J,11,K)+A(K,14,I)*B(J,5,K)+A(K,25,I)*B(J,4,K)
     $       +A(K,31,I)*B(J,1,K)
            OUT(J,61,I)=OUT(J,61,I)+A(K,6,I)*B(J,14,K)+A(K,11,I)
     $       *B(J,12,K)+A(K,12,I)*B(J,11,K)+A(K,14,I)*B(J,6,K)
     $       +A(K,26,I)*B(J,4,K)+A(K,31,I)*B(J,2,K)+A(K,32,I)*B(J,1,K)
            OUT(J,62,I)=OUT(J,62,I)+A(K,7,I)*B(J,14,K)+A(K,12,I)
     $       *B(J,12,K)+A(K,14,I)*B(J,7,K)+A(K,27,I)*B(J,4,K)
     $       +A(K,32,I)*B(J,2,K)
            OUT(J,63,I)=OUT(J,63,I)+A(K,8,I)*B(J,14,K)+A(K,11,I)
     $       *B(J,13,K)+A(K,13,I)*B(J,11,K)+A(K,14,I)*B(J,8,K)
     $       +A(K,28,I)*B(J,4,K)+A(K,31,I)*B(J,3,K)+A(K,33,I)*B(J,1,K)
            OUT(J,64,I)=OUT(J,64,I)+A(K,9,I)*B(J,14,K)+A(K,12,I)
     $       *B(J,13,K)+A(K,13,I)*B(J,12,K)+A(K,14,I)*B(J,9,K)
     $       +A(K,29,I)*B(J,4,K)+A(K,32,I)*B(J,3,K)+A(K,33,I)*B(J,2,K)
            OUT(J,65,I)=OUT(J,65,I)+A(K,10,I)*B(J,14,K)+A(K,13,I)
     $       *B(J,13,K)+A(K,14,I)*B(J,10,K)+A(K,30,I)*B(J,4,K)
     $       +A(K,33,I)*B(J,3,K)
            OUT(J,66,I)=OUT(J,66,I)+A(K,11,I)*B(J,14,K)+A(K,14,I)
     $       *B(J,11,K)+A(K,31,I)*B(J,4,K)+A(K,34,I)*B(J,1,K)
            OUT(J,67,I)=OUT(J,67,I)+A(K,12,I)*B(J,14,K)+A(K,14,I)
     $       *B(J,12,K)+A(K,32,I)*B(J,4,K)+A(K,34,I)*B(J,2,K)
            OUT(J,68,I)=OUT(J,68,I)+A(K,13,I)*B(J,14,K)+A(K,14,I)
     $       *B(J,13,K)+A(K,33,I)*B(J,4,K)+A(K,34,I)*B(J,3,K)
            OUT(J,69,I)=OUT(J,69,I)+A(K,14,I)*B(J,14,K)+A(K,34,I)
     $       *B(J,4,K)
            OUT(J,70,I)=OUT(J,70,I)+A(K,15,I)*B(J,5,K)
            OUT(J,71,I)=OUT(J,71,I)+A(K,15,I)*B(J,6,K)+A(K,16,I)
     $       *B(J,5,K)
            OUT(J,72,I)=OUT(J,72,I)+A(K,15,I)*B(J,7,K)+A(K,16,I)
     $       *B(J,6,K)+A(K,17,I)*B(J,5,K)
            OUT(J,73,I)=OUT(J,73,I)+A(K,16,I)*B(J,7,K)+A(K,17,I)
     $       *B(J,6,K)+A(K,18,I)*B(J,5,K)
            OUT(J,74,I)=OUT(J,74,I)+A(K,17,I)*B(J,7,K)+A(K,18,I)
     $       *B(J,6,K)
            OUT(J,75,I)=OUT(J,75,I)+A(K,18,I)*B(J,7,K)
            OUT(J,76,I)=OUT(J,76,I)+A(K,15,I)*B(J,8,K)+A(K,19,I)
     $       *B(J,5,K)
            OUT(J,77,I)=OUT(J,77,I)+A(K,15,I)*B(J,9,K)+A(K,16,I)
     $       *B(J,8,K)+A(K,19,I)*B(J,6,K)+A(K,20,I)*B(J,5,K)
            OUT(J,78,I)=OUT(J,78,I)+A(K,16,I)*B(J,9,K)+A(K,17,I)
     $       *B(J,8,K)+A(K,19,I)*B(J,7,K)+A(K,20,I)*B(J,6,K)+A(K,21,I)
     $       *B(J,5,K)
            OUT(J,79,I)=OUT(J,79,I)+A(K,17,I)*B(J,9,K)+A(K,18,I)
     $       *B(J,8,K)+A(K,20,I)*B(J,7,K)+A(K,21,I)*B(J,6,K)
            OUT(J,80,I)=OUT(J,80,I)+A(K,18,I)*B(J,9,K)+A(K,21,I)
     $       *B(J,7,K)
            OUT(J,81,I)=OUT(J,81,I)+A(K,15,I)*B(J,10,K)+A(K,19,I)
     $       *B(J,8,K)+A(K,22,I)*B(J,5,K)
            OUT(J,82,I)=OUT(J,82,I)+A(K,16,I)*B(J,10,K)+A(K,19,I)
     $       *B(J,9,K)+A(K,20,I)*B(J,8,K)+A(K,22,I)*B(J,6,K)+A(K,23,I)
     $       *B(J,5,K)
            OUT(J,83,I)=OUT(J,83,I)+A(K,17,I)*B(J,10,K)+A(K,20,I)
     $       *B(J,9,K)+A(K,21,I)*B(J,8,K)+A(K,22,I)*B(J,7,K)+A(K,23,I)
     $       *B(J,6,K)
            OUT(J,84,I)=OUT(J,84,I)+A(K,18,I)*B(J,10,K)+A(K,21,I)
     $       *B(J,9,K)+A(K,23,I)*B(J,7,K)
            OUT(J,85,I)=OUT(J,85,I)+A(K,19,I)*B(J,10,K)+A(K,22,I)
     $       *B(J,8,K)+A(K,24,I)*B(J,5,K)
            OUT(J,86,I)=OUT(J,86,I)+A(K,20,I)*B(J,10,K)+A(K,22,I)
     $       *B(J,9,K)+A(K,23,I)*B(J,8,K)+A(K,24,I)*B(J,6,K)
            OUT(J,87,I)=OUT(J,87,I)+A(K,21,I)*B(J,10,K)+A(K,23,I)
     $       *B(J,9,K)+A(K,24,I)*B(J,7,K)
            OUT(J,88,I)=OUT(J,88,I)+A(K,22,I)*B(J,10,K)+A(K,24,I)
     $       *B(J,8,K)
            OUT(J,89,I)=OUT(J,89,I)+A(K,23,I)*B(J,10,K)+A(K,24,I)
     $       *B(J,9,K)
            OUT(J,90,I)=OUT(J,90,I)+A(K,24,I)*B(J,10,K)
            OUT(J,91,I)=OUT(J,91,I)+A(K,15,I)*B(J,11,K)+A(K,25,I)
     $       *B(J,5,K)
            OUT(J,92,I)=OUT(J,92,I)+A(K,15,I)*B(J,12,K)+A(K,16,I)
     $       *B(J,11,K)+A(K,25,I)*B(J,6,K)+A(K,26,I)*B(J,5,K)
            OUT(J,93,I)=OUT(J,93,I)+A(K,16,I)*B(J,12,K)+A(K,17,I)
     $       *B(J,11,K)+A(K,25,I)*B(J,7,K)+A(K,26,I)*B(J,6,K)
     $       +A(K,27,I)*B(J,5,K)
            OUT(J,94,I)=OUT(J,94,I)+A(K,17,I)*B(J,12,K)+A(K,18,I)
     $       *B(J,11,K)+A(K,26,I)*B(J,7,K)+A(K,27,I)*B(J,6,K)
            OUT(J,95,I)=OUT(J,95,I)+A(K,18,I)*B(J,12,K)+A(K,27,I)
     $       *B(J,7,K)
            OUT(J,96,I)=OUT(J,96,I)+A(K,15,I)*B(J,13,K)+A(K,19,I)
     $       *B(J,11,K)+A(K,25,I)*B(J,8,K)+A(K,28,I)*B(J,5,K)
            OUT(J,97,I)=OUT(J,97,I)+A(K,16,I)*B(J,13,K)+A(K,19,I)
     $       *B(J,12,K)+A(K,20,I)*B(J,11,K)+A(K,25,I)*B(J,9,K)
     $       +A(K,26,I)*B(J,8,K)+A(K,28,I)*B(J,6,K)+A(K,29,I)*B(J,5,K)
            OUT(J,98,I)=OUT(J,98,I)+A(K,17,I)*B(J,13,K)+A(K,20,I)
     $       *B(J,12,K)+A(K,21,I)*B(J,11,K)+A(K,26,I)*B(J,9,K)
     $       +A(K,27,I)*B(J,8,K)+A(K,28,I)*B(J,7,K)+A(K,29,I)*B(J,6,K)
            OUT(J,99,I)=OUT(J,99,I)+A(K,18,I)*B(J,13,K)+A(K,21,I)
     $       *B(J,12,K)+A(K,27,I)*B(J,9,K)+A(K,29,I)*B(J,7,K)
            OUT(J,100,I)=OUT(J,100,I)+A(K,19,I)*B(J,13,K)+A(K,22,I)
     $       *B(J,11,K)+A(K,25,I)*B(J,10,K)+A(K,28,I)*B(J,8,K)
     $       +A(K,30,I)*B(J,5,K)
            OUT(J,101,I)=OUT(J,101,I)+A(K,20,I)*B(J,13,K)+A(K,22,I)
     $       *B(J,12,K)+A(K,23,I)*B(J,11,K)+A(K,26,I)*B(J,10,K)
     $       +A(K,28,I)*B(J,9,K)+A(K,29,I)*B(J,8,K)+A(K,30,I)*B(J,6,K)
            OUT(J,102,I)=OUT(J,102,I)+A(K,21,I)*B(J,13,K)+A(K,23,I)
     $       *B(J,12,K)+A(K,27,I)*B(J,10,K)+A(K,29,I)*B(J,9,K)
     $       +A(K,30,I)*B(J,7,K)
            OUT(J,103,I)=OUT(J,103,I)+A(K,22,I)*B(J,13,K)+A(K,24,I)
     $       *B(J,11,K)+A(K,28,I)*B(J,10,K)+A(K,30,I)*B(J,8,K)
            OUT(J,104,I)=OUT(J,104,I)+A(K,23,I)*B(J,13,K)+A(K,24,I)
     $       *B(J,12,K)+A(K,29,I)*B(J,10,K)+A(K,30,I)*B(J,9,K)
            OUT(J,105,I)=OUT(J,105,I)+A(K,24,I)*B(J,13,K)+A(K,30,I)
     $       *B(J,10,K)
            OUT(J,106,I)=OUT(J,106,I)+A(K,15,I)*B(J,14,K)+A(K,25,I)
     $       *B(J,11,K)+A(K,31,I)*B(J,5,K)
            OUT(J,107,I)=OUT(J,107,I)+A(K,16,I)*B(J,14,K)+A(K,25,I)
     $       *B(J,12,K)+A(K,26,I)*B(J,11,K)+A(K,31,I)*B(J,6,K)
     $       +A(K,32,I)*B(J,5,K)
            OUT(J,108,I)=OUT(J,108,I)+A(K,17,I)*B(J,14,K)+A(K,26,I)
     $       *B(J,12,K)+A(K,27,I)*B(J,11,K)+A(K,31,I)*B(J,7,K)
     $       +A(K,32,I)*B(J,6,K)
            OUT(J,109,I)=OUT(J,109,I)+A(K,18,I)*B(J,14,K)+A(K,27,I)
     $       *B(J,12,K)+A(K,32,I)*B(J,7,K)
            OUT(J,110,I)=OUT(J,110,I)+A(K,19,I)*B(J,14,K)+A(K,25,I)
     $       *B(J,13,K)+A(K,28,I)*B(J,11,K)+A(K,31,I)*B(J,8,K)
     $       +A(K,33,I)*B(J,5,K)
            OUT(J,111,I)=OUT(J,111,I)+A(K,20,I)*B(J,14,K)+A(K,26,I)
     $       *B(J,13,K)+A(K,28,I)*B(J,12,K)+A(K,29,I)*B(J,11,K)
     $       +A(K,31,I)*B(J,9,K)+A(K,32,I)*B(J,8,K)+A(K,33,I)*B(J,6,K)
            OUT(J,112,I)=OUT(J,112,I)+A(K,21,I)*B(J,14,K)+A(K,27,I)
     $       *B(J,13,K)+A(K,29,I)*B(J,12,K)+A(K,32,I)*B(J,9,K)
     $       +A(K,33,I)*B(J,7,K)
            OUT(J,113,I)=OUT(J,113,I)+A(K,22,I)*B(J,14,K)+A(K,28,I)
     $       *B(J,13,K)+A(K,30,I)*B(J,11,K)+A(K,31,I)*B(J,10,K)
     $       +A(K,33,I)*B(J,8,K)
            OUT(J,114,I)=OUT(J,114,I)+A(K,23,I)*B(J,14,K)+A(K,29,I)
     $       *B(J,13,K)+A(K,30,I)*B(J,12,K)+A(K,32,I)*B(J,10,K)
     $       +A(K,33,I)*B(J,9,K)
            OUT(J,115,I)=OUT(J,115,I)+A(K,24,I)*B(J,14,K)+A(K,30,I)
     $       *B(J,13,K)+A(K,33,I)*B(J,10,K)
            OUT(J,116,I)=OUT(J,116,I)+A(K,25,I)*B(J,14,K)+A(K,31,I)
     $       *B(J,11,K)+A(K,34,I)*B(J,5,K)
            OUT(J,117,I)=OUT(J,117,I)+A(K,26,I)*B(J,14,K)+A(K,31,I)
     $       *B(J,12,K)+A(K,32,I)*B(J,11,K)+A(K,34,I)*B(J,6,K)
            OUT(J,118,I)=OUT(J,118,I)+A(K,27,I)*B(J,14,K)+A(K,32,I)
     $       *B(J,12,K)+A(K,34,I)*B(J,7,K)
            OUT(J,119,I)=OUT(J,119,I)+A(K,28,I)*B(J,14,K)+A(K,31,I)
     $       *B(J,13,K)+A(K,33,I)*B(J,11,K)+A(K,34,I)*B(J,8,K)
            OUT(J,120,I)=OUT(J,120,I)+A(K,29,I)*B(J,14,K)+A(K,32,I)
     $       *B(J,13,K)+A(K,33,I)*B(J,12,K)+A(K,34,I)*B(J,9,K)
            OUT(J,121,I)=OUT(J,121,I)+A(K,30,I)*B(J,14,K)+A(K,33,I)
     $       *B(J,13,K)+A(K,34,I)*B(J,10,K)
            OUT(J,122,I)=OUT(J,122,I)+A(K,31,I)*B(J,14,K)+A(K,34,I)
     $       *B(J,11,K)
            OUT(J,123,I)=OUT(J,123,I)+A(K,32,I)*B(J,14,K)+A(K,34,I)
     $       *B(J,12,K)
            OUT(J,124,I)=OUT(J,124,I)+A(K,33,I)*B(J,14,K)+A(K,34,I)
     $       *B(J,13,K)
            OUT(J,125,I)=OUT(J,125,I)+A(K,34,I)*B(J,14,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE UPDATE_WL_0_0(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*16 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,0
            OUT(J,K,I)=(0.0D0,0.0D0)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_UPDATE_WL_0_0(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*32 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,0
            OUT(J,K,I)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE UPDATE_WL_3_0(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*16 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,34
            OUT(J,K,I)=(0.0D0,0.0D0)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,5,I)*B(J,0,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,6,I)*B(J,0,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,7,I)*B(J,0,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,8,I)*B(J,0,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,9,I)*B(J,0,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,10,I)*B(J,0,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,11,I)*B(J,0,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,12,I)*B(J,0,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,13,I)*B(J,0,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,14,I)*B(J,0,K)
            OUT(J,15,I)=OUT(J,15,I)+A(K,15,I)*B(J,0,K)
            OUT(J,16,I)=OUT(J,16,I)+A(K,16,I)*B(J,0,K)
            OUT(J,17,I)=OUT(J,17,I)+A(K,17,I)*B(J,0,K)
            OUT(J,18,I)=OUT(J,18,I)+A(K,18,I)*B(J,0,K)
            OUT(J,19,I)=OUT(J,19,I)+A(K,19,I)*B(J,0,K)
            OUT(J,20,I)=OUT(J,20,I)+A(K,20,I)*B(J,0,K)
            OUT(J,21,I)=OUT(J,21,I)+A(K,21,I)*B(J,0,K)
            OUT(J,22,I)=OUT(J,22,I)+A(K,22,I)*B(J,0,K)
            OUT(J,23,I)=OUT(J,23,I)+A(K,23,I)*B(J,0,K)
            OUT(J,24,I)=OUT(J,24,I)+A(K,24,I)*B(J,0,K)
            OUT(J,25,I)=OUT(J,25,I)+A(K,25,I)*B(J,0,K)
            OUT(J,26,I)=OUT(J,26,I)+A(K,26,I)*B(J,0,K)
            OUT(J,27,I)=OUT(J,27,I)+A(K,27,I)*B(J,0,K)
            OUT(J,28,I)=OUT(J,28,I)+A(K,28,I)*B(J,0,K)
            OUT(J,29,I)=OUT(J,29,I)+A(K,29,I)*B(J,0,K)
            OUT(J,30,I)=OUT(J,30,I)+A(K,30,I)*B(J,0,K)
            OUT(J,31,I)=OUT(J,31,I)+A(K,31,I)*B(J,0,K)
            OUT(J,32,I)=OUT(J,32,I)+A(K,32,I)*B(J,0,K)
            OUT(J,33,I)=OUT(J,33,I)+A(K,33,I)*B(J,0,K)
            OUT(J,34,I)=OUT(J,34,I)+A(K,34,I)*B(J,0,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_UPDATE_WL_3_0(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*32 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,34
            OUT(J,K,I)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,5,I)*B(J,0,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,6,I)*B(J,0,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,7,I)*B(J,0,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,8,I)*B(J,0,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,9,I)*B(J,0,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,10,I)*B(J,0,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,11,I)*B(J,0,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,12,I)*B(J,0,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,13,I)*B(J,0,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,14,I)*B(J,0,K)
            OUT(J,15,I)=OUT(J,15,I)+A(K,15,I)*B(J,0,K)
            OUT(J,16,I)=OUT(J,16,I)+A(K,16,I)*B(J,0,K)
            OUT(J,17,I)=OUT(J,17,I)+A(K,17,I)*B(J,0,K)
            OUT(J,18,I)=OUT(J,18,I)+A(K,18,I)*B(J,0,K)
            OUT(J,19,I)=OUT(J,19,I)+A(K,19,I)*B(J,0,K)
            OUT(J,20,I)=OUT(J,20,I)+A(K,20,I)*B(J,0,K)
            OUT(J,21,I)=OUT(J,21,I)+A(K,21,I)*B(J,0,K)
            OUT(J,22,I)=OUT(J,22,I)+A(K,22,I)*B(J,0,K)
            OUT(J,23,I)=OUT(J,23,I)+A(K,23,I)*B(J,0,K)
            OUT(J,24,I)=OUT(J,24,I)+A(K,24,I)*B(J,0,K)
            OUT(J,25,I)=OUT(J,25,I)+A(K,25,I)*B(J,0,K)
            OUT(J,26,I)=OUT(J,26,I)+A(K,26,I)*B(J,0,K)
            OUT(J,27,I)=OUT(J,27,I)+A(K,27,I)*B(J,0,K)
            OUT(J,28,I)=OUT(J,28,I)+A(K,28,I)*B(J,0,K)
            OUT(J,29,I)=OUT(J,29,I)+A(K,29,I)*B(J,0,K)
            OUT(J,30,I)=OUT(J,30,I)+A(K,30,I)*B(J,0,K)
            OUT(J,31,I)=OUT(J,31,I)+A(K,31,I)*B(J,0,K)
            OUT(J,32,I)=OUT(J,32,I)+A(K,32,I)*B(J,0,K)
            OUT(J,33,I)=OUT(J,33,I)+A(K,33,I)*B(J,0,K)
            OUT(J,34,I)=OUT(J,34,I)+A(K,34,I)*B(J,0,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE UPDATE_WL_3_1(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*16 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,69
            OUT(J,K,I)=(0.0D0,0.0D0)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,1,I)*B(J,1,K)+A(K,5,I)*B(J,0,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,1,I)*B(J,2,K)+A(K,2,I)*B(J,1,K)
     $       +A(K,6,I)*B(J,0,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,2,I)*B(J,2,K)+A(K,7,I)*B(J,0,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,1,I)*B(J,3,K)+A(K,3,I)*B(J,1,K)
     $       +A(K,8,I)*B(J,0,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,2,I)*B(J,3,K)+A(K,3,I)*B(J,2,K)
     $       +A(K,9,I)*B(J,0,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,3,I)*B(J,3,K)+A(K,10,I)
     $       *B(J,0,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,1,I)*B(J,4,K)+A(K,4,I)
     $       *B(J,1,K)+A(K,11,I)*B(J,0,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,2,I)*B(J,4,K)+A(K,4,I)
     $       *B(J,2,K)+A(K,12,I)*B(J,0,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,3,I)*B(J,4,K)+A(K,4,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,0,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,4,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,0,K)
            OUT(J,15,I)=OUT(J,15,I)+A(K,5,I)*B(J,1,K)+A(K,15,I)
     $       *B(J,0,K)
            OUT(J,16,I)=OUT(J,16,I)+A(K,5,I)*B(J,2,K)+A(K,6,I)
     $       *B(J,1,K)+A(K,16,I)*B(J,0,K)
            OUT(J,17,I)=OUT(J,17,I)+A(K,6,I)*B(J,2,K)+A(K,7,I)
     $       *B(J,1,K)+A(K,17,I)*B(J,0,K)
            OUT(J,18,I)=OUT(J,18,I)+A(K,7,I)*B(J,2,K)+A(K,18,I)
     $       *B(J,0,K)
            OUT(J,19,I)=OUT(J,19,I)+A(K,5,I)*B(J,3,K)+A(K,8,I)
     $       *B(J,1,K)+A(K,19,I)*B(J,0,K)
            OUT(J,20,I)=OUT(J,20,I)+A(K,6,I)*B(J,3,K)+A(K,8,I)
     $       *B(J,2,K)+A(K,9,I)*B(J,1,K)+A(K,20,I)*B(J,0,K)
            OUT(J,21,I)=OUT(J,21,I)+A(K,7,I)*B(J,3,K)+A(K,9,I)
     $       *B(J,2,K)+A(K,21,I)*B(J,0,K)
            OUT(J,22,I)=OUT(J,22,I)+A(K,8,I)*B(J,3,K)+A(K,10,I)
     $       *B(J,1,K)+A(K,22,I)*B(J,0,K)
            OUT(J,23,I)=OUT(J,23,I)+A(K,9,I)*B(J,3,K)+A(K,10,I)
     $       *B(J,2,K)+A(K,23,I)*B(J,0,K)
            OUT(J,24,I)=OUT(J,24,I)+A(K,10,I)*B(J,3,K)+A(K,24,I)
     $       *B(J,0,K)
            OUT(J,25,I)=OUT(J,25,I)+A(K,5,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,1,K)+A(K,25,I)*B(J,0,K)
            OUT(J,26,I)=OUT(J,26,I)+A(K,6,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,2,K)+A(K,12,I)*B(J,1,K)+A(K,26,I)*B(J,0,K)
            OUT(J,27,I)=OUT(J,27,I)+A(K,7,I)*B(J,4,K)+A(K,12,I)
     $       *B(J,2,K)+A(K,27,I)*B(J,0,K)
            OUT(J,28,I)=OUT(J,28,I)+A(K,8,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,1,K)+A(K,28,I)*B(J,0,K)
            OUT(J,29,I)=OUT(J,29,I)+A(K,9,I)*B(J,4,K)+A(K,12,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,2,K)+A(K,29,I)*B(J,0,K)
            OUT(J,30,I)=OUT(J,30,I)+A(K,10,I)*B(J,4,K)+A(K,13,I)
     $       *B(J,3,K)+A(K,30,I)*B(J,0,K)
            OUT(J,31,I)=OUT(J,31,I)+A(K,11,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,1,K)+A(K,31,I)*B(J,0,K)
            OUT(J,32,I)=OUT(J,32,I)+A(K,12,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,2,K)+A(K,32,I)*B(J,0,K)
            OUT(J,33,I)=OUT(J,33,I)+A(K,13,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,3,K)+A(K,33,I)*B(J,0,K)
            OUT(J,34,I)=OUT(J,34,I)+A(K,14,I)*B(J,4,K)+A(K,34,I)
     $       *B(J,0,K)
            OUT(J,35,I)=OUT(J,35,I)+A(K,15,I)*B(J,1,K)
            OUT(J,36,I)=OUT(J,36,I)+A(K,15,I)*B(J,2,K)+A(K,16,I)
     $       *B(J,1,K)
            OUT(J,37,I)=OUT(J,37,I)+A(K,16,I)*B(J,2,K)+A(K,17,I)
     $       *B(J,1,K)
            OUT(J,38,I)=OUT(J,38,I)+A(K,17,I)*B(J,2,K)+A(K,18,I)
     $       *B(J,1,K)
            OUT(J,39,I)=OUT(J,39,I)+A(K,18,I)*B(J,2,K)
            OUT(J,40,I)=OUT(J,40,I)+A(K,15,I)*B(J,3,K)+A(K,19,I)
     $       *B(J,1,K)
            OUT(J,41,I)=OUT(J,41,I)+A(K,16,I)*B(J,3,K)+A(K,19,I)
     $       *B(J,2,K)+A(K,20,I)*B(J,1,K)
            OUT(J,42,I)=OUT(J,42,I)+A(K,17,I)*B(J,3,K)+A(K,20,I)
     $       *B(J,2,K)+A(K,21,I)*B(J,1,K)
            OUT(J,43,I)=OUT(J,43,I)+A(K,18,I)*B(J,3,K)+A(K,21,I)
     $       *B(J,2,K)
            OUT(J,44,I)=OUT(J,44,I)+A(K,19,I)*B(J,3,K)+A(K,22,I)
     $       *B(J,1,K)
            OUT(J,45,I)=OUT(J,45,I)+A(K,20,I)*B(J,3,K)+A(K,22,I)
     $       *B(J,2,K)+A(K,23,I)*B(J,1,K)
            OUT(J,46,I)=OUT(J,46,I)+A(K,21,I)*B(J,3,K)+A(K,23,I)
     $       *B(J,2,K)
            OUT(J,47,I)=OUT(J,47,I)+A(K,22,I)*B(J,3,K)+A(K,24,I)
     $       *B(J,1,K)
            OUT(J,48,I)=OUT(J,48,I)+A(K,23,I)*B(J,3,K)+A(K,24,I)
     $       *B(J,2,K)
            OUT(J,49,I)=OUT(J,49,I)+A(K,24,I)*B(J,3,K)
            OUT(J,50,I)=OUT(J,50,I)+A(K,15,I)*B(J,4,K)+A(K,25,I)
     $       *B(J,1,K)
            OUT(J,51,I)=OUT(J,51,I)+A(K,16,I)*B(J,4,K)+A(K,25,I)
     $       *B(J,2,K)+A(K,26,I)*B(J,1,K)
            OUT(J,52,I)=OUT(J,52,I)+A(K,17,I)*B(J,4,K)+A(K,26,I)
     $       *B(J,2,K)+A(K,27,I)*B(J,1,K)
            OUT(J,53,I)=OUT(J,53,I)+A(K,18,I)*B(J,4,K)+A(K,27,I)
     $       *B(J,2,K)
            OUT(J,54,I)=OUT(J,54,I)+A(K,19,I)*B(J,4,K)+A(K,25,I)
     $       *B(J,3,K)+A(K,28,I)*B(J,1,K)
            OUT(J,55,I)=OUT(J,55,I)+A(K,20,I)*B(J,4,K)+A(K,26,I)
     $       *B(J,3,K)+A(K,28,I)*B(J,2,K)+A(K,29,I)*B(J,1,K)
            OUT(J,56,I)=OUT(J,56,I)+A(K,21,I)*B(J,4,K)+A(K,27,I)
     $       *B(J,3,K)+A(K,29,I)*B(J,2,K)
            OUT(J,57,I)=OUT(J,57,I)+A(K,22,I)*B(J,4,K)+A(K,28,I)
     $       *B(J,3,K)+A(K,30,I)*B(J,1,K)
            OUT(J,58,I)=OUT(J,58,I)+A(K,23,I)*B(J,4,K)+A(K,29,I)
     $       *B(J,3,K)+A(K,30,I)*B(J,2,K)
            OUT(J,59,I)=OUT(J,59,I)+A(K,24,I)*B(J,4,K)+A(K,30,I)
     $       *B(J,3,K)
            OUT(J,60,I)=OUT(J,60,I)+A(K,25,I)*B(J,4,K)+A(K,31,I)
     $       *B(J,1,K)
            OUT(J,61,I)=OUT(J,61,I)+A(K,26,I)*B(J,4,K)+A(K,31,I)
     $       *B(J,2,K)+A(K,32,I)*B(J,1,K)
            OUT(J,62,I)=OUT(J,62,I)+A(K,27,I)*B(J,4,K)+A(K,32,I)
     $       *B(J,2,K)
            OUT(J,63,I)=OUT(J,63,I)+A(K,28,I)*B(J,4,K)+A(K,31,I)
     $       *B(J,3,K)+A(K,33,I)*B(J,1,K)
            OUT(J,64,I)=OUT(J,64,I)+A(K,29,I)*B(J,4,K)+A(K,32,I)
     $       *B(J,3,K)+A(K,33,I)*B(J,2,K)
            OUT(J,65,I)=OUT(J,65,I)+A(K,30,I)*B(J,4,K)+A(K,33,I)
     $       *B(J,3,K)
            OUT(J,66,I)=OUT(J,66,I)+A(K,31,I)*B(J,4,K)+A(K,34,I)
     $       *B(J,1,K)
            OUT(J,67,I)=OUT(J,67,I)+A(K,32,I)*B(J,4,K)+A(K,34,I)
     $       *B(J,2,K)
            OUT(J,68,I)=OUT(J,68,I)+A(K,33,I)*B(J,4,K)+A(K,34,I)
     $       *B(J,3,K)
            OUT(J,69,I)=OUT(J,69,I)+A(K,34,I)*B(J,4,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_UPDATE_WL_3_1(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*32 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,69
            OUT(J,K,I)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,1,I)*B(J,1,K)+A(K,5,I)*B(J,0,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,1,I)*B(J,2,K)+A(K,2,I)*B(J,1,K)
     $       +A(K,6,I)*B(J,0,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,2,I)*B(J,2,K)+A(K,7,I)*B(J,0,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,1,I)*B(J,3,K)+A(K,3,I)*B(J,1,K)
     $       +A(K,8,I)*B(J,0,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,2,I)*B(J,3,K)+A(K,3,I)*B(J,2,K)
     $       +A(K,9,I)*B(J,0,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,3,I)*B(J,3,K)+A(K,10,I)
     $       *B(J,0,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,1,I)*B(J,4,K)+A(K,4,I)
     $       *B(J,1,K)+A(K,11,I)*B(J,0,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,2,I)*B(J,4,K)+A(K,4,I)
     $       *B(J,2,K)+A(K,12,I)*B(J,0,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,3,I)*B(J,4,K)+A(K,4,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,0,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,4,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,0,K)
            OUT(J,15,I)=OUT(J,15,I)+A(K,5,I)*B(J,1,K)+A(K,15,I)
     $       *B(J,0,K)
            OUT(J,16,I)=OUT(J,16,I)+A(K,5,I)*B(J,2,K)+A(K,6,I)
     $       *B(J,1,K)+A(K,16,I)*B(J,0,K)
            OUT(J,17,I)=OUT(J,17,I)+A(K,6,I)*B(J,2,K)+A(K,7,I)
     $       *B(J,1,K)+A(K,17,I)*B(J,0,K)
            OUT(J,18,I)=OUT(J,18,I)+A(K,7,I)*B(J,2,K)+A(K,18,I)
     $       *B(J,0,K)
            OUT(J,19,I)=OUT(J,19,I)+A(K,5,I)*B(J,3,K)+A(K,8,I)
     $       *B(J,1,K)+A(K,19,I)*B(J,0,K)
            OUT(J,20,I)=OUT(J,20,I)+A(K,6,I)*B(J,3,K)+A(K,8,I)
     $       *B(J,2,K)+A(K,9,I)*B(J,1,K)+A(K,20,I)*B(J,0,K)
            OUT(J,21,I)=OUT(J,21,I)+A(K,7,I)*B(J,3,K)+A(K,9,I)
     $       *B(J,2,K)+A(K,21,I)*B(J,0,K)
            OUT(J,22,I)=OUT(J,22,I)+A(K,8,I)*B(J,3,K)+A(K,10,I)
     $       *B(J,1,K)+A(K,22,I)*B(J,0,K)
            OUT(J,23,I)=OUT(J,23,I)+A(K,9,I)*B(J,3,K)+A(K,10,I)
     $       *B(J,2,K)+A(K,23,I)*B(J,0,K)
            OUT(J,24,I)=OUT(J,24,I)+A(K,10,I)*B(J,3,K)+A(K,24,I)
     $       *B(J,0,K)
            OUT(J,25,I)=OUT(J,25,I)+A(K,5,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,1,K)+A(K,25,I)*B(J,0,K)
            OUT(J,26,I)=OUT(J,26,I)+A(K,6,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,2,K)+A(K,12,I)*B(J,1,K)+A(K,26,I)*B(J,0,K)
            OUT(J,27,I)=OUT(J,27,I)+A(K,7,I)*B(J,4,K)+A(K,12,I)
     $       *B(J,2,K)+A(K,27,I)*B(J,0,K)
            OUT(J,28,I)=OUT(J,28,I)+A(K,8,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,1,K)+A(K,28,I)*B(J,0,K)
            OUT(J,29,I)=OUT(J,29,I)+A(K,9,I)*B(J,4,K)+A(K,12,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,2,K)+A(K,29,I)*B(J,0,K)
            OUT(J,30,I)=OUT(J,30,I)+A(K,10,I)*B(J,4,K)+A(K,13,I)
     $       *B(J,3,K)+A(K,30,I)*B(J,0,K)
            OUT(J,31,I)=OUT(J,31,I)+A(K,11,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,1,K)+A(K,31,I)*B(J,0,K)
            OUT(J,32,I)=OUT(J,32,I)+A(K,12,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,2,K)+A(K,32,I)*B(J,0,K)
            OUT(J,33,I)=OUT(J,33,I)+A(K,13,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,3,K)+A(K,33,I)*B(J,0,K)
            OUT(J,34,I)=OUT(J,34,I)+A(K,14,I)*B(J,4,K)+A(K,34,I)
     $       *B(J,0,K)
            OUT(J,35,I)=OUT(J,35,I)+A(K,15,I)*B(J,1,K)
            OUT(J,36,I)=OUT(J,36,I)+A(K,15,I)*B(J,2,K)+A(K,16,I)
     $       *B(J,1,K)
            OUT(J,37,I)=OUT(J,37,I)+A(K,16,I)*B(J,2,K)+A(K,17,I)
     $       *B(J,1,K)
            OUT(J,38,I)=OUT(J,38,I)+A(K,17,I)*B(J,2,K)+A(K,18,I)
     $       *B(J,1,K)
            OUT(J,39,I)=OUT(J,39,I)+A(K,18,I)*B(J,2,K)
            OUT(J,40,I)=OUT(J,40,I)+A(K,15,I)*B(J,3,K)+A(K,19,I)
     $       *B(J,1,K)
            OUT(J,41,I)=OUT(J,41,I)+A(K,16,I)*B(J,3,K)+A(K,19,I)
     $       *B(J,2,K)+A(K,20,I)*B(J,1,K)
            OUT(J,42,I)=OUT(J,42,I)+A(K,17,I)*B(J,3,K)+A(K,20,I)
     $       *B(J,2,K)+A(K,21,I)*B(J,1,K)
            OUT(J,43,I)=OUT(J,43,I)+A(K,18,I)*B(J,3,K)+A(K,21,I)
     $       *B(J,2,K)
            OUT(J,44,I)=OUT(J,44,I)+A(K,19,I)*B(J,3,K)+A(K,22,I)
     $       *B(J,1,K)
            OUT(J,45,I)=OUT(J,45,I)+A(K,20,I)*B(J,3,K)+A(K,22,I)
     $       *B(J,2,K)+A(K,23,I)*B(J,1,K)
            OUT(J,46,I)=OUT(J,46,I)+A(K,21,I)*B(J,3,K)+A(K,23,I)
     $       *B(J,2,K)
            OUT(J,47,I)=OUT(J,47,I)+A(K,22,I)*B(J,3,K)+A(K,24,I)
     $       *B(J,1,K)
            OUT(J,48,I)=OUT(J,48,I)+A(K,23,I)*B(J,3,K)+A(K,24,I)
     $       *B(J,2,K)
            OUT(J,49,I)=OUT(J,49,I)+A(K,24,I)*B(J,3,K)
            OUT(J,50,I)=OUT(J,50,I)+A(K,15,I)*B(J,4,K)+A(K,25,I)
     $       *B(J,1,K)
            OUT(J,51,I)=OUT(J,51,I)+A(K,16,I)*B(J,4,K)+A(K,25,I)
     $       *B(J,2,K)+A(K,26,I)*B(J,1,K)
            OUT(J,52,I)=OUT(J,52,I)+A(K,17,I)*B(J,4,K)+A(K,26,I)
     $       *B(J,2,K)+A(K,27,I)*B(J,1,K)
            OUT(J,53,I)=OUT(J,53,I)+A(K,18,I)*B(J,4,K)+A(K,27,I)
     $       *B(J,2,K)
            OUT(J,54,I)=OUT(J,54,I)+A(K,19,I)*B(J,4,K)+A(K,25,I)
     $       *B(J,3,K)+A(K,28,I)*B(J,1,K)
            OUT(J,55,I)=OUT(J,55,I)+A(K,20,I)*B(J,4,K)+A(K,26,I)
     $       *B(J,3,K)+A(K,28,I)*B(J,2,K)+A(K,29,I)*B(J,1,K)
            OUT(J,56,I)=OUT(J,56,I)+A(K,21,I)*B(J,4,K)+A(K,27,I)
     $       *B(J,3,K)+A(K,29,I)*B(J,2,K)
            OUT(J,57,I)=OUT(J,57,I)+A(K,22,I)*B(J,4,K)+A(K,28,I)
     $       *B(J,3,K)+A(K,30,I)*B(J,1,K)
            OUT(J,58,I)=OUT(J,58,I)+A(K,23,I)*B(J,4,K)+A(K,29,I)
     $       *B(J,3,K)+A(K,30,I)*B(J,2,K)
            OUT(J,59,I)=OUT(J,59,I)+A(K,24,I)*B(J,4,K)+A(K,30,I)
     $       *B(J,3,K)
            OUT(J,60,I)=OUT(J,60,I)+A(K,25,I)*B(J,4,K)+A(K,31,I)
     $       *B(J,1,K)
            OUT(J,61,I)=OUT(J,61,I)+A(K,26,I)*B(J,4,K)+A(K,31,I)
     $       *B(J,2,K)+A(K,32,I)*B(J,1,K)
            OUT(J,62,I)=OUT(J,62,I)+A(K,27,I)*B(J,4,K)+A(K,32,I)
     $       *B(J,2,K)
            OUT(J,63,I)=OUT(J,63,I)+A(K,28,I)*B(J,4,K)+A(K,31,I)
     $       *B(J,3,K)+A(K,33,I)*B(J,1,K)
            OUT(J,64,I)=OUT(J,64,I)+A(K,29,I)*B(J,4,K)+A(K,32,I)
     $       *B(J,3,K)+A(K,33,I)*B(J,2,K)
            OUT(J,65,I)=OUT(J,65,I)+A(K,30,I)*B(J,4,K)+A(K,33,I)
     $       *B(J,3,K)
            OUT(J,66,I)=OUT(J,66,I)+A(K,31,I)*B(J,4,K)+A(K,34,I)
     $       *B(J,1,K)
            OUT(J,67,I)=OUT(J,67,I)+A(K,32,I)*B(J,4,K)+A(K,34,I)
     $       *B(J,2,K)
            OUT(J,68,I)=OUT(J,68,I)+A(K,33,I)*B(J,4,K)+A(K,34,I)
     $       *B(J,3,K)
            OUT(J,69,I)=OUT(J,69,I)+A(K,34,I)*B(J,4,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE UPDATE_WL_2_1(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*16 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,34
            OUT(J,K,I)=(0.0D0,0.0D0)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,1,I)*B(J,1,K)+A(K,5,I)*B(J,0,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,1,I)*B(J,2,K)+A(K,2,I)*B(J,1,K)
     $       +A(K,6,I)*B(J,0,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,2,I)*B(J,2,K)+A(K,7,I)*B(J,0,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,1,I)*B(J,3,K)+A(K,3,I)*B(J,1,K)
     $       +A(K,8,I)*B(J,0,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,2,I)*B(J,3,K)+A(K,3,I)*B(J,2,K)
     $       +A(K,9,I)*B(J,0,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,3,I)*B(J,3,K)+A(K,10,I)
     $       *B(J,0,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,1,I)*B(J,4,K)+A(K,4,I)
     $       *B(J,1,K)+A(K,11,I)*B(J,0,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,2,I)*B(J,4,K)+A(K,4,I)
     $       *B(J,2,K)+A(K,12,I)*B(J,0,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,3,I)*B(J,4,K)+A(K,4,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,0,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,4,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,0,K)
            OUT(J,15,I)=OUT(J,15,I)+A(K,5,I)*B(J,1,K)
            OUT(J,16,I)=OUT(J,16,I)+A(K,5,I)*B(J,2,K)+A(K,6,I)*B(J,1,K)
            OUT(J,17,I)=OUT(J,17,I)+A(K,6,I)*B(J,2,K)+A(K,7,I)*B(J,1,K)
            OUT(J,18,I)=OUT(J,18,I)+A(K,7,I)*B(J,2,K)
            OUT(J,19,I)=OUT(J,19,I)+A(K,5,I)*B(J,3,K)+A(K,8,I)*B(J,1,K)
            OUT(J,20,I)=OUT(J,20,I)+A(K,6,I)*B(J,3,K)+A(K,8,I)
     $       *B(J,2,K)+A(K,9,I)*B(J,1,K)
            OUT(J,21,I)=OUT(J,21,I)+A(K,7,I)*B(J,3,K)+A(K,9,I)*B(J,2,K)
            OUT(J,22,I)=OUT(J,22,I)+A(K,8,I)*B(J,3,K)+A(K,10,I)
     $       *B(J,1,K)
            OUT(J,23,I)=OUT(J,23,I)+A(K,9,I)*B(J,3,K)+A(K,10,I)
     $       *B(J,2,K)
            OUT(J,24,I)=OUT(J,24,I)+A(K,10,I)*B(J,3,K)
            OUT(J,25,I)=OUT(J,25,I)+A(K,5,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,1,K)
            OUT(J,26,I)=OUT(J,26,I)+A(K,6,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,2,K)+A(K,12,I)*B(J,1,K)
            OUT(J,27,I)=OUT(J,27,I)+A(K,7,I)*B(J,4,K)+A(K,12,I)
     $       *B(J,2,K)
            OUT(J,28,I)=OUT(J,28,I)+A(K,8,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,1,K)
            OUT(J,29,I)=OUT(J,29,I)+A(K,9,I)*B(J,4,K)+A(K,12,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,2,K)
            OUT(J,30,I)=OUT(J,30,I)+A(K,10,I)*B(J,4,K)+A(K,13,I)
     $       *B(J,3,K)
            OUT(J,31,I)=OUT(J,31,I)+A(K,11,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,1,K)
            OUT(J,32,I)=OUT(J,32,I)+A(K,12,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,2,K)
            OUT(J,33,I)=OUT(J,33,I)+A(K,13,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,3,K)
            OUT(J,34,I)=OUT(J,34,I)+A(K,14,I)*B(J,4,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_UPDATE_WL_2_1(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*32 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,34
            OUT(J,K,I)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,1,I)*B(J,1,K)+A(K,5,I)*B(J,0,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,1,I)*B(J,2,K)+A(K,2,I)*B(J,1,K)
     $       +A(K,6,I)*B(J,0,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,2,I)*B(J,2,K)+A(K,7,I)*B(J,0,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,1,I)*B(J,3,K)+A(K,3,I)*B(J,1,K)
     $       +A(K,8,I)*B(J,0,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,2,I)*B(J,3,K)+A(K,3,I)*B(J,2,K)
     $       +A(K,9,I)*B(J,0,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,3,I)*B(J,3,K)+A(K,10,I)
     $       *B(J,0,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,1,I)*B(J,4,K)+A(K,4,I)
     $       *B(J,1,K)+A(K,11,I)*B(J,0,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,2,I)*B(J,4,K)+A(K,4,I)
     $       *B(J,2,K)+A(K,12,I)*B(J,0,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,3,I)*B(J,4,K)+A(K,4,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,0,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,4,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,0,K)
            OUT(J,15,I)=OUT(J,15,I)+A(K,5,I)*B(J,1,K)
            OUT(J,16,I)=OUT(J,16,I)+A(K,5,I)*B(J,2,K)+A(K,6,I)*B(J,1,K)
            OUT(J,17,I)=OUT(J,17,I)+A(K,6,I)*B(J,2,K)+A(K,7,I)*B(J,1,K)
            OUT(J,18,I)=OUT(J,18,I)+A(K,7,I)*B(J,2,K)
            OUT(J,19,I)=OUT(J,19,I)+A(K,5,I)*B(J,3,K)+A(K,8,I)*B(J,1,K)
            OUT(J,20,I)=OUT(J,20,I)+A(K,6,I)*B(J,3,K)+A(K,8,I)
     $       *B(J,2,K)+A(K,9,I)*B(J,1,K)
            OUT(J,21,I)=OUT(J,21,I)+A(K,7,I)*B(J,3,K)+A(K,9,I)*B(J,2,K)
            OUT(J,22,I)=OUT(J,22,I)+A(K,8,I)*B(J,3,K)+A(K,10,I)
     $       *B(J,1,K)
            OUT(J,23,I)=OUT(J,23,I)+A(K,9,I)*B(J,3,K)+A(K,10,I)
     $       *B(J,2,K)
            OUT(J,24,I)=OUT(J,24,I)+A(K,10,I)*B(J,3,K)
            OUT(J,25,I)=OUT(J,25,I)+A(K,5,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,1,K)
            OUT(J,26,I)=OUT(J,26,I)+A(K,6,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,2,K)+A(K,12,I)*B(J,1,K)
            OUT(J,27,I)=OUT(J,27,I)+A(K,7,I)*B(J,4,K)+A(K,12,I)
     $       *B(J,2,K)
            OUT(J,28,I)=OUT(J,28,I)+A(K,8,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,1,K)
            OUT(J,29,I)=OUT(J,29,I)+A(K,9,I)*B(J,4,K)+A(K,12,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,2,K)
            OUT(J,30,I)=OUT(J,30,I)+A(K,10,I)*B(J,4,K)+A(K,13,I)
     $       *B(J,3,K)
            OUT(J,31,I)=OUT(J,31,I)+A(K,11,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,1,K)
            OUT(J,32,I)=OUT(J,32,I)+A(K,12,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,2,K)
            OUT(J,33,I)=OUT(J,33,I)+A(K,13,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,3,K)
            OUT(J,34,I)=OUT(J,34,I)+A(K,14,I)*B(J,4,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE UPDATE_WL_0_2(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*16 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,14
            OUT(J,K,I)=(0.0D0,0.0D0)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,0,I)*B(J,5,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,0,I)*B(J,6,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,0,I)*B(J,7,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,0,I)*B(J,8,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,0,I)*B(J,9,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,0,I)*B(J,10,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,0,I)*B(J,11,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,0,I)*B(J,12,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,0,I)*B(J,13,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,0,I)*B(J,14,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_UPDATE_WL_0_2(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*32 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,14
            OUT(J,K,I)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,0,I)*B(J,5,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,0,I)*B(J,6,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,0,I)*B(J,7,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,0,I)*B(J,8,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,0,I)*B(J,9,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,0,I)*B(J,10,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,0,I)*B(J,11,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,0,I)*B(J,12,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,0,I)*B(J,13,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,0,I)*B(J,14,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE UPDATE_WL_2_0(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*16 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,14
            OUT(J,K,I)=(0.0D0,0.0D0)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,5,I)*B(J,0,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,6,I)*B(J,0,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,7,I)*B(J,0,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,8,I)*B(J,0,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,9,I)*B(J,0,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,10,I)*B(J,0,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,11,I)*B(J,0,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,12,I)*B(J,0,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,13,I)*B(J,0,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,14,I)*B(J,0,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_UPDATE_WL_2_0(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*32 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,14
            OUT(J,K,I)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,5,I)*B(J,0,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,6,I)*B(J,0,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,7,I)*B(J,0,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,8,I)*B(J,0,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,9,I)*B(J,0,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,10,I)*B(J,0,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,11,I)*B(J,0,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,12,I)*B(J,0,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,13,I)*B(J,0,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,14,I)*B(J,0,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE UPDATE_WL_2_2(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*16 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,69
            OUT(J,K,I)=(0.0D0,0.0D0)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,0,I)*B(J,5,K)+A(K,1,I)*B(J,1,K)
     $       +A(K,5,I)*B(J,0,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,0,I)*B(J,6,K)+A(K,1,I)*B(J,2,K)
     $       +A(K,2,I)*B(J,1,K)+A(K,6,I)*B(J,0,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,0,I)*B(J,7,K)+A(K,2,I)*B(J,2,K)
     $       +A(K,7,I)*B(J,0,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,0,I)*B(J,8,K)+A(K,1,I)*B(J,3,K)
     $       +A(K,3,I)*B(J,1,K)+A(K,8,I)*B(J,0,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,0,I)*B(J,9,K)+A(K,2,I)*B(J,3,K)
     $       +A(K,3,I)*B(J,2,K)+A(K,9,I)*B(J,0,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,0,I)*B(J,10,K)+A(K,3,I)
     $       *B(J,3,K)+A(K,10,I)*B(J,0,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,0,I)*B(J,11,K)+A(K,1,I)
     $       *B(J,4,K)+A(K,4,I)*B(J,1,K)+A(K,11,I)*B(J,0,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,0,I)*B(J,12,K)+A(K,2,I)
     $       *B(J,4,K)+A(K,4,I)*B(J,2,K)+A(K,12,I)*B(J,0,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,0,I)*B(J,13,K)+A(K,3,I)
     $       *B(J,4,K)+A(K,4,I)*B(J,3,K)+A(K,13,I)*B(J,0,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,0,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,4,K)+A(K,14,I)*B(J,0,K)
            OUT(J,15,I)=OUT(J,15,I)+A(K,1,I)*B(J,5,K)+A(K,5,I)*B(J,1,K)
            OUT(J,16,I)=OUT(J,16,I)+A(K,1,I)*B(J,6,K)+A(K,2,I)
     $       *B(J,5,K)+A(K,5,I)*B(J,2,K)+A(K,6,I)*B(J,1,K)
            OUT(J,17,I)=OUT(J,17,I)+A(K,1,I)*B(J,7,K)+A(K,2,I)
     $       *B(J,6,K)+A(K,6,I)*B(J,2,K)+A(K,7,I)*B(J,1,K)
            OUT(J,18,I)=OUT(J,18,I)+A(K,2,I)*B(J,7,K)+A(K,7,I)*B(J,2,K)
            OUT(J,19,I)=OUT(J,19,I)+A(K,1,I)*B(J,8,K)+A(K,3,I)
     $       *B(J,5,K)+A(K,5,I)*B(J,3,K)+A(K,8,I)*B(J,1,K)
            OUT(J,20,I)=OUT(J,20,I)+A(K,1,I)*B(J,9,K)+A(K,2,I)
     $       *B(J,8,K)+A(K,3,I)*B(J,6,K)+A(K,6,I)*B(J,3,K)+A(K,8,I)
     $       *B(J,2,K)+A(K,9,I)*B(J,1,K)
            OUT(J,21,I)=OUT(J,21,I)+A(K,2,I)*B(J,9,K)+A(K,3,I)
     $       *B(J,7,K)+A(K,7,I)*B(J,3,K)+A(K,9,I)*B(J,2,K)
            OUT(J,22,I)=OUT(J,22,I)+A(K,1,I)*B(J,10,K)+A(K,3,I)
     $       *B(J,8,K)+A(K,8,I)*B(J,3,K)+A(K,10,I)*B(J,1,K)
            OUT(J,23,I)=OUT(J,23,I)+A(K,2,I)*B(J,10,K)+A(K,3,I)
     $       *B(J,9,K)+A(K,9,I)*B(J,3,K)+A(K,10,I)*B(J,2,K)
            OUT(J,24,I)=OUT(J,24,I)+A(K,3,I)*B(J,10,K)+A(K,10,I)
     $       *B(J,3,K)
            OUT(J,25,I)=OUT(J,25,I)+A(K,1,I)*B(J,11,K)+A(K,4,I)
     $       *B(J,5,K)+A(K,5,I)*B(J,4,K)+A(K,11,I)*B(J,1,K)
            OUT(J,26,I)=OUT(J,26,I)+A(K,1,I)*B(J,12,K)+A(K,2,I)
     $       *B(J,11,K)+A(K,4,I)*B(J,6,K)+A(K,6,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,2,K)+A(K,12,I)*B(J,1,K)
            OUT(J,27,I)=OUT(J,27,I)+A(K,2,I)*B(J,12,K)+A(K,4,I)
     $       *B(J,7,K)+A(K,7,I)*B(J,4,K)+A(K,12,I)*B(J,2,K)
            OUT(J,28,I)=OUT(J,28,I)+A(K,1,I)*B(J,13,K)+A(K,3,I)
     $       *B(J,11,K)+A(K,4,I)*B(J,8,K)+A(K,8,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,1,K)
            OUT(J,29,I)=OUT(J,29,I)+A(K,2,I)*B(J,13,K)+A(K,3,I)
     $       *B(J,12,K)+A(K,4,I)*B(J,9,K)+A(K,9,I)*B(J,4,K)+A(K,12,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,2,K)
            OUT(J,30,I)=OUT(J,30,I)+A(K,3,I)*B(J,13,K)+A(K,4,I)
     $       *B(J,10,K)+A(K,10,I)*B(J,4,K)+A(K,13,I)*B(J,3,K)
            OUT(J,31,I)=OUT(J,31,I)+A(K,1,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,11,K)+A(K,11,I)*B(J,4,K)+A(K,14,I)*B(J,1,K)
            OUT(J,32,I)=OUT(J,32,I)+A(K,2,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,12,K)+A(K,12,I)*B(J,4,K)+A(K,14,I)*B(J,2,K)
            OUT(J,33,I)=OUT(J,33,I)+A(K,3,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,13,K)+A(K,13,I)*B(J,4,K)+A(K,14,I)*B(J,3,K)
            OUT(J,34,I)=OUT(J,34,I)+A(K,4,I)*B(J,14,K)+A(K,14,I)
     $       *B(J,4,K)
            OUT(J,35,I)=OUT(J,35,I)+A(K,5,I)*B(J,5,K)
            OUT(J,36,I)=OUT(J,36,I)+A(K,5,I)*B(J,6,K)+A(K,6,I)*B(J,5,K)
            OUT(J,37,I)=OUT(J,37,I)+A(K,5,I)*B(J,7,K)+A(K,6,I)
     $       *B(J,6,K)+A(K,7,I)*B(J,5,K)
            OUT(J,38,I)=OUT(J,38,I)+A(K,6,I)*B(J,7,K)+A(K,7,I)*B(J,6,K)
            OUT(J,39,I)=OUT(J,39,I)+A(K,7,I)*B(J,7,K)
            OUT(J,40,I)=OUT(J,40,I)+A(K,5,I)*B(J,8,K)+A(K,8,I)*B(J,5,K)
            OUT(J,41,I)=OUT(J,41,I)+A(K,5,I)*B(J,9,K)+A(K,6,I)
     $       *B(J,8,K)+A(K,8,I)*B(J,6,K)+A(K,9,I)*B(J,5,K)
            OUT(J,42,I)=OUT(J,42,I)+A(K,6,I)*B(J,9,K)+A(K,7,I)
     $       *B(J,8,K)+A(K,8,I)*B(J,7,K)+A(K,9,I)*B(J,6,K)
            OUT(J,43,I)=OUT(J,43,I)+A(K,7,I)*B(J,9,K)+A(K,9,I)*B(J,7,K)
            OUT(J,44,I)=OUT(J,44,I)+A(K,5,I)*B(J,10,K)+A(K,8,I)
     $       *B(J,8,K)+A(K,10,I)*B(J,5,K)
            OUT(J,45,I)=OUT(J,45,I)+A(K,6,I)*B(J,10,K)+A(K,8,I)
     $       *B(J,9,K)+A(K,9,I)*B(J,8,K)+A(K,10,I)*B(J,6,K)
            OUT(J,46,I)=OUT(J,46,I)+A(K,7,I)*B(J,10,K)+A(K,9,I)
     $       *B(J,9,K)+A(K,10,I)*B(J,7,K)
            OUT(J,47,I)=OUT(J,47,I)+A(K,8,I)*B(J,10,K)+A(K,10,I)
     $       *B(J,8,K)
            OUT(J,48,I)=OUT(J,48,I)+A(K,9,I)*B(J,10,K)+A(K,10,I)
     $       *B(J,9,K)
            OUT(J,49,I)=OUT(J,49,I)+A(K,10,I)*B(J,10,K)
            OUT(J,50,I)=OUT(J,50,I)+A(K,5,I)*B(J,11,K)+A(K,11,I)
     $       *B(J,5,K)
            OUT(J,51,I)=OUT(J,51,I)+A(K,5,I)*B(J,12,K)+A(K,6,I)
     $       *B(J,11,K)+A(K,11,I)*B(J,6,K)+A(K,12,I)*B(J,5,K)
            OUT(J,52,I)=OUT(J,52,I)+A(K,6,I)*B(J,12,K)+A(K,7,I)
     $       *B(J,11,K)+A(K,11,I)*B(J,7,K)+A(K,12,I)*B(J,6,K)
            OUT(J,53,I)=OUT(J,53,I)+A(K,7,I)*B(J,12,K)+A(K,12,I)
     $       *B(J,7,K)
            OUT(J,54,I)=OUT(J,54,I)+A(K,5,I)*B(J,13,K)+A(K,8,I)
     $       *B(J,11,K)+A(K,11,I)*B(J,8,K)+A(K,13,I)*B(J,5,K)
            OUT(J,55,I)=OUT(J,55,I)+A(K,6,I)*B(J,13,K)+A(K,8,I)
     $       *B(J,12,K)+A(K,9,I)*B(J,11,K)+A(K,11,I)*B(J,9,K)
     $       +A(K,12,I)*B(J,8,K)+A(K,13,I)*B(J,6,K)
            OUT(J,56,I)=OUT(J,56,I)+A(K,7,I)*B(J,13,K)+A(K,9,I)
     $       *B(J,12,K)+A(K,12,I)*B(J,9,K)+A(K,13,I)*B(J,7,K)
            OUT(J,57,I)=OUT(J,57,I)+A(K,8,I)*B(J,13,K)+A(K,10,I)
     $       *B(J,11,K)+A(K,11,I)*B(J,10,K)+A(K,13,I)*B(J,8,K)
            OUT(J,58,I)=OUT(J,58,I)+A(K,9,I)*B(J,13,K)+A(K,10,I)
     $       *B(J,12,K)+A(K,12,I)*B(J,10,K)+A(K,13,I)*B(J,9,K)
            OUT(J,59,I)=OUT(J,59,I)+A(K,10,I)*B(J,13,K)+A(K,13,I)
     $       *B(J,10,K)
            OUT(J,60,I)=OUT(J,60,I)+A(K,5,I)*B(J,14,K)+A(K,11,I)
     $       *B(J,11,K)+A(K,14,I)*B(J,5,K)
            OUT(J,61,I)=OUT(J,61,I)+A(K,6,I)*B(J,14,K)+A(K,11,I)
     $       *B(J,12,K)+A(K,12,I)*B(J,11,K)+A(K,14,I)*B(J,6,K)
            OUT(J,62,I)=OUT(J,62,I)+A(K,7,I)*B(J,14,K)+A(K,12,I)
     $       *B(J,12,K)+A(K,14,I)*B(J,7,K)
            OUT(J,63,I)=OUT(J,63,I)+A(K,8,I)*B(J,14,K)+A(K,11,I)
     $       *B(J,13,K)+A(K,13,I)*B(J,11,K)+A(K,14,I)*B(J,8,K)
            OUT(J,64,I)=OUT(J,64,I)+A(K,9,I)*B(J,14,K)+A(K,12,I)
     $       *B(J,13,K)+A(K,13,I)*B(J,12,K)+A(K,14,I)*B(J,9,K)
            OUT(J,65,I)=OUT(J,65,I)+A(K,10,I)*B(J,14,K)+A(K,13,I)
     $       *B(J,13,K)+A(K,14,I)*B(J,10,K)
            OUT(J,66,I)=OUT(J,66,I)+A(K,11,I)*B(J,14,K)+A(K,14,I)
     $       *B(J,11,K)
            OUT(J,67,I)=OUT(J,67,I)+A(K,12,I)*B(J,14,K)+A(K,14,I)
     $       *B(J,12,K)
            OUT(J,68,I)=OUT(J,68,I)+A(K,13,I)*B(J,14,K)+A(K,14,I)
     $       *B(J,13,K)
            OUT(J,69,I)=OUT(J,69,I)+A(K,14,I)*B(J,14,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_UPDATE_WL_2_2(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*32 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,69
            OUT(J,K,I)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,0,I)*B(J,5,K)+A(K,1,I)*B(J,1,K)
     $       +A(K,5,I)*B(J,0,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,0,I)*B(J,6,K)+A(K,1,I)*B(J,2,K)
     $       +A(K,2,I)*B(J,1,K)+A(K,6,I)*B(J,0,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,0,I)*B(J,7,K)+A(K,2,I)*B(J,2,K)
     $       +A(K,7,I)*B(J,0,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,0,I)*B(J,8,K)+A(K,1,I)*B(J,3,K)
     $       +A(K,3,I)*B(J,1,K)+A(K,8,I)*B(J,0,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,0,I)*B(J,9,K)+A(K,2,I)*B(J,3,K)
     $       +A(K,3,I)*B(J,2,K)+A(K,9,I)*B(J,0,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,0,I)*B(J,10,K)+A(K,3,I)
     $       *B(J,3,K)+A(K,10,I)*B(J,0,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,0,I)*B(J,11,K)+A(K,1,I)
     $       *B(J,4,K)+A(K,4,I)*B(J,1,K)+A(K,11,I)*B(J,0,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,0,I)*B(J,12,K)+A(K,2,I)
     $       *B(J,4,K)+A(K,4,I)*B(J,2,K)+A(K,12,I)*B(J,0,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,0,I)*B(J,13,K)+A(K,3,I)
     $       *B(J,4,K)+A(K,4,I)*B(J,3,K)+A(K,13,I)*B(J,0,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,0,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,4,K)+A(K,14,I)*B(J,0,K)
            OUT(J,15,I)=OUT(J,15,I)+A(K,1,I)*B(J,5,K)+A(K,5,I)*B(J,1,K)
            OUT(J,16,I)=OUT(J,16,I)+A(K,1,I)*B(J,6,K)+A(K,2,I)
     $       *B(J,5,K)+A(K,5,I)*B(J,2,K)+A(K,6,I)*B(J,1,K)
            OUT(J,17,I)=OUT(J,17,I)+A(K,1,I)*B(J,7,K)+A(K,2,I)
     $       *B(J,6,K)+A(K,6,I)*B(J,2,K)+A(K,7,I)*B(J,1,K)
            OUT(J,18,I)=OUT(J,18,I)+A(K,2,I)*B(J,7,K)+A(K,7,I)*B(J,2,K)
            OUT(J,19,I)=OUT(J,19,I)+A(K,1,I)*B(J,8,K)+A(K,3,I)
     $       *B(J,5,K)+A(K,5,I)*B(J,3,K)+A(K,8,I)*B(J,1,K)
            OUT(J,20,I)=OUT(J,20,I)+A(K,1,I)*B(J,9,K)+A(K,2,I)
     $       *B(J,8,K)+A(K,3,I)*B(J,6,K)+A(K,6,I)*B(J,3,K)+A(K,8,I)
     $       *B(J,2,K)+A(K,9,I)*B(J,1,K)
            OUT(J,21,I)=OUT(J,21,I)+A(K,2,I)*B(J,9,K)+A(K,3,I)
     $       *B(J,7,K)+A(K,7,I)*B(J,3,K)+A(K,9,I)*B(J,2,K)
            OUT(J,22,I)=OUT(J,22,I)+A(K,1,I)*B(J,10,K)+A(K,3,I)
     $       *B(J,8,K)+A(K,8,I)*B(J,3,K)+A(K,10,I)*B(J,1,K)
            OUT(J,23,I)=OUT(J,23,I)+A(K,2,I)*B(J,10,K)+A(K,3,I)
     $       *B(J,9,K)+A(K,9,I)*B(J,3,K)+A(K,10,I)*B(J,2,K)
            OUT(J,24,I)=OUT(J,24,I)+A(K,3,I)*B(J,10,K)+A(K,10,I)
     $       *B(J,3,K)
            OUT(J,25,I)=OUT(J,25,I)+A(K,1,I)*B(J,11,K)+A(K,4,I)
     $       *B(J,5,K)+A(K,5,I)*B(J,4,K)+A(K,11,I)*B(J,1,K)
            OUT(J,26,I)=OUT(J,26,I)+A(K,1,I)*B(J,12,K)+A(K,2,I)
     $       *B(J,11,K)+A(K,4,I)*B(J,6,K)+A(K,6,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,2,K)+A(K,12,I)*B(J,1,K)
            OUT(J,27,I)=OUT(J,27,I)+A(K,2,I)*B(J,12,K)+A(K,4,I)
     $       *B(J,7,K)+A(K,7,I)*B(J,4,K)+A(K,12,I)*B(J,2,K)
            OUT(J,28,I)=OUT(J,28,I)+A(K,1,I)*B(J,13,K)+A(K,3,I)
     $       *B(J,11,K)+A(K,4,I)*B(J,8,K)+A(K,8,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,1,K)
            OUT(J,29,I)=OUT(J,29,I)+A(K,2,I)*B(J,13,K)+A(K,3,I)
     $       *B(J,12,K)+A(K,4,I)*B(J,9,K)+A(K,9,I)*B(J,4,K)+A(K,12,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,2,K)
            OUT(J,30,I)=OUT(J,30,I)+A(K,3,I)*B(J,13,K)+A(K,4,I)
     $       *B(J,10,K)+A(K,10,I)*B(J,4,K)+A(K,13,I)*B(J,3,K)
            OUT(J,31,I)=OUT(J,31,I)+A(K,1,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,11,K)+A(K,11,I)*B(J,4,K)+A(K,14,I)*B(J,1,K)
            OUT(J,32,I)=OUT(J,32,I)+A(K,2,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,12,K)+A(K,12,I)*B(J,4,K)+A(K,14,I)*B(J,2,K)
            OUT(J,33,I)=OUT(J,33,I)+A(K,3,I)*B(J,14,K)+A(K,4,I)
     $       *B(J,13,K)+A(K,13,I)*B(J,4,K)+A(K,14,I)*B(J,3,K)
            OUT(J,34,I)=OUT(J,34,I)+A(K,4,I)*B(J,14,K)+A(K,14,I)
     $       *B(J,4,K)
            OUT(J,35,I)=OUT(J,35,I)+A(K,5,I)*B(J,5,K)
            OUT(J,36,I)=OUT(J,36,I)+A(K,5,I)*B(J,6,K)+A(K,6,I)*B(J,5,K)
            OUT(J,37,I)=OUT(J,37,I)+A(K,5,I)*B(J,7,K)+A(K,6,I)
     $       *B(J,6,K)+A(K,7,I)*B(J,5,K)
            OUT(J,38,I)=OUT(J,38,I)+A(K,6,I)*B(J,7,K)+A(K,7,I)*B(J,6,K)
            OUT(J,39,I)=OUT(J,39,I)+A(K,7,I)*B(J,7,K)
            OUT(J,40,I)=OUT(J,40,I)+A(K,5,I)*B(J,8,K)+A(K,8,I)*B(J,5,K)
            OUT(J,41,I)=OUT(J,41,I)+A(K,5,I)*B(J,9,K)+A(K,6,I)
     $       *B(J,8,K)+A(K,8,I)*B(J,6,K)+A(K,9,I)*B(J,5,K)
            OUT(J,42,I)=OUT(J,42,I)+A(K,6,I)*B(J,9,K)+A(K,7,I)
     $       *B(J,8,K)+A(K,8,I)*B(J,7,K)+A(K,9,I)*B(J,6,K)
            OUT(J,43,I)=OUT(J,43,I)+A(K,7,I)*B(J,9,K)+A(K,9,I)*B(J,7,K)
            OUT(J,44,I)=OUT(J,44,I)+A(K,5,I)*B(J,10,K)+A(K,8,I)
     $       *B(J,8,K)+A(K,10,I)*B(J,5,K)
            OUT(J,45,I)=OUT(J,45,I)+A(K,6,I)*B(J,10,K)+A(K,8,I)
     $       *B(J,9,K)+A(K,9,I)*B(J,8,K)+A(K,10,I)*B(J,6,K)
            OUT(J,46,I)=OUT(J,46,I)+A(K,7,I)*B(J,10,K)+A(K,9,I)
     $       *B(J,9,K)+A(K,10,I)*B(J,7,K)
            OUT(J,47,I)=OUT(J,47,I)+A(K,8,I)*B(J,10,K)+A(K,10,I)
     $       *B(J,8,K)
            OUT(J,48,I)=OUT(J,48,I)+A(K,9,I)*B(J,10,K)+A(K,10,I)
     $       *B(J,9,K)
            OUT(J,49,I)=OUT(J,49,I)+A(K,10,I)*B(J,10,K)
            OUT(J,50,I)=OUT(J,50,I)+A(K,5,I)*B(J,11,K)+A(K,11,I)
     $       *B(J,5,K)
            OUT(J,51,I)=OUT(J,51,I)+A(K,5,I)*B(J,12,K)+A(K,6,I)
     $       *B(J,11,K)+A(K,11,I)*B(J,6,K)+A(K,12,I)*B(J,5,K)
            OUT(J,52,I)=OUT(J,52,I)+A(K,6,I)*B(J,12,K)+A(K,7,I)
     $       *B(J,11,K)+A(K,11,I)*B(J,7,K)+A(K,12,I)*B(J,6,K)
            OUT(J,53,I)=OUT(J,53,I)+A(K,7,I)*B(J,12,K)+A(K,12,I)
     $       *B(J,7,K)
            OUT(J,54,I)=OUT(J,54,I)+A(K,5,I)*B(J,13,K)+A(K,8,I)
     $       *B(J,11,K)+A(K,11,I)*B(J,8,K)+A(K,13,I)*B(J,5,K)
            OUT(J,55,I)=OUT(J,55,I)+A(K,6,I)*B(J,13,K)+A(K,8,I)
     $       *B(J,12,K)+A(K,9,I)*B(J,11,K)+A(K,11,I)*B(J,9,K)
     $       +A(K,12,I)*B(J,8,K)+A(K,13,I)*B(J,6,K)
            OUT(J,56,I)=OUT(J,56,I)+A(K,7,I)*B(J,13,K)+A(K,9,I)
     $       *B(J,12,K)+A(K,12,I)*B(J,9,K)+A(K,13,I)*B(J,7,K)
            OUT(J,57,I)=OUT(J,57,I)+A(K,8,I)*B(J,13,K)+A(K,10,I)
     $       *B(J,11,K)+A(K,11,I)*B(J,10,K)+A(K,13,I)*B(J,8,K)
            OUT(J,58,I)=OUT(J,58,I)+A(K,9,I)*B(J,13,K)+A(K,10,I)
     $       *B(J,12,K)+A(K,12,I)*B(J,10,K)+A(K,13,I)*B(J,9,K)
            OUT(J,59,I)=OUT(J,59,I)+A(K,10,I)*B(J,13,K)+A(K,13,I)
     $       *B(J,10,K)
            OUT(J,60,I)=OUT(J,60,I)+A(K,5,I)*B(J,14,K)+A(K,11,I)
     $       *B(J,11,K)+A(K,14,I)*B(J,5,K)
            OUT(J,61,I)=OUT(J,61,I)+A(K,6,I)*B(J,14,K)+A(K,11,I)
     $       *B(J,12,K)+A(K,12,I)*B(J,11,K)+A(K,14,I)*B(J,6,K)
            OUT(J,62,I)=OUT(J,62,I)+A(K,7,I)*B(J,14,K)+A(K,12,I)
     $       *B(J,12,K)+A(K,14,I)*B(J,7,K)
            OUT(J,63,I)=OUT(J,63,I)+A(K,8,I)*B(J,14,K)+A(K,11,I)
     $       *B(J,13,K)+A(K,13,I)*B(J,11,K)+A(K,14,I)*B(J,8,K)
            OUT(J,64,I)=OUT(J,64,I)+A(K,9,I)*B(J,14,K)+A(K,12,I)
     $       *B(J,13,K)+A(K,13,I)*B(J,12,K)+A(K,14,I)*B(J,9,K)
            OUT(J,65,I)=OUT(J,65,I)+A(K,10,I)*B(J,14,K)+A(K,13,I)
     $       *B(J,13,K)+A(K,14,I)*B(J,10,K)
            OUT(J,66,I)=OUT(J,66,I)+A(K,11,I)*B(J,14,K)+A(K,14,I)
     $       *B(J,11,K)
            OUT(J,67,I)=OUT(J,67,I)+A(K,12,I)*B(J,14,K)+A(K,14,I)
     $       *B(J,12,K)
            OUT(J,68,I)=OUT(J,68,I)+A(K,13,I)*B(J,14,K)+A(K,14,I)
     $       *B(J,13,K)
            OUT(J,69,I)=OUT(J,69,I)+A(K,14,I)*B(J,14,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE UPDATE_WL_1_0(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*16 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,4
            OUT(J,K,I)=(0.0D0,0.0D0)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,4,I)*B(J,0,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_UPDATE_WL_1_0(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*32 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,4
            OUT(J,K,I)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,4,I)*B(J,0,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE UPDATE_WL_4_1(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*16 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,125
            OUT(J,K,I)=(0.0D0,0.0D0)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,1,I)*B(J,1,K)+A(K,5,I)*B(J,0,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,1,I)*B(J,2,K)+A(K,2,I)*B(J,1,K)
     $       +A(K,6,I)*B(J,0,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,2,I)*B(J,2,K)+A(K,7,I)*B(J,0,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,1,I)*B(J,3,K)+A(K,3,I)*B(J,1,K)
     $       +A(K,8,I)*B(J,0,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,2,I)*B(J,3,K)+A(K,3,I)*B(J,2,K)
     $       +A(K,9,I)*B(J,0,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,3,I)*B(J,3,K)+A(K,10,I)
     $       *B(J,0,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,1,I)*B(J,4,K)+A(K,4,I)
     $       *B(J,1,K)+A(K,11,I)*B(J,0,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,2,I)*B(J,4,K)+A(K,4,I)
     $       *B(J,2,K)+A(K,12,I)*B(J,0,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,3,I)*B(J,4,K)+A(K,4,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,0,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,4,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,0,K)
            OUT(J,15,I)=OUT(J,15,I)+A(K,5,I)*B(J,1,K)+A(K,15,I)
     $       *B(J,0,K)
            OUT(J,16,I)=OUT(J,16,I)+A(K,5,I)*B(J,2,K)+A(K,6,I)
     $       *B(J,1,K)+A(K,16,I)*B(J,0,K)
            OUT(J,17,I)=OUT(J,17,I)+A(K,6,I)*B(J,2,K)+A(K,7,I)
     $       *B(J,1,K)+A(K,17,I)*B(J,0,K)
            OUT(J,18,I)=OUT(J,18,I)+A(K,7,I)*B(J,2,K)+A(K,18,I)
     $       *B(J,0,K)
            OUT(J,19,I)=OUT(J,19,I)+A(K,5,I)*B(J,3,K)+A(K,8,I)
     $       *B(J,1,K)+A(K,19,I)*B(J,0,K)
            OUT(J,20,I)=OUT(J,20,I)+A(K,6,I)*B(J,3,K)+A(K,8,I)
     $       *B(J,2,K)+A(K,9,I)*B(J,1,K)+A(K,20,I)*B(J,0,K)
            OUT(J,21,I)=OUT(J,21,I)+A(K,7,I)*B(J,3,K)+A(K,9,I)
     $       *B(J,2,K)+A(K,21,I)*B(J,0,K)
            OUT(J,22,I)=OUT(J,22,I)+A(K,8,I)*B(J,3,K)+A(K,10,I)
     $       *B(J,1,K)+A(K,22,I)*B(J,0,K)
            OUT(J,23,I)=OUT(J,23,I)+A(K,9,I)*B(J,3,K)+A(K,10,I)
     $       *B(J,2,K)+A(K,23,I)*B(J,0,K)
            OUT(J,24,I)=OUT(J,24,I)+A(K,10,I)*B(J,3,K)+A(K,24,I)
     $       *B(J,0,K)
            OUT(J,25,I)=OUT(J,25,I)+A(K,5,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,1,K)+A(K,25,I)*B(J,0,K)
            OUT(J,26,I)=OUT(J,26,I)+A(K,6,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,2,K)+A(K,12,I)*B(J,1,K)+A(K,26,I)*B(J,0,K)
            OUT(J,27,I)=OUT(J,27,I)+A(K,7,I)*B(J,4,K)+A(K,12,I)
     $       *B(J,2,K)+A(K,27,I)*B(J,0,K)
            OUT(J,28,I)=OUT(J,28,I)+A(K,8,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,1,K)+A(K,28,I)*B(J,0,K)
            OUT(J,29,I)=OUT(J,29,I)+A(K,9,I)*B(J,4,K)+A(K,12,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,2,K)+A(K,29,I)*B(J,0,K)
            OUT(J,30,I)=OUT(J,30,I)+A(K,10,I)*B(J,4,K)+A(K,13,I)
     $       *B(J,3,K)+A(K,30,I)*B(J,0,K)
            OUT(J,31,I)=OUT(J,31,I)+A(K,11,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,1,K)+A(K,31,I)*B(J,0,K)
            OUT(J,32,I)=OUT(J,32,I)+A(K,12,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,2,K)+A(K,32,I)*B(J,0,K)
            OUT(J,33,I)=OUT(J,33,I)+A(K,13,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,3,K)+A(K,33,I)*B(J,0,K)
            OUT(J,34,I)=OUT(J,34,I)+A(K,14,I)*B(J,4,K)+A(K,34,I)
     $       *B(J,0,K)
            OUT(J,35,I)=OUT(J,35,I)+A(K,15,I)*B(J,1,K)+A(K,35,I)
     $       *B(J,0,K)
            OUT(J,36,I)=OUT(J,36,I)+A(K,15,I)*B(J,2,K)+A(K,16,I)
     $       *B(J,1,K)+A(K,36,I)*B(J,0,K)
            OUT(J,37,I)=OUT(J,37,I)+A(K,16,I)*B(J,2,K)+A(K,17,I)
     $       *B(J,1,K)+A(K,37,I)*B(J,0,K)
            OUT(J,38,I)=OUT(J,38,I)+A(K,17,I)*B(J,2,K)+A(K,18,I)
     $       *B(J,1,K)+A(K,38,I)*B(J,0,K)
            OUT(J,39,I)=OUT(J,39,I)+A(K,18,I)*B(J,2,K)+A(K,39,I)
     $       *B(J,0,K)
            OUT(J,40,I)=OUT(J,40,I)+A(K,15,I)*B(J,3,K)+A(K,19,I)
     $       *B(J,1,K)+A(K,40,I)*B(J,0,K)
            OUT(J,41,I)=OUT(J,41,I)+A(K,16,I)*B(J,3,K)+A(K,19,I)
     $       *B(J,2,K)+A(K,20,I)*B(J,1,K)+A(K,41,I)*B(J,0,K)
            OUT(J,42,I)=OUT(J,42,I)+A(K,17,I)*B(J,3,K)+A(K,20,I)
     $       *B(J,2,K)+A(K,21,I)*B(J,1,K)+A(K,42,I)*B(J,0,K)
            OUT(J,43,I)=OUT(J,43,I)+A(K,18,I)*B(J,3,K)+A(K,21,I)
     $       *B(J,2,K)+A(K,43,I)*B(J,0,K)
            OUT(J,44,I)=OUT(J,44,I)+A(K,19,I)*B(J,3,K)+A(K,22,I)
     $       *B(J,1,K)+A(K,44,I)*B(J,0,K)
            OUT(J,45,I)=OUT(J,45,I)+A(K,20,I)*B(J,3,K)+A(K,22,I)
     $       *B(J,2,K)+A(K,23,I)*B(J,1,K)+A(K,45,I)*B(J,0,K)
            OUT(J,46,I)=OUT(J,46,I)+A(K,21,I)*B(J,3,K)+A(K,23,I)
     $       *B(J,2,K)+A(K,46,I)*B(J,0,K)
            OUT(J,47,I)=OUT(J,47,I)+A(K,22,I)*B(J,3,K)+A(K,24,I)
     $       *B(J,1,K)+A(K,47,I)*B(J,0,K)
            OUT(J,48,I)=OUT(J,48,I)+A(K,23,I)*B(J,3,K)+A(K,24,I)
     $       *B(J,2,K)+A(K,48,I)*B(J,0,K)
            OUT(J,49,I)=OUT(J,49,I)+A(K,24,I)*B(J,3,K)+A(K,49,I)
     $       *B(J,0,K)
            OUT(J,50,I)=OUT(J,50,I)+A(K,15,I)*B(J,4,K)+A(K,25,I)
     $       *B(J,1,K)+A(K,50,I)*B(J,0,K)
            OUT(J,51,I)=OUT(J,51,I)+A(K,16,I)*B(J,4,K)+A(K,25,I)
     $       *B(J,2,K)+A(K,26,I)*B(J,1,K)+A(K,51,I)*B(J,0,K)
            OUT(J,52,I)=OUT(J,52,I)+A(K,17,I)*B(J,4,K)+A(K,26,I)
     $       *B(J,2,K)+A(K,27,I)*B(J,1,K)+A(K,52,I)*B(J,0,K)
            OUT(J,53,I)=OUT(J,53,I)+A(K,18,I)*B(J,4,K)+A(K,27,I)
     $       *B(J,2,K)+A(K,53,I)*B(J,0,K)
            OUT(J,54,I)=OUT(J,54,I)+A(K,19,I)*B(J,4,K)+A(K,25,I)
     $       *B(J,3,K)+A(K,28,I)*B(J,1,K)+A(K,54,I)*B(J,0,K)
            OUT(J,55,I)=OUT(J,55,I)+A(K,20,I)*B(J,4,K)+A(K,26,I)
     $       *B(J,3,K)+A(K,28,I)*B(J,2,K)+A(K,29,I)*B(J,1,K)+A(K,55,I)
     $       *B(J,0,K)
            OUT(J,56,I)=OUT(J,56,I)+A(K,21,I)*B(J,4,K)+A(K,27,I)
     $       *B(J,3,K)+A(K,29,I)*B(J,2,K)+A(K,56,I)*B(J,0,K)
            OUT(J,57,I)=OUT(J,57,I)+A(K,22,I)*B(J,4,K)+A(K,28,I)
     $       *B(J,3,K)+A(K,30,I)*B(J,1,K)+A(K,57,I)*B(J,0,K)
            OUT(J,58,I)=OUT(J,58,I)+A(K,23,I)*B(J,4,K)+A(K,29,I)
     $       *B(J,3,K)+A(K,30,I)*B(J,2,K)+A(K,58,I)*B(J,0,K)
            OUT(J,59,I)=OUT(J,59,I)+A(K,24,I)*B(J,4,K)+A(K,30,I)
     $       *B(J,3,K)+A(K,59,I)*B(J,0,K)
            OUT(J,60,I)=OUT(J,60,I)+A(K,25,I)*B(J,4,K)+A(K,31,I)
     $       *B(J,1,K)+A(K,60,I)*B(J,0,K)
            OUT(J,61,I)=OUT(J,61,I)+A(K,26,I)*B(J,4,K)+A(K,31,I)
     $       *B(J,2,K)+A(K,32,I)*B(J,1,K)+A(K,61,I)*B(J,0,K)
            OUT(J,62,I)=OUT(J,62,I)+A(K,27,I)*B(J,4,K)+A(K,32,I)
     $       *B(J,2,K)+A(K,62,I)*B(J,0,K)
            OUT(J,63,I)=OUT(J,63,I)+A(K,28,I)*B(J,4,K)+A(K,31,I)
     $       *B(J,3,K)+A(K,33,I)*B(J,1,K)+A(K,63,I)*B(J,0,K)
            OUT(J,64,I)=OUT(J,64,I)+A(K,29,I)*B(J,4,K)+A(K,32,I)
     $       *B(J,3,K)+A(K,33,I)*B(J,2,K)+A(K,64,I)*B(J,0,K)
            OUT(J,65,I)=OUT(J,65,I)+A(K,30,I)*B(J,4,K)+A(K,33,I)
     $       *B(J,3,K)+A(K,65,I)*B(J,0,K)
            OUT(J,66,I)=OUT(J,66,I)+A(K,31,I)*B(J,4,K)+A(K,34,I)
     $       *B(J,1,K)+A(K,66,I)*B(J,0,K)
            OUT(J,67,I)=OUT(J,67,I)+A(K,32,I)*B(J,4,K)+A(K,34,I)
     $       *B(J,2,K)+A(K,67,I)*B(J,0,K)
            OUT(J,68,I)=OUT(J,68,I)+A(K,33,I)*B(J,4,K)+A(K,34,I)
     $       *B(J,3,K)+A(K,68,I)*B(J,0,K)
            OUT(J,69,I)=OUT(J,69,I)+A(K,34,I)*B(J,4,K)+A(K,69,I)
     $       *B(J,0,K)
            OUT(J,70,I)=OUT(J,70,I)+A(K,35,I)*B(J,1,K)
            OUT(J,71,I)=OUT(J,71,I)+A(K,35,I)*B(J,2,K)+A(K,36,I)
     $       *B(J,1,K)
            OUT(J,72,I)=OUT(J,72,I)+A(K,36,I)*B(J,2,K)+A(K,37,I)
     $       *B(J,1,K)
            OUT(J,73,I)=OUT(J,73,I)+A(K,37,I)*B(J,2,K)+A(K,38,I)
     $       *B(J,1,K)
            OUT(J,74,I)=OUT(J,74,I)+A(K,38,I)*B(J,2,K)+A(K,39,I)
     $       *B(J,1,K)
            OUT(J,75,I)=OUT(J,75,I)+A(K,39,I)*B(J,2,K)
            OUT(J,76,I)=OUT(J,76,I)+A(K,35,I)*B(J,3,K)+A(K,40,I)
     $       *B(J,1,K)
            OUT(J,77,I)=OUT(J,77,I)+A(K,36,I)*B(J,3,K)+A(K,40,I)
     $       *B(J,2,K)+A(K,41,I)*B(J,1,K)
            OUT(J,78,I)=OUT(J,78,I)+A(K,37,I)*B(J,3,K)+A(K,41,I)
     $       *B(J,2,K)+A(K,42,I)*B(J,1,K)
            OUT(J,79,I)=OUT(J,79,I)+A(K,38,I)*B(J,3,K)+A(K,42,I)
     $       *B(J,2,K)+A(K,43,I)*B(J,1,K)
            OUT(J,80,I)=OUT(J,80,I)+A(K,39,I)*B(J,3,K)+A(K,43,I)
     $       *B(J,2,K)
            OUT(J,81,I)=OUT(J,81,I)+A(K,40,I)*B(J,3,K)+A(K,44,I)
     $       *B(J,1,K)
            OUT(J,82,I)=OUT(J,82,I)+A(K,41,I)*B(J,3,K)+A(K,44,I)
     $       *B(J,2,K)+A(K,45,I)*B(J,1,K)
            OUT(J,83,I)=OUT(J,83,I)+A(K,42,I)*B(J,3,K)+A(K,45,I)
     $       *B(J,2,K)+A(K,46,I)*B(J,1,K)
            OUT(J,84,I)=OUT(J,84,I)+A(K,43,I)*B(J,3,K)+A(K,46,I)
     $       *B(J,2,K)
            OUT(J,85,I)=OUT(J,85,I)+A(K,44,I)*B(J,3,K)+A(K,47,I)
     $       *B(J,1,K)
            OUT(J,86,I)=OUT(J,86,I)+A(K,45,I)*B(J,3,K)+A(K,47,I)
     $       *B(J,2,K)+A(K,48,I)*B(J,1,K)
            OUT(J,87,I)=OUT(J,87,I)+A(K,46,I)*B(J,3,K)+A(K,48,I)
     $       *B(J,2,K)
            OUT(J,88,I)=OUT(J,88,I)+A(K,47,I)*B(J,3,K)+A(K,49,I)
     $       *B(J,1,K)
            OUT(J,89,I)=OUT(J,89,I)+A(K,48,I)*B(J,3,K)+A(K,49,I)
     $       *B(J,2,K)
            OUT(J,90,I)=OUT(J,90,I)+A(K,49,I)*B(J,3,K)
            OUT(J,91,I)=OUT(J,91,I)+A(K,35,I)*B(J,4,K)+A(K,50,I)
     $       *B(J,1,K)
            OUT(J,92,I)=OUT(J,92,I)+A(K,36,I)*B(J,4,K)+A(K,50,I)
     $       *B(J,2,K)+A(K,51,I)*B(J,1,K)
            OUT(J,93,I)=OUT(J,93,I)+A(K,37,I)*B(J,4,K)+A(K,51,I)
     $       *B(J,2,K)+A(K,52,I)*B(J,1,K)
            OUT(J,94,I)=OUT(J,94,I)+A(K,38,I)*B(J,4,K)+A(K,52,I)
     $       *B(J,2,K)+A(K,53,I)*B(J,1,K)
            OUT(J,95,I)=OUT(J,95,I)+A(K,39,I)*B(J,4,K)+A(K,53,I)
     $       *B(J,2,K)
            OUT(J,96,I)=OUT(J,96,I)+A(K,40,I)*B(J,4,K)+A(K,50,I)
     $       *B(J,3,K)+A(K,54,I)*B(J,1,K)
            OUT(J,97,I)=OUT(J,97,I)+A(K,41,I)*B(J,4,K)+A(K,51,I)
     $       *B(J,3,K)+A(K,54,I)*B(J,2,K)+A(K,55,I)*B(J,1,K)
            OUT(J,98,I)=OUT(J,98,I)+A(K,42,I)*B(J,4,K)+A(K,52,I)
     $       *B(J,3,K)+A(K,55,I)*B(J,2,K)+A(K,56,I)*B(J,1,K)
            OUT(J,99,I)=OUT(J,99,I)+A(K,43,I)*B(J,4,K)+A(K,53,I)
     $       *B(J,3,K)+A(K,56,I)*B(J,2,K)
            OUT(J,100,I)=OUT(J,100,I)+A(K,44,I)*B(J,4,K)+A(K,54,I)
     $       *B(J,3,K)+A(K,57,I)*B(J,1,K)
            OUT(J,101,I)=OUT(J,101,I)+A(K,45,I)*B(J,4,K)+A(K,55,I)
     $       *B(J,3,K)+A(K,57,I)*B(J,2,K)+A(K,58,I)*B(J,1,K)
            OUT(J,102,I)=OUT(J,102,I)+A(K,46,I)*B(J,4,K)+A(K,56,I)
     $       *B(J,3,K)+A(K,58,I)*B(J,2,K)
            OUT(J,103,I)=OUT(J,103,I)+A(K,47,I)*B(J,4,K)+A(K,57,I)
     $       *B(J,3,K)+A(K,59,I)*B(J,1,K)
            OUT(J,104,I)=OUT(J,104,I)+A(K,48,I)*B(J,4,K)+A(K,58,I)
     $       *B(J,3,K)+A(K,59,I)*B(J,2,K)
            OUT(J,105,I)=OUT(J,105,I)+A(K,49,I)*B(J,4,K)+A(K,59,I)
     $       *B(J,3,K)
            OUT(J,106,I)=OUT(J,106,I)+A(K,50,I)*B(J,4,K)+A(K,60,I)
     $       *B(J,1,K)
            OUT(J,107,I)=OUT(J,107,I)+A(K,51,I)*B(J,4,K)+A(K,60,I)
     $       *B(J,2,K)+A(K,61,I)*B(J,1,K)
            OUT(J,108,I)=OUT(J,108,I)+A(K,52,I)*B(J,4,K)+A(K,61,I)
     $       *B(J,2,K)+A(K,62,I)*B(J,1,K)
            OUT(J,109,I)=OUT(J,109,I)+A(K,53,I)*B(J,4,K)+A(K,62,I)
     $       *B(J,2,K)
            OUT(J,110,I)=OUT(J,110,I)+A(K,54,I)*B(J,4,K)+A(K,60,I)
     $       *B(J,3,K)+A(K,63,I)*B(J,1,K)
            OUT(J,111,I)=OUT(J,111,I)+A(K,55,I)*B(J,4,K)+A(K,61,I)
     $       *B(J,3,K)+A(K,63,I)*B(J,2,K)+A(K,64,I)*B(J,1,K)
            OUT(J,112,I)=OUT(J,112,I)+A(K,56,I)*B(J,4,K)+A(K,62,I)
     $       *B(J,3,K)+A(K,64,I)*B(J,2,K)
            OUT(J,113,I)=OUT(J,113,I)+A(K,57,I)*B(J,4,K)+A(K,63,I)
     $       *B(J,3,K)+A(K,65,I)*B(J,1,K)
            OUT(J,114,I)=OUT(J,114,I)+A(K,58,I)*B(J,4,K)+A(K,64,I)
     $       *B(J,3,K)+A(K,65,I)*B(J,2,K)
            OUT(J,115,I)=OUT(J,115,I)+A(K,59,I)*B(J,4,K)+A(K,65,I)
     $       *B(J,3,K)
            OUT(J,116,I)=OUT(J,116,I)+A(K,60,I)*B(J,4,K)+A(K,66,I)
     $       *B(J,1,K)
            OUT(J,117,I)=OUT(J,117,I)+A(K,61,I)*B(J,4,K)+A(K,66,I)
     $       *B(J,2,K)+A(K,67,I)*B(J,1,K)
            OUT(J,118,I)=OUT(J,118,I)+A(K,62,I)*B(J,4,K)+A(K,67,I)
     $       *B(J,2,K)
            OUT(J,119,I)=OUT(J,119,I)+A(K,63,I)*B(J,4,K)+A(K,66,I)
     $       *B(J,3,K)+A(K,68,I)*B(J,1,K)
            OUT(J,120,I)=OUT(J,120,I)+A(K,64,I)*B(J,4,K)+A(K,67,I)
     $       *B(J,3,K)+A(K,68,I)*B(J,2,K)
            OUT(J,121,I)=OUT(J,121,I)+A(K,65,I)*B(J,4,K)+A(K,68,I)
     $       *B(J,3,K)
            OUT(J,122,I)=OUT(J,122,I)+A(K,66,I)*B(J,4,K)+A(K,69,I)
     $       *B(J,1,K)
            OUT(J,123,I)=OUT(J,123,I)+A(K,67,I)*B(J,4,K)+A(K,69,I)
     $       *B(J,2,K)
            OUT(J,124,I)=OUT(J,124,I)+A(K,68,I)*B(J,4,K)+A(K,69,I)
     $       *B(J,3,K)
            OUT(J,125,I)=OUT(J,125,I)+A(K,69,I)*B(J,4,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_UPDATE_WL_4_1(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*32 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,125
            OUT(J,K,I)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,1,I)*B(J,1,K)+A(K,5,I)*B(J,0,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,1,I)*B(J,2,K)+A(K,2,I)*B(J,1,K)
     $       +A(K,6,I)*B(J,0,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,2,I)*B(J,2,K)+A(K,7,I)*B(J,0,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,1,I)*B(J,3,K)+A(K,3,I)*B(J,1,K)
     $       +A(K,8,I)*B(J,0,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,2,I)*B(J,3,K)+A(K,3,I)*B(J,2,K)
     $       +A(K,9,I)*B(J,0,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,3,I)*B(J,3,K)+A(K,10,I)
     $       *B(J,0,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,1,I)*B(J,4,K)+A(K,4,I)
     $       *B(J,1,K)+A(K,11,I)*B(J,0,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,2,I)*B(J,4,K)+A(K,4,I)
     $       *B(J,2,K)+A(K,12,I)*B(J,0,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,3,I)*B(J,4,K)+A(K,4,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,0,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,4,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,0,K)
            OUT(J,15,I)=OUT(J,15,I)+A(K,5,I)*B(J,1,K)+A(K,15,I)
     $       *B(J,0,K)
            OUT(J,16,I)=OUT(J,16,I)+A(K,5,I)*B(J,2,K)+A(K,6,I)
     $       *B(J,1,K)+A(K,16,I)*B(J,0,K)
            OUT(J,17,I)=OUT(J,17,I)+A(K,6,I)*B(J,2,K)+A(K,7,I)
     $       *B(J,1,K)+A(K,17,I)*B(J,0,K)
            OUT(J,18,I)=OUT(J,18,I)+A(K,7,I)*B(J,2,K)+A(K,18,I)
     $       *B(J,0,K)
            OUT(J,19,I)=OUT(J,19,I)+A(K,5,I)*B(J,3,K)+A(K,8,I)
     $       *B(J,1,K)+A(K,19,I)*B(J,0,K)
            OUT(J,20,I)=OUT(J,20,I)+A(K,6,I)*B(J,3,K)+A(K,8,I)
     $       *B(J,2,K)+A(K,9,I)*B(J,1,K)+A(K,20,I)*B(J,0,K)
            OUT(J,21,I)=OUT(J,21,I)+A(K,7,I)*B(J,3,K)+A(K,9,I)
     $       *B(J,2,K)+A(K,21,I)*B(J,0,K)
            OUT(J,22,I)=OUT(J,22,I)+A(K,8,I)*B(J,3,K)+A(K,10,I)
     $       *B(J,1,K)+A(K,22,I)*B(J,0,K)
            OUT(J,23,I)=OUT(J,23,I)+A(K,9,I)*B(J,3,K)+A(K,10,I)
     $       *B(J,2,K)+A(K,23,I)*B(J,0,K)
            OUT(J,24,I)=OUT(J,24,I)+A(K,10,I)*B(J,3,K)+A(K,24,I)
     $       *B(J,0,K)
            OUT(J,25,I)=OUT(J,25,I)+A(K,5,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,1,K)+A(K,25,I)*B(J,0,K)
            OUT(J,26,I)=OUT(J,26,I)+A(K,6,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,2,K)+A(K,12,I)*B(J,1,K)+A(K,26,I)*B(J,0,K)
            OUT(J,27,I)=OUT(J,27,I)+A(K,7,I)*B(J,4,K)+A(K,12,I)
     $       *B(J,2,K)+A(K,27,I)*B(J,0,K)
            OUT(J,28,I)=OUT(J,28,I)+A(K,8,I)*B(J,4,K)+A(K,11,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,1,K)+A(K,28,I)*B(J,0,K)
            OUT(J,29,I)=OUT(J,29,I)+A(K,9,I)*B(J,4,K)+A(K,12,I)
     $       *B(J,3,K)+A(K,13,I)*B(J,2,K)+A(K,29,I)*B(J,0,K)
            OUT(J,30,I)=OUT(J,30,I)+A(K,10,I)*B(J,4,K)+A(K,13,I)
     $       *B(J,3,K)+A(K,30,I)*B(J,0,K)
            OUT(J,31,I)=OUT(J,31,I)+A(K,11,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,1,K)+A(K,31,I)*B(J,0,K)
            OUT(J,32,I)=OUT(J,32,I)+A(K,12,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,2,K)+A(K,32,I)*B(J,0,K)
            OUT(J,33,I)=OUT(J,33,I)+A(K,13,I)*B(J,4,K)+A(K,14,I)
     $       *B(J,3,K)+A(K,33,I)*B(J,0,K)
            OUT(J,34,I)=OUT(J,34,I)+A(K,14,I)*B(J,4,K)+A(K,34,I)
     $       *B(J,0,K)
            OUT(J,35,I)=OUT(J,35,I)+A(K,15,I)*B(J,1,K)+A(K,35,I)
     $       *B(J,0,K)
            OUT(J,36,I)=OUT(J,36,I)+A(K,15,I)*B(J,2,K)+A(K,16,I)
     $       *B(J,1,K)+A(K,36,I)*B(J,0,K)
            OUT(J,37,I)=OUT(J,37,I)+A(K,16,I)*B(J,2,K)+A(K,17,I)
     $       *B(J,1,K)+A(K,37,I)*B(J,0,K)
            OUT(J,38,I)=OUT(J,38,I)+A(K,17,I)*B(J,2,K)+A(K,18,I)
     $       *B(J,1,K)+A(K,38,I)*B(J,0,K)
            OUT(J,39,I)=OUT(J,39,I)+A(K,18,I)*B(J,2,K)+A(K,39,I)
     $       *B(J,0,K)
            OUT(J,40,I)=OUT(J,40,I)+A(K,15,I)*B(J,3,K)+A(K,19,I)
     $       *B(J,1,K)+A(K,40,I)*B(J,0,K)
            OUT(J,41,I)=OUT(J,41,I)+A(K,16,I)*B(J,3,K)+A(K,19,I)
     $       *B(J,2,K)+A(K,20,I)*B(J,1,K)+A(K,41,I)*B(J,0,K)
            OUT(J,42,I)=OUT(J,42,I)+A(K,17,I)*B(J,3,K)+A(K,20,I)
     $       *B(J,2,K)+A(K,21,I)*B(J,1,K)+A(K,42,I)*B(J,0,K)
            OUT(J,43,I)=OUT(J,43,I)+A(K,18,I)*B(J,3,K)+A(K,21,I)
     $       *B(J,2,K)+A(K,43,I)*B(J,0,K)
            OUT(J,44,I)=OUT(J,44,I)+A(K,19,I)*B(J,3,K)+A(K,22,I)
     $       *B(J,1,K)+A(K,44,I)*B(J,0,K)
            OUT(J,45,I)=OUT(J,45,I)+A(K,20,I)*B(J,3,K)+A(K,22,I)
     $       *B(J,2,K)+A(K,23,I)*B(J,1,K)+A(K,45,I)*B(J,0,K)
            OUT(J,46,I)=OUT(J,46,I)+A(K,21,I)*B(J,3,K)+A(K,23,I)
     $       *B(J,2,K)+A(K,46,I)*B(J,0,K)
            OUT(J,47,I)=OUT(J,47,I)+A(K,22,I)*B(J,3,K)+A(K,24,I)
     $       *B(J,1,K)+A(K,47,I)*B(J,0,K)
            OUT(J,48,I)=OUT(J,48,I)+A(K,23,I)*B(J,3,K)+A(K,24,I)
     $       *B(J,2,K)+A(K,48,I)*B(J,0,K)
            OUT(J,49,I)=OUT(J,49,I)+A(K,24,I)*B(J,3,K)+A(K,49,I)
     $       *B(J,0,K)
            OUT(J,50,I)=OUT(J,50,I)+A(K,15,I)*B(J,4,K)+A(K,25,I)
     $       *B(J,1,K)+A(K,50,I)*B(J,0,K)
            OUT(J,51,I)=OUT(J,51,I)+A(K,16,I)*B(J,4,K)+A(K,25,I)
     $       *B(J,2,K)+A(K,26,I)*B(J,1,K)+A(K,51,I)*B(J,0,K)
            OUT(J,52,I)=OUT(J,52,I)+A(K,17,I)*B(J,4,K)+A(K,26,I)
     $       *B(J,2,K)+A(K,27,I)*B(J,1,K)+A(K,52,I)*B(J,0,K)
            OUT(J,53,I)=OUT(J,53,I)+A(K,18,I)*B(J,4,K)+A(K,27,I)
     $       *B(J,2,K)+A(K,53,I)*B(J,0,K)
            OUT(J,54,I)=OUT(J,54,I)+A(K,19,I)*B(J,4,K)+A(K,25,I)
     $       *B(J,3,K)+A(K,28,I)*B(J,1,K)+A(K,54,I)*B(J,0,K)
            OUT(J,55,I)=OUT(J,55,I)+A(K,20,I)*B(J,4,K)+A(K,26,I)
     $       *B(J,3,K)+A(K,28,I)*B(J,2,K)+A(K,29,I)*B(J,1,K)+A(K,55,I)
     $       *B(J,0,K)
            OUT(J,56,I)=OUT(J,56,I)+A(K,21,I)*B(J,4,K)+A(K,27,I)
     $       *B(J,3,K)+A(K,29,I)*B(J,2,K)+A(K,56,I)*B(J,0,K)
            OUT(J,57,I)=OUT(J,57,I)+A(K,22,I)*B(J,4,K)+A(K,28,I)
     $       *B(J,3,K)+A(K,30,I)*B(J,1,K)+A(K,57,I)*B(J,0,K)
            OUT(J,58,I)=OUT(J,58,I)+A(K,23,I)*B(J,4,K)+A(K,29,I)
     $       *B(J,3,K)+A(K,30,I)*B(J,2,K)+A(K,58,I)*B(J,0,K)
            OUT(J,59,I)=OUT(J,59,I)+A(K,24,I)*B(J,4,K)+A(K,30,I)
     $       *B(J,3,K)+A(K,59,I)*B(J,0,K)
            OUT(J,60,I)=OUT(J,60,I)+A(K,25,I)*B(J,4,K)+A(K,31,I)
     $       *B(J,1,K)+A(K,60,I)*B(J,0,K)
            OUT(J,61,I)=OUT(J,61,I)+A(K,26,I)*B(J,4,K)+A(K,31,I)
     $       *B(J,2,K)+A(K,32,I)*B(J,1,K)+A(K,61,I)*B(J,0,K)
            OUT(J,62,I)=OUT(J,62,I)+A(K,27,I)*B(J,4,K)+A(K,32,I)
     $       *B(J,2,K)+A(K,62,I)*B(J,0,K)
            OUT(J,63,I)=OUT(J,63,I)+A(K,28,I)*B(J,4,K)+A(K,31,I)
     $       *B(J,3,K)+A(K,33,I)*B(J,1,K)+A(K,63,I)*B(J,0,K)
            OUT(J,64,I)=OUT(J,64,I)+A(K,29,I)*B(J,4,K)+A(K,32,I)
     $       *B(J,3,K)+A(K,33,I)*B(J,2,K)+A(K,64,I)*B(J,0,K)
            OUT(J,65,I)=OUT(J,65,I)+A(K,30,I)*B(J,4,K)+A(K,33,I)
     $       *B(J,3,K)+A(K,65,I)*B(J,0,K)
            OUT(J,66,I)=OUT(J,66,I)+A(K,31,I)*B(J,4,K)+A(K,34,I)
     $       *B(J,1,K)+A(K,66,I)*B(J,0,K)
            OUT(J,67,I)=OUT(J,67,I)+A(K,32,I)*B(J,4,K)+A(K,34,I)
     $       *B(J,2,K)+A(K,67,I)*B(J,0,K)
            OUT(J,68,I)=OUT(J,68,I)+A(K,33,I)*B(J,4,K)+A(K,34,I)
     $       *B(J,3,K)+A(K,68,I)*B(J,0,K)
            OUT(J,69,I)=OUT(J,69,I)+A(K,34,I)*B(J,4,K)+A(K,69,I)
     $       *B(J,0,K)
            OUT(J,70,I)=OUT(J,70,I)+A(K,35,I)*B(J,1,K)
            OUT(J,71,I)=OUT(J,71,I)+A(K,35,I)*B(J,2,K)+A(K,36,I)
     $       *B(J,1,K)
            OUT(J,72,I)=OUT(J,72,I)+A(K,36,I)*B(J,2,K)+A(K,37,I)
     $       *B(J,1,K)
            OUT(J,73,I)=OUT(J,73,I)+A(K,37,I)*B(J,2,K)+A(K,38,I)
     $       *B(J,1,K)
            OUT(J,74,I)=OUT(J,74,I)+A(K,38,I)*B(J,2,K)+A(K,39,I)
     $       *B(J,1,K)
            OUT(J,75,I)=OUT(J,75,I)+A(K,39,I)*B(J,2,K)
            OUT(J,76,I)=OUT(J,76,I)+A(K,35,I)*B(J,3,K)+A(K,40,I)
     $       *B(J,1,K)
            OUT(J,77,I)=OUT(J,77,I)+A(K,36,I)*B(J,3,K)+A(K,40,I)
     $       *B(J,2,K)+A(K,41,I)*B(J,1,K)
            OUT(J,78,I)=OUT(J,78,I)+A(K,37,I)*B(J,3,K)+A(K,41,I)
     $       *B(J,2,K)+A(K,42,I)*B(J,1,K)
            OUT(J,79,I)=OUT(J,79,I)+A(K,38,I)*B(J,3,K)+A(K,42,I)
     $       *B(J,2,K)+A(K,43,I)*B(J,1,K)
            OUT(J,80,I)=OUT(J,80,I)+A(K,39,I)*B(J,3,K)+A(K,43,I)
     $       *B(J,2,K)
            OUT(J,81,I)=OUT(J,81,I)+A(K,40,I)*B(J,3,K)+A(K,44,I)
     $       *B(J,1,K)
            OUT(J,82,I)=OUT(J,82,I)+A(K,41,I)*B(J,3,K)+A(K,44,I)
     $       *B(J,2,K)+A(K,45,I)*B(J,1,K)
            OUT(J,83,I)=OUT(J,83,I)+A(K,42,I)*B(J,3,K)+A(K,45,I)
     $       *B(J,2,K)+A(K,46,I)*B(J,1,K)
            OUT(J,84,I)=OUT(J,84,I)+A(K,43,I)*B(J,3,K)+A(K,46,I)
     $       *B(J,2,K)
            OUT(J,85,I)=OUT(J,85,I)+A(K,44,I)*B(J,3,K)+A(K,47,I)
     $       *B(J,1,K)
            OUT(J,86,I)=OUT(J,86,I)+A(K,45,I)*B(J,3,K)+A(K,47,I)
     $       *B(J,2,K)+A(K,48,I)*B(J,1,K)
            OUT(J,87,I)=OUT(J,87,I)+A(K,46,I)*B(J,3,K)+A(K,48,I)
     $       *B(J,2,K)
            OUT(J,88,I)=OUT(J,88,I)+A(K,47,I)*B(J,3,K)+A(K,49,I)
     $       *B(J,1,K)
            OUT(J,89,I)=OUT(J,89,I)+A(K,48,I)*B(J,3,K)+A(K,49,I)
     $       *B(J,2,K)
            OUT(J,90,I)=OUT(J,90,I)+A(K,49,I)*B(J,3,K)
            OUT(J,91,I)=OUT(J,91,I)+A(K,35,I)*B(J,4,K)+A(K,50,I)
     $       *B(J,1,K)
            OUT(J,92,I)=OUT(J,92,I)+A(K,36,I)*B(J,4,K)+A(K,50,I)
     $       *B(J,2,K)+A(K,51,I)*B(J,1,K)
            OUT(J,93,I)=OUT(J,93,I)+A(K,37,I)*B(J,4,K)+A(K,51,I)
     $       *B(J,2,K)+A(K,52,I)*B(J,1,K)
            OUT(J,94,I)=OUT(J,94,I)+A(K,38,I)*B(J,4,K)+A(K,52,I)
     $       *B(J,2,K)+A(K,53,I)*B(J,1,K)
            OUT(J,95,I)=OUT(J,95,I)+A(K,39,I)*B(J,4,K)+A(K,53,I)
     $       *B(J,2,K)
            OUT(J,96,I)=OUT(J,96,I)+A(K,40,I)*B(J,4,K)+A(K,50,I)
     $       *B(J,3,K)+A(K,54,I)*B(J,1,K)
            OUT(J,97,I)=OUT(J,97,I)+A(K,41,I)*B(J,4,K)+A(K,51,I)
     $       *B(J,3,K)+A(K,54,I)*B(J,2,K)+A(K,55,I)*B(J,1,K)
            OUT(J,98,I)=OUT(J,98,I)+A(K,42,I)*B(J,4,K)+A(K,52,I)
     $       *B(J,3,K)+A(K,55,I)*B(J,2,K)+A(K,56,I)*B(J,1,K)
            OUT(J,99,I)=OUT(J,99,I)+A(K,43,I)*B(J,4,K)+A(K,53,I)
     $       *B(J,3,K)+A(K,56,I)*B(J,2,K)
            OUT(J,100,I)=OUT(J,100,I)+A(K,44,I)*B(J,4,K)+A(K,54,I)
     $       *B(J,3,K)+A(K,57,I)*B(J,1,K)
            OUT(J,101,I)=OUT(J,101,I)+A(K,45,I)*B(J,4,K)+A(K,55,I)
     $       *B(J,3,K)+A(K,57,I)*B(J,2,K)+A(K,58,I)*B(J,1,K)
            OUT(J,102,I)=OUT(J,102,I)+A(K,46,I)*B(J,4,K)+A(K,56,I)
     $       *B(J,3,K)+A(K,58,I)*B(J,2,K)
            OUT(J,103,I)=OUT(J,103,I)+A(K,47,I)*B(J,4,K)+A(K,57,I)
     $       *B(J,3,K)+A(K,59,I)*B(J,1,K)
            OUT(J,104,I)=OUT(J,104,I)+A(K,48,I)*B(J,4,K)+A(K,58,I)
     $       *B(J,3,K)+A(K,59,I)*B(J,2,K)
            OUT(J,105,I)=OUT(J,105,I)+A(K,49,I)*B(J,4,K)+A(K,59,I)
     $       *B(J,3,K)
            OUT(J,106,I)=OUT(J,106,I)+A(K,50,I)*B(J,4,K)+A(K,60,I)
     $       *B(J,1,K)
            OUT(J,107,I)=OUT(J,107,I)+A(K,51,I)*B(J,4,K)+A(K,60,I)
     $       *B(J,2,K)+A(K,61,I)*B(J,1,K)
            OUT(J,108,I)=OUT(J,108,I)+A(K,52,I)*B(J,4,K)+A(K,61,I)
     $       *B(J,2,K)+A(K,62,I)*B(J,1,K)
            OUT(J,109,I)=OUT(J,109,I)+A(K,53,I)*B(J,4,K)+A(K,62,I)
     $       *B(J,2,K)
            OUT(J,110,I)=OUT(J,110,I)+A(K,54,I)*B(J,4,K)+A(K,60,I)
     $       *B(J,3,K)+A(K,63,I)*B(J,1,K)
            OUT(J,111,I)=OUT(J,111,I)+A(K,55,I)*B(J,4,K)+A(K,61,I)
     $       *B(J,3,K)+A(K,63,I)*B(J,2,K)+A(K,64,I)*B(J,1,K)
            OUT(J,112,I)=OUT(J,112,I)+A(K,56,I)*B(J,4,K)+A(K,62,I)
     $       *B(J,3,K)+A(K,64,I)*B(J,2,K)
            OUT(J,113,I)=OUT(J,113,I)+A(K,57,I)*B(J,4,K)+A(K,63,I)
     $       *B(J,3,K)+A(K,65,I)*B(J,1,K)
            OUT(J,114,I)=OUT(J,114,I)+A(K,58,I)*B(J,4,K)+A(K,64,I)
     $       *B(J,3,K)+A(K,65,I)*B(J,2,K)
            OUT(J,115,I)=OUT(J,115,I)+A(K,59,I)*B(J,4,K)+A(K,65,I)
     $       *B(J,3,K)
            OUT(J,116,I)=OUT(J,116,I)+A(K,60,I)*B(J,4,K)+A(K,66,I)
     $       *B(J,1,K)
            OUT(J,117,I)=OUT(J,117,I)+A(K,61,I)*B(J,4,K)+A(K,66,I)
     $       *B(J,2,K)+A(K,67,I)*B(J,1,K)
            OUT(J,118,I)=OUT(J,118,I)+A(K,62,I)*B(J,4,K)+A(K,67,I)
     $       *B(J,2,K)
            OUT(J,119,I)=OUT(J,119,I)+A(K,63,I)*B(J,4,K)+A(K,66,I)
     $       *B(J,3,K)+A(K,68,I)*B(J,1,K)
            OUT(J,120,I)=OUT(J,120,I)+A(K,64,I)*B(J,4,K)+A(K,67,I)
     $       *B(J,3,K)+A(K,68,I)*B(J,2,K)
            OUT(J,121,I)=OUT(J,121,I)+A(K,65,I)*B(J,4,K)+A(K,68,I)
     $       *B(J,3,K)
            OUT(J,122,I)=OUT(J,122,I)+A(K,66,I)*B(J,4,K)+A(K,69,I)
     $       *B(J,1,K)
            OUT(J,123,I)=OUT(J,123,I)+A(K,67,I)*B(J,4,K)+A(K,69,I)
     $       *B(J,2,K)
            OUT(J,124,I)=OUT(J,124,I)+A(K,68,I)*B(J,4,K)+A(K,69,I)
     $       *B(J,3,K)
            OUT(J,125,I)=OUT(J,125,I)+A(K,69,I)*B(J,4,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE UPDATE_WL_1_1(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*16 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,14
            OUT(J,K,I)=(0.0D0,0.0D0)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,1,I)*B(J,1,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,1,I)*B(J,2,K)+A(K,2,I)*B(J,1,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,2,I)*B(J,2,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,1,I)*B(J,3,K)+A(K,3,I)*B(J,1,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,2,I)*B(J,3,K)+A(K,3,I)*B(J,2,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,3,I)*B(J,3,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,1,I)*B(J,4,K)+A(K,4,I)*B(J,1,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,2,I)*B(J,4,K)+A(K,4,I)*B(J,2,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,3,I)*B(J,4,K)+A(K,4,I)*B(J,3,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,4,I)*B(J,4,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_UPDATE_WL_1_1(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE,OUT)
      INCLUDE 'coef_specs.inc'
      INTEGER I,J,K
      COMPLEX*32 A(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOP_MAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,14
            OUT(J,K,I)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,1,I)*B(J,1,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,1,I)*B(J,2,K)+A(K,2,I)*B(J,1,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,2,I)*B(J,2,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,1,I)*B(J,3,K)+A(K,3,I)*B(J,1,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,2,I)*B(J,3,K)+A(K,3,I)*B(J,2,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,3,I)*B(J,3,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,1,I)*B(J,4,K)+A(K,4,I)*B(J,1,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,2,I)*B(J,4,K)+A(K,4,I)*B(J,2,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,3,I)*B(J,4,K)+A(K,4,I)*B(J,3,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,4,I)*B(J,4,K)
          ENDDO
        ENDDO
      ENDDO
      END
