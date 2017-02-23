      SUBROUTINE SBORN_HEL(P,ANS_SUMMED)
C     
C     Return the sum of the split orders which are required in
C      orders.inc (BORN_ORDERS)
C     Also the values needed for the counterterms are stored in the
C      C_BORN_CNT common block
C     
C     
C     CONSTANTS
C     
      IMPLICIT NONE
      INCLUDE 'nexternal.inc'
      INTEGER NSQAMPSO
      PARAMETER (NSQAMPSO=1)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL), ANS_SUMMED
C     
C     VARIABLES
C     
      INTEGER I, J
      INCLUDE 'orders.inc'
      REAL*8 ANS(0:NSQAMPSO)
      LOGICAL KEEP_ORDER(NSQAMPSO), FIRSTTIME
      DATA KEEP_ORDER / NSQAMPSO * .TRUE. /
      DATA FIRSTTIME / .TRUE. /
      INCLUDE 'born_nhel.inc'
      DOUBLE PRECISION WGT_HEL(NSQAMPSO, MAX_BHEL)
      COMMON/C_BORN_HEL_SPLIT/WGT_HEL
      DOUBLE PRECISION WGT_HEL_SUMMED(MAX_BHEL)
      COMMON/C_BORN_HEL/WGT_HEL_SUMMED
C     
C     FUNCTIONS
C     
      INTEGER GETORDPOWFROMINDEX_B
C     
C     BEGIN CODE
C     
C     look for orders which match the born order constraint 

      IF (FIRSTTIME) THEN
        DO I = 1, NSQAMPSO
C         this is for the orders of the born to integrate
          DO J = 1, NSPLITORDERS
            IF(GETORDPOWFROMINDEX_B(J, I) .GT. BORN_ORDERS(J)) THEN
              KEEP_ORDER(I) = .FALSE.
              EXIT
            ENDIF
          ENDDO
        ENDDO
        FIRSTTIME = .FALSE.
      ENDIF

C     look for orders which match the born order constraint 
      CALL SBORN_HEL_SPLITORDERS(P,ANS)
      ANS_SUMMED = 0D0
      DO J = 1, MAX_BHEL
        WGT_HEL_SUMMED(J) = 0D0
      ENDDO
      DO I = 1, NSQAMPSO
        IF (KEEP_ORDER(I)) THEN
          ANS_SUMMED = ANS_SUMMED + ANS(I)
          DO J = 1, MAX_BHEL
            WGT_HEL_SUMMED(J) = WGT_HEL_SUMMED(J) + WGT_HEL(I,J)
          ENDDO
        ENDIF
      ENDDO

      RETURN
      END


      SUBROUTINE SBORN_HEL_SPLITORDERS(P1,ANS)
C     
C     Generated by MadGraph5_aMC@NLO v. %(version)s, %(date)s
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C     AND HELICITIES
C     FOR THE POINT IN PHASE SPACE P1(0:3,NEXTERNAL-1)
C     
C     Process: g g > t t~ [ real = QED QCD ] QCD^2=4 QED^2=2
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INCLUDE 'nexternal.inc'
      INCLUDE 'born_nhel.inc'
      INTEGER     NCOMB
      PARAMETER ( NCOMB=  16 )
      INTEGER NSQAMPSO
      PARAMETER (NSQAMPSO=1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*10)
      INTEGER NGRAPHS
      PARAMETER (NGRAPHS=   3)
C     
C     ARGUMENTS 
C     
      REAL*8 P1(0:3,NEXTERNAL-1),ANS(0:NSQAMPSO)
C     
C     LOCAL VARIABLES 
C     
      INTEGER IHEL,IDEN,I,J
      DOUBLE PRECISION T(NSQAMPSO)
      INTEGER IDEN_VALUES(10)
      DATA IDEN_VALUES /256, 256, 256, 256, 256, 256, 256, 256, 256,
     $  256/
C     
C     GLOBAL VARIABLES
C     
      LOGICAL GOODHEL(NCOMB,10)
      COMMON /C_GOODHEL/ GOODHEL
      DOUBLE PRECISION SAVEMOM(NEXTERNAL-1,2)
      COMMON/TO_SAVEMOM/SAVEMOM
      LOGICAL CALCULATEDBORN
      COMMON/CCALCULATEDBORN/CALCULATEDBORN
      INTEGER NFKSPROCESS
      COMMON/C_NFKSPROCESS/NFKSPROCESS
      DOUBLE PRECISION WGT_HEL(NSQAMPSO, MAX_BHEL)
      COMMON/C_BORN_HEL_SPLIT/WGT_HEL
C     ----------
C     BEGIN CODE
C     ----------
      IDEN=IDEN_VALUES(NFKSPROCESS)
      IF (CALCULATEDBORN) THEN
        DO J=1,NEXTERNAL-1
          IF (SAVEMOM(J,1).NE.P1(0,J) .OR. SAVEMOM(J,2).NE.P1(3,J))
     $      THEN
            CALCULATEDBORN=.FALSE.
            WRITE(*,*) 'Error in sborn_hel_splitorders: momenta not'
     $       //' the same in the born'
            STOP
          ENDIF
        ENDDO
      ELSE
        WRITE(*,*) 'Error in sborn_hel_splitorders: this should be'
     $   //' called only with calculatedborn = true'
        STOP
      ENDIF
      DO I=0,NSQAMPSO
        ANS(I) = 0D0
      ENDDO
      DO IHEL=1,NCOMB
        IF (GOODHEL(IHEL,NFKSPROCESS)) THEN
          CALL BORN_HEL_SPLITORDERS(P1,IHEL,T)
          DO I=1,NSQAMPSO
            WGT_HEL(I, IHEL) = T(I) / DBLE(IDEN)
            ANS(I)=ANS(I)+T(I)
          ENDDO
        ENDIF
      ENDDO
      DO I=1,NSQAMPSO
        ANS(I)=ANS(I)/DBLE(IDEN)
        ANS(0)=ANS(0)+ANS(I)
      ENDDO
      END


      SUBROUTINE BORN_HEL_SPLITORDERS(P,HELL,ANS)
C     
C     Generated by MadGraph5_aMC@NLO v. %(version)s, %(date)s
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C     FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL-1)

C     Process: g g > t t~ [ real = QED QCD ] QCD^2=4 QED^2=2
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER NAMPSO, NSQAMPSO
      PARAMETER (NAMPSO=1, NSQAMPSO=1)
      INTEGER     NGRAPHS
      PARAMETER ( NGRAPHS = 3 )
      INTEGER NCOLOR
      PARAMETER (NCOLOR=2)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1 = (0D0,1D0))
      INCLUDE 'nexternal.inc'
      INCLUDE 'born_nhel.inc'
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL-1)
      INTEGER HELL
      REAL*8 ANS(NSQAMPSO)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,M,N
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 ZTEMP, AMP(NGRAPHS), JAMP(NCOLOR,NAMPSO)
C     
C     GLOBAL VARIABLES
C     
      DOUBLE COMPLEX SAVEAMP(NGRAPHS,MAX_BHEL)
      COMMON/TO_SAVEAMP/SAVEAMP
      LOGICAL CALCULATEDBORN
      COMMON/CCALCULATEDBORN/CALCULATEDBORN
C     
C     FUNCTION
C     
      INTEGER SQSOINDEXB
C     
C     COLOR DATA
C     
      DATA DENOM(1)/3/
      DATA (CF(I,  1),I=  1,  2) /   16,   -2/
C     1 T(1,2,3,4)
      DATA DENOM(2)/3/
      DATA (CF(I,  2),I=  1,  2) /   -2,   16/
C     1 T(2,1,3,4)
C     ----------
C     BEGIN CODE
C     ----------
      IF (.NOT. CALCULATEDBORN) THEN
        WRITE(*,*) 'Error in b_sf: color_linked borns should be called'
     $   //' only with calculatedborn = true'
        STOP
      ELSEIF (CALCULATEDBORN) THEN
        DO I=1,NGRAPHS
          AMP(I)=SAVEAMP(I,HELL)
        ENDDO
      ENDIF
C     JAMPs contributing to orders QCD=2 QED=0
      JAMP(1,1)=+IMAG1*AMP(1)-AMP(2)
      JAMP(2,1)=-IMAG1*AMP(1)-AMP(3)
      DO I = 1, NSQAMPSO
        ANS(I) = 0D0
      ENDDO
      DO M = 1, NAMPSO
        DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
            ZTEMP = ZTEMP + CF(J,I)*JAMP(J,M)
          ENDDO
          DO N = 1, NAMPSO
            ANS(SQSOINDEXB(M,N))=ANS(SQSOINDEXB(M,N))+ZTEMP
     $       *DCONJG(JAMP(I,N))/DENOM(I)
          ENDDO
        ENDDO
      ENDDO
      END



