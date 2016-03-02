      SUBROUTINE ML5_0_TIRLOOP(I_SQSO,I_LOOPGROUP,I_LIB,NLOOPLINE,PL
     $ ,M2L,RANK,RES,STABLE)
C     
C     Generated by MadGraph5_aMC@NLO v. %(version)s, %(date)s
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     Interface between MG5 and TIR.
C     
C     Process: d u~ > m- vm~ g QED<=2 QCD<=1 [ virt = QCD ]
C     
C     
C     CONSTANTS 
C     
      INTEGER NLOOPGROUPS
      PARAMETER (NLOOPGROUPS=9)
C     These are constants related to the split orders
      INTEGER NSQUAREDSO
      PARAMETER (NSQUAREDSO=0)
      INCLUDE 'loop_max_coefs.inc'
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=5)
      LOGICAL CHECKPCONSERVATION
      PARAMETER (CHECKPCONSERVATION=.TRUE.)
      REAL*8 NORMALIZATION
      PARAMETER (NORMALIZATION = 1.D0/(16.D0*3.14159265358979323846D0*
     $ *2))
C     
C     ARGUMENTS 
C     
      INTEGER I_SQSO,I_LOOPGROUP,I_LIB
      INTEGER NLOOPLINE, RANK
      REAL*8 PL(0:3,NLOOPLINE)
      REAL*8 PCT(0:3,0:NLOOPLINE-1)
      REAL*8 PDEN(0:3,NLOOPLINE-1)
      COMPLEX*16 M2L(NLOOPLINE)
      COMPLEX*16 M2LCT(0:NLOOPLINE-1)
      COMPLEX*16 RES(3)
      LOGICAL STABLE
C     
C     LOCAL VARIABLES 
C     
      INTEGER I, J, K
      INTEGER NLOOPCOEFS
      LOGICAL CTINIT, TIRINIT, GOLEMINIT, SAMURAIINIT, NINJAINIT
      COMMON/REDUCTIONCODEINIT/CTINIT,TIRINIT,GOLEMINIT,SAMURAIINIT
     $ ,NINJAINIT

C     This variable will be used to detect changes in the TIR library
C      used so as to force the reset of the TIR filter.
      INTEGER LAST_LIB_USED
      DATA LAST_LIB_USED/-1/

      COMPLEX*16 TIRCOEFS(0:LOOPMAXCOEFS-1,3)
      COMPLEX*16 PJCOEFS(0:LOOPMAXCOEFS-1,3)
C     
C     EXTERNAL FUNCTIONS
C     
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'MadLoopParams.inc'
      INCLUDE 'coupl.inc'
      INTEGER CTMODE
      REAL*8 LSCALE
      COMMON/ML5_0_CT/LSCALE,CTMODE

C     The variables below are just for monitoring purposes. 
      INTEGER ID,SQSOINDEX,R
      COMMON/ML5_0_LOOP/ID,SQSOINDEX,R

C     The argument ILIB is the TIR library to be used for that
C      specific library.
      INTEGER LIBINDEX
      COMMON/ML5_0_I_LIB/LIBINDEX

      COMPLEX*16 LOOPCOEFS(0:LOOPMAXCOEFS-1,NSQUAREDSO,NLOOPGROUPS)
      COMMON/ML5_0_LCOEFS/LOOPCOEFS
C     ----------
C     BEGIN CODE
C     ----------

C     Initialize for the very first time ML is called LAST_ILIB with
C      the first ILIB used.
      IF(LAST_LIB_USED.EQ.-1) THEN
        LAST_LIB_USED = MLREDUCTIONLIB(LIBINDEX)
      ELSE
C       We changed the TIR library so we must refresh the cache.
        IF(MLREDUCTIONLIB(LIBINDEX).NE.LAST_LIB_USED) THEN
          LAST_LIB_USED = MLREDUCTIONLIB(LIBINDEX)
          CALL ML5_0_CLEAR_TIR_CACHE()
        ENDIF
      ENDIF

      IF (MLREDUCTIONLIB(I_LIB).EQ.4) THEN
C       Golem95 not available
        WRITE(*,*) 'ERROR:: Golem95 is not interfaced.'
        STOP
      ENDIF

C     INITIALIZE TIR IF NEEDED
      IF (TIRINIT) THEN
        TIRINIT=.FALSE.
        CALL ML5_0_INITTIR()
      ENDIF

C     CONVERT THE MASSES TO BE COMPLEX
      DO I=1,NLOOPLINE
        M2LCT(I-1)=M2L(I)
      ENDDO

C     CONVERT THE MOMENTA FLOWING IN THE LOOP LINES TO CT CONVENTIONS
      DO I=0,3
        DO J=0,(NLOOPLINE-1)
          PCT(I,J)=0.D0
        ENDDO
      ENDDO
      DO I=0,3
        DO J=1,NLOOPLINE
          PCT(I,0)=PCT(I,0)+PL(I,J)
        ENDDO
      ENDDO
      IF (CHECKPCONSERVATION) THEN
        IF (PCT(0,0).GT.1.D-6) THEN
          WRITE(*,*) 'energy is not conserved ',PCT(0,0)
          STOP 'energy is not conserved'
        ELSEIF (PCT(1,0).GT.1.D-6) THEN
          WRITE(*,*) 'px is not conserved ',PCT(1,0)
          STOP 'px is not conserved'
        ELSEIF (PCT(2,0).GT.1.D-6) THEN
          WRITE(*,*) 'py is not conserved ',PCT(2,0)
          STOP 'py is not conserved'
        ELSEIF (PCT(3,0).GT.1.D-6) THEN
          WRITE(*,*) 'pz is not conserved ',PCT(3,0)
          STOP 'pz is not conserved'
        ENDIF
      ENDIF
      DO I=0,3
        DO J=1,(NLOOPLINE-1)
          DO K=1,J
            PCT(I,J)=PCT(I,J)+PL(I,K)
          ENDDO
        ENDDO
      ENDDO

      DO I=0,3
        DO J=1,(NLOOPLINE-1)
          PDEN(I,J)=PCT(I,J)
        ENDDO
      ENDDO
C     NUMBER OF INDEPEDENT LOOPCOEFS FOR RANK=RANK
      NLOOPCOEFS=0
      DO I=0,RANK
        NLOOPCOEFS=NLOOPCOEFS+(3+I)*(2+I)*(1+I)/6
      ENDDO
      SELECT CASE(MLREDUCTIONLIB(I_LIB))
      CASE(2)
C     PJFry++
      WRITE(*,*) 'ERROR:: PJFRY++ is not interfaced.'
      STOP
      CASE(3)
C     IREGI
      WRITE(*,*) 'ERROR:: IREGI is not interfaced.'
      STOP
      END SELECT
      DO I=1,3
        RES(I)=(0.0D0,0.0D0)
        DO J=0,NLOOPCOEFS-1
          RES(I)=RES(I)+LOOPCOEFS(J,I_SQSO,I_LOOPGROUP)*TIRCOEFS(J,I)
        ENDDO
      ENDDO
      RES(1)=NORMALIZATION*2.0D0*DBLE(RES(1))
      RES(2)=NORMALIZATION*2.0D0*DBLE(RES(2))
      RES(3)=NORMALIZATION*2.0D0*DBLE(RES(3))
C     IF(MLReductionLib(I_LIB).EQ.2) THEN
C     WRITE(*,*) 'PJFry: Loop ID',ID,' =',RES(1),RES(2),RES(3)
C     ELSEIF(MLReductionLib(I_LIB).EQ.3) THEN
C     WRITE(*,*) 'Iregi: Loop ID',ID,' =',RES(1),RES(2),RES(3)
C     ENDIF
      END

      SUBROUTINE ML5_0_SWITCH_ORDER(CTMODE,NLOOPLINE,PL,PDEN,M2L)
      IMPLICIT NONE

      INTEGER CTMODE,NLOOPLINE

      REAL*8 PL(0:3,NLOOPLINE)
      REAL*8 PDEN(0:3,NLOOPLINE-1)
      COMPLEX*16 M2L(NLOOPLINE)
      REAL*8 NEW_PL(0:3,NLOOPLINE)
      REAL*8 NEW_PDEN(0:3,NLOOPLINE-1)
      COMPLEX*16 NEW_M2L(NLOOPLINE)

      INTEGER I,J,K

      IF (CTMODE.NE.2.AND.CTMODE.NE.5) THEN
        RETURN
      ENDIF

      IF (NLOOPLINE.LE.2) THEN
        RETURN
      ENDIF

      DO I=1,NLOOPLINE-1
        DO J=0,3
          NEW_PDEN(J,NLOOPLINE-I) = PDEN(J,I)
        ENDDO
      ENDDO
      DO I=1,NLOOPLINE-1
        DO J=0,3
          PDEN(J,I) = NEW_PDEN(J,I)
        ENDDO
      ENDDO

      DO I=2,NLOOPLINE
        NEW_M2L(I)=M2L(NLOOPLINE-I+2)
      ENDDO
      DO I=2,NLOOPLINE
        M2L(I)=NEW_M2L(I)
      ENDDO


      DO I=1,NLOOPLINE
        DO J=0,3
          NEW_PL(J,I) = -PL(J,NLOOPLINE+1-I)
        ENDDO
      ENDDO
      DO I=1,NLOOPLINE
        DO J=0,3
          PL(J,I) = NEW_PL(J,I)
        ENDDO
      ENDDO

      END

      SUBROUTINE ML5_0_MP_SWITCH_ORDER(CTMODE,NLOOPLINE,PL,PDEN,M2L)
      IMPLICIT NONE

      INTEGER CTMODE,NLOOPLINE

      REAL*16 PL(0:3,NLOOPLINE)
      REAL*16 PDEN(0:3,NLOOPLINE-1)
      COMPLEX*32 M2L(NLOOPLINE)
      REAL*16 NEW_PL(0:3,NLOOPLINE)
      REAL*16 NEW_PDEN(0:3,NLOOPLINE-1)
      COMPLEX*32 NEW_M2L(NLOOPLINE)

      INTEGER I,J,K

      IF (CTMODE.NE.2.AND.CTMODE.NE.5) THEN
        RETURN
      ENDIF

      IF (NLOOPLINE.LE.2) THEN
        RETURN
      ENDIF

      DO I=1,NLOOPLINE-1
        DO J=0,3
          NEW_PDEN(J,NLOOPLINE-I) = PDEN(J,I)
        ENDDO
      ENDDO
      DO I=1,NLOOPLINE-1
        DO J=0,3
          PDEN(J,I) = NEW_PDEN(J,I)
        ENDDO
      ENDDO

      DO I=2,NLOOPLINE
        NEW_M2L(I)=M2L(NLOOPLINE-I+2)
      ENDDO
      DO I=2,NLOOPLINE
        M2L(I)=NEW_M2L(I)
      ENDDO


      DO I=1,NLOOPLINE
        DO J=0,3
          NEW_PL(J,I) = -PL(J,NLOOPLINE+1-I)
        ENDDO
      ENDDO
      DO I=1,NLOOPLINE
        DO J=0,3
          PL(J,I) = NEW_PL(J,I)
        ENDDO
      ENDDO

      END

      SUBROUTINE ML5_0_INITTIR()
C     
C     INITIALISATION OF TIR
C     
C     LOCAL VARIABLES 
C     
      REAL*8 THRS
      LOGICAL EXT_NUM_FOR_R1
C     
C     GLOBAL VARIABLES 
C     
      INCLUDE 'MadLoopParams.inc'
      LOGICAL CTINIT, TIRINIT, GOLEMINIT, SAMURAIINIT, NINJAINIT
      COMMON/REDUCTIONCODEINIT/CTINIT,TIRINIT,GOLEMINIT,SAMURAIINIT
     $ ,NINJAINIT

C     ----------
C     BEGIN CODE
C     ----------

C     DEFAULT PARAMETERS FOR TIR
C     -------------------------------  
C     THRS1 IS THE PRECISION LIMIT BELOW WHICH THE MP ROUTINES
C      ACTIVATES
C     USE THE SAME MADLOOP PARAMETER IN CUTTOOLS AND TIR
C     IT IS NECESSARY TO INITIALIZE CT BECAUSE IREGI USES THE VERSION
C      OF ONELOOP
C     FROM CUTTOOLS LIBRARY
      THRS=CTSTABTHRES
C     LOOPLIB SET WHAT LIBRARY CT USES
C     1 -> LOOPTOOLS
C     2 -> AVH
C     3 -> QCDLOOP
      LOOPLIB=CTLOOPLIBRARY
C     The initialization below is for CT v1.9.2+
      IF (CTINIT) THEN
        CTINIT=.FALSE.
        CALL ML5_0_INITCT()
      ENDIF
      END

      SUBROUTINE ML5_0_CHOOSE_LOOPLIB(LIBINDEX,NLOOPLINE,RANK
     $ ,COMPLEX_MASS,LOOP_ID,DOING_QP,I_LIB)
C     
C     CHOOSE THE CORRECT LOOP LIB
C     Example:
C     MLReductionLib=3|2|1 and LIBINDEX=3
C     IF THE LOOP IS BEYOND THE SCOPE OF LOOP LIB MLReductionLib(3)=1
C     USE LIBINDEX=1, and LIBINDEX=2 ...
C     IF IT IS STILL NOT GOOD,STOP
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER NLOOPLIB
      PARAMETER (NLOOPLIB=4)
      INTEGER QP_NLOOPLIB
      PARAMETER (QP_NLOOPLIB=1)
      INTEGER NLOOPGROUPS
      PARAMETER (NLOOPGROUPS=9)
C     
C     ARGUMENTS
C     
      INTEGER LIBINDEX,NLOOPLINE,RANK,I_LIB,LOOP_ID
      LOGICAL COMPLEX_MASS,DOING_QP
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J_LIB,LIBNUM,SELECT_LIBINDEX
      LOGICAL LPASS
C     This list specifies what loop involve an Higgs effective vertex
C      so that CutTools limitations can be correctly implemented
      LOGICAL HAS_AN_HEFT_VERTEX(NLOOPGROUPS)
      DATA (HAS_AN_HEFT_VERTEX(I),I=     1,     9) /.FALSE.,.FALSE.
     $ ,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE./
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'MadLoopParams.inc'
      INCLUDE 'process_info.inc'
C     Change the list 'LOOPLIBS_QPAVAILABLE' in loop_matrix_standalone.
C     inc to change the list of QPTools availables
      LOGICAL QP_TOOLS_AVAILABLE
      INTEGER INDEX_QP_TOOLS(QP_NLOOPLIB+1)
      COMMON/ML5_0_LOOP_TOOLS/QP_TOOLS_AVAILABLE,INDEX_QP_TOOLS

C     ----------
C     BEGIN CODE
C     ----------

      IF(DOING_QP)THEN
C       QP EVALUATION, ONLY CUTTOOLS
        IF(.NOT.QP_TOOLS_AVAILABLE)THEN
          STOP 'No qp tools available, please make sure MLReductionLi'
     $     //'b is correct'
        ENDIF
        J_LIB=0
        SELECT_LIBINDEX=LIBINDEX
        DO WHILE(J_LIB.EQ.0)
          DO I=1,QP_NLOOPLIB
            IF(INDEX_QP_TOOLS(I).EQ.SELECT_LIBINDEX)THEN
              J_LIB=I
              EXIT
            ENDIF
          ENDDO
          IF(J_LIB.EQ.0)THEN
            SELECT_LIBINDEX=SELECT_LIBINDEX+1
            IF(SELECT_LIBINDEX.GT.NLOOPLIB.OR.MLREDUCTIONLIB(SELECT_LIB
     $       INDEX).EQ.0)SELECT_LIBINDEX=1
          ENDIF
        ENDDO
        I=J_LIB
        I_LIB=SELECT_LIBINDEX
        LIBNUM=MLREDUCTIONLIB(I_LIB)
        DO
        CALL DETECT_LOOPLIB(LIBNUM,NLOOPLINE,RANK,COMPLEX_MASS
     $   ,HAS_AN_HEFT_VERTEX(LOOP_ID),MAX_SPIN_CONNECTED_TO_LOOP,LPASS)
        IF(LPASS)EXIT
        I=I+1
        IF(I.GT.QP_NLOOPLIB.AND.INDEX_QP_TOOLS(I).EQ.0)THEN
          I=1
        ENDIF
        IF(I.EQ.J_LIB)THEN
          STOP 'No qp loop library can deal with this integral'
        ENDIF
        I_LIB=INDEX_QP_TOOLS(I)
        LIBNUM=MLREDUCTIONLIB(I_LIB)
        ENDDO
      ELSE
C       DP EVALUATION
        I_LIB=LIBINDEX
        LIBNUM=MLREDUCTIONLIB(I_LIB)
        DO
        CALL DETECT_LOOPLIB(LIBNUM,NLOOPLINE,RANK,COMPLEX_MASS
     $   ,HAS_AN_HEFT_VERTEX(LOOP_ID),MAX_SPIN_CONNECTED_TO_LOOP,LPASS)
        IF(LPASS)EXIT
        I_LIB=I_LIB+1
        IF(I_LIB.GT.NLOOPLIB.OR.MLREDUCTIONLIB(I_LIB).EQ.0)THEN
          I_LIB=1
        ENDIF
        IF(I_LIB.EQ.LIBINDEX)THEN
          STOP 'No dp loop library can deal with this integral'
        ENDIF
        LIBNUM=MLREDUCTIONLIB(I_LIB)
        ENDDO
      ENDIF
      RETURN
      END

      SUBROUTINE ML5_0_CLEAR_TIR_CACHE()
C     No TIR caching implemented, this is dummy. (The subroutine is
C      kept as it might be called by the MC).
      CONTINUE
      END SUBROUTINE


