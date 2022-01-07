      SUBROUTINE LOOP_CT_CALLS_1(P,NHEL,H,IC)
C     
C     Modules
C     
      USE POLYNOMIAL_CONSTANTS
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=4)
      INTEGER    NCOMB
      PARAMETER (NCOMB=16)
      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=1)
      INTEGER    NLOOPS, NLOOPGROUPS, NCTAMPS
      PARAMETER (NLOOPS=35, NLOOPGROUPS=25, NCTAMPS=15)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=50)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=6,NLOOPWAVEFUNCS=73)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      REAL*16     MP__ZERO
      PARAMETER (MP__ZERO=0.0E0_16)
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=2, NSQUAREDSO=1, NAMPSO=2)
C     
C     ARGUMENTS
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      INTEGER H
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K
      COMPLEX*16 COEFS(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)

      LOGICAL DUMMYFALSE
      DATA DUMMYFALSE/.FALSE./
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
      INCLUDE 'mp_coupl.inc'

      INTEGER HELOFFSET
      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/FILTERS/GOODAMP,GOODHEL,HELOFFSET

      LOGICAL CHECKPHASE
      LOGICAL HELDOUBLECHECKED
      COMMON/INIT/CHECKPHASE, HELDOUBLECHECKED

      INTEGER SQSO_TARGET
      COMMON/SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

      INTEGER I_SO
      COMMON/I_SO/I_SO
      INTEGER I_LIB
      COMMON/I_LIB/I_LIB

      COMPLEX*16 AMP(NBORNAMPS)
      COMMON/AMPS/AMP
      COMPLEX*16 W(20,NWAVEFUNCS)
      COMMON/W/W

      COMPLEX*16 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE,
     $ -1:NLOOPWAVEFUNCS)
      COMPLEX*16 PL(0:3,-1:NLOOPWAVEFUNCS)
      COMMON/WL/WL,PL

      COMPLEX*16 AMPL(3,NCTAMPS)
      COMMON/AMPL/AMPL

C     
C     ----------
C     BEGIN CODE
C     ----------

C     The target squared split order contribution is already reached
C      if true.
      IF (FILTER_SO.AND.CTCALL_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

C     CutTools call for loop numbers 1,29,30,2
      CALL LOOP_2(5,6,DCMPLX(CMASS_MDL_MW),DCMPLX(ZERO),2,I_SO,1)
C     CutTools call for loop numbers 3
      CALL LOOP_3(3,4,5,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(CMASS_MDL_MW)
     $ ,2,I_SO,2)
C     CutTools call for loop numbers 4
      CALL LOOP_3(1,2,6,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO
     $ ,3)
C     CutTools call for loop numbers 5
      CALL LOOP_3(1,2,6,DCMPLX(ZERO),DCMPLX(CMASS_MDL_MW),DCMPLX(ZERO)
     $ ,2,I_SO,4)
C     CutTools call for loop numbers 6
      CALL LOOP_4(1,2,3,4,DCMPLX(ZERO),DCMPLX(CMASS_MDL_MW)
     $ ,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,5)
C     CutTools call for loop numbers 7
      CALL LOOP_3(1,2,6,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(CMASS_MDL_MW)
     $ ,2,I_SO,6)
C     CutTools call for loop numbers 8
      CALL LOOP_4(1,2,4,3,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(CMASS_MDL_MW),2,I_SO,7)
C     CutTools call for loop numbers 9,25
      CALL LOOP_1_2(2,5,6,DCMPLX(CMASS_MDL_MZ),0,I_SO,8)
C     CutTools call for loop numbers 10,26,31,32,11
      CALL LOOP_2(5,6,DCMPLX(CMASS_MDL_MW),DCMPLX(CMASS_MDL_MZ),2,I_SO
     $ ,9)
C     CutTools call for loop numbers 12
      CALL LOOP_3(3,4,5,DCMPLX(ZERO),DCMPLX(CMASS_MDL_MW)
     $ ,DCMPLX(CMASS_MDL_MZ),2,I_SO,10)
C     CutTools call for loop numbers 13
      CALL LOOP_3(3,4,5,DCMPLX(CMASS_MDL_MZ),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,2,I_SO,11)
C     CutTools call for loop numbers 14
      CALL LOOP_3(3,4,5,DCMPLX(ZERO),DCMPLX(CMASS_MDL_MZ)
     $ ,DCMPLX(CMASS_MDL_MW),2,I_SO,12)
C     CutTools call for loop numbers 15
      CALL LOOP_3(1,2,6,DCMPLX(CMASS_MDL_MZ),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,2,I_SO,13)
C     CutTools call for loop numbers 16
      CALL LOOP_3(1,2,6,DCMPLX(ZERO),DCMPLX(CMASS_MDL_MW)
     $ ,DCMPLX(CMASS_MDL_MZ),2,I_SO,14)
C     CutTools call for loop numbers 17
      CALL LOOP_4(1,2,4,3,DCMPLX(ZERO),DCMPLX(CMASS_MDL_MW)
     $ ,DCMPLX(ZERO),DCMPLX(CMASS_MDL_MZ),2,I_SO,15)
C     CutTools call for loop numbers 18
      CALL LOOP_4(1,2,3,4,DCMPLX(ZERO),DCMPLX(CMASS_MDL_MW)
     $ ,DCMPLX(ZERO),DCMPLX(CMASS_MDL_MZ),2,I_SO,16)
C     CutTools call for loop numbers 19
      CALL LOOP_3(1,2,6,DCMPLX(ZERO),DCMPLX(CMASS_MDL_MZ)
     $ ,DCMPLX(CMASS_MDL_MW),2,I_SO,17)
C     CutTools call for loop numbers 20
      CALL LOOP_4(1,2,3,4,DCMPLX(ZERO),DCMPLX(CMASS_MDL_MZ)
     $ ,DCMPLX(ZERO),DCMPLX(CMASS_MDL_MW),2,I_SO,18)
C     CutTools call for loop numbers 21
      CALL LOOP_4(1,2,4,3,DCMPLX(ZERO),DCMPLX(CMASS_MDL_MZ)
     $ ,DCMPLX(ZERO),DCMPLX(CMASS_MDL_MW),2,I_SO,19)
C     CutTools call for loop numbers 22
      CALL LOOP_1_2(2,5,6,DCMPLX(CMASS_MDL_MH),0,I_SO,20)
C     CutTools call for loop numbers 23
      CALL LOOP_2(5,6,DCMPLX(CMASS_MDL_MW),DCMPLX(CMASS_MDL_MH),2,I_SO
     $ ,21)
C     CutTools call for loop numbers 24
      CALL LOOP_2(5,6,DCMPLX(CMASS_MDL_MH),DCMPLX(CMASS_MDL_MW),0,I_SO
     $ ,22)
C     CutTools call for loop numbers 27,28
      CALL LOOP_1_2(2,5,6,DCMPLX(CMASS_MDL_MW),0,I_SO,23)
C     CutTools call for loop numbers 33,35
      CALL LOOP_2(5,6,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,24)
C     CutTools call for loop numbers 34
      CALL LOOP_2(5,6,DCMPLX(CMASS_MDL_MT),DCMPLX(ZERO),2,I_SO,25)
C     At this point, all reductions needed for (QCD=0 QED=6), i.e. of
C      split order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 5000

      GOTO 1001
 5000 CONTINUE
      CTCALL_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

