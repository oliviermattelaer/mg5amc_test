      SUBROUTINE ML5_0_MP_HELAS_CALLS_UVCT_1(P,NHEL,H,IC)
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
      PARAMETER (NBORNAMPS=3)
      INTEGER    NLOOPS, NLOOPGROUPS, NCTAMPS
      PARAMETER (NLOOPS=44, NLOOPGROUPS=26, NCTAMPS=85)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=129)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=10,NLOOPWAVEFUNCS=93)
      INCLUDE 'loop_max_coefs.inc'
      INCLUDE 'coef_specs.inc'
      REAL*16     ZERO
      PARAMETER (ZERO=0.0E0_16)
      COMPLEX*32     IZERO
      PARAMETER (IZERO=CMPLX(0.0E0_16,0.0E0_16,KIND=16))
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=0, NSQUAREDSO=0, NAMPSO=0)
C     
C     ARGUMENTS
C     
      REAL*16 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      INTEGER H
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K
      COMPLEX*32 COEFS(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'mp_coupl_same_name.inc'

      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/ML5_0_FILTERS/GOODAMP,GOODHEL

      INTEGER SQSO_TARGET
      COMMON/ML5_0_SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/ML5_0_SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_S
     $ O_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

      COMPLEX*32 AMP(NBORNAMPS)
      COMMON/ML5_0_MP_AMPS/AMP
      COMPLEX*32 W(20,NWAVEFUNCS)
      COMMON/ML5_0_MP_W/W

      COMPLEX*32 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE,0:NLOOPWAVEF
     $ UNCS)
      COMPLEX*32 PL(0:3,0:NLOOPWAVEFUNCS)
      COMMON/ML5_0_MP_WL/WL,PL

      COMPLEX*32 AMPL(3,NCTAMPS)
      COMMON/ML5_0_MP_AMPL/AMPL

C     
C     ----------
C     BEGIN CODE
C     ----------

C     The target squared split order contribution is already reached
C      if true.
      IF (FILTER_SO.AND.MP_UVCT_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

C     Amplitude(s) for UVCT diagram with ID 40
      CALL MP_FFV1_0(W(1,4),W(1,3),W(1,5),GC_5,AMPL(1,80))
      AMPL(1,80)=AMPL(1,80)*(2.0D0*UVWFCT_G_2+2.0D0*UVWFCT_G_1
     $ +2.0D0*UVWFCT_T_0)
C     Amplitude(s) for UVCT diagram with ID 41
      CALL MP_FFV1_0(W(1,4),W(1,3),W(1,5),GC_5,AMPL(2,81))
      AMPL(2,81)=AMPL(2,81)*(2.0D0*UVWFCT_B_0_1EPS+4.0D0*UVWFCT_G_2_1EP
     $ S)
C     Amplitude(s) for UVCT diagram with ID 42
      CALL MP_FFV1_0(W(1,4),W(1,6),W(1,2),GC_5,AMPL(1,82))
      AMPL(1,82)=AMPL(1,82)*(2.0D0*UVWFCT_G_2+2.0D0*UVWFCT_G_1
     $ +2.0D0*UVWFCT_T_0)
C     Amplitude(s) for UVCT diagram with ID 43
      CALL MP_FFV1_0(W(1,4),W(1,6),W(1,2),GC_5,AMPL(2,83))
      AMPL(2,83)=AMPL(2,83)*(2.0D0*UVWFCT_B_0_1EPS+4.0D0*UVWFCT_G_2_1EP
     $ S)
C     Amplitude(s) for UVCT diagram with ID 44
      CALL MP_FFV1_0(W(1,7),W(1,3),W(1,2),GC_5,AMPL(1,84))
      AMPL(1,84)=AMPL(1,84)*(2.0D0*UVWFCT_G_2+2.0D0*UVWFCT_G_1
     $ +2.0D0*UVWFCT_T_0)
C     Amplitude(s) for UVCT diagram with ID 45
      CALL MP_FFV1_0(W(1,7),W(1,3),W(1,2),GC_5,AMPL(2,85))
      AMPL(2,85)=AMPL(2,85)*(2.0D0*UVWFCT_B_0_1EPS+4.0D0*UVWFCT_G_2_1EP
     $ S)

      GOTO 1001
 3000 CONTINUE
      MP_UVCT_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

