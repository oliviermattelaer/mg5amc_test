      SUBROUTINE HELAS_CALLS_AMPB_1(P,NHEL,H,IC)
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=4)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=4)
      INTEGER    NCOMB
      PARAMETER (NCOMB=8)
      INTEGER    NLOOPS, NLOOPGROUPS, NCTAMPS
      PARAMETER (NLOOPS=138, NLOOPGROUPS=16, NCTAMPS=86)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=10,NLOOPWAVEFUNCS=237)
      INTEGER MAXLWFSIZE
      PARAMETER (MAXLWFSIZE=4)
      INTEGER LOOPMAXCOEFS, VERTEXMAXCOEFS
      PARAMETER (LOOPMAXCOEFS=126, VERTEXMAXCOEFS=15)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      REAL*16     MP__ZERO
      PARAMETER (MP__ZERO=0.0E0_16)
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=1, NSQUAREDSO=1, NAMPSO=2)
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
      COMMON/SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DON
     $ E,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO

      INTEGER I_SO
      COMMON/I_SO/I_SO
      INTEGER I_LIB
      COMMON/I_LIB/I_LIB

      COMPLEX*16 AMP(NBORNAMPS)
      COMMON/AMPS/AMP
      COMPLEX*16 W(20,NWAVEFUNCS)
      COMMON/W/W

      COMPLEX*16 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE,0:NLOOPWAVEF
     $ UNCS)
      COMPLEX*16 PL(0:3,0:NLOOPWAVEFUNCS)
      COMMON/WL/WL,PL

      COMPLEX*16 LOOPCOEFS(0:LOOPMAXCOEFS-1,NSQUAREDSO,NLOOPGROUPS)
      COMMON/LCOEFS/LOOPCOEFS

      COMPLEX*16 AMPL(3,NCTAMPS)
      COMMON/AMPL/AMPL

      COMPLEX*16 LOOPRES(3,NSQUAREDSO,NLOOPGROUPS)
      LOGICAL S(NSQUAREDSO,NLOOPGROUPS)
      COMMON/LOOPRES/LOOPRES,S
C     
C     ----------
C     BEGIN CODE
C     ----------

C     The target squared split order contribution is already reached
C      if true.
      IF (FILTER_SO.AND.CT_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

      CALL VXXXXX(P(0,1),ZERO,NHEL(1),-1*IC(1),W(1,1))
      CALL VXXXXX(P(0,2),ZERO,NHEL(2),-1*IC(2),W(1,2))
      CALL SXXXXX(P(0,3),+1*IC(3),W(1,3))
      CALL VXXXXX(P(0,4),ZERO,NHEL(4),+1*IC(4),W(1,4))
C     Amplitude(s) for born diagram with ID 1
      CALL VVVS1_2_0(W(1,1),W(1,2),W(1,4),W(1,3),GC_3003H,GC_3003A
     $ ,AMP(1))
      CALL VVV1P0_1(W(1,1),W(1,2),GC_4,ZERO,ZERO,W(1,5))
C     Amplitude(s) for born diagram with ID 2
      CALL VVS2_10_0(W(1,4),W(1,5),W(1,3),GC_3002H,GC_3002A,AMP(2))
      CALL VVS2_10P0_1(W(1,1),W(1,3),GC_3002H,GC_3002A,ZERO,ZERO,W(1
     $ ,6))
C     Amplitude(s) for born diagram with ID 3
      CALL VVV1_0(W(1,2),W(1,4),W(1,6),GC_4,AMP(3))
      CALL VVV1P0_1(W(1,1),W(1,4),GC_4,ZERO,ZERO,W(1,7))
C     Amplitude(s) for born diagram with ID 4
      CALL VVS2_10_0(W(1,2),W(1,7),W(1,3),GC_3002H,GC_3002A,AMP(4))
      CALL VVS2_10P0_1(W(1,4),W(1,3),GC_3002H,GC_3002A,ZERO,ZERO,W(1
     $ ,8))
C     Counter-term amplitude(s) for loop diagram number 5
      CALL R2_GG_1_0(W(1,5),W(1,8),R2_GGQ,AMPL(1,1))
      CALL VVV1P0_1(W(1,2),W(1,4),GC_4,ZERO,ZERO,W(1,9))
C     Counter-term amplitude(s) for loop diagram number 6
      CALL R2_GG_1_0(W(1,6),W(1,9),R2_GGQ,AMPL(1,2))
C     Counter-term amplitude(s) for loop diagram number 7
      CALL VVV1_0(W(1,2),W(1,4),W(1,6),UV_3GB_1EPS,AMPL(2,3))
      CALL VVV1_0(W(1,2),W(1,4),W(1,6),UV_3GB_1EPS,AMPL(2,4))
      CALL VVV1_0(W(1,2),W(1,4),W(1,6),UV_3GB_1EPS,AMPL(2,5))
      CALL VVV1_0(W(1,2),W(1,4),W(1,6),UV_3GB_1EPS,AMPL(2,6))
      CALL VVV1_0(W(1,2),W(1,4),W(1,6),UV_3GB_1EPS,AMPL(2,7))
      CALL VVV1_0(W(1,2),W(1,4),W(1,6),UV_3GG_1EPS,AMPL(2,8))
      CALL VVV1_0(W(1,2),W(1,4),W(1,6),R2_3GQ,AMPL(1,9))
      CALL VVS2_10P0_1(W(1,2),W(1,3),GC_3002H,GC_3002A,ZERO,ZERO,W(1
     $ ,10))
C     Counter-term amplitude(s) for loop diagram number 9
      CALL R2_GG_1_0(W(1,10),W(1,7),R2_GGQ,AMPL(1,10))
C     Counter-term amplitude(s) for loop diagram number 10
      CALL VVV1_0(W(1,1),W(1,4),W(1,10),UV_3GB_1EPS,AMPL(2,11))
      CALL VVV1_0(W(1,1),W(1,4),W(1,10),UV_3GB_1EPS,AMPL(2,12))
      CALL VVV1_0(W(1,1),W(1,4),W(1,10),UV_3GB_1EPS,AMPL(2,13))
      CALL VVV1_0(W(1,1),W(1,4),W(1,10),UV_3GB_1EPS,AMPL(2,14))
      CALL VVV1_0(W(1,1),W(1,4),W(1,10),UV_3GB_1EPS,AMPL(2,15))
      CALL VVV1_0(W(1,1),W(1,4),W(1,10),UV_3GG_1EPS,AMPL(2,16))
      CALL VVV1_0(W(1,1),W(1,4),W(1,10),R2_3GQ,AMPL(1,17))
C     Counter-term amplitude(s) for loop diagram number 12
      CALL VVV1_0(W(1,1),W(1,2),W(1,8),UV_3GB_1EPS,AMPL(2,18))
      CALL VVV1_0(W(1,1),W(1,2),W(1,8),UV_3GB_1EPS,AMPL(2,19))
      CALL VVV1_0(W(1,1),W(1,2),W(1,8),UV_3GB_1EPS,AMPL(2,20))
      CALL VVV1_0(W(1,1),W(1,2),W(1,8),UV_3GB_1EPS,AMPL(2,21))
      CALL VVV1_0(W(1,1),W(1,2),W(1,8),UV_3GB_1EPS,AMPL(2,22))
      CALL VVV1_0(W(1,1),W(1,2),W(1,8),UV_3GG_1EPS,AMPL(2,23))
      CALL VVV1_0(W(1,1),W(1,2),W(1,8),R2_3GQ,AMPL(1,24))
C     Counter-term amplitude(s) for loop diagram number 14
      CALL R2_GG_1_0(W(1,5),W(1,8),R2_GGQ,AMPL(1,25))
C     Counter-term amplitude(s) for loop diagram number 15
      CALL R2_GG_1_0(W(1,6),W(1,9),R2_GGQ,AMPL(1,26))
C     Counter-term amplitude(s) for loop diagram number 16
      CALL VVV1_0(W(1,2),W(1,4),W(1,6),R2_3GQ,AMPL(1,27))
C     Counter-term amplitude(s) for loop diagram number 18
      CALL R2_GG_1_0(W(1,10),W(1,7),R2_GGQ,AMPL(1,28))
C     Counter-term amplitude(s) for loop diagram number 19
      CALL VVV1_0(W(1,1),W(1,4),W(1,10),R2_3GQ,AMPL(1,29))
C     Counter-term amplitude(s) for loop diagram number 21
      CALL VVV1_0(W(1,1),W(1,2),W(1,8),R2_3GQ,AMPL(1,30))
C     Counter-term amplitude(s) for loop diagram number 23
      CALL R2_GG_1_0(W(1,5),W(1,8),R2_GGQ,AMPL(1,31))
C     Counter-term amplitude(s) for loop diagram number 24
      CALL R2_GG_1_0(W(1,6),W(1,9),R2_GGQ,AMPL(1,32))
C     Counter-term amplitude(s) for loop diagram number 25
      CALL VVV1_0(W(1,2),W(1,4),W(1,6),R2_3GQ,AMPL(1,33))
C     Counter-term amplitude(s) for loop diagram number 27
      CALL R2_GG_1_0(W(1,10),W(1,7),R2_GGQ,AMPL(1,34))
C     Counter-term amplitude(s) for loop diagram number 28
      CALL VVV1_0(W(1,1),W(1,4),W(1,10),R2_3GQ,AMPL(1,35))
C     Counter-term amplitude(s) for loop diagram number 30
      CALL VVV1_0(W(1,1),W(1,2),W(1,8),R2_3GQ,AMPL(1,36))
C     Counter-term amplitude(s) for loop diagram number 32
      CALL R2_GG_1_0(W(1,5),W(1,8),R2_GGQ,AMPL(1,37))
C     Counter-term amplitude(s) for loop diagram number 33
      CALL R2_GG_1_0(W(1,6),W(1,9),R2_GGQ,AMPL(1,38))
C     Counter-term amplitude(s) for loop diagram number 34
      CALL VVV1_0(W(1,2),W(1,4),W(1,6),R2_3GQ,AMPL(1,39))
C     Counter-term amplitude(s) for loop diagram number 36
      CALL R2_GG_1_0(W(1,10),W(1,7),R2_GGQ,AMPL(1,40))
C     Counter-term amplitude(s) for loop diagram number 37
      CALL VVV1_0(W(1,1),W(1,4),W(1,10),R2_3GQ,AMPL(1,41))
C     Counter-term amplitude(s) for loop diagram number 39
      CALL VVV1_0(W(1,1),W(1,2),W(1,8),R2_3GQ,AMPL(1,42))
C     Counter-term amplitude(s) for loop diagram number 41
      CALL R2_GG_1_0(W(1,5),W(1,8),R2_GGQ,AMPL(1,43))
C     Counter-term amplitude(s) for loop diagram number 42
      CALL R2_GG_1_0(W(1,6),W(1,9),R2_GGQ,AMPL(1,44))
C     Counter-term amplitude(s) for loop diagram number 43
      CALL VVV1_0(W(1,2),W(1,4),W(1,6),R2_3GQ,AMPL(1,45))
C     Counter-term amplitude(s) for loop diagram number 45
      CALL R2_GG_1_0(W(1,10),W(1,7),R2_GGQ,AMPL(1,46))
C     Counter-term amplitude(s) for loop diagram number 46
      CALL VVV1_0(W(1,1),W(1,4),W(1,10),R2_3GQ,AMPL(1,47))
C     Counter-term amplitude(s) for loop diagram number 48
      CALL VVV1_0(W(1,1),W(1,2),W(1,8),R2_3GQ,AMPL(1,48))
C     Counter-term amplitude(s) for loop diagram number 50
      CALL R2_GG_1_R2_GG_2_0(W(1,5),W(1,8),R2_GGG_1,R2_GGG_2,AMPL(1
     $ ,49))
C     Counter-term amplitude(s) for loop diagram number 51
      CALL VVS2_10_0(W(1,5),W(1,4),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC_1E
     $ PS,AMPL(2,50))
      CALL VVS2_10_0(W(1,5),W(1,4),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC_1E
     $ PS,AMPL(2,51))
      CALL VVS2_10_0(W(1,5),W(1,4),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC_1E
     $ PS,AMPL(2,52))
      CALL VVS2_10_0(W(1,5),W(1,4),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC_1E
     $ PS,AMPL(2,53))
      CALL VVS2_10_0(W(1,5),W(1,4),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC_1E
     $ PS,AMPL(2,54))
      CALL VVS2_0(W(1,5),W(1,4),W(1,3),UV_GGX0_HG,AMPL(1,55))
      CALL VVS2_10_0(W(1,5),W(1,4),W(1,3),UV_GGX0_HG_1EPS,UV_GGX0_AG_1E
     $ PS,AMPL(2,56))
      CALL VVS3_10_0(W(1,5),W(1,4),W(1,3),R2_GGX0_H,R2_GGX0_A,AMPL(1
     $ ,57))
C     Counter-term amplitude(s) for loop diagram number 55
      CALL R2_GG_1_R2_GG_2_0(W(1,6),W(1,9),R2_GGG_1,R2_GGG_2,AMPL(1
     $ ,58))
C     Counter-term amplitude(s) for loop diagram number 56
      CALL VVV1_0(W(1,2),W(1,6),W(1,4),R2_3GG,AMPL(1,59))
C     Counter-term amplitude(s) for loop diagram number 60
      CALL R2_GG_1_R2_GG_2_0(W(1,10),W(1,7),R2_GGG_1,R2_GGG_2,AMPL(1
     $ ,60))
C     Counter-term amplitude(s) for loop diagram number 61
      CALL VVS2_10_0(W(1,2),W(1,7),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC_1E
     $ PS,AMPL(2,61))
      CALL VVS2_10_0(W(1,2),W(1,7),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC_1E
     $ PS,AMPL(2,62))
      CALL VVS2_10_0(W(1,2),W(1,7),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC_1E
     $ PS,AMPL(2,63))
      CALL VVS2_10_0(W(1,2),W(1,7),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC_1E
     $ PS,AMPL(2,64))
      CALL VVS2_10_0(W(1,2),W(1,7),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC_1E
     $ PS,AMPL(2,65))
      CALL VVS2_0(W(1,2),W(1,7),W(1,3),UV_GGX0_HG,AMPL(1,66))
      CALL VVS2_10_0(W(1,2),W(1,7),W(1,3),UV_GGX0_HG_1EPS,UV_GGX0_AG_1E
     $ PS,AMPL(2,67))
      CALL VVS3_10_0(W(1,2),W(1,7),W(1,3),R2_GGX0_H,R2_GGX0_A,AMPL(1
     $ ,68))
C     Counter-term amplitude(s) for loop diagram number 65
      CALL VVV1_0(W(1,1),W(1,4),W(1,10),R2_3GG,AMPL(1,69))
C     Counter-term amplitude(s) for loop diagram number 67
      CALL VVS2_10_0(W(1,1),W(1,9),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC_1E
     $ PS,AMPL(2,70))
      CALL VVS2_10_0(W(1,1),W(1,9),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC_1E
     $ PS,AMPL(2,71))
      CALL VVS2_10_0(W(1,1),W(1,9),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC_1E
     $ PS,AMPL(2,72))
      CALL VVS2_10_0(W(1,1),W(1,9),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC_1E
     $ PS,AMPL(2,73))
      CALL VVS2_10_0(W(1,1),W(1,9),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC_1E
     $ PS,AMPL(2,74))
      CALL VVS2_0(W(1,1),W(1,9),W(1,3),UV_GGX0_HG,AMPL(1,75))
      CALL VVS2_10_0(W(1,1),W(1,9),W(1,3),UV_GGX0_HG_1EPS,UV_GGX0_AG_1E
     $ PS,AMPL(2,76))
      CALL VVS3_10_0(W(1,1),W(1,9),W(1,3),R2_GGX0_H,R2_GGX0_A,AMPL(1
     $ ,77))
C     Counter-term amplitude(s) for loop diagram number 69
      CALL VVVS1_2_0(W(1,1),W(1,2),W(1,4),W(1,3),UV_3GX0_HB_1EPS
     $ ,UV_3GX0_AT_1EPS,AMPL(2,78))
      CALL VVVS1_2_0(W(1,1),W(1,2),W(1,4),W(1,3),UV_3GX0_HB_1EPS
     $ ,UV_3GX0_AT_1EPS,AMPL(2,79))
      CALL VVVS1_2_0(W(1,1),W(1,2),W(1,4),W(1,3),UV_3GX0_HB_1EPS
     $ ,UV_3GX0_AT_1EPS,AMPL(2,80))
      CALL VVVS1_2_0(W(1,1),W(1,2),W(1,4),W(1,3),UV_3GX0_HB_1EPS
     $ ,UV_3GX0_AT_1EPS,AMPL(2,81))
      CALL VVVS1_2_0(W(1,1),W(1,2),W(1,4),W(1,3),UV_3GX0_HB_1EPS
     $ ,UV_3GX0_AT_1EPS,AMPL(2,82))
      CALL VVVS1_0(W(1,1),W(1,2),W(1,4),W(1,3),UV_3GX0_HG,AMPL(1,83))
      CALL VVVS1_2_0(W(1,1),W(1,2),W(1,4),W(1,3),UV_3GX0_HG_1EPS
     $ ,UV_3GX0_AG_1EPS,AMPL(2,84))
      CALL VVVS1_2_0(W(1,1),W(1,2),W(1,4),W(1,3),R2_3GX0_H,R2_3GX0_A
     $ ,AMPL(1,85))
C     Counter-term amplitude(s) for loop diagram number 72
      CALL VVV1_0(W(1,1),W(1,2),W(1,8),R2_3GG,AMPL(1,86))
C     At this point, all CT amps needed for (QCD=4), i.e. of split
C      order ID=0, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 2000

      GOTO 1001
 2000 CONTINUE
      CT_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

