      SUBROUTINE MP_HELAS_CALLS_AMPB_1(P,NHEL,H,IC)
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
      REAL*16     ZERO
      PARAMETER (ZERO=0.0E0_16)
      COMPLEX*32     IZERO
      PARAMETER (IZERO=CMPLX(0.0E0_16,0.0E0_16,KIND=16))
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=1, NSQUAREDSO=1, NAMPSO=2)
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
      COMMON/FILTERS/GOODAMP,GOODHEL

      INTEGER SQSO_TARGET
      COMMON/SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DON
     $ E,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO

      COMPLEX*32 AMP(NBORNAMPS)
      COMMON/MP_AMPS/AMP
      COMPLEX*32 W(20,NWAVEFUNCS)
      COMMON/MP_W/W

      COMPLEX*32 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE,0:NLOOPWAVEF
     $ UNCS)
      COMPLEX*32 PL(0:3,0:NLOOPWAVEFUNCS)
      COMMON/MP_WL/WL,PL

      COMPLEX*32 LOOPCOEFS(0:LOOPMAXCOEFS-1,NSQUAREDSO,NLOOPGROUPS)
      COMMON/MP_LCOEFS/LOOPCOEFS

      COMPLEX*32 AMPL(3,NCTAMPS)
      COMMON/MP_AMPL/AMPL

      COMPLEX*16 LOOPRES(3,NSQUAREDSO,NLOOPGROUPS)
      LOGICAL S(NSQUAREDSO,NLOOPGROUPS)
      COMMON/LOOPRES/LOOPRES,S
C     
C     ----------
C     BEGIN CODE
C     ----------

C     The target squared split order contribution is already reached
C      if true.
      IF (FILTER_SO.AND.MP_CT_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

      CALL MP_VXXXXX(P(0,1),ZERO,NHEL(1),-1*IC(1),W(1,1))
      CALL MP_VXXXXX(P(0,2),ZERO,NHEL(2),-1*IC(2),W(1,2))
      CALL MP_SXXXXX(P(0,3),+1*IC(3),W(1,3))
      CALL MP_VXXXXX(P(0,4),ZERO,NHEL(4),+1*IC(4),W(1,4))
C     Amplitude(s) for born diagram with ID 1
      CALL MP_VVVS1_2_0(W(1,1),W(1,2),W(1,4),W(1,3),GC_3003H,GC_3003A
     $ ,AMP(1))
      CALL MP_VVV1P0_1(W(1,1),W(1,2),GC_4,ZERO,ZERO,W(1,5))
C     Amplitude(s) for born diagram with ID 2
      CALL MP_VVS2_10_0(W(1,4),W(1,5),W(1,3),GC_3002H,GC_3002A,AMP(2))
      CALL MP_VVS2_10P0_1(W(1,1),W(1,3),GC_3002H,GC_3002A,ZERO,ZERO
     $ ,W(1,6))
C     Amplitude(s) for born diagram with ID 3
      CALL MP_VVV1_0(W(1,2),W(1,4),W(1,6),GC_4,AMP(3))
      CALL MP_VVV1P0_1(W(1,1),W(1,4),GC_4,ZERO,ZERO,W(1,7))
C     Amplitude(s) for born diagram with ID 4
      CALL MP_VVS2_10_0(W(1,2),W(1,7),W(1,3),GC_3002H,GC_3002A,AMP(4))
      CALL MP_VVS2_10P0_1(W(1,4),W(1,3),GC_3002H,GC_3002A,ZERO,ZERO
     $ ,W(1,8))
C     Counter-term amplitude(s) for loop diagram number 5
      CALL MP_R2_GG_1_0(W(1,5),W(1,8),R2_GGQ,AMPL(1,1))
      CALL MP_VVV1P0_1(W(1,2),W(1,4),GC_4,ZERO,ZERO,W(1,9))
C     Counter-term amplitude(s) for loop diagram number 6
      CALL MP_R2_GG_1_0(W(1,6),W(1,9),R2_GGQ,AMPL(1,2))
C     Counter-term amplitude(s) for loop diagram number 7
      CALL MP_VVV1_0(W(1,2),W(1,4),W(1,6),UV_3GB_1EPS,AMPL(2,3))
      CALL MP_VVV1_0(W(1,2),W(1,4),W(1,6),UV_3GB_1EPS,AMPL(2,4))
      CALL MP_VVV1_0(W(1,2),W(1,4),W(1,6),UV_3GB_1EPS,AMPL(2,5))
      CALL MP_VVV1_0(W(1,2),W(1,4),W(1,6),UV_3GB_1EPS,AMPL(2,6))
      CALL MP_VVV1_0(W(1,2),W(1,4),W(1,6),UV_3GB_1EPS,AMPL(2,7))
      CALL MP_VVV1_0(W(1,2),W(1,4),W(1,6),UV_3GG_1EPS,AMPL(2,8))
      CALL MP_VVV1_0(W(1,2),W(1,4),W(1,6),R2_3GQ,AMPL(1,9))
      CALL MP_VVS2_10P0_1(W(1,2),W(1,3),GC_3002H,GC_3002A,ZERO,ZERO
     $ ,W(1,10))
C     Counter-term amplitude(s) for loop diagram number 9
      CALL MP_R2_GG_1_0(W(1,10),W(1,7),R2_GGQ,AMPL(1,10))
C     Counter-term amplitude(s) for loop diagram number 10
      CALL MP_VVV1_0(W(1,1),W(1,4),W(1,10),UV_3GB_1EPS,AMPL(2,11))
      CALL MP_VVV1_0(W(1,1),W(1,4),W(1,10),UV_3GB_1EPS,AMPL(2,12))
      CALL MP_VVV1_0(W(1,1),W(1,4),W(1,10),UV_3GB_1EPS,AMPL(2,13))
      CALL MP_VVV1_0(W(1,1),W(1,4),W(1,10),UV_3GB_1EPS,AMPL(2,14))
      CALL MP_VVV1_0(W(1,1),W(1,4),W(1,10),UV_3GB_1EPS,AMPL(2,15))
      CALL MP_VVV1_0(W(1,1),W(1,4),W(1,10),UV_3GG_1EPS,AMPL(2,16))
      CALL MP_VVV1_0(W(1,1),W(1,4),W(1,10),R2_3GQ,AMPL(1,17))
C     Counter-term amplitude(s) for loop diagram number 12
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),UV_3GB_1EPS,AMPL(2,18))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),UV_3GB_1EPS,AMPL(2,19))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),UV_3GB_1EPS,AMPL(2,20))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),UV_3GB_1EPS,AMPL(2,21))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),UV_3GB_1EPS,AMPL(2,22))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),UV_3GG_1EPS,AMPL(2,23))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),R2_3GQ,AMPL(1,24))
C     Counter-term amplitude(s) for loop diagram number 14
      CALL MP_R2_GG_1_0(W(1,5),W(1,8),R2_GGQ,AMPL(1,25))
C     Counter-term amplitude(s) for loop diagram number 15
      CALL MP_R2_GG_1_0(W(1,6),W(1,9),R2_GGQ,AMPL(1,26))
C     Counter-term amplitude(s) for loop diagram number 16
      CALL MP_VVV1_0(W(1,2),W(1,4),W(1,6),R2_3GQ,AMPL(1,27))
C     Counter-term amplitude(s) for loop diagram number 18
      CALL MP_R2_GG_1_0(W(1,10),W(1,7),R2_GGQ,AMPL(1,28))
C     Counter-term amplitude(s) for loop diagram number 19
      CALL MP_VVV1_0(W(1,1),W(1,4),W(1,10),R2_3GQ,AMPL(1,29))
C     Counter-term amplitude(s) for loop diagram number 21
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),R2_3GQ,AMPL(1,30))
C     Counter-term amplitude(s) for loop diagram number 23
      CALL MP_R2_GG_1_0(W(1,5),W(1,8),R2_GGQ,AMPL(1,31))
C     Counter-term amplitude(s) for loop diagram number 24
      CALL MP_R2_GG_1_0(W(1,6),W(1,9),R2_GGQ,AMPL(1,32))
C     Counter-term amplitude(s) for loop diagram number 25
      CALL MP_VVV1_0(W(1,2),W(1,4),W(1,6),R2_3GQ,AMPL(1,33))
C     Counter-term amplitude(s) for loop diagram number 27
      CALL MP_R2_GG_1_0(W(1,10),W(1,7),R2_GGQ,AMPL(1,34))
C     Counter-term amplitude(s) for loop diagram number 28
      CALL MP_VVV1_0(W(1,1),W(1,4),W(1,10),R2_3GQ,AMPL(1,35))
C     Counter-term amplitude(s) for loop diagram number 30
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),R2_3GQ,AMPL(1,36))
C     Counter-term amplitude(s) for loop diagram number 32
      CALL MP_R2_GG_1_0(W(1,5),W(1,8),R2_GGQ,AMPL(1,37))
C     Counter-term amplitude(s) for loop diagram number 33
      CALL MP_R2_GG_1_0(W(1,6),W(1,9),R2_GGQ,AMPL(1,38))
C     Counter-term amplitude(s) for loop diagram number 34
      CALL MP_VVV1_0(W(1,2),W(1,4),W(1,6),R2_3GQ,AMPL(1,39))
C     Counter-term amplitude(s) for loop diagram number 36
      CALL MP_R2_GG_1_0(W(1,10),W(1,7),R2_GGQ,AMPL(1,40))
C     Counter-term amplitude(s) for loop diagram number 37
      CALL MP_VVV1_0(W(1,1),W(1,4),W(1,10),R2_3GQ,AMPL(1,41))
C     Counter-term amplitude(s) for loop diagram number 39
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),R2_3GQ,AMPL(1,42))
C     Counter-term amplitude(s) for loop diagram number 41
      CALL MP_R2_GG_1_0(W(1,5),W(1,8),R2_GGQ,AMPL(1,43))
C     Counter-term amplitude(s) for loop diagram number 42
      CALL MP_R2_GG_1_0(W(1,6),W(1,9),R2_GGQ,AMPL(1,44))
C     Counter-term amplitude(s) for loop diagram number 43
      CALL MP_VVV1_0(W(1,2),W(1,4),W(1,6),R2_3GQ,AMPL(1,45))
C     Counter-term amplitude(s) for loop diagram number 45
      CALL MP_R2_GG_1_0(W(1,10),W(1,7),R2_GGQ,AMPL(1,46))
C     Counter-term amplitude(s) for loop diagram number 46
      CALL MP_VVV1_0(W(1,1),W(1,4),W(1,10),R2_3GQ,AMPL(1,47))
C     Counter-term amplitude(s) for loop diagram number 48
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),R2_3GQ,AMPL(1,48))
C     Counter-term amplitude(s) for loop diagram number 50
      CALL MP_R2_GG_1_R2_GG_2_0(W(1,5),W(1,8),R2_GGG_1,R2_GGG_2,AMPL(1
     $ ,49))
C     Counter-term amplitude(s) for loop diagram number 51
      CALL MP_VVS2_10_0(W(1,5),W(1,4),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC
     $ _1EPS,AMPL(2,50))
      CALL MP_VVS2_10_0(W(1,5),W(1,4),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC
     $ _1EPS,AMPL(2,51))
      CALL MP_VVS2_10_0(W(1,5),W(1,4),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC
     $ _1EPS,AMPL(2,52))
      CALL MP_VVS2_10_0(W(1,5),W(1,4),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC
     $ _1EPS,AMPL(2,53))
      CALL MP_VVS2_10_0(W(1,5),W(1,4),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC
     $ _1EPS,AMPL(2,54))
      CALL MP_VVS2_0(W(1,5),W(1,4),W(1,3),UV_GGX0_HG,AMPL(1,55))
      CALL MP_VVS2_10_0(W(1,5),W(1,4),W(1,3),UV_GGX0_HG_1EPS,UV_GGX0_AG
     $ _1EPS,AMPL(2,56))
      CALL MP_VVS3_10_0(W(1,5),W(1,4),W(1,3),R2_GGX0_H,R2_GGX0_A
     $ ,AMPL(1,57))
C     Counter-term amplitude(s) for loop diagram number 55
      CALL MP_R2_GG_1_R2_GG_2_0(W(1,6),W(1,9),R2_GGG_1,R2_GGG_2,AMPL(1
     $ ,58))
C     Counter-term amplitude(s) for loop diagram number 56
      CALL MP_VVV1_0(W(1,2),W(1,6),W(1,4),R2_3GG,AMPL(1,59))
C     Counter-term amplitude(s) for loop diagram number 60
      CALL MP_R2_GG_1_R2_GG_2_0(W(1,10),W(1,7),R2_GGG_1,R2_GGG_2
     $ ,AMPL(1,60))
C     Counter-term amplitude(s) for loop diagram number 61
      CALL MP_VVS2_10_0(W(1,2),W(1,7),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC
     $ _1EPS,AMPL(2,61))
      CALL MP_VVS2_10_0(W(1,2),W(1,7),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC
     $ _1EPS,AMPL(2,62))
      CALL MP_VVS2_10_0(W(1,2),W(1,7),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC
     $ _1EPS,AMPL(2,63))
      CALL MP_VVS2_10_0(W(1,2),W(1,7),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC
     $ _1EPS,AMPL(2,64))
      CALL MP_VVS2_10_0(W(1,2),W(1,7),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC
     $ _1EPS,AMPL(2,65))
      CALL MP_VVS2_0(W(1,2),W(1,7),W(1,3),UV_GGX0_HG,AMPL(1,66))
      CALL MP_VVS2_10_0(W(1,2),W(1,7),W(1,3),UV_GGX0_HG_1EPS,UV_GGX0_AG
     $ _1EPS,AMPL(2,67))
      CALL MP_VVS3_10_0(W(1,2),W(1,7),W(1,3),R2_GGX0_H,R2_GGX0_A
     $ ,AMPL(1,68))
C     Counter-term amplitude(s) for loop diagram number 65
      CALL MP_VVV1_0(W(1,1),W(1,4),W(1,10),R2_3GG,AMPL(1,69))
C     Counter-term amplitude(s) for loop diagram number 67
      CALL MP_VVS2_10_0(W(1,1),W(1,9),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC
     $ _1EPS,AMPL(2,70))
      CALL MP_VVS2_10_0(W(1,1),W(1,9),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC
     $ _1EPS,AMPL(2,71))
      CALL MP_VVS2_10_0(W(1,1),W(1,9),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC
     $ _1EPS,AMPL(2,72))
      CALL MP_VVS2_10_0(W(1,1),W(1,9),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC
     $ _1EPS,AMPL(2,73))
      CALL MP_VVS2_10_0(W(1,1),W(1,9),W(1,3),UV_GGX0_HB_1EPS,UV_GGX0_AC
     $ _1EPS,AMPL(2,74))
      CALL MP_VVS2_0(W(1,1),W(1,9),W(1,3),UV_GGX0_HG,AMPL(1,75))
      CALL MP_VVS2_10_0(W(1,1),W(1,9),W(1,3),UV_GGX0_HG_1EPS,UV_GGX0_AG
     $ _1EPS,AMPL(2,76))
      CALL MP_VVS3_10_0(W(1,1),W(1,9),W(1,3),R2_GGX0_H,R2_GGX0_A
     $ ,AMPL(1,77))
C     Counter-term amplitude(s) for loop diagram number 69
      CALL MP_VVVS1_2_0(W(1,1),W(1,2),W(1,4),W(1,3),UV_3GX0_HB_1EPS
     $ ,UV_3GX0_AT_1EPS,AMPL(2,78))
      CALL MP_VVVS1_2_0(W(1,1),W(1,2),W(1,4),W(1,3),UV_3GX0_HB_1EPS
     $ ,UV_3GX0_AT_1EPS,AMPL(2,79))
      CALL MP_VVVS1_2_0(W(1,1),W(1,2),W(1,4),W(1,3),UV_3GX0_HB_1EPS
     $ ,UV_3GX0_AT_1EPS,AMPL(2,80))
      CALL MP_VVVS1_2_0(W(1,1),W(1,2),W(1,4),W(1,3),UV_3GX0_HB_1EPS
     $ ,UV_3GX0_AT_1EPS,AMPL(2,81))
      CALL MP_VVVS1_2_0(W(1,1),W(1,2),W(1,4),W(1,3),UV_3GX0_HB_1EPS
     $ ,UV_3GX0_AT_1EPS,AMPL(2,82))
      CALL MP_VVVS1_0(W(1,1),W(1,2),W(1,4),W(1,3),UV_3GX0_HG,AMPL(1
     $ ,83))
      CALL MP_VVVS1_2_0(W(1,1),W(1,2),W(1,4),W(1,3),UV_3GX0_HG_1EPS
     $ ,UV_3GX0_AG_1EPS,AMPL(2,84))
      CALL MP_VVVS1_2_0(W(1,1),W(1,2),W(1,4),W(1,3),R2_3GX0_H
     $ ,R2_3GX0_A,AMPL(1,85))
C     Counter-term amplitude(s) for loop diagram number 72
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),R2_3GG,AMPL(1,86))
C     At this point, all CT amps needed for (QCD=4), i.e. of split
C      order ID=0, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 2000

      GOTO 1001
 2000 CONTINUE
      MP_CT_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

