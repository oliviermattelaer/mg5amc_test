      SUBROUTINE ML5_0_MP_BORN_AMPS_AND_WFS(P)
C     
C     Generated by MadGraph5_aMC@NLO v. %(version)s, %(date)s
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     Computes all the AMP and WFS in quadruple precision for the 
C     phase space point P(0:3,NEXTERNAL)
C     
C     Process: d d~ > t t~ QCD<=2 QED=0 [ virt = QCD ]
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=1)
      INTEGER    NLOOPAMPS, NCTAMPS
      PARAMETER (NLOOPAMPS=40, NCTAMPS=29)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=4)
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=6)
      INTEGER    NCOMB
      PARAMETER (NCOMB=16)
      REAL*16     ZERO
      PARAMETER (ZERO=0E0_16)
      COMPLEX*32 IMAG1
      PARAMETER (IMAG1=(0E0_16,1E0_16))

C     
C     ARGUMENTS 
C     
      REAL*16 P(0:3,NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,H
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      DATA IC/NEXTERNAL*1/
C     
C     FUNCTIONS
C     
      LOGICAL ML5_0_IS_HEL_SELECTED
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'mp_coupl_same_name.inc'

      INTEGER NTRY
      LOGICAL CHECKPHASE,HELDOUBLECHECKED
      REAL*8 REF
      COMMON/ML5_0_INIT/NTRY,CHECKPHASE,HELDOUBLECHECKED,REF

      LOGICAL GOODHEL(NCOMB)
      LOGICAL GOODAMP(NLOOPAMPS,NCOMB)
      COMMON/ML5_0_FILTERS/GOODAMP,GOODHEL

      INTEGER HELPICKED
      COMMON/ML5_0_HELCHOICE/HELPICKED

      COMPLEX*32 AMP(NBORNAMPS,NCOMB)
      COMMON/ML5_0_MP_AMPS/AMP
      COMPLEX*16 DPAMP(NBORNAMPS,NCOMB)
      COMMON/ML5_0_AMPS/DPAMP
      COMPLEX*32 W(20,NWAVEFUNCS,NCOMB)
      COMMON/ML5_0_MP_WFS/W

      COMPLEX*32 AMPL(3,NCTAMPS)
      COMMON/ML5_0_MP_AMPL/AMPL

      COMPLEX*16 DPW(20,NWAVEFUNCS,NCOMB)
      COMMON/ML5_0_WFCTS/DPW

      COMPLEX*16 DPAMPL(3,NLOOPAMPS)
      LOGICAL S(NLOOPAMPS)
      COMMON/ML5_0_AMPL/DPAMPL,S

      INTEGER HELC(NEXTERNAL,NCOMB)
      COMMON/ML5_0_HELCONFIGS/HELC

      LOGICAL MP_DONE_ONCE
      COMMON/ML5_0_MP_DONE_ONCE/MP_DONE_ONCE

C     This array specify potential special requirements on the
C      helicities to
C     consider. POLARIZATIONS(0,0) is -1 if there is not such
C      requirement.
      INTEGER POLARIZATIONS(0:NEXTERNAL,0:5)
      COMMON/ML5_0_BEAM_POL/POLARIZATIONS

C     ----------
C     BEGIN CODE
C     ---------

      MP_DONE_ONCE=.TRUE.

C     To be on the safe side, we always update the MP params here.
C     It can be redundant as this routine can be called a couple of
C      times for the same PS point during the stability checks.
C     But it is really not time consuming and I would rather be safe.
      CALL MP_UPDATE_AS_PARAM()

      DO H=1,NCOMB
        IF ((HELPICKED.EQ.H).OR.((HELPICKED.EQ.-1).AND.((CHECKPHASE.OR.
     $.NOT.HELDOUBLECHECKED).OR.GOODHEL(H)))) THEN
C         Handle the possible requirement of specific polarizations
          IF ((.NOT.CHECKPHASE).AND.HELDOUBLECHECKED.AND.POLARIZATIONS(
     $0,0).EQ.0.AND.(.NOT.ML5_0_IS_HEL_SELECTED(H))) THEN
            CYCLE
          ENDIF
          DO I=1,NEXTERNAL
            NHEL(I)=HELC(I,H)
          ENDDO
          CALL MP_IXXXXX(P(0,1),ZERO,NHEL(1),+1*IC(1),W(1,1,H))
          CALL MP_OXXXXX(P(0,2),ZERO,NHEL(2),-1*IC(2),W(1,2,H))
          CALL MP_OXXXXX(P(0,3),MDL_MT,NHEL(3),+1*IC(3),W(1,3,H))
          CALL MP_IXXXXX(P(0,4),MDL_MT,NHEL(4),-1*IC(4),W(1,4,H))
          CALL MP_FFV1P0_3(W(1,1,H),W(1,2,H),GC_5,ZERO,ZERO,W(1,5,H))
C         Amplitude(s) for born diagram with ID 1
          CALL MP_FFV1_0(W(1,4,H),W(1,3,H),W(1,5,H),GC_5,AMP(1,H))
          CALL MP_FFV1P0_3(W(1,4,H),W(1,3,H),GC_5,ZERO,ZERO,W(1,6,H))
C         Counter-term amplitude(s) for loop diagram number 2
          CALL MP_R2_GG_1_R2_GG_2_0(W(1,5,H),W(1,6,H),R2_GGG_1
     $     ,R2_GGG_2,AMPL(1,1))
C         Counter-term amplitude(s) for loop diagram number 3
          CALL MP_FFV1_0(W(1,4,H),W(1,3,H),W(1,5,H),R2_GQQ,AMPL(1,2))
          CALL MP_FFV1_0(W(1,4,H),W(1,3,H),W(1,5,H),UV_GQQG_1EPS
     $     ,AMPL(2,3))
          CALL MP_FFV1_0(W(1,4,H),W(1,3,H),W(1,5,H),UV_GQQQ_1EPS
     $     ,AMPL(2,4))
          CALL MP_FFV1_0(W(1,4,H),W(1,3,H),W(1,5,H),UV_GQQQ_1EPS
     $     ,AMPL(2,5))
          CALL MP_FFV1_0(W(1,4,H),W(1,3,H),W(1,5,H),UV_GQQQ_1EPS
     $     ,AMPL(2,6))
          CALL MP_FFV1_0(W(1,4,H),W(1,3,H),W(1,5,H),UV_GQQQ_1EPS
     $     ,AMPL(2,7))
          CALL MP_FFV1_0(W(1,4,H),W(1,3,H),W(1,5,H),UV_GQQQ_1EPS
     $     ,AMPL(2,8))
          CALL MP_FFV1_0(W(1,4,H),W(1,3,H),W(1,5,H),UV_GQQQ_1EPS
     $     ,AMPL(2,9))
          CALL MP_FFV1_0(W(1,4,H),W(1,3,H),W(1,5,H),UV_GQQB,AMPL(1,10))
          CALL MP_FFV1_0(W(1,4,H),W(1,3,H),W(1,5,H),UV_GQQT,AMPL(1,11))
C         Counter-term amplitude(s) for loop diagram number 5
          CALL MP_FFV1_0(W(1,1,H),W(1,2,H),W(1,6,H),R2_GQQ,AMPL(1,12))
          CALL MP_FFV1_0(W(1,1,H),W(1,2,H),W(1,6,H),UV_GQQG_1EPS
     $     ,AMPL(2,13))
          CALL MP_FFV1_0(W(1,1,H),W(1,2,H),W(1,6,H),UV_GQQQ_1EPS
     $     ,AMPL(2,14))
          CALL MP_FFV1_0(W(1,1,H),W(1,2,H),W(1,6,H),UV_GQQQ_1EPS
     $     ,AMPL(2,15))
          CALL MP_FFV1_0(W(1,1,H),W(1,2,H),W(1,6,H),UV_GQQQ_1EPS
     $     ,AMPL(2,16))
          CALL MP_FFV1_0(W(1,1,H),W(1,2,H),W(1,6,H),UV_GQQQ_1EPS
     $     ,AMPL(2,17))
          CALL MP_FFV1_0(W(1,1,H),W(1,2,H),W(1,6,H),UV_GQQQ_1EPS
     $     ,AMPL(2,18))
          CALL MP_FFV1_0(W(1,1,H),W(1,2,H),W(1,6,H),UV_GQQQ_1EPS
     $     ,AMPL(2,19))
          CALL MP_FFV1_0(W(1,1,H),W(1,2,H),W(1,6,H),UV_GQQB,AMPL(1,20))
          CALL MP_FFV1_0(W(1,1,H),W(1,2,H),W(1,6,H),UV_GQQT,AMPL(1,21))
C         Counter-term amplitude(s) for loop diagram number 10
          CALL MP_R2_GG_1_0(W(1,5,H),W(1,6,H),R2_GGQ,AMPL(1,22))
          CALL MP_R2_GG_1_0(W(1,5,H),W(1,6,H),R2_GGQ,AMPL(1,23))
          CALL MP_R2_GG_1_0(W(1,5,H),W(1,6,H),R2_GGQ,AMPL(1,24))
          CALL MP_R2_GG_1_0(W(1,5,H),W(1,6,H),R2_GGQ,AMPL(1,25))
C         Counter-term amplitude(s) for loop diagram number 11
          CALL MP_R2_GG_1_R2_GG_3_0(W(1,5,H),W(1,6,H),R2_GGQ,R2_GGB
     $     ,AMPL(1,26))
C         Counter-term amplitude(s) for loop diagram number 12
          CALL MP_R2_GG_1_R2_GG_3_0(W(1,5,H),W(1,6,H),R2_GGQ,R2_GGT
     $     ,AMPL(1,27))
C         Amplitude(s) for UVCT diagram with ID 13
          CALL MP_FFV1_0(W(1,4,H),W(1,3,H),W(1,5,H),GC_5,AMPL(1,28))
          AMPL(1,28)=AMPL(1,28)*(2.0D0*UVWFCT_T_0)
C         Amplitude(s) for UVCT diagram with ID 14
          CALL MP_FFV1_0(W(1,4,H),W(1,3,H),W(1,5,H),GC_5,AMPL(2,29))
          AMPL(2,29)=AMPL(2,29)*(2.0D0*UVWFCT_B_0_1EPS)
C         Copy the qp wfs to the dp ones as they are used to setup the
C          CT calls.
          DO I=1,NWAVEFUNCS
            DO J=1,20
              DPW(J,I,H)=W(J,I,H)
            ENDDO
          ENDDO
C         Same for the counterterms amplitudes
          DO I=1,NCTAMPS
            DO J=1,3
              DPAMPL(J,I)=AMPL(J,I)
              S(I)=.TRUE.
            ENDDO
          ENDDO
          DO I=1,NBORNAMPS
            DPAMP(I,H)=AMP(I,H)
          ENDDO
        ENDIF
      ENDDO

      END

