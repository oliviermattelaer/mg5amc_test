      SUBROUTINE ML5_0_MP_BORN_AMPS_AND_WFS(P)
C     
C     Generated by MadGraph5_aMC@NLO v. %(version)s, %(date)s
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     Computes all the AMP and WFS in quadruple precision for the 
C     phase space point P(0:3,NEXTERNAL)
C     
C     Process: d u~ > m- vm~ g QCD<=1 QED<=2 [ virt = QCD ]
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=2)
      INTEGER    NLOOPAMPS, NCTAMPS
      PARAMETER (NLOOPAMPS=39, NCTAMPS=28)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=5)
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=10)
      INTEGER    NCOMB
      PARAMETER (NCOMB=32)
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
        IF ((HELPICKED.EQ.H).OR.((HELPICKED.EQ.-1)
     $   .AND.((CHECKPHASE.OR..NOT.HELDOUBLECHECKED).OR.GOODHEL(H))))
     $    THEN
C         Handle the possible requirement of specific polarizations
          IF ((.NOT.CHECKPHASE)
     $     .AND.HELDOUBLECHECKED.AND.POLARIZATIONS(0,0)
     $     .EQ.0.AND.(.NOT.ML5_0_IS_HEL_SELECTED(H))) THEN
            CYCLE
          ENDIF
          DO I=1,NEXTERNAL
            NHEL(I)=HELC(I,H)
          ENDDO
          CALL MP_IXXXXX(P(0,1),ZERO,NHEL(1),+1*IC(1),W(1,1,H))
          CALL MP_OXXXXX(P(0,2),ZERO,NHEL(2),-1*IC(2),W(1,2,H))
          CALL MP_OXXXXX(P(0,3),ZERO,NHEL(3),+1*IC(3),W(1,3,H))
          CALL MP_IXXXXX(P(0,4),ZERO,NHEL(4),-1*IC(4),W(1,4,H))
          CALL MP_VXXXXX(P(0,5),ZERO,NHEL(5),+1*IC(5),W(1,5,H))
          CALL MP_FFV1_2(W(1,1,H),W(1,5,H),GC_5,ZERO,ZERO,W(1,6,H))
          CALL MP_FFV2_3(W(1,4,H),W(1,3,H),GC_47,MDL_MW,MDL_WW,W(1,7,H)
     $     )
C         Amplitude(s) for born diagram with ID 1
          CALL MP_FFV2_0(W(1,6,H),W(1,2,H),W(1,7,H),GC_47,AMP(1,H))
          CALL MP_FFV1_1(W(1,2,H),W(1,5,H),GC_5,ZERO,ZERO,W(1,8,H))
C         Amplitude(s) for born diagram with ID 2
          CALL MP_FFV2_0(W(1,1,H),W(1,8,H),W(1,7,H),GC_47,AMP(2,H))
          CALL MP_FFV2_1(W(1,2,H),W(1,7,H),GC_47,ZERO,ZERO,W(1,9,H))
C         Counter-term amplitude(s) for loop diagram number 3
          CALL MP_R2_QQ_1_0(W(1,6,H),W(1,9,H),R2_QQQ,AMPL(1,1))
C         Counter-term amplitude(s) for loop diagram number 4
          CALL MP_FFV2_0(W(1,6,H),W(1,2,H),W(1,7,H),R2_SXCW,AMPL(1,2))
C         Counter-term amplitude(s) for loop diagram number 5
          CALL MP_FFV2_0(W(1,1,H),W(1,8,H),W(1,7,H),R2_SXCW,AMPL(1,3))
C         Counter-term amplitude(s) for loop diagram number 7
          CALL MP_FFV1_0(W(1,1,H),W(1,9,H),W(1,5,H),UV_GQQQ_1EPS
     $     ,AMPL(2,4))
          CALL MP_FFV1_0(W(1,1,H),W(1,9,H),W(1,5,H),UV_GQQQ_1EPS
     $     ,AMPL(2,5))
          CALL MP_FFV1_0(W(1,1,H),W(1,9,H),W(1,5,H),UV_GQQQ_1EPS
     $     ,AMPL(2,6))
          CALL MP_FFV1_0(W(1,1,H),W(1,9,H),W(1,5,H),UV_GQQQ_1EPS
     $     ,AMPL(2,7))
          CALL MP_FFV1_0(W(1,1,H),W(1,9,H),W(1,5,H),UV_GQQB,AMPL(1,8))
          CALL MP_FFV1_0(W(1,1,H),W(1,9,H),W(1,5,H),UV_GQQQ_1EPS
     $     ,AMPL(2,9))
          CALL MP_FFV1_0(W(1,1,H),W(1,9,H),W(1,5,H),UV_GQQT,AMPL(1,10))
          CALL MP_FFV1_0(W(1,1,H),W(1,9,H),W(1,5,H),UV_GQQQ_1EPS
     $     ,AMPL(2,11))
          CALL MP_FFV1_0(W(1,1,H),W(1,9,H),W(1,5,H),UV_GQQG_1EPS
     $     ,AMPL(2,12))
          CALL MP_FFV1_0(W(1,1,H),W(1,9,H),W(1,5,H),R2_GQQ,AMPL(1,13))
          CALL MP_FFV2_2(W(1,1,H),W(1,7,H),GC_47,ZERO,ZERO,W(1,10,H))
C         Counter-term amplitude(s) for loop diagram number 11
          CALL MP_R2_QQ_1_0(W(1,10,H),W(1,8,H),R2_QQQ,AMPL(1,14))
C         Counter-term amplitude(s) for loop diagram number 12
          CALL MP_FFV1_0(W(1,10,H),W(1,2,H),W(1,5,H),UV_GQQQ_1EPS
     $     ,AMPL(2,15))
          CALL MP_FFV1_0(W(1,10,H),W(1,2,H),W(1,5,H),UV_GQQQ_1EPS
     $     ,AMPL(2,16))
          CALL MP_FFV1_0(W(1,10,H),W(1,2,H),W(1,5,H),UV_GQQQ_1EPS
     $     ,AMPL(2,17))
          CALL MP_FFV1_0(W(1,10,H),W(1,2,H),W(1,5,H),UV_GQQQ_1EPS
     $     ,AMPL(2,18))
          CALL MP_FFV1_0(W(1,10,H),W(1,2,H),W(1,5,H),UV_GQQB,AMPL(1,19)
     $     )
          CALL MP_FFV1_0(W(1,10,H),W(1,2,H),W(1,5,H),UV_GQQQ_1EPS
     $     ,AMPL(2,20))
          CALL MP_FFV1_0(W(1,10,H),W(1,2,H),W(1,5,H),UV_GQQT,AMPL(1,21)
     $     )
          CALL MP_FFV1_0(W(1,10,H),W(1,2,H),W(1,5,H),UV_GQQQ_1EPS
     $     ,AMPL(2,22))
          CALL MP_FFV1_0(W(1,10,H),W(1,2,H),W(1,5,H),UV_GQQG_1EPS
     $     ,AMPL(2,23))
          CALL MP_FFV1_0(W(1,10,H),W(1,2,H),W(1,5,H),R2_GQQ,AMPL(1,24))
C         Amplitude(s) for UVCT diagram with ID 14
          CALL MP_FFV2_0(W(1,6,H),W(1,2,H),W(1,7,H),GC_47,AMPL(1,25))
          AMPL(1,25)=AMPL(1,25)*(1.0D0*UVWFCT_G_2+1.0D0*UVWFCT_G_1)
C         Amplitude(s) for UVCT diagram with ID 15
          CALL MP_FFV2_0(W(1,6,H),W(1,2,H),W(1,7,H),GC_47,AMPL(2,26))
          AMPL(2,26)=AMPL(2,26)*(2.0D0*UVWFCT_G_2_1EPS)
C         Amplitude(s) for UVCT diagram with ID 16
          CALL MP_FFV2_0(W(1,1,H),W(1,8,H),W(1,7,H),GC_47,AMPL(1,27))
          AMPL(1,27)=AMPL(1,27)*(1.0D0*UVWFCT_G_2+1.0D0*UVWFCT_G_1)
C         Amplitude(s) for UVCT diagram with ID 17
          CALL MP_FFV2_0(W(1,1,H),W(1,8,H),W(1,7,H),GC_47,AMPL(2,28))
          AMPL(2,28)=AMPL(2,28)*(2.0D0*UVWFCT_G_2_1EPS)
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

