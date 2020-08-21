      SUBROUTINE MP_COMPUTE_LOOP_COEFS(PS,ANSDP)
C     
C     Generated by MadGraph5_aMC@NLO v. %(version)s, %(date)s
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     Returns amplitude squared summed/avg over colors
C     and helicities for the point in phase space P(0:3,NEXTERNAL)
C     and external lines W(0:6,NEXTERNAL)
C     
C     Process: u d~ > w+ [ all = QCD QED ] QCD^2=2 QED^2=2
C     Process: c s~ > w+ [ all = QCD QED ] QCD^2=2 QED^2=2
C     
C     Modules
C     
      USE POLYNOMIAL_CONSTANTS
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      CHARACTER*64 PARAMFILENAME
      PARAMETER ( PARAMFILENAME='MadLoopParams.dat')
      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=1)
      INTEGER    NLOOPS, NLOOPGROUPS, NCTAMPS
      PARAMETER (NLOOPS=1, NLOOPGROUPS=1, NCTAMPS=1)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=2)
      INTEGER    NCOLORROWS
      PARAMETER (NCOLORROWS=NLOOPAMPS)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=3)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=3,NLOOPWAVEFUNCS=3)
      INTEGER    NCOMB
      PARAMETER (NCOMB=12)
      REAL*16    ZERO
      PARAMETER (ZERO=0E0_16)
      COMPLEX*32 IMAG1
      PARAMETER (IMAG1=(0E0_16,1E0_16))
      COMPLEX*32 DP_IMAG1
      PARAMETER (DP_IMAG1=(0D0,1D0))
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=2, NSQUAREDSO=1, NAMPSO=2)


C     
C     ARGUMENTS 
C     
      REAL*16 PS(0:3,NEXTERNAL)
      REAL*8 ANSDP(3,0:NSQUAREDSO)
C     
C     LOCAL VARIABLES 
C     
      LOGICAL DPW_COPIED
      LOGICAL COMPUTE_INTEGRAND_IN_QP
      INTEGER I,J,K,H,HEL_MULT,ITEMP
      REAL*16 TEMP2
      REAL*8 DP_TEMP2
      COMPLEX*32 CTEMP
      COMPLEX*16 DP_CTEMP

      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      REAL*16 MP_P(0:3,NEXTERNAL)
      REAL*8 P(0:3,NEXTERNAL)

      DATA IC/NEXTERNAL*1/
      REAL*16 ANS(3,0:NSQUAREDSO)
      COMPLEX*32 COEFS(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 CFTOT
      COMPLEX*16 DP_CFTOT
C     
C     FUNCTIONS
C     
      LOGICAL IS_HEL_SELECTED
      INTEGER ML5SOINDEX_FOR_BORN_AMP
      INTEGER ML5SOINDEX_FOR_LOOP_AMP
      INTEGER ML5SQSOINDEX
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'mp_coupl_same_name.inc'

      INCLUDE 'MadLoopParams.inc'

      LOGICAL CHECKPHASE, HELDOUBLECHECKED
      COMMON/INIT/CHECKPHASE, HELDOUBLECHECKED

      INTEGER HELOFFSET
      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/FILTERS/GOODAMP,GOODHEL,HELOFFSET

      INTEGER HELPICKED
      COMMON/HELCHOICE/HELPICKED

      INTEGER USERHEL
      COMMON/USERCHOICE/USERHEL

      INTEGER SQSO_TARGET
      COMMON/SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

      COMPLEX*32 AMP(NBORNAMPS)
      COMMON/MP_AMPS/AMP
      COMPLEX*16 DP_AMP(NBORNAMPS)
      COMMON/AMPS/DP_AMP
      COMPLEX*32 W(20,NWAVEFUNCS)
      COMMON/MP_W/W

      COMPLEX*16 DPW(20,NWAVEFUNCS)
      COMMON/W/DPW

      COMPLEX*32 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE,
     $ -1:NLOOPWAVEFUNCS)
      COMPLEX*32 PL(0:3,-1:NLOOPWAVEFUNCS)
      COMMON/MP_WL/WL,PL

      COMPLEX*16 DP_WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE,
     $ -1:NLOOPWAVEFUNCS)
      COMPLEX*16 DP_PL(0:3,-1:NLOOPWAVEFUNCS)
      COMMON/WL/DP_WL,DP_PL

      COMPLEX*32 LOOPCOEFS(0:LOOPMAXCOEFS-1,NSQUAREDSO,NLOOPGROUPS)
      COMMON/MP_LCOEFS/LOOPCOEFS

      COMPLEX*16 DP_LOOPCOEFS(0:LOOPMAXCOEFS-1,NSQUAREDSO,NLOOPGROUPS)
      COMMON/LCOEFS/DP_LOOPCOEFS


      COMPLEX*32 AMPL(3,NCTAMPS)
      COMMON/MP_AMPL/AMPL

      COMPLEX*16 DP_AMPL(3,NCTAMPS)
      COMMON/AMPL/DP_AMPL


      INTEGER CF_D(NCOLORROWS,NBORNAMPS)
      INTEGER CF_N(NCOLORROWS,NBORNAMPS)
      COMMON/CF/CF_D,CF_N

      INTEGER HELC(NEXTERNAL,NCOMB)
      COMMON/HELCONFIGS/HELC

      LOGICAL MP_DONE_ONCE
      COMMON/MP_DONE_ONCE/MP_DONE_ONCE

      INTEGER LIBINDEX
      COMMON/I_LIB/LIBINDEX

C     This array specify potential special requirements on the
C      helicities to
C     consider. POLARIZATIONS(0,0) is -1 if there is not such
C      requirement.
      INTEGER POLARIZATIONS(0:NEXTERNAL,0:5)
      COMMON/BEAM_POL/POLARIZATIONS

C     ----------
C     BEGIN CODE
C     ----------

C     Decide whether to really compute the integrand in quadruple
C      precision or to fake it and copy the double precision
C      computation in the quadruple precision variables.
      COMPUTE_INTEGRAND_IN_QP = ((MLREDUCTIONLIB(LIBINDEX)
     $ .EQ.6.AND.USEQPINTEGRANDFORNINJA) .OR. (MLREDUCTIONLIB(LIBINDEX)
     $ .EQ.1.AND.USEQPINTEGRANDFORCUTTOOLS))

C     To be on the safe side, we always update the MP params here.
C     It can be redundant as this routine can be called a couple of
C      times for the same PS point during the stability checks.
C     But it is really not time consuming and I would rather be safe.
      CALL MP_UPDATE_AS_PARAM()

      MP_DONE_ONCE = .TRUE.

C     AS A SAFETY MEASURE WE FIRST COPY HERE THE PS POINT
      DO I=1,NEXTERNAL
        DO J=0,3
          MP_P(J,I)=PS(J,I)
          P(J,I) = REAL(PS(J,I),KIND=8)
        ENDDO
      ENDDO

      DO I=0,3
        PL(I,-1)=CMPLX(ZERO,ZERO,KIND=16)
        PL(I,0)=CMPLX(ZERO,ZERO,KIND=16)
        IF (.NOT.COMPUTE_INTEGRAND_IN_QP) THEN
          DP_PL(I,-1)=DCMPLX(0.0D0,0.0D0)
          DP_PL(I,0)=DCMPLX(0.0D0,0.0D0)
        ENDIF
      ENDDO

      DO I=1,MAXLWFSIZE
        DO J=0,LOOPMAXCOEFS-1
          DO K=1,MAXLWFSIZE
            WL(I,J,K,-1)=(ZERO,ZERO)
            DP_WL(I,J,K,-1)=(0.0D0,0.0D0)
            IF (I.EQ.K.AND.J.EQ.0) THEN
              WL(I,J,K,0)=(1.0E0_16,ZERO)
            ELSE
              WL(I,J,K,0)=(ZERO,ZERO)
            ENDIF
            IF (.NOT.COMPUTE_INTEGRAND_IN_QP) THEN
              IF (I.EQ.K.AND.J.EQ.0) THEN
                DP_WL(I,J,K,0)=(1.0D0,0.0D0)
              ELSE
                DP_WL(I,J,K,0)=(0.0D0,0.0D0)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

C     This is the chare conjugate version of the unit 4-currents in
C      the canonical cartesian basis.
C     This, for now, is only defined for 4-fermionic currents.
      WL(1,0,2,-1) = (-1.0E0_16,ZERO)
      WL(2,0,1,-1) = (1.0E0_16,ZERO)
      WL(3,0,4,-1) = (1.0E0_16,ZERO)
      WL(4,0,3,-1) = (-1.0E0_16,ZERO)
      DP_WL(1,0,2,-1) = DCMPLX(-1.0D0,0.0D0)
      DP_WL(2,0,1,-1) = DCMPLX(1.0D0,0.0D0)
      DP_WL(3,0,4,-1) = DCMPLX(1.0D0,0.0D0)
      DP_WL(4,0,3,-1) = DCMPLX(-1.0D0,0.0D0)


      DO K=1, 3
        DO I=1,NCTAMPS
          AMPL(K,I)=(ZERO,ZERO)
          IF (.NOT.COMPUTE_INTEGRAND_IN_QP) THEN
            DP_AMPL(K,I)=(0.0D0,0.0D0)
          ENDIF
        ENDDO
      ENDDO


      DO I=1, NBORNAMPS
        DP_AMP(I) = (0.0D0,0.0D0)
        AMP(I) = (ZERO, ZERO)
      ENDDO

      DO I=1,NLOOPGROUPS
        DO J=0,LOOPMAXCOEFS-1
          DO K=1,NSQUAREDSO
            LOOPCOEFS(J,K,I)=(ZERO,ZERO)
            IF (.NOT.COMPUTE_INTEGRAND_IN_QP) THEN
              DP_LOOPCOEFS(J,K,I)=(0.0D0,0.0D0)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DO K=1,3
        DO J=0,NSQUAREDSO
          ANSDP(K,J)=0.0D0
          ANS(K,J)=ZERO
        ENDDO
      ENDDO

      DPW_COPIED = .FALSE.
      DO H=1,NCOMB
        IF ((HELPICKED.EQ.H).OR.((HELPICKED.EQ.-1)
     $   .AND.(CHECKPHASE.OR.(.NOT.HELDOUBLECHECKED).OR.(GOODHEL(H)
     $   .GT.-HELOFFSET.AND.GOODHEL(H).NE.0)))) THEN

C         Handle the possible requirement of specific polarizations
          IF ((.NOT.CHECKPHASE)
     $     .AND.HELDOUBLECHECKED.AND.POLARIZATIONS(0,0)
     $     .EQ.0.AND.(.NOT.IS_HEL_SELECTED(H))) THEN
            CYCLE
          ENDIF

          DO I=1,NEXTERNAL
            NHEL(I)=HELC(I,H)
          ENDDO

          IF (COMPUTE_INTEGRAND_IN_QP) THEN
            MP_UVCT_REQ_SO_DONE=.FALSE.
            MP_CT_REQ_SO_DONE=.FALSE.
            MP_LOOP_REQ_SO_DONE=.FALSE.
          ELSE
            UVCT_REQ_SO_DONE=.FALSE.
            CT_REQ_SO_DONE=.FALSE.
            LOOP_REQ_SO_DONE=.FALSE.
          ENDIF

          IF (.NOT.CHECKPHASE.AND.HELDOUBLECHECKED.AND.HELPICKED.EQ.-1)
     $      THEN
            HEL_MULT=GOODHEL(H)
          ELSE
            HEL_MULT=1
          ENDIF


          IF (COMPUTE_INTEGRAND_IN_QP) THEN
            CALL MP_HELAS_CALLS_AMPB_1(MP_P,NHEL,H,IC)
            CONTINUE
          ELSE
            CALL HELAS_CALLS_AMPB_1(P,NHEL,H,IC)
            CONTINUE
          ENDIF

 2000     CONTINUE
          MP_CT_REQ_SO_DONE=.TRUE.

          IF (COMPUTE_INTEGRAND_IN_QP) THEN

            CONTINUE
          ELSE

            CONTINUE
          ENDIF

          IF (.NOT.COMPUTE_INTEGRAND_IN_QP) THEN
C           Copy back to the quantities computed in DP in the QP
C            containers (but only those needed)
            DO I=1,NBORNAMPS
              AMP(I)=CMPLX(DP_AMP(I),KIND=16)
            ENDDO
            DO I=1,NCTAMPS
              DO K=1,3
                AMPL(K,I)=CMPLX(DP_AMPL(K,I),KIND=16)
              ENDDO
            ENDDO
            DO I=1,NWAVEFUNCS
              DO J=1,MAXLWFSIZE+4
                W(J,I)=CMPLX(DPW(J,I),KIND=16)
              ENDDO
            ENDDO
          ENDIF

 3000     CONTINUE
          MP_UVCT_REQ_SO_DONE=.TRUE.

          IF (COMPUTE_INTEGRAND_IN_QP) THEN

            DO J=1,NBORNAMPS
              CTEMP = HEL_MULT*2.0E0_16*CONJG(AMP(J))
              DO I=1,NCTAMPS
                CFTOT=CMPLX(CF_N(I,J)/REAL(ABS(CF_D(I,J)),KIND=16)
     $           ,0.0E0_16,KIND=16)
                IF(CF_D(I,J).LT.0) CFTOT=CFTOT*IMAG1
                ITEMP = ML5SQSOINDEX(ML5SOINDEX_FOR_LOOP_AMP(I)
     $           ,ML5SOINDEX_FOR_BORN_AMP(J))
                IF (.NOT.FILTER_SO.OR.SQSO_TARGET.EQ.ITEMP) THEN
                  DO K=1,3
                    TEMP2 = REAL(CFTOT*AMPL(K,I)*CTEMP,KIND=16)
                    ANS(K,ITEMP)=ANS(K,ITEMP)+TEMP2
                    ANS(K,0)=ANS(K,0)+TEMP2
                  ENDDO
                ENDIF
              ENDDO
            ENDDO

          ELSE

            DO J=1,NBORNAMPS
              DP_CTEMP = HEL_MULT*2.0D0*DCONJG(DP_AMP(J))
              DO I=1,NCTAMPS
                DP_CFTOT=CMPLX(CF_N(I,J)/REAL(ABS(CF_D(I,J)),KIND=8)
     $           ,0.0D0,KIND=8)
                IF(CF_D(I,J).LT.0) DP_CFTOT=DP_CFTOT*DP_IMAG1
                ITEMP = ML5SQSOINDEX(ML5SOINDEX_FOR_LOOP_AMP(I)
     $           ,ML5SOINDEX_FOR_BORN_AMP(J))
                IF (.NOT.FILTER_SO.OR.SQSO_TARGET.EQ.ITEMP) THEN
                  DO K=1,3
                    DP_TEMP2 = REAL(DP_CFTOT*DP_AMPL(K,I)*DP_CTEMP
     $               ,KIND=8)
                    ANSDP(K,ITEMP)=ANSDP(K,ITEMP)+DP_TEMP2
                    ANSDP(K,0)=ANSDP(K,0)+DP_TEMP2
                  ENDDO
                ENDIF
              ENDDO
            ENDDO

          ENDIF


          IF (COMPUTE_INTEGRAND_IN_QP) THEN

            CALL MP_COEF_CONSTRUCTION_1(MP_P,NHEL,H,IC)

          ELSE

            CALL COEF_CONSTRUCTION_1(P,NHEL,H,IC)

C           Copy back to the coefficients computed in DP in the QP
C            containers
            DO I=0,LOOPMAXCOEFS-1
              DO K=1,NLOOPGROUPS
                DO J=1,NSQUAREDSO
                  LOOPCOEFS(I,J,K)=CMPLX(DP_LOOPCOEFS(I,J,K),KIND=16)
                ENDDO
              ENDDO
            ENDDO
          ENDIF

 4000     CONTINUE
          MP_LOOP_REQ_SO_DONE=.TRUE.


C         Copy the qp wfs to the dp ones as they are used to setup the
C          CT calls.
C         This needs to be done once since only the momenta of these
C          WF matters.
          IF(.NOT.DPW_COPIED.AND.COMPUTE_INTEGRAND_IN_QP) THEN
            DO I=1,NWAVEFUNCS
              DO J=1,MAXLWFSIZE+4
                DPW(J,I)=CMPLX(W(J,I),KIND=8)
              ENDDO
            ENDDO
            DPW_COPIED=.TRUE.
          ENDIF





        ENDIF
      ENDDO


C     If we were not computing the integrand in QP, then we were
C      already updating ANSDP all along, so that fetching it here from
C      the QP ANS(:,:) should not be done.
      IF (COMPUTE_INTEGRAND_IN_QP) THEN
        DO I=1,3
          DO J=0,NSQUAREDSO
            ANSDP(I,J)=REAL(ANS(I,J),KIND=8)
          ENDDO
        ENDDO
      ENDIF

C     Grouping of loop diagrams now done directly when creating the
C      LOOPCOEFS.
C     If some kind of coefficient merging was done above, do not
C      forget to copy back the LOOPCOEFS merged into DP_LOOPCOEFS if
C      COMPUTE_INTEGRAND_IN_QP is False.

      END

