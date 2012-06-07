C     -------------------------
      SUBROUTINE SMATRIX%(proc_id)s(P,ANS)
C     -------------------------
C  
%(info_lines)s
C 
C Color-ordered MadGraph for Madevent Version
C 
C Returns amplitude squared summed/avg over colors
c and helicities
c for the point in phase space P(0:3,NEXTERNAL)
C  
%(process_lines)s
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
    Include 'genps.inc'
    Include 'maxconfigs.inc'
    Include 'nexternal.inc'
    Include 'maxamps.inc'
      INTEGER                 NCOMB         
      PARAMETER (             NCOMB=%(ncomb)d)
    INTEGER    NGRAPHS
    PARAMETER (NGRAPHS=%(ngraphs)d) 
    INTEGER    NDIAGS
    PARAMETER (NDIAGS=%(ndiags)d) 
    INTEGER    THEL
    PARAMETER (THEL=NCOMB)
    INTEGER    NPERMS
    PARAMETER (NPERMS=%(nperms)d) 
    INTEGER    NFLOWS
    PARAMETER (NFLOWS=%(nflows)d) 
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL),ANS
C  
C LOCAL VARIABLES 
C  
    INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
    REAL*8 T, MATRIX%(proc_id)s
    REAL*8 R,SUMHEL,TS(NCOMB)
    INTEGER I,IDEN
    INTEGER IPROC,II
    LOGICAL GOODHEL(NCOMB)
    REAL*8 HWGT,XTOT,XTRY,XREJ,XR,YFRAC(0:NCOMB),FRAC
    INTEGER IDUM,IPERM,NGOOD,IGOOD(NCOMB),JHEL,J,JJ
    REAL*8 TOTFACT,FACTNUL
    INTEGER JHELBK,NUPPER
    LOGICAL FIRSTFLOW
    REAL*8   RVEC
    SAVE FRAC,FACTNUL
C  
C GLOBAL VARIABLES
C  
    DOUBLE PRECISION AMP2(MAXAMPS),JAMP2(0:MAXFLOW)
    COMMON/TO_AMPS/  AMP2,         JAMP2
    
    CHARACTER*101        HEL_BUFF
    COMMON/TO_HELICITY/  HEL_BUFF
    
    REAL*8 POL(2)
    COMMON/TO_POLARIZATION/ POL

    INTEGER          ISUM_HEL
    LOGICAL                    MULTI_CHANNEL
    COMMON/TO_MATRIX/ISUM_HEL, MULTI_CHANNEL

    INTEGER PERM(NEXTERNAL)    
    COMMON/TO_COPERM/PERM

%(define_iconfigs_lines)s
    DATA NTRY,IDUM /0,-1/
    DATA XTRY, XREJ, NGOOD /0,0,0/
    SAVE YFRAC, IGOOD, JHEL
    DATA GOODHEL/THEL*.FALSE./
%(helicity_lines)s
%(den_factor_line)s
C ----------
C BEGIN CODE
C ----------
    NTRY=NTRY+1

    IF (multi_channel) THEN
        DO I=1,NDIAGS
            AMP2(I)=0D0
        ENDDO
        JAMP2(0)=NFLOWS
        DO I=1,INT(JAMP2(0))
            JAMP2(I)=0D0
        ENDDO
    ENDIF
    ANS = 0D0
    WRITE(HEL_BUFF,'(20I5)') (0,I=1,NEXTERNAL)
    DO I=1,NCOMB
       TS(I)=0d0
    ENDDO
    IF(NTRY.EQ.1)THEN
c     Calculate summed permutation weight
      TOTFACT=NPERMS
      FRAC=MIN(NUMPERMS,NPERMS)*1D0/TOTFACT
      FACTNUL=1d0
C     Correction factor for events where no permutation is selected
      DO IPERM=1,NPERMS
         FACTNUL=FACTNUL*(1-FRAC)
      ENDDO
      FACTNUL=1-FACTNUL
    ENDIF
c   Decide between helicity all-sum case (ISUM_HEL=0) and partial sum case
      JHELBK=JHEL
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LE. MAXTRIES) THEN
         NUPPER=NCOMB
         HWGT=1.
      ELSE
         NUPPER=ISUM_HEL
         HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
      ENDIF
 10 FIRSTFLOW=.TRUE.
c     Start loop over permutations
      DO IPERM=1,NPERMS
        JHEL=JHELBK
C     Decide whether to use this permutation
        CALL RANMAR(RVEC)
        IF(RVEC.GT.FRAC) CYCLE
c     Get permutation
        CALL GETPERM%(proc_id)s(IPERM,PERM)
c     Start loop over helicities
        DO J=1,NUPPER
          IF (ISUM_HEL .EQ. 0 .OR. NTRY .LE. MAXTRIES) THEN
            I=J
            IF(.NOT.GOODHEL(I) .AND. NTRY .GT. MAXTRIES) CYCLE
          ELSE
            JHEL=JHEL+1
            IF (JHEL .GT. NGOOD) JHEL=1
            I = IGOOD(JHEL)
          ENDIF
          T=MATRIX%(proc_id)s(P,NHEL(1,I),PERM,FIRSTFLOW)
          DO JJ=1,NINCOMING
            IF(POL(JJ).NE.1D0.AND.NHEL(JJ,I).EQ.INT(SIGN(1D0,POL(JJ)))) THEN
              T=T*ABS(POL(JJ))
            ELSE IF(POL(JJ).NE.1D0)THEN
              T=T*(2D0-ABS(POL(JJ)))
            ENDIF
          ENDDO
          ANS=ANS+T*HWGT/FRAC*FACTNUL
          TS(I)=T*HWGT/FRAC*FACTNUL
        ENDDO
        FIRSTFLOW=.FALSE.
c        print *,'ANS: ',ANS
        IF(NTRY.LE.MAXTRIES)THEN
          DO I=1,NCOMB
            IF (.NOT.GOODHEL(I) .AND. (TS(I).GT.ANS*LIMHEL/NCOMB)) THEN
              GOODHEL(I)=.TRUE.
              NGOOD = NGOOD +1
              IGOOD(NGOOD) = I
              PRINT *,'Added good helicity ',I,TS(I)/ANS,' in event ',NTRY(IMIRROR)
            ENDIF
          ENDDO
        ENDIF
        IF(NTRY.EQ.MAXTRIES)THEN
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
        ENDIF
        IF (ISUM_HEL .EQ. 1) THEN
          WRITE(HEL_BUFF,'(20i5)')(NHEL(II,I),II=1,NEXTERNAL)
        ENDIF
      ENDDO
C     If no permutation was chosen, repeat
      IF(FIRSTFLOW) GOTO 10
      IF (ISUM_HEL .NE. 1) THEN
        CALL RANMAR(RVEC)
        R=RVEC*ANS
        SUMHEL=0d0
        DO I=1,NCOMB
          SUMHEL=SUMHEL+TS(I)
          IF(R.LT.SUMHEL)THEN
            WRITE(HEL_BUFF,'(20i5)')(NHEL(II,I),II=1,NEXTERNAL)
            GOTO 20
          ENDIF
        ENDDO
 20     CONTINUE
      ENDIF
      IF (MULTI_CHANNEL) THEN
        XTOT=0D0
        DO I=1,NDIAGS
          XTOT=XTOT+AMP2(I)
        ENDDO
c        print *,'XTOT,AMP2: ',XTOT,AMP2(SUBDIAG(1))
        IF (XTOT.NE.0D0) THEN
%(set_amp2_line)s
        ELSE
          ANS=0D0
        ENDIF
      ENDIF
      ANS=ANS/DBLE(IDEN)
      END

C     --------------------------------------------------------
      REAL*8 FUNCTION MATRIX%(proc_id)s(P,NHEL,PERM,FIRSTFLOW)
C     --------------------------------------------------------
C
C Returns amplitude squared summed/avg over colors
c for the point with external momenta P(0:3,NEXTERNAL)
C  
%(process_lines)s
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
    INTEGER    NGRAPHS
    PARAMETER (NGRAPHS=%(ngraphs)d) 
    include 'genps.inc'
    include 'nexternal.inc'
    include 'maxamps.inc'
    include 'coupl.inc'
    INTEGER    NWAVEFUNCS
    PARAMETER (NWAVEFUNCS=%(nwavefuncs)d)
    REAL*8     ZERO
    PARAMETER (ZERO=0D0)
    COMPLEX*16 IMAG1
    PARAMETER (IMAG1=(0D0,1D0))
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL),PERM(NEXTERNAL)
      LOGICAL FIRSTFLOW
C  
C GLOBAL VARIABLES
C  
    DOUBLE PRECISION AMP2(MAXAMPS), JAMP2(0:MAXFLOW)
    COMMON/TO_AMPS/  AMP2,       JAMP2
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J,IC(NEXTERNAL),IP(NEXTERNAL)
      COMPLEX*16 ZTEMP
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(%(wavefunctionsize)d,NWAVEFUNCS)
%(ic_data_line)s
%(ip_data_line)s
C
C EXTERNAL FUNCTIONS
C
      COMPLEX*16 ONEPERM%(proc_id)s
      EXTERNAL ONEPERM%(proc_id)s

c     Skip calculating amp2s for PS if secondary flows
      IF(.NOT.FIRSTFLOW) GOTO 10

%(helas_calls)s

      DO I=1,NGRAPHS
        AMP2(I)=AMP2(I)+AMP(I)*DCONJG(AMP(I))
      ENDDO

 10   ZTEMP = (0.D0,0.D0)
      ZTEMP=ZTEMP+ONEPERM%(proc_id)s(P,NHEL,PERM)
      MATRIX%(proc_id)s=REAL(ZTEMP)

      RETURN
      END
      
C     --------------------------------------
      COMPLEX*16 FUNCTION ONEPERM%(proc_id)s(P,NHEL,PM)
C     --------------------------------------
C
C Returns amplitude squared summed/avg over colors
c for the point with external lines W(0:6,NEXTERNAL)
C  
%(process_lines)s
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include 'nexternal.inc'
      include 'maxamps.inc'
      INTEGER    NJAMPS
      PARAMETER (NJAMPS=%(njamps)d) 
      INTEGER    NPERMS
      PARAMETER (NPERMS=%(nflowperms)d)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL)
      INTEGER PM(NEXTERNAL)
C  
C GLOBAL VARIABLES
C  
    DOUBLE PRECISION AMP2(MAXAMPS), JAMP2(0:MAXFLOW)
    COMMON/TO_AMPS/  AMP2,       JAMP2
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      COMPLEX*16 JAMP(NJAMPS)
      INTEGER PERMS(NEXTERNAL,NPERMS),IFERM(NPERMS),PERM(NEXTERNAL)
%(flow_perms_data_lines)s
%(flow_iferm_data_line)s
C
C EXTERNAL FUNCTIONS
C
%(flow_functions_lines)s
C  
C GLOBAL VARIABLES
C  
      include 'coupl.inc'
C ----------
C BEGIN CODE
C ----------
%(flow_call_lines)s
      ZTEMP = (0.D0,0.D0)
%(color_sum_lines)s
      ONEPERM%(proc_id)s=ZTEMP

%(jamp2_lines)s

      RETURN
      END

C     ------------------------------------
      SUBROUTINE GETPERM%(proc_id)s(IPERM,PERM)
C     ------------------------------------
C
C Gives permutation number IPERM. 
C Return value is the fermion factor due to PERM
C  
%(process_lines)s
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include 'nexternal.inc'
C  
C ARGUMENTS 
C  
      INTEGER IPERM,PERM(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J,IFLAG
      LOGICAL OK
      INTEGER COMP(NEXTERNAL)
%(comp_data_line)s
C ----------
C BEGIN CODE
C ----------
      DO I=1,NEXTERNAL
        PERM(I)=I
      ENDDO
      I=1
      DO WHILE(I.LT.IPERM)
         CALL IPNEXT(PERM,NEXTERNAL,IFLAG)
	 OK=.TRUE.
	 DO J=1,NEXTERNAL
           IF(COMP(PERM(J)).ne.COMP(J))THEN
              OK=.FALSE.
              EXIT
           ENDIF
	 ENDDO
         IF(OK) I=I+1
      ENDDO
      END

