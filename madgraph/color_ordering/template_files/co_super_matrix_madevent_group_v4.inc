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
    PARAMETER (THEL=2*NCOMB)
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
    INTEGER NHEL(NEXTERNAL,NCOMB),NTRY(2)
    INTEGER ISHEL(2)
    REAL*8 T, MATRIX%(proc_id)s
    REAL*8 R,SUMHEL,TS(NCOMB)
    INTEGER I,IDEN
    INTEGER IPROC,II
    LOGICAL GOODHEL(NCOMB,2)
    REAL*8 HWGT,XTOT,XTRY,XREJ,XR,YFRAC(0:NCOMB),FRAC
    INTEGER IDUM,IPERM,NGOOD(2),IGOOD(NCOMB,2)
    INTEGER JHEL(2),J,JJ
    REAL*8 TOTFACT,FACTNUL
    INTEGER JHELBK,NUPPER
    LOGICAL FIRSTFLOW
    REAL*8 COMBANS,PERMANS(NPERMS)
    INTEGER IPERMFACT(NPERMS),LPERM(NEXTERNAL),TPERM(NEXTERNAL)
    INTEGER PERMFACT%(proc_id)s
    EXTERNAL PERMFACT%(proc_id)s
    INTEGER PERMS(NEXTERNAL,LMAXCONFIGS)
    INCLUDE 'symperms.inc'
    REAL     XRAN1
    EXTERNAL XRAN1
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

    INTEGER IMIRROR
    COMMON/TO_MIRROR/ IMIRROR

%(define_iconfigs_lines)s
    DATA IDUM /0/
    DATA XTRY, XREJ /0,0/
    DATA NTRY /0,0/
    DATA NGOOD /0,0/
    DATA ISHEL /0,0/
    SAVE YFRAC, IGOOD, JHEL
    DATA GOODHEL/THEL*.FALSE./
%(helicity_lines)s
%(den_factor_line)s
C ----------
C BEGIN CODE
C ----------
    NTRY(IMIRROR)=NTRY(IMIRROR)+1

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
C     Reset permutation helper variables
      DO I=1,NEXTERNAL
        LPERM(I)=I
        TPERM(I)=I
      ENDDO
c     Get reversed config permutation
      CALL REVPERM(PERMS(1,MAPCONFIG(ICONFIG)),TPERM,NEXTERNAL)
      COMBANS=0d0
c     Calculate summed permutation weight
      TOTFACT=0
      DO IPERM=1,NPERMS
         IPERMFACT(IPERM)=PERMFACT%(proc_id)s(LPERM)
         TOTFACT=TOTFACT+IPERMFACT(IPERM)
C     Correction factor for events where no permutation is selected
         FRAC=MIN(NUMPERMS,NPERMS)*IPERMFACT(IPERM)*1D0/TOTFACT
         FACTNUL=FACTNUL*(1-FRAC)
      ENDDO
      FACTNUL=1-FACTNUL
c   Decide between helicity all-sum case (ISUM_HEL=0) and partial sum case
      JHELBK=JHEL(IMIRROR)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY(IMIRROR) .LE. MAXTRIES) THEN
         NUPPER=NCOMB
         HWGT=1.
      ELSE
         NUPPER=ISHEL(IMIRROR)
         HWGT = REAL(NGOOD(IMIRROR))/REAL(ISHEL(IMIRROR))
      ENDIF
 10 FIRSTFLOW=.TRUE.
c     Start loop over permutations
      DO IPERM=1,NPERMS
        JHEL(IMIRROR)=JHELBK
c     Set perm selection info
      FRAC=MIN(NUMPERMS,NPERMS)*IPERMFACT(IPERM)*1d0/TOTFACT
C     Decide whether to use this permutation
        IF(XRAN1(IDUM).GT.FRAC) CYCLE
c     Get permutation
        CALL GETPERM1(IPERM,LPERM)
c     Compensate for config permutation
        CALL SWITCHHEL(TPERM,PERM,LPERM,NEXTERNAL)
c     Start loop over helicities
        DO J=1,NUPPER
          IF (ISUM_HEL .EQ. 0 .OR. NTRY(IMIRROR) .LE. MAXTRIES) THEN
            I=J
            IF(.NOT.GOODHEL(I,IMIRROR) .AND. NTRY(IMIRROR) .GT. MAXTRIES) CYCLE
          ELSE
            JHEL(IMIRROR)=JHEL(IMIRROR)+1
            IF (JHEL(IMIRROR) .GT. NGOOD(IMIRROR)) JHEL(IMIRROR)=1
            I = IGOOD(JHEL(IMIRROR),IMIRROR)
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
        IF(NTRY(IMIRROR).LE.MAXTRIES)THEN
          DO I=1,NCOMB
            IF (.NOT.GOODHEL(I,IMIRROR) .AND. (TS(I).GT.ANS*LIMHEL/NCOMB)) THEN
              GOODHEL(I,IMIRROR)=.TRUE.
              NGOOD(IMIRROR) = NGOOD(IMIRROR) +1
              IGOOD(NGOOD(IMIRROR),IMIRROR) = I
              PRINT *,'Adding good helicity ',I,TS(I)/ANS
            ENDIF
          ENDDO
        ENDIF
        IF(NTRY(IMIRROR).EQ.MAXTRIES)THEN
          ISHEL(IMIRROR)=MIN(ISUM_HEL,NGOOD(IMIRROR))
        ENDIF
        IF (ISHEL(IMIRROR) .EQ. 1) THEN
          WRITE(HEL_BUFF,'(20i5)')(NHEL(II,I),II=1,NEXTERNAL)
        ENDIF
      ENDDO
C     If no permutation was chosen, repeat
      IF(FIRSTFLOW) GOTO 10
      IF (ISHEL(IMIRROR) .NE. 1) THEN
        R=XRAN1(IDUM)*ANS
        SUMHEL=0D0
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
          ANS=ANS*AMP2(SUBDIAG(1))/XTOT
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

C     ------------------------------------
      FUNCTION PERMFACT%(proc_id)s(PERM)
C     ------------------------------------
C     
C     Gives color factor for permutation PERM
C     
C     Process: g g > g g g g QCD=4 QED=0 WEIGHTED=4 singlet_QCD=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INCLUDE 'nexternal.inc'
      INCLUDE 'genps.inc'
      INTEGER NPAIRS
      PARAMETER(NPAIRS=%(npairs)d)
C     
C     ARGUMENTS 
C     
      INTEGER PERMFACT%(proc_id)s,PERM(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,IPAIRS
      LOGICAL OK
      INTEGER COMP(2,NPAIRS)
      DATA COMP/%(comp_pairs)s/
C     ----------
C     BEGIN CODE
C     ----------
      IPAIRS=0
      IF(NPAIRS.GT.0.AND.PERM(NEXTERNAL).EQ.COMP(2,1)) IPAIRS=IPAIRS+1
      DO I=1,NPAIRS
         DO J=1,NEXTERNAL-1
            IF(PERM(J).EQ.COMP(1,I).AND.PERM(J+1).EQ.COMP(2,I).OR.PERM(J).EQ.COMP(2,I).AND.PERM(J+1).EQ.COMP(1,I)) IPAIRS=IPAIRS+1
         ENDDO
      ENDDO
      PERMFACT%(proc_id)s=IFACT**IPAIRS
      RETURN
      END

