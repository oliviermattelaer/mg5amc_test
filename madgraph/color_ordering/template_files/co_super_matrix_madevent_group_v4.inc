C     -------------------------
      SUBROUTINE SMATRIX%(proc_id)s(P,ANS)
C     -------------------------
C  
%(info_lines)s
C 
C Color-ordered MadGraph for Madevent Version
C 
C Returns amplitude squared with sampling over colors and helicities
C for the point in phase space P(0:3,NEXTERNAL)
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
    INTEGER    ICO
    PARAMETER (ICO=%(color_order)d)
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
    REAL*8 R,SUMHEL,TS(NCOMB),TSCF(NCOMB,NPERMS)
    INTEGER I,IDEN,T_IDENT(NCOMB,2)
    INTEGER IPROC,II
    LOGICAL GOODHEL(NCOMB,2)
    REAL*8 HWGT,XTOT,XTRY,XREJ,XR,YFRAC(0:NCOMB),FRAC
    INTEGER IDUM,IPERM,NGOOD(2),IGOOD(NCOMB,2)
    INTEGER JHEL(2),J,JJ,ISEL
    REAL*8 TOTFACT
    INTEGER JHELBK,NUPPER
    INTEGER NSELPERM(NPERMS),NSORTPERM(NPERMS),NSELPERMS
    DOUBLE PRECISION JAC(NPERMS)
    LOGICAL FIRSTFLOW
    REAL*8   RVEC
    INTEGER ICOUNT
    SAVE FRAC,GOODHEL,T_IDENT,ICOUNT
    DATA ICOUNT/0/
    DOUBLE PRECISION ANSCF(NPERMS)
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

    LOGICAL CF_BASIS(NDIAGS,NPERMS)
    COMMON/TO_CF_BASIS/CF_BASIS

%(define_iconfigs_lines)s
    DATA IDUM /0/
    DATA XTRY, XREJ /0,0/
    DATA NTRY /0,0/
    DATA NGOOD /0,0/
    DATA ISHEL /0,0/
    SAVE YFRAC, IGOOD, JHEL
    DATA GOODHEL/THEL*.FALSE./
    DATA T_IDENT/THEL*-1/
%(helicity_lines)s
%(den_factor_line)s
C ----------
C BEGIN CODE
C ----------
    NTRY(IMIRROR)=NTRY(IMIRROR)+1

    IF (NTRY(IMIRROR).eq.1) THEN
        call generate_CF_integration_basis()
    ENDIF    
    
    IF (multi_channel) THEN
        DO I=1,NDIAGS
            AMP2(I)=0D0
        ENDDO
        JAMP2(0)=NFLOWS
        DO I=1,INT(JAMP2(0))
            JAMP2(I)=0D0
        ENDDO
    ENDIF
    DO I=1,NPERMS
        ANSCF(I) = 0d0
    ENDDO   
    ANS = 0D0
    WRITE(HEL_BUFF,'(20I5)') (0,I=1,NEXTERNAL)
    DO I=1,NCOMB
       TS(I)=0d0
       DO J=1,NPERMS
          TSCF(I,J)=0d0
       ENDDO 
    ENDDO
    NSELPERMS=MIN(NUMPERMS,NPERMS)
    IF(NSELPERMS.EQ.0) NSELPERMS=NPERMS
    IF(NTRY(IMIRROR).EQ.1)THEN
c     Calculate summed permutation weight
      TOTFACT=NPERMS
c      FRAC=NSELPERMS*1D0/TOTFACT
    ENDIF
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
c     Choose permutations completely randomly ! need to be improve!
      CALL SELECT_COLOR_FLOW(NSELPERMS, .true., SUBDIAG(1), NSELPERM, JAC)
c     Start loop over chosen permutations
      DO ISEL=1,NSELPERMS
c     Reset helicity info
        JHEL(IMIRROR)=JHELBK
c     Get permutation
        IPERM=NSELPERM(ISEL)
        CALL GETPERM%(proc_id)s(IPERM,PERM)
        IF(ICO.GT.1)THEN
C       Update ICOUNT to keep track of NLC contributions
           ICOUNT=ICOUNT+1
C       Reset if above highest relevant value
           IF(ICOUNT.EQ.NLCFACT**(ICO-1)+1) ICOUNT=1
        ENDIF
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
          IF(T_IDENT(I,IMIRROR).LT.0) THEN
C         Calculate matrix element if no identical helicity
             T=MATRIX%(proc_id)s(P,NHEL(1,I),PERM,FIRSTFLOW,1)
          ELSE IF(TS(T_IDENT(I,IMIRROR)).EQ.0)THEN
             T=MATRIX%(proc_id)s(P,NHEL(1,I),PERM,FIRSTFLOW,1)
          ENDIF
          IF(ICO.GT.1.AND.MOD(ICOUNT,NLCFACT).EQ.0) THEN
C         Sample NLC contributions
            T=T+NLCFACT*MATRIX%(proc_id)s(P,NHEL(1,I),PERM,FIRSTFLOW,2)
            IF(ICO.GT.2.AND.MOD(ICOUNT,NLCFACT**2).EQ.0) THEN
C         Sample NNLC and higher contributions
              DO JJ=3,ICO
                T=T+NLCFACT**2*MATRIX%(proc_id)s(P,NHEL(1,I),PERM,FIRSTFLOW,JJ)
              ENDDO
            ENDIF
          ENDIF
          DO JJ=1,NINCOMING
            IF(POL(JJ).NE.1D0.AND.NHEL(JJ,I).EQ.INT(SIGN(1D0,POL(JJ)))) THEN
              T=T*ABS(POL(JJ))
            ELSE IF(POL(JJ).NE.1D0)THEN
              T=T*(2D0-ABS(POL(JJ)))
            ENDIF
          ENDDO
          TSCF(I,IPERM) = TSCF(I,IPERM) + T*HWGT*JAC(ISEL)
        ENDDO
        FIRSTFLOW=.FALSE.
C       Calculate ANSCF 
c       First reset helicity info
c        JHEL(IMIRROR)=JHELBK
        DO J=1,NUPPER
          IF (ISUM_HEL .EQ. 0 .OR. NTRY(IMIRROR) .LE. MAXTRIES) THEN
            I=J
            IF(.NOT.GOODHEL(I,IMIRROR) .AND. NTRY(IMIRROR) .GT. MAXTRIES) CYCLE
          ELSE
            JHEL(IMIRROR)=JHEL(IMIRROR)+1
            IF (JHEL(IMIRROR) .GT. NGOOD(IMIRROR)) JHEL(IMIRROR)=1
            I = IGOOD(JHEL(IMIRROR),IMIRROR)
          ENDIF
          IF(T_IDENT(I,IMIRROR).GT.0.AND.TSCF(MAX(1,T_IDENT(I,IMIRROR)),IPERM).NE.0)THEN                                                                                                                    
C           If identical helicity, set TS(I) to that value                                                                                                                                                    
            TSCF(I,IPERM)=TSCF(T_IDENT(I,IMIRROR),IPERM)
          ENDIF
          ANSCF(IPERM)=ANSCF(IPERM)+DABS(TSCF(I,IPERM))
        ENDDO
C     Determine which helicities to include
      IF(NTRY(IMIRROR).LE.MAXTRIES)THEN
        DO I=1,NCOMB
          IF (.NOT.GOODHEL(I,IMIRROR) .AND. (DABS(TSCF(I,IPERM)).GT.ANSCF(IPERM)*LIMHEL/NCOMB)) THEN
            GOODHEL(I,IMIRROR)=.TRUE.
            NGOOD(IMIRROR) = NGOOD(IMIRROR) +1
            IGOOD(NGOOD(IMIRROR),IMIRROR) = I
C           Look for identical helicities
            DO J=1,I-1
              IF (TSCF(J,IPERM).ne.0d0)then
              IF(DABS(TSCF(I,IPERM)/TSCF(J,IPERM)-1D0).LT.1D-12)THEN
                PRINT *,'Found identical helicity ',I,J,TSCF(I,IPERM)/ANSCF(IPERM)
                T_IDENT(I,IMIRROR)=J
                EXIT
              ENDIF
              endif
            ENDDO
            IF(T_IDENT(I,IMIRROR).LT.0)THEN
              PRINT *,'Added good helicity ',I,TSCF(I,IPERM)/ANSCF(IPERM), ' in event ',NTRY(IMIRROR)
            ENDIF
          ENDIF
        ENDDO
      ENDIF        
C     Set how many helicities to sum over
      IF(NTRY(IMIRROR).EQ.MAXTRIES)THEN
        ISHEL(IMIRROR)=MIN(ISUM_HEL,NGOOD(IMIRROR))
      ENDIF


      IF (MULTI_CHANNEL) THEN
          IF (CF_BASIS(SUBDIAG(1),IPERM).AND.AMP2(SUBDIAG(1)).GT.0D0)THEN
            XTOT=0D0
            DO J=1,NDIAGS
              IF (CF_BASIS(J,IPERM))THEN
                XTOT = XTOT + AMP2(J)
              ENDIF
            ENDDO

            ANSCF(IPERM) = ANSCF(IPERM)*AMP2(SUBDIAG(1))/XTOT
          ELSE
            ANSCF(IPERM) = 0D0
          ENDIF
      ENDIF
      ENDDO        
      DO I = 1,NPERMS
        ANS = ANS +ANSCF(I)
      ENDDO
      do i = 1,NCOMB
         do j = 1,NPERMS
            TS(i) = TS(i) +TSCF(i,j)
         enddo
      enddo        
        

c     Choose a helicity to write out
      IF (ISHEL(IMIRROR) .EQ. 1) THEN
         WRITE(HEL_BUFF,'(20i5)')(NHEL(II,I),II=1,NEXTERNAL)
c         Set right sign for ANS, based on sign of chosen helicity
         ANS=DSIGN(ANS,TS(I))
      ELSE
        CALL RANMAR(RVEC)
        SUMHEL=0D0
        DO I=1,NCOMB
          SUMHEL=SUMHEL+DABS(TS(I))/ANS
          IF(RVEC.LT.SUMHEL)THEN
            WRITE(HEL_BUFF,'(20i5)')(NHEL(II,I),II=1,NEXTERNAL)
c           Set right sign for ANS, based on sign of chosen helicity
            ANS=DSIGN(ANS,TS(I))
            GOTO 20
          ENDIF
        ENDDO
 20     CONTINUE
      ENDIF
      ANS=ANS/DBLE(IDEN)
      END

C     ------------------------------------------------------------
      REAL*8 FUNCTION MATRIX%(proc_id)s(P,NHEL,PERM,FIRSTFLOW,ICO)
C     ------------------------------------------------------------
C
C Returns AMP2s and matrix element for a given helicity and permutation
c (color flow) at the color order specified by ICO
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
      INTEGER NHEL(NEXTERNAL),PERM(NEXTERNAL),ICO
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
      ZTEMP=ZTEMP+ONEPERM%(proc_id)s(P,NHEL,PERM,ICO)
      MATRIX%(proc_id)s=REAL(ZTEMP)

      RETURN
      END
      
C     -----------------------------------------------------
      COMPLEX*16 FUNCTION ONEPERM%(proc_id)s(P,NHEL,PM,ICO)
C     -----------------------------------------------------
C
C Returns the color flows for a given helicity and permutation, 
C at color order given by ICO
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
      INTEGER ICO
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

    SUBROUTINE SELECT_COLOR_FLOW(NB_TO_SELECT, FLAT, CHANNEL, NSELPERM, JAC)
    IMPLICIT NONE  
c     PARAMETER      
      INTEGER    NPERMS
      PARAMETER (NPERMS=%(nperms)d)
      INTEGER    NDIAGS
      PARAMETER (NDIAGS=%(ndiags)d)  
c     INPUT
      INTEGER NB_TO_SELECT ! can be change if impossible to reach taht number.
      LOGICAL FLAT
      INTEGER CHANNEL ! current diagram under consideration
c     OUTPUT
      INTEGER NSELPERM(*)
      DOUBLE PRECISION JAC(NPERMS)
c     LOCAL
      INTEGER ISEL,I,J,K
      INTEGER  TOT_PICK, ALLOWED
      DOUBLE PRECISION TOTAL,TMP
      DOUBLE PRECISION LIMIT(0:NPERMS)
      DOUBLE PRECISION RANDOM
      DOUBLE PRECISION LOC_JAC       
c     GLOBAL
      LOGICAL CF_BASIS(NDIAGS, NPERMS)
      COMMON/TO_CF_BASIS/CF_BASIS  
c
      INTEGER FACTORIAL
      include 'maxamps.inc' ! define MAXAMPS and MAXFLOW       
      DOUBLE PRECISION AMP2(MAXAMPS),JAMP2(0:MAXFLOW)
      COMMON/TO_AMPS/  AMP2, JAMP2       
       
      TOTAL=0d0
      JAC=1d0
c      
c    PREPARE GRID FOR PICKING COLOR FLOW RANDOMLY
c   
      LIMIT(0) = 0d0
      ALLOWED=0
      IF (FLAT) THEN
        DO I=1,NPERMS
           IF (CF_BASIS(CHANNEL,I)) THEN
                ALLOWED = ALLOWED +1
                TOTAL = TOTAL + 1d0
           ENDIF
           LIMIT(I) = TOTAL  
        ENDDO
c
c     WEIGHTED CASE
c        
      ELSE
        DO I=1,NPERMS
           TMP = 0d0
           IF (CF_BASIS(CHANNEL, I)) THEN
                DO J=1,NDIAGS
                    IF (CF_BASIS(J, I)) THEN
                        TMP = TMP + AMP2(J)
                    ENDIF
                ENDDO
           ENDIF
           IF (TMP.ne.0d0) THEN
               TOTAL = TOTAL + 1d0/TMP
               ALLOWED = ALLOWED + 1
           ENDIF     
           LIMIT(I) = TOTAL  
        ENDDO
      ENDIF
       
      
      IF (ALLOWED.le.NB_TO_SELECT)then
c        select all
         ISEL = 0
         DO I=1,NPERMS
           if (LIMIT(I).ne.LIMIT(I-1))then
              ISEL = ISEL +1
              NSELPERM(ISEL)=I
              JAC(ISEL) =1 
           endif        
        ENDDO
        DO I= ISEL+1, NB_TO_SELECT
            NSELPERM(I) = 0
            JAC(I) = 0
        ENDDO
        NB_TO_SELECT = ALLOWED
      ELSE
       
c     Choose permutations 
      LOC_JAC =1d0
      TOT_PICK = 0
      DO ISEL=1,NB_TO_SELECT
c     Pick a permutation
 20        TOT_PICK = TOT_PICK + 1
           CALL RANMAR(RANDOM)
           DO J=1,NPERMS
             IF ((LIMIT(J)/TOTAL)>RANDOM) THEN
                LOC_JAC =  TOTAL/(LIMIT(J)-LIMIT(J-1))
                do k=1,ISEL-1
                    if (NSELPERM(K).eq.J) then
                        JAC(K) = JAC(k) + LOC_JAC
                        goto 20
                    endif   
                enddo     
                NSELPERM(ISEL)=J
                JAC(ISEL) = LOC_JAC
c               FORBID THIS PICK TO BE TAKEN A SECOND TIME.                
c                TMP = LIMIT(J) - LIMIT(J-1)
c                TOTAL = TOTAL - TMP
c                DO K=J,NPERMS
c                    LIMIT(K) = LIMIT(K) - TMP
c                ENDDO
                GOTO 10
             ENDIF
         ENDDO
         IF (LIMIT(NPERMS).eq.0d0) THEN
             NSELPERM(ISEL)=0d0
             NB_TO_SELECT = ISEL - 1
             TOT_PICK = TOT_PICK - 1
             goto 30
         ELSE
             STOP 1
         ENDIF
 10    ENDDO 
c         
c     COMPUTE GLOBAL JAC FACTOR         
c
 30    DO ISEL=1,NB_TO_SELECT
         JAC(ISEL) = JAC(ISEL) / TOT_PICK
     ENDDO
     ENDIF
          
      RETURN
    END             
