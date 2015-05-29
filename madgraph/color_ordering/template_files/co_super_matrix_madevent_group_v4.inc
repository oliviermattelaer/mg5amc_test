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
      use DiscreteSampler
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
C       This is just to temporarily store the reference grid for helicity of the DiscreteSampler so as to obtain its number of entries with ref_helicity_grid%%n_tot_entries
    type(SampledDimension) ref_cf_grid
    
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
    
C     To be able to control when the matrix<i> subroutine can add                                                                                                                                           
C      entries to the grid for the MC over helicity configuration                                                                                                                                           
    LOGICAL ALLOW_HELICITY_GRID_ENTRIES
    COMMON/TO_ALLOW_HELICITY_GRID_ENTRIES/ALLOW_HELICITY_GRID_ENTRIES

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
    
!   If the CF grid status is 0, this means that it is not yet initialized.
    IF(DS_get_dim_status('color_flow').eq.0.or.cf_picked.eq.-1) THEN
        DO I=1, NPERMS ! loop over color-flow
           T=0
           DO J=1,NCOMB ! loop over helicity        
c           Get permutation/color flow
            CALL GETPERM%(proc_id)s(I, PERM)
            TS(J)=MATRIX%(proc_id)s(P,NHEL(1,J),PERM,.FALSE.,1)

            DO JJ=1,nincoming
               IF(POL(JJ).NE.1d0)THEN
                    IF(NHEL(JJ,J).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                        TS(J)=TS(J)*ABS(POL(JJ))
                    ELSE 
                        TS(J)=TS(J)*(2d0-ABS(POL(JJ)))
                    ENDIF
               ENDIF
            ENDDO
            if (ALLOW_HELICITY_GRID_ENTRIES.and.ISUM_HEL.ne.0) then
                call DS_add_entry('Helicity',J,TS(J))
                call DS_add_entry('color_flow',I,TS(J))
            endif
            T = T + TS(J)
           ENDDO
           if (ISUM_HEL.eq.0) then
               call DS_add_entry('color_flow',I, T)
c              check for irrrelevant helicity               
               do J=1,NCOMB
               IF (.NOT.GOODHEL(J,IMIRROR) .AND. (DABS(TS(J)).GT.T*LIMHEL/NCOMB)) THEN
                   GOODHEL(J,IMIRROR) = .True.
               endif
               enddo
c              Look for identical helicities
c               DO JJ = 1, NCOMB
c                    DO J=1,I-1
c                        IF (TS(J).ne.0d0.)then
c                            IF(DABS(TS(JJ)/TS(J)-1D0).LT.1D-12)THEN
c                            PRINT *,'Found identical helicity ',I,J,TS(JJ)/T
c                            T_IDENT(JJ,IMIRROR)=J
c                         ENDIF
c                     ENDDO
c               ENDDO   
           endif
        ENDDO
          CF_PICKED = -1
          HEL_PICKED = -1
          ANS=1d0
          RETURN
    ENDIF
    ref_cf_grid = DS_get_dimension(ref_grid,'color_flow')
    IF((DS_get_dim_status('color_flow').eq.1).and.(ref_cf_grid%%n_tot_entries.eq.0)) then 
!          If we finished the initialization we can update the grid so as to start sampling over it.
!           However the grid will now be filled by dsample with different kind of weights (including pdf, flux, etc...) so by setting the grid_mode of the reference grid to 'initialization' we make sure it will be overwritten (as opposed to 'combined') by the running grid at the next update.
       write(*,*) "calling the first(?) update of the grid"
       if (ISUM_HEL.ne.0)then
            CALL DS_UPDATE_GRID('Helicity')
            CALL DS_SET_GRID_MODE('Helicity','init')
       endif
       CALL DS_UPDATE_GRID('color_flow')
       CALL DS_SET_GRID_MODE('color_flow','init')
       call reset_cumulative_variable() ! avoid biais of the initialization
    endif 
    
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
    
    NSELPERMS = 1 ! For the moment only one CF by PS point
    
c    NSELPERMS=MIN(NUMPERMS,NPERMS)
c    IF(NSELPERMS.EQ.0) NSELPERMS=NPERMS
c    IF(NTRY(IMIRROR).EQ.1)THEN
c     Calculate summed permutation weight
c      TOTFACT=NPERMS
c      FRAC=NSELPERMS*1D0/TOTFACT
c    ENDIF
c   Decide between helicity all-sum case (ISUM_HEL=0) and partial sum case
c      JHELBK=JHEL(IMIRROR)
      IF (ISUM_HEL .EQ. 0) THEN
         NUPPER=NCOMB
         HWGT=1
      ELSE
         NUPPER=1
         HWGT = hel_jacobian
      ENDIF
 10 FIRSTFLOW=.TRUE.
c     Color-flow already selected!
c      CALL SELECT_COLOR_FLOW(NSELPERMS, .true., SUBDIAG(1), NSELPERM, JAC)
c     Start loop over chosen permutations ! only one for the moment.
      NSELPERM(1) = CF_PICKED
      JAC(1) = cf_jacobian
      DO ISEL=1,NSELPERMS
c     Reset helicity info
        JHEL(IMIRROR)= HEL_PICKED
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
          IF (ISUM_HEL .EQ. 0) THEN
            I= J
            IF(.NOT.GOODHEL(I,IMIRROR)) CYCLE
          ELSE
            I = HEL_PICKED
          ENDIF
C         Calculate matrix element
          T=MATRIX%(proc_id)s(P,NHEL(1,I),PERM,FIRSTFLOW,1)
C         partial probing of subleading color only 1/NLCFACT point used         
          IF(ICO.GT.1.AND.MOD(ICOUNT,NLCFACT).EQ.0) THEN
                T=T+NLCFACT*MATRIX%(proc_id)s(P,NHEL(1,I),PERM,FIRSTFLOW,2)
                IF(ICO.GT.2.AND.MOD(ICOUNT,NLCFACT**2).EQ.0) THEN
                DO JJ=3,ICO
                 T=T+NLCFACT**2*MATRIX%(proc_id)s(P,NHEL(1,I),PERM,FIRSTFLOW,JJ)
                ENDDO
                ENDIF
          ENDIF
C         END COMPUTE MATRIX ELEMENT
          DO JJ=1,NINCOMING
            IF(POL(JJ).NE.1D0.AND.NHEL(JJ,I).EQ.INT(SIGN(1D0,POL(JJ)))) THEN
              T=T*ABS(POL(JJ))
            ELSE IF(POL(JJ).NE.1D0)THEN
              T=T*(2D0-ABS(POL(JJ)))
            ENDIF
          ENDDO
          TSCF(I,IPERM) = TSCF(I,IPERM) + T*HWGT*JAC(ISEL)
          ANSCF(IPERM) = ANSCF(IPERM) + T*HWGT*JAC(ISEL)
        ENDDO ! end loop over helicities
        FIRSTFLOW=.FALSE.
        
c       First reset helicity info
c        JHEL(IMIRROR)=JHELBK
c        DO J=1,NUPPER
c          IF (ISUM_HEL .EQ. 0 .OR. NTRY(IMIRROR) .LE. MAXTRIES) THEN
c            I=J
c            IF(.NOT.GOODHEL(I,IMIRROR) .AND. NTRY(IMIRROR) .GT. MAXTRIES) CYCLE
c          ELSE
c            JHEL(IMIRROR)=JHEL(IMIRROR)+1
c            IF (JHEL(IMIRROR) .GT. NGOOD(IMIRROR)) JHEL(IMIRROR)=1
c            I = IGOOD(JHEL(IMIRROR),IMIRROR)
c          ENDIF
c          IF(T_IDENT(I,IMIRROR).GT.0.AND.TSCF(MAX(1,T_IDENT(I,IMIRROR)),IPERM).NE.0)THEN                                                                                                                    
C           If identical helicity, set TS(I) to that value                                                                                                                                                    
c            TSCF(I,IPERM)=TSCF(T_IDENT(I,IMIRROR),IPERM)
c          ENDIF
c          ANSCF(IPERM)=ANSCF(IPERM)+DABS(TSCF(I,IPERM))
c        ENDDO
        
C     Set how many helicities to sum over
c      IF(NTRY(IMIRROR).EQ.MAXTRIES)THEN
c        ISHEL(IMIRROR)=MIN(ISUM_HEL,NGOOD(IMIRROR))
c      ENDIF

      IF (MULTI_CHANNEL) THEN
          IF (AMP2(SUBDIAG(1)).GT.0D0)THEN
            XTOT=0D0
            DO J=1,NDIAGS
                XTOT = XTOT + AMP2(J)
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
      IF (HEL_PICKED.ne.-1) THEN
         WRITE(HEL_BUFF,'(20i5)')(NHEL(II,HEL_PICKED),II=1,NEXTERNAL)
c         Set right sign for ANS, based on sign of chosen helicity
         ANS=DSIGN(ANS,TS(HEL_PICKED))      
      ELSEIF (ISHEL(IMIRROR) .EQ. 1) THEN
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
      ANS = ANS
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
           
