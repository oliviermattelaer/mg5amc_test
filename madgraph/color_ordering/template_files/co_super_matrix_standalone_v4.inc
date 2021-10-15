C     -------------------------
      SUBROUTINE SMATRIX(P,ANS)
C     -------------------------
C  
%(info_lines)s
C 
C MadGraph StandAlone Version
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
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
      INTEGER    NCOMB         
      PARAMETER (NCOMB=%(ncomb)d)
      INTEGER    ICO
      PARAMETER (ICO=%(color_order)d)
      INTEGER    NPERMS
      PARAMETER (NPERMS=%(nperms)d) 
      INTEGER    NJAMPS
      PARAMETER (NJAMPS=%(njampsAL)d)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL),ANS
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY, JPERM(NJAMPS,NPERMS)
      REAL*8 T
      REAL*8 MATRIX
      COMPLEX*16 JAMP(NJAMPS)
      INTEGER IHEL,IDEN, I
      LOGICAL GOODHEL(NCOMB)
      DATA NTRY/0/
      DATA GOODHEL/NCOMB*.FALSE./
%(helicity_lines)s
%(den_factor_line)s
%(jamp_perm_lines)s


C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      ANS = 0D0
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL) .OR. NTRY .LT. 2) THEN
               T=0d0
               %(call_flow_lines)s
               DO I=1,ICO
                  T=T+MATRIX(P , NHEL(1,IHEL), I, JPERM, JAMP)
               ENDDO
               ANS=ANS+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL)) THEN
                   GOODHEL(IHEL)=.TRUE.
               ENDIF
             ENDIF
          ENDDO
      ANS=ANS/DBLE(IDEN)
      END
       
C     -----------------------------------
      REAL*8 FUNCTION MATRIX(P, NHEL, ICO, JPERM, JAMP)
C     -----------------------------------
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
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
      INTEGER    NPERMS
      PARAMETER (NPERMS=%(nperms)d) 
      INTEGER    NFPERMS
      PARAMETER (NFPERMS=%(nflowperms)d) 
      INTEGER    NJAMPS
      PARAMETER (NJAMPS=%(njampsAL)d)
      
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL),ICO, JPERM(NJAMPS,NPERMS)
      COMPLEX*16 JAMP(NJAMPS)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J,PERM(NEXTERNAL)
      COMPLEX*16 ZTEMP
      INTEGER PMFLOW(NFPERMS,NPERMS)

      ZTEMP = (0.D0,0.D0)
%(color_sum_lines)s
      MATRIX=REAL(ZTEMP)

      RETURN
      END
