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
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL),ANS
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T
      REAL*8 MATRIX
      COMPLEX*16 JAMP(NPERMS)
      INTEGER IHEL,IDEN, I
      LOGICAL GOODHEL(NCOMB)
      DATA NTRY/0/
      DATA GOODHEL/NCOMB*.FALSE./
%(helicity_lines)s
%(den_factor_line)s


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
                  T=T+MATRIX(P , NHEL(1,IHEL), I, JAMP)
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
      REAL*8 FUNCTION MATRIX(P, NHEL, ICO, JAMP)
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
      
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL),ICO
      COMPLEX*16 JAMP(NPERMS)
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
      
C     --------------------------------------
      COMPLEX*16 FUNCTION ONEPERM(P,NHEL,PM,ICO)
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
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
      INTEGER    NFLOWS
      PARAMETER (NFLOWS=%(nflows)d) 
      INTEGER    NFPERMS
      PARAMETER (NFPERMS=%(nflowperms)d) 
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL)
      INTEGER PM(NEXTERNAL)
      INTEGER ICO
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      COMPLEX*16 JAMP(NFLOWS)
      INTEGER PERMS(NEXTERNAL,NFPERMS),IFERM(NFPERMS),PERM(NEXTERNAL)
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
      ONEPERM=ZTEMP

      RETURN
      END
