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
      INTEGER                 NCOMB         
      PARAMETER (             NCOMB=%(ncomb)d)
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
                 T=MATRIX(P ,NHEL(1,IHEL))            
               ANS=ANS+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL)) THEN
                   GOODHEL(IHEL)=.TRUE.
               ENDIF
             ENDIF
          ENDDO
      ANS=ANS/DBLE(IDEN)
      END
       
C     ------------------------------
      REAL*8 FUNCTION MATRIX(P,NHEL)
C     ------------------------------
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
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J,PERM(NEXTERNAL)
      COMPLEX*16 ZTEMP
C
C EXTERNAL FUNCTIONS
C
      COMPLEX*16 ONEPERM
      EXTERNAL ONEPERM

      ZTEMP = (0.D0,0.D0)
      DO I=1,NPERMS
        CALL GETPERM(I,PERM)
	ZTEMP=ZTEMP+ONEPERM(P,NHEL,PERM)
      ENDDO
      MATRIX=REAL(ZTEMP)

      RETURN
      END
      
C     --------------------------------------
      COMPLEX*16 FUNCTION ONEPERM(P,NHEL,PM)
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
      INTEGER    NPERMS
      PARAMETER (NPERMS=%(nflowperms)d) 
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL)
      INTEGER PM(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      COMPLEX*16 JAMP(NFLOWS)
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
      ONEPERM=ZTEMP

      RETURN
      END

C     ------------------------------------
      SUBROUTINE GETPERM(IPERM,PERM)
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
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
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
