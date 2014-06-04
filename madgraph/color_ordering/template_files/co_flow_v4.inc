      COMPLEX*16 FUNCTION FLOW%(number)s(P,NHEL,IP)
C  
%(info_lines)s
C
C Returns color ordered amplitude for color flow %(number)s
c for the point with external momenta P(0:3, NEXTERNAL)
C  
%(process_lines)s
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=%(ngraphs)d) 
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=%(nwavefuncs)d) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1, ONE
      PARAMETER (IMAG1=(0D0,1D0),ONE=(1d0,0d0))
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL),IP(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J,IC(NEXTERNAL)
      COMPLEX*16 AMP(NGRAPHS), JAMP(1)
      COMPLEX*16 W(18,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0d0, 0d0), (1d0, 0d0)/
%(ic_data_line)s
C  
C GLOBAL VARIABLES
C  
      include 'coupl.inc'
C ----------
C BEGIN CODE
C ----------
%(helas_calls)s
%(jamp_lines)s

      FLOW%(number)s=JAMP(1)
      RETURN
      END