      SUBROUTINE FLOW%(number)s(P,NHEL,JAMP)
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
      INTEGER    NPERMS
      PARAMETER (NPERMS=%(nperms)d) 
      INTEGER    NJAMPS
      PARAMETER (NJAMPS=%(njampsAL)d)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1, ONE
      PARAMETER (IMAG1=(0D0,1D0),ONE=(1d0,0d0))
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL)
      COMPLEX*16 JAMP(NJAMPS)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J,IC(NEXTERNAL)
      COMPLEX*16 AMPN(NGRAPHS)
      COMPLEX*16 WN(18,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0d0, 0d0), (1d0, 0d0)/
%(fake_width_declaration)s
      logical first
      data first /.true./
      save first
%(ic_data_line)s
C  
C GLOBAL VARIABLES
C  
      include 'coupl.inc'

	double precision small_width_treatment
        common/narrow_width/small_width_treatment
C ----------
C BEGIN CODE
C ----------
if (first) then
    first=.false.
    %(fake_width_definitions)s
endif 
%(helas_calls)s
%(jamp_lines)s

      END
