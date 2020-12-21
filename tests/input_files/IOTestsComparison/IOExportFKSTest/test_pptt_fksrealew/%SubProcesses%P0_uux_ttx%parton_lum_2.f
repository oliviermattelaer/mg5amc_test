      SUBROUTINE DLUM_2(LUM)
C     ****************************************************            
C         
C     Generated by MadGraph5_aMC@NLO v. %(version)s, %(date)s
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     RETURNS PARTON LUMINOSITIES FOR MADFKS                          
C        
C     
C     Process: g u~ > t t~ u~ [ real = QCD QED ] QCD^2=4 QED^2=2
C     Process: g c~ > t t~ c~ [ real = QCD QED ] QCD^2=4 QED^2=2
C     
C     ****************************************************            
C         
      IMPLICIT NONE
C     
C     CONSTANTS                                                       
C         
C     
      INCLUDE 'genps.inc'
      INCLUDE 'nexternal.inc'
      DOUBLE PRECISION       CONV
      PARAMETER (CONV=389379660D0)  !CONV TO PICOBARNS             
C     
C     ARGUMENTS                                                       
C         
C     
      DOUBLE PRECISION LUM
C     
C     LOCAL VARIABLES                                                 
C         
C     
      INTEGER I, ICROSS,LP
      DOUBLE PRECISION G1
      DOUBLE PRECISION CX2,UX2
C     
C     EXTERNAL FUNCTIONS                                              
C         
C     
      DOUBLE PRECISION PDG2PDF
C     
C     GLOBAL VARIABLES                                                
C         
C     
      INTEGER              IPROC
      DOUBLE PRECISION PD(0:MAXPROC)
      COMMON /SUBPROC/ PD, IPROC
      INCLUDE 'coupl.inc'
      INCLUDE 'run.inc'
      INTEGER IMIRROR
      COMMON/CMIRROR/IMIRROR
C     
C     DATA                                                            
C         
C     
      DATA G1/1*1D0/
      DATA CX2,UX2/2*1D0/
      DATA ICROSS/1/
C     ----------                                                      
C         
C     BEGIN CODE                                                      
C         
C     ----------                                                      
C         
      LUM = 0D0
      IF (ABS(LPP(1)) .GE. 1) THEN
        LP=SIGN(1,LPP(1))
        G1=PDG2PDF(ABS(LPP(1)),0*LP,XBK(1),DSQRT(Q2FACT(1)))
      ENDIF
      IF (ABS(LPP(2)) .GE. 1) THEN
        LP=SIGN(1,LPP(2))
        CX2=PDG2PDF(ABS(LPP(2)),-4*LP,XBK(2),DSQRT(Q2FACT(2)))
        UX2=PDG2PDF(ABS(LPP(2)),-2*LP,XBK(2),DSQRT(Q2FACT(2)))
      ENDIF
      PD(0) = 0D0
      IPROC = 0
      IPROC=IPROC+1  ! g u~ > t t~ u~
      PD(IPROC) = G1*UX2
      IPROC=IPROC+1  ! g c~ > t t~ c~
      PD(IPROC) = G1*CX2
      DO I=1,IPROC
        IF (NINCOMING.EQ.2) THEN
          LUM = LUM + PD(I) * CONV
        ELSE
          LUM = LUM + PD(I)
        ENDIF
      ENDDO
      RETURN
      END

