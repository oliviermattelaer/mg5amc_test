      SUBROUTINE DLUM_1(LUM)
C     ****************************************************            
C         
C     Generated by MadGraph5_aMC@NLO v. %(version)s, %(date)s
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     RETURNS PARTON LUMINOSITIES FOR MADFKS                          
C        
C     
C     Process: u~ u > t t~ g [ real = QCD QED ] QCD^2<=6 QED^2<=0
C     Process: c~ c > t t~ g [ real = QCD QED ] QCD^2<=6 QED^2<=0
C     Process: d~ d > t t~ g [ real = QCD QED ] QCD^2<=6 QED^2<=0
C     Process: s~ s > t t~ g [ real = QCD QED ] QCD^2<=6 QED^2<=0
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
      DOUBLE PRECISION CX1,SX1,UX1,DX1
      DOUBLE PRECISION D2,U2,S2,C2
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
      DATA CX1,SX1,UX1,DX1/4*1D0/
      DATA D2,U2,S2,C2/4*1D0/
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
        CX1=PDG2PDF(ABS(LPP(1)),-4*LP,XBK(1),DSQRT(Q2FACT(1)))
        SX1=PDG2PDF(ABS(LPP(1)),-3*LP,XBK(1),DSQRT(Q2FACT(1)))
        UX1=PDG2PDF(ABS(LPP(1)),-2*LP,XBK(1),DSQRT(Q2FACT(1)))
        DX1=PDG2PDF(ABS(LPP(1)),-1*LP,XBK(1),DSQRT(Q2FACT(1)))
      ENDIF
      IF (ABS(LPP(2)) .GE. 1) THEN
        LP=SIGN(1,LPP(2))
        D2=PDG2PDF(ABS(LPP(2)),1*LP,XBK(2),DSQRT(Q2FACT(2)))
        U2=PDG2PDF(ABS(LPP(2)),2*LP,XBK(2),DSQRT(Q2FACT(2)))
        S2=PDG2PDF(ABS(LPP(2)),3*LP,XBK(2),DSQRT(Q2FACT(2)))
        C2=PDG2PDF(ABS(LPP(2)),4*LP,XBK(2),DSQRT(Q2FACT(2)))
      ENDIF
      PD(0) = 0D0
      IPROC = 0
      IPROC=IPROC+1  ! u~ u > t t~ g
      PD(IPROC) = UX1*U2
      IPROC=IPROC+1  ! c~ c > t t~ g
      PD(IPROC) = CX1*C2
      IPROC=IPROC+1  ! d~ d > t t~ g
      PD(IPROC) = DX1*D2
      IPROC=IPROC+1  ! s~ s > t t~ g
      PD(IPROC) = SX1*S2
      DO I=1,IPROC
        IF (NINCOMING.EQ.2) THEN
          LUM = LUM + PD(I) * CONV
        ELSE
          LUM = LUM + PD(I)
        ENDIF
      ENDDO
      RETURN
      END

