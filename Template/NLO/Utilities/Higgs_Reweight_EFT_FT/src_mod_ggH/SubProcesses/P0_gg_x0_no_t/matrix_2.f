      SUBROUTINE SMATRIX_2(P,ANS)
C     
C     Generated by MadGraph5_aMC@NLO v. 2.2.3, 2015-02-10
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     Returns amplitude squared summed/avg over colors
C     and helicities
C     for the point in phase space P(0:3,NEXTERNAL)
C     
C     Process: d~ g > x0 d~ XGLU=1 WEIGHTED=3 QNP=1 [ real = QCD ] / t
C     Process: u~ g > x0 u~ XGLU=1 WEIGHTED=3 QNP=1 [ real = QCD ] / t
C     Process: s~ g > x0 s~ XGLU=1 WEIGHTED=3 QNP=1 [ real = QCD ] / t
C     Process: c~ g > x0 c~ XGLU=1 WEIGHTED=3 QNP=1 [ real = QCD ] / t
C     Process: b~ g > x0 b~ XGLU=1 WEIGHTED=3 QNP=1 [ real = QCD ] / t
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INCLUDE 'nexternal.inc'
      INCLUDE '../../Source/MODEL/input.inc'
      INTEGER     NCOMB
      PARAMETER ( NCOMB=8)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL),ANS
      DOUBLE PRECISION       WGT_ME_BORN,WGT_ME_REAL
      COMMON /C_WGT_ME_TREE/ WGT_ME_BORN,WGT_ME_REAL
C     
C     LOCAL VARIABLES 
C     
      INTEGER IHEL,IDEN,I,T_IDENT(NCOMB)
      REAL*8 MATRIX_2
      REAL*8 T,T_SAVE(NCOMB)
      SAVE T_SAVE,T_IDENT
      INTEGER NHEL(NEXTERNAL,NCOMB)
      DATA (NHEL(I,   1),I=1,4) /-1,-1, 0, 1/
      DATA (NHEL(I,   2),I=1,4) /-1,-1, 0,-1/
      DATA (NHEL(I,   3),I=1,4) /-1, 1, 0, 1/
      DATA (NHEL(I,   4),I=1,4) /-1, 1, 0,-1/
      DATA (NHEL(I,   5),I=1,4) / 1,-1, 0, 1/
      DATA (NHEL(I,   6),I=1,4) / 1,-1, 0,-1/
      DATA (NHEL(I,   7),I=1,4) / 1, 1, 0, 1/
      DATA (NHEL(I,   8),I=1,4) / 1, 1, 0,-1/
      LOGICAL GOODHEL(NCOMB)
      DATA GOODHEL/NCOMB*.FALSE./
      INTEGER NTRY
      DATA NTRY/0/
      DATA IDEN/96/
      double precision AMPqg,AMPLO,sman,tman,uman,dot,pi
     f ,tpt1,tpt2,upt1,upt2,mh2 
     f ,AMPLOtotop 
      parameter (pi=3.1415926535897932385d0)
      INCLUDE 'coupl.inc'

c$$$C     ----------
c$$$C     BEGIN CODE
c$$$C     ----------
c$$$      NTRY=NTRY+1
c$$$      ANS = 0D0
c$$$      DO IHEL=1,NCOMB
c$$$        IF (GOODHEL(IHEL) .OR. NTRY .LT. 2) THEN
c$$$          IF (NTRY.LT.2) THEN
c$$$C           for the first ps-point, check for helicities that give
c$$$C           identical matrix elements
c$$$            T=MATRIX_2(P ,NHEL(1,IHEL))
c$$$            T_SAVE(IHEL)=T
c$$$            T_IDENT(IHEL)=-1
c$$$            DO I=1,IHEL-1
c$$$              IF (T.EQ.0D0) EXIT
c$$$              IF (T_SAVE(I).EQ.0D0) CYCLE
c$$$              IF (ABS(T/T_SAVE(I)-1D0) .LT. 1D-12) THEN
c$$$C               WRITE (*,*) 'FOUND IDENTICAL',T,IHEL,T_SAVE(I),I
c$$$                T_IDENT(IHEL) = I
c$$$              ENDIF
c$$$            ENDDO
c$$$          ELSE
c$$$            IF (T_IDENT(IHEL).GT.0) THEN
c$$$C             if two helicity states are identical, dont recompute
c$$$              T=T_SAVE(T_IDENT(IHEL))
c$$$              T_SAVE(IHEL)=T
c$$$            ELSE
c$$$              T=MATRIX_2(P ,NHEL(1,IHEL))
c$$$              T_SAVE(IHEL)=T
c$$$            ENDIF
c$$$          ENDIF
c$$$C         add to the sum of helicities
c$$$          ANS=ANS+T
c$$$          IF (T .NE. 0D0 .AND. .NOT. GOODHEL(IHEL)) THEN
c$$$            GOODHEL(IHEL)=.TRUE.
c$$$          ENDIF
c$$$        ENDIF
c$$$      ENDDO
c$$$      ANS=ANS/DBLE(IDEN) * AMPLOtotop()

      sman=(P(0,1)+P(0,2))**2-(P(1,1)+P(1,2))**2
     & -(P(2,1)+P(2,2))**2 -(P(3,1)+P(3,2))**2

c     CAVEAT to avoid rounding errors causing t to become 0d0.
c      t=(P(0,1)-P(0,3))**2-(P(1,1)-P(1,3))**2
c     #-(P(2,1)-P(2,3))**2 -(P(3,1)-P(3,3))**2
      upt1 = (P(0,1)-P(0,3))**2-(P(3,1)-P(3,3))**2
      upt2 = -(P(1,1)-P(1,3))**2-(P(2,1)-P(2,3))**2
      uman = upt1+upt2

c     CAVEAT to avoid rounding errors causing u/t to become 0d0.
c      u=(P(0,1)-P(0,4))**2-(P(1,1)-P(1,4))**2
c     #-(P(2,1)-P(2,4))**2 -(P(3,1)-P(3,4))**2
      tpt1 = (P(0,1)-P(0,4))**2-(P(3,1)-P(3,4))**2
      tpt2 = -(P(2,1)-P(2,4))**2 -(P(1,1)-P(1,4))**2
      tman = tpt1+tpt2

      mh2 = sman + tman + uman

      ANS = (g**2/4.d0/pi)**3*AMPqg(sman,tman,uman)
     & *dsqrt(2.d0)/12.d0/Pi*mh2*AMPLO(mh2) * 2d0/9d0*MDL_GF

c$$$      sman =   2d0 * dot(p(0,1), p(0,2))
c$$$      tman = - 2d0 * dot(p(0,1), p(0,4))
c$$$      uman = - 2d0 * dot(p(0,2), p(0,4))
c$$$
c$$$      ANS = AMPqg(sman,tman,uman)
c$$$      ANS = ANS * (g**2/(8d0*pi**2))**3*(32d0/3d0)/48d0
c$$$     f /0.60365910530599576d0 * AMPLOtotop()
c$$$!      print*, 'ours  : ',ANS,AMPqg(sman,tman,uman),g**2/(8d0*pi**2) 
c$$$!      print*, 'ratio : ',ANS/dum,sman,uman,tman
c$$$!      print*, 'born  : ',born_wgt_ao2pi
c$$$!      print*, '---------------------------------------------'
c$$$
c$$$!      ANS = dum
      WGT_ME_REAL=ANS

      END


      REAL*8 FUNCTION MATRIX_2(P,NHEL)
C     
C     Generated by MadGraph5_aMC@NLO v. 2.2.3, 2015-02-10
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: d~ g > x0 d~ XGLU=1 WEIGHTED=3 QNP=1 [ real = QCD ] / t
C     Process: u~ g > x0 u~ XGLU=1 WEIGHTED=3 QNP=1 [ real = QCD ] / t
C     Process: s~ g > x0 s~ XGLU=1 WEIGHTED=3 QNP=1 [ real = QCD ] / t
C     Process: c~ g > x0 c~ XGLU=1 WEIGHTED=3 QNP=1 [ real = QCD ] / t
C     Process: b~ g > x0 b~ XGLU=1 WEIGHTED=3 QNP=1 [ real = QCD ] / t
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=1)
      INTEGER    NWAVEFUNCS, NCOLOR
      PARAMETER (NWAVEFUNCS=5, NCOLOR=1)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(0D0,1D0))
      INCLUDE 'nexternal.inc'
      INCLUDE 'coupl.inc'
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J
      INTEGER IC(NEXTERNAL)
      DATA IC /NEXTERNAL*1/
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 ZTEMP, AMP(NGRAPHS), JAMP(NCOLOR), W(8,NWAVEFUNCS)
C     
C     COLOR DATA
C     
      DATA DENOM(1)/1/
      DATA (CF(I,  1),I=  1,  1) /    4/
C     1 T(2,1,4)
C     ----------
C     BEGIN CODE
C     ----------
      CALL OXXXXX(P(0,1),ZERO,NHEL(1),-1*IC(1),W(1,1))
      CALL VXXXXX(P(0,2),ZERO,NHEL(2),-1*IC(2),W(1,2))
      CALL SXXXXX(P(0,3),+1*IC(3),W(1,3))
      CALL IXXXXX(P(0,4),ZERO,NHEL(4),-1*IC(4),W(1,4))
      CALL FFV1P0_3(W(1,4),W(1,1),GC_5,ZERO,ZERO,W(1,5))
C     Amplitude(s) for diagram number 1
      CALL VVS2_10_0(W(1,2),W(1,5),W(1,3),GC_3002H,GC_3002A,AMP(1))
      JAMP(1)=-AMP(1)
      MATRIX_2 = 0.D0
      DO I = 1, NCOLOR
        ZTEMP = (0.D0,0.D0)
        DO J = 1, NCOLOR
          ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
        ENDDO
        MATRIX_2 = MATRIX_2+ZTEMP*DCONJG(JAMP(I))/DENOM(I)
      ENDDO
      END

