C This is zerox64.F of CERN libraries, stripped of pre-compiler
C instructions by SF, 26/7/2014.
C Error messages are also suppressed, and replaced by a return
C integer (IERR), which is set to zero when no errors are encountered. 
C When IERR#0, the function assumes a value equal to -1.d20. 
C These manipulations implied that lines #124 and 125
C in this file have been interchanged
*
* $Id: zerox64.F,v 1.1.1.1 1996/04/01 15:01:51 mclareni Exp $
*
* $Log: zerox64.F,v $
* Revision 1.1.1.1  1996/04/01 15:01:51  mclareni
* Mathlib gen
*
*
      FUNCTION DZEROX(A0,B0,EPS,MAXF,F,MODE,IERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     Based on
C
C        J.C.P. Bus and T.J. Dekker, Two Efficient Algorithms with
C        Guaranteed Convergence for Finding a Zero of a Function,
C        ACM Trans. Math. Software 1 (1975) 330-345.
C
C        (MODE = 1: Algorithm M;    MODE = 2: Algorithm R)
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'DZEROX')
      LOGICAL LMT

      DIMENSION IM1(2),IM2(2),LMT(2)

      PARAMETER (Z1 = 1, HALF = Z1/2)

      DATA IM1 /2,3/, IM2 /-1,3/

      IERR=0
      IF(MODE .NE. 1 .AND. MODE .NE. 2) THEN
        C=-1.D20
        IERR=1
c$$$       C=0
c$$$       WRITE(ERRTXT,101) MODE
c$$$       CALL MTLPRT(NAME,'C200.1',ERRTXT)
       GO TO 99
      ENDIF
      FA=F(B0)
      FB=F(A0)
      IF(FA*FB .GT. 0) THEN
       C=-1.D20
       IERR=2
c$$$       C=0
c$$$       WRITE(ERRTXT,102) A0,B0
c$$$       CALL MTLPRT(NAME,'C200.2',ERRTXT)
       GO TO 99
      ENDIF
      ATL=ABS(EPS)
      B=A0
      A=B0
      LMT(2)=.TRUE.
      MF=2
    1 C=A
      FC=FA
    2 IE=0
    3 IF(ABS(FC) .LT. ABS(FB)) THEN
       IF(C .NE. A) THEN
        D=A
        FD=FA
       END IF
       A=B
       B=C
       C=A
       FA=FB
       FB=FC
       FC=FA
      END IF
      TOL=ATL*(1+ABS(C))
      H=HALF*(C+B)
      HB=H-B
      IF(ABS(HB) .GT. TOL) THEN
       IF(IE .GT. IM1(MODE)) THEN
        W=HB
       ELSE
        TOL=TOL*SIGN(Z1,HB)
        P=(B-A)*FB
        LMT(1)=IE .LE. 1
        IF(LMT(MODE)) THEN
         Q=FA-FB
         LMT(2)=.FALSE.
        ELSE
         FDB=(FD-FB)/(D-B)
         FDA=(FD-FA)/(D-A)
         P=FDA*P
         Q=FDB*FA-FDA*FB
        END IF
        IF(P .LT. 0) THEN
         P=-P
         Q=-Q
        END IF
        IF(IE .EQ. IM2(MODE)) P=P+P
        IF(P .EQ. 0 .OR. P .LE. Q*TOL) THEN
         W=TOL
        ELSEIF(P .LT. HB*Q) THEN
         W=P/Q
        ELSE
         W=HB
        END IF
       END IF
       D=A
       A=B
       FD=FA
       FA=FB
       B=B+W
       MF=MF+1
       IF(MF .GT. MAXF) THEN
        C=-1.D20
        IERR=3
c$$$        CALL MTLPRT(NAME,'C200.3','TOO MANY FUNCTION CALLS')
        GO TO 99
       ENDIF
       FB=F(B)
       IF(FB .EQ. 0 .OR. SIGN(Z1,FC) .EQ. SIGN(Z1,FB)) GO TO 1
       IF(W .EQ. HB) GO TO 2
       IE=IE+1
       GO TO 3
      END IF
   99 CONTINUE
      DZEROX=C
      RETURN
  101 FORMAT('MODE = ',I3,' ILLEGAL')
  102 FORMAT('F(A) AND F(B) HAVE THE SAME SIGN, A = ',1P,D15.8,
     1       ', B = ',D15.8)
      END


*
* $Id: mtlprt.F,v 1.1.1.1 1996/04/01 15:02:52 mclareni Exp $
*
* $Log: mtlprt.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:52  mclareni
* Mathlib gen
*
*
c$$$#include "gen/pilot.h"
      SUBROUTINE MTLPRT(NAME,ERC,TEXT)
      CHARACTER*(*) NAME,ERC,TEXT
      LOGICAL LMF,LRF

      IF(ERC(5:6).NE.'.0') THEN
        CALL MTLMTR(ERC,MLG,LMF,LRF)
      ELSE
        LMF=.TRUE.
        LRF=.FALSE.
      ENDIF
      IF(LMF) THEN
        LT=LENOCC(TEXT)
        IF(MLG .LT. 1) WRITE(  *,100) ERC(1:4),NAME,ERC,TEXT(1:LT)
        IF(MLG .GE. 1) WRITE(MLG,100) ERC(1:4),NAME,ERC,TEXT(1:LT)
      ENDIF
      IF(.NOT.LRF) CALL ABEND
      RETURN
100   FORMAT(7X,'***** CERN ',A,1X,A,' ERROR ',A,': ',A)
      END


*
* $Id: mtlset.F,v 1.1.1.1 1996/04/01 15:02:53 mclareni Exp $
*
* $Log: mtlset.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:53  mclareni
* Mathlib gen
*
*
c$$$#include "gen/pilot.h"
      SUBROUTINE MTLSET(ERC,NLG,MXM,MXR)

      PARAMETER (KTE = 132)
      CHARACTER*6 ERC,CODE(KTE)
      LOGICAL LMF,LRF
      DIMENSION KNTM(KTE),KNTR(KTE)

      DATA ILG /0/

C     renumber the data statements after putting new codes in Unix with:
C     awk -F'[()]' '{ printf"%s(%s)%s(%s)%s(%s)%s\n",$1,NR,$3,NR,$5,NR,$7 }'
C     and modify KTE to the number of lines below

      DATA CODE(1),KNTM(1),KNTR(1) / 'B100.1', 255, 255 /
      DATA CODE(2),KNTM(2),KNTR(2) / 'B300.1', 255, 255 /
      DATA CODE(3),KNTM(3),KNTR(3) / 'B300.2', 255, 255 /
      DATA CODE(4),KNTM(4),KNTR(4) / 'C200.0', 255, 255 /
      DATA CODE(5),KNTM(5),KNTR(5) / 'C200.1', 255, 255 /
      DATA CODE(6),KNTM(6),KNTR(6) / 'C200.2', 255, 255 /
      DATA CODE(7),KNTM(7),KNTR(7) / 'C200.3', 255, 255 /
      DATA CODE(8),KNTM(8),KNTR(8) / 'C201.0', 255, 255 /
      DATA CODE(9),KNTM(9),KNTR(9) / 'C202.0', 255, 255 /
      DATA CODE(10),KNTM(10),KNTR(10) / 'C202.1', 255, 255 /
      DATA CODE(11),KNTM(11),KNTR(11) / 'C202.2', 255, 255 /
      DATA CODE(12),KNTM(12),KNTR(12) / 'C205.1', 255, 255 /
      DATA CODE(13),KNTM(13),KNTR(13) / 'C205.2', 255, 255 /
      DATA CODE(14),KNTM(14),KNTR(14) / 'C207.0', 255, 255 /
      DATA CODE(15),KNTM(15),KNTR(15) / 'C208.0', 255, 255 /
      DATA CODE(16),KNTM(16),KNTR(16) / 'C209.0', 255, 255 /
      DATA CODE(17),KNTM(17),KNTR(17) / 'C209.1', 255, 255 /
      DATA CODE(18),KNTM(18),KNTR(18) / 'C209.2', 255, 255 /
      DATA CODE(19),KNTM(19),KNTR(19) / 'C209.3', 255, 255 /
      DATA CODE(20),KNTM(20),KNTR(20) / 'C210.1', 255, 255 /
      DATA CODE(21),KNTM(21),KNTR(21) / 'C302.1', 255, 255 /
      DATA CODE(22),KNTM(22),KNTR(22) / 'C303.1', 255, 255 /
      DATA CODE(23),KNTM(23),KNTR(23) / 'C304.1', 255, 255 /
      DATA CODE(24),KNTM(24),KNTR(24) / 'C305.1', 255, 255 /
      DATA CODE(25),KNTM(25),KNTR(25) / 'C306.1', 255, 255 /
      DATA CODE(26),KNTM(26),KNTR(26) / 'C307.1', 255, 255 /
      DATA CODE(27),KNTM(27),KNTR(27) / 'C312.1', 255, 255 /
      DATA CODE(28),KNTM(28),KNTR(28) / 'C313.1', 255, 255 /
      DATA CODE(29),KNTM(29),KNTR(29) / 'C315.1', 255, 255 /
      DATA CODE(30),KNTM(30),KNTR(30) / 'C316.1', 255, 255 /
      DATA CODE(31),KNTM(31),KNTR(31) / 'C316.2', 255, 255 /
      DATA CODE(32),KNTM(32),KNTR(32) / 'C320.1', 255, 255 /
      DATA CODE(33),KNTM(33),KNTR(33) / 'C321.1', 255, 255 /
      DATA CODE(34),KNTM(34),KNTR(34) / 'C323.1', 255, 255 /
      DATA CODE(35),KNTM(35),KNTR(35) / 'C327.1', 255, 255 /
      DATA CODE(36),KNTM(36),KNTR(36) / 'C328.1', 255, 255 /
      DATA CODE(37),KNTM(37),KNTR(37) / 'C328.2', 255, 255 /
      DATA CODE(38),KNTM(38),KNTR(38) / 'C328.3', 255, 255 /
      DATA CODE(39),KNTM(39),KNTR(39) / 'C330.1', 255, 255 /
      DATA CODE(40),KNTM(40),KNTR(40) / 'C330.2', 255, 255 /
      DATA CODE(41),KNTM(41),KNTR(41) / 'C330.3', 255, 255 /
      DATA CODE(42),KNTM(42),KNTR(42) / 'C331.1', 255, 255 /
      DATA CODE(43),KNTM(43),KNTR(43) / 'C331.2', 255, 255 /
      DATA CODE(44),KNTM(44),KNTR(44) / 'C334.1', 255, 255 /
      DATA CODE(45),KNTM(45),KNTR(45) / 'C334.2', 255, 255 /
      DATA CODE(46),KNTM(46),KNTR(46) / 'C334.3', 255, 255 /
      DATA CODE(47),KNTM(47),KNTR(47) / 'C334.4', 255, 255 /
      DATA CODE(48),KNTM(48),KNTR(48) / 'C334.5', 255, 255 /
      DATA CODE(49),KNTM(49),KNTR(49) / 'C334.6', 255, 255 /
      DATA CODE(50),KNTM(50),KNTR(50) / 'C336.1', 255, 255 /
      DATA CODE(51),KNTM(51),KNTR(51) / 'C337.1', 255, 255 /
      DATA CODE(52),KNTM(52),KNTR(52) / 'C338.1', 255, 255 /
      DATA CODE(53),KNTM(53),KNTR(53) / 'C340.1', 255, 255 /
      DATA CODE(54),KNTM(54),KNTR(54) / 'C343.1', 255, 255 /
      DATA CODE(55),KNTM(55),KNTR(55) / 'C343.2', 255, 255 /
      DATA CODE(56),KNTM(56),KNTR(56) / 'C343.3', 255, 255 /
      DATA CODE(57),KNTM(57),KNTR(57) / 'C343.4', 255, 255 /
      DATA CODE(58),KNTM(58),KNTR(58) / 'C344.1', 255, 255 /
      DATA CODE(59),KNTM(59),KNTR(59) / 'C344.2', 255, 255 /
      DATA CODE(60),KNTM(60),KNTR(60) / 'C344.3', 255, 255 /
      DATA CODE(61),KNTM(61),KNTR(61) / 'C344.4', 255, 255 /
      DATA CODE(62),KNTM(62),KNTR(62) / 'C345.1', 255, 255 /
      DATA CODE(63),KNTM(63),KNTR(63) / 'C346.1', 255, 255 /
      DATA CODE(64),KNTM(64),KNTR(64) / 'C346.2', 255, 255 /
      DATA CODE(65),KNTM(65),KNTR(65) / 'C346.3', 255, 255 /
      DATA CODE(66),KNTM(66),KNTR(66) / 'C347.1', 255, 255 /
      DATA CODE(67),KNTM(67),KNTR(67) / 'C347.2', 255, 255 /
      DATA CODE(68),KNTM(68),KNTR(68) / 'C347.3', 255, 255 /
      DATA CODE(69),KNTM(69),KNTR(69) / 'C347.4', 255, 255 /
      DATA CODE(70),KNTM(70),KNTR(70) / 'C347.5', 255, 255 /
      DATA CODE(71),KNTM(71),KNTR(71) / 'C347.6', 255, 255 /
      DATA CODE(72),KNTM(72),KNTR(72) / 'C348.1', 255, 255 /
      DATA CODE(73),KNTM(73),KNTR(73) / 'C349.1', 255, 255 /
      DATA CODE(74),KNTM(74),KNTR(74) / 'C349.2', 255, 255 /
      DATA CODE(75),KNTM(75),KNTR(75) / 'C349.3', 255, 255 /
      DATA CODE(76),KNTM(76),KNTR(76) / 'D101.1', 255, 255 /
      DATA CODE(77),KNTM(77),KNTR(77) / 'D103.1', 255, 255 /
      DATA CODE(78),KNTM(78),KNTR(78) / 'D104.1', 255, 255 /
      DATA CODE(79),KNTM(79),KNTR(79) / 'D104.2', 255, 255 /
      DATA CODE(80),KNTM(80),KNTR(80) / 'D105.1', 255, 255 /
      DATA CODE(81),KNTM(81),KNTR(81) / 'D105.2', 255, 255 /
      DATA CODE(82),KNTM(82),KNTR(82) / 'D107.1', 255, 255 /
      DATA CODE(83),KNTM(83),KNTR(83) / 'D110.0', 255, 255 /
      DATA CODE(84),KNTM(84),KNTR(84) / 'D110.1', 255, 255 /
      DATA CODE(85),KNTM(85),KNTR(85) / 'D110.2', 255, 255 /
      DATA CODE(86),KNTM(86),KNTR(86) / 'D110.3', 255, 255 /
      DATA CODE(87),KNTM(87),KNTR(87) / 'D110.4', 255, 255 /
      DATA CODE(88),KNTM(88),KNTR(88) / 'D110.5', 255, 255 /
      DATA CODE(89),KNTM(89),KNTR(89) / 'D110.6', 255, 255 /
      DATA CODE(90),KNTM(90),KNTR(90) / 'D113.1', 255, 255 /
      DATA CODE(91),KNTM(91),KNTR(91) / 'D201.1', 255, 255 /
      DATA CODE(92),KNTM(92),KNTR(92) / 'D202.1', 255, 255 /
      DATA CODE(93),KNTM(93),KNTR(93) / 'D401.1', 255, 255 /
      DATA CODE(94),KNTM(94),KNTR(94) / 'D601.1', 255, 255 /
      DATA CODE(95),KNTM(95),KNTR(95) / 'E210.1', 255, 255 /
      DATA CODE(96),KNTM(96),KNTR(96) / 'E210.2', 255, 255 /
      DATA CODE(97),KNTM(97),KNTR(97) / 'E210.3', 255, 255 /
      DATA CODE(98),KNTM(98),KNTR(98) / 'E210.4', 255, 255 /
      DATA CODE(99),KNTM(99),KNTR(99) / 'E210.5', 255, 255 /
      DATA CODE(100),KNTM(100),KNTR(100) / 'E210.6', 255, 255 /
      DATA CODE(101),KNTM(101),KNTR(101) / 'E210.7', 255, 255 /
      DATA CODE(102),KNTM(102),KNTR(102) / 'E211.0', 255, 255 /
      DATA CODE(103),KNTM(103),KNTR(103) / 'E211.1', 255, 255 /
      DATA CODE(104),KNTM(104),KNTR(104) / 'E211.2', 255, 255 /
      DATA CODE(105),KNTM(105),KNTR(105) / 'E211.3', 255, 255 /
      DATA CODE(106),KNTM(106),KNTR(106) / 'E211.4', 255, 255 /
      DATA CODE(107),KNTM(107),KNTR(107) / 'E406.0', 255, 255 /
      DATA CODE(108),KNTM(108),KNTR(108) / 'E406.1', 255, 255 /
      DATA CODE(109),KNTM(109),KNTR(109) / 'E407.0', 255, 255 /
      DATA CODE(110),KNTM(110),KNTR(110) / 'E408.0', 255, 255 /
      DATA CODE(111),KNTM(111),KNTR(111) / 'E408.1', 255, 255 /
      DATA CODE(112),KNTM(112),KNTR(112) / 'F500.0', 255, 255 /
      DATA CODE(113),KNTM(113),KNTR(113) / 'F500.1', 255, 255 /
      DATA CODE(114),KNTM(114),KNTR(114) / 'F500.2', 255, 255 /
      DATA CODE(115),KNTM(115),KNTR(115) / 'F500.3', 255, 255 /
      DATA CODE(116),KNTM(116),KNTR(116) / 'G100.1', 255, 255 /
      DATA CODE(117),KNTM(117),KNTR(117) / 'G100.2', 255, 255 /
      DATA CODE(118),KNTM(118),KNTR(118) / 'G101.1', 255, 255 /
      DATA CODE(119),KNTM(119),KNTR(119) / 'G101.2', 255, 255 /
      DATA CODE(120),KNTM(120),KNTR(120) / 'G105.1', 255, 255 /
      DATA CODE(121),KNTM(121),KNTR(121) / 'G106.1', 255, 255 /
      DATA CODE(122),KNTM(122),KNTR(122) / 'G106.2', 255, 255 /
      DATA CODE(123),KNTM(123),KNTR(123) / 'G116.1', 255, 255 /
      DATA CODE(124),KNTM(124),KNTR(124) / 'G116.2', 255, 255 /
      DATA CODE(125),KNTM(125),KNTR(125) / 'H101.0', 255, 255 /
      DATA CODE(126),KNTM(126),KNTR(126) / 'H101.1', 255, 255 /
      DATA CODE(127),KNTM(127),KNTR(127) / 'H101.2', 255, 255 /
      DATA CODE(128),KNTM(128),KNTR(128) / 'H301.1', 255, 255 /
      DATA CODE(129),KNTM(129),KNTR(129) / 'U501.1', 255, 255 /
      DATA CODE(130),KNTM(130),KNTR(130) / 'V202.1', 255, 255 /
      DATA CODE(131),KNTM(131),KNTR(131) / 'V202.2', 255, 255 /
      DATA CODE(132),KNTM(132),KNTR(132) / 'V202.3', 255, 255 /

c$$$#if defined(CERNLIB_IBMVM)||defined(CERNLIB_IBMMVS)
c$$$      EXTERNAL C302ERR,C304ERR
c$$$
c$$$      LOGICAL LERSET
c$$$      DATA LERSET/.FALSE./
c$$$      SAVE LERSET
c$$$
c$$$      IF(.NOT.LERSET) THEN
c$$$        LERSET=.TRUE.
c$$$C---  AFB290I: Argument for GAMMA  is out of range
c$$$        CALL ERRSET(290,256,-1,1,C302ERR)
c$$$C---  AFB291I: Argument for ALGAMA is out of range
c$$$        CALL ERRSET(291,256,-1,1,C304ERR)
c$$$C---  AFB300I: Argument for DGAMMA is out of range
c$$$        CALL ERRSET(300,256,-1,1,C302ERR)
c$$$C---  AFB301I: Argument for DLGAMA is out of range
c$$$        CALL ERRSET(301,256,-1,1,C304ERR)
c$$$      ENDIF
c$$$#endif

      ILG=NLG
      L=0
      IF(ERC .NE. ' ') THEN
       DO 10 L = 1,6
       IF(ERC(1:L) .EQ. ERC) GOTO 12
   10  CONTINUE
   12  CONTINUE
      ENDIF
      DO 14 I = 1,KTE
      IF(L .EQ. 0 .OR. CODE(I)(1:L) .EQ. ERC(1:L)) THEN
       IF(MXM .GE. 0) KNTM(I)=MXM
       IF(MXR .GE. 0) KNTR(I)=MXR
      ENDIF
   14 CONTINUE
      RETURN

      ENTRY MTLMTR(ERC,MLG,LMF,LRF)

      MLG=ILG
      DO 20 I = 1,KTE
      IF(ERC .EQ. CODE(I))  GOTO 21
   20 CONTINUE
      WRITE(*,100) ERC
      CALL ABEND
      RETURN

   21 LMF=KNTM(I) .GE. 1
      LRF=KNTR(I) .GE. 1
      IF(LMF .AND. KNTM(I) .LT. 255)  KNTM(I)=KNTM(I)-1
      IF(LRF .AND. KNTR(I) .LT. 255)  KNTR(I)=KNTR(I)-1
      IF(.NOT.LRF) THEN
       IF(ILG .LT. 1) WRITE(  *,101) CODE(I)
       IF(ILG .GE. 1) WRITE(ILG,101) CODE(I)
      ENDIF
      RETURN
  100 FORMAT(7X,'***** CERN N002 MTLSET ... ERROR N002: ',
     1'ERROR CODE ',A6,' NOT RECOGNIZED BY ERROR MONITOR. RUN ABORTED.')
  101 FORMAT(7X,'***** CERN N002 MTLSET ... ERROR NOO2.1: ',
     1'RUN TERMINATED BY LIBRARY ERROR CONDITION ',A6)
      END


*
* $Id: lenocc.F,v 1.1.1.1 1996/02/15 17:49:49 mclareni Exp $
*
* $Log: lenocc.F,v $
* Revision 1.1.1.1  1996/02/15 17:49:49  mclareni
* Kernlib
*
*
c$$$#include "kerngen/pilot.h"
      FUNCTION LENOCC (CHV)
C
C CERN PROGLIB# M507    LENOCC          .VERSION KERNFOR  4.21  890323
C ORIG. March 85, A.Petrilli, re-write 21/02/89, JZ
C
C-    Find last non-blank character in CHV

      CHARACTER    CHV*(*)

      N = LEN(CHV)

      DO 17  JJ= N,1,-1
      IF (CHV(JJ:JJ).NE.' ') GO TO 99
   17 CONTINUE
      JJ = 0

   99 LENOCC = JJ
      RETURN
      END


*
* $Id: abend.F,v 1.1.1.1 1996/02/15 17:48:36 mclareni Exp $
*
* $Log: abend.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:36  mclareni
* Kernlib
*
*
c$$$#include "kernnumt/pilot.h"
          SUBROUTINE ABEND
c$$$#include "kernnumt/sysdat.inc"
          IF(LGFILE .EQ. 0)  WRITE(*,1000)
          IF(LGFILE .NE. 0)  WRITE(LGFILE,1000)
          RETURN
1000      FORMAT(31H ABEND ROUTINE HAS BEEN CALLED.)
          END
