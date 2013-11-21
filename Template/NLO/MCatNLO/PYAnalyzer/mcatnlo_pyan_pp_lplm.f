c
c Example analysis for "p p > l+ l- [QCD]" process.
c
C----------------------------------------------------------------------
      SUBROUTINE RCLOS()
C     DUMMY IF HBOOK IS USED
C----------------------------------------------------------------------
      END


C----------------------------------------------------------------------
      SUBROUTINE PYABEG
C     USER'S ROUTINE FOR INITIALIZATION
C----------------------------------------------------------------------
      implicit none
      include 'reweight0.inc'
      real * 8 bin,xmi,xms,pi
      PARAMETER (PI=3.14159265358979312D0)
      integer j,kk,l
      character*5 cc(2)
      data cc/'     ',' cuts'/
      integer nwgt,max_weight,nwgt_analysis
      common/cnwgt/nwgt
      common/c_analysis/nwgt_analysis
      parameter (max_weight=maxscales*maxscales+maxpdfs+1)
      character*15 weights_info(max_weight)
      common/cwgtsinfo/weights_info
      call inihist
      nwgt_analysis=nwgt
c
      xmi=50.d0
      xms=130.d0
      bin=0.8d0
      do kk=1,nwgt_analysis
      do j=1,2
      l=(kk-1)*42+(j-1)*21
      call mbook(l+ 1,'V pt      '//weights_info(kk)//cc(j)
     &     ,2.d0,0.d0,200.d0)
      call mbook(l+ 2,'V pt      '//weights_info(kk)//cc(j)
     &     ,10.d0,0.d0,1000.d0)
      call mbook(l+ 3,'V log[pt] '//weights_info(kk)//cc(j)
     &     ,0.05d0,0.1d0,5.d0)
      call mbook(l+ 4,'V y       '//weights_info(kk)//cc(j)
     &     ,0.2d0,-9.d0,9.d0)
      call mbook(l+ 5,'V eta     '//weights_info(kk)//cc(j)
     &     ,0.2d0,-9.d0,9.d0)
      call mbook(l+ 6,'mV        '//weights_info(kk)//cc(j)
     &     ,bin,xmi,xms)
c
      call mbook(l+ 7,'lm pt      '//weights_info(kk)//cc(j)
     &     ,2.d0,0.d0,200.d0)
      call mbook(l+ 8,'lm pt      '//weights_info(kk)//cc(j)
     &     ,10.d0,0.d0,1000.d0)
      call mbook(l+ 9,'lm log[pt] '//weights_info(kk)//cc(j)
     &     ,0.05d0,0.1d0,5.d0)
      call mbook(l+10,'lm eta     '//weights_info(kk)//cc(j)
     &     ,0.2d0,-9.d0,9.d0)
      call mbook(l+11,'lp pt      '//weights_info(kk)//cc(j)
     &     ,2.d0,0.d0,200.d0)
      call mbook(l+12,'lp pt      '//weights_info(kk)//cc(j)
     &     ,10.d0,0.d0,1000.d0)
      call mbook(l+13,'lp log[pt] '//weights_info(kk)//cc(j)
     &     ,0.05d0,0.1d0,5.d0)
      call mbook(l+14,'lp eta     '//weights_info(kk)//cc(j)
     &     ,0.2d0,-9.d0,9.d0)
c
      call mbook(l+15,'lmlp delta eta     '//weights_info(kk)//cc(j)
     $     ,0.2d0,-9.d0,9.d0)
      call mbook(l+16,'lmlp azimt         '//weights_info(kk)//cc(j)
     $     ,pi/20.d0,0.d0,pi)
      call mbook(l+17,'lmlp log[pi-azimt] '//weights_info(kk)//cc(j)
     $     ,0.05d0,-4.d0,0.1d0)
      call mbook(l+18,'lmlp inv m         '//weights_info(kk)//cc(j)
     $     ,bin,xmi,xms)
      call mbook(l+19,'lmlp pt            '//weights_info(kk)//cc(j)
     $     ,2.d0,0.d0,200.d0)
      call mbook(l+20,'lmlp log[pt]       '//weights_info(kk)//cc(j)
     $     ,0.05d0,0.1d0,5.d0)
c
      call mbook(l+21,'total'//weights_info(kk)//cc(j),1.d0,-1.d0,1.d0)
      enddo
      enddo
 999  END


C----------------------------------------------------------------------
      SUBROUTINE PYAEND(IEVT)
C     USER'S ROUTINE FOR TERMINAL CALCULATIONS, HISTOGRAM OUTPUT, ETC
C----------------------------------------------------------------------
      REAL*8 XNORM
      INTEGER I,J,KK,l,nwgt_analysis
      integer NPL
      parameter(NPL=15000)
      common/c_analysis/nwgt_analysis
      OPEN(UNIT=99,FILE='PYTLL.TOP',STATUS='UNKNOWN')
      XNORM=1.D0/IEVT
      DO I=1,NPL
 	CALL MFINAL3(I)             
        CALL MCOPY(I,I+NPL)
        CALL MOPERA(I+NPL,'F',I+NPL,I+NPL,(XNORM),0.D0)
 	CALL MFINAL3(I+NPL)             
      ENDDO                          
C
      do kk=1,nwgt_analysis
      do i=1,2
      l=(kk-1)*42+(i-1)*21
C
      call multitop(NPL+l+ 1,NPL-1,3,2,'V pt',' ','LOG')
      call multitop(NPL+l+ 2,NPL-1,3,2,'V pt',' ','LOG')
      call multitop(NPL+l+ 3,NPL-1,3,2,'V log[pt]',' ','LOG')
      call multitop(NPL+l+ 4,NPL-1,3,2,'V y',' ','LOG')
      call multitop(NPL+l+ 5,NPL-1,3,2,'V eta',' ','LOG')
      call multitop(NPL+l+ 6,NPL-1,3,2,'mV',' ','LOG')
c
      call multitop(NPL+l+ 7,NPL-1,3,2,'lm pt',' ','LOG')
      call multitop(NPL+l+ 8,NPL-1,3,2,'lm pt',' ','LOG')
      call multitop(NPL+l+ 9,NPL-1,3,2,'lm log[pt]',' ','LOG')
      call multitop(NPL+l+10,NPL-1,3,2,'lm eta',' ','LOG')
      call multitop(NPL+l+11,NPL-1,3,2,'lm pt',' ','LOG')
      call multitop(NPL+l+12,NPL-1,3,2,'lm pt',' ','LOG')
      call multitop(NPL+l+13,NPL-1,3,2,'lm log[pt]',' ','LOG')
      call multitop(NPL+l+14,NPL-1,3,2,'lm eta',' ','LOG')
c
      call multitop(NPL+l+15,NPL-1,3,2,'lmlp deta',' ','LOG')
      call multitop(NPL+l+16,NPL-1,3,2,'lmlp azi',' ','LOG')
      call multitop(NPL+l+17,NPL-1,3,2,'lmlp azi',' ','LOG')
      call multitop(NPL+l+18,NPL-1,3,2,'lmlp inv m',' ','LOG')
      call multitop(NPL+l+19,NPL-1,3,2,'lmlp pt',' ','LOG')
      call multitop(NPL+l+20,NPL-1,3,2,'lmlp pt',' ','LOG')
c
      call multitop(NPL+l+21,NPL-1,3,2,'total',' ','LOG')
      enddo
      enddo
c
      CLOSE(99)
      END

C----------------------------------------------------------------------
      SUBROUTINE PYANAL
C     USER'S ROUTINE TO ANALYSE DATA FROM EVENT
C----------------------------------------------------------------------
      implicit double precision(a-h, o-z)
      implicit integer(i-n)
      include 'reweight0.inc'
      DOUBLE PRECISION HWVDOT,PSUM(4),PPV(5),YCUT,XMV,PTV,YV,THV,ETAV,
     #  PPL(5),PPLB(5),PTL,YL,THL,ETAL,PLL,ENL,PTLB,YLB,THLB,ETALB,
     #  PLLB,ENLB,PTPAIR,DLL,CLL,AZI,AZINORM,XMLL,DETALLB,PPV0(5),
     #  PPL0(5),PPLB0(5),P1(4),P2(4),PIHEP(4)
      INTEGER ICHSUM,ICHINI,IHEP,IV,IFV,IST,ID,IJ,ID1,JPR,IDENT,
     #  ILL,ILLB,IHRD
      integer pychge
      external pydata
      common/pyjets/n,npad,k(4000,5),p(4000,5),v(4000,5)
      common/pydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/pydat2/kchg(500,4),pmas(500,4),parf(2000),vckm(4,4)
      common/pydat3/mdcy(500,3),mdme(8000,2),brat(8000),kfdp(8000,5)
      common/pysubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      common/pypars/mstp(200),parp(200),msti(200),pari(200)
      LOGICAL DIDSOF,TEST1,TEST2,TEST3,TEST4,TEST5,TEST6,TEST7,flag,ISLP,ISLM
      REAL*8 PI,wmass,wgamma,bwcutoff,getinvm,getdelphi,getrapidity,
     &getpseudorap
      PARAMETER (PI=3.14159265358979312D0)
      REAL*8 WWW0,TINY
      INTEGER KK,i,l
      DATA TINY/.1D-5/
      DOUBLE PRECISION EVWEIGHT
      COMMON/CEVWEIGHT/EVWEIGHT
      INTEGER IFAIL
      COMMON/CIFAIL/IFAIL
      SAVE INOBOSON,INOLEPTON,INOLEPTONB
      integer nwgt_analysis,max_weight
      common/c_analysis/nwgt_analysis
      parameter (max_weight=maxscales*maxscales+maxpdfs+1)
      double precision ww(max_weight),www(max_weight)
      common/cww/ww
      IF(IFAIL.EQ.1)RETURN
      IF (WW(1).EQ.0D0) THEN
         WRITE(*,*)'WW(1) = 0. Stopping'
         STOP
      ENDIF
      IDENT=23
C INCOMING PARTONS MAY TRAVEL IN THE SAME DIRECTION: IT'S A POWER-SUPPRESSED
C EFFECT, SO THROW THE EVENT AWAY
      IF(SIGN(1.D0,P(3,3)).EQ.SIGN(1.D0,P(4,3)))THEN
         WRITE(*,*)'WARNING 111 IN PYANAL'
         GOTO 999
      ENDIF
      DO I=1,nwgt_analysis
         WWW(I)=EVWEIGHT*ww(i)/ww(1)
      ENDDO
      DO I=1,4
         P1(I)=P(1,I)
         P2(I)=P(2,I)
      ENDDO
      CALL VVSUM(4,P1,P2,PSUM)
      CALL VSCA(4,-1D0,PSUM,PSUM)
      ICHSUM=0
      KF1=K(1,2)
      KF2=K(2,2)
      ICHINI=PYCHGE(KF1)+PYCHGE(KF2)
      IFV =0
      IFL =0
      IFLB=0
      IV0 =10000
      IL0 =10000
      ILB0=10000
C
      DO 100 IHEP=1,N
        DO J=1,4
          PIHEP(J)=P(IHEP,J)
        ENDDO
        IST =K(IHEP,1)
        ID1 =K(IHEP,2)
        IORI=K(IHEP,3)
        IF (IST.LE.10) THEN
          CALL VVSUM(4,PIHEP,PSUM,PSUM)
          ICHSUM=ICHSUM+PYCHGE(ID1)
        ENDIF
        TEST1=IORI.EQ.0
        TEST2=ID1.EQ.IDENT
        TEST3=ID1.GT.0.AND.ABS(ID1).GE.11.AND.ABS(ID1).LE.16
        TEST4=ID1.LT.0.AND.ABS(ID1).GE.11.AND.ABS(ID1).LE.16
        ISLP=ID1.EQ. 11.OR.ID1.EQ. 13.OR.ID1.EQ. 15
        ISLM=ID1.EQ.-11.OR.ID1.EQ.-13.OR.ID1.EQ.-15
        IF(TEST1.AND.TEST2)IV0=IHEP
        TEST5=IORI.EQ.IV0
        IF(TEST5)THEN
           IF(TEST2)THEN
              IV1=IHEP
              IFV=IFV+1
           ELSEIF(TEST3)THEN
              IL0 =IHEP
           ELSEIF(TEST4)THEN
              ILB0=IHEP
           ENDIF
        ENDIF
        TEST6=IORI.EQ.IL0
        TEST7=IORI.EQ.ILB0
        IF(TEST6.AND.TEST3)THEN
           IL =IHEP
           IFL=IFL+1
        ENDIF
        IF(TEST7.AND.TEST4)THEN
           ILB=IHEP
           IFLB=IFLB+1
        ENDIF
        IF((IST.GE.120.AND.IST.LE.125).AND.ISLM)ILL0=IHEP
        IF((IST.GE.120.AND.IST.LE.125).AND.ISLP)ILLB0=IHEP
 100  CONTINUE
C
      DO IJ=1,5
        PPV(IJ) =P(IV1,IJ)
        PPL(IJ) =P(IL, IJ)
        PPLB(IJ)=P(ILB,IJ)
      ENDDO
C CHECK MULTIPLICITIES
      IF(IFV.EQ.0) THEN
         INOBOSON=INOBOSON+1
         WRITE(*,*)'WARNING IN PYANAL: NO WEAK BOSON ',INOBOSON
         GOTO 999
      ELSEIF(IFV.GT.1) THEN
         WRITE(*,*)'WARNING IN PYANAL: TOO MANY WEAK BOSONS ',IFV
         GOTO 999
      ENDIF
      IF(IFL.EQ.0) THEN
         INOLEPTON=INOLEPTON+1
         WRITE(*,*)'WARNING IN PYANAL: NO LEPTON ',INOLEPTON
         GOTO 999
      ELSEIF(IFL.GT.1) THEN
         WRITE(*,*)'WARNING IN PYANAL: TOO MANY LEPTONS ',IFL
         GOTO 999
      ENDIF
      IF(IFLB.EQ.0) THEN
         INOLEPTONB=INOLEPTONB+1
         WRITE(*,*)'WARNING IN PYANAL: NO ANTILEPTON ',INOLEPTONB
         GOTO 999
      ELSEIF(IFLB.GT.1) THEN
         WRITE(*,*)'WARNING IN PYANAL: TOO MANY ANTILEPTONS ',IFLB
         GOTO 999
      ENDIF
C CHECK THAT THE LEPTONS ARE FINAL-STATE LEPTONS
      IF( (K(IL,1).NE.1).OR.K(ILB,1).NE.1 ) THEN
         WRITE(*,*)'WARNING 505 IN PYANAL'
         GOTO 999
      ENDIF
C CHECK MOMENTUM AND CHARGE CONSERVATION
      IF (VDOT(3,PSUM,PSUM).GT.1.E-4*P(1,4)**2) THEN
         WRITE(*,*)'WARNING 112 IN PYANAL'
         GOTO 999
      ENDIF
      IF (ICHSUM.NE.ICHINI) THEN
         WRITE(*,*)'WARNING 113 IN PYANAL'
         GOTO 999
      ENDIF
C FILL THE HISTOS
      YCUT=2.5D0
C Variables of the vector boson
      xmv=getinvm(ppv(4),ppv(1),ppv(2),ppv(3))
      ptv=sqrt(ppv(1)**2+ppv(2)**2)
      yv=getrapidity(ppv(4),ppv(3))
      etav=getpseudorap(ppv(4),ppv(1),ppv(2),ppv(3))
C Variables of the leptons
      ptl=sqrt(ppl(1)**2+ppl(2)**2)
      yl=getrapidity(ppl(4),ppl(3))
      etal=getpseudorap(ppl(4),ppl(1),ppl(2),ppl(3))
c
      ptlb=sqrt(pplb(1)**2+pplb(2)**2)
      ylb=getrapidity(pplb(4),pplb(3))
      etalb=getpseudorap(pplb(4),pplb(1),pplb(2),pplb(3))
c
      ptpair=ptv
      azi=getdelphi(ppl(1),pplb(1),ppl(2),pplb(2))
      azinorm=(pi-azi)/pi
      xmll=xmv
      detallb=etal-etalb
c
      do kk=1,nwgt_analysis
      l=(kk-1)*42
      call mfill(l+1,(ptv),(WWW(kk)))
      call mfill(l+2,(ptv),(WWW(kk)))
      if(ptv.gt.0.d0)call mfill(l+3,(log10(ptv)),(WWW(kk)))
      call mfill(l+4,(yv),(WWW(kk)))
      call mfill(l+5,(etav),(WWW(kk)))
      call mfill(l+6,(xmv),(WWW(kk)))
c
      call mfill(l+7,(ptl),(WWW(kk)))
      call mfill(l+8,(ptl),(WWW(kk)))
      if(ptl.gt.0.d0)call mfill(l+9,(log10(ptl)),(WWW(kk)))
      call mfill(l+10,(etal),(WWW(kk)))
      call mfill(l+11,(ptlb),(WWW(kk)))
      call mfill(l+12,(ptlb),(WWW(kk)))
      if(ptlb.gt.0.d0)call mfill(l+13,(log10(ptlb)),(WWW(kk)))
      call mfill(l+14,(etalb),(WWW(kk)))
c
      call mfill(l+15,(detallb),(WWW(kk)))
      call mfill(l+16,(azi),(WWW(kk)))
      if(azinorm.gt.0.d0)
     #  call mfill(l+17,(log10(azinorm)),(WWW(kk)))
      call mfill(l+18,(xmll),(WWW(kk)))
      call mfill(l+19,(ptpair),(WWW(kk)))
      if(ptpair.gt.0)call mfill(l+20,(log10(ptpair)),(WWW(kk)))
      call mfill(l+21,(0d0),(WWW(kk)))
c
      l=l+21

      if(abs(etav).lt.ycut)then
        call mfill(l+1,(ptv),(WWW(kk)))
        call mfill(l+2,(ptv),(WWW(kk)))
        if(ptv.gt.0.d0)call mfill(l+3,(log10(ptv)),(WWW(kk)))
      endif
      if(ptv.gt.20.d0)then
        call mfill(l+4,(yv),(WWW(kk)))
        call mfill(l+5,(etav),(WWW(kk)))
      endif
      if(abs(etav).lt.ycut.and.ptv.gt.20.d0)then
         call mfill(l+6,(xmv),(WWW(kk)))
         call mfill(l+21,(0d0),(WWW(kk)))
      endif
c
      if(abs(etal).lt.ycut)then
        call mfill(l+7,(ptl),(WWW(kk)))
        call mfill(l+8,(ptl),(WWW(kk)))
        if(ptl.gt.0.d0)call mfill(l+9,(log10(ptl)),(WWW(kk)))
      endif
      if(ptl.gt.20.d0)call mfill(l+10,(etal),(WWW(kk)))
      if(abs(etalb).lt.ycut)then
        call mfill(l+11,(ptlb),(WWW(kk)))
        call mfill(l+12,(ptlb),(WWW(kk)))
        if(ptlb.gt.0.d0)call mfill(l+13,(log10(ptlb)),(WWW(kk)))
      endif
      if(ptlb.gt.20.d0)call mfill(l+14,(etalb),(WWW(kk)))
c
      if( abs(etal).lt.ycut.and.abs(etalb).lt.ycut .and.
     #    ptl.gt.20.d0.and.ptlb.gt.20.d0)then
        call mfill(l+15,(detallb),(WWW(kk)))
        call mfill(l+16,(azi),(WWW(kk)))
        if(azinorm.gt.0.d0)
     #    call mfill(l+17,(log10(azinorm)),(WWW(kk)))
        call mfill(l+18,(xmll),(WWW(kk)))
        call mfill(l+19,(ptpair),(WWW(kk)))
        if(ptpair.gt.0) 
     #    call mfill(l+20,(log10(ptpair)),(WWW(kk)))
      endif

      enddo

 999  END

C-----------------------------------------------------------------------
      SUBROUTINE VVSUM(N,P,Q,R)
C-----------------------------------------------------------------------
C    VECTOR SUM
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N,I
      DOUBLE PRECISION P(N),Q(N),R(N)
      DO 10 I=1,N
   10 R(I)=P(I)+Q(I)
      END



C-----------------------------------------------------------------------
      SUBROUTINE VSCA(N,C,P,Q)
C-----------------------------------------------------------------------
C     VECTOR TIMES SCALAR
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N,I
      DOUBLE PRECISION C,P(N),Q(N)
      DO 10 I=1,N
   10 Q(I)=C*P(I)
      END

C-----------------------------------------------------------------------
      FUNCTION VDOT(N,P,Q)
C-----------------------------------------------------------------------
C     VECTOR DOT PRODUCT
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N,I
      DOUBLE PRECISION VDOT,PQ,P(N),Q(N)
      PQ=0.
      DO 10 I=1,N
   10 PQ=PQ+P(I)*Q(I)
      VDOT=PQ
      END

      function getrapidity(en,pl)
      implicit none
      real*8 getrapidity,en,pl,tiny,xplus,xminus,y
      parameter (tiny=1.d-8)
      xplus=en+pl
      xminus=en-pl
      if(xplus.gt.tiny.and.xminus.gt.tiny)then
         if( (xplus/xminus).gt.tiny.and.(xminus/xplus).gt.tiny)then
            y=0.5d0*log( xplus/xminus  )
         else
            y=sign(1.d0,pl)*1.d8
         endif
      else 
         y=sign(1.d0,pl)*1.d8
      endif
      getrapidity=y
      return
      end


      function getpseudorap(en,ptx,pty,pl)
      implicit none
      real*8 getpseudorap,en,ptx,pty,pl,tiny,pt,eta,th
      parameter (tiny=1.d-5)
c
      pt=sqrt(ptx**2+pty**2)
      if(pt.lt.tiny.and.abs(pl).lt.tiny)then
        eta=sign(1.d0,pl)*1.d8
      else
        th=atan2(pt,pl)
        eta=-log(tan(th/2.d0))
      endif
      getpseudorap=eta
      return
      end


      function getinvm(en,ptx,pty,pl)
      implicit none
      real*8 getinvm,en,ptx,pty,pl,tiny,tmp
      parameter (tiny=1.d-5)
c
      tmp=en**2-ptx**2-pty**2-pl**2
      if(tmp.gt.0.d0)then
        tmp=sqrt(tmp)
      elseif(tmp.gt.-tiny)then
        tmp=0.d0
      else
        write(*,*)'Attempt to compute a negative mass'
        stop
      endif
      getinvm=tmp
      return
      end


      function getdelphi(ptx1,pty1,ptx2,pty2)
      implicit none
      real*8 getdelphi,ptx1,pty1,ptx2,pty2,tiny,pt1,pt2,tmp
      parameter (tiny=1.d-5)
c
      pt1=sqrt(ptx1**2+pty1**2)
      pt2=sqrt(ptx2**2+pty2**2)
      if(pt1.ne.0.d0.and.pt2.ne.0.d0)then
        tmp=ptx1*ptx2+pty1*pty2
        tmp=tmp/(pt1*pt2)
        if(abs(tmp).gt.1.d0+tiny)then
          write(*,*)'Cosine larger than 1'
          stop
        elseif(abs(tmp).ge.1.d0)then
          tmp=sign(1.d0,tmp)
        endif
        tmp=acos(tmp)
      else
        tmp=1.d8
      endif
      getdelphi=tmp
      return
      end
