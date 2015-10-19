c
c Example analysis for "p p > t t~ [QCD]" process.
c
c It features the HwU format for histogram booking and output.
c The details of how to process/manipulate the resulting .HwU file,
c in particular how to plot it using gnuplot, I refer the reader to this
c FAQ:
c
c      https://answers.launchpad.net/mg5amcnlo/+faq/2671
c
c It mostly relies on using the following madgraph5 module in standalone
c
c  <MG5_aMC_install_dir>/madgraph/various/histograms.py
c
c You can learn about how to run it and what options are available with
c
c  python <MG5_aMC_install_dir>/madgraph/various/histograms.py --help
c
C----------------------------------------------------------------------
      SUBROUTINE RCLOS()
C     DUMMY IF HBOOK IS USED
C----------------------------------------------------------------------
      END


C----------------------------------------------------------------------
      SUBROUTINE PYABEG(nnn,wwwi)
C     USER''S ROUTINE FOR INITIALIZATION
C----------------------------------------------------------------------
      INCLUDE 'HEPMC.INC'
      include 'reweight0.inc'
      REAL*8 pi
      PARAMETER (PI=3.14159265358979312D0)
      integer j,kk,l,i,nnn
      character*5 Cuts(4)
      data Cuts/' inc ',' 2jet',' vbf1',' vbf2'/
      integer nwgt,max_weight,nwgt_analysis
      common/cnwgt/nwgt
      common/c_analysis/nwgt_analysis
      parameter (max_weight=maxscales*maxscales+maxpdfs+1)
      character*15 weights_info(max_weight),wwwi(max_weight)
      common/cwgtsinfo/weights_info
      
c
      weights_info(1)="central value  "
      do i=1,nnn+1
         weights_info(i+1)=wwwi(i)
      enddo
      nwgt=nnn+1
c Initialize histograms
      call HwU_inithist(nwgt,weights_info)
c Set method for error estimation to '0', i.e., use Poisson statistics
c for the uncertainty estimate
      call set_error_estimation(0)
      nwgt_analysis=nwgt
      call HwU_book(1,'cross section   ',10,-0.5d0,9.5d0)
c Loop over the 4 sets of cuts ("inclusive", "at least two jets", 
c "loose VBF cuts" and "tight VBF cuts")
      do i=1,4
      k=(i-1)*27+1
      call HwU_book(k+ 1,'Higgs pt        '//Cuts(i),40,0.d0,400.d0)
      call HwU_book(k+ 2,'Higgs y         '//Cuts(i),50,-5.d0,5.d0)
      call HwU_book(k+ 3,'== jets         '//Cuts(i),10,-0.5d0,9.5d0)
      call HwU_book(k+ 4,'>= jets         '//Cuts(i),10,-0.5d0,9.5d0)
      call HwU_book(k+ 5,'Higgs pt == 0j  '//Cuts(i),20,0.d0,400.d0)
      call HwU_book(k+ 6,'Higgs pt == 1j  '//Cuts(i),20,0.d0,400.d0)
      call HwU_book(k+ 7,'Higgs pt == 2j  '//Cuts(i),20,0.d0,400.d0)
      call HwU_book(k+ 8,'Higgs pt >= 3j  '//Cuts(i),20,0.d0,400.d0)
      call HwU_book(k+ 9,'d0              '//Cuts(i),40,0.d0,500.d0)
      call HwU_book(k+10,'log10[d0]       '//Cuts(i),40,-0.5d0,3.5d0)
      call HwU_book(k+11,'d1              '//Cuts(i),40,0.d0,500.d0)
      call HwU_book(k+12,'log10[d1]       '//Cuts(i),40,-0.5d0,3.5d0)
      call HwU_book(k+13,'d2              '//Cuts(i),40,0.d0,500.d0)
      call HwU_book(k+14,'log10[d2]       '//Cuts(i),40,-0.5d0,3.5d0)
      call HwU_book(k+15,'d3              '//Cuts(i),40,0.d0,500.d0)
      call HwU_book(k+16,'log10[d3]       '//Cuts(i),40,-0.5d0,3.5d0)
      call HwU_book(k+17,'pt[j1]          '//Cuts(i),30,0.d0,300.d0)
      call HwU_book(k+18,'pt[j2]          '//Cuts(i),30,0.d0,300.d0)
      call HwU_book(k+19,'pt[j3]          '//Cuts(i),30,0.d0,300.d0)
      call HwU_book(k+20,'pt[j4]          '//Cuts(i),30,0.d0,300.d0)
      call HwU_book(k+21,'eta[j1]         '//Cuts(i),20,-5.d0,5.d0)
      call HwU_book(k+22,'eta[j2]         '//Cuts(i),20,-5.d0,5.d0)
      call HwU_book(k+23,'eta[j3]         '//Cuts(i),20,-5.d0,5.d0)
      call HwU_book(k+24,'eta[j4]         '//Cuts(i),20,-5.d0,5.d0)
      call HwU_book(k+25,'mjj             '//Cuts(i),40,0.d0,1000.d0)
      call HwU_book(k+26,'Delta y[j1,j2]  '//Cuts(i),40,0.d0,10.d0)
      call HwU_book(k+27,'Delta phi[j1,j2]'//Cuts(i),20,0.d0,pi)
      enddo
 999  END

C----------------------------------------------------------------------
      SUBROUTINE PYAEND(IEVTTOT)
C     USER''S ROUTINE FOR TERMINAL CALCULATIONS, HISTOGRAM OUTPUT, ETC
C----------------------------------------------------------------------
      INCLUDE 'HEPMC.INC'
      REAL*8 XNORM,IEVTTOT
      INTEGER I,J,KK,l,nwgt_analysis
      common/c_analysis/nwgt_analysis
c Collect accumulated results. IEVTTOT is such that we need to multiply
c the results by this factor
      xnorm=ievttot
      call finalize_histograms(nevhep)
c Write the histograms to disk. 
      open (unit=99,file='MADatNLO.HwU',status='unknown')
      call HwU_output(99,xnorm)
      close (99)
      END

C----------------------------------------------------------------------
      SUBROUTINE PYANAL(nnn,xww)
C     USER''S ROUTINE TO ANALYSE DATA FROM EVENT
C----------------------------------------------------------------------
      INCLUDE 'HEPMC.INC'
      include 'reweight0.inc'
      DOUBLE PRECISION PSUM(4),PPH(5),XMH,PTH,YH
      INTEGER ICHSUM,ICHINI,IHEP,IFH,IST,ID,IJ,I,J,II
      INTEGER NMAX
      PARAMETER (NMAX=2000)
      INTEGER NN,NJET,JET(NMAX)
      DOUBLE PRECISION palg,rfj,sycut,PTCALC,ptj,PP(4,NMAX)
     &     ,PJET(4,NMAX),d01,d12,d23,d34,fastjetdmergemax,getrapidity
     &     ,getpseudorap,etaj,getdelphi,ptj1,ptj2,ptj3,dphi12,deta12
     &     ,invm,mjj
      REAL*8 WWW0,TINY
      INTEGER KK,KK1,ibin
      DATA TINY/.1D-5/
      integer nwgt_analysis,max_weight
      common/c_analysis/nwgt_analysis
      parameter (max_weight=maxscales*maxscales+maxpdfs+1)
      double precision ww(max_weight),www(max_weight),xww(max_weight)
      common/cww/ww
      integer icount
      data icount /0/
c
      ww(1)=xww(2)
      if(nnn.eq.0)ww(1)=1d0
      do i=2,nnn+1
         ww(i)=xww(i)
      enddo
c
      IF (WW(1).EQ.0D0) THEN
         WRITE(*,*)'WW(1) = 0. Stopping'
         STOP
      ENDIF
C INCOMING PARTONS MAY TRAVEL IN THE SAME DIRECTION: IT''S A POWER-SUPPRESSED
C EFFECT, SO THROW THE EVENT AWAY
      IF(SIGN(1.D0,PHEP(3,1)).EQ.SIGN(1.D0,PHEP(3,2)))THEN
         WRITE(*,*)'WARNING 111 IN PYANAL'
         GOTO 999
      ENDIF
      DO I=1,nwgt_analysis
         www(i)=ww(i)
      ENDDO

      NN=0
      IFH=0
      DO 100 IHEP=1,NHEP
        IST=ISTHEP(IHEP)      
        ID=IDHEP(IHEP)
        IF(ist.eq.1)THEN
          IF(ID.EQ.25)THEN
            IFH=IFH+1
            DO IJ=1,5
	      PPH(IJ)=PHEP(IJ,IHEP)
	    ENDDO
          ENDIF
        ENDIF
        IF (IST.EQ.1 .AND. (ABS(ID).GT.100.or.ID.eq.22)) THEN
          NN=NN+1
          IF (NN.GT.NMAX) STOP 'Too many particles!'
          DO I=1,4
             PP(I,NN)=PHEP(I,IHEP)
          ENDDO
        ENDIF
  100 CONTINUE
      IF(IFH.NE.1)THEN
         write (*,*) 'ERROR 501 in pyanal'
         stop
      ENDIF
      
      xmh=pph(5)
      pth=sqrt(pph(1)**2+pph(2)**2)
      yh=getrapidity(pph(4),pph(3))

c Call fastjet with kT-algo to get the jet scales
      palg=1.d0
      rfj=1.0d0
      sycut=0.d0
      call fastjetppgenkt(pp,nn,rfj,sycut,palg,pjet,njet,jet)
      if (NN.ge.1) then
         d01=dsqrt(fastjetdmergemax(0))
      else
         d01=-1d8
      endif
      if (NN.ge.2) then
         d12=dsqrt(fastjetdmergemax(1))
      else
         d12=-1d8
      endif
      if (NN.ge.3) then
         d23=dsqrt(fastjetdmergemax(2))
      else
         d23=-1d8
      endif
      if (NN.ge.4) then
         d34=dsqrt(fastjetdmergemax(3))
      else
         d34=-1d8
      endif

c Call fastjet with anti-kT to get the actual jets
      palg=-1.d0
      rfj=0.4d0
      sycut=30.d0
      njet=0
      call fastjetppgenkt(pp,nn,rfj,sycut,palg,pjet,njet,jet)
c Jets are already ordered in pT; remove jets that are outside of
c rapidity range
      do i=1,njet
         if (abs(getrapidity(pjet(4,i),pjet(3,i))).gt.4.4d0) then
            do j=i,njet-1
               do k=1,4
                  pjet(k,j)=pjet(k,j+1)
               enddo
            enddo
            njet=njet-1
         endif
      enddo

      do ii=1,4
         kk=(ii-1)*27+1
c Apply the cuts:
         if (ii.eq.2) then
            if (njet.lt.2) cycle
         endif
         if (ii.eq.3) then
            if (njet.lt.2) cycle
            if (mjj.lt.400d0 .or. deta12.lt.2.8d0) cycle
         endif
         if (ii.eq.4) then
            if (njet.lt.2) cycle
            if (mjj.lt.600d0 .or. deta12.lt.4.0d0) cycle
         endif

c Fill the cross section histogram         
         if (ii.eq.1) then
            call HwU_fill(1,0d0,WWW)
            if (njet.eq.2) call HwU_fill(1,1d0,WWW)
         elseif (ii.eq.2) then
            call HwU_fill(1,2d0,WWW)
            if (njet.eq.2) call HwU_fill(1,3d0,WWW)
         elseif (ii.eq.3) then
            call HwU_fill(1,4d0,WWW)
            if (njet.eq.2) call HwU_fill(1,5d0,WWW)
         elseif (ii.eq.4) then
            call HwU_fill(1,6d0,WWW)
            if (njet.eq.2) call HwU_fill(1,7d0,WWW)
         endif

c Higgs and number of jets plots
         call HwU_fill(kk+1,pth,WWW)
         call HwU_fill(kk+2,yh,WWW)
         call HwU_fill(kk+3,dble(njet),WWW)
         do i=0,njet
            call HwU_fill(kk+4,dble(i),WWW)
         enddo
         if (njet.eq.0) call HwU_fill(kk+5,pth,WWW)
         if (njet.eq.1) call HwU_fill(kk+6,pth,WWW)
         if (njet.eq.2) call HwU_fill(kk+7,pth,WWW)
         if (njet.ge.3) call HwU_fill(kk+8,pth,WWW)


c Diff. jet rates
         call HwU_fill(kk+9,d01,www)
         if (d01.ge.1d-3)
     &        call HwU_fill(kk+10,log10(d01),www)
         call HwU_fill(kk+11,d12,www)
         if (d12.ge.1d-3)
     &        call HwU_fill(kk+12,log10(d12),www)
         call HwU_fill(kk+13,d23,www)
         if (d23.ge.1d-3)
     &        call HwU_fill(kk+14,log10(d23),www)
         call HwU_fill(kk+15,d34,www)
         if (d34.ge.1d-3)
     &        call HwU_fill(kk+16,log10(d34),www)

         
c below are jet-observables. Need to have at least one jet
         if (njet.eq.0) cycle

         do i=1,njet
            if (i.eq.5) exit ! only fill 4 plots
            ptj=PTCALC(PJET(1,i))
            etaj=getrapidity(PJET(4,i),PJET(3,I))
            call HwU_fill(kk+16+i,ptj,WWW)
            call HwU_fill(kk+20+i,etaj,WWW)
         enddo

c below are 2-jet-observables. Need to have at least two jets
         if (njet.lt.2) cycle

         mjj=invm(pjet(1,1),pjet(1,2))
         deta12=abs(getrapidity(pjet(4,1),pjet(3,1))-
     $        getrapidity(pjet(4,2),pjet(3,2)))
         dphi12=getdelphi(PJET(1,1),PJET(2,1),PJET(1,2),PJET(2,2))
         call HwU_fill(kk+25,mjj,WWW)
         call HwU_fill(kk+26,deta12,WWW)
         call HwU_fill(kk+27,dphi12,WWW)
      enddo

      call HwU_add_points

 999  return
      end


      function getrapidity(en,pl)
      implicit none
      real*8 getrapidity,en,pl,tiny,xplus,xminus,y
      parameter (tiny=1.d-8)
c
      xplus=en+pl
      xminus=en-pl
      if(xplus.gt.tiny.and.xminus.gt.tiny)then
        if( (xplus/xminus).gt.tiny )then
          y=0.5d0*log( xplus/xminus )
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


      double precision function invm(p1,p2)
      implicit none
      integer i
      double precision p1(4),p2(4),tmp(4),getinvm
      do i=1,4
         tmp(i)=p1(i)+p2(i)
      enddo
      invm=getinvm(tmp(4),tmp(1),tmp(2),tmp(3))
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



      FUNCTION PTCALC(P)
      IMPLICIT NONE
      DOUBLE PRECISION PTCALC,P(4),PTSQ
      PTSQ=P(1)**2+P(2)**2
      IF (PTSQ.EQ.0D0) THEN
         PTCALC=0D0
      ELSE
         PTCALC=SQRT(PTSQ)
      ENDIF
      END



