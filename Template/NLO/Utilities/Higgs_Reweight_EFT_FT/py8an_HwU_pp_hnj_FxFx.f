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
c
c     The type suffix of the histogram title, with syntax 
c     |T@<type_name> is semantic in the HwU format. It allows for
c     various filtering when using the histogram.py module
c     (see comment at the beginning of this file).
c     It is in general a good idea to keep the same title for the
c     same observable (if they use the same range) and differentiate
c     them only using the type suffix.
c
      character*12 HwUtype
      data HwUtype/'|T@M30MFT012'/
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
      k=0
      call HwU_book(k+ 1,'Higgs pt        '//HwUtype,40,0.d0,400.d0)
      call HwU_book(k+ 2,'Higgs y         '//HwUtype,50,-5.d0,5.d0)
      call HwU_book(k+ 3,'== jets         '//HwUtype,10,-0.5d0,9.5d0)
      call HwU_book(k+ 4,'>= jets         '//HwUtype,10,-0.5d0,9.5d0)
      call HwU_book(k+ 5,'Higgs pt == 0j  '//HwUtype,20,0.d0,400.d0)
      call HwU_book(k+ 6,'Higgs pt == 1j  '//HwUtype,20,0.d0,400.d0)
      call HwU_book(k+ 7,'Higgs pt == 2j  '//HwUtype,20,0.d0,400.d0)
      call HwU_book(k+ 8,'Higgs pt >= 3j  '//HwUtype,20,0.d0,400.d0)
      call HwU_book(k+ 9,'pt[j1]          '//HwUtype,30,0.d0,300.d0)
      call HwU_book(k+10,'pt[j2]          '//HwUtype,30,0.d0,300.d0)
      call HwU_book(k+11,'pt[j3]          '//HwUtype,30,0.d0,300.d0)
      call HwU_book(k+12,'y[j1]           '//HwUtype,20,-5.d0,5.d0)
      call HwU_book(k+13,'y[j2]           '//HwUtype,20,-5.d0,5.d0)
      call HwU_book(k+14,'y[j3]           '//HwUtype,20,-5.d0,5.d0)
      call HwU_book(k+15,'mjj             '//HwUtype,40,0.d0,1000.d0)
      call HwU_book(k+16,'Delta y[j1,j2]  '//HwUtype,40,0.d0,10.d0)
      call HwU_book(k+17,'Delta phi[j1,j2]'//HwUtype,20,0.d0,pi)
      call HwU_book(k+18,'cross section   '//HwUtype,10,-0.5d0,9.5d0)
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
      INTEGER ITMP,NN,NJET,NJET30,JET(NMAX),IORJET(NMAX)
      DOUBLE PRECISION palg,rfj,sycut,PTU,PTL,PTCALC,ptj,yj, PP(4,NMAX)
     $     ,PJET(4,NMAX),d01,d12,d23,d34,fastjetdmergemax, getrapidity
     $     ,getpseudorap,etaj,getdelphi,ptj1,ptj2,ptj3,dphi12
     $     ,dphi13,dphi23,deta12,invm,mjj
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
c$$$         if (i.eq.1) then 
c$$$            WWW(I)=EVWGT*ww(2)/ww(1)
c$$$         else
c$$$            WWW(I)=EVWGT*ww(i)/ww(1)
c$$$         endif
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

C FILL THE HISTOS
      xmh=pph(5)
      pth=sqrt(pph(1)**2+pph(2)**2)
      yh=getrapidity(pph(4),pph(3))

c
      kk=0
         
      palg=-1.d0
      rfj=0.4d0
      sycut=30.d0
      call fastjetppgenkt(pp,nn,rfj,sycut,palg,pjet,njet,jet)
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

      call HwU_fill(kk+1,pth,WWW)
      call HwU_fill(kk+2,yh,WWW)
      call HwU_fill(kk+3,dble(njet),WWW)
      do i=njet,0,-1
         call HwU_fill(kk+4,dble(i),WWW)
      enddo
         
      if (njet.eq.0) call HwU_fill(kk+5,pth,WWW)
      if (njet.eq.1) call HwU_fill(kk+6,pth,WWW)
      if (njet.eq.2) call HwU_fill(kk+7,pth,WWW)
      if (njet.ge.3) call HwU_fill(kk+8,pth,WWW)

      do i=1,njet
         ptj=PTCALC(PJET(1,i))
         etaj=getrapidity(PJET(4,i),PJET(3,I))
         call HwU_fill(kk+8+min(i,3),ptj,WWW)
         call HwU_fill(kk+11+min(i,3),etaj,WWW)
      enddo
      if (njet.ge.2) then
         mjj=invm(pjet(1,1),pjet(1,2))
         deta12=abs(getrapidity(pjet(4,1),pjet(3,1))-
     $              getrapidity(pjet(4,2),pjet(3,2)))
         dphi12=getdelphi(PJET(1,1),PJET(2,1),PJET(1,2),PJET(2,2))

         call HwU_fill(kk+15,mjj,WWW)
         call HwU_fill(kk+16,deta12,WWW)
         call HwU_fill(kk+17,dphi12,WWW)
      endif
      
c total cross section
      call HwU_fill(kk+18,0d0,WWW)
      if (njet.ge.2) then
         if (mjj.gt.400d0 .and. deta12.ge.2.8d0) then
            call HwU_fill(kk+18,1d0,WWW)
            if (njet.eq.2) call HwU_fill(kk+18,2d0,WWW)
         endif
         if (mjj.gt.600d0 .and. deta12.ge.4.0d0) then
            call HwU_fill(kk+18,3d0,WWW)
            if (njet.eq.2) call HwU_fill(kk+18,4d0,WWW)
         endif
      endif
         
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
      double precision p1(4),p2(4),getinvm
      invm=getinvm(p1(4)+p2(4),p1(1)+p2(1),p1(2)+p2(2),p1(3)+p2(3))
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



