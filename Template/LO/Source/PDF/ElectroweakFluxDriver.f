c     /* ********************************************************* *
c                 Effective Vector Boson Approximation
c     /* ********************************************************* *      
c     File: ElectroweakFluxDriver.f
c     R. Ruiz (2021 February)
c     For details, see companion paper by Costantini, et al [arXiv:]
c     /* ********************************************************* *
c     function eva_get_pdf_by_PID:
c     - wrapper for eva_get_pdf_by_PID_evo
c     function eva_get_pdf_by_PID_evo
c     - set eva PDF couplings by PIDs
c     - call V_+,V_-,V_0 PDF by v polarization (vpol)
c     - call PDF for f_L,f_R by fL polarization (fLpol; fLpol=0.5 = unpolarized)
c     subroutine eva_get_mv2_by_PID
c     - assign mass by vPID
c     subroutine eva_get_mf2_by_PID      
c     - assign mass by fPID
c     subroutine eva_get_gg2_by_PID
c     - assign universal coupling strength by vPID
c     subroutine eva_get_gR2_by_PID
c     - assign right couplings of fermion by vPID and fPID      
c     subroutine eva_get_gL2_by_PID
c     - assign left couplings of fermion by vPID and fPID
c     /* ********************************************************* *
      double precision function eva_get_pdf_by_PID(vPID,fPID,vpol,fLpol,x,mu2)
      implicit none
      integer vPID,fPID,vpol
      double precision fLpol,x,mu2
      double precision eva_get_pdf_by_PID_evo
      double precision eva_get_pdf_photon_evo
      include 'ElectroweakFlux.inc'
      
      integer ievo
      double precision tiny,mu2min
      ievo  = 0 ! =0 for evolution by q^2 (!=0 for evolution by pT^2)
      tiny  = 1d-8
      mu2min = 1d2 ! (10 GeV)^2 reset mu2min by vPID  

c     do the following checks before calling PDF:
c     1. momentum fraction, x
c     2. fermion polarization fraction, fLpol
c     3. vector boson polarization by PID, vpol vPID
c     4. evolution scale, mu2
c     start checks
c     1. check momentum fraction
      if(x.lt.tiny.or.x.gt.(1d0-tiny)) then
         write(*,*) 'eva: x out of range',x
         eva_get_pdf_by_PID = 0d0
         return
      endif
c     2. check fermion polarization fraction
      if(fLpol.lt.0d0.or.fLpol.gt.1d0) then
         write(*,*) 'eva: fLpol out of range',fLpol
         stop
         eva_get_pdf_by_PID = 0d0
         return
      endif
c     3. check vector boson polarization by PID
c     also set by PID lower bound no mu2 scale evolution
      select case (iabs(vPID))
      case (23) ! z
         mu2min = eva_mz2
         if(iabs(vPol).ne.1.and.vPol.ne.0) then
            write(*,*) 'vPol out of range',vPol
            stop
            eva_get_pdf_by_PID = 0d0
            return
         endif
      case (24) ! w
         mu2min = eva_mw2
         if(iabs(vPol).ne.1.and.vPol.ne.0) then
            write(*,*) 'vPol out of range',vPol
            stop
            eva_get_pdf_by_PID = 0d0
            return
         endif
      case (7,22) ! photon (special treatment for mu2min)
         if(iabs(fPID).lt.5) then ! catch light quarks
            mu2min = (2d0)**2   ! set evolution boundary; 2 GeV is arbitrary
         else ! all others
            call eva_get_mf2_by_PID(mu2min,fPID) ! set scale to mass of parent fermion
         endif
         if(iabs(vPol).ne.1) then
            write(*,*) 'vPol out of range',vPol
            stop
            eva_get_pdf_by_PID = 0d0
            return
         endif         
c      case (32) (eva for bsm)
c         mu2min = eva_mx2
c         if(iabs(vPol).ne.1.and.vPol.ne.0) then
c            write(*,*) 'vPol out of range',vPol
c            stop
c            eva_get_pdf_by_PID = 0d0
c            return
c         endif
      case default         
         write(*,*) 'vPID out of range',vPID
         stop
         eva_get_pdf_by_PID = 0d0
         return
      end select
c     4. check evolution scale
      if(mu2.lt.mu2min) then
         write(*,*) 'mu2 too small. setting to mu2min',mu2,mu2min
         mu2 = mu2min
         eva_get_pdf_by_PID = 0d0
         return
      endif

c     celebrate by calling the PDF
      if(vPID.eq.22.or.vPID.eq.7) then
         eva_get_pdf_by_PID = eva_get_pdf_photon_evo(vPID,fPID,vpol,fLpol,x,mu2,ievo)
      else
         eva_get_pdf_by_PID = eva_get_pdf_by_PID_evo(vPID,fPID,vpol,fLpol,x,mu2,ievo)
      endif
      return
      end      
c     /* ********************************************************* *
c     /* ********************************************************* *
c     /* ********************************************************* *
c     /* ********************************************************* *
      double precision function eva_get_pdf_by_PID_evo(vPID,fPID,vpol,fLpol,x,mu2,ievo)
      implicit none
      integer vPID,fPID,vpol,ievo
      double precision fLpol,x,mu2
      double precision eva_fX_to_vm,eva_fX_to_v0,eva_fX_to_vp
      
      double precision gg2,gL2,gR2,mv2,tmpPDF
      call eva_get_mv2_by_PID(mv2,vPID)
      call eva_get_gg2_by_PID(gg2,vPID,fPID)
      if( fPID/iabs(fPID).gt.0 ) then ! particle
         call eva_get_gR2_by_PID(gR2,vPID,fPID)
         call eva_get_gL2_by_PID(gL2,vPID,fPID)
      else  ! antiparticle (invert parity)
         call eva_get_gR2_by_PID(gL2,vPID,fPID)
         call eva_get_gL2_by_PID(gR2,vPID,fPID)
      endif
      select case (vpol)
      case (-1)
         tmpPDF = eva_fX_to_vm(gg2,gL2,gR2,fLpol,mv2,x,mu2,ievo)
      case (0)
         tmpPDF = eva_fX_to_v0(gg2,gL2,gR2,fLpol,mv2,x,mu2,ievo)
      case (+1)
         tmpPDF = eva_fX_to_vp(gg2,gL2,gR2,fLpol,mv2,x,mu2,ievo)
      case default
         write(*,*) 'vPol out of range; should not be here',vPol
         stop
         tmpPDF = 0d0
      end select
      eva_get_pdf_by_PID_evo = tmpPDF
      return
      end
c     /* ********************************************************* *
      double precision function eva_get_pdf_photon_evo(vPID,fPID,vpol,fLpol,x,mu2,ievo)
      implicit none
      integer vPID,fPID,vpol,ievo
      double precision fLpol,x,mu2
      double precision eva_fX_to_vm,eva_fX_to_v0,eva_fX_to_vp
      
      double precision gg2,gL2,gR2,mf2,tmpPDF
      call eva_get_mf2_by_PID(mf2,fPID)
      call eva_get_gg2_by_PID(gg2,vPID,fPID)
      if( fPID/iabs(fPID).gt.0 ) then ! particle
         call eva_get_gR2_by_PID(gR2,vPID,fPID)
         call eva_get_gL2_by_PID(gL2,vPID,fPID)
      else  ! antiparticle (invert parity)
         call eva_get_gR2_by_PID(gL2,vPID,fPID)
         call eva_get_gL2_by_PID(gR2,vPID,fPID)
      endif
      select case (vpol)
      case (-1)
         tmpPDF = eva_fX_to_vm(gg2,gL2,gR2,fLpol,mf2,x,mu2,ievo)
      case (+1)
         tmpPDF = eva_fX_to_vp(gg2,gL2,gR2,fLpol,mf2,x,mu2,ievo)
      case default
         write(*,*) 'vPol out of range; should not be here',vPol
         stop
         tmpPDF = 0d0         
      end select
      eva_get_pdf_photon_evo = tmpPDF
      return
      end      
c     /* ********************************************************* *      
c     /* ********************************************************* *
c     /* ********************************************************* *
c     /* ********************************************************* *
c     /* ********************************************************* *                  
      subroutine eva_get_mv2_by_PID(mv2,vPID)
      implicit none
      integer vPID
      double precision mv2
      include 'ElectroweakFlux.inc'
      
      select case (iabs(vPID))
      case (7,22)
         mv2 = eva_zero
      case (23)
         mv2 = eva_mz2
      case (24)
         mv2 = eva_mw2
c      case (25)
c         mv2 = eva_mh2
c      case (32)
c         mv2 = eva_mx2
      case default
         write(*,*) 'eva: setting m_v to m_w. unknown vPID:', vPID
         mv2 = eva_mw2
      end select
      return
      end
c     /* ********************************************************* *
c     /* ********************************************************* *                  
      subroutine eva_get_mf2_by_PID(mf2,fPID)
      implicit none
      integer fPID
      double precision mf2
      include 'ElectroweakFlux.inc'
      
      select case (iabs(fPID))
      case (1)
         mf2 = eva_md2
      case (2)
         mf2 = eva_mu2
      case (3)
         mf2 = eva_ms2
      case (4)
         mf2 = eva_mc2
      case (5)
         mf2 = eva_mb2
      case (6)
         mf2 = eva_mt2
      case (11)
         mf2 = eva_me2
      case (13)
         mf2 = eva_mm2
      case (15)
         mf2 = eva_ml2  
      case default
         write(*,*) 'eva: setting m_f to m_e. unknown fPID:', fPID
         mf2 = eva_me2
      end select
      return
      end
c     /* ********************************************************* *            
c     /* ********************************************************* *      
      subroutine eva_get_gg2_by_PID(gg2,vPID,fPID)
      implicit none
      integer vPID,fPID
      double precision gg2
      include 'ElectroweakFlux.inc'
      
      select case (iabs(vPID))
      case (7,22)         
c     ******************************
         select case (iabs(fPID)) ! nested select case
         case (1)               ! down
            gg2 = eva_ee2*eva_qed2 ! = e^2 * (-1/3)^2
         case (2)               ! up
            gg2 = eva_ee2*eva_qeu2
         case (3)               ! strange
            gg2 = eva_ee2*eva_qed2
         case (4)               ! charm
            gg2 = eva_ee2*eva_qeu2
         case (5)               ! bottom
            gg2 = eva_ee2*eva_qed2
         case (6)               ! top
            gg2 = eva_ee2*eva_qeu2
         case (11)              ! electron
            gg2 = eva_ee2*eva_qee2
         case (12)              ! electron-neutrino
            write(*,*) 'eva: nu_e has zero QED charge.'
            gg2 = eva_zero
         case (13)              ! muon
            gg2 = eva_ee2*eva_qee2
         case (14)              ! muon-neutrino
            write(*,*) 'eva: nu_m has zero QED charge.'
            gg2 = eva_zero
         case (15)              ! tau
            gg2 = eva_ee2*eva_qee2
         case (16)              ! tau-neutrino
            write(*,*) 'eva: nu_tau has zero QED charge.'
            gg2 = eva_zero
         case default
            write(*,*) 'eva: setting QED coup to (e*Q_e). unknown fPID:', fPID
            gg2 = eva_ee2*eva_qee2
         end select
c     ******************************                     
      case (23)
         gg2 = eva_gz2
      case (24)
         gg2 = eva_gw2/2.d0
      case default
         gg2 = eva_gw2
      end select
      return
      end
c     /* ********************************************************* *      
c     /* ********************************************************* *      
      subroutine eva_get_gR2_by_PID(gR2,vPID,fPID)
      implicit none
      integer vPID,fPID
      double precision gR2
      include 'ElectroweakFlux.inc'
      
      select case (iabs(vPID))
      case (7,22)
         gR2 = eva_one
      case (23)
c     ******************************
         select case (iabs(fPID)) ! nested select case
         case (1)               ! down
            gR2 = eva_zRd**2
         case (2)               ! up
            gR2 = eva_zRu**2
         case (3)               ! strange
            gR2 = eva_zRd**2
         case (4)               ! charm
            gR2 = eva_zRu**2
         case (5)               ! bottom
            gR2 = eva_zRd**2
         case (6)               ! top
            gR2 = eva_zRu**2            
         case (11)              ! electron
            gR2 = eva_zRe**2
         case (12)              ! electron-neutrino
            gR2 = eva_zRv**2
         case (13)              ! muon
            gR2 = eva_zRe**2
         case (14)              ! muon-neutrino
            gR2 = eva_zRv**2
         case (15)              ! tau
            gR2 = eva_zRe**2
         case (16)              ! tau-neutrino
            gR2 = eva_zRv**2                        
         case default
            gR2 = eva_one**2
         end select
c     ******************************            
      case (24)
         gR2 = eva_zero
      case default
         gR2 = eva_one
      end select
      return
      end
c     /* ********************************************************* *      
c     /* ********************************************************* *      
      subroutine eva_get_gL2_by_PID(gL2,vPID,fPID)
      implicit none
      integer vPID,fPID
      double precision gL2
      include 'ElectroweakFlux.inc'
      
      select case (iabs(vPID))
      case (7,22)
         gL2 = eva_one
      case (23)
c     ******************************
         select case (iabs(fPID)) ! nested select case
         case (1)               ! down
            gL2 = eva_zLd**2
         case (2)               ! up
            gL2 = eva_zLu**2
         case (3)               ! strange
            gL2 = eva_zLd**2
         case (4)               ! charm
            gL2 = eva_zLu**2
         case (5)               ! bottom
            gL2 = eva_zLd**2
         case (6)               ! top
            gL2 = eva_zLu**2            
         case (11)              ! electron
            gL2 = eva_zLe**2
         case (12)              ! electron-neutrino
            gL2 = eva_zLv**2
         case (13)              ! muon
            gL2 = eva_zLe**2
         case (14)              ! muon-neutrino
            gL2 = eva_zLv**2
         case (15)              ! tau
            gL2 = eva_zLe**2
         case (16)              ! tau-neutrino
            gL2 = eva_zLv**2                        
         case default
            gL2 = eva_one**2
         end select
c     ******************************            
      case (24)
         gL2 = eva_one
      case default
         gL2 = eva_one
      end select
      return
      end
c     /* ********************************************************* *      
