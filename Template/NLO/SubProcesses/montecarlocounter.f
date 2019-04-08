      subroutine set_mc_matrices
      implicit none
      include "genps.inc"
      include 'nexternal.inc'
      include "born_nhel.inc"
      integer idup(nexternal-1,maxproc)
      integer mothup(2,nexternal-1,maxproc)
      integer icolup(2,nexternal-1,max_bcol)
c Nexternal is the number of legs (initial and final) al NLO, while max_bcol
c is the number of color flows at Born level
      integer i,j,k,l,k0,mothercol(2),i1(2)
      include "born_leshouche.inc"
      integer ipartners(0:nexternal-1),colorflow(nexternal-1,0:max_bcol)
      common /MC_info/ ipartners,colorflow
      integer i_fks,j_fks
      common/fks_indices/i_fks,j_fks
      integer fksfather
      common/cfksfather/fksfather
      logical notagluon,found
      common/cnotagluon/notagluon
      integer nglu,nsngl
      logical isspecial(max_bcol)
      common/cisspecial/isspecial
      logical spec_case
c
      logical is_leading_cflow(max_bcol)
      integer num_leading_cflows
      common/c_leading_cflows/is_leading_cflow,num_leading_cflows
      integer iforest(2,-max_branch:-1,lmaxconfigs)
      integer mapconfig(0:lmaxconfigs)
      integer sprop(-max_branch:-1,lmaxconfigs)
      integer tprid(-max_branch:-1,lmaxconfigs)
      include 'born_conf.inc'
      include 'coloramps.inc'
c
      do i=0,nexternal-1
         ipartners(i)=0
      enddo
      do i=1,nexternal-1
         do j=0,max_bcol
            colorflow(i,j)=0
         enddo
      enddo
c ipartners(0): number of particles that can be colour or anticolour partner 
c   of the father, the Born-level particle to which i_fks and j_fks are 
c   attached. If one given particle is the colour/anticolour partner of
c   the father in more than one colour flow, it is counted only once
c   in ipartners(0)
c ipartners(i), 1<=i<=nexternal-1: the label (according to Born-level
c   labelling) of the i^th colour partner of the father
c
c colorflow(i,0), 1<=i<=nexternal-1: number of colour flows in which
c   the particle ipartners(i) is a colour partner of the father
c colorflow(i,j): the actual label (according to born_leshouche.inc)
c   of the j^th colour flow in which the father and ipartners(i) are
c   colour partners
c
c Example: in the process q(1) qbar(2) -> g(3) g(4), the two color flows are
c
c j=1    i    icolup(1)    icolup(2)       j=2    i    icolup(1)    icolup(2)
c        1      500           0                   1      500           0
c        2       0           501                  2       0           501
c        3      500          502                  3      502          501
c        4      502          501                  4      500          502
c
c and if one fixes for example fksfather=3, then the situation is the following.
c
c fksfather = 3
c  
c ipartners(0) = 3
c ipartners(1,2,3) = 1, 4, 2
c  
c colorflow(1,0) = 1 = number of flows where ipartners(1) = 1 is connected to 3
c colorflow(2,0) = 2 = number of flows where ipartners(2) = 4 is connected to 3
c colorflow(3,0) = 1 = number of flows where ipartners(3) = 2 is connected to 3
c colorflow(1,1) = 1 = flow where ipartners(1) = 1 is connected to 3
c colorflow(1,2) = 0 -> no other flow connecting 1 and 3
c colorflow(2,1) = 1 = first flow where ipartners(2) = 4 is connected to 3
c colorflow(2,2) = 2 = second flow where ipartners(2) = 4 is connected to 3
c colorflow(3,1) = 2 = flow where ipartners(3) = 2 is connected to 3
c colorflow(3,2) = 0 -> no other flow connecting 2 and 3
c colorflow(4,1) = 0 -> there is no fourth partner of 3
c colorflow(4,2) = 0 -> there is no fourth partner of 3
c  
c Thus
c
c ipartners(0..3) = 3, 1, 4, 2
c  
c colorflow(1,0..2) = 1, 1, 0
c colorflow(2,0..2) = 2, 1, 2
c colorflow(3,0..2) = 1, 2, 0
c colorflow(4,0..2) = 0, 0, 0

      fksfather=min(i_fks,j_fks)

c isspecial will be set equal to .true. colour flow by colour flow only
c if the father is a gluon, and another gluon will be found which is
c connected to it by both colour and anticolour
      isspecial=.false.
c
c consider only leading colour flows
      num_leading_cflows=0
      do i=1,max_bcol
         is_leading_cflow(i)=.false.
         do j=1,mapconfig(0)
            if(icolamp(i,j,1))then
               is_leading_cflow(i)=.true.
               num_leading_cflows=num_leading_cflows+1
               exit
            endif
         enddo
      enddo
c
      do i=1,max_bcol
         if(.not.is_leading_cflow(i))cycle
c Loop over Born-level colour flows
c nglu and nsngl are the number of gluons (except for the father) and of 
c colour singlets in the Born process, according to the information 
c stored in ICOLUP
         nglu=0
         nsngl=0
         mothercol(1)=ICOLUP(1,fksfather,i)
         mothercol(2)=ICOLUP(2,fksfather,i)
         notagluon=(mothercol(1).eq.0 .or. mothercol(2).eq.0)
c
         do j=1,nexternal-1
c Loop over Born-level particles; j is the possible colour partner of father,
c and whether this is the case is determined inside this loop
            if (j.ne.fksfather) then
c Skip father (it cannot be its own colour partner)
               if(ICOLUP(1,j,i).eq.0.and.ICOLUP(2,j,i).eq.0)
     #           nsngl=nsngl+1
               if(ICOLUP(1,j,i).ne.0.and.ICOLUP(2,j,i).ne.0)
     #           nglu=nglu+1
               if ( (j.le.nincoming.and.fksfather.gt.nincoming) .or.
     #              (j.gt.nincoming.and.fksfather.le.nincoming) ) then
c father and j not both in the initial or in the final state -- connect
c colour (1) with colour (i1(1)), and anticolour (2) with anticolour (i1(2))
                  i1(1)=1
                  i1(2)=2
               else
c father and j both in the initial or in the final state -- connect
c colour (1) with anticolour (i1(2)), and anticolour (2) with colour (i1(1))
                  i1(1)=2
                  i1(2)=1
               endif
               do l=1,2
c Loop over colour and anticolour of father
                  found=.false.
                  if( ICOLUP(i1(l),j,i).eq.mothercol(l) .and.
     &                ICOLUP(i1(l),j,i).ne.0 ) then
c When ICOLUP(i1(l),j,i) = mothercol(l), the colour (if i1(l)=1) or
c the anticolour (if i1(l)=2) of particle j is connected to the
c colour (if l=1) or the anticolour (if l=2) of the father
                     k0=-1
                     do k=1,ipartners(0)
c Loop over previously-found colour/anticolour partners of father
                        if(ipartners(k).eq.j)then
                           if(found)then
c Safety measure: if this condition is met, it means that there exist
c k1 and k2 such that ipartners(k1)=ipartners(k2). This is thus a bug,
c since ipartners() is the list of possible partners of father, where each
c Born-level particle must appears at most once
                              write(*,*)'Error #1 in set_matrices'
                              write(*,*)i,j,l,k
                              stop
                           endif
                           found=.true.
                           k0=k
                        endif
                     enddo
                     if (.not.found) then
                        ipartners(0)=ipartners(0)+1
                        ipartners(ipartners(0))=j
                        k0=ipartners(0)
                     endif
c At this point, k0 is the k0^th colour/anticolour partner of father.
c Therefore, ipartners(k0)=j
                     if(k0.le.0.or.ipartners(k0).ne.j)then
                        write(*,*)'Error #2 in set_matrices'
                        write(*,*)i,j,l,k0,ipartners(k0)
                        stop
                     endif
                     spec_case=l.eq.2 .and. colorflow(k0,0).ge.1 .and.
     &                    colorflow(k0,colorflow(k0,0)).eq.i 
                     if (.not.spec_case)then
c Increase by one the number of colour flows in which the father is
c (anti)colour-connected with its k0^th partner (according to the
c list defined by ipartners)
                        colorflow(k0,0)=colorflow(k0,0)+1
c Store the label of the colour flow thus found
                        colorflow(k0,colorflow(k0,0))=i
                     elseif (spec_case)then
c Special case: father and ipartners(k0) are both gluons, connected
c by colour AND anticolour: the number of colour flows was overcounted
c by one unit, so decrease it
                         if( notagluon .or.
     &                       ICOLUP(i1(1),j,i).eq.0 .or.
     &                       ICOLUP(i1(2),j,i).eq.0 )then
                            write(*,*)'Error #3 in set_matrices'
                            write(*,*)i,j,l,k0,i1(1),i1(2)
                            stop
                         endif
                         colorflow(k0,colorflow(k0,0))=i
                         isspecial(i)=.true.
                     endif
                  endif
               enddo
            endif
         enddo
         if( ((nglu+nsngl).gt.(nexternal-2)) .or.
     #       (isspecial(i).and.(nglu+nsngl).ne.(nexternal-2)) )then
           write(*,*)'Error #4 in set_matrices'
           write(*,*)isspecial(i),nglu,nsngl
           stop
         endif
      enddo
      call check_mc_matrices
      return
      end



      subroutine check_mc_matrices
      implicit none
      include "nexternal.inc"
      include "born_nhel.inc"
c      include "fks.inc"
      integer fks_j_from_i(nexternal,0:nexternal)
     &     ,particle_type(nexternal),pdg_type(nexternal)
      common /c_fks_inc/fks_j_from_i,particle_type,pdg_type
      integer ipartners(0:nexternal-1),colorflow(nexternal-1,0:max_bcol)
      common /MC_info/ ipartners,colorflow
      integer i,j,ipart,iflow,ntot,ithere(1000)
      integer fksfather
      common/cfksfather/fksfather
      logical notagluon
      common/cnotagluon/notagluon
      logical isspecial(max_bcol)
      common/cisspecial/isspecial
      logical is_leading_cflow(max_bcol)
      integer num_leading_cflows
      common/c_leading_cflows/is_leading_cflow,num_leading_cflows
c
      if(ipartners(0).gt.nexternal-1)then
        write(*,*)'Error #1 in check_mc_matrices',ipartners(0)
        stop
      endif
c
      do i=1,ipartners(0)
        ipart=ipartners(i)
        if( ipart.eq.fksfather .or.
     #      ipart.le.0 .or. ipart.gt.nexternal-1 .or.
     #      ( abs(particle_type(ipart)).ne.3 .and.
     #        particle_type(ipart).ne.8 ) )then
          write(*,*)'Error #2 in check_mc_matrices',i,ipart,
     #particle_type(ipart)
          stop
        endif
      enddo
c
      do i=1,nexternal-1
        ithere(i)=1
      enddo
      do i=1,ipartners(0)
        ipart=ipartners(i)
        ithere(ipart)=ithere(ipart)-1
        if(ithere(ipart).lt.0)then
          write(*,*)'Error #3 in check_mc_matrices',i,ipart
          stop
        endif
      enddo
c
c ntot is the total number of colour plus anticolour partners of father
      ntot=0
      do i=1,ipartners(0)
        ntot=ntot+colorflow(i,0)
c
        if( colorflow(i,0).le.0 .or.
     #      colorflow(i,0).gt.max_bcol )then
          write(*,*)'Error #4 in check_mc_matrices',i,colorflow(i,0)
          stop
        endif
c
        do j=1,max_bcol
          ithere(j)=1
        enddo
        do j=1,colorflow(i,0)
          iflow=colorflow(i,j)
          ithere(iflow)=ithere(iflow)-1
          if(ithere(iflow).lt.0)then
            write(*,*)'Error #5 in check_mc_matrices',i,j,iflow
            stop
          endif
        enddo
c
      enddo
c
      if( (notagluon.and.ntot.ne.num_leading_cflows) .or.
     #    ( (.not.notagluon).and.
     #      ( (.not.isspecial(1)).and.ntot.ne.(2*num_leading_cflows) .or.
     #        (isspecial(1).and.ntot.ne.num_leading_cflows) ) ) )then
         write(*,*)'Error #6 in check_mc_matrices',
     #     notagluon,ntot,num_leading_cflows,max_bcol
         stop
      endif
c
      if(num_leading_cflows.gt.max_bcol)then
         write(*,*)'Error #7 in check_mc_matrices',
     #     num_leading_cflows,max_bcol
         stop
      endif
c
      return
      end



      subroutine xmcsubt_wrap(pp,xi_i_fks,y_ij_fks,wgt)
      implicit none
      include "nexternal.inc"
      include 'madfks_mcatnlo.inc'
      double precision pp(0:3,nexternal),wgt
      double precision xi_i_fks,y_ij_fks
      double precision xmc,xrealme,gfactsf,gfactcl,probne
      double precision xmcxsec(nexternal),z(nexternal)
      integer nofpartners
      logical lzone(nexternal),flagmc

      include "born_nhel.inc"
      double precision xmcxsec2(max_bcol)
      integer npartner,cflows
      integer ipartners(0:nexternal-1),colorflow(nexternal-1,0:max_bcol)
      common /MC_info/ ipartners,colorflow
      logical first_MCcnt_call,is_pt_hard
      common/cMCcall/first_MCcnt_call,is_pt_hard

      double precision xkern(2),xkernazi(2),factor,N_p
      double precision bornbars(max_bcol),bornbarstilde(max_bcol)
      double precision emscwgt(nexternal)
      double precision evnt_wgt
      integer i, j
      double precision mu_r
      double precision pb(0:4,-nexternal+3:2*nexternal-3)
      double precision p_read(0:4,2*nexternal-3), wgt_read
      integer npart
      double precision MCsec(nexternal-1,max_bcol)
      double precision dummy
      logical isspecial(max_bcol)
      common/cisspecial/isspecial
c
c True MC subtraction term
      first_MCcnt_call=.true.
      is_pt_hard=.false.
      xmcxsec=0d0
      xmcxsec2=0d0
      do cflows=1,max_bcol
         if(is_pt_hard)cycle
         N_p=1d0
         if(isspecial(cflows))N_p=2d0
         do npartner=1,ipartners(0)
            if(is_pt_hard)cycle
            factor=1d0
c
            call xmcsubt(pp,xi_i_fks,y_ij_fks,gfactsf,gfactcl,probne,
     &           nofpartners,lzone,flagmc,z,xkern,xkernazi,emscwgt,
     &           bornbars,bornbarstilde,npartner)
            if(dampMCsubt)factor=emscwgt(npartner)
            if(colorflow(npartner,cflows).eq.0)factor=0d0
            MCsec(npartner,cflows)=0d0
            if(colorflow(npartner,cflows).ne.0)then
               MCsec(npartner,cflows)=factor*
     &         (xkern(1)*N_p*bornbars(colorflow(npartner,cflows))+
     &         xkernazi(1)*N_p*bornbarstilde(colorflow(npartner,cflows)))
            endif
            xmcxsec(npartner)=xmcxsec(npartner)+MCsec(npartner,cflows)
            xmcxsec2(cflows)=xmcxsec2(cflows)+MCsec(npartner,cflows)
         enddo
      enddo
      if(.not.is_pt_hard)call complete_xmcsubt(pp,xmc,lzone,xmcxsec,
     #                                         xmcxsec2,MCsec,probne)
c G-function matrix element, to recover the real soft limit
      call xmcsubtME(pp,xi_i_fks,y_ij_fks,gfactsf,gfactcl,xrealme)

      wgt=xmc+xrealme

      return
      end



      subroutine xmcsubtME(pp,xi_i_fks,y_ij_fks,gfactsf,gfactcl,wgt)
      implicit none
      include "nexternal.inc"
      include "coupl.inc"
      double precision pp(0:3,nexternal),gfactsf,gfactcl,wgt,wgts,wgtc,wgtsc
      double precision xi_i_fks,y_ij_fks

      double precision zero,one
      parameter (zero=0d0)
      parameter (one=1d0)

      integer izero,ione,itwo
      parameter (izero=0)
      parameter (ione=1)
      parameter (itwo=2)

      double precision p1_cnt(0:3,nexternal,-2:2)
      double precision wgt_cnt(-2:2)
      double precision pswgt_cnt(-2:2)
      double precision jac_cnt(-2:2)
      common/counterevnts/p1_cnt,wgt_cnt,pswgt_cnt,jac_cnt

      double precision xi_i_fks_cnt(-2:2)
      common /cxiifkscnt/xi_i_fks_cnt

      integer i_fks,j_fks
      common/fks_indices/i_fks,j_fks

c Particle types (=colour) of i_fks, j_fks and fks_mother
      integer i_type,j_type,m_type
      common/cparticle_types/i_type,j_type,m_type

      double precision pmass(nexternal)
      include "pmass.inc"
c
      if (i_type.eq.8 .and. pmass(i_fks).eq.0d0)then
c i_fks is gluon
         call set_cms_stuff(izero)
         call sreal(p1_cnt(0,1,0),zero,y_ij_fks,wgts)
         call set_cms_stuff(ione)
         call sreal(p1_cnt(0,1,1),xi_i_fks,one,wgtc)
         call set_cms_stuff(itwo)
         call sreal(p1_cnt(0,1,2),zero,one,wgtsc)
         wgt=wgts+(1-gfactcl)*(wgtc-wgtsc)
         wgt=wgt*(1-gfactsf)
      elseif (abs(i_type).eq.3)then
c i_fks is (anti-)quark
         wgt=0d0
      else
         write(*,*) 'FATAL ERROR #1 in xmcsubtME',i_type,i_fks
         stop
      endif
c
      return
      end



c Main routine for MC counterterms. Now to be called inside a loop
c over colour partners
      subroutine xmcsubt(pp,xi_i_fks,y_ij_fks,gfactsf,gfactcl,probne,
     &     nofpartners,lzone,flagmc,z,xkern,xkernazi,emscwgt,
     &     bornbars,bornbarstilde,npartner)
      implicit none
      include "nexternal.inc"
      include "coupl.inc"
      include "born_nhel.inc"
      include "fks_powers.inc"
      include "madfks_mcatnlo.inc"
      include "run.inc"
      include "../../Source/MODEL/input.inc"
      include 'nFKSconfigs.inc'
      integer fks_j_from_i(nexternal,0:nexternal)
     &     ,particle_type(nexternal),pdg_type(nexternal)
      common /c_fks_inc/fks_j_from_i,particle_type,pdg_type

      double precision pp(0:3,nexternal),gfactsf,gfactcl,probne,wgt
      double precision xi_i_fks,y_ij_fks,xm12,xm22
      double precision xmcxsec(nexternal)
      integer nofpartners,cflows
      logical lzone(nexternal),flagmc,limit,non_limit

      double precision emsca_bare,ptresc,rrnd,ref_scale,
     & scalemin,scalemax,wgt1,qMC,emscainv,emscafun
      double precision emscwgt(nexternal),emscav(nexternal)
      double precision emscwgt_a(nexternal,nexternal),
     # emscav_a(nexternal,nexternal)
      double precision emscav_a2(nexternal,nexternal)
      integer jpartner,mpartner
      logical emscasharp
      double precision emscav_tmp(nexternal)
      double precision emscav_tmp_a(nexternal,nexternal)
      double precision emscav_tmp_a2(nexternal,nexternal)
      common/cemscav_tmp/emscav_tmp
      common/cemscav_tmp_a/emscav_tmp_a,emscav_tmp_a2

      double precision shattmp,dot,xkern(2),xkernazi(2),born_red,
     & born_red_tilde
      double precision bornbars(max_bcol), bornbarstilde(max_bcol)

      integer i,j,npartner,ileg,N_p
      double precision tk,uk,q1q,q2q,E0sq(nexternal),x,yi,yj,xij,ap,Q,
     & s,w1,w2,beta,xfact,prefact,kn,knbar,kn0,betae0,betad,betas,gfactazi,
     & gfunction,bogus_probne_fun,
     & z(nexternal),xi(nexternal),xjac(nexternal),ztmp,xitmp,xjactmp,
     & zHW6,xiHW6,xjacHW6_xiztoxy,zHWPP,xiHWPP,xjacHWPP_xiztoxy,zPY6Q,
     & xiPY6Q,xjacPY6Q_xiztoxy,zPY6PT,xiPY6PT,xjacPY6PT_xiztoxy,zPY8,
     & xiPY8,xjacPY8_xiztoxy,wcc
      common/cqMC/qMC

      double precision probne_bog
      common/cprobne_bog/probne_bog

      common/cscaleminmax/xm12,ileg
      double precision veckn_ev,veckbarn_ev,xp0jfks
      common/cgenps_fks/veckn_ev,veckbarn_ev,xp0jfks
      double precision p_born(0:3,nexternal-1)
      common/pborn/p_born
      integer i_fks,j_fks
      common/fks_indices/i_fks,j_fks
      double precision ybst_til_tolab,ybst_til_tocm,sqrtshat,shat
      common/parton_cms_stuff/ybst_til_tolab,ybst_til_tocm,
     #                        sqrtshat,shat

      integer ipartners(0:nexternal-1),colorflow(nexternal-1,0:max_bcol)
      common /MC_info/ ipartners,colorflow
      logical isspecial(max_bcol)
      common/cisspecial/isspecial

      integer fksfather
      common/cfksfather/fksfather

      logical softtest,colltest
      common/sctests/softtest,colltest

      double precision emsca
      common/cemsca/emsca,emsca_bare,emscasharp,scalemin,scalemax

      logical emscasharp_a(nexternal,nexternal)
      double precision emsca_a(nexternal,nexternal),
     #  emsca_bare_a(nexternal,nexternal)
      double precision emsca_bare_a2(nexternal,nexternal)
      double precision scalemin_a(nexternal,nexternal),
     #  scalemax_a(nexternal,nexternal)
      double precision ptresc_a(nexternal,nexternal)
      common/cemsca_a/emsca_a,emsca_bare_a,emsca_bare_a2,
     #                emscasharp_a,scalemin_a,scalemax_a

      double precision ran2,iseed
      external ran2
      logical extra

c Stuff to be written (depending on AddInfoLHE) onto the LHE file
      INTEGER NFKSPROCESS
      COMMON/C_NFKSPROCESS/NFKSPROCESS
      integer iSorH_lhe,ifks_lhe(fks_configs) ,jfks_lhe(fks_configs)
     &     ,fksfather_lhe(fks_configs) ,ipartner_lhe(fks_configs)
      double precision scale1_lhe(fks_configs),scale2_lhe(fks_configs)
      common/cto_LHE1/iSorH_lhe,ifks_lhe,jfks_lhe,
     #                fksfather_lhe,ipartner_lhe
      common/cto_LHE2/scale1_lhe,scale2_lhe

c Radiation hardness needed (pt_hardness) for the theta function
c Should be zero if there are no jets at the Born
      double precision shower_S_scale(fks_configs*2)
     &     ,shower_H_scale(fks_configs*2),ref_H_scale(fks_configs*2)
     &     ,pt_hardness
      common /cshowerscale2/shower_S_scale,shower_H_scale,ref_H_scale
     &     ,pt_hardness

      double precision becl,delta
c alsf and besf are the parameters that control gfunsoft
      double precision alsf,besf
      common/cgfunsfp/alsf,besf
c alazi and beazi are the parameters that control gfunazi
      double precision alazi,beazi
      common/cgfunazi/alazi,beazi

c Particle types (=color) of i_fks, j_fks and fks_mother
      integer i_type,j_type,m_type
      common/cparticle_types/i_type,j_type,m_type

      double precision zero,one,tiny,vtiny,ymin
      parameter (zero=0d0)
      parameter (one=1d0)
      parameter (vtiny=1d-10)
      parameter (ymin=0.9d0)

      double precision pi
      parameter(pi=3.1415926535897932384626433d0)

      double precision vcf,vtf,vca
      parameter (vcf=4d0/3d0)
      parameter (vtf=1d0/2d0)
      parameter (vca=3d0)

      logical first_MCcnt_call,is_pt_hard
      common/cMCcall/first_MCcnt_call,is_pt_hard

      double precision g_ew,charge,qi2,qj2
      double precision pmass(nexternal)
      save
      include "pmass.inc"

c Initialise if first time
      if(.not.first_MCcnt_call)goto 222
      flagmc   = .false.
      ztmp     = 0d0
      xitmp    = 0d0
      xjactmp  = 0d0
      gfactazi = 0d0
      xkern    = 0d0
      xkernazi = 0d0
      born_red = 0d0
      born_red_tilde = 0d0
      kn       = veckn_ev
      knbar    = veckbarn_ev
      kn0      = xp0jfks
      nofpartners = ipartners(0)
      g_ew=sqrt(4.d0*pi/aewm1)
      qi2=charge(pdg_type(i_fks))**2
      qj2=charge(pdg_type(j_fks))**2
      tiny = 1d-4
      if (softtest.or.colltest)tiny = 1d-6
c Logical variables to control the IR limits:
c one can remove any reference to xi_i_fks
      limit = 1-y_ij_fks.lt.tiny .and. xi_i_fks.ge.tiny
      non_limit = xi_i_fks.ge.tiny

c Discard if unphysical kinematics
      if(pp(0,1).le.0d0)return

c Determine invariants, ileg, and MC hardness qMC
      extra=dampMCsubt.or.AddInfoLHE.or.UseSudakov
      call kinematics_driver(xi_i_fks,y_ij_fks,shat,pp,ileg,
     &                       xm12,xm22,tk,uk,q1q,q2q,qMC,extra)
      w1=-q1q+q2q-tk
      w2=-q2q+q1q-uk
      if(extra.and.qMC.lt.0d0)then
         write(*,*)'Error in xmcsubt: qMC=',qMC
         stop
      endif

c Check ileg, and special case for PYTHIA6PT
      if(ileg.lt.0.or.ileg.gt.4)then
         write(*,*)'Error in xmcsubt: ileg=',ileg
         stop
      endif
      if(ileg.gt.2.and.shower_mc.eq.'PYTHIA6PT')then
         write(*,*)'FSR not allowed when matching PY6PT'
         stop
      endif

c New or standard MC@NLO formulation
      probne=bogus_probne_fun(qMC)
      if(.not.UseSudakov)probne=1.d0
      probne_bog=probne

c Call barred Born and assign shower scale
      call get_mbar(pp,y_ij_fks,ileg,bornbars,bornbarstilde)
      call assign_emsca(pp,xi_i_fks,y_ij_fks)
      call assign_emsca_array(pp,xi_i_fks,y_ij_fks)

c Distinguish ISR and FSR
      if(ileg.le.2)then
         delta=min(1d0,deltaI)
         yj=0d0
         yi=y_ij_fks
      elseif(ileg.ge.3)then
         delta=min(1d0,deltaO)
         yj=y_ij_fks
         yi=0d0
      endif
      x=1-xi_i_fks
      s=shat
      xij=2*(1-xm12/s-(1-x))/(2-(1-x)*(1-yj)) 

c G-function parameters 
      gfactsf=gfunction(x,alsf,besf,2d0)
      if(abs(i_type).eq.3)gfactsf=1d0
      becl=-(1d0-ymin)
      gfactcl=gfunction(y_ij_fks,alsf,becl,1d0)
      if(alazi.lt.0d0)gfactazi=1-gfunction(y_ij_fks,-alazi,beazi,delta)
c For processes that have jets at the Born level, we need to include a
c theta-function: The radiation from the shower should always be softer
c than the jets at the Born, hence no need to include the MC counter
c terms when the radiation is hard.
      if(pt_hardness.gt.shower_S_scale(nFKSprocess*2-1))then
         emsca=2d0*sqrt(ebeam(1)*ebeam(2))
         emsca_a=2d0*sqrt(ebeam(1)*ebeam(2))
         is_pt_hard=.true.
         return
      endif

c Shower variables
      if(shower_mc.eq.'HERWIGPP')then
         ztmp=zHWPP(ileg,xm12,xm22,shat,x,yi,yj,tk,uk,q1q,q2q)
         xitmp=xiHWPP(ileg,xm12,xm22,shat,x,yi,yj,tk,uk,q1q,q2q)
         xjactmp=xjacHWPP_xiztoxy(ileg,xm12,xm22,shat,x,yi,yj,tk,uk,q1q,q2q)
      elseif(shower_mc.eq.'PYTHIA6Q')then
         ztmp=zPY6Q(ileg,xm12,xm22,shat,x,yi,yj,tk,uk,q1q,q2q)
         xitmp=xiPY6Q(ileg,xm12,xm22,shat,x,yi,yj,tk,uk,q1q,q2q)
         xjactmp=xjacPY6Q_xiztoxy(ileg,xm12,xm22,shat,x,yi,yj,tk,uk,q1q,q2q)
      elseif(shower_mc.eq.'PYTHIA6PT')then
         ztmp=zPY6PT(ileg,xm12,xm22,shat,x,yi,yj,tk,uk,q1q,q2q)
         xitmp=xiPY6PT(ileg,xm12,xm22,shat,x,yi,yj,tk,uk,q1q,q2q)
         xjactmp=xjacPY6PT_xiztoxy(ileg,xm12,xm22,shat,x,yi,yj,tk,uk,q1q,q2q)
      elseif(shower_mc.eq.'PYTHIA8')then
         ztmp=zPY8(ileg,xm12,xm22,shat,x,yi,yj,tk,uk,q1q,q2q)
         xitmp=xiPY8(ileg,xm12,xm22,shat,x,yi,yj,tk,uk,q1q,q2q)
         xjactmp=xjacPY8_xiztoxy(ileg,xm12,xm22,shat,x,yi,yj,tk,uk,q1q,q2q)
      endif
      
      first_MCcnt_call=.false.
 222  continue
c Main loop over colour partners used to begin here
      E0sq(npartner)=dot(p_born(0,fksfather),p_born(0,ipartners(npartner)))
      if(E0sq(npartner).lt.0d0)then
         write(*,*)'Error in xmcsubt: negative E0sq'
         write(*,*)E0sq(npartner),ileg,npartner
         stop
      endif
      z(npartner)=ztmp
      xi(npartner)=xitmp
      xjac(npartner)=xjactmp
      if(shower_mc.eq.'HERWIG6')then
         z(npartner)=zHW6(ileg,E0sq(npartner),xm12,xm22,shat,
     &                    x,yi,yj,tk,uk,q1q,q2q)
         xi(npartner)=xiHW6(ileg,E0sq(npartner),xm12,xm22,shat,
     &                      x,yi,yj,tk,uk,q1q,q2q)
         xjac(npartner)=xjacHW6_xiztoxy(ileg,E0sq(npartner),xm12,xm22,
     &                                  shat,x,yi,yj,tk,uk,q1q,q2q)
      endif
c Compute dead zones
      call get_dead_zone(ileg,z(npartner),xi(npartner),s,x,yi,
     &  xm12,xm22,w1,w2,qMC,scalemax,ipartners(npartner),fksfather,
     &  lzone(npartner),wcc)
c Compute MC subtraction terms
      if(lzone(npartner))then
         if(.not.flagmc)flagmc=.true.
         if( (ileg.ge.3 .and. (m_type.eq.8.or.m_type.eq.0)) .or.
     &       (ileg.le.2 .and. (j_type.eq.8.or.j_type.eq.0)) )then
            if(i_type.eq.8)then
c g  --> g  g ( icode = 1 )
c go --> go g
               if(ileg.le.2)then
                  N_p=2
                  if(limit)then
                     xkern(1)=(g**2/N_p)*8*vca*(1-x*(1-x))**2/(s*x**2)
                     xkernazi(1)=-(g**2/N_p)*16*vca*(1-x)**2/(s*x**2)
                  elseif(non_limit)then
                     xfact=(1-yi)*(1-x)/x
                     prefact=4/(s*N_p)
                     call AP_reduced(m_type,i_type,one,z(npartner),ap)
                     ap=ap/(1-z(npartner))
                     xkern(1)=prefact*xfact*xjac(npartner)*ap/xi(npartner)
                     call Qterms_reduced_spacelike(m_type,i_type,one,z(npartner),Q)
                     Q=Q/(1-z(npartner))
                     xkernazi(1)=prefact*xfact*xjac(npartner)*Q/xi(npartner)
                  endif
c
               elseif(ileg.eq.3)then
                  N_p=2
                  if(non_limit)then
                     xfact=(2-(1-x)*(1-(kn0/kn)*yj))/kn*knbar*(1-x)*(1-yj)
                     prefact=2/(s*N_p)
                     call AP_reduced_SUSY(j_type,i_type,one,z(npartner),ap)
                     ap=ap/(1-z(npartner))
                     xkern(1)=prefact*xfact*xjac(npartner)*ap/xi(npartner)
                  endif
c
               elseif(ileg.eq.4)then
                  N_p=2
                  if(limit)then
                     xkern(1)=(g**2/N_p)*( 8*vca*
     &                    (s**2*(1-(1-x)*x)-s*(1+x)*xm12+xm12**2)**2 )/
     &                    ( s*(s-xm12)**2*(s*x-xm12)**2 )
                     xkernazi(1)=-(g**2/N_p)*(16*vca*s*(1-x)**2)/((s-xm12)**2)
                  elseif(non_limit)then
                     xfact=(2-(1-x)*(1-yj))/xij*(1-xm12/s)*(1-x)*(1-yj)
                     prefact=2/(s*N_p)
                     call AP_reduced(j_type,i_type,one,z(npartner),ap)
                     ap=ap/(1-z(npartner))
                     xkern(1)=prefact*xfact*xjac(npartner)*ap/xi(npartner)
                     call Qterms_reduced_timelike(j_type,i_type,one,z(npartner),Q)
                     Q=Q/(1-z(npartner))
                     xkernazi(1)=prefact*xfact*xjac(npartner)*Q/xi(npartner)
                  endif
               endif
            elseif(abs(i_type).eq.3)then
c g --> q q~ ( icode = 2 )
c a --> q q~
               if(ileg.le.2)then
                  N_p=1
                  if(limit)then
                     xkern(1)=(g**2/N_p)*4*vtf*(1-x)*((1-x)**2+x**2)/(s*x)
                     xkern(2)=xkern(1)*(g_ew**2/g**2)*(qi2*vca/vtf)
                  elseif(non_limit)then
                     xfact=(1-yi)*(1-x)/x
                     prefact=4/(s*N_p)
                     call AP_reduced(m_type,i_type,one,z(npartner),ap)
                     ap=ap/(1-z(npartner))
                     xkern(1)=prefact*xfact*xjac(npartner)*ap/xi(npartner)
                     xkern(2)=xkern(1)*(g_ew**2/g**2)*(qi2*vca/vtf)
                  endif
c
               elseif(ileg.eq.4)then
                  N_p=2
                  if(limit)then
                     xkern(1)=(g**2/N_p)*( 4*vtf*(1-x)*
     &                     (s**2*(1-2*(1-x)*x)-2*s*x*xm12+xm12**2) )/
     &                     ( (s-xm12)**2*(s*x-xm12) )
                     xkern(2)=xkern(1)*(g_ew**2/g**2)*(qi2*vca/vtf)
                     xkernazi(1)=(g**2/N_p)*(16*vtf*s*(1-x)**2)/((s-xm12)**2)
                     xkernazi(2)=xkernazi(1)*(g_ew**2/g**2)*(qi2*vca/vtf)
                  elseif(non_limit)then
                     xfact=(2-(1-x)*(1-yj))/xij*(1-xm12/s)*(1-x)*(1-yj)
                     prefact=2/(s*N_p)
                     call AP_reduced(j_type,i_type,one,z(npartner),ap)
                     ap=ap/(1-z(npartner))
                     xkern(1)=prefact*xfact*xjac(npartner)*ap/xi(npartner)
                     xkern(2)=xkern(1)*(g_ew**2/g**2)*(qi2*vca/vtf)
                     call Qterms_reduced_timelike(j_type,i_type,one,z(npartner),Q)
                     Q=Q/(1-z(npartner))
                     xkernazi(1)=prefact*xfact*xjac(npartner)*Q/xi(npartner)
                     xkernazi(2)=xkernazi(1)*(g_ew**2/g**2)*(qi2*vca/vtf)
                  endif
               endif
            else
               write(*,*)'Error 1 in xmcsubt: unknown particle type'
               write(*,*)i_type
               stop
            endif
         elseif( (ileg.ge.3 .and. abs(m_type).eq.3) .or.
     &           (ileg.le.2 .and. abs(j_type).eq.3) )then
            if(abs(i_type).eq.3)then
c q --> g q ( icode = 3 )
c a --> a q
               if(ileg.le.2)then
                  N_p=2
                  if(limit)then
                     xkern(1)=(g**2/N_p)*4*vcf*(1-x)*((1-x)**2+1)/(s*x**2)
                     xkern(2)=xkern(1)*(g_ew**2/g**2)*(qi2/vcf)
                     xkernazi(1)=-(g**2/N_p)*16*vcf*(1-x)**2/(s*x**2)
                     xkernazi(2)=xkernazi(1)*(g_ew**2/g**2)*(qi2/vcf)
                  elseif(non_limit)then
                     xfact=(1-yi)*(1-x)/x
                     prefact=4/(s*N_p)
                     call AP_reduced(m_type,i_type,one,z(npartner),ap)
                     ap=ap/(1-z(npartner))
                     xkern(1)=prefact*xfact*xjac(npartner)*ap/xi(npartner)
                     xkern(2)=xkern(1)*(g_ew**2/g**2)*(qi2/vcf)
                     call Qterms_reduced_spacelike(m_type,i_type,one,z(npartner),Q)
                     Q=Q/(1-z(npartner))
                     xkernazi(1)=prefact*xfact*xjac(npartner)*Q/xi(npartner)
                     xkernazi(2)=xkernazi(1)*(g_ew**2/g**2)*(qi2/vcf)
                  endif
c
               elseif(ileg.eq.3)then
                  N_p=1
                  if(non_limit)then
                     xfact=(2-(1-x)*(1-(kn0/kn)*yj))/kn*knbar*(1-x)*(1-yj)
                     prefact=2/(s*N_p)
                     call AP_reduced(j_type,i_type,one,z(npartner),ap)
                     ap=ap/(1-z(npartner))
                     xkern(1)=prefact*xfact*xjac(npartner)*ap/xi(npartner)
                     xkern(2)=xkern(1)*(g_ew**2/g**2)*(qi2/vcf)
                  endif
c
               elseif(ileg.eq.4)then
                  N_p=1
                  if(limit)then
                     xkern(1)=(g**2/N_p)*
     &                    ( 4*vcf*(1-x)*(s**2*(1-x)**2+(s-xm12)**2) )/
     &                    ( (s-xm12)*(s*x-xm12)**2 )
                     xkern(2)=xkern(1)*(g_ew**2/g**2)*(qi2/vcf)
                  elseif(non_limit)then
                     xfact=(2-(1-x)*(1-yj))/xij*(1-xm12/s)*(1-x)*(1-yj)
                     prefact=2/(s*N_p)
                     call AP_reduced(j_type,i_type,one,z(npartner),ap)
                     ap=ap/(1-z(npartner))
                     xkern(1)=prefact*xfact*xjac(npartner)*ap/xi(npartner)
                     xkern(2)=xkern(1)*(g_ew**2/g**2)*(qi2/vcf)
                  endif
               endif
            elseif(i_type.eq.8)then
c q  --> q  g ( icode = 4 )
c sq --> sq g
               if(ileg.le.2)then
                  N_p=1
                  if(limit)then
                     xkern(1)=(g**2/N_p)*4*vcf*(1+x**2)/(s*x)
                  elseif(non_limit)then
                     xfact=(1-yi)*(1-x)/x
                     prefact=4/(s*N_p)
                     call AP_reduced(m_type,i_type,one,z(npartner),ap)
                     ap=ap/(1-z(npartner))
                     xkern(1)=prefact*xfact*xjac(npartner)*ap/xi(npartner)
                  endif
c
               elseif(ileg.eq.3)then
                  N_p=1
                  if(non_limit)then
                     xfact=(2-(1-x)*(1-(kn0/kn)*yj))/kn*knbar*(1-x)*(1-yj)
                     prefact=2/(s*N_p)
                     if(abs(PDG_type(j_fks)).le.6)then
                        if(shower_mc.ne.'HERWIGPP')
     &                  call AP_reduced(j_type,i_type,one,z(npartner),ap)
                        if(shower_mc.eq.'HERWIGPP')
     &                  call AP_reduced_massive(j_type,i_type,one,z(npartner),
     &                                              xi(npartner),xm12,ap)
                     else
                        call AP_reduced_SUSY(j_type,i_type,one,z(npartner),ap)
                     endif
                     ap=ap/(1-z(npartner))
                     xkern(1)=prefact*xfact*xjac(npartner)*ap/xi(npartner)
                  endif
c
               elseif(ileg.eq.4)then
                  N_p=1
                  if(limit)then
                     xkern(1)=(g**2/N_p)*4*vcf*
     &                     ( s**2*(1+x**2)-2*xm12*(s*(1+x)-xm12) )/
     &                     ( s*(s-xm12)*(s*x-xm12) )
                  elseif(non_limit)then
                     xfact=(2-(1-x)*(1-yj))/xij*(1-xm12/s)*(1-x)*(1-yj)
                     prefact=2/(s*N_p)
                     call AP_reduced(j_type,i_type,one,z(npartner),ap)
                     ap=ap/(1-z(npartner))
                     xkern(1)=prefact*xfact*xjac(npartner)*ap/xi(npartner)
                  endif
               endif
            elseif(i_type.eq.0)then
c q  --> q  a ( icode = 4 )
c sq --> sq a
               if(ileg.le.2)then
                  N_p=1
                  if(limit)then
                     xkern(2)=(g_ew**2/N_p)*4*qj2*(1+x**2)/(s*x)
                  else
                     xfact=(1-yi)*(1-x)/x
                     prefact=4/(s*N_p)
                     call AP_reduced(m_type,i_type,one,z(npartner),ap)
                     ap=ap/(1-z(npartner))
                     xkern(2)=prefact*xfact*xjac(npartner)*ap/xi(npartner)
                     xkern(2)=xkern(2)*(g_ew**2/g**2)*(qj2/vcf)
                  endif
c
               elseif(ileg.eq.3)then
                  N_p=1
                  if(non_limit)then
                     xfact=(2-(1-x)*(1-(kn0/kn)*yj))/kn*knbar*(1-x)*(1-yj)
                     prefact=2/(s*N_p)
                     if(abs(PDG_type(j_fks)).le.6)then
                        if(shower_mc.ne.'HERWIGPP')
     &                  call AP_reduced(j_type,i_type,one,z(npartner),ap)
                        if(shower_mc.eq.'HERWIGPP')
     &                    call AP_reduced_massive(j_type,i_type,one,
     &                      z(npartner),xi(npartner),xm12,ap)
                     else
                        call AP_reduced_SUSY(j_type,i_type,one,z(npartner),ap)
                     endif
                     ap=ap/(1-z(npartner))
                     xkern(2)=prefact*xfact*xjac(npartner)*ap/xi(npartner)
                     xkern(2)=xkern(2)*(g_ew**2/g**2)*(qj2/vcf)
                  endif
c
               elseif(ileg.eq.4)then
                  N_p=1
                  if(limit)then
                     xkern(2)=(g_ew**2/N_p)*4*qj2*
     &                    ( s**2*(1+x**2)-2*xm12*(s*(1+x)-xm12) )/
     &                    ( s*(s-xm12)*(s*x-xm12) )
                  else
                     xfact=(2-(1-x)*(1-yj))/xij*(1-xm12/s)*(1-x)*(1-yj)
                     prefact=2/(s*N_p)
                     call AP_reduced(j_type,i_type,one,z(npartner),ap)
                     ap=ap/(1-z(npartner))
                     xkern(2)=prefact*xfact*xjac(npartner)*ap/xi(npartner)
                     xkern(2)=xkern(2)*(g_ew**2/g**2)*(qj2/vcf)
                  endif
               endif
            else
               write(*,*)'Error 2 in xmcsubt: unknown particle type'
               write(*,*)i_type
               stop
            endif
         else
            write(*,*)'Error 3 in xmcsubt: unknown particle type'
            write(*,*)j_type,i_type
            stop
         endif

      else
c Dead zone
        xkern=0d0
        xkernazi=0d0
      endif
c
      xkern=xkern*gfactsf*wcc
      xkernazi=xkernazi*gfactazi*gfactsf*wcc

c Emsca stuff
      if(dampMCsubt)then
         if(emscasharp)then
            if(qMC.le.scalemax)then
               emscwgt(npartner)=1d0
               emscav(npartner)=emsca_bare
            else
               emscwgt(npartner)=0d0
               emscav(npartner)=scalemax
            endif
         else
            ptresc=(qMC-scalemin)/(scalemax-scalemin)
            if(ptresc.le.0d0)then
               emscwgt(npartner)=1d0
               emscav(npartner)=emsca_bare
            elseif(ptresc.lt.1d0)then 
               emscwgt(npartner)=1-emscafun(ptresc,one)
               emscav(npartner)=emsca_bare
            else
               emscwgt(npartner)=0d0
               emscav(npartner)=scalemax
            endif
         endif
      endif
      emscav_tmp(npartner)=emscav(npartner)
c Emsca stuff
      do i=1,nexternal-2
         do j=i+1,nexternal-1 
            if(dampMCsubt)then
               if(emscasharp_a(i,j))then
                  if(qMC.le.scalemax_a(i,j))then
                     emscwgt_a(i,j)=1d0
                     emscav_a(i,j)=emsca_bare_a(i,j)
                     emscav_a2(i,j)=emsca_bare_a2(i,j)
                  else
                     emscwgt_a(i,j)=0d0
                     emscav_a(i,j)=scalemax_a(i,j)
                     emscav_a2(i,j)=scalemax_a(i,j)
                  endif
               else
                  ptresc_a(i,j)=(qMC-scalemin_a(i,j))/
     &                          (scalemax_a(i,j)-scalemin_a(i,j))
                  if(ptresc_a(i,j).le.0d0)then
                     emscwgt_a(i,j)=1d0
                     emscav_a(i,j)=emsca_bare_a(i,j)
                     emscav_a2(i,j)=emsca_bare_a2(i,j)
                  elseif(ptresc_a(i,j).lt.1d0)then 
                     emscwgt_a(i,j)=1-emscafun(ptresc_a(i,j),one)
                     emscav_a(i,j)=emsca_bare_a(i,j)
                     emscav_a2(i,j)=emsca_bare_a2(i,j)
                  else
                     emscwgt_a(i,j)=0d0
                     emscav_a(i,j)=scalemax_a(i,j)
                     emscav_a2(i,j)=scalemax_a(i,j)
                  endif
               endif
            endif
            emscav_tmp_a(i,j)=emscav_a(i,j)
            emscav_tmp_a2(i,j)=emscav_a2(i,j)
c
            ptresc_a(j,i)=ptresc_a(i,j)
            emscwgt_a(j,i)=emscwgt_a(i,j)
            emscav_a(j,i)=emscav_a(i,j)
            emscav_a2(j,i)=emscav_a2(i,j)
            emscav_tmp_a(j,i)=emscav_tmp_a(i,j)
            emscav_tmp_a2(j,i)=emscav_tmp_a2(i,j)
         enddo
      enddo
c Main loop over colour partners used to end here
      return
      end


c Finalises the MC counterterm computations performed in xmcsubt(),
c fills arrays relevant to shower scales, and computes Delta
      subroutine complete_xmcsubt(p,wgt,lzone,xmcxsec,xmcxsec2,MCsec,probne)
      implicit none
      include "born_nhel.inc"
      include 'nFKSconfigs.inc'
      include 'nexternal.inc'
      include 'madfks_mcatnlo.inc'
      include 'run.inc'

      integer i_fks,j_fks
      common/fks_indices/i_fks,j_fks

      double precision emsca_bare,ptresc,rrnd,ref_scale,
     & scalemin,scalemax,wgt11,qMC,emscainv,emscafun
      double precision emscwgt(nexternal),emscav(nexternal)
      double precision emscwgt_a(nexternal,nexternal),
     & emscav_a(nexternal,nexternal)
      double precision emscav_a2(nexternal,nexternal)
      integer jpartner,mpartner,cflows,jflow
      common/c_colour_flow/jflow
      logical emscasharp

      double precision emsca
      common/cemsca/emsca,emsca_bare,emscasharp,scalemin,scalemax

      double precision emsca_a(nexternal,nexternal)
      double precision emsca_bare_a(nexternal,nexternal),
     & emsca_bare_a2(nexternal,nexternal)
      logical emscasharp_a(nexternal,nexternal)
      double precision scalemin_a(nexternal,nexternal),
     & scalemax_a(nexternal,nexternal)
      common/cemsca_a/emsca_a,emsca_bare_a,emsca_bare_a2,
     #                emscasharp_a,scalemin_a,scalemax_a
      common/cqMC/qMC

      integer ipartners(0:nexternal-1),colorflow(nexternal-1,0:max_bcol)
      common /MC_info/ ipartners,colorflow

      integer fksfather
      common/cfksfather/fksfather

      double precision ran2,iseed
      external ran2
      logical extra

c Controls assignments of scales in H events in LHE file.
c Set iHscale=0 for scale=target_scale
c     iHscale=1 for scale=dipole_mass
      integer iHscale,jbar,ifksscl(2)
      parameter (iHscale=0)
      double precision dipole_mass,fksscales(3)
      external dipole_mass

c Maps real labels onto Born labels, by excluding i_fks
c  1<=iRtoB(k)<=nexternal-1,  1<=k<=nexternal
      integer iRtoB(nexternal)
c Maps Born labels onto real labels, by excluding i_fks
c  1<=iBtoR(k)<=nexternal,  1<=k<=nexternal-1
      integer iBtoR(nexternal-1)

c Stuff to be written (depending on AddInfoLHE) onto the LHE file
      INTEGER NFKSPROCESS
      COMMON/C_NFKSPROCESS/NFKSPROCESS
      integer iSorH_lhe,ifks_lhe(fks_configs) ,jfks_lhe(fks_configs)
     &     ,fksfather_lhe(fks_configs) ,ipartner_lhe(fks_configs)
      double precision scale1_lhe(fks_configs),scale2_lhe(fks_configs)
      common/cto_LHE1/iSorH_lhe,ifks_lhe,jfks_lhe,
     #                fksfather_lhe,ipartner_lhe
      common/cto_LHE2/scale1_lhe,scale2_lhe

      integer npartner
      double precision emscav_tmp(nexternal)
      double precision emscav_tmp_a(nexternal,nexternal)
      double precision emscav_tmp_a2(nexternal,nexternal)
      common/cemscav_tmp/emscav_tmp
      common/cemscav_tmp_a/emscav_tmp_a,emscav_tmp_a2

      double precision xmcxsec(nexternal),xmcxsec2(max_bcol),probne,wgt,wgt2
      logical lzone(nexternal)

      integer i,j,k
      double precision tiny
      parameter(tiny=1d-7)

      double precision p(0:3,nexternal)
      double precision xkern(2),xkernazi(2),factor
      double precision bornbars(max_bcol),bornbarstilde(max_bcol)
      double precision MCsec(nexternal-1,max_bcol)
      include "genps.inc"
      integer idup(nexternal-1,maxproc)
      integer mothup(2,nexternal-1,maxproc)
      integer icolup(2,nexternal-1,max_bcol)
      integer idup_s(nexternal-1)
      integer mothup_s(2,nexternal-1)
      integer icolup_s(2,nexternal-1)
      integer idup_h(nexternal)
      integer mothup_h(2,nexternal)
      integer icolup_h(2,nexternal)
      integer spinup_local(nexternal)
      integer istup_local(nexternal)
      double precision wgt_sudakov
      double precision scales(0:99)
      common /colour_connections/ icolup_s,icolup_h

      integer itmp,jtmp
      double precision scltarget,sclstart,sudpdffact
      common/cscalprob/scltarget,sclstart,sudpdffact

c To access Pythia8 control variables
      include 'pythia8_control.inc'
      include "born_leshouche.inc"
      integer jpart(7,-nexternal+3:2*nexternal-3),lc,iflow
      logical firsttime1
      data firsttime1 /.true./
      include 'leshouche_decl.inc'
      save mothup_d, icolup_d, niprocs_d

      logical         Hevents
      common/SHevents/Hevents
      integer nexternal_now
      double precision sumMCsec(max_bcol)
c SCALUP_tmp_S = m_ij scales that determine S-event scales written onto LHE
      double precision SCALUP_tmp_S(nexternal,nexternal)
c SCALUP_tmp_S2 = m_ij starting scales for Delta
      double precision SCALUP_tmp_S2(nexternal,nexternal)
c SCALUP_tmp_H = t_ij scales that determine H-event scales written onto LHE
      double precision SCALUP_tmp_H(nexternal,nexternal)
c SCALUP_tmp_H2 = t_ij target scales for Delta
      double precision SCALUP_tmp_H2(nexternal,nexternal)
      common/c_SCALUP_tmp/SCALUP_tmp_S,SCALUP_tmp_H

c Lower and upper limits of fitted st and xm ranges.
c Require one prior call to pysudakov() to be set,
c here done in the firsttime1 clause
      real*8 cstlow,cstupp,cxmlow,cxmupp
      common/cstxmbds/cstlow,cstupp,cxmlow,cxmupp

c Set Delta(pt,..)=0 for pt<smallptlow, and interpolate
c between 0 and Delta(smallptupp,..) for smallptlow<pt<smallptupp
c For things to work properly, one must have:
c               cstlow <= smallptupp
      real*8 smallptlow,smallptupp,get_to_zero
      parameter (smallptlow=0.5d0)
      parameter (smallptupp=1.01d0)

      integer iii,jjj,LP
      double precision xscales(0:99,0:99)
      double precision xmasses(0:99,0:99)
      double precision xscales2(0:99,0:99)
      double precision xmasses2(0:99,0:99)
      logical*1 dzones(0:99,0:99)
      logical*1 dzones2(0:99,0:99)

      integer id,type,icount
cSF ARE noemProb AND mDipole USEFUL?
      double precision noemProb, startingScale, stoppingScale, mDipole
      double precision mcmass(21)
      double precision pysudakov,deltanum,deltaden,deltarat
      integer nG_S,nQ_S,i_dipole_counter,isudtype
      integer i_dipole_dead_counter
c
      integer fks_j_from_i(nexternal,0:nexternal)
     &     ,particle_type(nexternal),pdg_type(nexternal)
      common /c_fks_inc/fks_j_from_i,particle_type,pdg_type

      double precision xbjrk_ev(2),xbjrk_cnt(2,-2:2)
      common/cbjorkenx/xbjrk_ev,xbjrk_cnt

      double precision pdg2pdf,pdffnum,pdffden
      external pdg2pdf
c
      LOGICAL  IS_A_J(NEXTERNAL),IS_A_LP(NEXTERNAL),IS_A_LM(NEXTERNAL)
      LOGICAL  IS_A_PH(NEXTERNAL)
      COMMON /TO_SPECISA/IS_A_J,IS_A_LP,IS_A_LM,IS_A_PH
      integer idIn1, idIn2
      integer idOut(0:9)
      double precision tBefore,tAfter
      logical do_time_profiling
      parameter (do_time_profiling=.false.)
      double precision masses_to_MC(0:25)
      double precision pi
      parameter(pi=3.1415926535897932384626433d0)
      logical are_col_conn_S(nexternal-1,nexternal-1)
      logical are_col_conn_H(nexternal,nexternal)
      integer get_mass_from_id
      external get_mass_from_id
      double precision p_born(0:3,nexternal-1)
      common/pborn/p_born
      double complex cdummy(2)
c Jamp amplitudes of the Born (to be filled with a call the sborn())
      double Precision amp2(ngraphs), jamp2(0:ncolor)
      common/to_amps/  amp2,       jamp2
c
      mcmass=0d0
      masses_to_MC=0d0
      include 'MCmasses_PYTHIA8.inc'
c
      do i=1,2
        istup_local(i) = -1
      enddo
      do i=3,nexternal
        istup_local(i) = 1
      enddo
      do i=1,nexternal
        spinup_local(i) = -9
      enddo
      pythia_cmd_file=''

c Input check
      do npartner=1,ipartners(0)
         if(xmcxsec(npartner).lt.0d0)then
            write(*,*)'Fatal error 1 in complete_xmcsubt'
            write(*,*)npartner,xmcxsec(npartner)
            stop
         endif
      enddo
      do cflows=1,max_bcol
         if(xmcxsec2(cflows).lt.0d0)then
            write(*,*)'Fatal error 2 in complete_xmcsubt'
            write(*,*)cflows,xmcxsec2(cflows)
            stop
         endif
      enddo

c Compute MC cross section
      wgt=0d0
      wgt2=0d0
      do i=1,max_bcol
         sumMCsec(i)=0d0
      enddo
      do npartner=1,ipartners(0)
         wgt=wgt+xmcxsec(npartner)
      enddo
      do cflows=1,max_bcol
         wgt2=wgt2+xmcxsec2(cflows)
      enddo
      do cflows=1,max_bcol
         do npartner=1,ipartners(0)
            sumMCsec(cflows)=sumMCsec(cflows)+MCsec(npartner,cflows)
         enddo
      enddo
c
      if( (abs(wgt).gt.1.d-10 .and.abs(wgt-wgt2)/abs(wgt).gt.tiny) .or.
     &    (abs(wgt).le.1.d-10 .and.abs(wgt-wgt2).gt.tiny) )then
         write(*,*)'Fatal error 3 in complete_xmcsubt'
         write(*,*)wgt,wgt2
         stop
      endif
      do cflows=1,max_bcol
         if( (abs(sumMCsec(cflows)).gt.1.d-10 .and.
     &        abs(sumMCsec(cflows)-xmcxsec2(cflows))/
     &        abs(sumMCsec(cflows)).gt.tiny) .or.
     &       (abs(sumMCsec(cflows)).le.1.d-10 .and.
     &        abs(sumMCsec(cflows)-xmcxsec2(cflows)).gt.tiny) )then
            write(*,*)'Fatal error 3 in complete_xmcsubt'
            write(*,*)sumMCsec(cflows),xmcxsec2(cflows)
            stop
         endif
      enddo

c Assign flow on statistical basis
      if (wgt2.gt.0d0) then
         ! use born-bars times kernels
         rrnd=ran2()
         wgt11=0d0
         jflow=0
         cflows=0
         do while(jflow.eq.0.and.cflows.lt.max_bcol)
            cflows=cflows+1
            wgt11=wgt11+xmcxsec2(cflows)
            if(wgt11.ge.rrnd*wgt2)jflow=cflows
         enddo
         if(jflow.eq.0)then
            write(*,*)'Error in xmcsubt: flow unweighting failed'
            stop
         endif
      else
         ! use the born-bars
         call sborn(p_born,cdummy)
         wgt11=0.d0
         do i=1,max_bcol
            wgt11=wgt11+jamp2(i)
         enddo
         wgt2=ran2()*wgt11
         jflow=0
         wgt11=0d0
         do while (wgt11 .lt. wgt2)
            jflow=jflow+1
            wgt11=wgt11+jamp2(jflow)
         enddo
      endif
         

c Assign emsca (scalar) on statistical basis -- ensure backward compatibility
      if(dampMCsubt.and.wgt.gt.1d-30)then
        rrnd=ran2()
        wgt11=0d0
        jpartner=0
        do npartner=1,ipartners(0)
           if(lzone(npartner).and.jpartner.eq.0)then
              wgt11=wgt11+MCsec(npartner,jflow)
              if(wgt11.ge.rrnd*xmcxsec2(jflow))then
                 jpartner=ipartners(npartner)
                 mpartner=npartner
              endif
           endif
        enddo
        if(jpartner.eq.0)then
           write(*,*)'Error in xmcsubt: emsca unweighting failed'
           stop
        else
           emsca=emscav_tmp(mpartner)
        endif
      endif
      if(dampMCsubt.and.wgt.lt.1d-30)emsca=scalemax

c Additional information for LHE
      if(AddInfoLHE)then
         ifks_lhe(nFKSprocess)=i_fks
         jfks_lhe(nFKSprocess)=j_fks
         fksfather_lhe(nFKSprocess)=fksfather
         if(jpartner.ne.0)then
            ipartner_lhe(nFKSprocess)=jpartner
         else
c min() avoids troubles if ran2()=1
            ipartner_lhe(nFKSprocess)=min(int(ran2()*ipartners(0))+1,
     #                                    ipartners(0) )
            ipartner_lhe(nFKSprocess)=ipartners(ipartner_lhe(nFKSprocess))
         endif
         scale1_lhe(nFKSprocess)=qMC
      endif

      if(dampMCsubt)then
         if(emsca.lt.scalemin)then
            write(*,*)'Error in xmcsubt: emsca too small'
            write(*,*)emsca,jpartner,lzone(jpartner)
            stop
         endif
      endif

c S-event information:
c id's and mothers read from born_leshouche.inc;
c colour configuration read from born_leshouche.inc and jflow 
      do i=1,nexternal-1
        IDUP_S(i)=IDUP(i,1)
        MOTHUP_S(1,i)=MOTHUP(1,i,1)
        MOTHUP_S(2,i)=MOTHUP(2,i,1)
        ICOLUP_S(1,i)=ICOLUP(1,i,jflow)
        ICOLUP_S(2,i)=ICOLUP(2,i,jflow)
      enddo
      are_col_conn_S=.false.
      do i=1,nexternal-1
         do j=1,nexternal-1
            if(i.ne.j)
     &        are_col_conn_S(i,j)=
     &      (ICOLUP_S(1,i).ne.0.and.ICOLUP_S(1,i).eq.ICOLUP_S(1,j)).or.
     &      (ICOLUP_S(1,i).ne.0.and.ICOLUP_S(1,i).eq.ICOLUP_S(2,j)).or.
     &      (ICOLUP_S(2,i).ne.0.and.ICOLUP_S(2,i).eq.ICOLUP_S(1,j)).or.
     &      (ICOLUP_S(2,i).ne.0.and.ICOLUP_S(2,i).eq.ICOLUP_S(2,j))
         enddo
      enddo
c SCALUP_tmp_S* are the m_ij scales, ie the starting scales (as determined
c by the D(mu) function) for extra radiation; they are copies of the
c emscav_tmp_a* arrays, originally filled by xmcsubt(). Only the (i,j) 
c entries associated with a colour line that belongs to jflow have
c meaningful values; the others are set equal to -1.
c SCALUP_tmp_S and SCALUP_tmp_S2 are chosen in exactly the same way, except
c for the random numbers that enter their definitions. The former will
c help determine the S-event shower scales written onto the LHE file, 
c the latter is employed in the computation of Delta
      SCALUP_tmp_S=-1d0
      SCALUP_tmp_S2=-1d0
      do i=1,nexternal-2
         do j=i+1,nexternal-1
            if(are_col_conn_S(i,j))then
               SCALUP_tmp_S(i,j)=emscav_tmp_a(i,j)
               SCALUP_tmp_S(j,i)=emscav_tmp_a(j,i)
               SCALUP_tmp_S2(i,j)=emscav_tmp_a2(i,j)
               SCALUP_tmp_S2(j,i)=emscav_tmp_a2(j,i)
            endif
         enddo
      enddo
c
c H-event information.
c First write ids, mothers and all colours.
cSF NOTE: reconsider how much H-event information is actually needed
cSF by Pythia. For example, the colour flow is used only to reconstruct
cSF the underlying S-event flow, which is already available here, and
cSF thus can be directly passed rather than reconstructed
      if (firsttime1)then
        firsttime1=.false.
        call read_leshouche_info2(idup_d,mothup_d,icolup_d,niprocs_d)
c Fake call for initialisation
        deltanum=pysudakov(1.d2,2.d2,1,1,mcmass)
        if(cstlow.gt.smallptupp)then
          write(*,*)'Error in xmcsubt: cstlow,smallptupp',
     &               cstlow,smallptupp
          stop
        endif
      endif
      do i=1,nexternal
         IDUP_H(i)=IDUP_D(nFKSprocess,i,1)
         MOTHUP_H(1,i)=MOTHUP_D(nFKSprocess,1,i,1)
         MOTHUP_H(2,i)=MOTHUP_D(nFKSprocess,2,i,1)
      enddo
c Fill selected color configuration into jpart array. 
      call fill_icolor_H(jflow,jpart)
      do i=1,nexternal
        ICOLUP_H(1,i)=jpart(4,i)
        ICOLUP_H(2,i)=jpart(5,i)
      enddo
      are_col_conn_H=.false.
      do i=1,nexternal
         do j=1,nexternal
            if(i.ne.j)
     &        are_col_conn_H(i,j)=
     &      (ICOLUP_H(1,i).ne.0.and.ICOLUP_H(1,i).eq.ICOLUP_H(1,j)).or.
     &      (ICOLUP_H(1,i).ne.0.and.ICOLUP_H(1,i).eq.ICOLUP_H(2,j)).or.
     &      (ICOLUP_H(2,i).ne.0.and.ICOLUP_H(2,i).eq.ICOLUP_H(1,j)).or.
     &      (ICOLUP_H(2,i).ne.0.and.ICOLUP_H(2,i).eq.ICOLUP_H(2,j))
         enddo
      enddo
      do i=1,nexternal
        if(i.lt.i_fks)then
          iRtoB(i)=i
          iBtoR(i)=i
        elseif(i.eq.i_fks)then
          iRtoB(i)=-1
          if(i.lt.nexternal)iBtoR(i)=i+1
        elseif(i.gt.i_fks)then
          iRtoB(i)=i-1
          if(i.lt.nexternal)iBtoR(i)=i+1
        endif
      enddo
c
      nexternal_now=nexternal
      call clear_HEPEUP_event()
      call fill_HEPEUP_event_2(p, wgt, nexternal_now, idup_h,
     &       istup_local, mothup_h, icolup_h, spinup_local,
     &       emsca, emscav_tmp, emscav_tmp_a)
      xscales=-1d0
      xmasses=-1d0
      dzones=.true.
      xscales2=-1d0
      xmasses2=-1d0
      dzones2=.true.
      if(do_time_profiling)then
         if (is_pythia_active.eq.0) then
c Fill masses
            do i=7,20
               if(i.le.10.or.i.ge.17)masses_to_MC(i)=-1d0
            enddo
            masses_to_MC(5) =get_mass_from_id(5)
            masses_to_MC(6) =get_mass_from_id(6)
            masses_to_MC(15)=get_mass_from_id(15)
            masses_to_MC(23)=get_mass_from_id(23)
            masses_to_MC(24)=get_mass_from_id(24)
            masses_to_MC(25)=get_mass_from_id(25)
c
            call cpu_time(tBefore)
            idOut=0
            do i=3,nexternal-1
              idOut(i-3) = IDUP_S(i)
              if ( is_a_j(i) ) idOut(i-3)=2212
            enddo
            idIn1 = idup_s(1)
            idIn2 = idup_s(2)
            if ( abs(idIn1) < 10 .or. idIn1 .eq. 21) idIn1=2212
            if ( abs(idIn2) < 10 .or. idIn2 .eq. 21) idIn2=2212
            call dire_init_default(idIn1, idIn2, idOut, masses_to_MC)
            call cpu_time(tAfter)
            write(*,*)'time elapsed in dire_init_default',tAfter-tBefore
         endif
         call cpu_time(tBefore)
         call dire_setevent()
         call cpu_time(tAfter)
         write(*,*)'time elapsed in dire_setevent',tAfter-tBefore
         call cpu_time(tBefore)
         call dire_next()
         call cpu_time(tAfter)
         write(*,*)'time elapsed in dire_next',tAfter-tBefore
         call cpu_time(tBefore)
         call dire_get_stopping_info(xscales,xmasses)
         call cpu_time(tAfter)
         write(*,*)'time elapsed in dire_get_stopping_info',tAfter-tBefore
         call cpu_time(tBefore)
         call dire_get_dead_zones(dzones)
         call cpu_time(tAfter)
         write(*,*)'time elapsed in dire_get_dead_zones',tAfter-tBefore
         write(*,*)
      else
         if (is_pythia_active.eq.0) then
c Fill masses
            do i=7,20
               if(i.le.10.or.i.ge.17)masses_to_MC(i)=-1d0
            enddo
            masses_to_MC(5) =get_mass_from_id(5)
            masses_to_MC(6) =get_mass_from_id(6)
            masses_to_MC(15)=get_mass_from_id(15)
            masses_to_MC(23)=get_mass_from_id(23)
            masses_to_MC(24)=get_mass_from_id(24)
            masses_to_MC(25)=get_mass_from_id(25)
c
            idOut=0
            do i=3,nexternal-1
              idOut(i-3) = IDUP_S(i)
              if ( is_a_j(i) ) idOut(i-3)=2212
            enddo
            idIn1 = idup_s(1)
            idIn2 = idup_s(2)
            if ( abs(idIn1) < 10 .or. idIn1 .eq. 21) idIn1=2212
            if ( abs(idIn2) < 10 .or. idIn2 .eq. 21) idIn2=2212
            call dire_init_default(idIn1, idIn2, idOut, masses_to_MC)
         endif
         call dire_setevent()
         call dire_next()
         call dire_get_stopping_info(xscales,xmasses)
         call dire_get_dead_zones(dzones)
         call dire_clear()
      endif
c After the calls above, we have
c   xscales(i,j)=t_ij
c with t_ij == scale(Pythia)_{emitter,recoiler}, and the particle being
c emitted equal to the FKS parton. Although both emitter and recoiler
c are Born-level quantities, their labellings follow the real-process
c conventions. Thus, in the matrix xscales(i,j) one has 1<=i,j<=nexternal, 
c with xscales(i_fks,*)=xscales(*,i_fks)=-1.
c The same labelling conventions apply to xmasses(i,j) (which is the
c dipole mass associated with the colour line that connects i and j)
c and dzones(i,j) (which is the dead zone relevant to the emission
c from parton i colour-connected with recoiler j).
c
c Since any the pair of indices (i,j) associated with sensible
c entries in the arrays returned by Pythia is in one-to-one correspondence
c with Born-level quantities, it is convenient to define relabelled
c copies of such arrays (which we call xscales2, xmasses2, and dzones2),
c for which 1<=i,j<=nexternal-1
c
c By construction, t_ij are the target scales. For notational consistency
c with the case of SCALUP_tmp_S2, a copy of xscales2 is created and called
c SCALUP_tmp_H2, meant to be used in the computation of Delta. 
      SCALUP_tmp_H2=-1d0
      do i=1,nexternal
         if(i.eq.i_fks)cycle
         do j=1,nexternal
            if(j.eq.i_fks)cycle
            xscales2(iRtoB(i),iRtoB(j))=xscales(i,j)
            xmasses2(iRtoB(i),iRtoB(j))=xmasses(i,j)
            dzones2(iRtoB(i),iRtoB(j))=dzones(i,j)
            scalup_tmp_H2(iRtoB(i),iRtoB(j))=
     &        xscales2(iRtoB(i),iRtoB(j))
         enddo
      enddo
c Checks
      do i=1,nexternal-1
         do j=1,nexternal-1
            if((xscales2(i,j).ne.-1d0.and.xmasses2(i,j).eq.-1d0).or.
     &         (xscales2(i,j).eq.-1d0.and.xmasses2(i,j).ne.-1d0))then
               write(*,*)'Error in xmcsubt: xscales, xmasses',
     &                   i,j,xscales2(i,j),xmasses2(i,j)
               stop
            endif
         enddo
      enddo
c
c Scales written onto the LHE file for H events.
c For notational consistency with the case of SCALUP_tmp_S, the array that
c contains these scales is called SCALUP_tmp_H. If the scales are identified
c with the target scales, SCALUP_tmp_H is a copy of xscales, but this is
c only one of the possible options. Furthermore, at this point of the code
c xscales is correctly labelled with 1<=i,j<=nexternal, but lacks the entries 
c relevant to i and j equal to i_fks; these will be provided later
      SCALUP_tmp_H=-1d0
      do i=1,nexternal
         icount=0
         do j=1,nexternal
            if(i.eq.j.or.i.eq.i_fks.or.i.eq.j_fks)cycle
            if(iHscale.eq.0)then
              if(are_col_conn_H(i,j))then
                if(iRtoB(j).gt.0.and.are_col_conn_S(iRtoB(i),iRtoB(j)))then
c If the H-event colour line corresponds to an S-event colour line that
c connects the two partons whose labels are obtained from mapping the 
c labels at the real-emission level, the content of xscales() must be right
                  SCALUP_tmp_H(i,j)=xscales(i,j)
                else
c Otherwise, an S-event colour line must have been broken, which can happen
c only in the splitting mother->i_fks+j_fks. If j=j_fks, the previous if clause
c must be fulfilled, which leaves one with only the case j=i_fks to deal with
                  if(j.ne.i_fks)then
                    write(*,*)'Error #6 in complete_xmcsubt:',
     #                i,j,i_fks,j_fks
                    stop
                  else
                    SCALUP_tmp_H(i,j)=xscales(i,j_fks)
                  endif
                endif
              endif
            elseif(iHscale.eq.1)then
              if(are_col_conn_H(i,j))then
                SCALUP_tmp_H(i,j)=dipole_mass(p,i,j)
                icount=icount+1
              endif
            else
              write(*,*)'Error in complete_xmcsubt:'
              write(*,*)' unknown iHscale:',iHscale
              stop
            endif
         enddo
         if(iHscale.eq.1)then
           if( (pdg_type(i).eq.21.and.icount.ne.2) .or.
     #         (abs(pdg_type(i)).le.6.and.icount.ne.1) )then
              write(*,*)'Error #1 in complete_xmcsubt:',
     #                  i,pdg_type(i),icount
              stop
           endif
         endif
      enddo
c
c Assignment of shower scale for i_fks and j_fks. The scales for i_fks are not
c present in xscales, while those for j_fks get overwritten below, so to
c have it identical to the scales for i_fks. First clean up
      do i=1,nexternal
        SCALUP_tmp_H(i_fks,i)=-1.d0
        SCALUP_tmp_H(j_fks,i)=-1.d0
      enddo
      jbar=iRtoB(j_fks)
      if(iBtoR(jbar).ne.j_fks)then
        write(*,*)'Error #2 in complete_xmcsubt:',
     #            j_fks,jbar,iRtoB(j_fks),iBtoR(jbar)
      endif
      icount=0
      fksscales(1)=-1.d0
      fksscales(2)=-1.d0
      ifksscl(1)=0
      ifksscl(2)=0
      do i=1,nexternal-1
        if(are_col_conn_S(jbar,i).and.i.ne.jbar)then
          icount=icount+1
          fksscales(icount)=xscales(iBtoR(jbar),iBtoR(i))
          if(abs(IDUP_S(jbar)).le.6)then
c Here and just below: in the case of more sophisticated scale
c assignments (in assign_ifks_Hscale), this could be replaced by:
c  ifksscl(icount)=iBtoR(i)
            ifksscl(icount)=1
          elseif(IDUP_S(jbar).eq.21)then
            if(are_col_conn_H(i_fks,iBtoR(i)))
     #        ifksscl(icount)=1
          endif
        endif
      enddo
      if( (IDUP_S(jbar).eq.21.and.icount.ne.2) .or.
     #    (abs(IDUP_S(jbar)).le.6.and.icount.ne.1) )then
         write(*,*)'Error #3 in complete_xmcsubt:',
     #             jbar,IDUP_S(jbar),icount
         stop
      endif
c Associate a single scale with i_fks and j_fks;
c what follows can be generalised if need be
      call assign_ifks_Hscale(IDUP_S(jbar),ifksscl,fksscales)
      do i=1,nexternal
        if(are_col_conn_H(i_fks,i))SCALUP_tmp_H(i_fks,i)=fksscales(3)
c The assignments for (j_fks,*) stem from MC considerations; in particular,
c they guarantee a symmetric treatment of i_fks and j_fks in the case of
c final-state branchings, and a shower from an initial-state j_fks
c that respects the ordering in t. Note that assigments driven by
c matrix-element considerations could be significantly different
        if(are_col_conn_H(j_fks,i))SCALUP_tmp_H(j_fks,i)=fksscales(3)
      enddo
c Final consistency checks
      do i=1,nexternal
         icount=0
         do j=1,nexternal
           if(are_col_conn_H(i,j))then
             if(SCALUP_tmp_H(i,j).ne.-1.d0)icount=icount+1
           else
             if(SCALUP_tmp_H(i,j).ne.-1.d0)then
               write(*,*)'Error #4 in complete_xmcsubt:',
     #          i,j,pdg_type(i),pdg_type(j),i_fks,j_fks,SCALUP_tmp_H(i,j)
               do k=1,nexternal
                 write(*,*)k,ICOLUP_H(1,k),ICOLUP_H(2,k)
               enddo
               stop
             endif
           endif
         enddo
         if( (pdg_type(i).eq.21.and.icount.ne.2) .or.
     #       (abs(pdg_type(i)).le.6.and.icount.ne.1) )then
           write(*,*)'Error #5 in complete_xmcsubt:',
     #                i,pdg_type(i),icount
           stop
         endif
      enddo
c
c Computation of Delta = wgt_sudakov as the product of Sudakovs between
c starting scales (SCALUP_tmp_S2) and target scales (SCALUP_tmp_H2).
c For initial-state legs, Delta contains a PDF ratio with S-event Bjorken
c fraction and SCALUP_tmp_S2, SCALUP_tmp_H2 scales, see also formula (5.62)
c in Ellis-Stirling-Webber
      wgt_sudakov=1d0
      i_dipole_counter=0
      i_dipole_dead_counter=0
      nG_S=0
      nQ_S=0
      scltarget=1.d10
      sudpdffact=1.d0
      itmp=-1
      jtmp=-1
c
      do i=1,nexternal-1
         if(idup_s(i).eq.21)nG_S=nG_S+1
         if(abs(idup_s(i)).le.6)nQ_S=nQ_S+1
         do j=1,nexternal-1
            if(j.eq.i)cycle
            if(xscales2(i,j).eq.-1d0)cycle
cSF The following definition of startingScale is unprotected:
cSF cstupp must be sufficiently large
            startingScale = min(SCALUP_tmp_S2(i,j),cstupp)
            stoppingScale = SCALUP_tmp_H2(i,j)
c Passing the following if clause must be exceedingly rare
            if(startingScale.le.smallptupp)then
              write(*,*)'Warning in xmcsubt: startingScale, smallptupp'
              write(*,*)startingScale,smallptupp
c$$$              stop
              startingScale=smallptupp
            endif
c only colour-connected partons livezone and sensible scales
c contribute to Delta
            if(xscales2(i,j).eq.-1d0.or.dzones2(i,j))then
               i_dipole_dead_counter=i_dipole_dead_counter+1
               cycle
            endif
            if(stoppingScale.lt.scltarget)then
              scltarget=stoppingScale
              sclstart=startingScale
              itmp=i
              jtmp=j
            endif
            if(stoppingScale.gt.startingScale)then
               i_dipole_dead_counter=i_dipole_dead_counter+1
               cycle
            endif
            if(i.le.2.and.j.le.2)then
               isudtype=1
            elseif(i.gt.2.and.j.gt.2)then
               isudtype=2
            elseif(i.le.2.and.j.gt.2)then
               isudtype=3
            elseif(i.gt.2.and.j.le.2)then
               isudtype=4
            endif
c
            if(stoppingScale.le.smallptlow)then
              deltanum=0.d0
            elseif( stoppingScale.gt.smallptlow .and.
     #              stoppingScale.le.smallptupp )then
              deltanum=pysudakov(smallptupp,xmasses2(i,j),
     &                           idup_s(i),isudtype,mcmass)
              deltanum=deltanum*
     #                 get_to_zero(stoppingScale,smallptlow,smallptupp)
            else
              deltanum=pysudakov(stoppingScale,xmasses2(i,j),
     &                           idup_s(i),isudtype,mcmass)
            endif
            deltaden=pysudakov(startingScale,xmasses2(i,j),
     &                         idup_s(i),isudtype,mcmass)
            if(i.le.nincoming)then
               LP=SIGN(1,LPP(i))
               if (idup_s(i).le.6) then ! (anti-)quark 
                  id=LP*idup_s(i)
               elseif (idup_s(i).eq.21) then ! gluon
                  id=0
               elseif (idup_s(i).eq.22) then ! photon
                  id=7
               endif
               pdffnum=pdg2pdf(abs(lpp(i)),id,
     &                         xbjrk_cnt(i,0),stoppingScale)
               deltanum=deltanum*pdffnum
               pdffden=pdg2pdf(abs(lpp(i)),id,
     &                         xbjrk_cnt(i,0),startingScale)
               deltaden=deltaden*pdffden
               sudpdffact=sudpdffact*pdffnum/pdffden
            endif
            if(deltaden.eq.0.d0)then
              deltarat=1.d0
            else
              deltarat=deltanum/deltaden
            endif
            if(deltarat.ge.1.d0)deltarat=1.d0
            if(deltarat.le.0.d0)deltarat=0.d0
            wgt_sudakov=wgt_sudakov*deltarat
            i_dipole_counter=i_dipole_counter+1
         enddo
      enddo
      if(i_dipole_counter+i_dipole_dead_counter.ne.nQ_S+2*nG_S)then
         write(*,*)'Mismatch in number of dipole ends and Delta factors'
         write(*,*)i_dipole_counter+i_dipole_dead_counter,nQ_S,nG_S
         stop
      endif
      if(itmp.eq.-1.or.jtmp.eq.-1)then
        scltarget=-1.d8
        sclstart=-1.d8
      endif

cSF THE FOLLOWING MIGHT READ
c        if(UseDelta)probne = wgt_sudakov
cSF IF A UseDelta parameter is introduced
      probne = wgt_sudakov
cSF ELIMINATE WHAT FOLLOWS AFTER CHECKS
      if(probne.lt.0.d0)then
        write(*,*)'SFWARNING1',probne
        probne=0.d0
      endif
      if(probne.gt.1.d0)then
        write(*,*)'SFWARNING2',probne
        probne=1.d0
      endif
c
      do i=1,nexternal
         if(i.le.ipartners(0))xmcxsec(i)=xmcxsec(i)*probne
         if(i.gt.ipartners(0))xmcxsec(i)=0d0
      enddo

      return
      end


      function get_to_zero(sc,xlow,xupp)
      implicit none
      double precision get_to_zero,sc,xlow,xupp
      double precision x,tmp,emscafun
c
      if(sc.le.xlow)then
        tmp=0d0
      elseif(sc.le.xupp)then
        x=(xupp-sc)/(xupp-xlow)
        tmp=1-emscafun(x,2d0)
      else
        tmp=1d0
      endif
      get_to_zero=tmp
      return
      end


      function dipole_mass(p,i,j)
      implicit none
      include 'nexternal.inc'
      double precision dipole_mass,sign,tmp
      double precision p(0:3,nexternal)
      integer i,j,k
c
      sign=1.d0
      if(i.le.2)sign=-sign
      if(j.le.2)sign=-sign
      tmp=(p(0,i)+sign*p(0,j))**2
      do k=1,3
        tmp=tmp-(p(k,i)+sign*p(k,j))**2
      enddo
      tmp=sqrt(max(0.d0,tmp))
      dipole_mass=tmp
      return
      end


      subroutine assign_ifks_Hscale(ipdg,ifksscl,fksscales)
      implicit none
      double precision fksscales(3)
      integer ipdg,ifksscl(2),i,icount
      integer i_fks,j_fks
      common/fks_indices/i_fks,j_fks
      logical wrong
c Set itype=0 to set scale according to the colour line to which i_fks belongs
c     itype=1 to take the minimum of the two scales in the case of mother=gluon
      integer itype
      parameter (itype=0)
c
      fksscales(3)=1.d10
      if(itype.eq.0)then
        wrong=.not.((ifksscl(1).eq.1.and.ifksscl(2).eq.0).or.
     &              (ifksscl(1).eq.0.and.ifksscl(2).eq.1))
        if(wrong)then
          write(*,*)'Something wrong in assign_ifks_Hscale (0):'
          write(*,*)ipdg,icount,i_fks,j_fks
          write(*,*)ifksscl(1),ifksscl(2),fksscales(1),fksscales(2)
          stop
        endif
        if(ifksscl(1).ne.0)fksscales(3)=fksscales(1)
        if(ifksscl(2).ne.0)fksscales(3)=fksscales(2)
      elseif(itype.eq.1)then
        do i=1,2
          if(fksscales(i).gt.0d0)
     #      fksscales(3)=min(fksscales(3),fksscales(i))
        enddo
      endif
      if(fksscales(3).eq.1.d10)then
        write(*,*)'Could not assign scale in assign_ifks_Hscale:'
        write(*,*)ipdg,ifksscl(1),ifksscl(2),fksscales(1),fksscales(2)
        stop
      endif
      return
      end


      subroutine get_mbar(p,y_ij_fks,ileg,bornbars,bornbarstilde)
c Computes barred amplitudes (bornbars) squared according
c to Odagiri's prescription (hep-ph/9806531).
c Computes barred azimuthal amplitudes (bornbarstilde) with
c the same method 
      implicit none

      include "genps.inc"
      include 'nexternal.inc'
      include "born_nhel.inc"

      double precision p(0:3,nexternal)
      double precision y_ij_fks,bornbars(max_bcol),bornbarstilde(max_bcol)

      double precision zero
      parameter (zero=0.d0)
      double precision p_born_rot(0:3,nexternal-1)

      integer imother_fks,ileg

      double precision p_born(0:3,nexternal-1)
      common/pborn/p_born

      double Precision amp2(ngraphs), jamp2(0:ncolor)
      common/to_amps/  amp2,       jamp2

      integer i_fks,j_fks
      common/fks_indices/i_fks,j_fks

      double complex wgt1(2),W1(6),W2(6),W3(6),W4(6),Wij_angle,Wij_recta
      double complex azifact

      double complex xij_aor
      common/cxij_aor/xij_aor

      double precision born,sumborn,borntilde
      integer i

      double precision vtiny,pi(0:3),pj(0:3),cphi_mother,sphi_mother
      parameter (vtiny=1d-8)
      double complex ximag
      parameter (ximag=(0.d0,1.d0))

      double precision xi_i_fks_ev,y_ij_fks_ev,t
      double precision p_i_fks_ev(0:3),p_i_fks_cnt(0:3,-2:2)
      common/fksvariables/xi_i_fks_ev,y_ij_fks_ev,p_i_fks_ev,p_i_fks_cnt

      double precision cthbe,sthbe,cphibe,sphibe
      common/cbeangles/cthbe,sthbe,cphibe,sphibe

      logical calculatedBorn
      common/ccalculatedBorn/calculatedBorn
      double precision iden_comp
      common /c_iden_comp/iden_comp

c Particle types (=color) of i_fks, j_fks and fks_mother
      integer i_type,j_type,m_type
      common/cparticle_types/i_type,j_type,m_type
c
      logical is_leading_cflow(max_bcol)
      integer num_leading_cflows
      common/c_leading_cflows/is_leading_cflow,num_leading_cflows
      
c
c BORN
      call sborn(p_born,wgt1)
      born=dble(wgt1(1))

c born is the total born amplitude squared
      sumborn=0.d0
      do i=1,max_bcol
         if(is_leading_cflow(i))sumborn=sumborn+jamp2(i)
c sumborn is the sum of the leading-color amplitudes squared
      enddo
      

c BORN TILDE
      if(ileg.eq.1.or.ileg.eq.2)then
         if(j_fks.eq.2 .and. nexternal-1.ne.3)then
c Rotation according to innerpin.m. Use rotate_invar() if a more 
c general rotation is needed.
c Exclude 2->1 (at the Born level) processes: matrix elements are
c independent of the PS point, but non-zero helicity configurations
c might flip when rotating the momenta.
            do i=1,nexternal-1
               p_born_rot(0,i)=p_born(0,i)
               p_born_rot(1,i)=-p_born(1,i)
               p_born_rot(2,i)=p_born(2,i)
               p_born_rot(3,i)=-p_born(3,i)
            enddo
            calculatedBorn=.false.
            call sborn(p_born_rot,wgt1)
            calculatedBorn=.false.
         else
            call sborn(p_born,wgt1)
         endif
         if (abs(m_type).eq.3) then
            wgt1(2)=0d0
         else
c Insert <ij>/[ij] which is not included by sborn()
            if (1d0-y_ij_fks.lt.vtiny)then
               azifact=xij_aor
            else
               do i=0,3
                  pi(i)=p_i_fks_ev(i)
                  pj(i)=p(i,j_fks)
               enddo
               if(j_fks.eq.2)then
c Rotation according to innerpin.m. Use rotate_invar() if a more 
c general rotation is needed
                  pi(1)=-pi(1)
                  pi(3)=-pi(3)
                  pj(1)=-pj(1)
                  pj(3)=-pj(3)
               endif
               CALL IXXXSO(pi ,ZERO ,+1,+1,W1)        
               CALL OXXXSO(pj ,ZERO ,-1,+1,W2)        
               CALL IXXXSO(pi ,ZERO ,-1,+1,W3)        
               CALL OXXXSO(pj ,ZERO ,+1,+1,W4)        
               Wij_angle=(0d0,0d0)
               Wij_recta=(0d0,0d0)
               do i=1,4
                  Wij_angle = Wij_angle + W1(i)*W2(i)
                  Wij_recta = Wij_recta + W3(i)*W4(i)
               enddo
               azifact=Wij_angle/Wij_recta
            endif
c Insert the extra factor due to Madgraph convention for polarization vectors
            if(j_fks.eq.2)then
               cphi_mother=-1.d0
               sphi_mother=0.d0
            else
               cphi_mother=1.d0
               sphi_mother=0.d0
            endif
            wgt1(2) = -(cphi_mother+ximag*sphi_mother)**2 *
     #                wgt1(2) * dconjg(azifact)
         endif
      elseif(ileg.eq.3.or.ileg.eq.4)then
         if(abs(j_type).eq.3.and.i_type.eq.8)then
            wgt1(2)=0.d0
         elseif(m_type.eq.8)then
c Insert <ij>/[ij] which is not included by sborn()
            if(1.d0-y_ij_fks.lt.vtiny)then
               azifact=xij_aor
            else
               do i=0,3
                  pi(i)=p_i_fks_ev(i)
                  pj(i)=p(i,j_fks)
               enddo
               CALL IXXXSO(pi ,ZERO ,+1,+1,W1)        
               CALL OXXXSO(pj ,ZERO ,-1,+1,W2)        
               CALL IXXXSO(pi ,ZERO ,-1,+1,W3)        
               CALL OXXXSO(pj ,ZERO ,+1,+1,W4)        
               Wij_angle=(0d0,0d0)
               Wij_recta=(0d0,0d0)
               do i=1,4
                  Wij_angle = Wij_angle + W1(i)*W2(i)
                  Wij_recta = Wij_recta + W3(i)*W4(i)
               enddo
               azifact=Wij_angle/Wij_recta
            endif
c Insert the extra factor due to Madgraph convention for polarization vectors
            imother_fks=min(i_fks,j_fks)
            call getaziangles(p_born(0,imother_fks),
     #                        cphi_mother,sphi_mother)
            wgt1(2) = -(cphi_mother-ximag*sphi_mother)**2 *
     #                  wgt1(2) * azifact
         else
            write(*,*)'FATAL ERROR in get_mbar',
     #           i_type,j_type,i_fks,j_fks
            stop
         endif
      else
         write(*,*)'unknown ileg in get_mbar'
         stop
      endif

      borntilde=dble(wgt1(2))


c BARRED AMPLITUDES
      do i=1,max_bcol
         if (sumborn.ne.0d0.and.is_leading_cflow(i)) then
            bornbars(i)=jamp2(i)/sumborn * born *iden_comp
         elseif (born.eq.0d0 .or. jamp2(i).eq.0d0
     &           .or..not.is_leading_cflow(i)) then
            bornbars(i)=0d0
         else
            write (*,*) 'ERROR #1, dividing by zero'
            stop
         endif
         if (sumborn.ne.0d0.and.is_leading_cflow(i)) then
            bornbarstilde(i)=jamp2(i)/sumborn * borntilde *iden_comp
         elseif (borntilde.eq.0d0 .or. jamp2(i).eq.0d0
     &           .or..not.is_leading_cflow(i)) then
            bornbarstilde(i)=0d0
         else
            write (*,*) 'ERROR #2, dividing by zero'
            stop
         endif      
c bornbars(i) is the i-th leading-color amplitude squared re-weighted
c in such a way that the sum of bornbars(i) is born rather than sumborn.
c the same holds for bornbarstilde(i).
      enddo

      return
      end



      function gfunction(w,alpha,beta,delta)
c Gets smoothly to 0 as w goes to 1.
c Call with
c   alpha > 1, or alpha < 0; if alpha < 0, gfunction = 1;
c   0 < |beta| <= 1;
c   0 < delta <= 2.
      implicit none
      double precision tiny
      parameter (tiny=1.d-5)
      double precision gfunction,alpha,beta,delta,w,wmin,wg,tt,tmp
      logical firsttime
      save firsttime
      data firsttime /.true./
      double precision cutoff,cutoff2
      parameter(cutoff=1d0)
      parameter(cutoff2=0.99d0)
c
c     set cutoff < 1 and cutoff2 = cutoff in the final version
c
      if(firsttime)then
        firsttime=.false.
        if(alpha.ge.0d0.and.alpha.lt.1d0)then
          write(*,*)'Incorrect alpha in gfunction',alpha
          stop
        endif
        if(abs(beta).gt.1d0)then
          write(*,*)'Incorrect beta in gfunction',beta
          stop
        endif
        if(delta.gt.2d0.or.delta.le.0d0)then
          write(*,*)'Incorrect delta in gfunction',delta
          stop
        endif
      endif
c
      tmp=1d0
      if(alpha.gt.0d0)then
        if(beta.lt.0d0)then
          wmin=0d0
        else
          wmin=max(0d0,1d0-delta)
        endif
        wg=min(1d0-(1-wmin)*abs(beta),cutoff-tiny)
        if(abs(w).gt.wg.and.abs(w).lt.cutoff2)then
          tt=(abs(w)-wg)/(cutoff-wg)
          if(tt.gt.1d0)then
            write(*,*)'Fatal error in gfunction',tt
            stop
          endif
          tmp=(1-tt)**(2*alpha)/(tt**(2*alpha)+(1-tt)**(2*alpha))
        elseif(abs(w).ge.cutoff2)then
          tmp=0d0
        endif
      endif
      gfunction=tmp
      return
      end



      subroutine kinematics_driver(xi_i_fks,y_ij_fks,sh,pp,ileg,
     &                             xm12,xm22,xtk,xuk,xq1q,xq2q,qMC,extra)
c Determines Mandelstam invariants and assigns ileg and shower-damping
c variable qMC
      implicit none
      include "nexternal.inc"
      include "coupl.inc"
      include "run.inc"
      double precision pp(0:3,nexternal),pp_rec(0:3)
      double precision xi_i_fks,y_ij_fks,xij

      integer ileg,j,i,nfinal
      double precision xp1(0:3),xp2(0:3),xk1(0:3),xk2(0:3),xk3(0:3)
      common/cpkmomenta/xp1,xp2,xk1,xk2,xk3
      double precision sh,xtk,xuk,w1,w2,xq1q,xq2q,xm12,xm22
      double precision qMC,zPY8,zeta1,zeta2,get_zeta,z,qMCarg,dot
      logical extra
      double precision p_born(0:3,nexternal-1)
      common/pborn/p_born
      integer fksfather
      common/cfksfather/fksfather
      integer i_fks,j_fks
      common/fks_indices/i_fks,j_fks
      double precision tiny
      parameter(tiny=1d-5)
      double precision zero
      parameter(zero=0d0)

      integer isqrtneg
      save isqrtneg

      double precision pmass(nexternal)
      include "pmass.inc"

c Initialise
      do i=0,3
         pp_rec(i)=0d0
         xp1(i)=0d0
         xp2(i)=0d0
         xk1(i)=0d0
         xk2(i)=0d0
         xk3(i)=0d0
      enddo
      nfinal=nexternal-2
      xm12=0d0
      xm22=0d0
      xq1q=0d0
      xq2q=0d0
      qMC=-1d0

c Discard if unphysical FKS variables
      if(xi_i_fks.lt.0d0.or.xi_i_fks.gt.1d0.or.
     &   abs(y_ij_fks).gt.1d0)then
         write(*,*)'Error 0 in kinematics_driver: fks variables'
         write(*,*)xi_i_fks,y_ij_fks
         stop
      endif

c Determine ileg
c ileg = 1 ==> emission from left     incoming parton
c ileg = 2 ==> emission from right    incoming parton
c ileg = 3 ==> emission from massive  outgoing parton
c ileg = 4 ==> emission from massless outgoing parton
c Instead of pmass(j_fks), one should use pmass(fksfather), but the
c kernels where pmass(fksfather) != pmass(j_fks) are non-singular
      if(fksfather.le.2)then
        ileg=fksfather
      elseif(pmass(j_fks).ne.0d0)then
        ileg=3
      elseif(pmass(j_fks).eq.0d0)then
        ileg=4
      else
        write(*,*)'Error 1 in kinematics_driver: unknown ileg'
        write(*,*)ileg,fksfather,pmass(j_fks)
        stop
      endif

c Determine and assign momenta:
c xp1 = incoming left parton  (emitter (recoiler) if ileg = 1 (2))
c xp2 = incoming right parton (emitter (recoiler) if ileg = 2 (1))
c xk1 = outgoing parton       (emitter (recoiler) if ileg = 3 (4))
c xk2 = outgoing parton       (emitter (recoiler) if ileg = 4 (3))
c xk3 = extra parton          (FKS parton)
      do j=0,3
c xk1 and xk2 are never used for ISR
         xp1(j)=pp(j,1)
         xp2(j)=pp(j,2)
         xk3(j)=pp(j,i_fks)
         if(ileg.gt.2)pp_rec(j)=pp(j,1)+pp(j,2)-pp(j,i_fks)-pp(j,j_fks)
         if(ileg.eq.3)then
            xk1(j)=pp(j,j_fks)
            xk2(j)=pp_rec(j)
         elseif(ileg.eq.4)then
            xk1(j)=pp_rec(j)
            xk2(j)=pp(j,j_fks)
         endif
      enddo

c Determine the Mandelstam invariants needed in the MC functions in terms
c of FKS variables: the argument of MC functions are (p+k)^2, NOT 2 p.k
c
c Definitions of invariants in terms of momenta
c
c xm12 =     xk1 . xk1
c xm22 =     xk2 . xk2
c xtk  = - 2 xp1 . xk3
c xuk  = - 2 xp2 . xk3
c xq1q = - 2 xp1 . xk1 + xm12
c xq2q = - 2 xp2 . xk2 + xm22
c w1   = + 2 xk1 . xk3        = - xq1q + xq2q - xtk
c w2   = + 2 xk2 . xk3        = - xq2q + xq1q - xuk
c xq1c = - 2 xp1 . xk2        = - s - xtk - xq1q + xm12
c xq2c = - 2 xp2 . xk1        = - s - xuk - xq2q + xm22
c
c Parametrisation of invariants in terms of FKS variables
c
c ileg = 1
c xp1  =  sqrt(s)/2 * ( 1 , 0 , 0 , 1 )
c xp2  =  sqrt(s)/2 * ( 1 , 0 , 0 , -1 )
c xk3  =  B * ( 1 , 0 , sqrt(1-yi**2) , yi )
c xk1  =  irrelevant
c xk2  =  irrelevant
c yi = y_ij_fks
c x = 1 - xi_i_fks
c B = sqrt(s)/2*(1-x)
c
c ileg = 2
c xp1  =  sqrt(s)/2 * ( 1 , 0 , 0 , 1 )
c xp2  =  sqrt(s)/2 * ( 1 , 0 , 0 , -1 )
c xk3  =  B * ( 1 , 0 , sqrt(1-yi**2) , -yi )
c xk1  =  irrelevant
c xk2  =  irrelevant
c yi = y_ij_fks
c x = 1 - xi_i_fks
c B = sqrt(s)/2*(1-x)
c
c ileg = 3
c xp1  =  sqrt(s)/2 * ( 1 , 0 , sqrt(1-yi**2) , yi )
c xp2  =  sqrt(s)/2 * ( 1 , 0 , -sqrt(1-yi**2) , -yi )
c xk1  =  ( sqrt(veckn_ev**2+xm12) , 0 , 0 , veckn_ev )
c xk2  =  xp1 + xp2 - xk1 - xk3
c xk3  =  B * ( 1 , 0 , sqrt(1-yj**2) , yj )
c yj = y_ij_fks
c yi = irrelevant
c x = 1 - xi_i_fks
c veckn_ev is such that xk2**2 = xm22
c B = sqrt(s)/2*(1-x)
c azimuth = irrelevant (hence set = 0)
c
c ileg = 4
c xp1  =  sqrt(s)/2 * ( 1 , 0 , sqrt(1-yi**2) , yi )
c xp2  =  sqrt(s)/2 * ( 1 , 0 , -sqrt(1-yi**2) , -yi )
c xk1  =  xp1 + xp2 - xk2 - xk3
c xk2  =  A * ( 1 , 0 , 0 , 1 )
c xk3  =  B * ( 1 , 0 , sqrt(1-yj**2) , yj )
c yj = y_ij_fks
c yi = irrelevant
c x = 1 - xi_i_fks
c A = (s*x-xm12)/(sqrt(s)*(2-(1-x)*(1-yj)))
c B = sqrt(s)/2*(1-x)
c azimuth = irrelevant (hence set = 0)

      if(ileg.eq.1)then
         xtk=-sh*xi_i_fks*(1-y_ij_fks)/2
         xuk=-sh*xi_i_fks*(1+y_ij_fks)/2
         if(extra)then
            if(shower_mc.eq.'HERWIG6'  .or.
     &         shower_mc.eq.'HERWIGPP' )qMC=xi_i_fks/2*sqrt(sh*(1-y_ij_fks**2))
            if(shower_mc.eq.'PYTHIA6Q' )qMC=sqrt(-xtk)
            if(shower_mc.eq.'PYTHIA6PT'.or.
     &         shower_mc.eq.'PYTHIA8'  )qMC=sqrt(-xtk*xi_i_fks)
         endif
      elseif(ileg.eq.2)then
         xtk=-sh*xi_i_fks*(1+y_ij_fks)/2
         xuk=-sh*xi_i_fks*(1-y_ij_fks)/2
         if(extra)then
            if(shower_mc.eq.'HERWIG6'  .or.
     &         shower_mc.eq.'HERWIGPP' )qMC=xi_i_fks/2*sqrt(sh*(1-y_ij_fks**2))
            if(shower_mc.eq.'PYTHIA6Q' )qMC=sqrt(-xuk)
            if(shower_mc.eq.'PYTHIA6PT'.or.
     &         shower_mc.eq.'PYTHIA8'  )qMC=sqrt(-xuk*xi_i_fks)
         endif
      elseif(ileg.eq.3)then
         xm12=pmass(j_fks)**2
         xm22=dot(pp_rec,pp_rec)
         xtk=-2*dot(xp1,xk3)
         xuk=-2*dot(xp2,xk3)
         xq1q=-2*dot(xp1,xk1)+xm12
         xq2q=-2*dot(xp2,xk2)+xm22
         w1=-xq1q+xq2q-xtk
         w2=-xq2q+xq1q-xuk
         if(extra)then
            if(shower_mc.eq.'HERWIG6'.or.
     &         shower_mc.eq.'HERWIGPP')then
               zeta1=get_zeta(sh,w1,w2,xm12,xm22)
               qMCarg=zeta1*((1-zeta1)*w1-zeta1*xm12)
               if(qMCarg.lt.0d0.and.qMCarg.ge.-tiny)qMCarg=0d0
               if(qMCarg.lt.-tiny) then
                  isqrtneg=isqrtneg+1
                  write(*,*)'Error 2 in kinematics_driver: negtive sqrt'
                  write(*,*)qMCarg,isqrtneg
                  if(isqrtneg.ge.100)stop
               endif
               qMC=sqrt(qMCarg)
            elseif(shower_mc.eq.'PYTHIA6Q')then
               qMC=sqrt(w1+xm12)
            elseif(shower_mc.eq.'PYTHIA6PT')then
               write(*,*)'PYTHIA6PT not available for FSR'
               stop
            elseif(shower_mc.eq.'PYTHIA8')then
               z=zPY8(ileg,xm12,xm22,sh,1d0-xi_i_fks,0d0,y_ij_fks,xtk,xuk,xq1q,xq2q)
               qMC=sqrt(z*(1-z)*w1)
            endif
         endif
      elseif(ileg.eq.4)then
         xm12=dot(pp_rec,pp_rec)
         xm22=0d0
         xtk=-2*dot(xp1,xk3)
         xuk=-2*dot(xp2,xk3)
         xij=2*(1-xm12/sh-xi_i_fks)/(2-xi_i_fks*(1-y_ij_fks))
         w2=sh*xi_i_fks*xij*(1-y_ij_fks)/2d0
         xq2q=-sh*xij*(2-dot(xp1,xk2)*4d0/(sh*xij))/2d0
         xq1q=xuk+xq2q+w2
         w1=-xq1q+xq2q-xtk
         if(extra)then
            if(shower_mc.eq.'HERWIG6'.or.
     &         shower_mc.eq.'HERWIGPP')then
               zeta2=get_zeta(sh,w2,w1,xm22,xm12)
               qMCarg=zeta2*(1-zeta2)*w2
               if(qMCarg.lt.0d0.and.qMCarg.ge.-tiny)qMCarg=0d0
               if(qMCarg.lt.-tiny)then
                  isqrtneg=isqrtneg+1
                  write(*,*)'Error 3 in kinematics_driver: negtive sqrt'
                  write(*,*)qMCarg,isqrtneg
                  if(isqrtneg.ge.100)stop
               endif
               qMC=sqrt(qMCarg)
            elseif(shower_mc.eq.'PYTHIA6Q')then
               qMC=sqrt(w2)
            elseif(shower_mc.eq.'PYTHIA6PT')then
               write(*,*)'PYTHIA6PT not available for FSR'
               stop
            elseif(shower_mc.eq.'PYTHIA8')then
               z=zPY8(ileg,xm12,xm22,sh,1d0-xi_i_fks,0d0,y_ij_fks,xtk,xuk,xq1q,xq2q)
               qMC=sqrt(z*(1-z)*w2)
            endif
         endif
      else
         write(*,*)'Error 4 in kinematics_driver: assigned wrong ileg'
         stop
      endif

c Checks on invariants
      call check_invariants(ileg,sh,xtk,xuk,w1,w2,xq1q,xq2q,xm12,xm22)

      return
      end



      subroutine check_invariants(ileg,sh,xtk,xuk,w1,w2,xq1q,xq2q,xm12,xm22)
      implicit none
      integer ileg

      double precision xp1(0:3),xp2(0:3),xk1(0:3),xk2(0:3),xk3(0:3)
      common/cpkmomenta/xp1,xp2,xk1,xk2,xk3
      double precision sh,xtk,xuk,w1,w2,xq1q,xq2q,xm12,xm22

      double precision tiny,dot
      parameter(tiny=1d-5)

      if(ileg.le.2)then
         if((abs(xtk+2*dot(xp1,xk3))/sh.ge.tiny).or.
     &      (abs(xuk+2*dot(xp2,xk3))/sh.ge.tiny))then
            write(*,*)'Imprecision 1 in check_invariants'
            write(*,*)abs(xtk+2*dot(xp1,xk3))/sh,
     &                abs(xuk+2*dot(xp2,xk3))/sh
            stop
         endif
c
      elseif(ileg.eq.3)then
         if(sqrt(w1+xm12).ge.sqrt(sh)-sqrt(xm22))then
            write(*,*)'Imprecision 2a in check_invariants'
            write(*,*)sqrt(w1),sqrt(sh),xm22
            stop
         endif
         if(((abs(w1-2*dot(xk1,xk3))/sh.ge.tiny)).or.
     &      ((abs(w2-2*dot(xk2,xk3))/sh.ge.tiny)))then
            write(*,*)'Imprecision 2b in check_invariants'
            write(*,*)abs(w1-2*dot(xk1,xk3))/sh,
     &                abs(w2-2*dot(xk2,xk3))/sh
            stop
         endif
         if(xm12.eq.0d0)then
            write(*,*)'Error 2c in check_invariants'
            stop
         endif
c
      elseif(ileg.eq.4)then
         if(sqrt(w2).ge.sqrt(sh)-sqrt(xm12))then
            write(*,*)'Imprecision 3a in check_invariants'
            write(*,*)sqrt(w2),sqrt(sh),xm12
            stop
         endif
         if(((abs(w2-2*dot(xk2,xk3))/sh.ge.tiny)).or.
     &      ((abs(xq2q+2*dot(xp2,xk2))/sh.ge.tiny)).or.
     &      ((abs(xq1q+2*dot(xp1,xk1)-xm12)/sh.ge.tiny)))then
            write(*,*)'Imprecision 3b in check_invariants'
            write(*,*)abs(w2-2*dot(xk2,xk3))/sh,
     &                abs(xq2q+2*dot(xp2,xk2))/sh,
     &                abs(xq1q+2*dot(xp1,xk1)-xm12)/sh
            stop
         endif
         if(xm22.ne.0d0)then
            write(*,*)'Error 3c in check_invariants'
            stop
         endif

      endif

      return
      end



c Monte Carlo functions
c
c The invariants given in input to these routines follow FNR conventions
c (i.e., are defined as (p+k)^2, NOT 2 p.k). 
c The invariants used inside these routines follow MNR conventions
c (i.e., are defined as -2p.k, NOT (p+k)^2)

c Herwig6

      function zHW6(ileg,e0sq,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
c Shower energy variable
      implicit none
      integer ileg
      double precision zHW6,e0sq,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q,tiny,
     &ss,w1,w2,tbeta1,zeta1,tbeta2,zeta2,get_zeta,beta,betae0,betad,betas
      parameter (tiny=1d-5)

      w1=-xq1q+xq2q-xtk
      w2=-xq2q+xq1q-xuk
c
      if(ileg.eq.1)then
         if(1-x.lt.tiny)then
            zHW6=1-(1-x)*(s*(1-yi)+4*e0sq*(1+yi))/(8*e0sq)
         elseif(1-yi.lt.tiny)then
            zHW6=x-(1-yi)*(1-x)*(s*x**2-4*e0sq)/(8*e0sq)
         else
            ss=1-(1+xuk/s)/(e0sq/xtk)
            if(ss.lt.0d0)goto 999
            zHW6=2*(e0sq/xtk)*(1-sqrt(ss))
         endif
c
      elseif(ileg.eq.2)then
         if(1-x.lt.tiny)then
            zHW6=1-(1-x)*(s*(1-yi)+4*e0sq*(1+yi))/(8*e0sq)
         elseif(1-yi.lt.tiny)then
            zHW6=x-(1-yi)*(1-x)*(s*x**2-4*e0sq)/(8*e0sq)
         else
            ss=1-(1+xtk/s)/(e0sq/xuk)
            if(ss.lt.0d0)goto 999
            zHW6=2*(e0sq/xuk)*(1-sqrt(ss))
         endif
c
      elseif(ileg.eq.3)then
         if(e0sq.le.(w1+xm12))goto 999
         if(1-x.lt.tiny)then
            beta=1-xm12/s
            betae0=sqrt(1-xm12/e0sq)
            betad=sqrt((1-(xm12-xm22)/s)**2-(4*xm22/s))
            betas=1+(xm12-xm22)/s
            zHW6=1+(1-x)*( s*(yj*betad-betas)/(4*e0sq*(1+betae0))-
     &                     betae0*(xm12-xm22+s*(1+(1+yj)*betad-betas))/
     &                     (betad*(xm12-xm22+s*(1+betad))) )
         else
            tbeta1=sqrt(1-(w1+xm12)/e0sq)
            zeta1=get_zeta(s,w1,w2,xm12,xm22)
            zHW6=1-tbeta1*zeta1-w1/(2*(1+tbeta1)*e0sq)
         endif
c
      elseif(ileg.eq.4)then
         if(e0sq.le.w2)goto 999
         if(1-x.lt.tiny)then
            zHW6=1-(1-x)*( (s-xm12)*(1-yj)/(8*e0sq)+
     &                     s*(1+yj)/(2*(s-xm12)) )
         elseif(1-yj.lt.tiny)then
            zHW6=(s*x-xm12)/(s-xm12)+(1-yj)*(1-x)*(s*x-xm12)*
     &           ( (s-xm12)**2*(s*(1-2*x)+xm12)+
     &             4*e0sq*s*(s*x-xm12*(2-x)) )/
     &           ( 8*e0sq*(s-xm12)**3 )
         else
            tbeta2=sqrt(1-w2/e0sq)
            zeta2=get_zeta(s,w2,w1,xm22,xm12)
            zHW6=1-tbeta2*zeta2-w2/(2*(1+tbeta2)*e0sq)
         endif
c
      else
         write(*,*)'zHW6: unknown ileg'
         stop
      endif

      if(zHW6.lt.0d0.or.zHW6.gt.1d0)goto 999

      return
 999  continue
      zHW6=-1d0
      return
      end



      function xiHW6(ileg,e0sq,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
c Shower evolution variable
      implicit none
      integer ileg
      double precision xiHW6,e0sq,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q,tiny,z,
     &zHW6,w1,w2,beta,betae0,betad,betas
      parameter (tiny=1d-5)

      z=zHW6(ileg,e0sq,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
      if(z.lt.0d0)goto 999
      w1=-xq1q+xq2q-xtk
      w2=-xq2q+xq1q-xuk
c
      if(ileg.eq.1)then
         if(1-x.lt.tiny)then
            xiHW6=2*s*(1-yi)/(s*(1-yi)+4*e0sq*(1+yi))
         elseif(1-yi.lt.tiny)then
            xiHW6=(1-yi)*s*x**2/(4*e0sq)
         else
            xiHW6=2*(1+xuk/(s*(1-z)))
         endif
c
      elseif(ileg.eq.2)then
         if(1-x.lt.tiny)then
            xiHW6=2*s*(1-yi)/(s*(1-yi)+4*e0sq*(1+yi))
         elseif(1-yi.lt.tiny)then
            xiHW6=(1-yi)*s*x**2/(4*e0sq)
         else
            xiHW6=2*(1+xtk/(s*(1-z)))
         endif
c
      elseif(ileg.eq.3)then
         if(e0sq.le.(w1+xm12))goto 999
         if(1-x.lt.tiny)then
            beta=1-xm12/s
            betae0=sqrt(1-xm12/e0sq)
            betad=sqrt((1-(xm12-xm22)/s)**2-(4*xm22/s))
            betas=1+(xm12-xm22)/s
            xiHW6=( s*(1+betae0)*betad*(xm12-xm22+s*(1+betad))*(yj*betad-betas) )/
     &            ( -4*e0sq*betae0*(1+betae0)*(xm12-xm22+s*(1+(1+yj)*betad-betas))+
     &              (s*betad*(xm12-xm22+s*(1+betad))*(yj*betad-betas)) )
         else
            xiHW6=w1/(2*z*(1-z)*e0sq)
         endif
c
      elseif(ileg.eq.4)then
         if(e0sq.le.w2)goto 999
         if(1-x.lt.tiny)then
            xiHW6=2*(s-xm12)**2*(1-yj)/( (s-xm12)**2*(1-yj)+4*e0sq*s*(1+yj) )
         elseif(1-yj.lt.tiny)then
            xiHW6=(s-xm12)**2*(1-yj)/(4*e0sq*s)
         else
            xiHW6=w2/(2*z*(1-z)*e0sq)
         endif
c
      else
         write(*,*)'xiHW6: unknown ileg'
         stop
      endif

      if(xiHW6.lt.0d0)goto 999

      return
 999  continue
      xiHW6=-1d0
      return
      end



      function xjacHW6_xiztoxy(ileg,e0sq,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
c Returns the jacobian d(z,xi)/d(x,y), where z and xi are the shower 
c variables, and x and y are FKS variables
      implicit none
      integer ileg
      double precision xjacHW6_xiztoxy,e0sq,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,
     &xq2q,tiny,tmp,z,zHW6,xi,xiHW6,w1,w2,tbeta1,zeta1,dw1dx,dw2dx,dw1dy,dw2dy,
     &tbeta2,get_zeta,beta,betae0,betad,betas,eps1,eps2,beta1,beta2,zmo
      parameter (tiny=1d-5)

      tmp=0d0
      z=zHW6(ileg,e0sq,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
      xi=xiHW6(ileg,e0sq,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
      if(z.lt.0d0.or.xi.lt.0d0)goto 999
      w1=-xq1q+xq2q-xtk
      w2=-xq2q+xq1q-xuk
c
      if(ileg.eq.1)then
         if(1-x.lt.tiny)then
            tmp=-2*s/(s*(1-yi)+4*(1+yi)*e0sq)
         elseif(1-yi.lt.tiny)then
            tmp=-s*x**2/(4*e0sq)
         else
            tmp=-s*(1-x)*z**3/(4*e0sq*(1-z)*(xi*(1-z)+z))
         endif
c
      elseif(ileg.eq.2)then
         if(1-x.lt.tiny)then
            tmp=-2*s/(s*(1-yi)+4*(1+yi)*e0sq)
         elseif(1-yi.lt.tiny)then
            tmp=-s*x**2/(4*e0sq)
         else
            tmp=-s*(1-x)*z**3/(4*e0sq*(1-z)*(xi*(1-z)+z))
         endif
c
      elseif(ileg.eq.3)then
         if(e0sq.le.(w1+xm12))goto 999
         if(1-x.lt.tiny)then
            beta=1-xm12/s
            betae0=sqrt(1-xm12/e0sq)
            betad=sqrt((1-(xm12-xm22)/s)**2-(4*xm22/s))
            betas=1+(xm12-xm22)/s
            tmp=( s*betae0*(1+betae0)*betad*(xm12-xm22+s*(1+betad)) )/
     &          ( (-4*e0sq*(1+betae0)*(xm12-xm22+s*(1+betad*(1+yj)-betas)))+
     &            (xm12-xm22+s*(1+betad))*(xm12*(4+yj*betad-betas)-
     &            (xm22-s)*(yj*betad-betas)) )
         else
            eps2=1-(xm12-xm22)/(s-w1)
            beta2=sqrt(eps2**2-4*s*xm22/(s-w1)**2)
            tbeta1=sqrt(1-(w1+xm12)/e0sq)
            call dinvariants_dFKS(ileg,s,x,yi,yj,xm12,xm22,dw1dx,dw1dy,dw2dx,dw2dy)
            tmp=-(dw1dy*dw2dx-dw1dx*dw2dy)*tbeta1/(2*e0sq*z*(1-z)*(s-w1)*beta2)
         endif
c
      elseif(ileg.eq.4)then
         if(e0sq.le.w2)goto 999
         if(1-x.lt.tiny)then
            zmo=(s-xm12)*(1-yj)/(8*e0sq)+s*(1+yj)/(2*(s-xm12))
            tmp=-s/(4*e0sq*zmo)
         elseif(1-yj.lt.tiny)then
            tmp=-(s-xm12)/(4*e0sq)
         else
            eps1=1+xm12/(s-w2)
            beta1=sqrt(eps1**2-4*s*xm12/(s-w2)**2)
            tbeta2=sqrt(1-w2/e0sq)
            call dinvariants_dFKS(ileg,s,x,yi,yj,xm12,xm22,dw1dx,dw1dy,dw2dx,dw2dy)
            tmp=-(dw1dy*dw2dx-dw1dx*dw2dy)*tbeta2/(2*e0sq*z*(1-z)*(s-w2)*beta1)
         endif
c
      else
         write(*,*)'xjacHW6_xiztoxy: unknown ileg'
         stop
      endif
      xjacHW6_xiztoxy=abs(tmp)

      return
 999  continue
      xjacHW6_xiztoxy=0d0
      return
      end



c Hewrig++

      function zHWPP(ileg,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
c Shower energy variable
      implicit none
      integer ileg
      double precision zHWPP,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q,tiny,
     &w1,w2,zeta1,zeta2,get_zeta,betad,betas
      parameter (tiny=1d-5)

      w1=-xq1q+xq2q-xtk
      w2=-xq2q+xq1q-xuk
c
      if(ileg.eq.1)then
         zHWPP=1-(1-x)*(1+yi)/2d0
c
      elseif(ileg.eq.2)then
         zHWPP=1-(1-x)*(1+yi)/2d0
c
      elseif(ileg.eq.3)then
         if(1-x.lt.tiny)then
            betad=sqrt((1-(xm12-xm22)/s)**2-(4*xm22/s))
            betas=1+(xm12-xm22)/s
            zHWPP=1-(1-x)*(1+yj)/(betad+betas)
         else
            zeta1=get_zeta(s,w1,w2,xm12,xm22)
            zHWPP=1-zeta1
         endif
c
      elseif(ileg.eq.4)then
         if(1-x.lt.tiny)then
            zHWPP=1-(1-x)*(1+yj)*s/(2*(s-xm12))
         elseif(1-yj.lt.tiny)then
            zHWPP=(s*x-xm12)/(s-xm12)+(1-yj)*(1-x)*s*(s*x+xm12*(x-2))*
     &                                (s*x-xm12)/(2*(s-xm12)**3)
         else
            zeta2=get_zeta(s,w2,w1,xm22,xm12)
            zHWPP=1-zeta2 
         endif
c
      else
         write(*,*)'zHWPP: unknown ileg'
         stop
      endif

      if(zHWPP.lt.0d0.or.zHWPP.gt.1d0)goto 999

      return
 999  continue
      zHWPP=-1d0
      return
      end



      function xiHWPP(ileg,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
c Shower evolution variable
      implicit none
      integer ileg
      real*8 xiHWPP,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q,tiny,w1,w2,
     &betad,betas,z,zHWPP
      parameter (tiny=1d-5)

      z=zHWPP(ileg,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
      if(z.lt.0d0)goto 999
      w1=-xq1q+xq2q-xtk
      w2=-xq2q+xq1q-xuk
c 
      if(ileg.eq.1)then
         xiHWPP=s*(1-yi)/(1+yi)
c
      elseif(ileg.eq.2)then
         xiHWPP=s*(1-yi)/(1+yi)
c
      elseif(ileg.eq.3)then
         if(1-x.lt.tiny)then
            betad=sqrt((1-(xm12-xm22)/s)**2-(4*xm22/s))
            betas=1+(xm12-xm22)/s
            xiHWPP=-s*(betad+betas)*(yj*betad-betas)/(2*(1+yj))
         else
            xiHWPP=w1/(z*(1-z))
         endif
c
      elseif(ileg.eq.4)then
         if(1-x.lt.tiny)then
            xiHWPP=(1-yj)*(s-xm12)**2/(s*(1+yj))
         elseif(1-yj.lt.tiny)then
            xiHWPP=(1-yj)*(s-xm12)**2/(2*s)
         else
            xiHWPP=w2/(z*(1-z))
         endif
c
      else
         write(*,*)'xiHWPP: unknown ileg'
         stop
      endif

      if(xiHWPP.lt.0d0)goto 999

      return
 999  continue
      xiHWPP=-1d0
      return
      end



      function xjacHWPP_xiztoxy(ileg,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
c Returns the jacobian d(z,xi)/d(x,y), where z and xi are the shower 
c variables, and x and y are FKS variables
      implicit none
      integer ileg
      double precision xjacHWPP_xiztoxy,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,
     &xq2q,tiny,tmp,z,zHWPP,w1,w2,zeta1,dw1dx,dw2dx,dw1dy,dw2dy,get_zeta,
     &betad,betas,eps1,eps2,beta1,beta2
      parameter (tiny=1d-5)

      tmp=0d0
      z=zHWPP(ileg,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
      if(z.lt.0d0)goto 999
      w1=-xq1q+xq2q-xtk
      w2=-xq2q+xq1q-xuk 
c
      if(ileg.eq.1)then
         tmp=-s/(1+yi)
c
      elseif(ileg.eq.2)then
         tmp=-s/(1+yi)
c
      elseif(ileg.eq.3)then
         if(1-x.lt.tiny)then
            betad=sqrt((1-(xm12-xm22)/s)**2-(4*xm22/s))
            betas=1+(xm12-xm22)/s
            tmp=-s*(betad+betas)/(2*(1+yj))
         else
            eps2=1-(xm12-xm22)/(s-w1)
            beta2=sqrt(eps2**2-4*s*xm22/(s-w1)**2)
            call dinvariants_dFKS(ileg,s,x,yi,yj,xm12,xm22,dw1dx,dw1dy,dw2dx,dw2dy)
            tmp=-(dw1dy*dw2dx-dw1dx*dw2dy)/(z*(1-z))/((s-w1)*beta2)
         endif
c
      elseif(ileg.eq.4)then
         if(1-x.lt.tiny)then
            tmp=-(s-xm12)/(1+yj)
         elseif(1-yj.lt.tiny)then
            tmp=-(s-xm12)/2
         else
            eps1=1+xm12/(s-w2)
            beta1=sqrt(eps1**2-4*s*xm12/(s-w2)**2)
            call dinvariants_dFKS(ileg,s,x,yi,yj,xm12,xm22,dw1dx,dw1dy,dw2dx,dw2dy)
            tmp=-(dw1dy*dw2dx-dw1dx*dw2dy)/(z*(1-z))/((s-w2)*beta1)
         endif
c
      else
         write(*,*)'xjacHWPP_xiztoxy: unknown ileg'
         stop
      endif
      xjacHWPP_xiztoxy=abs(tmp)

      return
 999  continue
      xjacHWPP_xiztoxy=0d0
      return
      end



c Pythia6Q

      function zPY6Q(ileg,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
c Shower energy variable
      implicit none
      integer ileg
      double precision zPY6Q,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q,tiny,
     &w1,w2,betad,betas
      parameter(tiny=1d-5)

      w1=-xq1q+xq2q-xtk
      w2=-xq2q+xq1q-xuk
c
      if(ileg.eq.1)then
         zPY6Q=x
c
      elseif(ileg.eq.2)then
         zPY6Q=x
c
      elseif(ileg.eq.3)then
         if(1-x.lt.tiny)then
            betad=sqrt((1-(xm12-xm22)/s)**2-(4*xm22/s))
            betas=1+(xm12-xm22)/s
            zPY6Q=1-(2*xm12)/(s*betas*(betas-betad*yj))
         else
            zPY6Q=1-s*(1-x)*(xm12+w1)/w1/(s+w1+xm12-xm22)
c This is equation (3.10) of hep-ph/1102.3795. In the partonic
c CM frame it is equal to (xk1(0)+xk3(0)*f)/(xk1(0)+xk3(0)),
c where f = xm12/( s+xm12-xm22-2*sqrt(s)*(xk1(0)+xk3(0)) )
         endif
c
      elseif(ileg.eq.4)then
         if(1-x.lt.tiny)then
            zPY6Q=1-s*(1-x)/(s-xm12)
         elseif(1-yj.lt.tiny)then
            zPY6Q=(s*x-xm12)/(s-xm12)+(1-yj)*(1-x)**2*s*(s*x-xm12)/
     &                                ( 2*(s-xm12)**2 )
         else
            zPY6Q=1-s*(1-x)/(s+w2-xm12)
         endif
c
      else
         write(*,*)'zPY6Q: unknown ileg'
         stop
      endif

      if(zPY6Q.lt.0d0.or.zPY6Q.gt.1d0)goto 999

      return
 999  continue
      zPY6Q=-1d0
      return
      end



      function xiPY6Q(ileg,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
c Shower evolution variable
      implicit none
      integer ileg
      double precision xiPY6Q,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q,tiny,z,
     &zPY6Q,w1,w2,betad,betas
      parameter(tiny=1d-5)

      w1=-xq1q+xq2q-xtk
      w2=-xq2q+xq1q-xuk
c
      if(ileg.eq.1)then
         xiPY6Q=s*(1-x)*(1-yi)/2
c
      elseif(ileg.eq.2)then
         xiPY6Q=s*(1-x)*(1-yi)/2
c
      elseif(ileg.eq.3)then
         if(1-x.lt.tiny)then
            betad=sqrt((1-(xm12-xm22)/s)**2-(4*xm22/s))
            betas=1+(xm12-xm22)/s
            xiPY6Q=s*(1-x)*(betas-betad*yj)/2
         else
            xiPY6Q=w1
         endif
c
      elseif(ileg.eq.4)then
         if(1-x.lt.tiny)then
            xiPY6Q=(1-yj)*(1-x)*(s-xm12)/2
         elseif(1-yj.lt.tiny)then
            xiPY6Q=(1-yj)*(1-x)*(s*x-xm12)/2
         else
            xiPY6Q=w2
         endif
c
      else
        write(*,*)'xiPY6Q: unknown ileg'
        stop
      endif

      if(xiPY6Q.lt.0d0)goto 999

      return
 999  continue
      xiPY6Q=-1d0
      return
      end



      function xjacPY6Q_xiztoxy(ileg,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
c Returns the jacobian d(z,xi)/d(x,y), where z and xi are the shower 
c variables, and x and y are FKS variables
      implicit none
      integer ileg
      double precision xjacPY6Q_xiztoxy,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,
     &xq2q,tiny,tmp,z,zPY6Q,w1,w2,dw1dx,dw1dy,dw2dx,dw2dy,betad,betas
      parameter (tiny=1d-5)

      tmp=0d0
      z=zPY6Q(ileg,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
      if(z.lt.0d0)goto 999
      w1=-xq1q+xq2q-xtk
      w2=-xq2q+xq1q-xuk
c
      if(ileg.eq.1)then
         tmp=-s*(1-x)/2
c
      elseif(ileg.eq.2)then
         tmp=-s*(1-x)/2
c
      elseif(ileg.eq.3)then
         if(1-x.lt.tiny)then
            betad=sqrt((1-(xm12-xm22)/s)**2-(4*xm22/s))
            betas=1+(xm12-xm22)/s
            tmp=xm12*betad/betas/(betas-betad*yj)
         else
            call dinvariants_dFKS(ileg,s,x,yi,yj,xm12,xm22,dw1dx,dw1dy,dw2dx,dw2dy)
            tmp=s*(xm12+w1)/w1/(s+w1+xm12-xm22)*dw1dy
         endif
c
      elseif(ileg.eq.4)then
         if(1-x.lt.tiny)then
            tmp=s*(1-x)/2
         elseif(1-yj.lt.tiny)then
            tmp=-s*(1-x)*(s*x-xm12)/( 2*(s-xm12) )
         else
            call dinvariants_dFKS(ileg,s,x,yi,yj,xm12,xm22,dw1dx,dw1dy,dw2dx,dw2dy) 
            tmp=s/(s+w2-xm12)*dw2dy
         endif
c
      else
         write(*,*)'xjacPY6Q_xiztoxy: unknown ileg'
         stop
      endif
      xjacPY6Q_xiztoxy=abs(tmp)

      return
 999  continue
      xjacPY6Q_xiztoxy=0d0
      return
      end



c Pythia6PT

      function zPY6PT(ileg,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
c Shower energy variable
      implicit none
      integer ileg
      double precision zPY6PT,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q,tiny
      parameter(tiny=1d-5)

      if(ileg.eq.1)then
         zPY6PT=x
c
      elseif(ileg.eq.2)then
         zPY6PT=x
c
      elseif(ileg.eq.3)then
         write(*,*)'PYTHIA6PT not available for FSR'
         stop
c
      elseif(ileg.eq.4)then
         write(*,*)'PYTHIA6PT not available for FSR'
         stop
c
      else
         write(*,*)'zPY6PT: unknown ileg'
         stop
      endif

      if(zPY6PT.lt.0d0.or.zPY6PT.gt.1d0)goto 999

      return
 999  continue
      zPY6PT=-1d0
      return
      end



      function xiPY6PT(ileg,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
c Shower evolution variable
      implicit none
      integer ileg
      double precision xiPY6PT,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q,tiny,z
      parameter(tiny=1d-5)

      if(ileg.eq.1)then
         xiPY6PT=s*(1-x)**2*(1-yi)/2
c
      elseif(ileg.eq.2)then
         xiPY6PT=s*(1-x)**2*(1-yi)/2
c
      elseif(ileg.eq.3)then
         write(*,*)'PYTHIA6PT not available for FSR'
         stop
c
      elseif(ileg.eq.4)then
         write(*,*)'PYTHIA6PT not available for FSR'
         stop
c
      else
         write(*,*)'xiPY6PT: unknown ileg'
         stop
      endif

      if(xiPY6PT.lt.0d0)goto 999

      return
 999  continue
      xiPY6PT=-1d0
      return
      end



      function xjacPY6PT_xiztoxy(ileg,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
c Returns the jacobian d(z,xi)/d(x,y), where z and xi are the shower 
c variables, and x and y are FKS variables
      implicit none
      integer ileg
      double precision xjacPY6PT_xiztoxy,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,
     &xq2q,tiny,tmp
      parameter(tiny=1d-5)

      if(ileg.eq.1)then
         tmp=-s*(1-x)**2/2
c
      elseif(ileg.eq.2)then
         tmp=-s*(1-x)**2/2
c
      elseif(ileg.eq.3)then
         write(*,*)'PYTHIA6PT not available for FSR'
         stop
c
      elseif(ileg.eq.4)then
         write(*,*)'PYTHIA6PT not available for FSR'
         stop
c
      else
         write(*,*)'xjacPY6PT_xiztoxy: unknown ileg'
         stop
      endif
      xjacPY6PT_xiztoxy=abs(tmp)

      return
 999  continue
      xjacPY6PT_xiztoxy=0d0
      return
      end



c Pythia8

      function zPY8(ileg,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
c Shower energy variable
      implicit none
      integer ileg
      double precision zPY8,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q,tiny,
     &w1,w2,betad,betas
      parameter(tiny=1d-5)

      w1=-xq1q+xq2q-xtk
      w2=-xq2q+xq1q-xuk
c
      if(ileg.eq.1)then
         zPY8=x
c
      elseif(ileg.eq.2)then
         zPY8=x
c
      elseif(ileg.eq.3)then
         if(1-x.lt.tiny)then
            betad=sqrt((1-(xm12-xm22)/s)**2-(4*xm22/s))
            betas=1+(xm12-xm22)/s
            zPY8=1-(2*xm12)/(s*betas*(betas-betad*yj))
         else
            zPY8=1-s*(1-x)*(xm12+w1)/w1/(s+w1+xm12-xm22)
c This is equation (3.10) of hep-ph/1102.3795. In the partonic
c CM frame it is equal to (xk1(0)+xk3(0)*f)/(xk1(0)+xk3(0)),
c where f = xm12/( s+xm12-xm22-2*sqrt(s)*(xk1(0)+xk3(0)) )
         endif
c
      elseif(ileg.eq.4)then
         if(1-x.lt.tiny)then
            zPY8=1-s*(1-x)/(s-xm12)
         elseif(1-yj.lt.tiny)then
            zPY8=(s*x-xm12)/(s-xm12)+(1-yj)*(1-x)**2*s*(s*x-xm12)/
     &                               ( 2*(s-xm12)**2 )
         else
            zPY8=1-s*(1-x)/(s+w2-xm12)
         endif
c
      else
         write(*,*)'zPY8: unknown ileg'
         stop
      endif

      if(zPY8.lt.0d0.or.zPY8.gt.1d0)goto 999

      return
 999  continue
      zPY8=-1d0
      return
      end



      function xiPY8(ileg,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
c Shower evolution variable
      implicit none
      integer ileg
      double precision xiPY8,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q,tiny,z,
     &zPY8,w1,w2,betas,betad,z0
      parameter(tiny=1d-5)

      z=zPY8(ileg,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
      if(z.lt.0d0)goto 999
      w1=-xq1q+xq2q-xtk
      w2=-xq2q+xq1q-xuk
c
      if(ileg.eq.1)then
         xiPY8=s*(1-x)**2*(1-yi)/2
c
      elseif(ileg.eq.2)then
         xiPY8=s*(1-x)**2*(1-yi)/2
c
      elseif(ileg.eq.3)then
         if(1-x.lt.tiny)then
            betad=sqrt((1-(xm12-xm22)/s)**2-(4*xm22/s))
            betas=1+(xm12-xm22)/s
            z0=1-(2*xm12)/(s*betas*(betas-betad*yj))
            xiPY8=s*(1-x)*(betas-betad*yj)*z0*(1-z0)/2
         else
            xiPY8=z*(1-z)*w1
         endif
c
      elseif(ileg.eq.4)then
         if(1-x.lt.tiny)then
            xiPY8=s*(1-x)**2*(1-yj)/2
         elseif(1-yj.lt.tiny)then
            xiPY8=s*(1-x)**2*(1-yj)*(s*x-xm12)**2/(2*(s-xm12)**2)
         else
            xiPY8=z*(1-z)*w2
         endif
c
      else
        write(*,*)'xiPY8: unknown ileg'
        stop
      endif

      if(xiPY8.lt.0d0)goto 999

      return
 999  continue
      xiPY8=-1d0
      return
      end



      function xjacPY8_xiztoxy(ileg,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
c Returns the jacobian d(z,xi)/d(x,y), where z and xi are the shower 
c variables, and x and y are FKS variables
      implicit none
      integer ileg
      double precision xjacPY8_xiztoxy,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,
     &xq2q,tiny,tmp,z,zPY8,w1,w2,dw1dx,dw1dy,dw2dx,dw2dy,betad,betas,z0
      parameter (tiny=1d-5)

      tmp=0d0
      z=zPY8(ileg,xm12,xm22,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
      if(z.lt.0d0)goto 999
      w1=-xq1q+xq2q-xtk
      w2=-xq2q+xq1q-xuk
c
      if(ileg.eq.1)then
         tmp=-s*(1-x)**2/2
c
      elseif(ileg.eq.2)then
         tmp=-s*(1-x)**2/2
c
      elseif(ileg.eq.3)then
         if(1-x.lt.tiny)then
            betad=sqrt((1-(xm12-xm22)/s)**2-(4*xm22/s))
            betas=1+(xm12-xm22)/s
            z0=1-(2*xm12)/(s*betas*(betas-betad*yj))
            tmp=xm12*betad/betas/(betas-betad*yj)*z0*(1-z0)
         else
            call dinvariants_dFKS(ileg,s,x,yi,yj,xm12,xm22,dw1dx,dw1dy,dw2dx,dw2dy)
            tmp=s*(xm12+w1)/w1/(s+w1+xm12-xm22)*dw1dy*z*(1-z)
         endif
c
      elseif(ileg.eq.4)then
         if(1-x.lt.tiny)then
            tmp=s**2*(1-x)**2/( 2*(s-xm12) )
         elseif(1-yj.le.tiny)then
            tmp=4*s**2*(1-x)**2*(s*x-xm12)**2/( 2*(s-xm12) )**3
         else
            call dinvariants_dFKS(ileg,s,x,yi,yj,xm12,xm22,dw1dx,dw1dy,dw2dx,dw2dy)
            tmp=s/(s+w2-xm12)*dw2dy*z*(1-z)
         endif
c
      else
         write(*,*)'xjacPY8_xiztoxy: unknown ileg'
         stop
      endif
      xjacPY8_xiztoxy=abs(tmp)

      return
 999  continue
      xjacPY8_xiztoxy=0d0
      return
      end

c End of Monte Carlo functions



      function get_zeta(xs,xw1,xw2,xxm12,xxm22)
      implicit none
      double precision get_zeta,xs,xw1,xw2,xxm12,xxm22
      double precision eps,beta
c
      eps=1-(xxm12-xxm22)/(xs-xw1)
      beta=sqrt(eps**2-4*xs*xxm22/(xs-xw1)**2)
      get_zeta=( (2*xs-(xs-xw1)*eps)*xw2+(xs-xw1)*((xw1+xw2)*beta-eps*xw1) )/
     &         ( (xs-xw1)*beta*(2*xs-(xs-xw1)*eps+(xs-xw1)*beta) )
c
      return
      end



      function emscafun(x,alpha)
      implicit none
      double precision emscafun,x,alpha
c
      if(x.lt.0d0.or.x.gt.1d0)then
         write(*,*)'Fatal error in emscafun'
         stop
      endif
      emscafun=x**(2*alpha)/(x**(2*alpha)+(1-x)**(2*alpha))
      return
      end



      function emscainv(r,alpha)
c Inverse of emscafun, implemented only for alpha=1 for the moment
      implicit none
      double precision emscainv,r,alpha
c
      if(r.lt.0d0.or.r.gt.1d0.or.alpha.ne.1d0)then
         write(*,*)'Fatal error in emscafun'
         stop
      endif
      emscainv=sqrt(r)/(sqrt(r)+sqrt(1d0-r))
      return
      end



      function bogus_probne_fun(qMC)
      implicit none
      double precision bogus_probne_fun,qMC
      double precision x,tmp,emscafun
      integer itype
      data itype/2/
c
      if(itype.eq.1)then
c Theta function
         tmp=1d0
         if(qMC.le.2d0)tmp=0d0
      elseif(itype.eq.2)then
c Smooth function
        if(qMC.le.0.5d0)then
          tmp=0d0
        elseif(qMC.le.1d1)then
          x=(1d1-qMC)/(1d1-0.5d0)
          tmp=1-emscafun(x,2d0)
        else
          tmp=1d0
        endif
      else
        write(*,*)'Error in bogus_probne_fun: unknown option',itype
        stop
      endif
      bogus_probne_fun=tmp
      return
      end



      function get_angle(p1,p2)
      implicit none
      double precision get_angle,p1(0:3),p2(0:3)
      double precision tiny,mod1,mod2,cosine
      parameter (tiny=1d-5)
c
      mod1=sqrt(p1(1)**2+p1(2)**2+p1(3)**2)
      mod2=sqrt(p2(1)**2+p2(2)**2+p2(3)**2)

      if(mod1.eq.0d0.or.mod2.eq.0d0)then
         write(*,*)'Undefined angle in get_angle',mod1,mod2
         stop
      endif
c
      cosine=p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3)
      cosine=cosine/(mod1*mod2)
c
      if(abs(cosine).gt.1d0+tiny)then
         write(*,*)'cosine larger than 1 in get_angle',cosine,p1,p2
         stop
      elseif(abs(cosine).ge.1d0)then
         cosine=sign(1d0,cosine)
      endif
c
      get_angle=acos(cosine)

      return
      end



c Shower scale

      subroutine assign_emsca(pp,xi_i_fks,y_ij_fks)
      implicit none
      include "nexternal.inc"
      include "madfks_mcatnlo.inc"
      include "run.inc"

      double precision pp(0:3,nexternal),xi_i_fks,y_ij_fks
      double precision shattmp,dot,emsca_bare,ref_scale,scalemin,
     &scalemax,rrnd,ran2,emscainv,dum(5),xm12,qMC,ptresc
      integer ileg
      double precision p_born(0:3,nexternal-1)
      common/pborn/p_born

      logical emscasharp
      double precision emsca
      common/cemsca/emsca,emsca_bare,emscasharp,scalemin,scalemax

      double precision ybst_til_tolab,ybst_til_tocm,sqrtshat,shat
      common/parton_cms_stuff/ybst_til_tolab,ybst_til_tocm,sqrtshat,shat

c Consistency check
      shattmp=2d0*dot(pp(0,1),pp(0,2))
      if(abs(shattmp/shat-1d0).gt.1d-5)then
         write(*,*)'Error in assign_emsca: inconsistent shat'
         write(*,*)shattmp,shat
         stop
      endif

      call kinematics_driver(xi_i_fks,y_ij_fks,shat,pp,ileg,
     &                       xm12,dum(1),dum(2),dum(3),dum(4),dum(5),qMC,.true.)

      emsca=2d0*sqrt(ebeam(1)*ebeam(2))
      if(dampMCsubt)then
         call assign_scaleminmax(shat,xi_i_fks,scalemin,scalemax,ileg,xm12)
         emscasharp=(scalemax-scalemin).lt.(1d-3*scalemax)
         if(emscasharp)then
            emsca_bare=scalemax
            emsca=emsca_bare
         else
            rrnd=ran2()
            rrnd=emscainv(rrnd,1d0)
            emsca_bare=scalemin+rrnd*(scalemax-scalemin)
            ptresc=(qMC-scalemin)/(scalemax-scalemin)
            if(ptresc.lt.1d0)emsca=emsca_bare
            if(ptresc.ge.1d0)emsca=scalemax
         endif
      endif

      return
      end


      subroutine assign_emsca_array(pp,xi_i_fks,y_ij_fks)
      implicit none
      include "nexternal.inc"
      include "madfks_mcatnlo.inc"
      include "run.inc"
      include "born_nhel.inc"
      double precision pp(0:3,nexternal),xi_i_fks,y_ij_fks,shattmp,dot
      double precision emsca_bare_a(nexternal,nexternal),emsca_bare_a2(nexternal,nexternal),
     &scalemin_a(nexternal,nexternal),scalemax_a(nexternal,nexternal),rrnd,ran2,emscainv,
     &dum(5),xm12,qMC,ptresc_a(nexternal,nexternal),ref_scale_a(nexternal,nexternal)
      integer ileg,npartner,i,j
      double precision p_born(0:3,nexternal-1)
      common/pborn/p_born
      integer ipartners(0:nexternal-1),colorflow(nexternal-1,0:max_bcol)
      common /MC_info/ ipartners,colorflow

      logical emscasharp_a(nexternal,nexternal)
      double precision emsca_a(nexternal,nexternal)
      common/cemsca_a/emsca_a,emsca_bare_a,emsca_bare_a2,
     #                emscasharp_a,scalemin_a,scalemax_a

      double precision ybst_til_tolab,ybst_til_tocm,sqrtshat,shat
      common/parton_cms_stuff/ybst_til_tolab,ybst_til_tocm,sqrtshat,shat

c Consistency check
      shattmp=2d0*dot(pp(0,1),pp(0,2))
      if(abs(shattmp/shat-1d0).gt.1d-5)then
         write(*,*)'Error in assign_emsca_array: inconsistent shat'
         write(*,*)shattmp,shat
         stop
      endif

      call kinematics_driver(xi_i_fks,y_ij_fks,shat,pp,ileg,
     &                       xm12,dum(1),dum(2),dum(3),dum(4),dum(5),qMC,.true.)
      call assign_scaleminmax_array(shat,xi_i_fks,scalemin_a,scalemax_a,ileg,xm12)
      emsca_a=-1d0

      do i=1,nexternal-2
         do j=i+1,nexternal-1
            if(dampMCsubt)then
               emscasharp_a(i,j)=(scalemax_a(i,j)-scalemin_a(i,j)).lt.(1d-3*scalemax_a(i,j))
               if(emscasharp_a(i,j))then
                  emsca_bare_a(i,j)=scalemax_a(i,j)
                  emsca_bare_a2(i,j)=scalemax_a(i,j)
                  emsca_a(i,j)=emsca_bare_a(i,j)
               else
                  rrnd=ran2()
                  rrnd=emscainv(rrnd,1d0)
                  emsca_bare_a(i,j)=scalemin_a(i,j)+rrnd*(scalemax_a(i,j)-scalemin_a(i,j))
                  rrnd=ran2()
                  rrnd=emscainv(rrnd,1d0)
                  emsca_bare_a2(i,j)=scalemin_a(i,j)+rrnd*(scalemax_a(i,j)-scalemin_a(i,j))
                  ptresc_a(i,j)=(qMC-scalemin_a(i,j))/(scalemax_a(i,j)-scalemin_a(i,j))
                  if(ptresc_a(i,j).lt.1d0)emsca_a(i,j)=emsca_bare_a(i,j)
                  if(ptresc_a(i,j).ge.1d0)emsca_a(i,j)=scalemax_a(i,j)
               endif
            endif
            emsca_a(j,i)=emsca_a(i,j)
            emsca_bare_a(j,i)=emsca_bare_a(i,j)
            emsca_bare_a2(j,i)=emsca_bare_a2(i,j)
            emscasharp_a(j,i)=emscasharp_a(i,j)
            ptresc_a(j,i)=ptresc_a(i,j)
         enddo
      enddo

      return
      end



      subroutine assign_scaleminmax(shat,xi,xscalemin,xscalemax,ileg,xm12)
      implicit none
      include "nexternal.inc"
      include "run.inc"
      include "madfks_mcatnlo.inc"
      integer i,ileg
      double precision shat,xi,ref_scale,xscalemax,xscalemin,xm12
      character*4 abrv
      common/to_abrv/abrv
      double precision p_born(0:3,nexternal-1)
      common/pborn/p_born

      call assign_ref_scale(p_born,xi,shat,ref_scale)
      xscalemin=max(shower_scale_factor*frac_low*ref_scale,scaleMClow)
      xscalemax=max(shower_scale_factor*frac_upp*ref_scale,
     &              xscalemin+scaleMCdelta)
      xscalemax=min(xscalemax,2d0*sqrt(ebeam(1)*ebeam(2)))
      xscalemin=min(xscalemin,xscalemax)
c
      if(abrv.ne.'born'.and.shower_mc(1:7).eq.'PYTHIA6'.and.ileg.eq.3)then
         xscalemin=max(xscalemin,sqrt(xm12))
         xscalemax=max(xscalemin,xscalemax)
      endif

      return
      end


      subroutine assign_scaleminmax_array(shat,xi,xscalemin_a,xscalemax_a,ileg,xm12)
      implicit none
      include "nexternal.inc"
      include "run.inc"
      include "madfks_mcatnlo.inc"
      integer i,j,ileg
      double precision shat,xi,ref_scale_a(nexternal,nexternal),xm12
      double precision xscalemax_a(nexternal,nexternal),xscalemin_a(nexternal,nexternal)
      character*4 abrv
      common/to_abrv/abrv
      double precision p_born(0:3,nexternal-1)
      common/pborn/p_born

      call assign_ref_scale_array(p_born,xi,shat,ref_scale_a)
      do i=1,nexternal-2
         do j=i+1,nexternal-1
            xscalemin_a(i,j)=max(shower_scale_factor*frac_low*ref_scale_a(i,j),scaleMClow)
            xscalemax_a(i,j)=max(shower_scale_factor*frac_upp*ref_scale_a(i,j),xscalemin_a(i,j)+scaleMCdelta)
            xscalemax_a(i,j)=min(xscalemax_a(i,j),2d0*sqrt(ebeam(1)*ebeam(2)))
            xscalemin_a(i,j)=min(xscalemin_a(i,j),xscalemax_a(i,j))
c
            if(abrv.ne.'born'.and.shower_mc(1:7).eq.'PYTHIA6'.and.ileg.eq.3)then
               xscalemin_a(i,j)=max(xscalemin_a(i,j),sqrt(xm12))
               xscalemax_a(i,j)=max(xscalemin_a(i,j),xscalemax_a(i,j))
            endif
c
            xscalemin_a(j,i)=xscalemin_a(i,j)
            xscalemax_a(j,i)=xscalemax_a(i,j)
         enddo
      enddo
c
      return
      end


      subroutine assign_ref_scale(p,xii,sh,ref_sc)
      implicit none
      include "nexternal.inc"
      include "madfks_mcatnlo.inc"
      double precision p(0:3,nexternal-1),xii,sh,ref_sc
      integer i_scale,i
      parameter(i_scale=1)

      ref_sc=0d0
      if(i_scale.eq.0)then
c Born-level CM energy squared
         ref_sc=dsqrt(max(0d0,(1-xii)*sh))
      elseif(i_scale.eq.1)then
c Sum of final-state transverse masses
         do i=3,nexternal-1
            ref_sc=ref_sc+dsqrt(max(0d0,(p(0,i)+p(3,i))*(p(0,i)-p(3,i))))
         enddo
         ref_sc=ref_sc/2d0
      else
         write(*,*)'Wrong i_scale in assign_ref_scale',i_scale
         stop
      endif
c Safety threshold for the reference scale
      ref_sc=max(ref_sc,scaleMClow+scaleMCdelta)

      return
      end


      subroutine assign_ref_scale_array(p,xii,sh,ref_sc_a)
      implicit none
      include "nexternal.inc"
      include "madfks_mcatnlo.inc"
      double precision p(0:3,nexternal-1),xii,sh
      double precision ref_sc_a(nexternal,nexternal)
      double precision ref_sc
      integer i,j
c
      call assign_ref_scale(p,xii,sh,ref_sc)
      do i=1,nexternal-2
         do j=i+1,nexternal-1
            ref_sc_a(i,j)=sqrt( max(0d0,(p(0,i)+p(0,j))**2-(p(1,i)+p(1,j))**2
     &                                 -(p(2,i)+p(2,j))**2-(p(3,i)+p(3,j))**2) )
            ref_sc_a(i,j)=min(ref_sc,ref_sc_a(i,j))
            ref_sc_a(i,j)=max(ref_sc_a(i,j),scaleMClow+scaleMCdelta)
            ref_sc_a(j,i)=ref_sc_a(i,j)
         enddo
      enddo
c
      return
      end


      subroutine dinvariants_dFKS(ileg,s,x,yi,yj,xm12,xm22,dw1dx,dw1dy,dw2dx,dw2dy)
c Returns derivatives of Mandelstam invariants with respect to FKS variables
      implicit none
      integer ileg
      double precision s,x,yi,yj,xm12,xm22,dw1dx,dw2dx,dw1dy,dw2dy
      double precision afun,bfun,cfun,mom_fks_sister_p,mom_fks_sister_m,
     &diff_p,diff_m,signfac,dadx,dady,dbdx,dbdy,dcdx,dcdy,mom_fks_sister,
     &dmomfkssisdx,dmomfkssisdy,en_fks,en_fks_sister,dq1cdx,dq2qdx,dq1cdy,
     &dq2qdy
      double precision veckn_ev,veckbarn_ev,xp0jfks
      common/cgenps_fks/veckn_ev,veckbarn_ev,xp0jfks
      double precision tiny
      parameter(tiny=1d-5)

      if(ileg.eq.1)then
         write(*,*)'dinvariants_dFKS should not be called for ileg = 1'
         stop
c
      elseif(ileg.eq.2)then
         write(*,*)'dinvariants_dFKS should not be called for ileg = 2'
         stop
c
      elseif(ileg.eq.3)then
c For ileg = 3, the mother 3-momentum is [afun +- sqrt(bfun) ] / cfun
         afun=sqrt(s)*(1-x)*(xm12-xm22+s*x)*yj
         bfun=s*( (1+x)**2*(xm12**2+(xm22-s*x)**2-
     &        xm12*(2*xm22+s*(1+x**2)))+xm12*s*(1-x**2)**2*yj**2 )
         cfun=s*(-(1+x)**2+(1-x)**2*yj**2)
         dadx=sqrt(s)*yj*(xm22-xm12+s*(1-2*x))
         dady=sqrt(s)*(1-x)*(xm12-xm22+s*x)
         dbdx=2*s*(1+x)*( xm12**2+(xm22-s*x)*(xm22-s*(1+2*x))
     &        -xm12*(2*xm22+s*(1+x+2*(x**2)+2*(1-x)*x*(yj**2))) )
         dbdy=2*xm12*(s**2)*((1-x**2)**2)*yj
         dcdx=-2*s*(1+x+(yj**2)*(1-x))
         dcdy=2*s*((1-x)**2)*yj
c Determine correct sign
         mom_fks_sister_p=(afun+sqrt(bfun))/cfun
         mom_fks_sister_m=(afun-sqrt(bfun))/cfun
         diff_p=abs(mom_fks_sister_p-veckn_ev)
         diff_m=abs(mom_fks_sister_m-veckn_ev)
         if(min(diff_p,diff_m)/max(abs(veckn_ev),1d0).ge.1d-3)then
            write(*,*)'Fatal error 1 in dinvariants_dFKS'
            write(*,*)mom_fks_sister_p,mom_fks_sister_m,veckn_ev
            stop
         elseif(min(diff_p,diff_m)/max(abs(veckn_ev),1d0).ge.tiny)then
            write(*,*)'Numerical imprecision 1 in dinvariants_dFKS'
         endif
         signfac=1d0
         if(diff_p.ge.diff_m)signfac=-1d0
         mom_fks_sister=veckn_ev
         en_fks=sqrt(s)*(1-x)/2
         en_fks_sister=sqrt(mom_fks_sister**2+xm12)
         dmomfkssisdx=(dadx+signfac*dbdx/(2*sqrt(bfun))-dcdx*mom_fks_sister)/cfun
         dmomfkssisdy=(dady+signfac*dbdy/(2*sqrt(bfun))-dcdy*mom_fks_sister)/cfun
         dw1dx=sqrt(s)*( yj*mom_fks_sister-en_fks_sister+(1-x)*
     &                   (mom_fks_sister/en_fks_sister-yj)*dmomfkssisdx )
         dw1dy=-sqrt(s)*(1-x)*( mom_fks_sister+
     &                   (yj-mom_fks_sister/en_fks_sister)*dmomfkssisdy )
         dw2dx=-dw1dx-s
         dw2dy=-dw1dy
c
      elseif(ileg.eq.4)then
         dq1cdx=-(1-yi)*(s*(1+yj)+xm12*(1-yj))/(1+yj+x*(1-yj))**2
         dq2qdx=-(1+yi)*(s*(1+yj)+xm12*(1-yj))/(1+yj+x*(1-yj))**2
         dw2dx=(1-yj)*(s*(1+yj-x*(2*(1+yj)+x*(1-yj)))+2*xm12)/(1+yj+x*(1-yj))**2
         dw1dx=dq1cdx+dq2qdx
         dq1cdy=(1-x)*(1-yi)*(s*x-xm12)/(1+yj+x*(1-yj))**2
         dq2qdy=(1-x)*(1+yi)*(s*x-xm12)/(1+yj+x*(1-yj))**2
         dw2dy=-2*(1-x)*(s*x-xm12)/(1+yj+x*(1-yj))**2
         dw1dy=dq1cdy+dq2qdy
c
      else
         write(*,*)'Error in dinvariants_dFKS: unknown ileg',ileg
         stop
      endif

      return
      end



      subroutine get_dead_zone(ileg,z,xi,s,x,yi,xm12,xm22,w1,w2,qMC,
     &                         scalemax,ip,ifat,lzone,wcc)
      implicit none
      include "run.inc"
      include "nexternal.inc"
      include 'nFKSconfigs.inc'
      include "madfks_mcatnlo.inc"

      integer ileg,ip,ifat,i
      double precision z,xi,s,x,yi,xm12,xm22,w1,w2,qMC,scalemax,wcc
      logical lzone

      double precision max_scale,upscale,upscale2,xmp2,xmm2,xmr2,ww,Q2,
     &lambda,dot,e0sq,beta,dum,ycc,mdip,mdip_g,zp1,zm1,zp2,zm2,zp3,zm3
      external dot

      double precision p_born(0:3,nexternal-1)
      common/pborn/p_born
      double precision pip(0:3),pifat(0:3),psum(0:3)

      INTEGER NFKSPROCESS
      COMMON/C_NFKSPROCESS/NFKSPROCESS
      double precision shower_S_scale(fks_configs*2)
     &     ,shower_H_scale(fks_configs*2),ref_H_scale(fks_configs*2)
     &     ,pt_hardness
      common /cshowerscale2/shower_S_scale,shower_H_scale,ref_H_scale
     &     ,pt_hardness

      integer mstj50,mstp67
      double precision parp67
      parameter (mstj50=2,mstp67=2,parp67=1d0)
      double precision get_angle,theta2p

c Stop if wrong ileg
      if(ileg.lt.1.or.ileg.gt.4)then
         write(*,*)'Subroutine get_dead_zone: unknown ileg ',ileg
         stop
      endif

c Check compatibility among masses and ileg
      if( (ileg.eq.3.and.xm12.eq.0d0) .or.
     &    (ileg.eq.4.and.xm22.ne.0d0) .or.
     &    (ileg.le.2.and.(xm12.ne.0d0.or.xm22.ne.0d0)) )then
         write(*,*)'Subroutine get_dead_zone: wrong masses '
         write(*,*)ileg,sqrt(xm12),sqrt(xm22)
         stop
      endif

c Skip if unphysical shower variables
      if(z.lt.0d0.or.xi.lt.0d0)goto 999

c Definition and initialisation of variables
      lzone=.true.
      do i=0,3
         pifat(i)=p_born(i,ifat)
         pip(i)  =p_born(i,ip)
         psum(i) =pifat(i)+pip(i) 
      enddo
      max_scale=scalemax
      xmp2=dot(pip,pip)
      e0sq=dot(pip,pifat)
      theta2p=get_angle(pip,pifat)
      theta2p=theta2p**2
      xmm2=xm12*(4-ileg)
      xmr2=xm22*(4-ileg)-xm12*(3-ileg)
      ww=w1*(4-ileg)-w2*(3-ileg)
      Q2=dot(psum,psum)
      lambda=sqrt((Q2+xmm2-xmp2)**2-4*Q2*xmm2)
      beta=sqrt(1-4*s*(xmm2+ww)/(s-xmr2+xmm2+ww)**2)
      wcc=1d0
      ycc=1-parp67*x/(1-x)**2/2
      mdip  =sqrt((sqrt(xmp2+xmm2+2*e0sq)-sqrt(xmp2))**2-xmm2)
      mdip_g=sqrt((sqrt(s)-sqrt(xmr2))**2-xmm2)
      zp1=(1+(xmm2+beta*ww)/(xmm2+ww))/2
      zm1=(1+(xmm2-beta*ww)/(xmm2+ww))/2
      zp2=(1+beta)/2
      zm2=(1-beta)/2
      zp3=(1+sqrt(1-4*xi/mdip_g**2))/2
      zm3=(1-sqrt(1-4*xi/mdip_g**2))/2

c Dead zones
c IMPLEMENT QED DZ's!
      if(shower_mc.eq.'HERWIG6')then
         lzone=.false.
         if(ileg.le.2.and.z**2.ge.xi)lzone=.true.
         if(ileg.gt.2.and.e0sq*xi*z**2.ge.xmm2
     &               .and.xi.le.1d0)lzone=.true.
         if(e0sq.eq.0d0)lzone=.false.
c
      elseif(shower_mc.eq.'HERWIGPP')then
         lzone=.false.
         if(ileg.le.2)upscale2=2*e0sq
         if(ileg.gt.2)then
            upscale2=2*e0sq+xmm2
            if(ip.gt.2)upscale2=(Q2+xmm2-xmp2+lambda)/2
         endif
         if(xi.lt.upscale2)lzone=.true.
c
      elseif(shower_mc.eq.'PYTHIA6Q')then
         if(ileg.le.2)then
            if(mstp67.eq.2.and.ip.gt.2.and.
     &         4*xi/s/(1-z).ge.theta2p)lzone=.false.
         elseif(ileg.gt.2)then
            if(mstj50.eq.2.and.ip.le.2.and.
c around line 71636 of pythia6428: V(IEP(1),5)=virtuality, P(IM,4)=sqrt(s)
     &         max(z/(1-z),(1-z)/z)*4*(xi+xmm2)/s.ge.theta2p)lzone=.false.
            if(z.gt.zp1.or.z.lt.zm1)lzone=.false.
         endif
         if(.not.dampMCsubt)then
            call assign_scaleminmax(s,1-x,dum,upscale,ileg,xmm2)
            if(ileg.gt.2)upscale=max(upscale,sqrt(xmm2))
            if(qMC.gt.upscale)lzone=.false.
         endif
c
      elseif(shower_mc.eq.'PYTHIA6PT')then
         if(mstp67.eq.1.and.yi.lt.ycc)lzone=.false.
         if(mstp67.eq.2)wcc=min(1d0,(1-ycc)/(1-yi))
         if(.not.dampMCsubt)then
            call assign_scaleminmax(s,1-x,dum,upscale,ileg,xmm2)
            if(qMC.gt.upscale)lzone=.false.
         endif
c
      elseif(shower_mc.eq.'PYTHIA8')then
         if(ileg.le.2.and.z.gt.1-sqrt(xi/z/s)*
     &      (sqrt(1+xi/4/z/s)-sqrt(xi/4/z/s)))lzone=.false.
         if(ileg.gt.2)then
            max_scale=min(min(scalemax,mdip/2),mdip_g/2)
            if(z.gt.min(zp2,zp3).or.z.lt.max(zm2,zm3))lzone=.false.
         endif
         if(.not.dampMCsubt)then
            call assign_scaleminmax(s,1-x,dum,upscale,ileg,xmm2)
            if(qMC.gt.upscale)lzone=.false.
         endif

      endif
 
      max_scale=min(max_scale,shower_S_scale(nFKSprocess*2-1))
      max_scale=max(max_scale,scaleMCcut)
      if(qMC.gt.max_scale)lzone=.false.

      return
 999  continue
      lzone=.false.
      return
      end



      function charge(ipdg)
c computes the electric charge given the pdg code
      implicit none
      integer ipdg
      double precision charge,tmp,dipdg

      dipdg=dble(ipdg)
c quarks
      if(abs(dipdg).eq.1) tmp=-1d0/3d0*sign(1d0,dipdg)
      if(abs(dipdg).eq.2) tmp= 2d0/3d0*sign(1d0,dipdg)
      if(abs(dipdg).eq.3) tmp=-1d0/3d0*sign(1d0,dipdg)
      if(abs(dipdg).eq.4) tmp= 2d0/3d0*sign(1d0,dipdg)
      if(abs(dipdg).eq.5) tmp=-1d0/3d0*sign(1d0,dipdg)
      if(abs(dipdg).eq.6) tmp= 2d0/3d0*sign(1d0,dipdg)
c leptons
      if(abs(dipdg).eq.11)tmp=-1d0*sign(1d0,dipdg)
      if(abs(dipdg).eq.12)tmp= 0d0
      if(abs(dipdg).eq.13)tmp=-1d0*sign(1d0,dipdg)
      if(abs(dipdg).eq.14)tmp= 0d0
      if(abs(dipdg).eq.15)tmp=-1d0*sign(1d0,dipdg)
      if(abs(dipdg).eq.16)tmp= 0d0
c bosons
      if(dipdg.eq.21)     tmp= 0d0
      if(dipdg.eq.22)     tmp= 0d0
      if(dipdg.eq.23)     tmp= 0d0
      if(abs(dipdg).eq.24)tmp= 1d0*sign(1d0,dipdg)
      if(dipdg.eq.25)     tmp= 0d0
c
      charge=tmp

      return
      end
