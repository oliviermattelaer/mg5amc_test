      SUBROUTINE GENERATE_CF_INTEGRATION_BASIS()
      IMPLICIT NONE
      INTEGER NDIAGS
      PARAMETER (NDIAGS=12)
      INTEGER    NPERMS
      PARAMETER (NPERMS=24)
      
      LOGICAL CF_BASIS(NDIAGS,NPERMS)
      COMMON/TO_CF_BASIS/CF_BASIS
      
C     CONSTANTS
C     GLOBAL VARIABLES
      include 'maxconfigs.inc'
      include 'maxparticles.inc'
C     ARGUMENTS
C     LOCAL VARIABLES 
      integer i,j,k
      integer natm,max_ch,max_chcf,ierr
      parameter (natm=int((max_particles+1)/2d0))
      parameter (max_ch=10000,max_chcf=100)
      integer atm(2,natm),cf(max_particles,max_ch)
      integer cfnum(max_chcf,max_ch)
      integer ndiag,ngluons,ncf,nchcf(max_ch)
      integer ch_cf(max_particles,max_chcf,max_ch)
      logical KEEP(max_ch,10000)
C     EXTERNAL FUNCTIONS
      integer fact
      external fact
C     ----------
C     BEGIN CODE
C     ----------
c      ndiag = lmaxconfigs
      ngluons = max_particles
c      ncf = fact(ngluons -1)
      do i = 1,ndiags
         do j = 1,nperms
            CF_BASIS(i,j) = .false.
         enddo
      enddo         
      
      do i = 1,ndiags
         call make_atm(i,atm,natm)
         call find_perchcf(ngluons,natm,atm,ch_cf(1,1,i),cfnum(1,i),nchcf(i))

         do j = 1,nchcf(i)
            do k = 1,nperms
               if (cfnum(j,i).eq.k) then
                  CF_BASIS(i,k) = .true.
                  exit
               endif
            enddo
         enddo
      enddo

c      open(1,file="chcf_table.dat",status='replace')
c      do i = 1,nperms
c         write(1,*) i,(KEEP(j,i),j=1,ndiags)
c      enddo
c      close(1)
c      call write_output2(ngluons,ch_cf,nchcf,ndiag,nperms,KEEP)

      end


      subroutine make_atm(ich,atm_out,max_natm)
C     ****************************************************
C     By Yoshitaro Takaesu @U.tokyo Jun.25 2014
C     ****************************************************
      implicitnone
C     CONSTANTS     
      include 'maxconfigs.inc'
      include 'genps.inc'
      include 'maxamps.inc'
C     ARGUMENTS 
      integer max_natm,ich
      integer atm_out(2,max_natm)
C     GLOBAL VARIABLES
      integer iforest(2,-max_branch:-1,lmaxconfigs)
      common/to_forest/ iforest
      integer sprop(maxsproc,-max_branch:-1,lmaxconfigs)
      integer tprid(-max_branch:-1,lmaxconfigs)
      common/to_sprop/sprop,tprid
      integer            mapconfig(0:lmaxconfigs), this_config
      common/to_mconfigs/mapconfig, this_config
C     LOCAL VARIABLES
      integer i
c      include 'configs.inc'
      integer atm(2,20),iprop(max_natm),nt,ns
      integer iused(max_particles),ipost,iposs,idiag,next,nfin,ipos
      integer maxnprop,ipos_tmp,nn,ipos_tmp2,atm1(10),iatm1(10)
      integer prop2_1,prop2_2
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
C Initialization
      idiag = ich
      next = max_particles
      nfin = next -2
      maxnprop = -1*nfin ! the numbering of the external line of particle 2
      ipos_tmp = 10  ! temporary atm number for s-splitting atm attached to particle 2
      atm(1,1) = 1
      atm(1,2) = 2

C inerpriting iforest information into atm information
      prop2_1 = iforest(1,maxnprop,idiag) 
      prop2_2 = iforest(2,maxnprop,idiag) 
C store a splitting leg attached to the particle 2
      if ((prop2_1.lt.0).and.(SPROP(1,MIN(prop2_1,-1),idiag).ne.0)) then
         atm(1,ipos_tmp) = iforest(1,prop2_1,idiag) 
         atm(2,ipos_tmp) = iforest(2,prop2_1,idiag) 
         atm(2,2) = 0
      elseif ((prop2_2.lt.0).and.(SPROP(1,MIN(prop2_2,-1),idiag).ne.0)) then
         write(*,*) ipos_tmp, prop2_2, idiag
         atm(1,ipos_tmp) = iforest(1,prop2_2,idiag) 
         atm(2,ipos_tmp) = iforest(2,prop2_2,idiag) 
         atm(2,2) = 0
C store an external attached to the particle 2
      else
         if (iforest(2,maxnprop,idiag).gt.0) then
            atm(2,2) = iforest(2,maxnprop,idiag)
         else
            atm(2,2) = 0
         endif
      endif

      ipos = 2
      ipos_tmp2 = 10
      do i = maxnprop+1,-1,1
C check atm1
         if (iforest(1,i,idiag).eq.1) then 
C store an external leg attarched to the particle 1
            if (iforest(2,i,idiag).gt.0) then
               atm(2,1) = iforest(2,i,idiag)
               iprop(1) = 0
C store a splitting leg attached to the particle 1
            else
               ipos = ipos +1
               atm(2,1) = 0
               iprop(1) = iforest(2,i,idiag)
               atm(1,ipos) = iforest(1,iprop(1),idiag) 
               atm(2,ipos) = iforest(2,iprop(1),idiag) 
            endif
C check the rest of external legs
         elseif (iforest(2,i,idiag).gt.0) then 
C store an external leg attached to a propagator
            if (iforest(1,i,idiag).lt.0) then 
               ipos_tmp2 = ipos_tmp2 +1
               atm(1,ipos_tmp2) = iforest(2,i,idiag)
               atm(2,ipos_tmp2) = 0
C store a splitting leg attached to a propagator
            elseif (iforest(1,i,idiag).gt.0) then 
               ipos_tmp2 = ipos_tmp2 +1
               atm(1,ipos_tmp2) = iforest(1,i,idiag)
               atm(2,ipos_tmp2) = iforest(2,i,idiag)
            endif
         endif         
      enddo

C sort the atm's ascendent order w.r.t the first element, atm(1,*)
      nn = ipos_tmp2 -10
      do i = 1,nn
         atm1(i) = atm(1,i+10)
      enddo
      call isort_asc(nn,atm1,iatm1)
      do i = 1,nn
         ipos = ipos +1
         atm(1,ipos) = atm(1,iatm1(i)+10)
         atm(2,ipos) = atm(2,iatm1(i)+10)
      enddo         
      ipos = ipos +1
      atm(1,ipos) = atm(1,ipos_tmp)
      atm(2,ipos) = atm(2,ipos_tmp)

      do i = 1,max_natm
         atm_out(1,i) = atm(1,i)
         atm_out(2,i) = atm(2,i)
      enddo

      return
      end


      subroutine isort_asc(n, a,m)
c by Yoshitaro Takaesu -Jun/25/2014 @U.Toyo
      implicit none
      integer i,j
      integer n
      integer m(n)
      integer a(n)
      integer mini,tmp
      integer min
      do i = 1,n
         m(i) = i
      enddo
      do i = 1,n-1
         mini = i
         min = a(i)
         do j = i+1,n
            if( a(j) < min ) then
               min = a(j)
               mini = j
            endif
         enddo
         if( mini.ne.i ) then
            a(mini) = a(i)
            a(i) = min
            tmp = m(mini)
            m(mini) = m(i)
            m(i) = tmp
         endif
      enddo
      end

      subroutine find_chcf(next,natm,atm,cflow,ncf)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Nov.19 2012
C     ****************************************************
      implicitnone
C     CONSTANTS
C     ARGUMENTS 
      integer next,natm
      integer atm(2,natm),cflow(next,100000),ncf
C     GLOBAL VARIABLES
C     LOCAL VARIABLES 
      integer i,j,k
      integer iflag,nperm,pos2,ierr,a(natm),nncf,maxk(natm)
      integer cflowtmp
C     EXTERNAL FUNCTIONS
      integer fact
      external fact
C     ----------
C     BEGIN CODE
C     ----------
      iflag = 0
      do i = 1,natm
         a(i) = i
      enddo
      do i = 1,natm
         maxk(i) = 1
         if (atm(2,i).ne.0) maxk(i) = 2
      enddo
       
      ncf = 0
      nperm = fact(natm-1)
      do i = 1,nperm
         if (i.ne.1) then
            call ipnext(a(2),natm-1,iflag)
         endif
         do j = 2,natm
            if (a(j).eq.2) pos2 = j
         enddo
         ierr = 0
         do j = 2,pos2-2
            do k = j+1,pos2-1
               if (a(j).lt.a(k)) ierr = 1
            enddo
         enddo
         if (ierr.ne.0) cycle
         ierr = 0
         do j = pos2+1,natm-1
            do k = j+1,natm
               if (a(j).gt.a(k)) ierr = 1
            enddo
         enddo
         if (ierr.ne.0) cycle

         call get_cflow(next,natm,atm,a,maxk,nncf,cflow(1,ncf+1))
         ncf = ncf +nncf
      enddo

      do i = 1,ncf
         if (cflow(1,i).ne.1) then
            cflowtmp = cflow(1,i)
            do j = 1,next-1
               cflow(j,i) = cflow(j+1,i)
            enddo
            cflow(next,i) = cflowtmp
         endif
      enddo

      return
      end


      subroutine find_perchcf(next,natm,atm,cflow,cflow_num,ncf)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN.14 2013
C     ****************************************************
      implicitnone
C     CONSTANTS
C     ARGUMENTS 
      integer next,natm,ierr
      integer atm(2,natm),cflow(next,10000),ncf,cflow_num(10000)
C     GLOBAL VARIABLES
C     LOCAL VARIABLES 
      integer i,j,k
      integer iflag,nperm,a(natm),nonzero,zero,nonzero_atm(2,natm)
      integer zero_atm(2,natm),incf,iatm,aatm(2,natm),aaatm(2,natm)
C     EXTERNAL FUNCTIONS
      integer fact
      external fact
C     ----------
C     BEGIN CODE
C     ----------
      iflag = 0
      ncf = 0
      ierr = 0
      do i = 3,natm
         a(i) = i
      enddo
      do i = 1,natm
         aatm(1,i) = atm(1,i)
         aatm(2,i) = atm(2,i)
      enddo

      if ((atm(2,1)*atm(2,2)).ne.0) then
         nperm = fact(natm-2)
         do i = 1,nperm
            if (i.ne.1) then
               call ipnext(a(3),natm-2,iflag)
            endif
            do j = 3,natm
               aatm(1,j) = atm(1,a(j))
               aatm(2,j) = atm(2,a(j))
            enddo
            call find_chcf(next,natm,aatm,cflow(1,ncf+1),incf)
            ncf = ncf +incf
         enddo
      elseif ((atm(2,1).eq.0).and.(atm(2,2).ne.0)) then
         nperm = fact(natm-3)
         do i = 1,nperm
            if (i.ne.1) then
               call ipnext(a(3),natm-3,iflag)
            endif
            do j = 3,natm-1
               aatm(1,j) = atm(1,a(j))
               aatm(2,j) = atm(2,a(j))
            enddo
            call find_chcf(next,natm,aatm,cflow(1,ncf+1),incf)
            ncf = ncf +incf
         enddo
      elseif ((atm(2,1).ne.0).and.(atm(2,2).eq.0)) then
         nperm = fact(natm-3)
         do i = 1,nperm
            if (i.ne.1) then
               call ipnext(a(4),natm-3,iflag)
            endif
            do j = 4,natm
               aatm(1,j) = atm(1,a(j))
               aatm(2,j) = atm(2,a(j))
            enddo
            call find_chcf(next,natm,aatm,cflow(1,ncf+1),incf)
            ncf = ncf +incf
         enddo
      elseif ((atm(2,1).eq.0).and.(atm(2,2).eq.0)) then
         nperm = fact(natm-4)
         do i = 1,nperm
            if (i.ne.1) then
               call ipnext(a(4),natm-4,iflag)
            endif
            do j = 4,natm-1
               aatm(1,j) = atm(1,a(j))
               aatm(2,j) = atm(2,a(j))
            enddo
            call find_chcf(next,natm,aatm,cflow(1,ncf+1),incf)
            ncf = ncf +incf
         enddo
         do i = 1,natm
            aatm(1,i) = atm(1,i)
            aatm(2,i) = atm(2,i)
         enddo
         aatm(1,3) = atm(1,natm)
         aatm(2,3) = atm(2,natm)
         aatm(1,natm) = atm(1,3)
         aatm(2,natm) = atm(2,3)

         do i = 1,natm
            aaatm(1,i) = aatm(1,i)
            aaatm(2,i) = aatm(2,i)
         enddo
         iflag = 0
         do i = 1,natm
            a(i) = i
         enddo
         do i = 1,nperm
            if (i.ne.1) then
               call ipnext(a(4),natm-4,iflag)
            endif
            do j = 4,natm-1
               aaatm(1,j) = aatm(1,a(j))
               aaatm(2,j) = aatm(2,a(j))
            enddo
            call find_chcf(next,natm,aaatm,cflow(1,ncf+1),incf)
            ncf = ncf +incf
         enddo
      endif

C     Get color-flow numbers of each channel
      do i = 1,ncf
         call find_cfnumber(next,cflow(1,i),cflow_num(i))
         if (cflow_num(i).eq.-1) then
            write(*,*) "ERROR: Wrong color-flow is generated."
            return
         endif
      enddo

      return
      end


      subroutine find_cfnumber(next,cflow,cflow_num)
C     ****************************************************
C     By Yoshitaro Takaesu @ U.Tokyo JUN.27 2014
C     ****************************************************
      implicitnone
C     CONSTANTS
C     ARGUMENTS 
      integer next,cflow(next),cflow_num
C     GLOBAL VARIABLES
C     LOCAL VARIABLES 
      integer i
      integer floworder(next),iflag1,iflag2,ncf
C     EXTERNAL FUNCTIONS
      integer fact
      external fact
C     ----------
C     BEGIN CODE
C     ----------
C     Initialization
      ncf = fact(next-1)
      do i = 1,next
         floworder(i) = i
      enddo
C     Finding the color-flow number of cflow
      do i = 1,ncf
         if (i.ne.1) then
            call ipnext(floworder(2),next-1,iflag1)
         endif
         call check_cf(next,floworder,cflow,iflag2)
C     If the color-flow number is found, return the number
         if (iflag2.eq.1) then
            cflow_num = i 
            return
         endif
      enddo
C     If the color-flow number is not found, return -1
      cflow_num = -1

      return
      end



      subroutine check_cf(next,cf1,cf2,iflag)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Dec.8 2012
C     ****************************************************
      implicitnone
C     ARGUMENTS
      integer next
      integer cf1(next),cf2(next),iflag
C     LOCAL VARIABLES 
      integer i
C     ----------
C     BEGIN CODE
C     ----------
      iflag = 0
      do i = 1,next
         if (cf1(i).ne.cf2(i)) return
      enddo
      iflag = 1
      return
      end

      integer function fact(n)
      implicit none
      integer n,i      
      fact = 1
      if (n.eq.0) then
         fact = 1
      else
         do i = 1, n
            fact = fact*i
         enddo
      endif       
      return
      end


      subroutine get_cflow(next,natm,atm,a,maxk,ncf,cflow)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS AUG 25 2012
C     ****************************************************
      implicit none

C     GLOBAL VARIABLES

C     CONSTANTS

C     ARGUMENTS 
      integer natm,next
      integer a(natm),atm(2,natm),maxk(natm),ncf
      integer cflow(next,1000)
C     LOCAL VARIABLES 
      integer i,j,cfpos,k1,k2,k3,k4,k5,k6,k7,k8
      integer aatm(2,natm)
C     EXTERNAL FUNCTIONS

C     ----------
C     BEGIN CODE
C     ----------
      ncf = 0

      if (natm.eq.3) then
         do k1 = 1,maxk(1)
            do k2 = 1,maxk(2)
               do k3 = 1,maxk(3)
                  do j = 1,2
                     aatm(j,1) = atm(mod(j*k1,3),1)
                     aatm(j,2) = atm(mod(j*k2,3),2)
                     aatm(j,3) = atm(mod(j*k3,3),3)
                  enddo
                  ncf = ncf +1
                  call sort_cflow(aatm,a,natm,next,cflow(1,ncf))
               enddo 
            enddo
         enddo
      elseif (natm.eq.4) then
         do k1 = 1,maxk(1)
            do k2 = 1,maxk(2)
               do k3 = 1,maxk(3)
                  do k4 = 1,maxk(4)
                     do j = 1,2
                        aatm(j,1) = atm(mod(j*k1,3),1)
                        aatm(j,2) = atm(mod(j*k2,3),2)
                        aatm(j,3) = atm(mod(j*k3,3),3)
                        aatm(j,4) = atm(mod(j*k4,3),4)
                     enddo
                     ncf = ncf +1
                     call sort_cflow(aatm,a,natm,next,cflow(1,ncf))
                  enddo 
               enddo
            enddo
         enddo
      elseif (natm.eq.5) then
         do k1 = 1,maxk(1)
            do k2 = 1,maxk(2)
               do k3 = 1,maxk(3)
                  do k4 = 1,maxk(4)
                     do k5 = 1,maxk(5)
                        do j = 1,2
                           aatm(j,1) = atm(mod(j*k1,3),1)
                           aatm(j,2) = atm(mod(j*k2,3),2)
                           aatm(j,3) = atm(mod(j*k3,3),3)
                           aatm(j,4) = atm(mod(j*k4,3),4)
                           aatm(j,5) = atm(mod(j*k5,3),5)
                        enddo
                        ncf = ncf +1
                        call sort_cflow(aatm,a,natm,next,cflow(1,ncf))
                     enddo
                  enddo 
               enddo
            enddo
         enddo      
      elseif (natm.eq.6) then
         do k1 = 1,maxk(1)
            do k2 = 1,maxk(2)
               do k3 = 1,maxk(3)
                  do k4 = 1,maxk(4)
                     do k5 = 1,maxk(5)
                        do k6 = 1,maxk(6)
                           do j = 1,2
                              aatm(j,1) = atm(mod(j*k1,3),1)
                              aatm(j,2) = atm(mod(j*k2,3),2)
                              aatm(j,3) = atm(mod(j*k3,3),3)
                              aatm(j,4) = atm(mod(j*k4,3),4)
                              aatm(j,5) = atm(mod(j*k5,3),5)
                              aatm(j,6) = atm(mod(j*k6,3),6)
                           enddo
                           ncf = ncf +1
                          call sort_cflow(aatm,a,natm,next,cflow(1,ncf))
                        enddo
                     enddo
                  enddo 
               enddo
            enddo
         enddo      
      elseif (natm.eq.7) then
         do k1 = 1,maxk(1)
            do k2 = 1,maxk(2)
               do k3 = 1,maxk(3)
                  do k4 = 1,maxk(4)
                     do k5 = 1,maxk(5)
                        do k6 = 1,maxk(6)
                           do k7 = 1,maxk(7)
                              do j = 1,2
                                 aatm(j,1) = atm(mod(j*k1,3),1)
                                 aatm(j,2) = atm(mod(j*k2,3),2)
                                 aatm(j,3) = atm(mod(j*k3,3),3)
                                 aatm(j,4) = atm(mod(j*k4,3),4)
                                 aatm(j,5) = atm(mod(j*k5,3),5)
                                 aatm(j,6) = atm(mod(j*k6,3),6)
                                 aatm(j,7) = atm(mod(j*k7,3),7)
                              enddo
                              ncf = ncf +1
                          call sort_cflow(aatm,a,natm,next,cflow(1,ncf))
                       enddo
                    enddo
                 enddo
              enddo 
           enddo
        enddo
      enddo      
      elseif (natm.eq.8) then
         do k1 = 1,maxk(1)
            do k2 = 1,maxk(2)
               do k3 = 1,maxk(3)
                  do k4 = 1,maxk(4)
                     do k5 = 1,maxk(5)
                        do k6 = 1,maxk(6)
                           do k7 = 1,maxk(7)
                              do k8 = 1,maxk(8)
                                 do j = 1,2
                                    aatm(j,1) = atm(mod(j*k1,3),1)
                                    aatm(j,2) = atm(mod(j*k2,3),2)
                                    aatm(j,3) = atm(mod(j*k3,3),3)
                                    aatm(j,4) = atm(mod(j*k4,3),4)
                                    aatm(j,5) = atm(mod(j*k5,3),5)
                                    aatm(j,6) = atm(mod(j*k6,3),6)
                                    aatm(j,7) = atm(mod(j*k7,3),7)
                                    aatm(j,8) = atm(mod(j*k8,3),8)
                                 enddo
                                 ncf = ncf +1
                          call sort_cflow(aatm,a,natm,next,cflow(1,ncf))
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo 
               enddo
            enddo
         enddo
      endif

      return
      end
      

      subroutine sort_cflow(atm,ia,natm,next,cflow)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS AUG 25 2012
C
C     
C     ****************************************************
      implicit none

C     GLOBAL VARIABLES

C     CONSTANTS

C     ARGUMENTS 
      integer natm,next
      integer ia(natm),cflow(next),atm(2,natm)
C     LOCAL VARIABLES 
      integer i,j,cfpos
C     EXTERNAL FUNCTIONS

C     ----------
C     BEGIN CODE
C     ----------
      cfpos = 0
      do i = 1,natm
         do j = 1,2
            if (atm(j,ia(i)).ne.0) then
               cfpos = cfpos +1
               cflow(cfpos) = atm(j,ia(i))
            endif
         enddo
      enddo

      return
      end

