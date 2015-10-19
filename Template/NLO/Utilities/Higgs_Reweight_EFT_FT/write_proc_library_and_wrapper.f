      program write_proc_labels_virtual
      implicit none
      character*200 procid,str_j,str_mult,str_temp,library,line
     &     ,str_permut,str_i,incoming,incoming1,incoming2
     &     ,str_init1,str_init2
      character*1 s1
      character*2 s2
      character*3 s3
      character*60 buff,model
      integer maxproc,maxpart,maxchannels
      parameter (maxproc=999,maxpart=20,maxchannels=999)
      integer total_proc,i,j,proc_number(maxproc),npart(maxproc),
     &     partid(maxpart,maxproc),lstr(maxproc),l,iproc,mult
      character*140 str(maxproc),initial_states(maxchannels)
      character*2 ls
      logical need_switching,found,firsttime(maxchannels)
      integer date(3),time(3),len_line
      integer flst_ml5(maxproc),k,alpha(maxproc),alphas(maxproc)
      integer n,permut(10),ident,ident_max
      logical nextp
      external nextp

      model='sm-no_b_mass'
      library='ML5lib_reweight'

      open (unit=1,file='processes.dat',status='old',err=33)
      i = 0
      iproc = 0
 100  i = i + 1
      read(1,*,END=101) line
      if ((line(1:1).eq.' ').or.(line(1:1).eq.'#')) then
         GOTO 100
      endif
      iproc = iproc+1
      backspace(1)
      read(1,*,err=33) proc_number(iproc),npart(iproc),
     &     (partid(j,iproc),j=1,npart(iproc)),alpha(iproc),alphas(iproc)
      total_proc=proc_number(iproc)
      GOTO 100
 101  close(1)

com-- this is to split the wrapper function with respect to the 
com-- different initial states
      ident = 0
      
      do j=1,total_proc
         write (incoming1,*) partid(1,j)
         write (incoming2,*) partid(2,j)
         do i=1,200
            if (incoming1(i:i) .eq. "-") incoming1(i:i) = "m"
            if (incoming2(i:i) .eq. "-") incoming2(i:i) = "m"
         enddo
        incoming=trim(adjustl(incoming1))//"_"//trim(adjustl(incoming2))
!         print*, trim(incoming)
         if ( .not. ANY( initial_states 
     &       .eq. trim(adjustl(incoming)) ) ) then
            ident = ident+1
            initial_states(ident)=trim(adjustl(incoming))
            firsttime(ident)=.true.
         endif
      enddo
      ident_max = ident

c$$$      do ident=10,ident_max
c$$$         print*, ident,trim(adjustl(initial_states(ident)))
c$$$      enddo

      open(unit=3,file='MadLoop.mg5',err=35)
      write(*,*) '************************************************'//
     $     '******'
      write(*,*) 'Writing the orderfile for MadLoop5...'
      write(*,*) 'and a wrapper function to access the library...'
      write(*,*) ''
       

      write(3,*) '# MadLoop.mg5'
      write(3,*) '# Created for heft-mergin + reweighting'
      write(3,*)
      write(3,*) 'import model loop_'//model

com-- this creates a wrapper function of the wrapper functions
      open(unit=2,file='loop_matrix_lib.f',status='unknown')
      write(2,*) '     subroutine loop_matrix_lib(procid,P,WGT)'
      write(2,*) '     implicit none'
      write(2,*) '     character*200 procid'
      write(2,*) '     integer init1, init2'
      write(2,*) '     double precision P(0:3,*),WGT'
      write(2,*) '     '
      write(2,*) '     read(procid,*) init1 , init2'
com-- this creates now separate files for each subchannel
com-- each of these files contain a different library
com-- with the name of file and subroutine containing the subchannel 
      do i=1,ident_max
       j=i+10
        open(unit=j,file='loop_matrix_lib_'
     f        //trim(adjustl(initial_states(i)))
     f        //'.f',status='unknown')
        write(j,*) '     subroutine loop_matrix_lib_'
     f      //trim(adjustl(initial_states(i)))//'(procid,P,WGT)'
        write(j,*) '     implicit none'
        write(j,*) '     character*200 procid'
        write(j,*) '     integer i,j,ans_dim'
        write(j,*) '     double precision P3(0:3,3),P4(0:3,4),P5(0:3,5)'
        write(j,*) '    f     ,P6(0:3,6),P7(0:3,7),P8(0:3,8),P9(0:3,9)'
        write(j,*) '    f     ,P(0:3,*),ANS(0:3,0:1),WGT'
        write(j,*) '     '
      enddo

      do j=1,total_proc

         if (npart(j).gt.9) then
            write (*,*) 'npart(j) is too large',npart(j)
            stop 1
         endif

         do k=1,npart(j)
            flst_ml5(k)=partid(k,j)
            if (flst_ml5(k).eq.0) flst_ml5(k)=21
         enddo
         write(3,*) 'add process '
     $        ,(flst_ml5(k),k=1,2),' >', (flst_ml5(k),k=3,npart(j))
     $        ,' QED=',alpha(j),' QCD=',alphas(j)-1,' [virt=QCD] @ ',j

com-- create array with incremental numbers to be permuted
         do i=1,npart(j)-2
            permut(i)=i
         enddo
com-- 1110 goto-loop loops over all permutations of array permut
 1110    procid=""
com-- do-loop loops over initial state particles to apply permutations
         do k=1,2
            write(str_temp,*) flst_ml5(k)
            procid = trim(adjustl(procid))//' '//trim(adjustl(str_temp))
         enddo
com-- do-loop loops over all final state particles to apply permutations
         do i=1,npart(j)-2
            write(str_temp,*) flst_ml5(2+permut(i))
            procid = trim(adjustl(procid))//' '//trim(adjustl(str_temp))
         enddo

com-- find out to which file to write the process (ie, which subchannel)
         write (incoming1,*) partid(1,j)
         write (incoming2,*) partid(2,j)
         do i=1,200
            if (incoming1(i:i) .eq. "-") incoming1(i:i) = "m"
            if (incoming2(i:i) .eq. "-") incoming2(i:i) = "m"
         enddo
        incoming=trim(adjustl(incoming1))//"_"//trim(adjustl(incoming2))
         found = .false.
         do i = 1,200
            if (initial_states(i) .eq. incoming) then
               ident = i
               found = .true.
               exit
            endif
         end do
         if (.not. found) then
            write(*,*) "ERROR could not find channel in initial_states"
            stop
         endif

         write(str_j,*) j
         write(str_mult,*) npart(j)
com-- determine first write for each ident value
         if (firsttime(ident)) then
            write(ident+10,*) '     if(trim(procid) .eq. "',trim(procid)
     &           ,'") then'
            firsttime(ident)=.false.
com-- use this loop also to create the wrapper of the wrappers
            write(str_init1,*) partid(1,j)
            write(str_init2,*) partid(2,j)            
            if(j.eq.1) then
               write (2,*)
     &       '     if(init1 .eq. ',trim(adjustl(str_init1))
     &      ,' .and. init2 .eq. ',trim(adjustl(str_init2)),') then'
               write(2,*) '        call loop_matrix_lib_'
     f              //trim(adjustl(initial_states(i)))//'(procid,P,WGT)'
            else
               write (2,*)
     &       '     elseif(init1 .eq. ',trim(adjustl(str_init1))
     &      ,' .and. init2 .eq. ',trim(adjustl(str_init2)),') then'
               write(2,*) '        call loop_matrix_lib_'
     f              //trim(adjustl(initial_states(i)))//'(procid,P,WGT)'
            endif
         else
            write(ident+10,*)'     elseif(trim(procid) .eq. "'
     &           ,trim(procid),'") then'
         endif
         write(ident+10,*) '         do j=0,3'
         write(ident+10,*) '            P',trim(adjustl(str_mult))
     &        ,'(j,1)=P(j,1)'
         write(ident+10,*) '            P',trim(adjustl(str_mult))
     &        ,'(j,2)=P(j,2)'
         do i=1,npart(j)-2
         write(str_permut,*) 2+permut(i)
         write(str_i,*) 2+i
         write(ident+10,*) '            P',trim(adjustl(str_mult))
     &        ,'(j,',trim(adjustl(str_permut)),')=P(j,'
     &        ,trim(adjustl(str_i)),')'
         enddo
         write(ident+10,*) '         enddo'
         write(ident+10,*) '         call ML5_',trim(adjustl(str_j))
     &        ,'_GET_ANSWER_DIMENSION(ans_dim)'
         write(ident+10,*) '         if (ans_dim.ne.1) then'
       write(ident+10,*)
     &        '            write (*,*) "ERROR ans_dim not zero",'
     &        //'ans_dim'
         write(ident+10,*) '            stop 1'
         write(ident+10,*) '         endif'
         write(ident+10,*) '         call ML5_',trim(adjustl(str_j))
     f        ,'_SLOOPMATRIX(P',trim(adjustl(str_mult)),',ANS)'
         if(nextp(npart(j)-2,permut) .and. npart(j) .gt. 3) go to 1110

         if (j.eq. total_proc) then
           do i=1,ident_max
             write(i+10,*) '     else'
            write(i+10,*)'        write (*,*) "ERROR procid not found",'
     &            //' procid'
             write(i+10,*) '        stop 1'
             write(i+10,*) '     endif'
           enddo
           write(2,*) '     else'
           write(2,*)
     &'        write (*,*)"ERROR procid does not match initial states",'
     &            //' procid'
           write(2,*) '        stop 1'
           write(2,*) '     endif'
         endif
         
      enddo

      write(3,*) 'output ',trim(adjustl(library))
      write(3,*) 'launch -f'
      close(3)
      write(2,*) '     end'
      write(2,*) '     '
      close(2)
      do i=1,ident_max
         write(i+10,*) '     '
         write(i+10,*) '     WGT=ANS(1,0)'
         write(i+10,*) '     '
         write(i+10,*) '     return'
         write(i+10,*) '     end'
         write(i+10,*) '     '
         close(i+10)
      enddo

      call system("../../bin/mg5_aMC MadLoop.mg5")
c$$$      call system("mv loop_matrix_lib.f "
c$$$     f     //trim(adjustl(library))//"/SubProcesses/")
c$$$      call system("mv "//trim(adjustl(library))//
c$$$     f     "/SubProcesses/makefile "//trim(adjustl(library))//
c$$$     f     "/SubProcesses/makefile_orig")
c$$$c$$$      call system("cp "//trim(adjustl(library))//
c$$$c$$$     f "/SubProcesses/P1_gg_h/nsquaredSO.inc "//trim(adjustl(library))//
c$$$c$$$     f  "/SubProcesses/")
c$$$
c$$$      open (unit=20,file=trim(adjustl(library))//
c$$$     f     "/SubProcesses/makefile_orig",status='old',err=43)
c$$$      open(unit=21,file=trim(adjustl(library))//
c$$$     f     "/SubProcesses/makefile",err=44)
c$$$      i = 0
c$$$      iproc = 0
c$$$ 200  i = i + 1
c$$$      read(20,10,END=201) line
c$$$      if (line(1:11).eq.'OLP_PROCESS') then
c$$$         len_line=LEN(TRIM(line))
c$$$         line=line(1:len_line-2)//" loop_matrix_lib.o \ "
c$$$      endif
c$$$      write(21,10) trim(line)
c$$$      GOTO 200
c$$$ 201  close(20)
c$$$      close(21)

      call system("cd "//trim(adjustl(library))//"/SubProcesses/;"
     & //"make OLP_static > create_lib.sh;"
     & //"chmod +x create_lib.sh;"
     & //"./create_lib.sh;"
     & //"cp libMadLoop.* ../lib/;"
     & //"cd ../../")
      write(*,*) ""
      write(*,*) "Ignore previous ERROR! library should be created now."
      write(*,*) ""


      return
 32   write (*,*)
     &     "ERROR: file 'proc_number' not found or not correct format."
      return
 33   write (*,*)
     &    "ERROR: file 'processes.dat' not found or not correct format."
 34   write (*,*)
     &  "ERROR: file 'proc_couplings' not found or not correct format."
      return
 35   write(*,*) '******************************************'
      write(*,*) 'Problem in writing the order file'
      write(*,*) 'CANNOT PROCEED. POWHEG execution abort'
      write(*,*) '******************************************'
      return
 36   write (*,*) "ERROR: file 'proc_couplings'"/
     $     /" not found or not correct format."
      return
 37   write (*,*) "ERROR: file '../Cards/proc_card.dat'"/
     $     /" not found or not correct format."
      return
 43   write (*,*)
     &    "ERROR: file '"//trim(adjustl(library))//"/SubProcesses/"//
     f     "makefile_orig' not found or not correct format."
 44   write (*,*)
     &    "ERROR: file '"//trim(adjustl(library))//"/SubProcesses/"//
     f     "makefile' not found or not correct format."
 10   format(A)
      end

      function nextp(n,a)
      integer n,a,i,j,k,t
      logical nextp
      dimension a(n)
      i=n-1
   10 if(a(i).lt.a(i+1)) go to 20
      i=i-1
      if(i.eq.0) go to 20
      go to 10
   20 j=i+1
      k=n
   30 t=a(j)
      a(j)=a(k)
      a(k)=t
      j=j+1
      k=k-1
      if(j.lt.k) go to 30
      j=i
      if(j.ne.0) go to 40
      nextp=.false.
      return
   40 j=j+1
      if(a(j).lt.a(i)) go to 40
      t=a(i)
      a(i)=a(j)
      a(j)=t
      nextp=.true.
      end
