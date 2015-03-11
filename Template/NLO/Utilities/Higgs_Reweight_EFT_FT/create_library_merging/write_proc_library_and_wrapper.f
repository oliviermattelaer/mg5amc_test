      program write_proc_labels_virtual
      implicit none
      character*200 procid,str_j,str_mult,str_temp,library,line
     &     ,str_permut,str_i
      character*1 s1
      character*2 s2
      character*3 s3
      character*60 buff,model
      integer maxproc,maxpart
      parameter (maxproc=999,maxpart=20)
      integer total_proc,i,j,proc_number(maxproc),npart(maxproc),
     &     partid(maxpart,maxproc),lstr(maxproc),l,iproc,mult
      character*140 str(maxproc)
      character*2 ls
      logical need_switching
      integer date(3),time(3),len_line
      integer flst_ml5(maxproc),k,alpha(maxproc),alphas(maxproc)
      integer n,permut(10)
      logical nextp
      external nextp

      model='sm-no_b_mass'
      library='ML5lib_reweight'

      open (unit=11,file='processes.dat',status='old',err=33)
      i = 0
      iproc = 0
 100  i = i + 1
      read(11,*,END=101) line
      if ((line(1:1).eq.' ').or.(line(1:1).eq.'#')) then
         GOTO 100
      endif
      iproc = iproc+1
      backspace(11)
      read(11,*,err=33) proc_number(iproc),npart(iproc),
     &     (partid(j,iproc),j=1,npart(iproc)),alpha(iproc),alphas(iproc)
      total_proc=proc_number(iproc)
      GOTO 100
 101  close(11)

      open(unit=18,file='MadLoop.mg5',err=35)
      open(unit=12,file='loop_matrix_lib.f',status='unknown')
      write(*,*) '************************************************'//
     $     '******'
      write(*,*) 'Writing the orderfile for MadLoop5...'
      write(*,*) 'and a wrapper function to access the library...'
      write(*,*) ''
       

      write(18,*) '# MadLoop.mg5'
      write(18,*) '# Created for heft-mergin + reweighting'
      write(18,*)
      write(18,*) 'import model loop_'//model
      write (12,*) 
     f '     subroutine loop_matrix_lib(procid,P,WGT)'
      write (12,*) '     implicit none'
      write (12,*) '     character*200 procid'
      write (12,*) '     integer i,j,ans_dim'
      write (12,*) '     double precision P3(0:3,3),P4(0:3,4),P5(0:3,5)'
      write (12,*) '    f     ,P6(0:3,6),P7(0:3,7),P8(0:3,8),P9(0:3,9)'
      write (12,*) '    f     ,P(0:3,*),ANS(0:3,0:1),WGT'
      write (12,*) '     '
      do j=1,total_proc

         if (npart(j).gt.9) then
            write (*,*) 'npart(j) is too large',npart(j)
            stop 1
         endif

         do k=1,npart(j)
            flst_ml5(k)=partid(k,j)
            if (flst_ml5(k).eq.0) flst_ml5(k)=21
         enddo
         write(18,*) 'add process '
     $        ,(flst_ml5(k),k=1,2),' >', (flst_ml5(k),k=3,npart(j))
     $        ,' QED=',alpha(j),' QCD=',alphas(j)-1,' [virt=QCD] @ ',j


com-- create array with incremental numbers to be permuted
         do i=1,npart(j)-2
            permut(i)=i
         enddo
com-- 1110 goto-loop loops over all permutations of array permute
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

         write(str_j,*) j
         write(str_mult,*) npart(j)
         if (j.eq.1) then
            write (12,*) '     if(trim(procid) .eq. "',trim(procid)
     &           ,'") then'
         else
            write (12,*) '     elseif(trim(procid) .eq. "',trim(procid)
     &           ,'") then'
         endif
         write(12,*) '         do j=0,3'
         write(12,*) '            P',trim(adjustl(str_mult))
     &        ,'(j,1)=P(j,1)'
         write(12,*) '            P',trim(adjustl(str_mult))
     &        ,'(j,2)=P(j,2)'
         do i=1,npart(j)-2
         write(str_permut,*) 2+permut(i)
         write(str_i,*) 2+i
         write(12,*) '            P',trim(adjustl(str_mult))
     &        ,'(j,',trim(adjustl(str_i)),')=P(j,'
     &        ,trim(adjustl(str_permut)),')'
         enddo
         write(12,*) '         enddo'
         write(12,*) '         call ML5_',trim(adjustl(str_j))
     &        ,'_GET_ANSWER_DIMENSION(ans_dim)'
         write(12,*) '         if (ans_dim.ne.1) then'
         write(12,*) '            write (*,*) "ERROR ans_dim not zero",'
     &        //'ans_dim'
         write(12,*) '            stop 1'
         write(12,*) '         endif'
         write(12,*) '         call ML5_',trim(adjustl(str_j))
     f        ,'_SLOOPMATRIX(P',trim(adjustl(str_mult)),',ANS)'
         if(nextp(npart(j)-2,permut) .and. npart(j) .gt. 3) go to 1110

         if (j.eq. total_proc) then
            write(12,*) '     else'
            write(12,*) '        write (*,*) "ERROR procid not found",'
     &           //' procid'
            write(12,*) '        stop 1'
            write(12,*) '     endif'
         endif
         
      enddo

      write(18,*) 'output ',trim(adjustl(library))
      write(18,*) 'launch -f'
      write(12,*) '     '
      write(12,*) '     WGT=ANS(1,0)'
      write(12,*) '     '
      write(12,*) '     return'
      write(12,*) '     end'
      write(12,*) '     '
      close(18)
      close(12)

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

      call system("cd "//trim(adjustl(library))//"/SubProcesses/ "/
     &     /"; make OLP_static ; cd ../../" )
      

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
