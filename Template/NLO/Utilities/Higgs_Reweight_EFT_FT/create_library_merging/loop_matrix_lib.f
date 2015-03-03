      subroutine loop_matrix_lib(procid,P,WGT)
      implicit none
      INCLUDE "nexternal.inc"
      character*200 procid
      integer i,j,ans_dim
      double precision P3(0:3,3),P4(0:3,4),P5(0:3,5)
     f     ,P6(0:3,6),P7(0:3,7),P8(0:3,8),P9(0:3,9)
     f     ,P(0:3,nexternal),ANS(0:3,0:0),WGT
      
      if(trim(procid) .eq. "21 21 25") then
          do j=0,3
             P3(j,1)=P(j,1)
             P3(j,2)=P(j,2)
             P3(j,3)=P(j,3)
          enddo
          call ML5_1_GET_ANSWER_DIMENSION(ans_dim)
          if (ans_dim.ne.0) then
             write (*,*) "ERROR ans_dim not zero",ans_dim
             stop 1
          endif
          call ML5_1_SLOOPMATRIX(P3,ANS)
      elseif(trim(procid) .eq. "21 21 21 25") then
          do j=0,3
             P4(j,1)=P(j,1)
             P4(j,2)=P(j,2)
             P4(j,3)=P(j,3)
             P4(j,4)=P(j,4)
          enddo
          call ML5_2_GET_ANSWER_DIMENSION(ans_dim)
          if (ans_dim.ne.0) then
             write (*,*) "ERROR ans_dim not zero",ans_dim
             stop 1
          endif
          call ML5_2_SLOOPMATRIX(P4,ANS)
      elseif(trim(procid) .eq. "21 21 25 21") then
          do j=0,3
             P4(j,1)=P(j,1)
             P4(j,2)=P(j,2)
             P4(j,3)=P(j,4)
             P4(j,4)=P(j,3)
          enddo
          call ML5_2_GET_ANSWER_DIMENSION(ans_dim)
          if (ans_dim.ne.0) then
             write (*,*) "ERROR ans_dim not zero",ans_dim
             stop 1
          endif
          call ML5_2_SLOOPMATRIX(P4,ANS)
      elseif(trim(procid) .eq. "21 1 25 1") then
          do j=0,3
             P4(j,1)=P(j,1)
             P4(j,2)=P(j,2)
             P4(j,3)=P(j,3)
             P4(j,4)=P(j,4)
          enddo
          call ML5_3_GET_ANSWER_DIMENSION(ans_dim)
          if (ans_dim.ne.0) then
             write (*,*) "ERROR ans_dim not zero",ans_dim
             stop 1
          endif
          call ML5_3_SLOOPMATRIX(P4,ANS)
      elseif(trim(procid) .eq. "21 1 1 25") then
          do j=0,3
             P4(j,1)=P(j,1)
             P4(j,2)=P(j,2)
             P4(j,3)=P(j,4)
             P4(j,4)=P(j,3)
          enddo
          call ML5_3_GET_ANSWER_DIMENSION(ans_dim)
          if (ans_dim.ne.0) then
             write (*,*) "ERROR ans_dim not zero",ans_dim
             stop 1
          endif
          call ML5_3_SLOOPMATRIX(P4,ANS)
      elseif(trim(procid) .eq. "1 -1 21 25") then
          do j=0,3
             P4(j,1)=P(j,1)
             P4(j,2)=P(j,2)
             P4(j,3)=P(j,3)
             P4(j,4)=P(j,4)
          enddo
          call ML5_4_GET_ANSWER_DIMENSION(ans_dim)
          if (ans_dim.ne.0) then
             write (*,*) "ERROR ans_dim not zero",ans_dim
             stop 1
          endif
          call ML5_4_SLOOPMATRIX(P4,ANS)
      elseif(trim(procid) .eq. "1 -1 25 21") then
          do j=0,3
             P4(j,1)=P(j,1)
             P4(j,2)=P(j,2)
             P4(j,3)=P(j,4)
             P4(j,4)=P(j,3)
          enddo
          call ML5_4_GET_ANSWER_DIMENSION(ans_dim)
          if (ans_dim.ne.0) then
             write (*,*) "ERROR ans_dim not zero",ans_dim
             stop 1
          endif
          call ML5_4_SLOOPMATRIX(P4,ANS)
      elseif(trim(procid) .eq. "21 -5 25 21 -5") then
          do j=0,3
             P5(j,1)=P(j,1)
             P5(j,2)=P(j,2)
             P5(j,3)=P(j,3)
             P5(j,4)=P(j,4)
             P5(j,5)=P(j,5)
          enddo
          call ML5_5_GET_ANSWER_DIMENSION(ans_dim)
          if (ans_dim.ne.0) then
             write (*,*) "ERROR ans_dim not zero",ans_dim
             stop 1
          endif
          call ML5_5_SLOOPMATRIX(P5,ANS)
      elseif(trim(procid) .eq. "21 -5 25 -5 21") then
          do j=0,3
             P5(j,1)=P(j,1)
             P5(j,2)=P(j,2)
             P5(j,3)=P(j,3)
             P5(j,4)=P(j,5)
             P5(j,5)=P(j,4)
          enddo
          call ML5_5_GET_ANSWER_DIMENSION(ans_dim)
          if (ans_dim.ne.0) then
             write (*,*) "ERROR ans_dim not zero",ans_dim
             stop 1
          endif
          call ML5_5_SLOOPMATRIX(P5,ANS)
      elseif(trim(procid) .eq. "21 -5 21 25 -5") then
          do j=0,3
             P5(j,1)=P(j,1)
             P5(j,2)=P(j,2)
             P5(j,3)=P(j,4)
             P5(j,4)=P(j,3)
             P5(j,5)=P(j,5)
          enddo
          call ML5_5_GET_ANSWER_DIMENSION(ans_dim)
          if (ans_dim.ne.0) then
             write (*,*) "ERROR ans_dim not zero",ans_dim
             stop 1
          endif
          call ML5_5_SLOOPMATRIX(P5,ANS)
      elseif(trim(procid) .eq. "21 -5 21 -5 25") then
          do j=0,3
             P5(j,1)=P(j,1)
             P5(j,2)=P(j,2)
             P5(j,3)=P(j,4)
             P5(j,4)=P(j,5)
             P5(j,5)=P(j,3)
          enddo
          call ML5_5_GET_ANSWER_DIMENSION(ans_dim)
          if (ans_dim.ne.0) then
             write (*,*) "ERROR ans_dim not zero",ans_dim
             stop 1
          endif
          call ML5_5_SLOOPMATRIX(P5,ANS)
      elseif(trim(procid) .eq. "21 -5 -5 25 21") then
          do j=0,3
             P5(j,1)=P(j,1)
             P5(j,2)=P(j,2)
             P5(j,3)=P(j,5)
             P5(j,4)=P(j,3)
             P5(j,5)=P(j,4)
          enddo
          call ML5_5_GET_ANSWER_DIMENSION(ans_dim)
          if (ans_dim.ne.0) then
             write (*,*) "ERROR ans_dim not zero",ans_dim
             stop 1
          endif
          call ML5_5_SLOOPMATRIX(P5,ANS)
      elseif(trim(procid) .eq. "21 -5 -5 21 25") then
          do j=0,3
             P5(j,1)=P(j,1)
             P5(j,2)=P(j,2)
             P5(j,3)=P(j,5)
             P5(j,4)=P(j,4)
             P5(j,5)=P(j,3)
          enddo
          call ML5_5_GET_ANSWER_DIMENSION(ans_dim)
          if (ans_dim.ne.0) then
             write (*,*) "ERROR ans_dim not zero",ans_dim
             stop 1
          endif
          call ML5_5_SLOOPMATRIX(P5,ANS)
      elseif(trim(procid) .eq. "2 21 25 21 2") then
          do j=0,3
             P5(j,1)=P(j,1)
             P5(j,2)=P(j,2)
             P5(j,3)=P(j,3)
             P5(j,4)=P(j,4)
             P5(j,5)=P(j,5)
          enddo
          call ML5_6_GET_ANSWER_DIMENSION(ans_dim)
          if (ans_dim.ne.0) then
             write (*,*) "ERROR ans_dim not zero",ans_dim
             stop 1
          endif
          call ML5_6_SLOOPMATRIX(P5,ANS)
      elseif(trim(procid) .eq. "2 21 25 2 21") then
          do j=0,3
             P5(j,1)=P(j,1)
             P5(j,2)=P(j,2)
             P5(j,3)=P(j,3)
             P5(j,4)=P(j,5)
             P5(j,5)=P(j,4)
          enddo
          call ML5_6_GET_ANSWER_DIMENSION(ans_dim)
          if (ans_dim.ne.0) then
             write (*,*) "ERROR ans_dim not zero",ans_dim
             stop 1
          endif
          call ML5_6_SLOOPMATRIX(P5,ANS)
      elseif(trim(procid) .eq. "2 21 21 25 2") then
          do j=0,3
             P5(j,1)=P(j,1)
             P5(j,2)=P(j,2)
             P5(j,3)=P(j,4)
             P5(j,4)=P(j,3)
             P5(j,5)=P(j,5)
          enddo
          call ML5_6_GET_ANSWER_DIMENSION(ans_dim)
          if (ans_dim.ne.0) then
             write (*,*) "ERROR ans_dim not zero",ans_dim
             stop 1
          endif
          call ML5_6_SLOOPMATRIX(P5,ANS)
      elseif(trim(procid) .eq. "2 21 21 2 25") then
          do j=0,3
             P5(j,1)=P(j,1)
             P5(j,2)=P(j,2)
             P5(j,3)=P(j,4)
             P5(j,4)=P(j,5)
             P5(j,5)=P(j,3)
          enddo
          call ML5_6_GET_ANSWER_DIMENSION(ans_dim)
          if (ans_dim.ne.0) then
             write (*,*) "ERROR ans_dim not zero",ans_dim
             stop 1
          endif
          call ML5_6_SLOOPMATRIX(P5,ANS)
      elseif(trim(procid) .eq. "2 21 2 25 21") then
          do j=0,3
             P5(j,1)=P(j,1)
             P5(j,2)=P(j,2)
             P5(j,3)=P(j,5)
             P5(j,4)=P(j,3)
             P5(j,5)=P(j,4)
          enddo
          call ML5_6_GET_ANSWER_DIMENSION(ans_dim)
          if (ans_dim.ne.0) then
             write (*,*) "ERROR ans_dim not zero",ans_dim
             stop 1
          endif
          call ML5_6_SLOOPMATRIX(P5,ANS)
      elseif(trim(procid) .eq. "2 21 2 21 25") then
          do j=0,3
             P5(j,1)=P(j,1)
             P5(j,2)=P(j,2)
             P5(j,3)=P(j,5)
             P5(j,4)=P(j,4)
             P5(j,5)=P(j,3)
          enddo
          call ML5_6_GET_ANSWER_DIMENSION(ans_dim)
          if (ans_dim.ne.0) then
             write (*,*) "ERROR ans_dim not zero",ans_dim
             stop 1
          endif
          call ML5_6_SLOOPMATRIX(P5,ANS)
      else
         write (*,*) "ERROR procid not found", procid
         stop 1
      endif
      
      WGT=ANS(1,0)
      
      return
      end
      
