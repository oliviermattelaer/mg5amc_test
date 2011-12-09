C
C ipnext gives the next permutation given a permutation with n elements
C
      subroutine ipnext(ia, n, flag)

      implicit none

      integer ia(*)
      integer n
      integer flag
      
      integer i,j
      integer itemp


      i = n-1

      do while (i.ge.1 .and. ia(i).gt.ia(i+1))
         i = i-1
      enddo

      if (i<1) then
c      if (i.lt.1) then
         flag = -1
         return
      endif

      j = n

      do while (ia(i).gt.ia(j))
         j = j-1
      enddo

      itemp = ia(i)
      ia(i) = ia(j)
      ia(j) = itemp

      i = i+1
      j = n

      do while (i.lt.j)
         itemp = ia(i)
         ia(i) = ia(j)
         ia(j) = itemp
         i = i+1
         j = j-1
      enddo

      flag = 1

      return
      end
