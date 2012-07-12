C
C IPNEXT gives the next permutation given a permutation with N elements
C
      SUBROUTINE IPNEXT(IA, N, FLAG)

      IMPLICIT NONE

      INTEGER IA(*)
      INTEGER N
      INTEGER FLAG
      
      INTEGER I,J
      INTEGER ITEMP


      I = N-1

      DO WHILE (I.GE.1 .AND. IA(I).GT.IA(I+1))
         I = I-1
      ENDDO

      IF (I<1) THEN
C      IF (I.LT.1) THEN
         FLAG = -1
         RETURN
      ENDIF

      J = N

      DO WHILE (IA(I).GT.IA(J))
         J = J-1
      ENDDO

      ITEMP = IA(I)
      IA(I) = IA(J)
      IA(J) = ITEMP

      I = I+1
      J = N

      DO WHILE (I.LT.J)
         ITEMP = IA(I)
         IA(I) = IA(J)
         IA(J) = ITEMP
         I = I+1
         J = J-1
      ENDDO

      FLAG = 1

      RETURN
      END


C
C     INSORT inserts INEW in the sorted array IA with N entries
C
      SUBROUTINE INSORT(INEW,N,IA)
      IMPLICIT NONE
      INTEGER INEW,N,IA(*),I,ITMP
      
      DO I=N,2,-1
         IF(INEW.GE.IA(I))THEN
            IA(I+1)=INEW
            EXIT
         ELSE
            IA(I+1)=IA(I)
         ENDIF
      ENDDO
      END

