      SUBROUTINE GENERATE_CF_INTEGRATION_BASIS()
      IMPLICIT NONE
      INTEGER NDIAGS
      PARAMETER (NDIAGS=%(ndiags)i)
      INTEGER    NPERMS
      PARAMETER (NPERMS=%(nperms)i)

      LOGICAL CF_BASIS(NDIAGS,NPERMS)
      COMMON/TO_CF_BASIS/CF_BASIS

C     CONSTANTS
C     GLOBAL VARIABLES
      INCLUDE 'maxconfigs.inc'
      INCLUDE 'maxparticles.inc'
C     ARGUMENTS
C     LOCAL VARIABLES 
      INTEGER I,J,K
      INTEGER MAX_NATM,NATM,MAX_CH,MAX_CHCF,IERR
      PARAMETER (MAX_NATM=MAX_PARTICLES -2)
      PARAMETER (MAX_CH=10000,MAX_CHCF=%(nperms)i)
      INTEGER ATM(2,MAX_NATM),CF(MAX_PARTICLES,MAX_CH)
      INTEGER CFNUM(MAX_CHCF,MAX_CH)
      INTEGER NDIAG,NGLUONS,NCF,NCHCF(MAX_CH)
      INTEGER CH_CF(MAX_PARTICLES,MAX_CHCF,MAX_CH)
      LOGICAL KEEP(MAX_CH,10000)
C     EXTERNAL FUNCTIONS
      INTEGER FACT
      EXTERNAL FACT
C     ----------
C     BEGIN CODE
C     ----------
      NGLUONS = MAX_PARTICLES
      DO I = 1,NDIAGS
        DO J = 1,NPERMS
          CF_BASIS(I,J) = .TRUE.
        ENDDO
      ENDDO
      RETURN
      DO I = 1,NDIAGS
        CALL MAKE_ATM(I,MAX_NATM,ATM,NATM)
        CALL FIND_PERCHCF(NGLUONS,NATM,ATM,CH_CF(1,1,I),CFNUM(1,I),NCHCF(I))

        DO J = 1,NCHCF(I)
          DO K = 1,NPERMS
            IF (CFNUM(J,I).EQ.K) THEN
              CF_BASIS(I,K) = .TRUE.
              EXIT
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      END


      SUBROUTINE MAKE_ATM(ICH,MAX_NATM,ATM_OUT,NATM)
C     ****************************************************
C     By Yoshitaro Takaesu @U.tokyo Jun.25 2014
C     ****************************************************
      IMPLICITNONE
C     CONSTANTS     
      INCLUDE 'maxconfigs.inc'
      INCLUDE 'genps.inc'
      INCLUDE 'maxamps.inc'
C     ARGUMENTS 
      INTEGER MAX_NATM,ICH,NATM
      INTEGER ATM_OUT(2,MAX_NATM)
C     GLOBAL VARIABLES
      INTEGER IFOREST(2,-MAX_BRANCH:-1,LMAXCONFIGS)
      COMMON/TO_FOREST/ IFOREST
      INTEGER SPROP(MAXSPROC,-MAX_BRANCH:-1,LMAXCONFIGS)
      INTEGER TPRID(-MAX_BRANCH:-1,LMAXCONFIGS)
      COMMON/TO_SPROP/SPROP,TPRID
      INTEGER            MAPCONFIG(0:LMAXCONFIGS), THIS_CONFIG
      COMMON/TO_MCONFIGS/MAPCONFIG, THIS_CONFIG
C     LOCAL VARIABLES
      INTEGER I
c      include 'configs.inc'
      integer MAXNPROP,NEXT,NFIN
      parameter (NEXT=MAX_PARTICLES,NFIN=MAX_PARTICLES-2)
      parameter (MAXNPROP=-1*NFIN+1)
      INTEGER ATM(2,20),IPROP(MAX_NATM),NT,NS
      INTEGER IUSED(MAX_PARTICLES),IPOST,IPOSS,IDIAG,IPOS
      INTEGER IPOS_TMP,NN,IPOS_TMP2,ATM1(10),IATM1(10)
      INTEGER PROP2_1,PROP2_2,iused_prop(-1*nfin+1:-1)
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
C     Initialization
      IDIAG = ICH
      IPOS_TMP = 9  ! temporary atm number for s-splitting atm attached to particle 2
      do i = 1,max_natm
         atm(1,i) = 0
         atm(2,i) = 0
      enddo
      do i = MAXNPROP,-1,1
         iused_prop(i) = 0
      enddo
      ATM(1,1) = 1
      ATM(1,2) = 2
      IPOS = 2
C     inerpriting iforest information into atm information
      PROP2_1 = IFOREST(1,MAXNPROP-1,IDIAG)
      PROP2_2 = IFOREST(2,MAXNPROP-1,IDIAG)
C     store a splitting leg attached to the particle 2
      IF ((PROP2_1.LT.0).AND.(SPROP(1,MIN(PROP2_1,-1),IDIAG).NE.0)) THEN
        IPOS = IPOS +1
        ATM(1,IPOS) = IFOREST(1,PROP2_1,IDIAG)
        ATM(2,IPOS) = IFOREST(2,PROP2_1,IDIAG)
        ATM(2,2) = 0
        iused_prop(prop2_1) = 1
      ELSEIF ((PROP2_2.LT.0).AND.(SPROP(1,MIN(PROP2_2,-1),IDIAG).NE.0)) THEN
        IPOS = IPOS +1
        ATM(1,IPOS) = IFOREST(1,PROP2_2,IDIAG)
        ATM(2,IPOS) = IFOREST(2,PROP2_2,IDIAG)
        ATM(2,2) = 0
        iused_prop(prop2_2) = 1
C       store an external attached to the particle 2
      ELSE
        IF (IFOREST(2,MAXNPROP-1,IDIAG).GT.0) THEN
          ATM(2,2) = IFOREST(2,MAXNPROP-1,IDIAG)
        ELSE
         write(*,*) "ERROR: atm for particle 2 is not set appropretely."
           return
       ENDIF
      ENDIF

      IPOS_TMP2 = 10
      DO I = MAXNPROP,-1,1
        IF (iused_prop(i).ne.1) THEN
C         check atm1
          IF (IFOREST(1,I,IDIAG).EQ.1) THEN
C           store an external leg attarched to the particle 1
            IF (IFOREST(2,I,IDIAG).GT.0) THEN
              ATM(2,1) = IFOREST(2,I,IDIAG)
C             store a splitting leg attached to the particle 1
            ELSE
              IPOS_TMP = IPOS_TMP +1
              ATM(2,1) = 0
              IPROP(1) = IFOREST(2,I,IDIAG)
              ATM(1,IPOS_TMP) = IFOREST(1,IPROP(1),IDIAG)
              ATM(2,IPOS_TMP) = IFOREST(2,IPROP(1),IDIAG)
              iused_prop(iprop(1)) = 1
            ENDIF
C           check the rest of external legs
          ELSEIF (IFOREST(2,I,IDIAG).GT.0) THEN
C           store an external leg attached to a propagator
            IF (IFOREST(1,I,IDIAG).LT.0) THEN
              IPOS_TMP2 = IPOS_TMP2 +1
              ATM(1,IPOS_TMP2) = IFOREST(2,I,IDIAG)
              ATM(2,IPOS_TMP2) = 0
C             store a splitting leg attached to a propagator
            ELSEIF (IFOREST(1,I,IDIAG).GT.0) THEN
              IPOS_TMP2 = IPOS_TMP2 +1
              ATM(1,IPOS_TMP2) = IFOREST(1,I,IDIAG)
              ATM(2,IPOS_TMP2) = IFOREST(2,I,IDIAG)
            ENDIF
          ENDIF
        ENDIF
      ENDDO

C     sort the atm's ascendent order w.r.t the first element, atm(1,*)
      NN = IPOS_TMP2 -10
      DO I = 1,NN
        ATM1(I) = ATM(1,I+10)
      ENDDO
      CALL ISORT_ASC(NN,ATM1,IATM1)
      DO I = 1,NN
        IPOS = IPOS +1
        ATM(1,IPOS) = ATM(1,IATM1(I)+10)
        ATM(2,IPOS) = ATM(2,IATM1(I)+10)
      ENDDO
      IF (IPOS_TMP.EQ.10) THEN
        IPOS = IPOS +1
        ATM(1,IPOS) = ATM(1,IPOS_TMP)
        ATM(2,IPOS) = ATM(2,IPOS_TMP)
      ENDIF

      natm = max_natm
      DO I = 1,MAX_NATM
         if (atm(1,i).eq.0) then
            natm = natm -1
            ATM_OUT(1,I) = 0
            ATM_OUT(2,I) = 0
         ELSE
            ATM_OUT(1,I) = ATM(1,I)
            ATM_OUT(2,I) = ATM(2,I)
         ENDIF
      ENDDO

      RETURN
      END


      SUBROUTINE ISORT_ASC(N, A,M)
C     by Yoshitaro Takaesu -Jun/25/2014 @U.Toyo
      IMPLICIT NONE
      INTEGER I,J
      INTEGER N
      INTEGER M(N)
      INTEGER A(N)
      INTEGER MINI,TMP
      INTEGER MIN
      DO I = 1,N
        M(I) = I
      ENDDO
      DO I = 1,N-1
        MINI = I
        MIN = A(I)
        DO J = I+1,N
          IF( A(J) < MIN ) THEN
            MIN = A(J)
            MINI = J
          ENDIF
        ENDDO
        IF( MINI.NE.I ) THEN
          A(MINI) = A(I)
          A(I) = MIN
          TMP = M(MINI)
          M(MINI) = M(I)
          M(I) = TMP
        ENDIF
      ENDDO
      END

      SUBROUTINE FIND_CHCF(NEXT,NATM,ATM,CFLOW,NCF)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Nov.19 2012
C     ****************************************************
      IMPLICITNONE
C     CONSTANTS
C     ARGUMENTS 
      INTEGER NEXT,NATM
      INTEGER ATM(2,NATM),CFLOW(NEXT,100000),NCF
C     GLOBAL VARIABLES
C     LOCAL VARIABLES 
      INTEGER I,J,K
      INTEGER IFLAG,NPERM,POS2,IERR,A(NATM+1),NNCF,MAXK(NATM)
      INTEGER CFLOWTMP
C     EXTERNAL FUNCTIONS
      INTEGER FACT
      EXTERNAL FACT
C     ----------
C     BEGIN CODE
C     ----------
      IFLAG = 0
      DO I = 1,NATM
        A(I) = I
      ENDDO
      A(NATM+1) = 0
      DO I = 1,NATM
        MAXK(I) = 1
        IF (ATM(2,I).NE.0) MAXK(I) = 2
      ENDDO

      NCF = 0
      NPERM = FACT(NATM-1)
      DO I = 1,NPERM
        IF (I.NE.1) THEN
          CALL IPNEXT(A(2),NATM-1,IFLAG)
        ENDIF
        DO J = 2,NATM
          IF (A(J).EQ.2) POS2 = J
        ENDDO

        IERR = 1
        if (atm(2,1)+atm(2,2).eq.0) then
           if ((a(pos2+1).eq.3).or.(a(pos2-1).eq.3)) then
              if ((a(natm).eq.natm).or.(a(2).eq.natm)) ierr = 0
           endif
           if ((a(natm).eq.3).or.(a(2).eq.3)) then
              if ((a(pos2+1).eq.natm).or.(a(pos2-1).eq.natm)) ierr = 0
           ENDIF
        elseif (atm(2,1).eq.0) then
           if ((a(2).eq.natm).or.(a(natm).eq.natm)) ierr = 0
        elseif (atm(2,2).eq.0) then
           if ((a(pos2-1).eq.3).or.(a(pos2+1).eq.3)) ierr = 0
        ELSE
           ierr = 0
        ENDIF
        IF (IERR.NE.0) CYCLE

        CALL GET_CFLOW(NEXT,NATM,ATM,A,MAXK,NNCF,CFLOW(1,NCF+1))
        NCF = NCF +NNCF
      ENDDO

c     Put the first element of CFLOW to the end if CFLOW(1) =! 1
c     ex) 31245 -> 12453
      DO I = 1,NCF
        IF (CFLOW(1,I).NE.1) THEN
          CFLOWTMP = CFLOW(1,I)
          DO J = 1,NEXT-1
            CFLOW(J,I) = CFLOW(J+1,I)
          ENDDO
          CFLOW(NEXT,I) = CFLOWTMP
        ENDIF
      ENDDO

      RETURN
      END


      SUBROUTINE FIND_PERCHCF(NEXT,NATM,ATM,CFLOW,CFLOW_NUM,NCF)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN.14 2013
C     ****************************************************
      IMPLICITNONE
C     CONSTANTS
C     ARGUMENTS 
      INTEGER NEXT,NATM,IERR
      INTEGER ATM(2,NATM),CFLOW(NEXT,10000),NCF,CFLOW_NUM(10000)
C     GLOBAL VARIABLES
C     LOCAL VARIABLES 
      INTEGER I,J,K
      INTEGER IFLAG,NPERM,A(NATM),NONZERO,ZERO,NONZERO_ATM(2,NATM)
      INTEGER ZERO_ATM(2,NATM),INCF,IATM,AATM(2,NATM),AAATM(2,NATM)
C     EXTERNAL FUNCTIONS
      INTEGER FACT
      EXTERNAL FACT
C     ----------
C     BEGIN CODE
C     ----------
      NCF = 0
      DO I = 1,NATM
        AATM(1,I) = ATM(1,I)
        AATM(2,I) = ATM(2,I)
      ENDDO

      CALL FIND_CHCF(NEXT,NATM,AATM,CFLOW(1,NCF+1),INCF)
      NCF = NCF +INCF

C     Get color-flow numbers of each channel
      DO I = 1,NCF
        CALL FIND_CFNUMBER(NEXT,CFLOW(1,I),CFLOW_NUM(I))
        IF (CFLOW_NUM(I).EQ.-1) THEN
          WRITE(*,*) 'ERROR: Wrong color-flow is generated.'
          RETURN
        ENDIF
      ENDDO

      RETURN
      END


      SUBROUTINE FIND_CFNUMBER(NEXT,CFLOW,CFLOW_NUM)
C     ****************************************************
C     By Yoshitaro Takaesu @ U.Tokyo JUN.27 2014
C     ****************************************************
      IMPLICITNONE
C     CONSTANTS
C     ARGUMENTS 
      INTEGER NEXT,CFLOW(NEXT),CFLOW_NUM
C     GLOBAL VARIABLES
C     LOCAL VARIABLES 
      INTEGER I
      INTEGER FLOWORDER(NEXT),IFLAG1,IFLAG2,NCF
C     EXTERNAL FUNCTIONS
      INTEGER FACT
      EXTERNAL FACT
C     ----------
C     BEGIN CODE
C     ----------
C     Initialization
      NCF = FACT(NEXT-1)
      DO I = 1,NEXT
        FLOWORDER(I) = I
      ENDDO
C     Finding the color-flow number of cflow
      DO I = 1,NCF
        IF (I.NE.1) THEN
          CALL IPNEXT(FLOWORDER(2),NEXT-1,IFLAG1)
        ENDIF
        CALL CHECK_CF(NEXT,FLOWORDER,CFLOW,IFLAG2)
C       If the color-flow number is found, return the number
        IF (IFLAG2.EQ.1) THEN
          CFLOW_NUM = I
          RETURN
        ENDIF
      ENDDO
C     If the color-flow number is not found, return -1
      CFLOW_NUM = -1

      RETURN
      END



      SUBROUTINE CHECK_CF(NEXT,CF1,CF2,IFLAG)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Dec.8 2012
C     ****************************************************
      IMPLICITNONE
C     ARGUMENTS
      INTEGER NEXT
      INTEGER CF1(NEXT),CF2(NEXT),IFLAG
C     LOCAL VARIABLES 
      INTEGER I
C     ----------
C     BEGIN CODE
C     ----------
      IFLAG = 0
      DO I = 1,NEXT
        IF (CF1(I).NE.CF2(I)) RETURN
      ENDDO
      IFLAG = 1
      RETURN
      END

      INTEGER FUNCTION FACT(N)
      IMPLICIT NONE
      INTEGER N,I
      FACT = 1
      IF (N.EQ.0) THEN
        FACT = 1
      ELSE
        DO I = 1, N
          FACT = FACT*I
        ENDDO
      ENDIF
      RETURN
      END


      SUBROUTINE GET_CFLOW(NEXT,NATM,ATM,A,MAXK,NCF,CFLOW)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS AUG 25 2012
C     ****************************************************
      IMPLICIT NONE

C     GLOBAL VARIABLES

C     CONSTANTS

C     ARGUMENTS 
      INTEGER NATM,NEXT
      INTEGER A(NATM),ATM(2,NATM),MAXK(NATM),NCF
      INTEGER CFLOW(NEXT,1000)
C     LOCAL VARIABLES 
      INTEGER I,J,CFPOS,K1,K2,K3,K4,K5,K6,K7,K8
      INTEGER AATM(2,NATM)
C     EXTERNAL FUNCTIONS

C     ----------
C     BEGIN CODE
C     ----------
      NCF = 0

      IF (NATM.EQ.3) THEN
        DO K1 = 1,MAXK(1)
          DO K2 = 1,MAXK(2)
            DO K3 = 1,MAXK(3)
              DO J = 1,2
                AATM(J,1) = ATM(MOD(J*K1,3),1)
                AATM(J,2) = ATM(MOD(J*K2,3),2)
                AATM(J,3) = ATM(MOD(J*K3,3),3)
              ENDDO
              NCF = NCF +1
              CALL SORT_CFLOW(AATM,A,NATM,NEXT,CFLOW(1,NCF))
            ENDDO
          ENDDO
        ENDDO
      ELSEIF (NATM.EQ.4) THEN
        DO K1 = 1,MAXK(1)
          DO K2 = 1,MAXK(2)
            DO K3 = 1,MAXK(3)
              DO K4 = 1,MAXK(4)
                DO J = 1,2
                  AATM(J,1) = ATM(MOD(J*K1,3),1)
                  AATM(J,2) = ATM(MOD(J*K2,3),2)
                  AATM(J,3) = ATM(MOD(J*K3,3),3)
                  AATM(J,4) = ATM(MOD(J*K4,3),4)
                ENDDO
                NCF = NCF +1
                CALL SORT_CFLOW(AATM,A,NATM,NEXT,CFLOW(1,NCF))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSEIF (NATM.EQ.5) THEN
        DO K1 = 1,MAXK(1)
          DO K2 = 1,MAXK(2)
            DO K3 = 1,MAXK(3)
              DO K4 = 1,MAXK(4)
                DO K5 = 1,MAXK(5)
                  DO J = 1,2
                    AATM(J,1) = ATM(MOD(J*K1,3),1)
                    AATM(J,2) = ATM(MOD(J*K2,3),2)
                    AATM(J,3) = ATM(MOD(J*K3,3),3)
                    AATM(J,4) = ATM(MOD(J*K4,3),4)
                    AATM(J,5) = ATM(MOD(J*K5,3),5)
                  ENDDO
                  NCF = NCF +1
                  CALL SORT_CFLOW(AATM,A,NATM,NEXT,CFLOW(1,NCF))
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSEIF (NATM.EQ.6) THEN
        DO K1 = 1,MAXK(1)
          DO K2 = 1,MAXK(2)
            DO K3 = 1,MAXK(3)
              DO K4 = 1,MAXK(4)
                DO K5 = 1,MAXK(5)
                  DO K6 = 1,MAXK(6)
                    DO J = 1,2
                      AATM(J,1) = ATM(MOD(J*K1,3),1)
                      AATM(J,2) = ATM(MOD(J*K2,3),2)
                      AATM(J,3) = ATM(MOD(J*K3,3),3)
                      AATM(J,4) = ATM(MOD(J*K4,3),4)
                      AATM(J,5) = ATM(MOD(J*K5,3),5)
                      AATM(J,6) = ATM(MOD(J*K6,3),6)
                    ENDDO
                    NCF = NCF +1
                    CALL SORT_CFLOW(AATM,A,NATM,NEXT,CFLOW(1,NCF))
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSEIF (NATM.EQ.7) THEN
        DO K1 = 1,MAXK(1)
          DO K2 = 1,MAXK(2)
            DO K3 = 1,MAXK(3)
              DO K4 = 1,MAXK(4)
                DO K5 = 1,MAXK(5)
                  DO K6 = 1,MAXK(6)
                    DO K7 = 1,MAXK(7)
                      DO J = 1,2
                        AATM(J,1) = ATM(MOD(J*K1,3),1)
                        AATM(J,2) = ATM(MOD(J*K2,3),2)
                        AATM(J,3) = ATM(MOD(J*K3,3),3)
                        AATM(J,4) = ATM(MOD(J*K4,3),4)
                        AATM(J,5) = ATM(MOD(J*K5,3),5)
                        AATM(J,6) = ATM(MOD(J*K6,3),6)
                        AATM(J,7) = ATM(MOD(J*K7,3),7)
                      ENDDO
                      NCF = NCF +1
                      CALL SORT_CFLOW(AATM,A,NATM,NEXT,CFLOW(1,NCF))
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSEIF (NATM.EQ.8) THEN
        DO K1 = 1,MAXK(1)
          DO K2 = 1,MAXK(2)
            DO K3 = 1,MAXK(3)
              DO K4 = 1,MAXK(4)
                DO K5 = 1,MAXK(5)
                  DO K6 = 1,MAXK(6)
                    DO K7 = 1,MAXK(7)
                      DO K8 = 1,MAXK(8)
                        DO J = 1,2
                          AATM(J,1) = ATM(MOD(J*K1,3),1)
                          AATM(J,2) = ATM(MOD(J*K2,3),2)
                          AATM(J,3) = ATM(MOD(J*K3,3),3)
                          AATM(J,4) = ATM(MOD(J*K4,3),4)
                          AATM(J,5) = ATM(MOD(J*K5,3),5)
                          AATM(J,6) = ATM(MOD(J*K6,3),6)
                          AATM(J,7) = ATM(MOD(J*K7,3),7)
                          AATM(J,8) = ATM(MOD(J*K8,3),8)
                        ENDDO
                        NCF = NCF +1
                        CALL SORT_CFLOW(AATM,A,NATM,NEXT,CFLOW(1,NCF))
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END


      SUBROUTINE SORT_CFLOW(ATM,IA,NATM,NEXT,CFLOW)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS AUG 25 2012
C     ****************************************************
      IMPLICIT NONE
C     GLOBAL VARIABLES
C     CONSTANTS
C     ARGUMENTS 
      INTEGER NATM,NEXT
      INTEGER IA(*),CFLOW(*),ATM(2,*)
C     LOCAL VARIABLES 
      INTEGER I,J,CFPOS
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      CFPOS = 0
      DO I = 1,NATM
        DO J = 1,2
          IF (ATM(J,IA(I)).NE.0) THEN
            CFPOS = CFPOS +1
            CFLOW(CFPOS) = ATM(J,IA(I))
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END
