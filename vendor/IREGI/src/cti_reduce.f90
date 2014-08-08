MODULE cti_reduce
  USE global
  USE csi_reduce
  USE cpave_reduce
  USE funlib
  USE ti_reduce
  IMPLICIT NONE
CONTAINS
  SUBROUTINE comp_tensor_integral_reduce(IMODE,NLOOPLINE,MAXRANK,NCOEFS,PDEN_IN,M2L_IN,MU,TICOEFS,TSTABLE,CTMODE)
    ! CTMODE=1 ---> DO NOTHING
    ! CTMODE=2 ---> ROTATE THE LOOP PROPAGATOR DIRECTION,i.e D0...DN-1 to DN-1...D0
    ! IMODE=0, IBP reduction
    ! IMODE=1, PaVe reduction
    ! IMODE=2, PaVe reduction with stablility improved by IBP reduction 
    IMPLICIT NONE
    INTEGER,INTENT(IN)::NLOOPLINE,MAXRANK,IMODE,NCOEFS
    INTEGER,INTENT(IN),OPTIONAL::CTMODE
    REAL(KIND(1d0)),INTENT(IN)::MU
    REAL(KIND(1d0)),DIMENSION(NLOOPLINE,0:3),INTENT(IN)::PDEN_IN
    REAL(KIND(1d0)),DIMENSION(NLOOPLINE,0:3)::PDEN
    REAL(KIND(1d0)),DIMENSION(0:3)::PDEN_N
    COMPLEX(KIND(1d0)),DIMENSION(NLOOPLINE),INTENT(IN)::M2L_IN
    COMPLEX(KIND(1d0)),DIMENSION(NLOOPLINE)::M2L
    COMPLEX(KIND(1d0)),DIMENSION(0:NCOEFS-1,1:4),INTENT(OUT)::TICOEFS
    COMPLEX(KIND(1d0)),DIMENSION(0:NCOEFS-1,1:4)::TICOEFS_SAVE
    LOGICAL,INTENT(OUT)::TSTABLE
    INTEGER,DIMENSION(1,1)::sol11
    REAL(KIND(1d0)),DIMENSION(1)::factor1
    INTEGER::first=0,i,j,numzerp,n,r,nr,init,ntot
    INTEGER,DIMENSION(NLOOPLINE)::zerp
    REAL(KIND(1d0)),DIMENSION(xiarraymax2)::syfactor
    INTEGER,DIMENSION(xiarraymax2,-1:NLOOPLINE)::sy
!    REAL(KIND(1d0)),PARAMETER::EPS=1d-10
    REAL(KIND(1d0))::temp
    INTEGER::ctmode_local
    LOGICAL::do_shift
    SAVE first
    IF(PRESENT(CTMODE))THEN
       ctmode_local=CTMODE
    ELSE
       ctmode_local=1
    ENDIF
    MU_R_IREGI=MU
    IF(first.EQ.0)THEN
       WRITE(*,*)"#####################################################################################"
       WRITE(*,*)"#                                                                                   #"
       WRITE(*,*)"#                           IREGI-alpha-1.0.0                                       #"
       WRITE(*,*)"#    package for one-loop tensor Integral REduction with General propagator Indices #"
       WRITE(*,*)"#    Author: Hua-Sheng Shao (huasheng.shao@cern.ch/erdissshaw@gmail.com)            #"
       WRITE(*,*)"#    a) Physics School, Peking University, Beijing, China                           #"
       WRITE(*,*)"#    b) PH Department, TH Unit, CERN, Geneva, Switzerland                           #"
       WRITE(*,*)"#                                                                                   #"
       WRITE(*,*)"#####################################################################################"
       ! initialization xiarray and metric
       CALL all_Integers(1,1,1,sol11,factor1)
       DO i=0,3
          DO j=0,3
             IF(i.NE.j)THEN
                metric(i,j)=0d0
             ELSEIF(i.EQ.0)THEN
                metric(i,j)=1d0
             ELSE
                metric(i,j)=-1d0
             ENDIF
          ENDDO
       END DO
       first=1
    ENDIF
    IF(ctmode_local.EQ.2)THEN
       PDEN_N(0:3)=PDEN_IN(NLOOPLINE,0:3)
       DO I=1,NLOOPLINE
          PDEN(NLOOPLINE-I+1,0:3)=PDEN_IN(I,0:3)-PDEN_N(0:3)
          M2L(NLOOPLINE-I+1)=M2L_IN(I)
       ENDDO
    ELSE
       PDEN_N(0:3)=PDEN_IN(1,0:3)
       DO I=1,NLOOPLINE
          PDEN(I,0:3)=PDEN_IN(I,0:3)-PDEN_N(0:3)
          M2L(1:NLOOPLINE)=M2L_IN(1:NLOOPLINE)
       ENDDO
    ENDIF
    numzerp=0
    temp=(ABS(PDEN_N(0))+ABS(PDEN_N(1))+ABS(PDEN_N(2))+ABS(PDEN_N(3)))/4d0
    IF(temp.LT.EPS)THEN
       do_shift=.FALSE.
    ELSE
       do_shift=.TRUE.
    ENDIF
    DO i=1,NLOOPLINE
       temp=(ABS(PDEN(i,0))+ABS(PDEN(i,1))+ABS(PDEN(i,2))+ABS(PDEN(i,3)))/4d0
       IF(temp.LT.EPS)THEN
          zerp(i)=0
          numzerp=numzerp+1
       ELSE
          zerp(i)=1
       ENDIF
    ENDDO
    n=NLOOPLINE-numzerp
    IF(n.GE.6.OR.MAXRANK.GE.6)THEN
       WRITE(*,*)"ERROR: out of range of comp_tensor_integral_reduce (N<8,R<7)"
       STOP
    ENDIF
    CALL sytensor(n,MAXRANK,ntot,sy(1:xiarraymax2,-1:n),syfactor(1:xiarraymax2))
    CALL comp_symmetry(IMODE,NLOOPLINE,ntot,sy(1:xiarraymax2,-1:n),syfactor(1:xiarraymax2),&
         numzerp,zerp,PDEN,M2L,NCOEFS,TICOEFS)
    TSTABLE=STABLE_IREGI
    IF(TSTABLE.AND.do_shift)THEN
       ! do momentum shift
       TICOEFS_SAVE(0:NCOEFS-1,1:4)=TICOEFS(0:NCOEFS-1,1:4)
       PDEN_N(0:3)=-PDEN_N(0:3)
       CALL SHIFT_MOM(PDEN_N,NCOEFS,TICOEFS_SAVE,TICOEFS)
    ENDIF
    RETURN
  END SUBROUTINE comp_tensor_integral_reduce

  SUBROUTINE comp_symmetry(IMODE,NLOOPLINE,ntot,sy,syfactor,numzerp,zerp,PDEN,M2L,NCOEFS,coefs)
    ! IMODE=0, IBP reduction  
    ! IMODE=1, PaVe reduction
    ! IMODE=2, PaVe reduction with stablility improved by IBP reduction 
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ntot,numzerp,NLOOPLINE,IMODE,NCOEFS
    INTEGER,DIMENSION(xiarraymax2,-1:NLOOPLINE-numzerp),INTENT(IN)::sy
    REAL(KIND(1d0)),DIMENSION(xiarraymax2),INTENT(IN)::syfactor
    INTEGER,DIMENSION(NLOOPLINE),INTENT(IN)::zerp
    REAL(KIND(1d0)),DIMENSION(NLOOPLINE,0:3),INTENT(IN)::PDEN
    COMPLEX(KIND(1d0)),DIMENSION(NLOOPLINE),INTENT(IN)::M2L
    INTEGER::idim
    INTEGER,DIMENSION(NLOOPLINE)::indices
    INTEGER,DIMENSION(0:NLOOPLINE)::paveindices
    COMPLEX(KIND(1d0)),DIMENSION(0:NCOEFS-1,1:4),INTENT(OUT)::coefs
    INTEGER,DIMENSION(xiarraymax,4)::sol
    REAL(KIND(1d0)),DIMENSION(xiarraymax)::factor
    REAL(KIND(1d0)),DIMENSION(NLOOPLINE-numzerp,0:3)::PCLL
    INTEGER::i,j,k,init,r,nntot,nt
    LOGICAL::lzero
    INTEGER,DIMENSION(0:NLOOPLINE)::nj
    INTEGER,DIMENSION(0:NLOOPLINE,1:5)::jlist
    REAL(KIND(1d0))::res
    COMPLEX(KIND(1d0)),DIMENSION(1:4)::scalar
    REAL(KIND(1d0)),DIMENSION(xiarraymax3)::coco
    INTEGER,DIMENSION(5)::lor
    coefs(0:NCOEFS-1,1:4)=DCMPLX(0d0)
    j=1
    DO i=1,NLOOPLINE
       IF(zerp(i).NE.0)THEN
          PCLL(j,0:3)=PDEN(i,0:3)
          j=j+1
       ENDIF
    ENDDO
    nt=NLOOPLINE-numzerp
    DO i=1,ntot
       r=sy(i,-1)
       nntot=ntot_xiarray(r,4)
       CALL all_integers(4,nntot,r,sol(1:nntot,1:4),factor(1:nntot))
       lzero=.TRUE.
       coco(1:nntot)=0d0
       DO j=1,nntot
          res=0d0
          init=1
          IF(sol(j,1).GE.1)THEN
             lor(init:init+sol(j,1)-1)=0
             init=init+sol(j,1)
          ENDIF
          IF(sol(j,2).GE.1)THEN
             lor(init:init+sol(j,2)-1)=1
             init=init+sol(j,2)
          ENDIF
          IF(sol(j,3).GE.1)THEN
             lor(init:init+sol(j,3)-1)=2
             init=init+sol(j,3)
          ENDIF
          IF(sol(j,4).GE.1)THEN
             lor(init:init+sol(j,4)-1)=3
             init=init+sol(j,4)
          ENDIF
          nj(0:nt)=0
          jlist(0:nt,1:5)=0
          CALL recursive_symmetry(nt,1,r,PCLL(1:nt,0:3),sy(i,-1:nt),lor(1:r),&
               nj(0:nt),jlist(0:nt,1:5),res)
          lzero=lzero.AND.(ABS(res).LT.EPS)
          IF(.NOT.ML5_CONVENTION)THEN
             coco(j)=res*factor(j)
          ELSE
             coco(j)=res
          ENDIF
          !coco(j)=res*factor(j)
       ENDDO
       IF(.NOT.lzero)THEN
          IF(IMODE.EQ.0)THEN
             ! IBP reduction
             indices(1:NLOOPLINE)=1 ! canonical form
             idim=2*(r-sy(i,0))
             init=1
             DO j=1,NLOOPLINE
                IF(zerp(j).NE.1)CYCLE
                indices(j)=indices(j)+sy(i,init)
                init=init+1
             ENDDO
             scalar(1:4)=comp_scalar_integral_reduce(NLOOPLINE,idim,indices,PDEN,M2L)
             IF(.NOT.STABLE_IREGI)THEN
                WRITE(*,*)"IREGI:WARNING, it detects unstable case, some integrals may set to be 0."
                RETURN
             ENDIF
             DO j=1,nntot
                init=calc_pos(sol(j,1:4))
                coefs(init,1:4)=coefs(init,1:4)+coco(j)*syfactor(i)*scalar(1:4)
             ENDDO
          ELSE
             ! PaVe reduction
             paveindices(0)=2*sy(i,0)
             init=1
             DO j=1,NLOOPLINE
                IF(zerp(j).NE.1)THEN
                   paveindices(j)=0
                ELSE
                   paveindices(j)=sy(i,init)
                   init=init+1
                ENDIF
             ENDDO
             IF(IMODE.EQ.1)THEN
                ! pure PaVe reduction
                scalar(1:4)=comp_pavefun_reduce(NLOOPLINE,paveindices,PDEN,M2L)
             ELSE
                ! optimized PaVe reduction
                ! the stability has been improved by IBP reduction
                scalar(1:4)=comp_pave_opt_reduce(NLOOPLINE,paveindices,PDEN,M2L)
             ENDIF
             IF(.NOT.STABLE_IREGI)THEN
                WRITE(*,*)"IREGI:WARNING, it detects unstable case, some integrals may set to be 0."
                RETURN
             ENDIF
             DO j=1,nntot
                init=calc_pos(sol(j,1:4))
                coefs(init,1:4)=coefs(init,1:4)+coco(j)*scalar(1:4)
             ENDDO
          ENDIF
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE comp_symmetry

END MODULE cti_reduce
