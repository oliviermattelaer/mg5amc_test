C     This File Automatically generated by MadGraph 5/FeynRules HELAS
C      writer 
C     The process calculated in this file is: 
C     number of lorentz index :0
      NUMBER OF SPIN INDEX :0
      OTHER INFO SET(['S1', 'V3', 'V4'])
      (0,) --> ( S1_1 * ( 1J * ( V4_1 * V3_1 ) + -1J * ( V4_2 
     $ * V3_2 ) + -1J * ( V4_3 * V3_3 ) + -1J * ( V4_4 * V3_4 ) ) )
      / NUMBER OF LORENTZ INDEX :0
      NUMBER OF SPIN INDEX :0
      OTHER INFO SET(['P2', 'W2', 'M2'])
      (0,) --> ( ( M2 * ( -1 * M2 + 1J * W2 ) ) + ( ( P2_0**2 ) 
     $ + -1 * ( P2_1**2 ) + -1 * ( P2_2**2 ) + -1 * ( P2_3**2 ) ) )
      
      SUBROUTINE /USERS/OMATT/DOCUMENTS/ECLIPSE2/EXPORT_MODEL/MODELS
     $ /SM2/FORTRAN/L_21_1011(C,V3,V4,S1,M2,W2,V4)
      IMPLICIT NONE
      DOUBLE PRECISION C
      DOUBLE COMPLEX S1(3)
      DOUBLE COMPLEX S2(3)
      DOUBLE COMPLEX V3(6)
      DOUBLE COMPLEX V4(6)
      DOUBLE COMPLEX DENOM
      DOUBLE PRECISION M2,W2
      DOUBLE PRECISION P2(0:3)
      V4(5)=-(-S1(2)-S2(2)-V3(5))
      V4(6)=-(-S1(3)-S2(3)-V3(6))
      P2(0) =  DBLE(S2(2))
      P2(1) =  DBLE(S2(3))
      P2(2) =  DIMAG(S2(3))
      P2(3) =  DIMAG(S2(2))
      DENOM =1D0/(((M2*(-M2+(0.000000000D0, 1.000000000D0)*W2))
     $ +((P2(0)**2)-(P2(1)**2)-(P2(2)**2)-(P2(3)**2))))
      S2(1)= C*DENOM*(S1(1)*((0.000000000D0, 1.000000000D0)*(V4(1)
     $ *V3(1))+(0.000000000D0, -1.000000000D0)*(V4(2)*V3(2))+(0.0000000
     $ 00D0, -1.000000000D0)*(V4(3)*V3(3))+(0.000000000D0, -1.000000000
     $ D0)*(V4(4)*V3(4))))
      END
