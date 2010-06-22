C     This File Automatically generated by MadGraph 5/FeynRules HELAS
C      writer 
C     The process calculated in this file is: 
C     number of lorentz index :0
      NUMBER OF SPIN INDEX :0
      OTHER INFO SET(['P2', 'P3', 'V3', 'S2'])
      (0,) --> ( S2_1 * ( ( V3_1 * ( 1J * P2_0 + 1J * P3_0 ) ) 
     $ + ( ( V3_2 * ( -1J * P2_1 + -1J * P3_1 ) ) + ( ( V3_3 * ( 
     $ -1J * P2_2 + -1J * P3_2 ) ) + ( V3_4 * ( -1J * P2_3 + 
     $ -1J * P3_3 ) ) ) ) ) )
      / NUMBER OF LORENTZ INDEX :0
      NUMBER OF SPIN INDEX :0
      OTHER INFO SET(['P1', 'M1', 'W1'])
      (0,) --> ( ( M1 * ( -1 * M1 + 1J * W1 ) ) + ( ( P1_0**2 ) 
     $ + -1 * ( P1_1**2 ) + -1 * ( P1_2**2 ) + -1 * ( P1_3**2 ) ) )
      
      SUBROUTINE /USERS/OMATT/DOCUMENTS/ECLIPSE2/EXPORT_MODEL/MODELS
     $ /SM2/FORTRAN/L_3_011(C,V3,S1,M1,W1,V3)
      IMPLICIT NONE
      DOUBLE PRECISION C
      DOUBLE COMPLEX S1(3)
      DOUBLE COMPLEX S2(3)
      DOUBLE COMPLEX V3(6)
      DOUBLE COMPLEX DENOM
      DOUBLE PRECISION M1,W1
      DOUBLE PRECISION P2(0:3),P3(0:3),P1(0:3)
      V3(5)=-(-S1(2)-S2(2))
      V3(6)=-(-S1(3)-S2(3))
      P2(0) =  DBLE(S2(2))
      P2(1) =  DBLE(S2(3))
      P2(2) =  DIMAG(S2(3))
      P2(3) =  DIMAG(S2(2))
      P3(0) = - DBLE(V3(5))
      P3(1) = - DBLE(V3(6))
      P3(2) = - DIMAG(V3(6))
      P3(3) = - DIMAG(V3(5))
      P1(0) =  DBLE(S1(2))
      P1(1) =  DBLE(S1(3))
      P1(2) =  DIMAG(S1(3))
      P1(3) =  DIMAG(S1(2))
      DENOM =1D0/(((M1*(-M1+(0.000000000D0, 1.000000000D0)*W1))
     $ +((P1(0)**2)-(P1(1)**2)-(P1(2)**2)-(P1(3)**2))))
      S1(1)= C*DENOM*(S2(1)*((V3(1)*((0.000000000D0, 1.000000000D0)
     $ *P2(0)+(0.000000000D0, 1.000000000D0)*P3(0)))+((V3(2)*((0.000000
     $ 000D0, -1.000000000D0)*P2(1)+(0.000000000D0, -1.000000000D0)
     $ *P3(1)))+((V3(3)*((0.000000000D0, -1.000000000D0)*P2(2)
     $ +(0.000000000D0, -1.000000000D0)*P3(2)))+(V3(4)*((0.000000000D0
     $ , -1.000000000D0)*P2(3)+(0.000000000D0, -1.000000000D0)
     $ *P3(3)))))))
      END
