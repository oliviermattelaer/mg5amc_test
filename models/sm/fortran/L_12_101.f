C     This File Automatically generated by MadGraph 5/FeynRules HELAS
C      writer 
C     The process calculated in this file is: 
C     number of lorentz index :0
      NUMBER OF SPIN INDEX :0
      OTHER INFO SET(['P2', 'V1', 'S3', 'P3'])
      (0,) --> ( S3_1 * ( ( V1_1 * ( 1J * P2_0 + -1J * P3_0 ) ) 
     $ + ( ( V1_2 * ( -1J * P2_1 + 1J * P3_1 ) ) + ( ( V1_3 * ( 
     $ -1J * P2_2 + 1J * P3_2 ) ) + ( V1_4 * ( -1J * P2_3 + 1J 
     $ * P3_3 ) ) ) ) ) )
      / NUMBER OF LORENTZ INDEX :0
      NUMBER OF SPIN INDEX :0
      OTHER INFO SET(['P2', 'W2', 'M2'])
      (0,) --> ( ( M2 * ( -1 * M2 + 1J * W2 ) ) + ( ( P2_0**2 ) 
     $ + -1 * ( P2_1**2 ) + -1 * ( P2_2**2 ) + -1 * ( P2_3**2 ) ) )
      
      SUBROUTINE /USERS/OMATT/DOCUMENTS/ECLIPSE2/EXPORT_MODEL/MODELS
     $ /SM2/FORTRAN/L_12_101(C,V1,S2,M2,W2,S3)
      IMPLICIT NONE
      DOUBLE PRECISION C
      DOUBLE COMPLEX V1(6)
      DOUBLE COMPLEX S2(3)
      DOUBLE COMPLEX S3(3)
      DOUBLE COMPLEX DENOM
      DOUBLE PRECISION M2,W2
      DOUBLE PRECISION P2(0:3),P3(0:3)
      S3(2)=-(-V1(5)-S2(2))
      S3(3)=-(-V1(6)-S2(3))
      P2(0) =  DBLE(S2(2))
      P2(1) =  DBLE(S2(3))
      P2(2) =  DIMAG(S2(3))
      P2(3) =  DIMAG(S2(2))
      P3(0) =  DBLE(S3(2))
      P3(1) =  DBLE(S3(3))
      P3(2) =  DIMAG(S3(3))
      P3(3) =  DIMAG(S3(2))
      DENOM =1D0/(((M2*(-M2+(0.000000000D0, 1.000000000D0)*W2))
     $ +((P2(0)**2)-(P2(1)**2)-(P2(2)**2)-(P2(3)**2))))
      S2(1)= C*DENOM*(S3(1)*((V1(1)*((0.000000000D0, 1.000000000D0)
     $ *P2(0)+(0.000000000D0, -1.000000000D0)*P3(0)))+((V1(2)
     $ *((0.000000000D0, -1.000000000D0)*P2(1)+(0.000000000D0
     $ , 1.000000000D0)*P3(1)))+((V1(3)*((0.000000000D0, -1.000000000D0
     $ )*P2(2)+(0.000000000D0, 1.000000000D0)*P3(2)))+(V1(4)*((0.000000
     $ 000D0, -1.000000000D0)*P2(3)+(0.000000000D0, 1.000000000D0)
     $ *P3(3)))))))
      END
