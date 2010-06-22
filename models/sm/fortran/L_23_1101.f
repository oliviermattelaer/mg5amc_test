C     This File Automatically generated by MadGraph 5/FeynRules HELAS
C      writer 
C     The process calculated in this file is: 
C     number of lorentz index :1
      NUMBER OF SPIN INDEX :0
      OTHER INFO SET(['P1', 'mass1', 'S3', 'V2', 'S4', 'OM1'])
      (0,) --> ( S3_1 * ( ( OM1 * ( 1J * ( V2_1 * P1_0 ) + -1J 
     $ * ( V2_2 * P1_1 ) + -1J * ( V2_3 * P1_2 ) + -1J * ( V2_4 
     $ * P1_3 ) ) * P1_0 ) + ( -1J * V2_1 ) ) * S4_1 )
      (1,) --> ( S3_1 * ( ( OM1 * ( 1J * ( V2_1 * P1_0 ) + -1J 
     $ * ( V2_2 * P1_1 ) + -1J * ( V2_3 * P1_2 ) + -1J * ( V2_4 
     $ * P1_3 ) ) * P1_1 ) + ( -1J * V2_2 ) ) * S4_1 )
      (2,) --> ( S3_1 * ( ( OM1 * ( 1J * ( V2_1 * P1_0 ) + -1J 
     $ * ( V2_2 * P1_1 ) + -1J * ( V2_3 * P1_2 ) + -1J * ( V2_4 
     $ * P1_3 ) ) * P1_2 ) + ( -1J * V2_3 ) ) * S4_1 )
      (3,) --> ( S3_1 * ( ( OM1 * ( 1J * ( V2_1 * P1_0 ) + -1J 
     $ * ( V2_2 * P1_1 ) + -1J * ( V2_3 * P1_2 ) + -1J * ( V2_4 
     $ * P1_3 ) ) * P1_3 ) + ( -1J * V2_4 ) ) * S4_1 )
      / NUMBER OF LORENTZ INDEX :0
      NUMBER OF SPIN INDEX :0
      OTHER INFO SET(['P1', 'M1', 'W1'])
      (0,) --> ( ( M1 * ( -1 * M1 + 1J * W1 ) ) + ( ( P1_0**2 ) 
     $ + -1 * ( P1_1**2 ) + -1 * ( P1_2**2 ) + -1 * ( P1_3**2 ) ) )
      
      SUBROUTINE /USERS/OMATT/DOCUMENTS/ECLIPSE2/EXPORT_MODEL/MODELS
     $ /SM2/FORTRAN/L_23_1101(C,V1,V2,S3,M1,W1,S4)
      IMPLICIT NONE
      DOUBLE PRECISION C
      DOUBLE COMPLEX V1(6)
      DOUBLE COMPLEX V2(6)
      DOUBLE COMPLEX S3(3)
      DOUBLE COMPLEX S4(3)
      DOUBLE COMPLEX DENOM
      DOUBLE PRECISION M1,W1
      DOUBLE COMPLEX OM1
      DOUBLE PRECISION P1(0:3)
      S4(2)=-(-V1(5)-V2(5)-S3(2))
      S4(3)=-(-V1(6)-V2(6)-S3(3))
      P1(0) =  DBLE(V1(5))
      P1(1) =  DBLE(V1(6))
      P1(2) =  DIMAG(V1(6))
      P1(3) =  DIMAG(V1(5))
      OM1 = 0D0
      IF (M1 .NE. 0D0) OM1=1D0/DCMPLX(M1**2,-W1*M1)
      DENOM =1D0/(((M1*(-M1+(0.000000000D0, 1.000000000D0)*W1))
     $ +((P1(0)**2)-(P1(1)**2)-(P1(2)**2)-(P1(3)**2))))
      S3(1)= C*DENOM*(S3(1)*((OM1*((0.000000000D0, 1.000000000D0)
     $ *(V2(1)*P1(0))+(0.000000000D0, -1.000000000D0)*(V2(2)*P1(1))
     $ +(0.000000000D0, -1.000000000D0)*(V2(3)*P1(2))+(0.000000000D0, 
     $ -1.000000000D0)*(V2(4)*P1(3)))*P1(0))+((0.000000000D0, 
     $ -1.000000000D0)*V2(1)))*S4(1))
      S3(2)= C*DENOM*(S3(1)*((OM1*((0.000000000D0, 1.000000000D0)
     $ *(V2(1)*P1(0))+(0.000000000D0, -1.000000000D0)*(V2(2)*P1(1))
     $ +(0.000000000D0, -1.000000000D0)*(V2(3)*P1(2))+(0.000000000D0, 
     $ -1.000000000D0)*(V2(4)*P1(3)))*P1(1))+((0.000000000D0, 
     $ -1.000000000D0)*V2(2)))*S4(1))
      S3(3)= C*DENOM*(S3(1)*((OM1*((0.000000000D0, 1.000000000D0)
     $ *(V2(1)*P1(0))+(0.000000000D0, -1.000000000D0)*(V2(2)*P1(1))
     $ +(0.000000000D0, -1.000000000D0)*(V2(3)*P1(2))+(0.000000000D0, 
     $ -1.000000000D0)*(V2(4)*P1(3)))*P1(2))+((0.000000000D0, 
     $ -1.000000000D0)*V2(3)))*S4(1))
      S3(4)= C*DENOM*(S3(1)*((OM1*((0.000000000D0, 1.000000000D0)
     $ *(V2(1)*P1(0))+(0.000000000D0, -1.000000000D0)*(V2(2)*P1(1))
     $ +(0.000000000D0, -1.000000000D0)*(V2(3)*P1(2))+(0.000000000D0, 
     $ -1.000000000D0)*(V2(4)*P1(3)))*P1(3))+((0.000000000D0, 
     $ -1.000000000D0)*V2(4)))*S4(1))
      END
