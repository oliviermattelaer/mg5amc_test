C     This File Automatically generated by MadGraph 5/FeynRules HELAS
C      writer 
C     The process calculated in this file is: 
C     number of lorentz index :1
      NUMBER OF SPIN INDEX :0
      OTHER INFO SET(['F2', 'F3', 'mass1', 'OM1', 'P1'])
      (0,) --> ( ( OM1 * ( ( F2_3 * ( ( F3_2 * ( 1J * P1_1 + P1_2 ) ) 
     $ + ( F3_1 * ( 1J * P1_0 + 1J * P1_3 ) ) ) ) + ( F2_4 * ( ( F3_1 
     $ * ( 1J * P1_1 + (-1+0J) * P1_2 ) ) + ( F3_2 * ( 1J * P1_0 
     $ + -1J * P1_3 ) ) ) ) ) * P1_0 ) + ( -1J * ( F3_1 * F2_3 ) 
     $ + -1J * ( F3_2 * F2_4 ) ) )
      (1,) --> ( ( OM1 * ( ( F2_3 * ( ( F3_2 * ( 1J * P1_1 + P1_2 ) ) 
     $ + ( F3_1 * ( 1J * P1_0 + 1J * P1_3 ) ) ) ) + ( F2_4 * ( ( F3_1 
     $ * ( 1J * P1_1 + (-1+0J) * P1_2 ) ) + ( F3_2 * ( 1J * P1_0 
     $ + -1J * P1_3 ) ) ) ) ) * P1_1 ) + ( 1J * ( F3_2 * F2_3 ) 
     $ + 1J * ( F3_1 * F2_4 ) ) )
      (2,) --> ( ( OM1 * ( ( F2_3 * ( ( F3_2 * ( 1J * P1_1 + P1_2 ) ) 
     $ + ( F3_1 * ( 1J * P1_0 + 1J * P1_3 ) ) ) ) + ( F2_4 * ( ( F3_1 
     $ * ( 1J * P1_1 + (-1+0J) * P1_2 ) ) + ( F3_2 * ( 1J * P1_0 
     $ + -1J * P1_3 ) ) ) ) ) * P1_2 ) + ( ( F3_2 * F2_3 ) + (-1
     $ +0J) * ( F3_1 * F2_4 ) ) )
      (3,) --> ( ( OM1 * ( ( F2_3 * ( ( F3_2 * ( 1J * P1_1 + P1_2 ) ) 
     $ + ( F3_1 * ( 1J * P1_0 + 1J * P1_3 ) ) ) ) + ( F2_4 * ( ( F3_1 
     $ * ( 1J * P1_1 + (-1+0J) * P1_2 ) ) + ( F3_2 * ( 1J * P1_0 
     $ + -1J * P1_3 ) ) ) ) ) * P1_3 ) + ( 1J * ( F3_1 * F2_3 ) 
     $ + -1J * ( F3_2 * F2_4 ) ) )
      / NUMBER OF LORENTZ INDEX :0
      NUMBER OF SPIN INDEX :0
      OTHER INFO SET(['P1', 'M1', 'W1'])
      (0,) --> ( ( M1 * ( -1 * M1 + 1J * W1 ) ) + ( ( P1_0**2 ) 
     $ + -1 * ( P1_1**2 ) + -1 * ( P1_2**2 ) + -1 * ( P1_3**2 ) ) )
      
      SUBROUTINE /USERS/OMATT/DOCUMENTS/ECLIPSE2/EXPORT_MODEL/MODELS
     $ /SM2/FORTRAN/L_15_011(C,F2,F3,M1,W1,F3)
      IMPLICIT NONE
      DOUBLE PRECISION C
      DOUBLE COMPLEX V1(6)
      DOUBLE COMPLEX F2(6)
      DOUBLE COMPLEX F3(6)
      DOUBLE COMPLEX DENOM
      DOUBLE PRECISION M1,W1
      DOUBLE COMPLEX OM1
      DOUBLE PRECISION P1(0:3)
      F3(5)=-(-V1(5)-F2(5))
      F3(6)=-(-V1(6)-F2(6))
      P1(0) =  DBLE(V1(5))
      P1(1) =  DBLE(V1(6))
      P1(2) =  DIMAG(V1(6))
      P1(3) =  DIMAG(V1(5))
      OM1 = 0D0
      IF (M1 .NE. 0D0) OM1=1D0/DCMPLX(M1**2,-W1*M1)
      DENOM =1D0/(((M1*(-M1+(0.000000000D0, 1.000000000D0)*W1))
     $ +((P1(0)**2)-(P1(1)**2)-(P1(2)**2)-(P1(3)**2))))
      V1(1)= C*DENOM*((OM1*((F2(3)*((F3(2)*((0.000000000D0, 1.000000000
     $ D0)*P1(1)+P1(2)))+(F3(1)*((0.000000000D0, 1.000000000D0)*P1(0)
     $ +(0.000000000D0, 1.000000000D0)*P1(3)))))+(F2(4)*((F3(1)
     $ *((0.000000000D0, 1.000000000D0)*P1(1)+-P1(2)))+(F3(2)
     $ *((0.000000000D0, 1.000000000D0)*P1(0)+(0.000000000D0, 
     $ -1.000000000D0)*P1(3))))))*P1(0))+((0.000000000D0, -1.000000000D
     $ 0)*(F3(1)*F2(3))+(0.000000000D0, -1.000000000D0)*(F3(2)*F2(4))))
      V1(2)= C*DENOM*((OM1*((F2(3)*((F3(2)*((0.000000000D0, 1.000000000
     $ D0)*P1(1)+P1(2)))+(F3(1)*((0.000000000D0, 1.000000000D0)*P1(0)
     $ +(0.000000000D0, 1.000000000D0)*P1(3)))))+(F2(4)*((F3(1)
     $ *((0.000000000D0, 1.000000000D0)*P1(1)+-P1(2)))+(F3(2)
     $ *((0.000000000D0, 1.000000000D0)*P1(0)+(0.000000000D0, 
     $ -1.000000000D0)*P1(3))))))*P1(1))+((0.000000000D0, 1.000000000D0
     $ )*(F3(2)*F2(3))+(0.000000000D0, 1.000000000D0)*(F3(1)*F2(4))))
      V1(3)= C*DENOM*((OM1*((F2(3)*((F3(2)*((0.000000000D0, 1.000000000
     $ D0)*P1(1)+P1(2)))+(F3(1)*((0.000000000D0, 1.000000000D0)*P1(0)
     $ +(0.000000000D0, 1.000000000D0)*P1(3)))))+(F2(4)*((F3(1)
     $ *((0.000000000D0, 1.000000000D0)*P1(1)+-P1(2)))+(F3(2)
     $ *((0.000000000D0, 1.000000000D0)*P1(0)+(0.000000000D0, 
     $ -1.000000000D0)*P1(3))))))*P1(2))+((F3(2)*F2(3))+-(F3(1)
     $ *F2(4))))
      V1(4)= C*DENOM*((OM1*((F2(3)*((F3(2)*((0.000000000D0, 1.000000000
     $ D0)*P1(1)+P1(2)))+(F3(1)*((0.000000000D0, 1.000000000D0)*P1(0)
     $ +(0.000000000D0, 1.000000000D0)*P1(3)))))+(F2(4)*((F3(1)
     $ *((0.000000000D0, 1.000000000D0)*P1(1)+-P1(2)))+(F3(2)
     $ *((0.000000000D0, 1.000000000D0)*P1(0)+(0.000000000D0, 
     $ -1.000000000D0)*P1(3))))))*P1(3))+((0.000000000D0, 1.000000000D0
     $ )*(F3(1)*F2(3))+(0.000000000D0, -1.000000000D0)*(F3(2)*F2(4))))
      END
