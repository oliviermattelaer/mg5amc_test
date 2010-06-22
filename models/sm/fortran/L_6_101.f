C     This File Automatically generated by MadGraph 5/FeynRules HELAS
C      writer 
C     The process calculated in this file is: 
C     number of lorentz index :1
      NUMBER OF SPIN INDEX :1
      OTHER INFO SET(['P2', 'S1', 'F3', 'M2'])
      (0, 0) --> ( S1_1 * ( ( F3_2 * ( 1J * P2_1 + (-1+0J) * P2_2 ) ) 
     $ + ( ( F3_1 * ( 1J * P2_0 + 1J * P2_3 ) ) + ( 1J * ( F3_3 
     $ * M2 ) ) ) ) )
      (1, 0) --> ( S1_1 * ( ( F3_1 * ( -1J * P2_1 + P2_2 ) ) 
     $ + ( ( F3_2 * ( -1J * P2_0 + -1J * P2_3 ) ) + ( 1J * ( F3_4 
     $ * M2 ) ) ) ) )
      (2, 0) --> ( S1_1 * ( ( F3_1 * ( P2_1 + 1J * P2_2 ) ) + ( ( F3_2
     $  * ( (-1+0J) * P2_0 + (-1+0J) * P2_3 ) ) + ( ( F3_4 * M2 ) ) ) 
     $ ) )
      (3, 0) --> ( S1_1 * ( ( F3_2 * ( 1J * P2_1 + (-1+0J) * P2_2 ) ) 
     $ + ( ( F3_1 * ( -1J * P2_0 + -1J * P2_3 ) ) + ( 1J * ( F3_3 
     $ * M2 ) ) ) ) )
      (0, 1) --> ( S1_1 * ( ( F3_1 * ( 1J * P2_1 + P2_2 ) ) + ( ( F3_2
     $  * ( 1J * P2_0 + -1J * P2_3 ) ) + ( 1J * ( F3_4 * M2 ) ) ) ) )
      (1, 1) --> ( S1_1 * ( ( F3_2 * ( -1J * P2_1 + (-1+0J) * P2_2 ) 
     $ ) + ( ( F3_1 * ( -1J * P2_0 + 1J * P2_3 ) ) + ( 1J * ( F3_3 
     $ * M2 ) ) ) ) )
      (2, 1) --> ( S1_1 * ( ( F3_2 * ( (-1+0J) * P2_1 + 1J * P2_2 ) ) 
     $ + ( ( F3_1 * ( P2_0 + (-1+0J) * P2_3 ) ) + ( (-1+0J) * ( F3_3 
     $ * M2 ) ) ) ) )
      (3, 1) --> ( S1_1 * ( ( F3_1 * ( -1J * P2_1 + (-1+0J) * P2_2 ) 
     $ ) + ( ( F3_2 * ( 1J * P2_0 + -1J * P2_3 ) ) + ( -1J * ( F3_4 
     $ * M2 ) ) ) ) )
      (0, 2) --> ( S1_1 * ( ( F3_4 * ( -1J * P2_1 + P2_2 ) ) 
     $ + ( ( F3_3 * ( 1J * P2_0 + -1J * P2_3 ) ) + ( 1J * ( F3_1 
     $ * M2 ) ) ) ) )
      (1, 2) --> ( S1_1 * ( ( F3_3 * ( -1J * P2_1 + P2_2 ) ) 
     $ + ( ( F3_4 * ( 1J * P2_0 + -1J * P2_3 ) ) + ( -1J * ( F3_2 
     $ * M2 ) ) ) ) )
      (2, 2) --> ( S1_1 * ( ( F3_3 * ( P2_1 + 1J * P2_2 ) ) + ( ( F3_4
     $  * ( P2_0 + (-1+0J) * P2_3 ) ) + ( (-1+0J) * ( F3_2 * M2 ) ) ) 
     $ ) )
      (3, 2) --> ( S1_1 * ( ( F3_4 * ( 1J * P2_1 + (-1+0J) * P2_2 ) ) 
     $ + ( ( F3_3 * ( 1J * P2_0 + -1J * P2_3 ) ) + ( -1J * ( F3_1 
     $ * M2 ) ) ) ) )
      (0, 3) --> ( S1_1 * ( ( F3_3 * ( -1J * P2_1 + (-1+0J) * P2_2 ) 
     $ ) + ( ( F3_4 * ( 1J * P2_0 + 1J * P2_3 ) ) + ( 1J * ( F3_2 
     $ * M2 ) ) ) ) )
      (1, 3) --> ( S1_1 * ( ( F3_4 * ( -1J * P2_1 + (-1+0J) * P2_2 ) 
     $ ) + ( ( F3_3 * ( 1J * P2_0 + 1J * P2_3 ) ) + ( -1J * ( F3_1 
     $ * M2 ) ) ) ) )
      (2, 3) --> ( S1_1 * ( ( F3_4 * ( (-1+0J) * P2_1 + 1J * P2_2 ) ) 
     $ + ( ( F3_3 * ( (-1+0J) * P2_0 + (-1+0J) * P2_3 ) ) + ( ( F3_1 
     $ * M2 ) ) ) ) )
      (3, 3) --> ( S1_1 * ( ( F3_3 * ( -1J * P2_1 + (-1+0J) * P2_2 ) 
     $ ) + ( ( F3_4 * ( -1J * P2_0 + -1J * P2_3 ) ) + ( 1J * ( F3_2 
     $ * M2 ) ) ) ) )
      / NUMBER OF LORENTZ INDEX :0
      NUMBER OF SPIN INDEX :0
      OTHER INFO SET(['P2', 'W2', 'M2'])
      (0,) --> ( ( M2 * ( -1 * M2 + 1J * W2 ) ) + ( ( P2_0**2 ) 
     $ + -1 * ( P2_1**2 ) + -1 * ( P2_2**2 ) + -1 * ( P2_3**2 ) ) )
      
      SUBROUTINE /USERS/OMATT/DOCUMENTS/ECLIPSE2/EXPORT_MODEL/MODELS
     $ /SM2/FORTRAN/L_6_101(C,F2,F3,M2,W2,F3)
      IMPLICIT NONE
      DOUBLE PRECISION C
      DOUBLE COMPLEX S1(3)
      DOUBLE COMPLEX F2(6)
      DOUBLE COMPLEX F3(6)
      DOUBLE COMPLEX DENOM
      DOUBLE PRECISION M2,W2
      DOUBLE PRECISION P2(0:3)
      F3(5)=-(-S1(2)-F2(5))
      F3(6)=-(-S1(3)-F2(6))
      P2(0) =  DBLE(F2(5))
      P2(1) =  DBLE(F2(6))
      P2(2) =  DIMAG(F2(6))
      P2(3) =  DIMAG(F2(5))
      DENOM =1D0/(((M2*(-M2+(0.000000000D0, 1.000000000D0)*W2))
     $ +((P2(0)**2)-(P2(1)**2)-(P2(2)**2)-(P2(3)**2))))
      F2(1)= C*DENOM*(S1(1)*((F3(2)*((0.000000000D0, 1.000000000D0)
     $ *P2(1)+-P2(2)))+((F3(1)*((0.000000000D0, 1.000000000D0)*P2(0)
     $ +(0.000000000D0, 1.000000000D0)*P2(3)))+((0.000000000D0
     $ , 1.000000000D0)*(F3(3)*M2)))))
      F2(2)= C*DENOM*(S1(1)*((F3(1)*((0.000000000D0, -1.000000000D0)
     $ *P2(1)+P2(2)))+((F3(2)*((0.000000000D0, -1.000000000D0)*P2(0)
     $ +(0.000000000D0, -1.000000000D0)*P2(3)))+((0.000000000D0
     $ , 1.000000000D0)*(F3(4)*M2)))))
      F2(3)= C*DENOM*(S1(1)*((F3(1)*(P2(1)+(0.000000000D0, 1.000000000D
     $ 0)*P2(2)))+((F3(2)*(-P2(0)+-P2(3)))+((F3(4)*M2)))))
      F2(4)= C*DENOM*(S1(1)*((F3(2)*((0.000000000D0, 1.000000000D0)
     $ *P2(1)+-P2(2)))+((F3(1)*((0.000000000D0, -1.000000000D0)*P2(0)
     $ +(0.000000000D0, -1.000000000D0)*P2(3)))+((0.000000000D0
     $ , 1.000000000D0)*(F3(3)*M2)))))
      F2(5)= C*DENOM*(S1(1)*((F3(1)*((0.000000000D0, 1.000000000D0)
     $ *P2(1)+P2(2)))+((F3(2)*((0.000000000D0, 1.000000000D0)*P2(0)
     $ +(0.000000000D0, -1.000000000D0)*P2(3)))+((0.000000000D0
     $ , 1.000000000D0)*(F3(4)*M2)))))
      F2(6)= C*DENOM*(S1(1)*((F3(2)*((0.000000000D0, -1.000000000D0)
     $ *P2(1)+-P2(2)))+((F3(1)*((0.000000000D0, -1.000000000D0)*P2(0)
     $ +(0.000000000D0, 1.000000000D0)*P2(3)))+((0.000000000D0
     $ , 1.000000000D0)*(F3(3)*M2)))))
      F2(7)= C*DENOM*(S1(1)*((F3(2)*(-P2(1)+(0.000000000D0, 1.000000000
     $ D0)*P2(2)))+((F3(1)*(P2(0)+-P2(3)))+(-(F3(3)*M2)))))
      F2(8)= C*DENOM*(S1(1)*((F3(1)*((0.000000000D0, -1.000000000D0)
     $ *P2(1)+-P2(2)))+((F3(2)*((0.000000000D0, 1.000000000D0)*P2(0)
     $ +(0.000000000D0, -1.000000000D0)*P2(3)))+((0.000000000D0, 
     $ -1.000000000D0)*(F3(4)*M2)))))
      F2(9)= C*DENOM*(S1(1)*((F3(4)*((0.000000000D0, -1.000000000D0)
     $ *P2(1)+P2(2)))+((F3(3)*((0.000000000D0, 1.000000000D0)*P2(0)
     $ +(0.000000000D0, -1.000000000D0)*P2(3)))+((0.000000000D0
     $ , 1.000000000D0)*(F3(1)*M2)))))
      F2(10)= C*DENOM*(S1(1)*((F3(3)*((0.000000000D0, -1.000000000D0)
     $ *P2(1)+P2(2)))+((F3(4)*((0.000000000D0, 1.000000000D0)*P2(0)
     $ +(0.000000000D0, -1.000000000D0)*P2(3)))+((0.000000000D0, 
     $ -1.000000000D0)*(F3(2)*M2)))))
      F2(11)= C*DENOM*(S1(1)*((F3(3)*(P2(1)+(0.000000000D0, 1.000000000
     $ D0)*P2(2)))+((F3(4)*(P2(0)+-P2(3)))+(-(F3(2)*M2)))))
      F2(12)= C*DENOM*(S1(1)*((F3(4)*((0.000000000D0, 1.000000000D0)
     $ *P2(1)+-P2(2)))+((F3(3)*((0.000000000D0, 1.000000000D0)*P2(0)
     $ +(0.000000000D0, -1.000000000D0)*P2(3)))+((0.000000000D0, 
     $ -1.000000000D0)*(F3(1)*M2)))))
      F2(13)= C*DENOM*(S1(1)*((F3(3)*((0.000000000D0, -1.000000000D0)
     $ *P2(1)+-P2(2)))+((F3(4)*((0.000000000D0, 1.000000000D0)*P2(0)
     $ +(0.000000000D0, 1.000000000D0)*P2(3)))+((0.000000000D0
     $ , 1.000000000D0)*(F3(2)*M2)))))
      F2(14)= C*DENOM*(S1(1)*((F3(4)*((0.000000000D0, -1.000000000D0)
     $ *P2(1)+-P2(2)))+((F3(3)*((0.000000000D0, 1.000000000D0)*P2(0)
     $ +(0.000000000D0, 1.000000000D0)*P2(3)))+((0.000000000D0, 
     $ -1.000000000D0)*(F3(1)*M2)))))
      F2(15)= C*DENOM*(S1(1)*((F3(4)*(-P2(1)+(0.000000000D0, 1.00000000
     $ 0D0)*P2(2)))+((F3(3)*(-P2(0)+-P2(3)))+((F3(1)*M2)))))
      F2(16)= C*DENOM*(S1(1)*((F3(3)*((0.000000000D0, -1.000000000D0)
     $ *P2(1)+-P2(2)))+((F3(4)*((0.000000000D0, -1.000000000D0)*P2(0)
     $ +(0.000000000D0, -1.000000000D0)*P2(3)))+((0.000000000D0
     $ , 1.000000000D0)*(F3(2)*M2)))))
      END
