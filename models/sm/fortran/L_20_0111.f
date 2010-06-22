C     This File Automatically generated by MadGraph 5/FeynRules HELAS
C      writer 
C     The process calculated in this file is: 
C     number of lorentz index :0
      NUMBER OF SPIN INDEX :0
      OTHER INFO SET(['S3', 'S2', 'S4'])
      (0,) --> 1J * ( S4_1 * S3_1 * S2_1 )
      / NUMBER OF LORENTZ INDEX :0
      NUMBER OF SPIN INDEX :0
      OTHER INFO SET(['P1', 'M1', 'W1'])
      (0,) --> ( ( M1 * ( -1 * M1 + 1J * W1 ) ) + ( ( P1_0**2 ) 
     $ + -1 * ( P1_1**2 ) + -1 * ( P1_2**2 ) + -1 * ( P1_3**2 ) ) )
      
      SUBROUTINE /USERS/OMATT/DOCUMENTS/ECLIPSE2/EXPORT_MODEL/MODELS
     $ /SM2/FORTRAN/L_20_0111(C,S1,S2,S3,M1,W1,S4)
      IMPLICIT NONE
      DOUBLE PRECISION C
      DOUBLE COMPLEX S1(3)
      DOUBLE COMPLEX S2(3)
      DOUBLE COMPLEX S3(3)
      DOUBLE COMPLEX S4(3)
      DOUBLE COMPLEX DENOM
      DOUBLE PRECISION M1,W1
      DOUBLE PRECISION P1(0:3)
      S4(2)=-(-S1(2)-S2(2)-S3(2))
      S4(3)=-(-S1(3)-S2(3)-S3(3))
      P1(0) =  DBLE(S1(2))
      P1(1) =  DBLE(S1(3))
      P1(2) =  DIMAG(S1(3))
      P1(3) =  DIMAG(S1(2))
      DENOM =1D0/(((M1*(-M1+(0.000000000D0, 1.000000000D0)*W1))
     $ +((P1(0)**2)-(P1(1)**2)-(P1(2)**2)-(P1(3)**2))))
      S1(1)= C*DENOM*(0.000000000D0, 1.000000000D0)*(S4(1)*S3(1)*S2(1))
      END
