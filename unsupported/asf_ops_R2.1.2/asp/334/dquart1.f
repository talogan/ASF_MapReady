C Alaska SAR Processor (ASP) %W% %E% %U%
      SUBROUTINE DQUART1 (*,*,COEFF,N12,N34,XR,XI)
C/*   SUBROUTINE DQUART1 (*,*,COEFF,N12,N34,XR,XI) ------------
C
C              ***** DOUBLE PRECISION SUBROUTINE *****
C     SOLVES A QUARTIC EQUATION.
C     COEFFICIENTS ARE IN COEFF (COEFF(1)=COEFFICIENT OF X**4)
C     REAL PARTS OF SOLUTIONS IN XR
C     IMAGINARY PARTS OF SOLUTIONS IN XI
C     XR(1).LE.XR(2)
C     XR(3).LE.XR(4)
C     XI(1),XI(3) .LE. 0.D0
C     XI(2),XI(4) .GE. 0.D0
C     N12=DESCRIPTOR OF ROOTS 1,2 (SEE SUBROUTINE QUAD)
C     N34=DESCRIPTOR OF ROOTS 3,4 (SEE SUBROUTINE QUAD)
C     RETURN 1   DEGREE.LT.4 
C     RETURN 2   ERROR IN RELATED CUBIC
C*/
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CB(4),ZR(3),ZI(3),COEFF(5),XR(4),XI(4)
      IF (COEFF(1).EQ.0.D0) RETURN1   
      A3=COEFF(2)/COEFF(1)  
      A2=COEFF(3)/COEFF(1)  
      A1=COEFF(4)/COEFF(1)  
      A0=COEFF(5)/COEFF(1)  
      CB(1)=1.D0
      CB(2)=-A2
      CB(3)=A1*A3-4.D0*A0
      CB(4)=-A1**2-A0*A3**2+4.D0*A0*A2
      CALL DCUBE2 (*9,CB,ZR,ZI)
    9 HAFA3=.5D0*A3
      BC=-A0
      DO 1 I=1,3
      IF (ZI(I).NE.0.D0) GO TO 1
      BB=ZR(I)-A2
      CALL DPRECIS1 (*1,B1,HAFA3,BB,1)
      HAFZR=.5D0*ZR(I)
      CALL DPRECIS1 (*1,C1,HAFZR,BC,1)
      CALL DPRECIS1 (*1,B2,HAFA3,BB,2)
      CALL DPRECIS1 (*1,C2,HAFZR,BC,2)
      D=B1*C2+B2*C1
      DD=B1*C1+B2*C2
      IF (DABS(D-A1).LE.DABS(DD-A1)) GO TO 2
      D=C1
      C1=C2
      C2=D
      GO TO 2
    1 CONTINUE
      RETURN 2
    2 CALL DQUAD (1.D0,B1,C1,N12,XR(1),XI(1),XR(2),XI(2))
      CALL DQUAD (1.D0,B2,C2,N34,XR(3),XI(3),XR(4),XI(4))
      RETURN
      END   
