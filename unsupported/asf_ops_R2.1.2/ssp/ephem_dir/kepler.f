c SccsId = @(#)kepler.f	2.41 3/24/98
      SUBROUTINE KEPLER(AM,E,EA,SE,CE)  
C/*   SUBROUTINE KEPLER(AM,E,EA,SE,CE)  -----------------
C
C  VERSION OF 1/27/87
C  PURPOSE   
C    SOLVES KEPLER EQUATION
C  INPUT 
C    AM     = MEAN ANOMALY (RAD)
C    E      = ECCENTRICITY 
C  OUTPUT
C    EA     = ECCENTRIC ANOMALY (RAD)
C    SE     = SIN(EA)  
C    CE     = COS(EA)  
C  CALL SUBROUTINES  
C    NONE  
C  REFERENCES
C    JPL EM 312/87-153, 20 APRIL 1987
C    JPL IOM 312/79.4-426, FINAL REPORT JPL CONTRACT #955140   
C    AND USER'S GUIDE FOR MULCON, 1 JAN 79.
C  ANALYSIS  
C    J. H. KWOK - JPL  
C  PROGRAMMER
C    J. H. KWOK - JPL  
C  PROGRAM MODIFICATION  
C    NONE  
C  COMMENTS  
C    THIS PROGRAM IS A STRIPED VERSION OF THE ORGINAL VERSION FOR  
C    MULCON.  THIS VERSION DOES NOT WORK FOR HYPERBOLIC CASES.  
C    THE TOLERANCES ARE SET FOR DOUBLE PRECISION. CHANGE TOL1 AND  
C    TOL2 FOR OTHER PRECISION.   SEE REFERENCE 2 FOR ALGORITHM 
C    EXPLANATION.  TOL1 SHOULD BE SET TO 1/3 OF MACHINE PRECISION.
C    TOL2 SHOULD BE SET TO 1/2 OF MACHINE PRECISION.  FOR EXAMPLE,
C    TOL1=1.D-6 AND TOL2=1.D-9 FOR 18 DIGIT MACHINE PRECISION.
C*/
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DATA TOL1,TOL2/1.D-5,1.D-8/
      DATA HALF,ZERO,ONE/0.5D0,0.D0,1.D0/
      DATA PI,TPI/3.141592653589793D0,6.283185307179586D0/  
      K=0   
   10 CONTINUE  
      ABM=DABS(AM)  
      IF (ABM.LE.PI) GO TO 20   
      IF (AM.GT.PI) AM=AM-TPI   
      IF (AM.LT.-PI) AM=AM+TPI  
      GO TO 10  
   20 SM=DSIN(ABM)  
      EA=ABM+E*SM/(ONE-DSIN(ABM+E)+SM) 
   30 SE=DSIN(EA)   
      CE=DCOS(EA)   
      K=K+1 
      IF (K.GT.10) STOP 
   50 ESE=E*SE  
      F=ABM+ESE-EA  
      D=ONE-E*CE   
      C=F/(D+HALF*F*ESE/D) 
      EA=EA+C   
      IF (DABS(C).GT.TOL1) GO TO 30 
      IF (DABS(C).LT.TOL2) GO TO 40 
      A=ONE-HALF*C*C  
      SEN=A*SE+C*CE 
      CEN=A*CE-C*SE 
      SE=SEN
      CE=CEN
      GO TO 50  
   40 SEN=SE+C*CE   
      CEN=CE-C*SE   
      SE=SEN
      CE=CEN
      IF (AM.GT.ZERO) RETURN
      SE=-SE
      EA=-EA
      RETURN
      END
