C--  Copyright (c)1996, California Institute of Technology.
C--  U.S. Government Sponsorship acknowledged.
C-- ==========================================================================
C--
C--  Fortran Filename:	delm.for
C--
C--  Description:	
C--	
C--  Notes:
C--
C-- ==========================================================================
      SUBROUTINE DELM(RE,GE,AJ2,X,DX)

      character*100 SccsFileID
     -/'@(#)delm.for	5.1 98/01/08 APS/ASF\0'/

C  VERSION 4/17/87
C    COMPUTES OSCULATING - MEAN ELEMENTS WHEN GIVEN MEAN ELEMENTS
C  INPUT
C    RE      = RADIUS OF PLANET (KM)
C    GE      = GRAVITATIONAL CONSTANT * MASS OF PLANET (KM**3/SEC**2)
C    AJ2     = J2, = -C20
C    X       = MEAN ELEMENTS
C    X(1)    = A, SEMI-MAJOR AXIS (KM)
C    X(2)    = E, ECCENTRICITY
C    X(3)    = I, INCLINATION (RAD)
C    X(4)    = NODE, LONGITUDE OF ASCENDING NODE (RAD)
C    X(5)    = W, ARGUMENT OF PERIAPSIS (RAD)
C    X(6)    = M, MEAN ANOMALY (RAD)
C  OUTPUT
C    DX      = ARRAY OF 6 VALUES = OSCULATING - MEAN ELEMENTS
C  CALL SUBROUTINES
C    ANOMLY
C  REFERENCES
C    JPL IOM 312/85.2-927, 23 JANUARY 1985, BY C. UPHOFF
C    JPL EM 312/87-153, 20 APRIL 1987
C  ANALYSIS
C    JOHNNY H. KWOK - JPL
C  PROGRAMMER
C    JOHNNY H. KWOK - JPL
C  MODIFICATIONS
C    NONE
C  COMMENTS
C    THIS PROGRAM USES AN ALGORITHM DERIVED BY C. UPHOFF (REF) WHICH
C    USES A COMBINATION OF KOZAI'S AND IZSAK'S THEORY
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(6),DX(6)
      DATA HALF,ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN
     1  /0.5D0,1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0/
      DATA F14,F34,F32,F18/0.25D0,0.75D0,1.5D0,0.125D0/
      DATA TPI/6.283185307179586D0/
      A=X(1)
      E=X(2)
      AI=X(3)
      AN=X(4)
      W=X(5)
      AM=X(6)
      CALL ANOMLY(2,AM,E,F)
C
C *** AM AND F HAVE TO BE IN THE SAME QUADRANT
C
      F=DMOD(F+TPI,TPI)
      AM=DMOD(AM+TPI,TPI)
      SI=DSIN(AI)
      CI=DCOS(AI)
      TI=SI/CI
      SI2=SI*SI
      CI2=CI*CI
      SF=DSIN(F)
      CF=DCOS(F)
      S2F=DSIN(TWO*F)
      U=F+W
      E2=E*E
      ESF=E*SF
      D1=ONE-E2
      D2=DSQRT(D1)
      D3=E*CF
      D4=ONE+D3
      D42=D4*D4
      D5=ONE+D2
      D6=(THREE*CI2-ONE)/D5
      P=A*D1
      D7=DSQRT(GE/P)
      R=P/D4
      RDOT=D7*ESF
      TWOU=TWO*U
      TWOW=TWO*W
      S2U=DSIN(TWOU)
      C2U=DCOS(TWOU)
      SF2W=DSIN(F+TWOW)
      D8=THREE*F+TWOW
      S3F2W=DSIN(D8)
      CF2W=DCOS(F+TWOW)
      C3F2W=DCOS(D8)
      Q1=AJ2*(RE/P)**2
      DI=F34*Q1*SI*CI*(C2U+E*CF2W+E/THREE*C3F2W)
      DP=TWO*P*TI*DI
      DUMMY1=F-AM+ESF-HALF*S2U-HALF*E*SF2W-E*S3F2W/SIX
      DN=-F32*Q1*CI*DUMMY1
      DR=-F14*P*Q1*((THREE*CI2-ONE)*(TWO*D2/D4+D3/D5+ONE)-SI2*C2U)
      DRDOT=F14*D7*Q1*(D6*ESF*(D2*D5+D42)-TWO*SI2*D42*S2U)
      DU=-F18*Q1*(SIX*(ONE-FIVE*CI2)*(F-AM)+FOUR*ESF*((ONE-SIX*CI2)
     1  -D6)-D6*E2*S2F+TWO*(FIVE*CI2-TWO)*E*SF2W+(SEVEN*CI2-ONE)*S2U
     2  +TWO*CI2*E*S3F2W)
      PNW=P+DP
      AINW=AI+DI
      ANNW=AN+DN
      RNW=R+DR
      RDOTNW=RDOT+DRDOT
      UNW=U+DU
      AA=PNW/RNW-ONE
      BB=DSQRT(PNW/GE)*RDOTNW
      ENW2=AA*AA+BB*BB
      ENW=DSQRT(ENW2)
      FNW=DATAN2(BB,AA)
      ANW=PNW/(ONE-ENW2)
      WNW=UNW-FNW
      CALL ANOMLY(5,FNW,ENW,AMNW)
      DX(1)=ANW-A
      DX(2)=ENW-E
      DX(3)=AINW-AI
      DX(4)=ANNW-AN
      DX(5)=WNW-W
      DX(6)=AMNW-AM
      RETURN
      END
