C--  Copyright (c)1996, California Institute of Technology.
C--  U.S. Government Sponsorship acknowledged.
C-- ==========================================================================
C--
C--  Fortran Filename:	kozsak.for
C--
C--  Description:	
C--	
C--  Notes:
C--
C-- ==========================================================================
      SUBROUTINE KOZSAK(IFLAG,GE,RE,AJ2,X,XN)

      character*100 SccsFileID
     -/'@(#)kozsak.for	5.1 98/01/08 APS/ASF\0'/

C  VERSION 4/28/88
C  PURPOSE
C    CONVERTS BETWEEN MEAN AND OSCULATING ORBITAL ELEMENTS
C  INPUT
C    IFLAG   = 1, MEAN TO OSCULATING
C            = 2, OSCULATING TO MEAN
C    GE      = PRODUCT OF GRAVITATIONAL CONSTANT * MASS OF
C              PLANET (KM**3/SEC**2)
C    RE      = RADIUS OF PLANET (KM)
C    AJ2     = J2 = -C20
C    X       = AN ARRAY OF 6 ORBITAL ELEMENTS, A, E, I, NODE, W, AND M
C              (KM, RAD)
C  OUTPUT ARGUMENTS
C    XN      = AN ARRAY OF 6 ORBITAL ELEMENTS AFTER CONVERSION
C  CALL SUBROUTINES
C    DELM
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
C    USES A COMBINATION OF KOZAI'S AND IZSAK'S THEORY.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(6),XN(6),XI(6),DX(6)
      DATA TPI/6.283185307179586D0/
      DATA ZERO/0.D0/
      IMAX=5
      IF (X(2).LT.1.D-1) IMAX=10
      IF (X(2).LT.1.D-2) IMAX=20
      IF (X(2).LT.1.D-3) IMAX=30
      DO 40 I=1,6
   40 XI(I)=X(I)
      ICOUNT=0
  100 CONTINUE
      CALL DELM(RE,GE,AJ2,X,DX)
      IF (IFLAG.EQ.2) THEN
      ICOUNT=ICOUNT+1
      DO 20 I=1,6
   20 XN(I)=XI(I)-DX(I)
      IF (XN(2).LT.ZERO) XN(2)=ZERO
      DO 30 I=3,6
   30 XN(I)=DMOD(XN(I)+TPI,TPI)
      IF (ICOUNT.EQ.IMAX) GO TO 900
      DO 50 I=1,6
   50 X(I)=XN(I)
      GO TO 100
      ELSE
      DO 60 I=1,6
   60 XN(I)=XI(I)+DX(I)
      DO 61 I=3,6
   61 XN(I)=DMOD(XN(I)+TPI,TPI)
      ENDIF
  900 CONTINUE
      RETURN
      END
