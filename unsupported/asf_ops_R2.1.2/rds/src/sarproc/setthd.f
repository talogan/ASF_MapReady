      SUBROUTINE SETTHD(X,Y)  
C/*   SUBROUTINE SETTHD(X,Y)  ------------------------------
C
C  VERSION OF 1/27/85
C  PURPOSE
C    COMPUTES ROTATION MATRIX FOR CONVERSION FROM ORBITAL  
C    ELEMENTS TO CARTESIAN COORDINATES 
C  INPUT 
C    X      = 7-D ORBITAL ELEMENT SET OF THIRD BODY  
C     (1)   = SEMI-MAJOR AXIS (KM)
C     (2)   = ECCENTRICITY
C     (3)   = INCLINATION (RAD)
C     (4)   = LONGGITUDE OF ASCENDING NODE (RAD)
C     (5)   = ARGUMENT OF PERIAPSIS (RAD)
C     (6)   = MEAN ANOMALY OF THE THIRD BODY AT TR (RAD)
C     (7)   = MEAN MOTION (RAD/SEC)
C    Y      = 7-D ELEMENTS TO TRANFORM ORBITAL ELEMENTS TO CART.
C             COORD
C      (1)  = L1, COS(NODE)*COS(W) - SIN(NODE)*SIN(W)*COS(I)  
C      (2)  = M1, SIN(NODE)*COS(W) + COS(NODE)*SIN(W)*COS(I)  
C      (3)  = N1, SIN(W)*SIN(I)   
C      (4)  = L2, -COS(NODE)*SIN(W) - SIN(NODE)*COS(W)*COS(I) 
C      (5)  = M2, -SIN(NODE)*SIN(W) + COS(NODE)*COS(W)*COS(I) 
C      (6)  = N2, COS(W)*SIN(I)   
C      (7)  = SEMI-MINOR AXIS, A*SQRT(1-E*E)  
C  REFENENCES
C    JPL EM 312/87-153, 20 APRIL 1987
C    ORBITAL MOTION, A. E. ROY, 1978, P.102
C  ANALYSIS
C    J. H. KWOK - JPL  
C  PROGRAMMER
C    J. H. KWOK - JPL  
C  PROGRAM MODIFICATION  
C    NONE  
C  COMMENTS  
C    NONE  
C*/
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DIMENSION X(7),Y(7)   
      character*80 sccsid
      data sccsid /'@(#)setthd.f	1.3 96/04/09 22:51:57\0'/
      CI=DCOS(X(3)) 
      SI=DSIN(X(3)) 
      CC=DCOS(X(4)) 
      SC=DSIN(X(4)) 
      CW=DCOS(X(5)) 
      SW=DSIN(X(5)) 
      Y(1)=CC*CW-SC*SW*CI   
      Y(2)=SC*CW+CC*SW*CI   
      Y(3)=SW*SI
      Y(4)=-CC*SW-SC*CW*CI  
      Y(5)=-SC*SW+CC*CW*CI  
      Y(6)=CW*SI
      Y(7)=DSQRT(1.D0-X(2)**2)*X(1) 
      RETURN
      END   
