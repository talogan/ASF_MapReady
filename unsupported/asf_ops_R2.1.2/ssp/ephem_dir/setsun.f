c SccsId = @(#)setsun.f	2.41 3/24/98
      SUBROUTINE SETSUN(T,ES)
C/*   SUBROUTINE SETSUN(T,ES) -----------------------------------
C
C  VERSION 1/31/87
C  PURPOSE
C    SETS UP SOME SOLAR ORBITAL ELEMENTS OF THE BUILT-IN LUNI-SOLAR
C    EPHEMERIDES
C  INPUT ARGUMENTS
C    T      = CURRENT JULIAN TIME (SEC)
C  OUTPUT ARGUMENTS
C    ES     = MEAN ORBITAL ELEMENTS OF THE SUN IN EARTH MEAN EQUATOR
C             AND EQUINOX OF DATE
C      (1)  = A, SEMI-MAJOR AXIS (KM)
C      (2)  = E, ECCENTRICITY
C      (3)  = I, INCLINATION (RAD)
C      (4)  = CAPW, LONGITUDE OF ASCENDING NODE, =0 BY DEFINITION
C      (5)  = W, ARGUMENT OF PERIAPSIS (RAD)
C  CALL SUBROUTINES
C    NONE
C  REFERENCES
C    JPL EM 312/86-153, 20 APRIL 1987
C    EXPLANATORY SUPPLEMENT TO THE ASTRONOMICAL EPHEMERIS AND THE
C      AMERICAN EPHEMERIS AND NAUTICAL ALMANAC
C  ANALYSIS
C    JOHNNY H. KWOK - JPL
C  PROGRAMMER
C    JOHNNY H. KWOK - JPL
C  MODIFICATIONS
C    NONE
C  COMMENTS
C    NONE
C*/
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ES(7)
      DATA RJD/2.41502D6/
      DATA DTS/8.64D4/
      DT=T/DTS-RJD
      DT1=DT*1.D-4
      DT2=DT1*DT1
      DT3=DT1*DT2
      ES(1)=1.496D8
      ES(2)=1.675104D-2-1.1444D-5*DT1-9.4D-9*DT2
      ES(3)=.4093197474D0-6.217910D-5*DT1-2.1468D-9*DT2+1.7977D-10*DT3
      ES(4)=0.D0
      ES(5)=4.908229653D0+8.2149855D-7*DT+5.9167D-7*DT2+1.22D-9*DT3
      RETURN
      END
