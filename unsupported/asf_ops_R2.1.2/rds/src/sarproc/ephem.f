      SUBROUTINE EPHEM(ISUN,IMOON,ISRP,T,ES,EM)
C/*   SUBROUTINE EPHEM(ISUN,IMOON,ISRP,T,ES,EM)  -------------
C
C  VERSION OF 1/31/87
C  PURPOSE
C    COMPUTES THE ORBITAL ELEMENTS OF A SUN AND/OR A MOON
C  INPUT ARGUMENTS
C    ISUN   = SEE SUBROUTINE ASAP
C    IMOON  = SEE SUBROUTINE ASAP
C    ISRP   = SEE SUBROUTINE ASAP
C    T      = CURRENT JULIAN TIME (SEC)
C    ES(3)  = INCLINATION OF THE SUN (RAD)
C  OUTPUT ARGUMENTS
C    ES(6)  = MEAN ANOMALY OF THE SUN AT T (RAD)
C    ES(7)  = MEAN MOTION OF THE SUN (RAD/SEC)
C    EM     = SEVEN ORBITAL ELEMENTS OF THE MOON AT T IN EARTH MEAN
C             EQUATOR AND EQUINOX OF EPOCH
C      (1)  = SEMI-MAJOR AXIS (KM)
C      (2)  = ECCENTRICITY
C      (3)  = INCLINATION (RAD)
C      (4)  = LONGITUDE OF ASCENDING NODE (RAD)
C      (5)  = ARGUMENT OF PERIAPSIS (RAD)
C      (6)  = MEAN ANOMALY AT T (RAD)
C      (7)  = MEAN MOTION (RAD/SEC)
C  CALL SUBROUTINES
C    NONE
C  REFERENCES
C    JPL EM 312/87-153, 20 APRIL 1987
C    EXPLANATORY SUPPLEMENT TO THE ASTRONOMICAL EPHEMERIS AND THE
C      AMERICAN EPHEMERIS AND NAUTICAL ALMANAC
C  ANALYSIS
C    JOHNNY H. KWOK - JPL
C  PROGRAMMER
C    JOHNNY H. KWOK - JPL
C  MODIFICATIONS
C    NONE
C  COMMENTS
C    FOR MERCURY, VENUS, AND MARS, THERE IS NO LUNI-PERTURBATION, USER
C    SHOULD INPUT THE SUN'S ORBITAL ELEMENTS RELATIVE TO PLANET EQUATOR
C    OF EPOCH.  THE INERTIAL X-AXIS COULD BE PLANET EQUINOX OR THE
C    DIRECTION OF THE PRIME MERIDIAN AT SOME EPOCH SUCH AS DEFINED BY
C    THE IAU.  FOR EARTH, THE ORBITAL ELEMENTS OF THE MOON VARIES TOO
C    MUCH TO USE TWO-BODY EPHEMERIS.  USER MAY USE THE MEAN LUNAR
C    EPHEMERIS BUILT INTO THE PROGRAM.  THE SHORT PERIOD VARIATIONS OF
C    THE MOON MAY INTRODUCE AN ERROR UP TO TWO DEGREES.  AS FOR THE SUN,
C    THE MEAN ORBITAL ELEMENTS OF THE SUN RELATIVE TO EARTH EQUATOR AND
C    EQUINOX ARE BUILT IN ALSO.  BUT SINCE THEY DON'T CHANGE AS RAPIDLY
C    AS THOSE OF THE MOON, THEY ARE INITIALIZE ELSEWHERE AND THEN HELD
C    FIXED FOR THE DURATION OF THE TRAJECTORY PROPAGATION
C*/
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ES(7),EM(7)
      DATA RJD/2.41502D6/
      DATA PI,TPI/3.141592653589793D0,6.283185307179586D0/
      DATA ANS,ANM/1.720196977D-2,.22997150294101D0/
      DATA DTS/8.64D4/
      character*80 sccsid
      data sccsid /'@(#)ephem.f	1.3 96/04/09 22:51:47\0'/
      DT=T/DTS-RJD
      DT1=DT*1.D-4
      DT2=DT1*DT1
      DT3=DT1*DT2
C
C  SUN RELATIVE TO EARTH IN MEAN EQUATOR OF DATE
C
      IF (ISUN.EQ.1.OR.ISRP.EQ.1) THEN
      ES(6)=DMOD(6.256583575D0+ANS*DT
     1  -1.9548D-7*DT2-1.22D-9*DT3,TPI)
      ES(7)=ANS/DTS
      ENDIF
C
C  MOON IN EARTH MEAN EQUINOX OF DATE
C
      IF (IMOON.EQ.1) THEN
      EM(1)=3.844D5
      EM(2)=5.4900489D-2
      EM(3)=8.980410805D-2
      EM(4)=4.523601515D0-9.242202942D-4*DT+2.71748D-6*DT2+8.73D-10*DT3
      EM(5)=5.835151539D0+1.944368001D-3*DT-1.35071D-5*DT2-4.538D-9*DT3
     1  -EM(4)
      EM(6)=4.719966572D0+ANM*DT-1.4835D-6*DT2
     1  +6.80678D-10*DT3-EM(4)-EM(5)
      DO 130 I=4,6
  130 EM(I)=DMOD(EM(I),TPI)
      EM(7)=ANM/DTS
C
C  CONVERT TO EARTH MEAN EQUATOR OF DATE
C
      SE=DSIN(ES(3))
      CE=DCOS(ES(3))
      SI=DSIN(EM(3))
      CI=DCOS(EM(3))
      SO=DSIN(EM(4))
      CO=DCOS(EM(4))
      CA=SI*SE*CO-CE*CI
      ALPHA=DACOS(CA)
      EM(3)=PI-ALPHA
      SA=DSIN(ALPHA)
      SW=SO*SI/SA
      SD=SO*SE/SA
      SSO=DSIGN(1.D0,SO)
      SWSO=SW*SO*SSO
      CWSO=(SW*CE*CO+SD*CI)*SSO
      EM(4)=DATAN2(SWSO,CWSO)
      SDSO=SD*SO*SSO
      CDSO=(SD*CO*CI+SW*CE)*SSO
      DEL=DATAN2(SDSO,CDSO)
      EM(5)=EM(5)+DEL
      ENDIF
      RETURN
      END
