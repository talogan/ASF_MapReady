c SccsId = @(#)COMMON_ZERO_ONE_CEOS.f	2.41 3/24/98
        SUBROUTINE CROSS(X1,Y1,Z1,X2,Y2,Z2,U1,U2,U3)
        IMPLICIT REAL*8 (A-H,O-Z)

        U1=Y1*Z2-Z1*Y2
        U2=Z1*X2-X1*Z2
        U3=X1*Y2-Y1*X2
        R=(U1**2+U2**2+U3**2)**.5
        U1=U1/R
        U2=U2/R
        U3=U3/R
        RETURN
        END
        subroutine EME_to_EBF(omega0_e,second0,second,x,y,z,x1,y1,z1)
        implicit real*8 (a-h,o-z)
        include 'ssp2_const.inc'
cc      parameter (omega_e = 360.9856296d0 / (24.0d0 * 3600.0d0 ))

        omega_e_day = omega0_e+omega_e*(second-second0)*180.d0/pi
        x1=cosd(omega_e_day)*x+sind(omega_e_day)*y
        y1=-1*sind(omega_e_day)*x+cosd(omega_e_day)*y
        z1 = z

        return
        end
	function dis(x1,x2,x3,y1,y2,y3)
	real*8 x1,x2,x3,y1,y2,y3,dis

	dis = sqrt((x1-y1)**2+(x2-y2)**2+(x3-y3)**2)

	return
	end
	function fd(lambda,x1,x2,x3,v1,v2,v3,t1,t2,t3,w1,w2,w3)
        real*8 lambda,x1,x2,x3,v1,v2,v3,t1,t2,t3,w1,w2,w3,fd
	real*8 vdotr,rs

	vdotr = (v1-w1)*(x1-t1)+(v2-w2)*(x2-t2)+(v3-w3)*(x3-t3)
	rs = sqrt((x1-t1)**2+(x2-t2)**2+(x3-t3)**2)
	fd = -2*vdotr/(lambda*rs)

	return
	end
c****************************************************************
c       Jet Propulsion Laboratory
c       Project: ASF RadarSAT
c       Task:    ScanSAR Processor
c       Code:    Fortran 77
c------------------------------------------------------------------
c
c       8/31/94 A.Chu Initial code & utilize the modified Johnny Kwok's 
c	              ASAP code.
c------------------------------------------------------------------
	subroutine propagation(r_x1,r_y1,r_z1,v_x1,v_y1,v_z1,
     *                         start_time,delta_t,
     *                         r_x2,r_y2,r_z2,v_x2,v_y2,v_z2)
	implicit none

        character*128 SccsId_COMMON_ZERO_ONE_CEOS
        data SccsId_COMMON_ZERO_ONE_CEOS
     +  /'@(#)PPCOMMON_ZERO_ONE_CEOS.f:2.41'/


c       ----------------
c       Abstract:
c               Propagate one state vector based on starting state vector,
c		start time in second(convert from GMT) and time interval.
c
c               Input:
c		     1).State Vector:
c		        S/C Position & velocity.
c                    2).GMT Time in second:
c	                Start time & time interval
c               Output:
c	             State Vector of the start time + time interval.
c
c               Subroutines included:
c                    propag_init.f
c                    eme_kepler.f
c                    asap.f (modified version by A.Chu)
c	-----------------
c	INPUT PARAMETERS PASSING
c	-----------------
	real*8	r_x1,r_y1,r_z1	!Position of S/C
	real*8	v_x1,v_y1,v_z1	!Velocity of S/C
	real*8	start_time	!Start_time in second
	real*8	delta_t		!Time interval in sec.

c	-----------------
c	OUTPUT PARAMETERS PASSING
c	-----------------
	real*8	r_x2,r_y2,r_z2	!New Position of S/C of delta_t
	real*8	v_x2,v_y2,v_z2	!New velocity of S/C of delta_t

c	-----------------
c	OUTPUT PARAMETERS PASSING
c	-----------------
	real*8	kepler(6)	!Kepler 6 elements
	real*8	state_vec(6)	!New state vectors

c	------------	
c	1. Setup flags & coefficients
c	2. Convert State Vector to Kepler Elements
c	------------	
	call propag_init(r_x1,r_y1,r_z1,v_x1,v_y1,v_z1,kepler)

c	   write(6,*)'kepler r',kepler(1),kepler(2),kepler(3)
c	   write(6,*)'kepler v',kepler(4),kepler(5),kepler(6)

c	------------	
c	State Vector propagation
c	------------	
	call asap(kepler,start_time,delta_t,state_vec)

	r_x2=state_vec(1)
	r_y2=state_vec(2)
	r_z2=state_vec(3)
	v_x2=state_vec(4)
	v_y2=state_vec(5)
	v_z2=state_vec(6)

D	write(6,*)' state R:',state_vec(1),state_vec(2),state_vec(3)
D	write(6,*)' state  V:',state_vec(4),state_vec(5),state_vec(6)
	end

c****************************************************************
c       Jet Propulsion Laboratory
c       Project: ASF RadarSAT
c       Task:    ScanSAR
c       Code:    Fortran 77
c------------------------------------------------------------------
c
c       8/31/94 A.Chu Modified for ScanSAR
c------------------------------------------------------------------
	subroutine propag_init(r_x,r_y,r_z,v_x,v_y,v_z,kepler)

c     *                         GE,RE,RATE,PM,ELLIP,RATM)
c	---------------------
c	Abstract:
c	        1. Convert state vector to Kepler coordinate
c	        2. Set up flags & coefficients.
c	---------------------
c	
C    ALL THE UNIT CONVERSIONS SHOULD BE MADE IN THIS DRIVER
C    PROGRAM SO THAT ALL SUBROUTINES WOULD BE CONSISTENT IN UNITS. 
C   
      	IMPLICIT DOUBLE PRECISION (A-H,O-Z)   

	real*8 r_x,r_y,r_z,v_x,v_y,v_z
	real*8 kepler(6)
	integer ii,jj

      	DIMENSION ORB(6),Y(6),X(6)
      	DIMENSION TINT(2),TFIN(2),TREF(2)
      	COMMON/OPTION/L,M,IRES,ISUN,IMOON,IEPHEM,IDRAG,IDENS,ISRP,IORB
     1  ,IPRINT,INODE,IPLOT
      	COMMON/ATIME/TI,TF,TR
      	COMMON/PLTCON/GE,RE,RATE,PM,AJ2,ELLIP,RATM
      	COMMON/ATMCON/RDENS,RHT,SHT,ALTMAX,WT
      	COMMON/SPCCON/AREAD,AREAS,SCMASS,CDRAG,CSRP
      	COMMON/SUNCON/GS,ES(7),ET(7)
      	COMMON/MUNCON/GM,EM(7),EN(7)
      	COMMON/HARMON/C(41,41),S(41,41)
      	COMMON/accurcy_err/relerr,abserr
      	EXTERNAL DER
      	DATA NEQ/6/
      	DATA DTR,DTS/.1745329251994330D-1,8.64D4/
      	DATA ZERO/0.D0/
      	DATA HSTART,HLARGE/60.D0,1.D99/

	call eme_kepler (r_x,r_y,r_z,v_x,v_y,v_z,kepler)

	L=20	!Mean value, (2)
	M=20	!Mean value, (0)
	IRES=0	!1=Resonance effect in Tesserals,not affect in zonal term
	ISUN=0  !No Solar Gravity
	IMOON=0 !No Lunar gravity
	IEPHEM=1!Earth orbiting spacecraft using earth mean equator
	IDRAG=0 !No Air Drag
	IDENS=1	!Use static 1977 Earth Model
	ISRP=0	!No Solar Radiation Pressure
	IORB=0	!Orbit element are osculating value
	IPRINT=0!Print at constant step as specified by step,1=periapsis
	INODE=0 !No Print nodal crossing info
	IPLOT=0 !No plotting
	RELERR=1.d-12	!Relative accuracy of the integration
	ABSERR=1.d-12	!Absolute accuracy of the integration
        GE = 3.9860045D5 !EARTH gmass KM**3/SEC**2 (JPL DE118)
cc	GE=GE/1.d9	!KM**3/SEC**2
        RE= 6378.140D0 !FOR EARTH(KM)
        RATE= 4.178074216D-3 !FOR EARTH (deg/sec)
	PM=99.652865509d0 !Prime meridian relative to the inertial
			  !x-axis at tref (deg)
        ELLIP = .8182D-1  !Ellipticity of the reference ellipsoid to
			  !compute geodetic altitude
        RATM= RE+90 	  !Atmospheric blckage (KM)
        RDENS=0.d0
        RHT=0.d0
        SHT=0.d0
        ALTMAX=1.d0	  !Maximum altitude for Drag perturbation
        WT=1.d0		  !Weight factor to be applied the density
        AREAD=10.d-6	  !s/c area for Drag (km**2)
        AREAS=10.d-6	  !s/c area for solar wind (km**2)
        SCMASS=2000.d0    !Effective S/C mass for the length
			  !propagation (KG)
        CDRAG=2.d0	  !Drag coefficient
        CSRP=6.6d-3       !SRP Constant (kg/km-sec**2)
        GS=.13271244d12	  !Sun gmass
	ES(1)=0.d0        !Sun Ephemeris
	ES(2)=0.d0
	ES(3)=0.d0
	ES(4)=0.d0
	ES(5)=0.d0
	ES(6)=0.d0
	ES(7)=0.d0
        GM=.490279d4	  !Moon gmass
        EM(1)=0.d0
        EM(2)=0.d0
        EM(3)=0.d0
        EM(4)=0.d0
        EM(5)=0.d0
        EM(6)=0.d0
        EM(7)=0.d0
	do ii=1,20
	   do jj=1,20
	   c(ii,jj)=0.d0
	   s(ii,jj)=0.d0
	   enddo
	enddo
	C(3,1)=-.10826271d-2
	S(3,1)=-.10826271d-2
C

	return
	end
         SUBROUTINE ASAP(orb,tints,delta_t,x)
ccc         SUBROUTINE ASAP(orb,tints,delta_t,trefs,x)
ccc         SUBROUTINE ASAP(orb,tints,tfins,trefs,x)
ccc         SUBROUTINE ASAP(orb,tint,tfin,tref,x)

c	----------------
c	Abstract:
c	        Propagate one state vector based on one Kepler coordinate
c	        point and GMT Time Interval.
c	        
c               Input: 
c	            1).ORB(6): [Kepler 6 elements] 
c	               a,ecc,inc,ascending node,argument periapsis &
c		       mean anomly.
c	            2).GMT Time :
c	               Tint: starting time
c	               Tfin: ending time
c	               Tref: reference time
c               Output: 
c	             x(6): [state vectors]
c	             r_x,r_y,r_z, & v_x,v_y,v_z
c	----------------
C
C/*	SUBROUTINE ASAP  --------------------------
C
C  VERSION OF 2/29/88
C  PURPOSE   
C    MAIN DRIVER FOR THE ARTIFICIAL SATELLITE ANALYSIS PROGRAM (ASAP).
C    ASAP IS AN EPHEMERIS PROPAGATION PROGRAM FOR ORBITING PLANETARY
C    SPACECRAFTS.  IT USES COWELL'S METHOD OF SPECIAL PERTURBATION.  IT
C    INCLUDES HARMONICS OF UP TO 40X40 FIELD, LUNI-SOLAR GRAVITY, DRAG,
C    AND SOLAR RADIATION PRESSURE.  IT USES A RUNGE-KUTTA 8TH ORDER
C    METHOD WITH STEP SIZE CONTROL FOR NUMERICAL INTEGRATION.  THE
C    PROGRAM IS IN MODULAR FORM SO THAT USERS CAN MODIFY IT TO INCLUDE
C    DIFFERENT I/O MODULES, OR DENSITY MODELS, OR DIFFERENT INTEGRATORS,
C    ETC.  IT ASSUMES A PLANET MEAN  EQUATOR OF EPOCH SYSTEM AND IGNORES
C    POLAR MOTION.  ALL INPUTS ARE EITHER I30 FOR INTEGERS OR D30. FOR
C    DOUBLE PRECISION WITH ONE VALUE TO EACH RECORD.  THE VALUE MUST BE
C    PLACED WITHIN THE FIRST 30 COLUMNS BUT THERE IS NO NEED TO RIGHT
C    JUSTIFY.  COLUMNS 31 TO 80 CAN BE USED FOR COMMENTS.  THE EXCEPTION
C    IS THE INPUT OF THE SPHERICAL HARMONIC COEFFICIENTS.  THIS PROGRAM
C    AND ITS DOCUMENTATION ARE AVAILABLE FROM COSMIC, CAT. #NPO-16731.
C    THERE IS A COMPANION PROGRAM CALLED LONG-TERM PREDICTION (LOP) THAT
C    USES AN AVERAGING METHOD.  LOP IS AVAILABLE FROM COSMIC, CAT.
C    #NPO-17052.
C  INPUT
C    L      = DEGREE OF GRAVITY HARMONICS TO BE INCLUDED
C    M      = ORDER OF GRAVITY HARMONICS TO BE INCLUDED.  MAXIMUM FOR L
C             AND M IS 40.  L CAN BE BIGGER THAN M.
C    IRES   = AN OPTION TO ONLY INCLUDE TESSERALS OF ORDER IRES.  CAN BE
C             USED TO STUDY RESONANCE EFFECTS WHEREBY NON-RESONANT
C             TERRESALS ARE NOT INCLUDED IN THE CALCULATION AND THUS
C             CUTTING PROPAGATION TIME.  DOES NOT AFFECT ZONAL TERMS.
C           = 0, TO TURN OFF THE OPTION.
C    ISUN   = 0, NO SOLAR GRAVITY
C           = 1, WITH SOLAR GRAVITY
C    IMOON  = 0, NO LUNAR GRAVITY
C           = 1, WITH LUNAR GRAVITY
C    IEPHEM = 0, USE TWO-BODY EPHEMERIDES FOR SUN AND MOON
C           = 1, FOR EARTH ORBITING SPACECRAFT ONLY, USE BUILT-IN
C                LUNI-SOLAR EPHEMERIDES.  MUST USE EARTH MEAN EQUATOR
C                AND EQUINOX OF EPOCH SYSTEM
C    IDRAG  = 0, NO DRAG
C           = 1, WITH DRAG
C    IDENS  = 0, USE EXPONENTIAL MODEL
C           = 1, USE STATIC 1977 EARTH MODEL
C    ISRP   = 0, NO SOLAR RADIATION PRESSURE
C           = 1, WITH SOLAR RADIATION PRESSURE
C    IORB     FLAG FOR INPUTING EITHER MEAN OR OSCULATING ORBITAL
C             ELEMENTS.  USES THE VALUE OF C20 TO COMPUTE SHORT PERIOD
C             EFFECTS
C           = 0, INPUT ORBITAL ELEMENTS ORB ARE OSCULATING VALUES
C           = 1, INPUT ORBITAL ELEMENTS ORB ARE MEAN VALUES
C    IPRINT = 0, PRINT AT CONSTANT STEP AS SPECIFIED BY STEP
C           = 1, ALSO PRINT AT PERIAPSIS AND APOAPSIS
C    INODE  = 0, NO NODAL CROSSING PRINT
C           = 1, NODAL CROSSING TIME AND INFO PRINT
C    IPLOT  = 0, NO OUTPUT ASCII FILE FOR PLOTTING
C           = 1, WRITE OUTPUT AT EVERY 'STEP', APSIS AND NODAL CROSSINGS
C                TO FILE 8
C    ORB      OSCULATING OR MEAN ORBITAL ELEMENTS OF THE SPACECRAFT AT
C             TINT.  SEE IORB
C       (1) = A, SEMI-MAJOR AXIS (KM)
C       (2) = E, ECCENTRICITY
C       (3) = I, INCLINATION (DEG)
C       (4) = CAPW, LONGITUDE OF ASCENDING NODE (DEG)
C       (5) = W, ARGUMENT OF PERIAPSIS (DEG)
C       (6) = M, MEAN ANOMALY (DEG)
C    RELERR = RELATIVE ACCURACY OF THE INTEGRATOR.  RECOMMENDED VALUES
C             BETWEEN 1.D-6 TO 1.D-12
C    ABSERR = ABSOLUTE ACCURACY OF THE INTEGRATOR.  RECOMMENDED VALUES
C             BETWEEN 1.D-6 TO 1.D-12
C    STEP   = TIME STEP TO PRINT (SEC).
C    TINT(1)= INITIAL CALENDAR DATE OF RUN (YYYYMMDD.).  FOR EXAMPLE,
C             19880726.D0, FOR 26 JULY 1988.  ALL TIME USED IN THIS
C             PROGRAM ASSUMES EPHEMERIS TIME AND NOT UNIVERSAL TIME.
C        (2)= INITIAL TIME OF DAY OF RUN (HHMMSS.SS...D0).  FOR EXAMPLE,
C             130723.1234D0, FOR 13 HR 7 MIN 23.1234 SEC
C    TFIN   = SAME AS TINT EXCEPT FOR END OF RUN
C    TREF   = SAME AS TINT EXCEPT THIS IS THE TIME CORRESPONDING TO THE
C             POSITION OF THE LUNI-SOLAR EPHEMERIDES (ES AND EM) AND THE
C             PRIME MERIDIAN (PM) INPUT.  IF IEPHEM=1 TO USE BUILT-IN
C             LUNI-SOLAR EPHEMERIDES, THEN THIS IS THE TIME
C             CORRESPONDING TO THE PRIME MERIDIAN INPUT ONLY
C    GE     = PRODUCT OF GRAVITATIONAL CONSTANT AND MASS OF PLANET
C             (KM**3/SEC**2).  RECOMMENDED VALUES (JPL DE118),
C           = 3.9860045D5, FOR EARTH
C           = 3.2485877D5, FOR VENUS
C           = 4.2828287D4, FOR MARS
C    RE     = RADIUS OF PLANET (KM).  RECOMMENDED VALUES (IAU 1982),
C           = 6378.140D0, FOR EARTH
C           = 6051.D0, FOR VENUS
C           = 3393.4D0, FOR MARS
C    RATE   = ROTATION RATE OF THE PLANET (DEG/SEC).  RECOMMENDED VALUES
C             (IAU 1982),
C           = 4.178074216D-3, FOR EARTH
C           = -1.71460706D-5, FOR VENUS
C           = 4.061249803D-3, FOR MARS
C    PM     = LOCATION OF THE PRIME MERIDIAN RELATIVE TO THE INERTIAL
C             X-AXIS AT TREF (DEG).
C    ELLIP  = ELLIPTICITY OF THE REFERENCE ELLIPSOID.  USED BY DRAG
C             TO COMPUTE GEODETIC ALTITUDE FOR ATMOSPHERE DENSITY
C             EVALUATION AND BY OUTPUT ROUTINE TO COMPUTE GEODETIC
C             ALTITUDE.  RECOMMENDED VALUES (IAU 1982),
C           = 0.D0, TO USE A SPHERE
C           = .8182D-1, FOR EARTH
C           = 0.D0, FOR VENUS
C           = .1017D0, FOR MARS
C    RATM     RADIUS OF THE PLANET INCLUDING ATMOSPHERIC BLOCKAGE (KM).
C             IT IS USED TO COMPUTE SHADOW ENTRY AND EXIT IN SOLAR
C             RADIATION PRESSURE COMPUTATION.  RECOMMENDED VALUES,
C           = RE+90 KM FOR VENUS, EARTH, AND MARS
C           = RE FOR MERCURY OR THE MOON
C    RDENS    REFERENCE DENSITY AT REFERENCE HEIGHT RHT TO BE USED BY
C             THE EXPONENTIAL DENSITY MODEL (KG/KM**3)
C    RHT      REFERENCE HEIGHT FOR THE EXPONENTIAL DENSITY MODEL (KM).
C             USE PERIAPSIS ALTITUDE IF POSSIBLE.
C    SHT      SCALE HEIGHT OF THE EXPONENTIAL DENSITY MODEL (KM)
C    ALTMAX   MAXIMUM ALTITUDE TO INCLUDE DRAG PERTURBATION (KM)
C    WT       WEIGHT FACTOR TO BE APPLIED TO THE DENSITY, 1.D0 FOR
C             NOMINAL, 2.D0 FOR TWICE DENSER, ETC.
C    AREAD    EFFECTIVE SPACECRAFT AREA FOR DRAG (KM**2)
C    AREAS    EFFECTIVE SPACECRAFT AREA FOR SOLAR RADIATION PRESSURE
C             (KM**2)
C    SCMASS   EFFECTIVE SPACECRAFT MASS FOR THE LENGTH OF PROPOGATION
C             (KG)
C    CDRAG    DRAG COEFFICIENT, RECOMMENDED VALUES BETWEEN 2.D0 TO 2.2D0
C    CSRP     SOLAR RADIATION PRESSURE CONSTANT (KG/KM-SEC**2).
C             RECOMMENDED VALUES,
C           = G * 4.4D-3, FOR EARTH,
C           = G * 8.4D-3, FOR VENUS,
C           = G * 1.9D-3, FOR MARS,
C             WHERE G < 1 FOR TRANSLUCENT MATERIAL
C                     = 1 FOR BLACK BODY
C                     = 2 FOR PERFECTLY REFLECTIVE MATERIAL
C    GS       PRODUCT OF GRAVITATIONAL CONSTANT AND MASS OF SUN
C             (KM**3/SEC**2).  RECOMMENDED VALUE (JPL DE118),
C           = .13271244D12
C    ES     = ORBITAL ELEMENTS OF THE SUN IN PLANET EQUATOR OF EPOCH,
C             USED IN CALCULATING POINT MASS PERTURBATION DUE TO THE
C             SUN AND SOLAR RADIATION PRESSURE AS WELL.  SEE IEPHEM FOR
C             BUILT-IN LUNI-SOLAR EPHEMERIDES FOR THE EARTH,
C      (1)  = SEMI-MAJOR AXIS (KM)
C      (2)  = ECCENTRICITY
C      (3)  = INCLINATION (DEG)
C      (4)  = LONGITUDE OF THE ASCENDING NODE (DEG), EQUAL TO ZERO IF
C             X-AXIS IS EQUINOX OF EPOCH
C      (5)  = ARGUMENT OF PERIAPSIS (DEG)
C      (6)  = MEAN ANOMALY AT TREF (DEG)
C      (7)  = MEAN MOTION (DEG/SEC)
C    GM       PRODUCT OF GRAVITATIONAL CONSTANT AND MASS OF THE MOON
C             (KM**3/SEC**2).  RECOMMENDED VALUE,
C           = .490279D4, FOR EARTH'S MOON
C    EM       AN ARRAY OF 7 ORBITAL ELEMENTS OF A MOON IN PLANET EQUATOR
C             OF EPOCH SIMILAR TO ES.  SEE IEPHEM FOR BUILT-IN
C             LUNI-SOLAR EPHEMERIDES FOR THE EARTH.
C    N, M, C, S
C             DEGREE, ORDER, OF THE SPHERICAL HARMONIC COEFFICIENTS.
C             CONTINUE TO AS MANY SPHERICAL HARMONICS AS THE FIELD
C             REQUIRES, UP TO 40X40 FIELD.  N SHOULD BE WITHIN COLUMNS
C             1 TO 5, M SHOULD BE WITHIN COLUMNS 6 TO 10, CNM SHOULD BE
C             WITHIN COLUMNS 11 TO 40, AND SNM SHOULD BE WITHIN COLUMNS
C             41 TO 70.  GRAVITY COEFFICIENTS GREATER THAN L AND M ARE
C             NOT COMPUTED IN THE FORCE COMPUTATION.  ONE MAY SET UP A
C             40X40 FIELD HERE ANYWAY EVEN THOUGH A SMALLER FIELD IS
C             RUN.  THE PENALTY IS LONGER DISK READ TIME.
C  CALL SUBROUTINES  
C    JULIAN, KOZSAK, KEPLER, COORD, SETSM, SETTHD, RK78CN, POUT, RK78
C  REFERENCES
C    JPL EM 312/87-153, 20 APRIL 1987
C  ANALYSIS
C    J. H. KWOK - JPL  
C  PROGRAMMER
C    J. H. KWOK - JPL  
C  PROGRAM MODIFICATIONS 
C    NONE
C  COMMENTS
C    NOTE THAT THE SPHERICAL HARMONIC COEFFICIENTS ARE STORED
C    AS C22=C(3,3), ETC.  COEFS. ARE DIMENSIONED FOR 40X40 FIELD.   
C    INPUTS ASSUME KM, SEC, KG, AND DEGREES.
C
C    ALL THE UNIT CONVERSIONS SHOULD BE MADE IN THIS DRIVER
C    PROGRAM SO THAT ALL SUBROUTINES WOULD BE CONSISTENT IN UNITS. 
C   
C*/
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DIMENSION ORB(6),Y(6),X(6)
      DIMENSION TINT(2),TFIN(2),TREF(2)
      COMMON/OPTION/L,M,IRES,ISUN,IMOON,IEPHEM,IDRAG,IDENS,ISRP,IORB
     1  ,IPRINT,INODE,IPLOT
      COMMON/ATIME/TI,TF,TR
      COMMON/PLTCON/GE,RE,RATE,PM,AJ2,ELLIP,RATM
      COMMON/ATMCON/RDENS,RHT,SHT,ALTMAX,WT
      COMMON/SPCCON/AREAD,AREAS,SCMASS,CDRAG,CSRP
      COMMON/SUNCON/GS,ES(7),ET(7)
      COMMON/MUNCON/GM,EM(7),EN(7)
      COMMON/HARMON/C(41,41),S(41,41)
      COMMON/accurcy_err/relerr,abserr	!A.Chu
      EXTERNAL DER
      DATA NEQ/6/
      DATA DTR,DTS/.1745329251994330D-1,8.64D4/
      DATA ZERO/0.D0/
      DATA HSTART,HLARGE/60.D0,1.D99/
C
C  BEGIN READING INPUT DATA.  THIS BLOCK CAN BE REPLACED
C  IF A NAMELIST ROUTINE IS AVAILABLE
C
cAC      OPEN(5,FILE='5')
cAC      READ(5,4000)L,M,IRES,ISUN,IMOON,IEPHEM,IDRAG,IDENS,ISRP,IORB
cAc     1  ,IPRINT,INODE,IPLOT
cAC      READ(5,3000)(ORB(I),I=1,6),RELERR,ABSERR,STEP
cAC      READ(5,3000)(TINT(I),I=1,2),(TFIN(I),I=1,2),(TREF(I),I=1,2)
cAC      READ(5,3000)GE,RE,RATE,PM,ELLIP,RATM
cAC      READ(5,3000)RDENS,RHT,SHT,ALTMAX,WT
cAC      READ(5,3000)AREAD,AREAS,SCMASS,CDRAG,CSRP
cAC      READ(5,3000)GS,(ES(I),I=1,7)
cAC      READ(5,3000)GM,(EM(I),I=1,7)
	 orb(1)=orb(1)/1000.d0	!km A.Chu
C
C  BEGIN OUTPUT INPUT DATA
C
c	write(6,*)' write i/p data'
c     	WRITE((6,*)4000)L,M,IRES,ISUN,IMOON,IEPHEM,IDRAG,IDENS,ISRP,
c     1       iorb,IPRINT,INODE,IPLOT
c     	WRITE((6,*)3000)(ORB(I),I=1,6),RELERR,ABSERR,STEP
c     	WRITE((6,*)3000)TINT,TFIN,TREF
c     	WRITE((6,*)3000)GE,RE,RATE,PM,ELLIP,RATM
c     	WRITE((6,*)3000)RDENS,RHT,SHT,ALTMAX,WT
c     	WRITE((6,*)3000)AREAD,AREAS,SCMASS,CDRAG,CSRP
c    	WRITE((6,*)3000)GS,(ES(I),I=1,7)
c     	WRITE((6,*)3000)GM,(EM(I),I=1,7)

C  INPUT AND OUTPUT GRAVITY FIELD
C
      DO 10 I=1,41
      DO 10 J=1,41
      C(I,J)=ZERO
   10 S(I,J)=ZERO
      DO 20 N=1,2000
cAC      READ(5,5000,END=30)I,J,C(I+1,J+1),S(I+1,J+1)
   20 CONTINUE
   30 CONTINUE
      IF (L.EQ.0.AND.M.EQ.0) GO TO 50
      DO 40 I=2,L
      DO 40 J=0,I
      IF (J.GT.M) GO TO 40
c      WRITE((6,*)5000)I,J,C(I+1,J+1),S(I+1,J+1)
   40 CONTINUE
   50 CONTINUE
c qdn 4/3/97 modified to add the missing input
      C(3,1) = -0.10826271D-2
      AJ2=-C(3,1)
      IF (IRES.EQ.0) IRES=1
C
C  CONVERT CALENDAR DATE AND TIME TO JULIAN DATE
C
cAC	write(6,*)'gmt:',tint,tfin,tref
ccc      CALL JULIAN(TINT,TID)
ccc      CALL JULIAN(TFIN,TFD)
ccc      CALL JULIAN(TREF,TRD)
ccc	write(6,*)'julian:',tid,tfd,trd
C     WRITE(7,6000)TID,TFD,TRD
C
C  UNIT CONVERSION FROM DEG TO RADIAN AND DAY TO SECONDS
C
      PM=PM*DTR
      RATE=RATE*DTR
      DO 60 I=3,6
      ORB(I)=ORB(I)*DTR
      ES(I)=ES(I)*DTR
      EM(I)=EM(I)*DTR
   60 CONTINUE
      ES(7)=ES(7)*DTR
      EM(7)=EM(7)*DTR
cAC      TI=TID*DTS
cAC      TF=TFD*DTS
cAC      TR=TRD*DTS
	ti=tints	!Input passing A.Chu 
	tf=ti+delta_t
	tr=ti		!reference time in sec
C
C  CONVERT MEAN ELEMENTS TO OSCULATING ELEMENTS IF NECESSARY
C
      IF (IORB.EQ.1) THEN
        CALL KOZSAK(1,GE,RE,AJ2,ORB,X)
        DO 70 I=1,NEQ
   70   ORB(I)=X(I)
      ENDIF
C
C  CHANGE INPUT ORBITAL ELEMENTS TO CARTESIAN COORDINATES
C
      DO 80 I=1,NEQ
   80 Y(I)=ORB(I)
      CALL KEPLER(Y(6),Y(2),EA,SE,CE)
      Y(6)=EA
      CALL COORD(Y,GE,X)
C
C  SET UP ROTATIONAL MATRIX FOR ANALYTICAL EPHEMERIDES FOR THE SUN
C  AND THE MOON.  BESIDES SUN AND SRP, SETSUN IS REQUIRED FOR MOON
C  PERTURBATION BECAUSE THE ECLIPTIC ANGLE IS NEEDED FOR TRANSFORMATION
C  FROM EMO OF DATE TO EME OF DATE
C
      IF (IEPHEM.EQ.1) THEN
        CALL SETSUN(TI,ES)
      ENDIF
      IF (ISUN.EQ.1.OR.ISRP.EQ.1) CALL SETTHD(ES,ET)      
      IF (IMOON.EQ.1.AND.IEPHEM.EQ.0) CALL SETTHD(EM,EN)
C
C  SET UP INTEGRATOR PARAMETERS
C
      H=HSTART
      T=TI
      CALL RK78CN
C
C  THE VARIABLE SP GIVES AN APPROX. TIME (SEC) WHEN THE ROUTINE USEROP
C  IS CALLED.  USEROP IS A USER OPTION ROUTINE.  IT IS USED TO COMPUTE
C  CERTAIN STATE OF THE TRAJECTORY TO DETECT PERIAPSIS OR APOAPSIS
C  PASSAGE.  WHEN SP IS SET TO A LARGE NUMBER, IT WON'T BE CALLED.  WHEN
C  SET TO ZERO, IT WILL BE CALLED EVERY INTEGRATION STEP.  USEROP CAN
C  BE USED TO COMPUTE OTHER QUATITIES DURING THE PROPAGATION SUCH AS
C  SHADOWING.  YOU CAN ALSO USE USEROP TO PERFORM MANEUVERS AT A
C  PARTICULAR TIME OR DO IT AUTOMATICALLY WHEN CERTAIN ORBITAL STATE
C  IS REACHED SUCH AS ECCENTRICITY EXCEEDS CERTAIN VALUES.
C
      IF (IPRINT.EQ.0.AND.INODE.EQ.0) THEN
        SP=HLARGE
      ELSE
        SP=ZERO
      ENDIF
C
C  SET UP PLOTTING FILE
C
      IF (IPLOT.EQ.1) OPEN (8,FILE='8')
C
C  THIS IS SET UP TO PRINT AT FIXED STEP INTERVALS, SEE SP EXPLANATION
C  IF YOU WANT TO DO ADVANCE PROGRAMMING
C
c	A.Chu 8/31/94
cAC  500 CONTINUE
      CALL POUT(T,X)
c      TOUT=T+STEP   
cAC      IF (TOUT.GT.TF) TOUT=TF   
cAC1      CALL RK78(DER,T,TOUT,NEQ,X,H,RELERR,ABSERR,SP)
        CALL RK78(DER,T,TF,NEQ,X,H,RELERR,ABSERR,SP)
	do i=1,6
	   x(i)=x(i)*1000.d0
	enddo
	
cAC      IF (TOUT.EQ.TF) GO TO 900 
cAC      GO TO 500 
cAC  900 CONTINUE  
      CALL POUT(T,X)
 3000 FORMAT(1P,BN,D30.16)
 4000 FORMAT(BN,I5)
 5000 FORMAT(1P,BN,2I5,2D30.16)
 6000 FORMAT(1P,1H1,/
     1      ,5X,'RUN STARTS ON JULIAN DATE             = ',D25.16,/
     2      ,5X,'RUN ENDS ON JULIAN DATE               = ',D25.16,/
     3      ,5X,'REFERENCE JULIAN DATE OF PM AND EPHEM = ',D25.16)
cAC      CLOSE (5)
      CLOSE (8)
      CLOSE (7)
      END
       subroutine elevation(lat,lon,h_update,avg_terrain,
     *  dem_ellip,dlon_topo,dlat_topo,topo )

	implicit real*8 (a-h, o-z)
c	include 'key_const.inc'
c	include 'key_pp.inc'
C*****WAYNE******
       include 'ssp2_const.inc'
       integer     dem_ellip
       real*8     dlon_topo
       real*8     dlat_topo
       integer*2     topo(ns_ew,ns_ns)
C*****WAYNE******
	real*8 lat,lon
        real*8 avg_terrain

	if(dem_ellip.eq.ellip) then
	h_update = avg_terrain
	else

        ai = mod(nint(360.d0+lon),360)/dlon_topo+1
c	ai = lon/dlon_topo+1
	aj = (90.d0-lat)/dlat_topo+1

	i = ai
	j = aj
	i1 = i + 1
	j1 = j + 1
	if(i1.gt.ns_ew) i1 = i1-ns_ew
	if(j1.gt.ns_ns) j1 = j1-ns_ns

	d1 = ai - i
	d2 = aj - j

	value = topo(i,j)*(1-d2)*(1-d1)
     *	      + topo(i1,j)*d2*(1-d1)
     *	      + topo(i,j1)*(1-d2)*d1
     *	      + topo(i1,j1)*d2*d1

	h_update = value
	end if

	return
	end
