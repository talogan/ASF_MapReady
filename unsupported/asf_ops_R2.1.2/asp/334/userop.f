C Alaska SAR Processor (ASP) %W% %E% %U%
      SUBROUTINE USEROP(T,X,N,H,NREJ,NSTP,ISTART,IFLAG)
C/*   SUBROUTINE USEROP(T,X,N,H,NREJ,NSTP,ISTART,IFLAG) ----------
C
C  VERSION OF 3/31/87
C  PURPOSE
C    SPECIAL INPUT/OUTPUT ROUTINE
C  INPUT
C    T      = CURRENT TIME (SEC)
C    X      = 6-D CARTESIAN COORD X,Y,Z,XD,YD,ZD (KM,KM/SEC)
C    DX     = 6-D DERIVATIVES OF POSITION AND VELOCITY (KM/SEC,KM/SEC**2)
C    N      = NUMBER OF EQUATIONS TO BE INTEGRATED
C    H      = CURRENT STEP-SIZE (SEC)
C    NREJ   = NUMBER OF STEPS REJECTED
C    NSTP   = NUMBER OF STEPS OF INTEGRATION TAKEN
C    ISTART = 0, FLAG TO REJECT LAST STEP TAKEN BY SETTING TO 1
C    IFLAG  = 0, FLAG TO ACCEPT LAST STEP TAKEN BY SETTING TO 1
C  OUTPUT
C    T      = TIME TO RESTART INTEGRATION
C    X      = STATE TO RESTART INTEGRATION
C    H      = STEP-SIZE TO RESTART INTEGRATION
C    ISTART = RESET TO 1 ONLY TO REJECT LAST STEP TAKEN
C    IFLAG  = RESET TO 1 ONLY TO ACCEPT LAST STEP TAKEN UNCONDITIONALLY
C             (OTHERWISE, THERE MAY A CHANCE THAT THE STEP WILL BE
C             REJECTED BECAUSE OF NOT MEETING ERROR TOLERANCE)
C  CALLED BY SUBROUTINES
C    RK
C  CALL SUBROUTINES
C    EQNOX, POUT
C  REFERENCES
C    JPL EM 312/87-153, 20 APRIL 1987
C  ANALYSIS
C    J. H. KWOK - JPL
C  PROGRAMMER
C    J. H. KWOK - JPL
C  PROGRAM MODIFICATIONS
C    NONE
C  COMMENTS
C    THIS ROUTINE IS USED TO PRINT ORBIT INFORMATION AT PERIAPSIS AND
C    APOAPSIS, AND NODAL CROSSING.  AT EVERY INTEGRATION STEP, THE DOT
C    PRODUCT OF RADIUS AND VELOCITY IS CHECKED TO DETECT EITHER A
C    PERIAPSIS OR APOAPSIS PASSAGE.  IF TRUE THEN THE TIME OF PERIASIS
C    IS COMPUTED AND THE INTEGRATION TRYS TO INTEGRATE TO THAT TIME.
C    FOR NODAL CROSSING, THE TIME OF NODAL CROSSING IS COMPUTED AND THE
C    INTEGRATION TRYS TO INTEGRATE TO THAT TIME.  THEN FLAGS ARE SET TO
C    PRINT APSES PASSAGE INFORMATION OR NODAL CROSSING INFORMATION.
C    THIS ROUTINE CAN BE MODIFIED TO PRINT AT OTHER CRITERIA
C    SUCH AS MAXIMUM AND MINIMUM ALTITUDE POINTS.  IT CAN ALSO BE
C    MODIFIED TO INCLUDE AUTOMATIC MANEUVERS BY CHECKING ORBIT
C    PARAMETERS, OR TO INCLUDE SOLAR RADIATION PRESUURE BY CHECKING
C    SHADOWING EFFECT, ETC.  OBVIOUSLY, THE TIME OF APSIS PASSAGE DOES
C    NOT WORK TOO WELL WITH CIRCULAR ORBITS.
C*/
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(6)
      DIMENSION XSAV(6),Y(18)
      COMMON/OPTION/L,M,IRES,ISUN,IMOON,IEPHEM,IDRAG,IDENS,ISRP,IORB
     1,IPRINT,INODE,IPLOT
      COMMON/TIME/TI,TF,TR
      COMMON/PLTCON/GE,RE,RATE,PM,AJ2,ELLIP,RATM
      COMMON/COUNT/NODE,NREV
      COMMON/CASE/ICASE
      DATA PI,TPI/3.141592653589793D0,6.283185307179586D0/
      DATA ZERO,HALF,ONE,TWO/0.D0,0.5D0,1.D0,2.D0/
      IF (IPRINT.EQ.0.AND.INODE.EQ.0) GO TO 990
C
C *** INITIALIZATION
C
      IF (NREJ.EQ.0.AND.NSTP.EQ.0) THEN
        IZ1=0
        ISN2=0
        ICASE=0
        IPNODE=0
        IPAPSE=0
        NSAV=0
      ENDIF
C
C *** IF THIS IS A REJECTED STEP, RETURN
C
      IF (NSTP.EQ.NSAV) GO TO 990
      IF (INODE.EQ.0) GO TO 100
C
C *** IF THIS IS A PRINT NODE STEP, THEN
C
      IF (IPNODE.EQ.1) THEN
C
C *** PRINT INFORMATION AT ASCENDING OR DESCENDING NODE
C
        IF (XSAV(3).LT.ZERO) THEN
          ICASE=1
        ELSE
          ICASE=2
        ENDIF
        CALL POUT(T,X)
        IPNODE=0
      ELSE
C
C *** CHECK FOR NODAL CROSSING
C
        IZ2=DSIGN(ONE,X(3))
        IF (IZ1*IZ2.EQ.-1) THEN
C
C *** IF NODAL CROSSING, COMPUTE NODAL CROSSING TIME
C
          CALL EQNOX(XSAV,GE,Y)
          D2=DSQRT((ONE-Y(7))/(ONE+Y(7)))
            IF (X(3).GT.ZERO) THEN
              F2=TPI-Y(10)
            ELSE
              F2=PI-Y(10)
            ENDIF
          EA2=TWO*DATAN(D2*DTAN(HALF*F2))
          AM2=EA2-Y(7)*DSIN(EA2)
          Y(11)=DMOD(Y(11)+TPI,TPI)
          AM2=DMOD(AM2+TPI,TPI)
          DM=DMOD(AM2-Y(11)+TPI,TPI)
          HNODE=DSQRT(Y(1)**3/GE)*(DM)
          IPNODE=1
        ENDIF
      ENDIF
  100 CONTINUE
      IF (IPRINT.EQ.0) GO TO 200
C
C *** IF THIS IS A PRINT APSIS STEP, THEN
C
      IF (IPAPSE.EQ.1) THEN
C
C *** PRINT INFORMATION AT PERIAPSIS OR APOAPSIS
C
        IF (ISN2.EQ.1) THEN
          ICASE=3
        ELSE
          ICASE=4
        ENDIF
        CALL POUT(T,X)
        IPAPSE=0
      ELSE
C
C *** IF WANTS APSIS INFO, THEN CHECK APSIS PASSAGE
C
        SN2=X(1)*X(4)+X(2)*X(5)+X(3)*X(6)
        ISN2=DSIGN(ONE,SN2)
        IF (ISN1*ISN2.EQ.-1) THEN
C
C *** COMPUTE APSIS CROSSING TIME
C
          HAPSE=(-SN1)*(T-TSAV)/(SN2-SN1)
          IPAPSE=1
        ENDIF
      ENDIF
  200 CONTINUE
C
C *** IF EVENT WILL OCCUR, PICK THE RIGHT STEP
C
      IF (IPNODE.EQ.1.OR.IPAPSE.EQ.1) THEN
        ISTART=1
        IFLAG=1
        IF (IPNODE.EQ.1) H=HNODE
        IF (IPAPSE.EQ.1) H=HAPSE
        IF (IPNODE.EQ.1.AND.IPAPSE.EQ.1) THEN
          H=DMIN1(HNODE,HAPSE)
          IF (H.EQ.HNODE) THEN
            IPAPSE=0
          ELSE
            IPNODE=0
          ENDIF
        ENDIF
      ENDIF
C
C *** IF ISTART=0 SAVE STEP
C
      IF (ISTART.EQ.0) THEN
        NSAV=NSTP
        SN1=SN2
        ISN1=ISN2
        IZ1=IZ2
        TSAV=T
        DO 910 I=1,N
  910   XSAV(I)=X(I)
      ENDIF
  990 CONTINUE
      RETURN
      END
