C--  Copyright (c)1996, California Institute of Technology.
C--  U.S. Government Sponsorship acknowledged.
C-- ==========================================================================
C--
C--  Fortran Filename:	polgrid.for
C--
C--  Description:	
C--	
C--  Notes:
C--
C-- ==========================================================================

C-----------------------------------------------------------------------
C SUBROUTINE POLGRID
C
C PURPOSE
C	DRAW THE GRID LINES FOR THE POLAR STEREOGRAPHIC PROJECTION
C
C $Logfile:   ACS003:[BLD.MPS.MAPR.SRC]POLGRID.FOV  $
C
C INPUT
C	PROJN		PROJECTION NO.
C	MINLAT,MAXLAT	MIN/MAX LAT OF WINDOW (INNER/OUTER LAT)
C	MINLON,MAXLON	MIN/MAX LON OF WINDOW (CLOCKWISE DIRECTION)
C	OBSLAT,OBSLON	CENTER LAT/LON
C
C OUTPUT
C	LATINCR,LONINCR	LAT/LON INCREMENTERS
C	STRTLT,STOPLT,
C	STRTLN,STOPLN	START/STOP LAT/LON
C
C INTERNAL
C	SIGN		SIGN -/+
C	N		LOOP COUNTER
C	LATDIST,LONDIT	LAT/LON DISTANCE
C
C SUBROUTINES CALLED
C
C WRITTEN BY CRAIG K. FUJIMOTO - SEP 89
C
C MODIFICATIONS
C $Date$ $Revision$ $Author$
C 4/24/95  Nadia Adhami  GSTXCI(1) -> (0) ; white -> black
C-----------------------------------------------------------------------
       SUBROUTINE POLGRID(PROJN,GRIDLN,OBSLAT,OBSLON,
     1                     MINLAT,MAXLAT,MINLON,MAXLON,
     1                     LATINCR,LONINCR,
     1                     STRTLT,STOPLT,STRTLN,STOPLN)

      character*100 SccsFileID
     -/'@(#)polgrid.for	5.1 98/01/08 APS/ASF\0'/

       IMPLICIT NONE

c--port       INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'


       INCLUDE 'GKS_ROOT:include/fgksenum.inc'

C INPUT
      INTEGER PROJN,GRIDLN
      REAL OBSLAT,OBSLON
      REAL MINLAT,MAXLAT,MINLON,MAXLON
C OUTPUT
      REAL STRTLT,STOPLT,STRTLN,STOPLN
      REAL LATINCR,LONINCR
C INTERNAL
      REAL SIGN
      REAL N
      REAL LATDIST,LONDIST
      REAL TICKLAT,TICKLON
      REAL LATTICK_SIZE,LONTICK_SIZE
C-----------------------------------------------------------------------

      IF (PROJN .EQ. 5) THEN
        SIGN = 1.0      
      ELSE 
        SIGN = -1.0
      ENDIF

      LATDIST = ABS(MAXLAT-MINLAT)

      CALL STLOND(PROJN,MINLON,MAXLON,LONDIST)

C DEFINE LATINCR
C IF >=90 DEGREES, DISPLAY GRID LINE EVERY 10 DEGREES
      IF (LATDIST .GE. 90.0) THEN
        LATINCR = 10.0
C IF >10 DEGREES, DISPLAY GRID LINE EVERY 5 DEGREES
      ELSE IF (LATDIST .GT. 10.0) THEN
        LATINCR = 5.0
C IF <=10 DEGREES, DISPLAY GRID LINE EVERY 1 DEGREE
      ELSE        
        LATINCR = 1.0
      END IF
                                     
C DEFINE LONTICK_SIZE
      IF (LATDIST .GE. 90.0) THEN
        LONTICK_SIZE = 3.0
      ELSE IF (LATDIST .GT. 50.0) THEN
        LONTICK_SIZE = 2.0
      ELSE IF (LATDIST .GT. 10.0) THEN
        LONTICK_SIZE = 2.0
      ELSE        
        LONTICK_SIZE = 0.5
      END IF

C DEFINE LONINCR
C IF >=90 DEGREES, DISPLAY GRID LINE EVERY 10 DEGREES
      IF (LONDIST .GT. 90.0) THEN
        LONINCR = 10.0
C IF >10 DEGREES, DISPLAY GRID LINE EVERY 5 DEGREES
      ELSE IF (LONDIST .GT. 10.0) THEN
        LONINCR = 5.0
C IF <=10 DEGREES, DISPLAY GRID LINE EVERY 1 DEGREE
      ELSE
        LONINCR = 1.0
      END IF

C DEFINE LATTICK_SIZE
      IF (LONDIST .GT. 90.0) THEN
        LATTICK_SIZE = 3.0
      ELSE IF (LONDIST .GT. 50.0) THEN
        LATTICK_SIZE = 2.0
      ELSE IF (LONDIST .GT. 10.0) THEN
        LATTICK_SIZE = 1.0
      ELSE
        LATTICK_SIZE = 0.5
      END IF


C DETERMINE THE START LATITUDE INDEX
      DO 1000 N = SIGN*90.0,0.0,-SIGN*LATINCR

        IF ((PROJN.EQ.5 .AND. N.LE.MINLAT) .OR. 
     1      (PROJN.EQ.6 .AND. N.GE.MINLAT)) THEN

          STRTLT = N
          GO TO 1050

        END IF

 1000 CONTINUE
 1050 CONTINUE

C DETERMINE THE STOP LATITUDE INDEX
      DO 1100 N = 0.0,SIGN*90.0,SIGN*LATINCR
        IF ((PROJN.EQ.5 .AND. N.GE.MAXLAT) .OR.
     1      (PROJN.EQ.6 .AND. N.LE.MAXLAT)) THEN
          STOPLT = N
          GO TO 1150
        END IF
 1100 CONTINUE
 1150 CONTINUE

C DRAW A GRID LINE AT THE SELECTED MINLAT
      CALL POL_ARC(PROJN,OBSLAT,OBSLON,MINLAT,MINLON,MAXLON)

C DRAW A GRID LINE AT THE SELECTED MAXLAT
      CALL POL_ARC(PROJN,OBSLAT,OBSLON,MAXLAT,MINLON,MAXLON)


C DETERMINE THE START LONGITUDE INDEX
      DO 2000 N = SIGN*180.0,-SIGN*180.0,-SIGN*LONINCR
        IF ((PROJN.EQ.5 .AND. N.LE.MINLON) .OR.
     1      (PROJN.EQ.6 .AND. N.GE.MINLON)) THEN
          STRTLN = N
          GO TO 2050
        END IF
 2000 CONTINUE
 2050 CONTINUE
                             
C DETERMINE THE STOP LONGITUDE INDEX
      DO 2100 N = -SIGN*180.0,SIGN*180.0,SIGN*LONINCR
        IF ((PROJN.EQ.5 .AND. N.GE.MAXLON) .OR.
     1      (PROJN.EQ.6 .AND. N.LE.MAXLON)) THEN
          STOPLN = N
          GO TO 2150
        END IF
 2100 CONTINUE
 2150 CONTINUE

C DRAW A GRID LINE AT THE SELECTED MINLON
      CALL POL_LINE(PROJN,OBSLAT,OBSLON,MINLAT,MINLON,MAXLAT,MINLON)

C DRAW A GRID LINE AT THE SELECTED MAXLON
      CALL POL_LINE(PROJN,OBSLAT,OBSLON,MINLAT,MAXLON,MAXLAT,MAXLON)



C WITH GRID LINES
      IF (GRIDLN .EQ. 1) THEN


C LATITUDE GRID LINES - CONCENTRIC CIRCLES

C DRAW THE ARCS INBETWEEN
        DO 3000 N = STRTLT,STOPLT,-SIGN*LATINCR

          CALL POL_ARC(PROJN,OBSLAT,OBSLON,N,MINLON,MAXLON)
 3000   CONTINUE


C LONGITUDE GRID LINES - RADIAL LINES

C FOR FULL CIRCLE
        IF (STRTLN .EQ. STOPLN) THEN

          DO 4000 N = -180.0,180.0,LONINCR

            CALL POL_LINE(PROJN,OBSLAT,OBSLON,MINLAT,N,MAXLAT,N)
 4000     CONTINUE

        ELSE

C DRAW THE GRID LINES INBETWEEN

          DO 5000 N = STRTLN,-SIGN*180.0,-SIGN*LONINCR
         
            IF (N.EQ.STOPLN) GO TO 5200

            CALL POL_LINE(PROJN,OBSLAT,OBSLON,MINLAT,N,MAXLAT,N)
 5000     CONTINUE

          DO 5100 N = SIGN*180.0,STOPLN,-SIGN*LONINCR

            CALL POL_LINE(PROJN,OBSLAT,OBSLON,MINLAT,N,MAXLAT,N)
 5100     CONTINUE           
      
 5200     CONTINUE

        END IF

C------------------------------------------------------------------
C NO GRID - TICK MARKS ONLY
C------------------------------------------------------------------

      ELSE

C LATITUDE TICK MARKS

C SET LON OF TICK MARK
        TICKLON = MINLON - SIGN * LONTICK_SIZE
        IF (TICKLON .GT. 180.0) THEN
          TICKLON = TICKLON - 360.0
        ELSE IF (TICKLON .LT. -180.0) THEN
          TICKLON = TICKLON + 360.0
        END IF

C DRAW THE TICK MARKS ALONG MINLON LINE
        DO 6000 N = STRTLT,STOPLT,-SIGN*LATINCR

          CALL POL_ARC(PROJN,OBSLAT,OBSLON,N,MINLON,TICKLON)

 6000   CONTINUE

C SET LON OF TICK ARC
        TICKLON = MAXLON + SIGN * LONTICK_SIZE
        IF (TICKLON .GT. 180.0) THEN
          TICKLON = TICKLON - 360.0
        ELSE IF (TICKLON .LT. -180.0) THEN
          TICKLON = TICKLON + 360.0
        END IF

C DRAW THE TICK MARKS ALONG MAXLON LINE
        DO 7000 N = STRTLT,STOPLT,-SIGN*LATINCR

          CALL POL_ARC(PROJN,OBSLAT,OBSLON,N,TICKLON,MAXLON)

 7000   CONTINUE


C LONGITUDE TICK MARKS


C DRAW FULL CIRCLE OF LONGITUDE TICKS
        IF (STRTLN .EQ. STOPLN) THEN

C SET TICK LATITUDE
          TICKLAT = MAXLAT + SIGN * LONTICK_SIZE

          DO 8000 N = -180.0,180.0,LONINCR

            CALL POL_LINE(PROJN,OBSLAT,OBSLON,TICKLAT,N,MAXLAT,N)

 8000     CONTINUE


        ELSE

C DO LONGITUDE TICK MARKS ALONG MAX LAT LINE

C SET TICK LATITUDE
          TICKLAT = MAXLAT + SIGN * LONTICK_SIZE

          DO 9000 N = STRTLN,-SIGN*180.0,-SIGN*LONINCR
                          
            IF (N.EQ.STOPLN) GO TO 9500

            CALL POL_LINE(PROJN,OBSLAT,OBSLON,TICKLAT,N,MAXLAT,N)

 9000     CONTINUE

          DO 9100 N = SIGN*180.0,STOPLN,-SIGN*LONINCR

            CALL POL_LINE(PROJN,OBSLAT,OBSLON,TICKLAT,N,MAXLAT,N)

 9100     CONTINUE           
      
 9500     CONTINUE


C DO LONGITUDE TICK MARKS ALONG MIN LAT LINE

C DON'T DO TICK MARKS AT 90/-90 LATITUDE
          IF (MINLAT .EQ. (SIGN*90.0)) GO TO 9600

C SET TICK LATITUDE
          TICKLAT = MINLAT - SIGN * LONTICK_SIZE

          DO 9200 N = STRTLN,-SIGN*180.0,-SIGN*LONINCR

            IF (N.EQ.STOPLN) GO TO 9600

            CALL POL_LINE(PROJN,OBSLAT,OBSLON,MINLAT,N,TICKLAT,N)

 9200     CONTINUE

          DO 9300 N = SIGN*180.0,STOPLN,-SIGN*LONINCR

            CALL POL_LINE(PROJN,OBSLAT,OBSLON,MINLAT,N,TICKLAT,N)

 9300     CONTINUE           
      
 9600     CONTINUE

        END IF

      END IF

 9999 CONTINUE
      RETURN
      END

                                     
C-----------------------------------------------------------------------
C SUBROUTINE POLBORDER
C
C PURPOSE
C	DRAW THE TEXT ALONG THE BORDER FOR THE POLAR STEREO PROJECTION
C
C INPUT
C	WSID		WORKSTATION ID
C	PROJN		PROJECTION NUMBER
C	OBSLAT,OBSLON	CENTER LAT/LON
C INTERNAL
C	LEN
C	DR
C	I
C	LATTXT,LONTXT	LAT/LON TEXT
C	DEGS,RADS
C	ONCX,ONCY
C	X,Y
C	CHARX,CHARY
C	DUMR
C	RADIUS
C	COUNTER
C	CENTLAT,CENTLON
C	PTLAT,PTLON
C	SIGN
C	HIDDEN,BRNCUT
C
C WRITTEN BY CRAIG K. FUJIMOTO
C-----------------------------------------------------------------------
       SUBROUTINE POLBORDER(WSID,PROJN,OBSLAT,OBSLON)

       IMPLICIT NONE

c--port       INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'

       INCLUDE 'GKS_ROOT:include/fgksenum.inc'
       INCLUDE 'mapper_port.inc'

C INPUT
      INTEGER WSID,PROJN
      REAL OBSLAT,OBSLON
C INTERNAL
      INTEGER LEN
      INTEGER DR,I

      CHARACTER*4 LATTXT,LONTXT

      REAL PI /3.14159265359/
      REAL X(2),Y(2)
      REAL CHARX,CHARY

      REAL CENTLAT,CENTLON
      REAL PTLAT,PTLON,SIGN
      REAL  GET_DEF_CH

      LOGICAL HIDDEN,BRNCUT

C POLY LINE COLOR INDEX
      CALL GSPLCI (5)
C TEXT COLOR INDEX WHITE -> BLACK
C     CALL GSTXCI(1)
      CALL GSTXCI(0)
C TEXT CHARACTER HEIGHT
      CALL GSCHH(GET_DEF_CH())
C TEXT ALIGNMENT
      CALL GSTXAL(GACENT,GAHALF)

      IF (PROJN .EQ. 5) THEN
        SIGN = 1.0
      ELSE IF (PROJN .EQ. 6) THEN
        SIGN = -1.0
      ENDIF

      CENTLAT = SIGN * 90.0
      CENTLON = 0.0

      CALL LTRANS(PROJN,OBSLAT,OBSLON,
     1            CENTLON,CENTLAT,
     2            X(1),Y(1),
     3            HIDDEN,BRNCUT)
      
C GRID LINES - LATITUDES - CONCENTRIC CIRCLES

      PTLON = 0.0

      DO 100 I = 1,8

        PTLAT = SIGN * I * 10.0

        CALL LTRANS(PROJN,OBSLAT,OBSLON,
     1              PTLON,PTLAT,
     1              X(2),Y(2),
     1              HIDDEN,BRNCUT)

        CALL GGDP (2,X,Y,GGCCP,DR,0)
 
  100 CONTINUE                
                        

C GRID LINE - LONGITUDES - RADIAL LINES

      PTLAT = 0.0

      DO 200 I = 0,35

        PTLON = SIGN * I * 10.0
      
        CALL LTRANS(PROJN,OBSLAT,OBSLON,
     1              PTLON,PTLAT,
     1              X(2),Y(2),
     1              HIDDEN,BRNCUT)

        CALL GPL(2,X,Y)

  200 CONTINUE


C LONGITUDE LINE NUMBER TEXT

C PLACE LONGITUDE NUMBERS 2 DEGREES FROM EDGE OF MAP
      PTLAT = -1.0 * SIGN * 2.0

      DO 300 I = -16,18,2

        PTLON = 10.0 * I

        CALL LTRANS(PROJN,OBSLAT,OBSLON,
     1              PTLON,PTLAT,
     1              CHARX,CHARY,
     1              HIDDEN,BRNCUT)

        ENCODE(4,'(I4)',LONTXT) INT(PTLON)
        CALL TRIM(LONTXT,4,LEN)

        CALL GTXS(CHARX,CHARY,LEN,LONTXT)

  300 CONTINUE

C LATITUDE LINE NUMBER TEXT

C PLACE LATITUDE LINES ALONG 0 LONGITUDE
      PTLON = 360.0

      DO 400 I = 1,8
                                
        PTLAT = SIGN * I * 10.0

        CALL LTRANS(PROJN,OBSLAT,OBSLON,
     1              PTLON,PTLAT,
     1              CHARX,CHARY,
     1              HIDDEN,BRNCUT)

        ENCODE(4,'(I4)',LATTXT) INT(PTLAT)
        CALL TRIM(LATTXT,4,LEN)

        CALL GTXS(CHARX,CHARY,LEN,LATTXT)

  400 CONTINUE


 9999 CONTINUE
      RETURN
      END
