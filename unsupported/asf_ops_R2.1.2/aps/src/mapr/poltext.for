C--  Copyright (c)1996, California Institute of Technology.
C--  U.S. Government Sponsorship acknowledged.
C-- ==========================================================================
C--
C--  Fortran Filename:	poltext.for
C--
C--  Description:	
C--	
C--  Notes:
C--
C-- ==========================================================================

C-----------------------------------------------------------------
C SUBROUTINE POLTEXT
C
C PURPOSE
C	DISPLAY THE LAT/LON NUMBERS ON THE POLAR STEREOGRAPHIC PROJECTION
C
C $Logfile:   ACS003:[BLD.MPS.MAPR.SRC]POLTEXT.FOV  $
C
C INPUT
C	PROJN		PROJECTION NO.
C	MINLAT,MAXLAT	MIN/MAX LAT OF WINDOW
C	MINLON,MAXLON	MIN/MAX LON OF WINDOW
C	OBSLAT,OBSLON	CENTER LAT/LON
C	LATINCR,LONINCR	LAT/LON INCREMENTERS
C	STRTLT,STOPLT,
C	STRTLN,STOPLN	START/STOP LAT/LON
C INTERNAL
C	SIGN
C	N		LOOP INDEX
C	DIFLAT,DIFLON	HEIGHT/WIDTH OF WINDOW
C	DX,DY
C	LATLINE,LONLINE
C
C WRITTEN BY CRAIG K. FUJIMOTO
C
C MODIFICATION
C $Date$ $Revision$ $Author$
C-----------------------------------------------------------------
      SUBROUTINE POLTEXT(PROJN,OBSLAT,OBSLON,
     1                   MINLAT,MAXLAT,MINLON,MAXLON,
     2                   LATINCR,LONINCR,
     3                   STRTLT,STOPLT,STRTLN,STOPLN)

      character*100 SccsFileID
     -/'@(#)poltext.for	5.1 98/01/08 APS/ASF\0'/

      IMPLICIT NONE

C INPUT
c--port      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'

       INCLUDE 'GKS_ROOT:include/fgksenum.inc'

      INTEGER PROJN
      REAL OBSLAT,OBSLON
      REAL MINLAT,MAXLAT,MINLON,MAXLON
      REAL LATINCR,LONINCR
      REAL STRTLT,STOPLT,STRTLN,STOPLN
C INTERNAL
      REAL SIGN,N
      REAL DIFLAT,DIFLON
      REAL DX,DY
      REAL LATLINE,LONLINE
C-----------------------------------------------------------------

      IF (PROJN .EQ. 5) THEN
        SIGN = 1.0
      ELSE
        SIGN = -1.0           
      END IF

C DO LAT LINE NUMBERS

C SET UP ALIGNMENT
      CALL GSTXAL(GACENT,GAHALF)

C SET DISPLACEMENTS
      DX = 0.0
      DY = 0.0

C DETERMINE WHAT LON LINE LAT NUMBERS WILL BE DISPLAYED ALONG
      CALL F_LONLINE(PROJN,OBSLAT,OBSLON,
     1             MINLON,MAXLON,LONLINE)

C IF THE DISTANCE BETWEEN THE MINLAT AND STRTLT IS LARGE ENOUGH,
C  DISPLAY USER SPECIFIED MINLAT TEXT

      DIFLAT = ABS(MINLAT-STRTLT)

      IF (DIFLAT .GT. (LATINCR / 3.0)) THEN
        CALL P_TXT(PROJN,OBSLAT,OBSLON,
     1             MINLAT,LONLINE,0,DX,DY)
      END IF

C IF THE DISTANCE BETWEEN THE MAXLAT AND STOPLT IS LARGE ENOUGH,
C DISPLAY USER SPECIFIED MAXLAT TEXT

      DIFLAT = ABS(MAXLAT-STOPLT)

      IF (DIFLAT .GT. (LATINCR / 3.0)) THEN
        CALL P_TXT(PROJN,OBSLAT,OBSLON,
     1             MAXLAT,LONLINE,0,DX,DY)
      END IF

C DISPLAY TEXT FOR LAT LINES INBETWEEN

      DO 1200 N = STRTLT,STOPLT,-SIGN*LATINCR

        IF (N .NE. 0) THEN
          CALL P_TXT(PROJN,OBSLAT,OBSLON,
     1             N,LONLINE,0,DX,DY)
        END IF

 1200 CONTINUE
      

C D0 LON LINE NUMBERS

C SET UP ALIGNMENT
C      CALL GSTXAL(GACENT,GAHALF)

C DETERMINE DISPLACEMENTS
C      CALL TXTDISP(DIFLAT,X,Y)
      DX = 0.0
      DY = 0.0

C DETERMINE WHICH LAT LINE THE LON NUMBERS WILL BE DISPLAYED AT
      LATLINE = MAXLAT

C IF DISTANCE BETWEEN MINLON AND STRTLON IS LARGE ENOUGH,
C DISPLAY THE USER SPECIFIED MINLAT TEXT

      DIFLON = ABS(MINLON - STRTLN)

      IF (DIFLON .GT. (LONINCR / 3.0)) THEN
        CALL P_TXT(PROJN,OBSLAT,OBSLON,
     1             LATLINE,MINLON,1,DX,DY)
      END IF

C IF DISTANCE BETWEEN MINLON AND STRTLON IS LARGE ENOUGH,
C DISPLAY THE USER SPECIFIED MINLAT TEXT

      DIFLON = ABS(MAXLON - STOPLN)

      IF (DIFLON .GT. (LONINCR / 3.0)) THEN
        CALL P_TXT(PROJN,OBSLAT,OBSLON,
     1             LATLINE,MINLON,1,DX,DY)
      END IF

C FOR FULL CIRCLE
      IF ((STRTLN .EQ. STOPLN) .OR.
     1    (STRTLN .EQ. -180.0 .AND. STOPLN .EQ. 180.0) .OR.
     1    (STRTLN .EQ. 180.0 .AND. STOPLN .EQ. -180.0)) THEN

        DO 3000 N = -170.0,180.0,LONINCR

          CALL P_TXT(PROJN,OBSLAT,OBSLON,
     1               LATLINE,N,1,DX,DY)

 3000   CONTINUE

      ELSE

C FOR SEMI CIRCLES

        DO 4000 N = STRTLN,-SIGN*180.0,-SIGN*LONINCR
             
          CALL P_TXT(PROJN,OBSLAT,OBSLON,
     1               LATLINE,N,1,DX,DY)

          IF (N.EQ.STOPLN) THEN
            GO TO 4500
          END IF

 4000   CONTINUE

        DO 4100 N = SIGN*180.0,STOPLN,-SIGN*LONINCR

C          IF (N .EQ. STOPLN) THEN
C            GO TO 4500
C          END IF
              
          CALL P_TXT(PROJN,OBSLAT,OBSLON,
     1               LATLINE,N,1,DX,DY)

 4100   CONTINUE           
      
 4500   CONTINUE

      END IF

 9999 CONTINUE
      RETURN
      END


C-----------------------------------------------------------------
C SUBROUTINE F_LONLINE
C
C PURPOSE
C	DETERMINE WHICH LONGITUDE LINE TO USE TO DISPLAY
C	THE GRID TEXT LATITUDE NUMBERS AT.
C
C INPUT
C	PROJN		PROJECTION NO.
C	OBSLAT,OBSLON	CENTER LAT/LON
C	MINLON,MAXLON	MIN/MAX LON OF WINDOW
C	LONLINE
C INTERNAL
C	X,Y
C	HIDDEN,BRNCUT	LTRANS RETURN FLAGS
C
C WRITTEN BY CRAIG K. FUJIMOTO
C-----------------------------------------------------------------
      SUBROUTINE F_LONLINE(PROJN,OBSLAT,OBSLON,
     1             MINLON,MAXLON,LONLINE)

      IMPLICIT NONE

      INTEGER PROJN
      REAL OBSLAT,OBSLON
      REAL MINLON,MAXLON
      REAL LONLINE

C-----------------------------------------------------------------

      IF ((MINLON .EQ. MAXLON) .OR.
     1    (MINLON .EQ. -180.0 .AND. MAXLON .EQ. 180.0) .OR.
     1    (MINLON .EQ. 180.0 .AND. MAXLON .EQ. -180.0)) THEN

        LONLINE = 0.0

      ELSE
             
        LONLINE = MINLON

      END IF

C        CALL LTRANS(PROJN,OBSLAT,OBSLON,
C     1              MINLON,0.0,X(1),Y(1),
C     1              HIDDEN,BRNCUT)
C
C        CALL LTRANS(PROJN,OBSLAT,OBSLON,
C     1              MAXLON,0.0,X(2),Y(2),
C     1              HIDDEN,BRNCUT)
C
C        IF (X(1) .LT. X(2)) THEN
C          LONLINE = MINLON
C        ELSE
C          LONLINE = MAXLON
C        END IF
C
C      END IF

 9999 CONTINUE
      RETURN
      END

C-----------------------------------------------------------------
C SUBROUTINE P_TXT
C
C PURPOSE
C	PRINT TEXT AT THE SPECIFIED LAT AND LON BUT DISPLACE
C	IT BY GIVEN DX AND DY
C
C INPUT
C	PROJN		PROJECTION NUMBER
C	OBSLAT,OBSLON	CENTER LAT/LON
C	NUM		NUMBER TO BE DISPLAYED
C	LAT,LON		LAT/LON OF POINT TO DISPLAY TEXT AT
C	LLFLAG		 0 = LAT TEXT / 1 = LON TEXT
C	DX,DY
C INTERNAL
C	LEN
C	NUM
C	X,Y
C	CX,CY
C	NUMTXT
C	QD
C	HIDDEN,BRNCUT	LTRANS RETURNS
C
C WRITTEN BY CRAIG K. FUJIMOTO
C-----------------------------------------------------------------
      SUBROUTINE P_TXT(PROJN,OBSLAT,OBSLON,
     1                 LAT,LON,LLFLAG,DX,DY)

      IMPLICIT NONE

C INPUT
c--port      INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'

       INCLUDE 'GKS_ROOT:include/fgksenum.inc'

      INTEGER PROJN
      INTEGER LLFLAG

      REAL OBSLAT,OBSLON
      REAL LAT,LON
      REAL DX,DY

C OUTPUT
C INTERNAL
      INTEGER LEN
      REAL NUM
      REAL X(1),Y(1)
      REAL CX,CY,ABX,ABY

      CHARACTER*5 NUMTXT
      LOGICAL QD(1)
      LOGICAL HIDDEN,BRNCUT
C-----------------------------------------------------------------

C DETERMINE WHICH NUMBER IS TO BE DISPLAYED

C FOR LATS
      IF (LLFLAG .EQ. 0) THEN
        NUM = LAT
C FOR LONS
      ELSE
        NUM = LON
      END IF

C CONVERT THE CENTER POINT TO X/Y
      CALL LTRANS(PROJN,OBSLAT,OBSLON,
     1            OBSLON,OBSLAT,CX,CY,
     1            HIDDEN,BRNCUT)

C CONVERT THE POINT TO X/Y
      CALL LTRANS(PROJN,OBSLAT,OBSLON,
     1            LON,LAT,X(1),Y(1),
     1            HIDDEN,BRNCUT)

C CONVERT THE NUMBER TO ASCII
      ENCODE(4,'(I4)',NUMTXT) INT(NUM)
      CALL TRIM(NUMTXT,4,LEN)


      CALL QUADRANT(CX,CY,X,Y,1,QD)

C FOR LON POINTS, DETERMINE WHICH QUADRANT THE POINT IS IN
      IF (LLFLAG .EQ. 1) THEN

        IF (QD(1) .EQ. 1) THEN
          CALL GSTXAL(GALEFT,GABOTT)
        ELSE IF (QD(1) .EQ. 2) THEN
          CALL GSTXAL(GARITE,GABOTT)
        ELSE IF (QD(1) .EQ. 3) THEN
          CALL GSTXAL(GARITE,GATOP)
        ELSE
          CALL GSTXAL(GALEFT,GATOP)
        END IF

      ELSE

        ABX = ABS (CX - X(1))
        ABY = ABS (CY - Y(1))

        IF (QD(1) .EQ. 1) THEN
          IF (ABX .EQ. 0.0) THEN
            CALL GSTXAL(GARITE,GAHALF)
          ELSE IF (ABY .LT. ABX) THEN
            CALL GSTXAL(GARITE,GABOTT)
          ELSE
            CALL GSTXAL(GARITE,GABASE)
          END IF
        ELSE IF (QD(1) .EQ. 2) THEN
          IF (ABY .EQ. 0.0) THEN
            CALL GSTXAL(GACENT,GATOP)
          ELSE IF (ABY .LT. ABX) THEN
            CALL GSTXAL(GACENT,GATOP)
          ELSE
            CALL GSTXAL(GARITE,GACAP)
          END IF
        ELSE IF (QD(1) .EQ. 3) THEN
          IF (ABX .EQ. 0.0) THEN
            CALL GSTXAL(GALEFT,GAHALF)
          ELSE IF (ABY .LT. ABX) THEN     
            CALL GSTXAL(GALEFT,GATOP)
          ELSE
            CALL GSTXAL(GALEFT,GACAP)
          END IF
        ELSE IF (QD(1) .EQ. 4) THEN
          IF (ABY .EQ. 0.0) THEN
            CALL GSTXAL(GACENT,GABOTT)
          ELSE IF (ABY .LT. ABX) THEN
            CALL GSTXAL(GACENT,GABOTT)
          ELSE
            CALL GSTXAL(GALEFT,GABASE)
          END IF
        END IF

      END IF

      CALL GTXS(X(1)+DX,Y(1)+DY,LEN,NUMTXT)
                          
 9999 CONTINUE
      RETURN
      END
