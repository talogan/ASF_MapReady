C--  Copyright (c)1996, California Institute of Technology.
C--  U.S. Government Sponsorship acknowledged.
C-- ==========================================================================
C--
C--  Fortran Filename:	zmseg.for
C--
C--  Description:	
C--	
C--  Notes:
C--
C-- ==========================================================================

C-----------------------------------------------------------------
C SUBROUTINE ZMSEG
C
C PURPOSE
C	REDISPLAY SEGMENTS ON ZOOMED MAP
C
C $LOGFILE$
C
C INPUT
C	DISP		DISPLAY 1=WORLD MAP / 2=METAFILE MAP
C	WSID		WORKSTATION ID
C	WISS		WORKSTATION IND. SEG. STORAGE ID
C	PROJN		PROJECTION NO.
C	GRIDLN		GRID LINE FLAG
C	WNSIZE		WINDOW SIZE
C	MINLAT,MAXLAT	MIN/MAX LAT OF WINDOW
C	MINLON,MAXLON	MIN/MAX LON OF WINDOW
C	OBSLAT,OBSLON	CENTER LAT/LON
C	OBSMIN,OBSMAX	MIN/MAX LON OF WINDOW
C	LONPT		LONGITUDE DIVISIONS
C	NSEG		OVERLAY COUNT
C	CRSEGN		OVERLAY STATUS ARRAY
C	CRSEGT		OVERLAY NAME ARRAY
C	MINX,MAXX,
C	MINY,MAXY	MIN/MAX X/Y COORDS
C	MINXWN,MAXXWN,
C	MINYWN,MAXYWN	MIN/MAX X/Y OF ENTIRE WINDOW
C	STMNLT,STMXLT,
C	STMNLN,STMXLN	MIN/MAX WINDOW LAT/LON
C
C OUTPUT
C	MINX,MAXX,
C	MINY,MAXY	MIN/MAX X/Y COORDS
C	MINXWN,MAXXWN,
C	MINYWN,MAXYWN	MIN/MAX X/Y OF ENTIRE WINDOW
C
C INTERNAL
C	XEDGE,YEDGE	ARRAY OF X/Y COORDS ON THE UNIT CIRCLE
C
C SUBROUTINES CALLED
C
C ORIGINALLY WRITTEN BY RICHARD P. LEE
C MODIFIED FOR ASF BY CRAIG K. FUJIMOTO - SEP 89
C
C MODIFICATIONS
C $DATE$ $REVISION$ $AUTHOR$
C-----------------------------------------------------------------
       SUBROUTINE ZMSEG (DISP,WSID,WISS,PROJN,GRIDLN,WNSIZE,
     1                   MINLAT,MAXLAT,MINLON,MAXLON,
     2                   OBSLAT,OBSLON,OBSMIN,OBSMAX,LONPT,
     4                   NSEG,CRSEGN,
     5                   STMNLT,STMXLT,STMNLN,STMXLN,
     6                   MINX,MAXX,MINY,MAXY,
     7                   MINXWN,MAXXWN,MINYWN,MAXYWN)

      character*100 SccsFileID
     -/'@(#)zmseg.for	5.1 98/01/08 APS/ASF\0'/

       IMPLICIT NONE

C INPUT: 
c--port       INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'

       INTEGER PROJN,DISP,WSID,WISS,GRIDLN
       INTEGER NSEG,CRSEGN(*)

       REAL OBSLAT,OBSLON
       REAL MINLAT,MAXLAT,MINLON,MAXLON
       REAL OBSMIN,OBSMAX
       REAL STMNLT,STMXLT,STMNLN,STMXLN
       REAL WNSIZE,LONPT(*)

C OUTPUT:
       REAL MINX,MAXX,MINY,MAXY
       REAL MINXWN,MAXXWN,MINYWN,MAXYWN

C INTERNAL:
       REAL DIFLAT,DIFLON

C-----------------------------------------------------------------------

C CALCULATE WINDOW DIMENSIONS

       CALL DIFDEG (OBSLON,MINLAT,MAXLAT,MINLON,MAXLON,
     1              DIFLAT,DIFLON)

C  SET THE MINIMUM AND MAXIMUM WINDOW LONGITUDE

       OBSMIN = MINLON
       OBSMAX = MAXLON
       IF (DIFLON .GE. 359.97) THEN
         OBSMIN = OBSLON + 180.01
         OBSMAX = OBSLON + 179.99
         MINLON = STMNLN
         MAXLON = STMXLN
       END IF
  950  IF (OBSMIN .GT. 360.0) THEN 
         OBSMIN = OBSMIN - 360.0
         GO TO 950
       END IF
  955  IF (OBSMAX .GT. 360.0) THEN 
         OBSMAX = OBSMAX - 360.0
         GO TO 955            
       END IF
  960  IF (OBSLON .GT. 360.0) THEN 
         OBSLON = OBSLON - 360.0
         GO TO 960
       END IF
            
C DISPLAY THE SEGMENTS

       CALL DISSEG (WSID,PROJN,GRIDLN,WNSIZE,
     1              MINLAT,MAXLAT,MINLON,MAXLON,
     2              OBSLAT,OBSLON,OBSMIN,OBSMAX,LONPT,
     4              NSEG,CRSEGN,
     5              MINX,MAXX,MINY,MAXY,
     6              MINXWN,MAXXWN,MINYWN,MAXYWN)

C END SUBROUTINE

 9999  CONTINUE
       RETURN
       END
