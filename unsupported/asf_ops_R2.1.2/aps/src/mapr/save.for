C--  Copyright (c)1996, California Institute of Technology.
C--  U.S. Government Sponsorship acknowledged.
C-- ==========================================================================
C--
C--  Fortran Filename:	save.for
C--
C--  Description:	
C--	
C--  Notes:
C--
C-- ==========================================================================

C-----------------------------------------------------------------------
C SUBROUTINE SAVE
C
C PURPOSE
C	STORES THE CURRENT DISPLAY INTO A MAP WINDOW 
C	FILE. A HEADER FILE IS ALSO CREATED CONTAINING INFORMATION ABOUT
C	THE STORED DISPLAY.
C
C       MAP WINDOW FILE FILE NAME:  NAME.GKS
C       MAP WINDOW HEADER FILE NAME:  NAME.HDR
C
C $Logfile:   ACS003:[BLD.MPS.MAPR.SRC]SAVE.FOV  $
C
C INPUT   :  
C	MINLAT		MINIMUM ZOOM WINDOW LATITUDE (DEG)
C	MAXLAT		MAXIMUM ZOOM WINDOW LATITUDE (DEG)
C	MINLON		MINIMUM ZOOM WINDOW LONGITUDE (DEG)
C	MAXLON		MAXIMUM ZOOM WINDOW LONGITUDE (DEG)
C	PROJN		PROJECTION TYPE TO DISPLAY MAP
C	OBSLAT		CENTERING LATITUDE OF THE WINDOW (DEG)
C	OBSLON		CENTERING LONGITUDE OF THE WINDOW (DEG)
C	GRIDLN		DISPLAY GRID LINES
C	CRSEGN		STATUS ARRAY OF OVERLAYS
C	CRSEGT		NAME ARRAY OF OVERLAYS
C	NSEG 		NUMBER OF SEGMENTS IN THE WORKSTATION
C
C INTERNAL:
C	IOS		STATUS
C	UINIT		USER'S INITIALS
C	MWFIL		MAP WINDOW FILE NAME (.GKS)
C	MWHDR		MAP WINDOW HEADER FILENAME (.HDR)
C	COMMNT		COMMENTS
C	MWFM 		MAP WINDOW FILE ID NUMBER
C	MWFH 		MAP WINDOW HEADER FILE ID NUMBER
C	FN   		FILE NAME LENGTH
C	I    		DO LOOP INDEX                   
C
C                                    
C WRITTEN BY CRAIG K. FUJIMOTO
C
C MODIFICATION
C $Date$ $Revision$ $Author$
C 8/8/94 Nadia Adhami  replaced GKS$OPEN_WS with Prior GKS calls
C
C-----------------------------------------------------------------------
       SUBROUTINE SAVE (MWID,PROJN,GRIDLN,
     1                  MINLAT,MAXLAT,MINLON,MAXLON,
     2                  STMNLT,STMXLT,STMNLN,STMXLN,
     3                  OBSLAT,OBSLON,
     4                  NSEG,CRSEGN,CRSEGT,CRDARID)

      character*100 SccsFileID
     -/'@(#)save.for	5.1 98/01/08 APS/ASF\0'/

       IMPLICIT NONE

C INPUT:
c--port       INCLUDE 'SYS$LIBRARY:GKSDEFS.FOR'
c--port       INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'

       INCLUDE 'GKS_ROOT:include/fgksenum.inc'


       CHARACTER*(*) CRSEGT(*)

       INTEGER MWID,PROJN,GRIDLN
       INTEGER NSEG,CRSEGN(*),CRDARID(*)

       REAL MINLAT,MAXLAT,MINLON,MAXLON
       REAL OBSLAT,OBSLON
       REAL STMNLT,STMXLT,STMNLN,STMXLN

C INTERNAL:
       CHARACTER*23 DATE
       CHARACTER*30 MWFIL,MWHDR
       CHARACTER*76 COMMNT

       INTEGER IOS,I

c--port declare external c routine: gettime 

C---      external gettime          !$pragma c(gettime)
      INCLUDE 'APS_HOME:include/local/timeconv.inc'

C-----------------------------------------------------------------------

      WRITE (6,100)
  100 FORMAT (/,' SAVE DISPLAY')

 1000 CONTINUE

C ASK FOR MAP WINDOW FILE AND HEADER FILE NAME
      CALL ASK_MFNAME(MWFIL,MWHDR,1,IOS)
      IF (IOS .NE. 0) GO TO 9999

C ASK FOR COMMENTS
      CALL ASK_COMM(COMMNT)

C GET THE SYSTEM TIME - ASF FORMAT
C---      CALL GETTIME(DATE)
      IOS = tc_systime2asf(DATE)

C MESSAGE TO USER
      CALL DISMSGW ('Creating Header file...')

C OPEN THE HEADER FILE
      OPEN (UNIT=10,FILE=MWHDR,STATUS='NEW',RECL=195,IOSTAT=IOS)

      IF (IOS .NE. 0) THEN
        CALL DISMSG('Error opening header file.')
      END IF

C WRITE HEADER INFORMATION

      WRITE (10,200) COMMNT,DATE,PROJN,GRIDLN,
     1               MINLAT,MAXLAT,MINLON,MAXLON,
     2               OBSLAT,OBSLON,
     3               STMNLT,STMXLT,STMNLN,STMXLN
  200 FORMAT (A76,A21,2(I2),10(F9.3))

C WRITE NAME OF GKS METAFILE
C      WRITE (10,220) MWFIL
C  220 FORMAT (1X,A57) 

C WRITE NAMES OF ALL SEGMENTS BEING STORED INCLUDING DARID OF SSC SEGS
      DO 2000 I = 1,NSEG

        IF (CRSEGN(I) .EQ. 1) THEN

          WRITE (10,230) CRSEGT(I),CRDARID(I)
  230     FORMAT (1X,A20,I5)

        END IF

 2000 CONTINUE


C CLOSE HEADER FILE

       CLOSE (10)

C MESSAGE TO USER
       CALL DISMSGW('Creating GKS Metafile...')

C OPEN AND ACTIVATE GKS METAFILE AND ASSOCIATE ALL SEGMENTS CURRENTLY 
C   IN THE WINDOW TO THE GKS METAFILE

c--port       CALL GKS$OPEN_WS (MWID,MWFIL,GKS$K_GKSM_OUTPUT)

       OPEN (UNIT=9, FILE=MWFIL , STATUS='UNKNOWN')
       CALL GOPWK( MWID, 9, 3 )

       CALL GACWK (MWID)

C ASSOCIATE GLOBE,WORLD AND GRID SEGMENTS TO THE METAFILE
       CALL DISMSGW('Saving world...')
       CALL GASGWK (MWID,500)
       CALL GASGWK (MWID,501)
       CALL GASGWK (MWID,502)
       CALL GASGWK (MWID,503)

       CALL DISMSGW('Saving grid...')
       CALL GASGWK (MWID,504)
       CALL GASGWK (MWID,505)
       CALL GASGWK (MWID,506)

C ASSOCIATE SEGMENTS TO THE GKS METAFILE
       IF (NSEG .GT. 0) THEN

         CALL DISMSGW('Saving overlays...')
         DO 3000 I = 1,NSEG

           IF (CRSEGN(I) .EQ. 1) THEN

             CALL GASGWK (MWID,I)
                                
           END IF

 3000    CONTINUE

       END IF

C DEACTIVATE AND CLOSE GKS METAFILE

       CALL GDAWK (MWID)
       CALL GCLWK (MWID)

       CALL DISMSG('Saving has completed.')

 9999  CONTINUE
       RETURN 
       END
