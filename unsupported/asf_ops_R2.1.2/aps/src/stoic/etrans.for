C--  Copyright (c)1996, California Institute of Technology.
C--  U.S. Government Sponsorship acknowledged.
C-- ==========================================================================
C--
C--  Fortran Filename:  etrans.for
C--
C--  Description:
C--
C--  Notes:
C--
C-- ==========================================================================

      DOUBLE PRECISION FUNCTION TC_ET2UTC(ET)
      character*100 SccsFileID
     -/'@(#)etrans.for	5.1 98/01/08 APS/ASF\0'/
 
********************************************************************
*  Name: TC_ET2UTC
*  Module Type: FUNCTION	Language: FORTRAN 77
*  Purpose: 	CONVERT EPHEMERIS(TDT) TIME TO UTC.  REAL JULIAN DAYS.
*
*		TDT = TERRESTRIAL DYNAMICAL TIME - REFERENCE FRAME IS THE EARTH.
*		      DIFFERENCE BETWEEN TDT AND UTC IS ALWAYS AN INTEGER 
*		      NUMBER OF SECONDS + 0.184. (except during a positive 
*		      leap second, when UTC stops for 1 second.  Therefore, the
*		      difference increases uniformly by one second during the 
*		      leap second, i.e., a 45 degree ramp on a graph of seconds
*		      difference vs seconds TDT.  )
*
*		ET FROM VECTOR LIBRARY HAS A REFERENCE FRAME DIFFERENT FROM THE 
*		EARTH.  (IT IS TDB, BARYCENTRIC DYNAMICAL TIME - AT THE 
*		BARICENTER OF THE SOLAR SYSTEM.  )
*		RELATIVISTIC EFFECTS CAUSE THE DIFFERENCE BETWEEN TDB AND TDT 
*		TO VARY.  
*		SEE THE ASTRONOMICAL ALMANAC 1989 (U.S.Govt Printing Office) 
*		SECTION B "Relationships between time-scales" 
*		FOR THE DIFFERENCE EXPRESSED AS AN EQUATION WITH THE MEAN 
*		ANOMALY OF THE EARTH AS THE THE INDEPENDENT VARIABLE.
c  Calls:
c     Etmtim - VECTOR LIBRARY
c  Programmer:
c     L. STEVENS
*  Input Parameters:
*  Name         Type    Definition
*  ET		REAL	REAL JULIAN DAYS.
*  Output Parameters:
*  TC_ET2UTC	REAL	UTC.  EXPRESSED AS REAL JULIAN DAYS.
*  Modification History:                                            
*                                                                   
*  Date			Revision	Author
*  $Date$ $Revision$ $Author$
c**********************************************************************
      IMPLICIT NONE
      REAL*8 ETMTIM, ET, XSEC, XSEC1, XSEC2, XDIF1, XDIF2, XFRAC
      INTEGER ISECS, IDAY


      XSEC = ETMTIM(ET,'UTCP')
C---	SAMPLE AT .0005 SEC BEFORE ET:
      XSEC1 = ETMTIM( (ET - .0005D0/86400.0D0),'UTCP')
C---	SAMPLE AT .0005 SEC AFTER ET:
      XSEC2 = ETMTIM( (ET + .0005D0/86400.0D0),'UTCP')
C---	NOW CHECK TO SEE IF WE ARE IN A LEAP SECOND AT ET:
      XDIF1 = ABS(XSEC1 - XSEC)
      XDIF2 = ABS(XSEC2 - XSEC)
      XFRAC = ( ET - INT(ET) )
      IF(XDIF1 .GE. .00049D0/86400.0D0 .AND. 
     ?   XDIF2 .GE. .00049D0/86400.0D0 .AND.
     ?   ABS(XFRAC - 0.50D0) .LT. 0.00138D0     ) THEN
C
C---	ET IS DURING LEAP SECOND:  ==> AT MIDNIGHT UTC <== 
C---		WE DEDUCE THIS FROM THE FACT THAT WE ARE:
C---		1) IN A RAMP OF INCREASING OR DECREASING SECONDS DIFFERENCE 
C---			BETWEEN ET AND UTC.  THE SLOPE OF THE RAMP IS ALMOST 1
C---			OR GREATER.  THE SLOPE EXISTS ON BOTH SIDES OF ET.  
C---			THIS IS NOT JUST A FLUCTUATION BETWEEN TDT AND TDB.
C---			(ETMTIM GIVES TDB-UTC, NOT TDT-UTC)
C---	  AND	2) ET IS WITHIN 2 MINUTES (0.00138 DAYS) OF MIDNIGHT (0.5 DAYS) 
C---			ET, WHEN LEAP EVENTS ARE SCHEDULED.  
C
C---		SINCE ET IS WITHIN A LEAP SECOND, WE JUST TAKE THE DAY INTEGER 
C---		AND FORCE MIDNIGHT (.5) ON THE RESULTING UTC:  LEAP EVENTS ARE 
C---		A STOPPAGE OF UTC AT EXACTLY MIDNIGHT UTC.  
         IDAY = ET
         TC_ET2UTC = (IDAY + 0.50D0)
C
      ELSE
C---		NOT DURING A LEAP SECOND.  
C---		HERE, ETMTIM GIVES BDT - UTC.  WHAT WE WANT IS TDT - UTC, 
C---		WHICH IS ALWAYS AN INTEGER SECONDS WITH A .184 DECIMAL.  
C---		BDT - UTC IS VERY NEAR THIS VALUE, WITHIN 0.001658 SECONDS.
C---		SO WE CAN TAKE THE INTEGER SECONDS FROM ETMTIM AND FORCE THE 
C---		0.184 ON THE DIFFERENCE.  
         ISECS = XSEC
         TC_ET2UTC = ET - (ISECS + 0.184D0) / 86400.0D0
C
      ENDIF
C
      RETURN
      END


      DOUBLE PRECISION FUNCTION TC_UTC2ET(UTC)
********************************************************************
*  Name: TC_UTC2ET
*  Module Type: FUNCTION	Language: FORTRAN 77
*  Purpose: CONVERT UTC TO EPHEMERIS TIME (TDT).  
*
*		TDT = TERRESTRIAL DYNAMICAL TIME - REFERENCE FRAME IS THE EARTH.
*		      DIFFERENCE BETWEEN TDT AND UTC IS ALWAYS AN INTEGER 
*		      NUMBER OF SECONDS + 0.184.
*		ET FROM VECTOR LIBRARY HAS A REFERENCE FRAME DIFFERENT FROM THE 
*		EARTH.  RELATIVISTIC EFFECTS CAUSE THE FRACTIONAL SECONDS 
*		DIFFERENCE BETWEEN ET AND UTC TO VARY FROM 0.184.  SEE THE 
*		ASTRONOMICAL ALMANAC 1989 SECTION B "Relationships between 
*		time-scales" FOR THE DIFFERENCE EXPRESSED AS AN EQUATION WITH
*		THE MEAN ANOMALY OF THE EARTH AS THE THE INDEPENDENT VARIABLE.
c  Calls:
c     VECTOR LIBRARY
*  Input Parameters:
*  Name         Type    Definition
*  UTC		REAL	JULIAN DAYS UTC
*  Output Parameters:
*  TC_UTC2ET	REAL	JULIAN DAYS EPHEMERIS TIME (TDT).
*  Modification History:                                            
*  Date			Revision	Author
*  $Date$ $Revision$ $Author$
c**********************************************************************
C---      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT NONE
      REAL*8 SPERD, UTC
      REAL*8 XET, XUTC, XET1, XET2, XUTC1, XUTC2
      REAL*8 TC_ET2UTC, DIFSEC
      DATA SPERD /86400.0D0/

C---	FIRST GUESS (XET) AT EPHEMERIS TIME.
      XET = UTC + ( UTC - TC_ET2UTC(UTC) )
C---
 1000 CONTINUE
C---	USE XET TO GET A RELATED UTC (XUTC).
      XUTC = TC_ET2UTC(XET)
C---	NOW XET AND XUTC ARE A LEGITIMATE ET/UTC PAIR.  
C---	SEE IF THIS XUTC IS WITHIN .0001 SECONDS FROM THE GIVEN UTC.  
C---	IF SO, THEN THE XET IS GOOD ENOUGH AND IS RATIFIED.  
      IF ( ABS(XUTC-UTC) .LT. (.0001D0/SPERD) ) GO TO 9000
C---	THE OBVIOUS TRY DID NOT WORK.  
C---	NOT CLOSE ENOUGH YET.  
C---	NOW MAKE 3 TRIES.  1) RAISE XET BY 1 SECOND.
C---			   2) RAISE XET BY THE DIFFERENCE OF (UTC-XUTC)
C---			   3) LOWER XET BY 1 SECOND.  (PROXIMITY TO A NEGATIVE 
C---				LEAP SECOND)
C---			   4) IF UTC IS IN THE LAST SECOND BEFORE MIDNIGHT 0.5,
C---				THEN SET XUTC TO THE SAME INTEGER DAY + 0.5 
C---				AND USE THE RESULTING ET.  (UTC IS WITHIN THE 
C---				SKIPPED SECOND OF UTC.)
C---	ONE OF THESE TRIES SHOULD WORK; IF NOT, BOMB THE PROGRAM.  
C---	COMPUTE TRY 1
      XET1 = XET + 1.0D0 / SPERD
      XUTC1 = TC_ET2UTC(XET1)
      IF ( ABS(XUTC1-UTC) .LT. (.0001D0/SPERD) ) XET = XET1 
      IF ( ABS(XUTC1-UTC) .LT. (.0001D0/SPERD) ) GO TO 9000
C---      WRITE(*,*) ' TRY 1 DIDN''T WORK.  '
C---	COMPUTE TRY 2.  
      XET2 = XET + (UTC - XUTC)
      XUTC2 = TC_ET2UTC(XET2)
      IF ( ABS(XUTC2-UTC) .LT. (.0001D0/SPERD) ) XET = XET2
      IF ( ABS(XUTC2-UTC) .LT. (.0001D0/SPERD) ) GO TO 9000
C---	COMPUTE TRY 3
      XET1 = XET - 1.0D0 / SPERD
      XUTC1 = TC_ET2UTC(XET1)
      IF ( ABS(XUTC1-UTC) .LT. (.0001D0/SPERD) ) XET = XET1 
      IF ( ABS(XUTC1-UTC) .LT. (.0001D0/SPERD) ) GO TO 9000
C---	CHECK FOR LAST SECOND.  
      DIFSEC = (UTC - INT(UTC)) * SPERD
      IF( DIFSEC .GT. 1.0D0 .OR. DIFSEC .LT. 0.0D0) GO TO 8888
C---	COMPUTE TRY 4
      XUTC = INT(UTC) + 0.5D0
      XET = XUTC + ( XUTC - TC_ET2UTC(XUTC) ) - 1.0D0/SPERD
      GO TO 9000
 8888 CONTINUE
      WRITE(*,*)' TC_UTC2ET:  ERROR:  TOO MANY ITERATIONS.  '
      WRITE(*,*)' TC_UTC2ET:  COULD NOT FIND THE ET FOR THIS UTC.  '
      WRITE(*,*)' TC_UTC2ET:  THIS CODE SHOULD BE REVIEWED SO THAT THIS'
      WRITE(*,*)'          CASE WILL WORK.'
      WRITE(*,*)' TC_UTC2ET:  UTC   = ', UTC
C
      WRITE(*,*)' TC_UTC2ET:  XUTC  = ', XUTC
      WRITE(*,*)' TC_UTC2ET:  XET   = ', XET
C
      WRITE(*,*)' TC_UTC2ET:  XUTC1 = ', XUTC1
      WRITE(*,*)' TC_UTC2ET:  XET1  = ', XET1
C
      WRITE(*,*)' TC_UTC2ET:  XUTC2 = ', XUTC2
      WRITE(*,*)' TC_UTC2ET:  XET2  = ', XET2
C
      WRITE(*,*)' TC_UTC2ET:  MPS SYSTEM DETECTED ERROR.  '
      WRITE(*,*)' TC_UTC2ET:  MPS SYSTEM DETECTED ERROR.  '
      WRITE(*,*)' TC_UTC2ET:  MPS SYSTEM DETECTED ERROR.  '
      WRITE(*,*)' TC_UTC2ET:  NOW BOMBING THE PROGRAM.  '
      CALL GCBOMB
 9000 CONTINUE
      TC_UTC2ET = XET
      RETURN
      END
*
      SUBROUTINE TC_JULIAN(TIN,TOUT)
C  VERSION OF 4/1/85
C  PURPOSE
C    COMPUTES JULIAN DATE WHEN GIVEN CALENDER DATE AND TIME
C  INPUT
C    TIN(1) = CALENDER DATE, YYYYMMDD.
C       (2) = CALENDER TIME, HHMMSS.SSSS
C  OUTPUT
C    TOUT   = JULIAN DATE (DAYS)
C  CALL SUBROUTINES
C    NONE
C  REFERENCES
C    JPL EM 312/87-153, 20 APRIL 1987
C  ANALYSIS
C    J. H. KWOK
C  PROGRAMMER
C    J. H. KWOK
C  PROGRAM MODIFICATIONS
C    NONE
C  COMMENTS
C    NONE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TIN(2)
      DATA D0,D1,D2,D3,D4,D5/2415020.5D0,1.D-4,1.D-2,24.D0,1.44D3
     1,8.64D4/
C
C  THIS IS THE JULIAN DATE FOR 19000101
C
      TOUT=D0
C
      IY=INT(TIN(1)*D1-1.9D3)
      IM=INT(TIN(1)*D2-1.9D5)-IY*100
      ID=INT(TIN(1)-1.9D7)-IY*10000-IM*100
      IHOUR=INT(TIN(2)*D1)
      IMIN=INT(TIN(2)*D2)-IHOUR*100
      SEC=TIN(2)-IHOUR*10000-IMIN*100
      JD=IY*365+(IY-1)/4
      IM1=IM-1
      IF (IM1.EQ.0) GO TO 12
      GO TO (1,2,3,4,5,6,7,8,9,10,11),IM1
    1 JD=JD+31
      GO TO 12
    2 JD=JD+59
      GO TO 12
    3 JD=JD+90
      GO TO 12
    4 JD=JD+120
      GO TO 12
    5 JD=JD+151
      GO TO 12
    6 JD=JD+181
      GO TO 12
    7 JD=JD+212
      GO TO 12
    8 JD=JD+243
      GO TO 12
    9 JD=JD+273
      GO TO 12
   10 JD=JD+304
      GO TO 12
   11 JD=JD+334
   12 CONTINUE
      IF (IY/4*4-IY.EQ.0.AND.IM.GT.2) JD=JD+1
      JD=JD+ID-1
      TOUT=TOUT+JD+IHOUR/D3+IMIN/D4+SEC/D5
      RETURN
      END




