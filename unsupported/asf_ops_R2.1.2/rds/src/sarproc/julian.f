      SUBROUTINE JULIAN(TIN,TOUT)
C/*   SUBROUTINE JULIAN(TIN,TOUT) ---------------
C
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
C*/
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TIN(2)
      DATA D0,D1,D2,D3,D4,D5/2415020.5D0,1.D-4,1.D-2,24.D0,1.44D3
     1,8.64D4/
      character*80 sccsid
      data sccsid /'@(#)julian.f	1.3 96/04/09 22:51:51\0'/
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
