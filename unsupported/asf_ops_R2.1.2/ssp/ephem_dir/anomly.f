c SccsId = @(#)anomly.f	2.41 3/24/98
      SUBROUTINE ANOMLY(NFLG,AIN,E,AOUT)
C
C/*    SUBROUTINE ANOMLY(NFLG,AIN,E,AOUT)  ----------------
C
C  VERSION 6/30/86
C  PURPOSE
C    TO CONVERT AMONG MEAN, ECCENTRIC AND TRUE ANOMALY
C  INPUT ARGUMENTS
C    NFLG    = 1, MEAN TO ECCENTRIC
C              2, MEAN TO TRUE
C              3, ECCENTRIC TO MEAN
C              4, ECCENTRIC TO TRUE
C              5, TRUE TO MEAN
C              6, TRUE TO ECCENTRIC
C    AIN     = INPUT ANGLE (RAD)
C    E       = ECCENTRICITY
C  OUTPUT ARGUMENT
C    AOUT    = OUTPUT ANGLE (RAD)
C  CALL SUBROUTINES
C    KEPLER
C  REFERENCE
C    JPL EM 312/87-153, 20 APRIL 1987
C  ANALYSIS
C    JOHNNY H. KWOK
C  PROGRAMMER
C    JOHNNY H. KWOK
C  MODIFICATIONS
C    NONE
C  COMMENTS
C    NONE
C
C*/
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA ONE,TWO/1.D0,2.D0/
      GO TO (11,11,12,12,13,13) NFLG
   11 AM=AIN
      GO TO 14
   12 EA=AIN
      GO TO 14
   13 F=AIN
   14 CONTINUE
      IF (NFLG.EQ.1.OR.NFLG.EQ.2) CALL KEPLER(AM,E,EA,SE,CE)
      IF (NFLG.EQ.1) GO TO 900
      IF (NFLG.EQ.2.OR.NFLG.EQ.4) THEN
        F=TWO*DATAN(DSQRT((ONE+E)/(ONE-E))*DTAN(EA/TWO))
        GO TO 900
      ENDIF
      IF (NFLG.EQ.5.OR.NFLG.EQ.6)
     1  EA=TWO*DATAN(DSQRT((ONE-E)/(ONE+E))*DTAN(F/TWO))
      IF (NFLG.EQ.6) GO TO 900
      IF (NFLG.EQ.5.OR.NFLG.EQ.3) THEN
        AM=EA-E*DSIN(EA)
        GO TO 900
      ENDIF
  900 CONTINUE
      GO TO (21,22,23,22,23,21) NFLG
   21 AOUT=EA
      GO TO 24
   22 AOUT=F
      GO TO 24
   23 AOUT=AM
   24 CONTINUE
      RETURN
      END
