C--  Copyright (c)1996, California Institute of Technology.
C--  U.S. Government Sponsorship acknowledged.

*************************************************************************
*  Name:	SWENTQ
*  Module Type: SUBROUTINE	Language: FORTRAN
*  $Logfile:   ACS003:[BLD.MPS.LIB.SRC]SWENTQ.FOV  $
*  Purpose:	GIVEN A QUADRILATERAL SITE AND SWATH POINTS AT TIMES 
*		1 AND 2, COMPUTE THE FRACTIONAL DISTANCE FROM 1 TO 2 
*		THE SWATH MUST TRAVEL TO START COVERING THE SITE.  
*  Subroutines called:
*  GCBOMB, GCDIST, GCENTC
*  MUST LINK VECTOR LIBRARY, DOUBLE PRECISION:  DOT
*  Input Parameters:
*  Name         Type    Definition
*  Q1(3)	REAL*8	1ST POINT OF THE QUADRILATERAL
*  Q2(3)	REAL*8	2ND POINT OF THE QUADRILATERAL
*  Q3(3)	REAL*8	3RD POINT OF THE QUADRILATERAL
*  Q4(3)	REAL*8	4TH POINT OF THE QUADRILATERAL
*		NOTE THAT THESE POINTS ARE ENTERED CLOCKWISE.
*  SL1(3)	REAL*8	AT TIME 1, THE LEFT POINT OF THE SWATH.
*  SR1(3)	REAL*8	AT TIME 1, THE RIGHT POINT OF THE SWATH.
*  SL2(3)	REAL*8	AT TIME 2, THE LEFT POINT OF THE SWATH.
*  SR2(3)	REAL*8	AT TIME 2, THE RIGHT POINT OF THE SWATH.
*  Output Parameters:
*  Name         Type    Definition
*  F		REAL*8	THE FRACTIONAL DISTANCE FROM 1 TO 2 THAT THE SWATH 
*			MUST TRAVEL TO START COVERING THE SITE.
*			0.0 <= F < 1.0
*			IF THE SITE DOESN'T OVERLAP THE SWATH RECTANGLE, 
*			A VALUE OF -1.0D00 IS RETURNED.  
*  Variables:
*  Locals :
*  Externals :
*  Modification History:                                            
*  Date			Revision	Author
*  $Date$ $Revision$ $Author$
*                                                                   
*********************************************************************/
      SUBROUTINE SWENTQ (Q1,Q2,Q3,Q4,SL1,SR1,SL2,SR2,F)
      character*100 SccsFileID
     -/'@(#)swentq.for	5.1 98/01/08 APS/ASF\0'/

      IMPLICIT NONE
      REAL*8 Q1(3), Q2(3), Q3(3), Q4(3)
      REAL*8 SL1(3), SR1(3), SL2(3), SR2(3), F
      REAL*8 FTRY(6)
      INTEGER IFLAG, J
      F = 0.0D0
C---	SEE IF THE SWATH COVERS THE QUADRILATERAL FROM THE START.  
      CALL SWINQ(SL1,SR1,SL2,SR2,Q1,Q2,Q3,Q4,IFLAG)
      IF (IFLAG .EQ. 1) GO TO 9999
      F = -1.0D00
C---	THE SWATH COVERAGE DOES NOT START RIGHT AWAY; VALUES OF F = 0.0
C---	WILL BE REJECTED AFTER THIS.  
      DO 1100 J = 1,6
 1100 FTRY(J) = -1.0D0
C
C---	GET F FOR THE LEFT BORDER OF THE SWATH.  
      CALL GCLINQ(F,SL1,SL2,Q1,Q2,Q3,Q4)
      IF(F .EQ. 0.0D00) GO TO 9999
      FTRY(1) = F
C---	GET F FOR THE RIGHT BORDER OF THE SWATH.  
      CALL GCRINQ(F,SR1,SR2,Q1,Q2,Q3,Q4)
      IF(F .EQ. 0.0D00) GO TO 9999
      FTRY(2) = F
C---	GET F FOR ANY QUADRILATERAL POINTS IN THE MIDDLE OF THE SWATH.
C---	TAKE THE SMALLEST OF THE F VALUES.  
C---	NOTE THAT THE SWATH IS INPUT TO GCPINQ AS A QUADRILATERAL WITH 
C---	THE POINTS IN CLOCKWISE ORDER.  
      CALL GCPINQ(IFLAG,SL2,SR2,SR1,SL1,Q1)
C---	Q1 MUST BE INSIDE, NOT ON AN EDGE.  
      IF(IFLAG.NE.1) GO TO 2000
C---	USING THE SWATH ENTERS A CIRCLE ROUTINE WITH THE POINT AND RADIUS = 0.
      CALL SWENTC(Q1,0.0D00,SL1,SR1,SL2,SR2,F)
      IF(F .EQ. 0.0D00) GO TO 9999
      FTRY(3) = F
 2000 CONTINUE
      CALL GCPINQ(IFLAG,SL2,SR2,SR1,SL1,Q2)
      IF(IFLAG.NE.1) GO TO 3000
      CALL SWENTC(Q2,0.0D00,SL1,SR1,SL2,SR2,F)
      IF(F .EQ. 0.0D00) GO TO 9999
      FTRY(4) = F
 3000 CONTINUE
      CALL GCPINQ(IFLAG,SL2,SR2,SR1,SL1,Q3)
      IF(IFLAG.NE.1) GO TO 4000
      CALL SWENTC(Q3,0.0D00,SL1,SR1,SL2,SR2,F)
      IF(F .EQ. 0.0D00) GO TO 9999
      FTRY(5) = F
 4000 CONTINUE
      CALL GCPINQ(IFLAG,SL2,SR2,SR1,SL1,Q4)
      IF(IFLAG.NE.1) GO TO 5000
      CALL SWENTC(Q4,0.0D00,SL1,SR1,SL2,SR2,F)
      IF(F .EQ. 0.0D00) GO TO 9999
      FTRY(6) = F
 5000 CONTINUE
C---	NOW TAKE THE SMALLEST FTRY .GT. 0.  IF ALL ARE -1, NO COVERAGE.
      F = 100
      DO 6000 J = 1,6
      IF(FTRY(J) .LT. 0.0D00) GO TO 6000
      IF(FTRY(J) .LT. F) F = FTRY(J)
 6000 CONTINUE
      IF(F .GT. 1.0D00) F = -1.0D0
 9999 CONTINUE
      RETURN
      END
