C--  Copyright (c)1996, California Institute of Technology.
C--  U.S. Government Sponsorship acknowledged.

********************************************************************
*  Name:	GCROVR
*  Module Type: SUBROUTINE	Language: FORTRAN
*  $Logfile:   ACS003:[BLD.MPS.LIB.SRC]GCROVR.FOV  $
*  Purpose:	COMPUTE AT WHAT POINT THE SEGMENT SR1 SR2 INTERSECTS THE 
*		SEGMENT Q WHILE CROSSING OVER IT FROM THE NEAR SIDE
*		TO THE FAR SIDE.  SR1 SR2 IS CONSIDERED TO BE THE RIGHT
*		BOUNDARY OF A SENSOR SWATH WHEN APPROPRIATE FOR SPECIAL CASES.
*  Subroutines called:
*  GCBTWN, GCLOSE, GCBOMB
*  VECTOR LIBRARY: UCROSS, VUNIT, VUNITN
*  Input Parameters:
*  4 POINTS ON A UNIT SPHERE, (X,Y,Z)
*  Name         Type    Definition
*  SR1		REAL*8	FIRST POINT OF SEGMENT 
*  SR2		REAL*8	SECOND POINT OF SEGMENT 
*  QL		REAL*8	LEFT POINT OF SEGMENT Q.  
*  QR		REAL*8	RIGHT POINT OF SEGMENT Q.  
*  Output Parameters:
*  Name         Type    Definition
*  IFLAG	INTEGER	CONDITION FLAG
*			= 2 	SR1 IS ALREADY ON THE FAR SIDE OF THE Q PLANE.
*			= 1 	SR1 SR2 STARTS ON THE NEAR SIDE OR ON Q ITSELF 
*				(BUT NOT QL)
*				AND CROSSES OVER AND BEYOND THE SEGMENT Q, 
*				INTERSECTING THE SEGMENT AT A POINT NOT QL. 
*			= 0 	ALL OTHER CASES.  
*  X(3)		REAL*8	VALUE IS BASED ON IFLAG.  
*			IFLAG = 2:  X = SR1.
*			IFLAG = 1:  X = THE INTERSECTION OF S AND Q.  
*			IFLAG = 0:  THE 3 X VALUES ARE = 2.0
*  Variables:
*  Locals :
*  Externals :
*  Modification History:                                            
*  Date			Revision	Author
*  $Date$ $Revision$ $Author$
*                                                                   
*********************************************************************/
      SUBROUTINE GCROVR (IFLAG,X,SR1,SR2,QL,QR)
      character*100 SccsFileID
     -/'@(#)gcrovr.for	5.1 98/01/08 APS/ASF\0'/

      IMPLICIT NONE
      REAL*8 X(3), SR1(3), SR2(3), QL(3), QR(3)
      REAL*8 QP(3), SRP(3), X1(3)
C---	FUNCTION FROM VECTOR LIBRARY
      REAL*8 DOT
      INTEGER IFLAG, IFLAGS, IFLAGQ
      X(1) = 2.
      X(2) = 2.
      X(3) = 2.
      IFLAG = 0
C--	SPECIAL CASE:  SR1 = -SR2 OR QL = -QR:  BOMB THE PROGRAM
      CALL GCLOSE(IFLAGS,SR1,SR2)
      IF(IFLAGS.EQ.-1) CALL GCBOMB
      CALL GCLOSE(IFLAGQ,QL,QR)
      IF(IFLAGQ.EQ.-1) CALL GCBOMB
C--	SPECIAL CASE:  QL = QR:  RETURN NO INTERSECTION.
      IF(IFLAGQ.EQ.1) GO TO 9999
C---	QL .NE. QR.  
C---		QP = QL X QR	QP IS PERPENDICULAR TO PLANE Q
      CALL UCROSS(QP,QL,QR)
C--	SPECIAL CASE:  CHECK IF SR1 IS ALREADY BEYOND THE PLANE Q.  
C---			IF SO, RETURN IFLAG = 2.  
      IF(DOT(SR1,QP,3) .GT. 0.0D00) GO TO 9002
C--	SPECIAL CASE:  CHECK IF SR2 IS ON THIS SIDE OF THE PLANE Q OR ON THE 
C---	PLANE Q.  IF SO, RETURN IFLAG = 0.  NO INTERSECTION.
      IF(DOT(SR2,QP,3) .LE. 0.0D00) GO TO 9999
C---	NOW COMPUTE INTERSECTION OF SL AND Q.  IF THERE IS ONE, IFLAG = 1.
C---	ON A ZERO LENGTH SR1 SR2, NO CROSSOVER.
      IF(IFLAGS.EQ.1) GO TO 9999
C---		SRP = SR1 X SR2	SRP IS PERPENDICULAR TO PLANE SL
      CALL UCROSS(SRP,SR1,SR2)
      CALL GCLOSE(IFLAGS,QP,SRP)
C---	SPECIAL CASE:  IF SL AND Q ARE IN THE SAME PLANES, THEN NO CROSSOVER.
      IF(IFLAGS.NE.0) GO TO 9999
C---	FIRST DETERMINE A UNIT VECTOR X1 THAT LIES IN BOTH THE Q PLANE AND SR
C---	PLANE.  (THE Q PLANE IS DEFINED BY THE THREE POINTS QL, QR, AND THE 
C---	ORIGIN.)
C---		QP = QL X QR	QP IS PERPENDICULAR TO PLANE Q
C---		SRP = SR1 X SR2	SRP IS PERPENDICULAR TO PLANE SR
      CALL UCROSS(X1,SRP,QP)
C---		X1 = SRP X QP	X1 IS PERPENDICULAR TO BOTH SRP AND QP AND 
C---				THEREFORE LIES IN BOTH PLANES SR AND Q.
C---				THE ROUTINE UCROSS YIELDS UNIT VECTORS.
C---	NOW SEE IF EITHER X1 OR ITS NEGATIVE IS 
C---	BETWEEN SR1 AND SR2 IN THE PLANE.
C---	IFLAGS WILL BE 0 FOR NOT BETWEEN, 1 FOR X1 BETWEEN, -1 FOR -X1 BETWEEN.
      CALL GCBTWN(IFLAGS,SR1,SR2,X1)
      IF (IFLAGS.EQ.0 ) GO TO 9999
C---	OK.  X OR -X LIES IN SR1 SR2.
C---	NOW SEE IF EITHER X1 OR -X1 IS BETWEEN QL AND QR IN THE PLANE.
      CALL GCBTWN(IFLAGQ,QL,QR,X1)
      IF (IFLAGQ .EQ. 0 ) GO TO 9999
C---	OK.  X OR -X LIES IN QL QR.
C---	IS THE SAME POINT BETWEEN THEM BOTH?
      IF (IFLAGS .NE. IFLAGQ) GO TO 9999
C---	YES. EITHER X1 -X1 IS BETWEEN THEM BOTH AND LIES IN BOTH PLANES.
      IF (IFLAGS .EQ. 1)  CALL VUNIT(X,X1,3)
      IF (IFLAGS .EQ. -1) CALL VUNITN(X,X1,3)
      IFLAG = 1
C---	NOW CHECK FOR = QL.  
      CALL GCLOSE(IFLAGQ,X,QL)
      IF(IFLAGQ .NE. 1) GO TO 9999
C---	BACK TO IFLAG = 0
      X(1) = 2.
      X(2) = 2.
      X(3) = 2.
      IFLAG = 0
      GO TO 9999
 9002 CONTINUE
C---	CASE 2.   
      IFLAG = 2
      X(1) = SR1(1)
      X(2) = SR1(2)
      X(3) = SR1(3)
 9999 CONTINUE
      RETURN
      END
