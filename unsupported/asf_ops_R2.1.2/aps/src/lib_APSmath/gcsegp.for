C--  Copyright (c)1996, California Institute of Technology.
C--  U.S. Government Sponsorship acknowledged.

*************************************************************************
*  Name:	GCSEGP
*  Module Type: SUBROUTINE	Language: FORTRAN
*  $Logfile:   ACS003:[BLD.MPS.LIB.SRC]GCSEGP.FOV  $
*  Purpose:	DETERMINE THE CLOSEST POINT ON SEGMENT PL PR TO A POINT Q.
*		ALSO COMPUTE THE DISTANCE.  PL AND PR ARE CONSIDERED AS 
*		BEING ON THE LEFT AND ON THE RIGHT RESPECTIVELY.  YOU ARE 
*		STANDING IN THE MIDDLE OF THE LINE SEGMENT.  YOUR LEFT SIDE IS 
*		TOWARDS PL AND YOUR RIGHT SIDE IS TOWARDS PR.  IMAGINE A PLANE 
*		THAT SLICES THE EARTH INTO TWO EQUAL PARTS AND THE SEGMENT PL 
*		PR AND THE EARTH'S CENTER LIES IN THAT PLANE.  THE COMPUTED 
*		DISTANCE (D) IS POSITIVE IF Q IS IN THE HALF OF THE EARTH THAT 
*		IS IN FRONT OF YOU AND THE LINE SEGMENT.  IT IS ALSO POSITIVE 
*		IF IT IS ANYWHERE IN THE SAME PLANE AS PL PR.  
*		D IS NEGATIVE IF Q IS IN THE HALF OF THE EARTH THAT IS BEHIND 
*		YOU AND THE SEGMENT PL PR.  
*  Subroutines called:
*  GCBOMB, GCBTWN, GCDIST, GCLOSE, GCPERP
*  VECTOR LIBRARY: UCROSS, DOT
*  Input Parameters:
*  3 POINTS ON A UNIT SPHERE, (X,Y,Z)
*  Name         Type    Definition
*  PL(3)	REAL*8	LEFT POINT OF SEGMENT
*  PR(3)	REAL*8	RIGHT POINT OF SEGMENT
*  Q(3)		REAL*8	TEST POINT. 
*  Output Parameters:
*  Name         Type    Definition
*  D		REAL*8	DISTANCE FROM PL PR TO Q. 
*  PC		REAL*8	POINT ON PL PR THAT IS CLOSEST TO Q.  
*  Variables:
*  Locals :
*  Externals :
*  Modification History:                                            
*  Date			Revision	Author
*  $Date$ $Revision$ $Author$
*                                                                   
*********************************************************************/
      SUBROUTINE GCSEGP (D,PC,PL,PR,Q)
      character*100 SccsFileID
     -/'@(#)gcsegp.for	5.1 98/01/08 APS/ASF\0'/

      IMPLICIT NONE
      REAL*8 PL(3), PR(3), Q(3)
      REAL*8 D, PC(3)
      REAL*8 D1, D2, PP(3), C
C---	DOT IS A VECTOR LIBRARY FUNCTION.  
      REAL*8 DOT
      INTEGER IFLAG
C---	SPECIAL CASE:  PL=PR
      CALL GCLOSE(IFLAG,PL,PR)
      IF(IFLAG.NE.1) GO TO 1000
C---	IFLAG = 1 AND PL = PR.
      CALL GCDIST(D,PL,Q)
      PC(1) = PL(1)
      PC(2) = PL(2)
      PC(3) = PL(3)
      GO TO 9999
 1000 CONTINUE
C---	CHECK FOR UNDEFINED SEGMENT:  
      IF(IFLAG.EQ.-1) CALL GCBOMB
C---	OK. NORMAL CASE FOR PL PR.  
C---	CHECK FOR Q NORMAL TO THE P PLANE.  IF SO, THEN ALL PONTS ON PL PR ARE
C---	CLOSEST TO Q.  
C---	VECTOR PP IS NORMAL TO PL PR AND POINTING IN THE FORWARD DIRECTION.
      CALL UCROSS(PP,PL,PR)
      CALL GCLOSE(IFLAG,PP,Q)
      IF (IFLAG.EQ.0) GO TO 2000 
C---	Q IS NORMAL TO THE P PLANE.  
      CALL GCDIST(D,PL,Q)
      PC(1) = PL(1)
      PC(2) = PL(2)
      PC(3) = PL(3)
      GO TO 9000
 2000 CONTINUE
C---	O.K.  REGULAR CASE.
      CALL GCPERP(PC,PL,PR,Q)
      CALL GCBTWN(IFLAG,PL,PR,PC)
      IF(IFLAG.NE.1) GO TO 3000
C---	PC WAS BETWEEN PL AND PR.  PC IS OUR POINT, CLOSEST TO Q.
      CALL GCDIST(D,PC,Q)
      GO TO 9000
 3000 CONTINUE
C---	THE CLOSEST POINT TO Q IS EITHER PL OR PR.  
      CALL GCDIST(D1,PL,Q)
      CALL GCDIST(D2,PR,Q)
      IF(D1.GT.D2) GO TO 4000 
C---	D1 .LE. D2 HERE.  
      PC(1) = PL(1)
      PC(2) = PL(2)
      PC(3) = PL(3)
      D = D1
      GO TO 9000
 4000 CONTINUE
C---	D1 .GT. D2 HERE.  
      PC(1) = PR(1)
      PC(2) = PR(2)
      PC(3) = PR(3)
      D = D2
 9000 CONTINUE
C---	GET POSITIVE OR NEGATIVE DIRECTION:  
C---	USE PP, THE NORMAL TO PLANE PL PR IN THE POSITIVE DIRECTION.  
      C = DOT(PP,Q,3)
C---	IF THE COSINE OF THE ANGLE BETWEEN Q AND PP IS MORE THAN 90 DEGREES, 
C---	THEN Q IS BEHIND US AND THE DISTANCE IS NEGATIVE.  
      IF(C.LT.0.0D00) D = - D
 9999 CONTINUE
      RETURN
      END
