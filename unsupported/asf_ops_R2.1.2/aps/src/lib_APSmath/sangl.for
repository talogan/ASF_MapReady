C--  Copyright (c)1996, California Institute of Technology.
C--  U.S. Government Sponsorship acknowledged.

********************************************************************
*  Name:	SANGL
*  Module Type: SUBROUTINE	Language: FORTRAN
*  $Logfile:   ACS003:[BLD.MPS.LIB.SRC]SANGL.FOV  $
*  Purpose:	COMPUTE SPHERICAL ANGLE IN RADIANS FROM THREE POINTS ON THE 
*		UNIT SPHERE
*  Subroutines called:
*  GCBOMB, GCLOSE
*  VECTOR LIBRARY: UCROSS, UCROSSM, DOT
*  Input Parameters:
*  3 POINTS ON A UNIT SPHERE, (X,Y,Z)
*  Name         Type    Definition
*  P1		REAL*8	FIRST POINT 
*  P2		REAL*8	SECOND POINT
*  P3		REAL*8	THIRD POINT
*  Output Parameters:
*  Name         Type    Definition
*  ARAD		REAL*8	ANGLE IN RADIANS:  -PI < ARAD <= +PI
*  			ANGLE WITH P2 AT THE VERTEX.  
*			AS YOU GO FROM P1 TO P2 AND THEN TO P3, IF YOU TURN
*			RIGHT, THEN ARAD IS POSITIVE;  LEFT:  ARAD IS NEGATIVE.
*			THINK OF P2 AS THE ORIGIN ON THE XY PLANE.  P1 IS OUT
*			ON THE X-AXIS IN THE POSITIVE X DIRECTION SOMEWHERE.
*			THE P2 P3 FORMS A RAY IN THE XY PLANE AT AN ANGLE
*			USUALLY KNOWN AS THETA.  THE POSITIVE AND NEGATIVE 
*			ANGLE CONVENTION COMES FROM THERE.  
*  Variables:
*  Locals :
*  Externals :
*  Modification History:                                            
*  Date			Revision	Author
*  $Date$ $Revision$ $Author$
*                                                                   
*********************************************************************/
      SUBROUTINE SANGL (ARAD,P1,P2,P3)
      character*100 SccsFileID
     -/'@(#)sangl.for	5.1 98/01/08 APS/ASF\0'/

      IMPLICIT NONE
      REAL*8 P1(3),P2(3),P3(3),U12(3),U23(3),XVEC(3), D, ARAD
      REAL*8 DOTMAG, XMAG, PI8
C---	FUNCTION FROM VECTOR LIBRARY
      REAL*8 DOT
      REAL*8 ABS
      INCLUDE 'APS_HOME:include/local/mps_const_math.inc'
      INTEGER I
      ARAD = 0
C---	PI IN REAL*8
      PI8 = PI
C---	SPECIAL CASES:  P1=P2 OR P1=-P2
C---			P3=P2 OR P3=-P2
C---	IN THESE CASES, WHERE ARAD IS UNDEFINED, BOMB THE PROGRAM.
      CALL GCLOSE(I,P1,P2)
      IF (I.NE.0) CALL GCBOMB
      CALL GCLOSE(I,P3,P2)
      IF (I.NE.0) CALL GCBOMB
C---	O.K.  THE VECTORS WEREN'T CLOSE AND THE OPPOSITES WEREN'T CLOSE.
C---	FIRST WE ARE JUST TRYING TO SEE IF THE ANGLE IS POSITIVE OR NEGATIVE, 
C---	WHICH SIDE OF THE P1 P2 PLANE THAT P3 IS ON.
C---	COMPUTE A UNIT NORM (U12) OF PLANE P1,P2, BY P2XP1.  THEN TAKE 
C---	P3 DOT U12 TO GET THE PROJECTION OF P3 ONTO THE NORM AND DIVIDE BY 
C---	THE MAGNITUDE OF P3(=1.0).  
C---	THE RESULT IS THE COSINE OF THE ANGLE 
C---	P3 MAKES WITH THE NORMAL U12.  U12 WAS COMPUTED SUCH THAT IF THE 
C---	COSINE IS POSITIVE, THEN ARAD WILL BE POSITIVE, AND IF THE COSINE 
C---	IS NEGATIVE, THEN ARAD IS NEGATIVE.  
C---	THEN COMPUTE A UNIT NORM (U23) OF P2,P3, BY P2XP3.  THEN COMPUTE THE 
C---	CROSS PRODUCT OF THE 2 UNIT NORMS, YIELDING THE SINE OF THE ANGLE 
C---	BETWEEN THE TWO PLANES(XMAG).  THIS IS THE ANGLE AT P2 ON THE SPHERE.  
C---	THE COSINE OF THE ANGLE BETWEEN THE TWO PLANES IS THE DOT PRODUCT.
C---	KNOWING THE SINE AND COSINE, THE ANGLE CAN BE DETERMINED BETWEEN 
C---	0 TO 180, WITH THE SIGN ALREADY DETERMINED.  
      CALL UCROSS(U12,P2,P1)
      D = DOT(P3,U12,3)
      IF (D.GE.0.0) I =  1
      IF (D.LT.0.0) I = -1
      CALL UCROSS(U23,P2,P3)
      CALL UCROSM(XMAG,XVEC,U12,U23)
      DOTMAG = DOT(U12,U23,3)
      ARAD = I * ABS( ATAN2(XMAG,DOTMAG) )
      IF (ARAD .LE. -PI8) ARAD = ARAD + 2.0D00 * PI8
      RETURN
      END
