      SUBROUTINE DITSCT3 (AX,RCENTR,RHO,RZERO,R1,R2,N)
C/*   SUBROUTINE DITSCT3 (AX,RCENTR,RHO,RZERO,R1,R2,N) ---------
C
C              ***** DOUBLE PRECISION SUBROUTINE *****
C     DERIVED FROM ITSCT2
C     GIVEN AX, THE THREE SEMI-AXES OF A CENTERED ELLIPSOID (AXES
C     COINCIDING WITH CARTESIAN FRAME), A CIRCLE CENTERED AT 
C     RCENTR WITH UNIT VECTOR RHO BEING THE NORMAL TO ITS PLANE
C     AND RZERO BEING A RADIUS VECTOR, THE SUBROUTINE COMPUTES
C     THE CIRCLE-ELLIPSOID INTERSECTION POINTS R1,R2.
C     THE SUBROUTINE ALSO OUTPUTS (IN LABELLED COMMON/CTSCT2/) 
C     THE CORRESPONDING ANGLES PHI1,PHI2 (IN DEGREES). PHI IS 
C     DEFINED AS THE ANGLE THROUGH WHICH RZERO HAS TO BE ROTATED 
C     (ABOUT RHO) TO BRING IT TO AN INTERSECTION POINT.
C     R1,R2 ARE CHOSEN SUCH THAT ABS(PHI1) .LE. ABS(PHI2)
C     THE SUBROUTINE IS VALID ONLY WHEN ABS(PHI1) IS SMALL
C     N IS AN OUTPUT CODE DESCRIBING THE SITUATION AS FOLLOWS:
C     N=1  TWO INTERSECTIONS
C     N=2  ONE INTERSECTION (TANGENCY)
C     N=3  NO INTERSECTIONS
C*/
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AX(3),RCENTR(3),RHO(3),RZERO(3),G(3),R1(3)
      DIMENSION RDIF(3),R2(3),DM1(3),RSUM(3)
      COMMON/CDTSCT2/PHI1,PHI2
      DATA DGPRAD/57.2957795130823208768D0/
      character*80 sccsid
      data sccsid /'@(#)ditsct3.f	1.3 96/04/09 22:51:43\0'/
      CALL DSCALAR (3,RHO,RCENTR,H)
      CALL DVECTP (RHO,RZERO,G)
      CALL DSUMV (RSUM,RCENTR,RZERO,3,2)
      CALL DSUMV (RDIF,RCENTR,RZERO,3,1)
      A=-1.D0
      B=0.D0
      CQ=-1.D0
      DO 3 I=1,3
      A=A+(RSUM(I)*RDIF(I)+2.D0*G(I)**2)/AX(I)**2
      B=B+(G(I)*RSUM(I))/AX(I)**2
    3 CQ=CQ+(RSUM(I)/AX(I))**2
      A=.5D0*A
      B=2.D0*B
      CALL DQUAD (A,B,CQ,NQ,PSI1,PSI1I,PSI2,PSI2I)
      N=3
      IF (NQ.GT.2) RETURN
      N=NQ
      D=PSI1+PSI2
      IF (D.GE.0.D0) GO TO 2
      D=PSI1
      PSI1=PSI2
      PSI2=D
    2 CALL DROTV (2,PSI1,D,D1,RHO,RZERO,DM1)
      CALL DSUMV2 (R1,1.D0,RCENTR,1.D0,DM1,3)
      CALL DROTV (2,PSI2,D,D1,RHO,RZERO,DM1)
      CALL DSUMV2 (R2,1.D0,RCENTR,1.D0,DM1,3)
      PHI1=DGPRAD*PSI1
      PHI2=DGPRAD*PSI2
      RETURN
      END
