C--  Copyright (c)1996, California Institute of Technology.
C--  U.S. Government Sponsorship acknowledged.
C-- ==========================================================================
C--
C--  Fortran Filename:	rk78cn.for
C--
C--  Description:	
C--	
C--  Notes:
C--
C-- ==========================================================================
      SUBROUTINE RK78CN

      character*100 SccsFileID
     -/'@(#)rk78cn.for	5.1 98/01/08 APS/ASF\0'/

C  VERSION OF 4/1/85
C  PURPOSE
C    COMPUTES THE FEHLBERG COEFFICIENTS FOR A RUNGE-KUTTA 78 INTEGRATOR
C  INPUT
C    NONE
C  COMMON BLOCK OUTPUT
C    CH, ALPH, BETA, NORDER, ORDRCP, NTIMES
C  REFERENCES
C    JPL EM 312/85-153, 20 APRIL 1987
C    NASA TR R-287, OCTOBER 1968
C  ANALYSIS
C    J. H. KWOK - JPL
C  PROGRAMMER
C    J. H. KWOK - JPL
C  PROGRAM MODIFICATIONS
C    NONE
C  COMMENTS
C    THIS ROUTINE MUST BE CALLED BEFORE CALLING ROUTINE RK78
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      COMMON/FELCON/CH(13),ALPH(13),BETA(13,12),ORDRCP,NORDER,NTIMES
      DATA ZERO,ONE,TWO,THREE/0.D0,1.D0,2.D0,3.D0/
      NORDER=8
      ORDRCP=ONE/NORDER
      NTIMES=13 
      DO 30 I=1,13  
      DO 40 J=1,12  
   40 BETA(I,J)=ZERO
      ALPH(I)=ZERO
   30 CH(I)=ZERO
      CH(6)=34.D0/105.D0
      CH(7)=9.D0/35.D0  
      CH(8)=CH(7)   
      CH(9)=9.D0/280.D0 
      CH(10)=CH(9)  
      CH(12)=41.D0/840.D0   
      CH(13)=CH(12) 
      ALPH(2)=TWO/27.D0
      ALPH(3)=ONE/9.D0 
      ALPH(4)=ONE/6.D0 
      ALPH(5)=5.D0/12.D0
      ALPH(6)=.5D0  
      ALPH(7)=5.D0/6.D0 
      ALPH(8)=ONE/6.D0 
      ALPH(9)=TWO/THREE 
      ALPH(10)=ONE/THREE
      ALPH(11)=ONE
      ALPH(13)=ONE 
      BETA(2,1)=TWO/27.D0  
      BETA(3,1)=ONE/36.D0  
      BETA(4,1)=ONE/24.D0  
      BETA(5,1)=5.D0/12.D0  
      BETA(6,1)=.5D-1   
      BETA(7,1)=-25.D0/108.D0   
      BETA(8,1)=31.D0/300.D0
      BETA(9,1)=TWO
      BETA(10,1)=-91.D0/108.D0  
      BETA(11,1)=2383.D0/4100.D0
      BETA(12,1)=THREE/205.D0
      BETA(13,1)=-1777.D0/4100.D0   
      BETA(3,2)=ONE/12.D0  
      BETA(4,3)=ONE/8.D0   
      BETA(5,3)=-25.D0/16.D0
      BETA(5,4)=-BETA(5,3)  
      BETA(6,4)=.25D0   
      BETA(7,4)=125.D0/108.D0   
      BETA(9,4)=-53.D0/6.D0 
      BETA(10,4)=23.D0/108.D0   
      BETA(11,4)=-341.D0/164.D0 
      BETA(13,4)=BETA(11,4) 
      BETA(6,5)=.2D0
      BETA(7,5)=-65.D0/27.D0
      BETA(8,5)=61.D0/225.D0
      BETA(9,5)=704.D0/45.D0
      BETA(10,5)=-976.D0/135.D0 
      BETA(11,5)=4496.D0/1025.D0
      BETA(13,5)=BETA(11,5) 
      BETA(7,6)=125.D0/54.D0
      BETA(8,6)=-TWO/9.D0  
      BETA(9,6)=-107.D0/9.D0
      BETA(10,6)=311.D0/54.D0   
      BETA(11,6)=-301.D0/82.D0  
      BETA(12,6)=-6.D0/41.D0
      BETA(13,6)=-289.D0/82.D0  
      BETA(8,7)=13.D0/900.D0
      BETA(9,7)=67.D0/90.D0 
      BETA(10,7)=-19.D0/60.D0   
      BETA(11,7)=2133.D0/4100.D0
      BETA(12,7)=-THREE/205.D0   
      BETA(13,7)=2193.D0/4100.D0
      BETA(9,8)=THREE
      BETA(10,8)=17.D0/6.D0 
      BETA(11,8)=45.D0/82.D0
      BETA(12,8)=-THREE/41.D0
      BETA(13,8)=51.D0/82.D0
      BETA(10,9)=-ONE/12.D0
      BETA(11,9)=45.D0/164.D0   
      BETA(12,9)=THREE/41.D0 
      BETA(13,9)=33.D0/164.D0   
      BETA(11,10)=18.D0/41.D0   
      BETA(12,10)=6.D0/41.D0
      BETA(13,10)=12.D0/41.D0   
      BETA(13,12)=ONE  
      RETURN
      END
