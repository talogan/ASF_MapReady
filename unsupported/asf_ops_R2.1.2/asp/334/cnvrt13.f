C Alaska SAR Processor (ASP) %W% %E% %U%
      SUBROUTINE CNVRT13 (P,JP,NM1,SCLNG)
C/*   SUBROUTINE CNVRT13 (P,JP,NM1,SCLNG)  ----------------
C
C     WE ARE GIVEN A REAL ARRAY (P) OF NON-NEGATIVE NUMBERS.
C     THE SUBROUTINE PLACES A FIXED-POINT SCALED VERSION OF 
C     THIS ARRAY IN INTEGER*2 ARRAY JP.
C     THE NUMBER REPRESENTATION IN ARRAY JP IS "UNSIGNED 
C     MAGNITUDE" IN WHICH THE MOST SIGNIFICANT BIT REPRESENTS 
C     2**(-1).
C     THE SCALING IS SUCH THAT THE MAXIMAL VALUE OF P (PMAX)
C     SCALES TO (1-2**(-16)=65535./65536. IN ARRAY JP.
C     NOTE: FORTRAN INTERPRETS THE SCALED PMAX AS A NEGATIVE 
C     NUMBER.
C     OUTPUT VARIABLE SCLNG IS A SCALING PARAMETER SATISFYING
C          PMAX = SCLNG*(65535./65536.)
C	   SCLNG IS ALSO AN INPUT TO COMPUTE PMAX.
C     BOTH ARRAYS ARE DIMENSIONED (0:NM1)
C*/
      INTEGER*2 JP(0:NM1)
      INTEGER TWO15,TWO16
      REAL P(0:NM1)
      DATA TWO15/32768/ ! =(2**15)
      DATA TWO16/65536/ ! =(2**16)
      DATA FLIM/65535./ ! =(2**16)-1.
      CALL MAX5 (P,0,NM1,PMAX)
C Ming, 11/18/92
      IF (PMAX .LT. SCLNG) THEN
	PMAX = SCLNG*(65535./65536.)
      ELSE
	WRITE(6,*) 'WARMING, PMAX = ', PMAX, ' TOO MUCH'
      END IF		
      FMLT=FLIM/PMAX
      DO 1 I=0,NM1
      MP=NINT(FMLT*P(I))
      IF (MP.GE.TWO15) MP=MP-TWO16
    1 JP(I)=MP
      SCLNG=FLOAT(TWO16)/FMLT
      RETURN
      END
