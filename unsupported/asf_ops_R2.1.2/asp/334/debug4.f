C Alaska SAR Processor (ASP) %W% %E% %U%
      SUBROUTINE DEBUG4(NP,L,J,V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,
     *                         V11,V12)
C/*   SUBROUTINE DEBUG4(NP,L,J,V1,V2,V3,V4,V5,V6,V7,V8,V9,V10, -----
C    *                         V11,V12)
C
C     USED FOR DEBUGGING
C     PRINTS NP (LOCATION IDENTIFIER) AND FIRST J OF ARGUMENTS 
C     FOLLOWING J. 
C     J.GT.0 (BUT SEE ENTRY DEBUG44(NP) )
C     L=0  HEX PRINTOUT OF 16-BIT WORDS (J.LE.12)
C          (EXACT MACHINE REPRESENTATION)
C     L=1  HEX PRINTOUT OF 32-BIT WORDS (J.LE.4)  
C          (EXACT MACHINE REPRESENTATION)
C     L=2  INTEGER PRINTOUT (J.LE.4)  
C     L=3  FLOATING POINT PRINTOUT (SINGLE PRECISION WORDS; 
C          8-DIGIT MANTISSAS; J.LE.4)   
C     L=4  FLOATING POINT PRINTOUT (DOUBLE PRECISION WORDS; 
C          8-DIGIT MANTISSAS; J.LE.4)  
C     L=5  FLOATING POINT PRINTOUT (DOUBLE PRECISION WORDS; 
C          17-DIGIT MANTISSAS; J.LE.2)
C     L=6  COMPLEX*8 PRINTOUT (3-DIGIT MANTISSAS; J.LE.3)
C     L=7  COMPLEX*8 PRINTOUT (8-DIGIT MANTISSAS; J.LE.2)
C     L=8  CHARACTER PRINTOUT (FOUR CHARACTERS PER WORD; 
C          J.LE.12)
C     L=9  HEX PRINTOUT OF DOUBLE-PRECISION WORDS (J.LE.2)
C          (EXACT MACHINE REPRESENTATION)
C     ENTRY DEBUG44(NP)
C     PRINTS NP (LOCATION IDENTIFIER) ONLY
C     NP.LT.1000
C*/
  101 FORMAT(' DEBUG4(',I1,') POINT',I4.3,':',4(7X,A8))   
  102 FORMAT(' DEBUG4(',I1,') POINT',I4.3,':',4I15)   
  103 FORMAT(' DEBUG4(',I1,') POINT',I4.3,':',1P4E15.7)  
  104 FORMAT(' DEBUG4(',I1,') POINT',I4.3,':',1P4D15.7)
  105 FORMAT(' DEBUG4(',I1,') POINT',I4.3,':',1P2D30.16)
  106 FORMAT(' DEBUG4(',I1,') POINT',I4.3,':',1P6E10.2)
  108 FORMAT(' DEBUG4(',I1,') POINT',I4.3,':',12(1X,A4))
  109 FORMAT(' DEBUG4(',I1,') POINT',I4.3,':',2(14X,A16))
  114 FORMAT(' DEBUG44   POINT',I4.3,':')
      INTEGER IARG(12),JL(0:9)
      INTEGER*2 MSHORT(2,12)
      REAL V1(1),V2(1),V3(1),V4(1),V5(1),V6(1)
      REAL V7(1),V8(1),V9(1),V10(1),V11(1),V12(1)
      REAL ARG(12),V(12),VV(2,4)   
      DOUBLE PRECISION DARG(4)
      COMPLEX CARG(3)
      CHARACTER C(12)*4,N(4)*9,D(2)*17,M(12)*5
      EQUIVALENCE (IARG,ARG,V,C,MSHORT),(DARG,CARG,VV)
      DATA JL/12,4*4,2,3,2,12,2/
      GO TO (3,3,3,3,1,1,1,1,3,1),L+1
    3 V(1)=V1(1)
      IF (J.LT.2) GO TO 9
      V(2)=V2(1)
      IF (J.LT.3) GO TO 9
      V(3)=V3(1)
      IF (J.LT.4) GO TO 9
      V(4)=V4(1)
      IF (J.LT.5) GO TO 9
      V(5)=V5(1)
      IF (J.LT.6) GO TO 9
      V(6)=V6(1)
      IF (J.LT.7) GO TO 9
      V(7)=V7(1)
      IF (J.LT.8) GO TO 9
      V(8)=V8(1)
      IF (J.LT.9) GO TO 9
      V(9)=V9(1)
      IF (J.LT.10) GO TO 9
      V(10)=V10(1)
      IF (J.LT.11) GO TO 9
      V(11)=V11(1)
      IF (J.LT.12) GO TO 9
      V(12)=V12(1)
    9 IF (L.EQ.1) THEN
         DO 4 I=1,J
    4    CALL ITOHEX (IARG(I),N(I))
      ELSE IF (L.EQ.0) THEN
	      DO 6 I=1,J
    6         CALL SITOHEX (MSHORT(1,I),M(I))
      END IF
      GO TO 2
    1 VV(1,1)=V1(1)
      VV(2,1)=V1(2)
      IF (J.LT.2) GO TO 8
      VV(1,2)=V2(1)
      VV(2,2)=V2(2)
      IF (J.LT.3) GO TO 8
      VV(1,3)=V3(1)
      VV(2,3)=V3(2)
      IF (J.LT.4) GO TO 8
      VV(1,4)=V4(1)
      VV(2,4)=V4(2)
    8 IF (L.NE.9) GO TO 2
      DO 5 I=1,J
    5 CALL DTOHEX (DARG(I),D(I))
    2 JP=J
      IF (J.GT.JL(L)) JP=JL(L)
      GO TO (20,21,22,23,24,25,26,27,28,29),L+1  
   20 WRITE(6,108) L,NP,(M(K),K=1,JP)   
      RETURN  
   21 WRITE(6,101) L,NP,(N(K),K=1,JP)   
      RETURN  
   22 WRITE(6,102) L,NP,(IARG(K),K=1,JP)   
      RETURN  
   23 WRITE(6,103) L,NP,(ARG(K),K=1,JP)   
      RETURN  
   24 WRITE(6,104) L,NP,(DARG(K),K=1,JP)   
      RETURN  
   25 WRITE(6,105) L,NP,(DARG(K),K=1,JP)   
      RETURN  
   26 WRITE(6,106) L,NP,(CARG(K),K=1,JP)   
      RETURN  
   27 WRITE(6,104) L,NP,(CARG(K),K=1,JP)   
      RETURN  
   28 WRITE(6,108) L,NP,(C(K),K=1,JP)   
      RETURN  
   29 WRITE(6,109) L,NP,(D(K),K=1,JP)   
      RETURN 
      ENTRY DEBUG44(NP) 
      WRITE (6,114) NP  
      RETURN  
      END   
