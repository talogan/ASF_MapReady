      SUBROUTINE INTRP81 (*,Y,F,YT,N,IFIRST)
C/*   SUBROUTINE INTRP81 (*,Y,F,YT,N,IFIRST)  --------------
C
C     INTERPOLTES LINEARLY TO FIND Y(F*DX) FROM THE N-TERM TABLE YT 
C     ORGANIZED AS FOLLOWS:
C         YT(1)=Y(IFIRST*DX)
C         YT(N)=Y((N-1+IFIRST)*DX)
C     RETURN1:  INPUT ARGUMENT (F*DX) FALLS OUTSIDE TABLE
C*/
      DIMENSION YT(1)
      character*80 sccsid
      data sccsid /'@(#)intrp81.f	1.3 96/04/09 22:51:51\0'/
      NF=INT(F)
      IF (F) 1,2,2
    1 NF=NF-1
    2 NFT=NF-IFIRST+1
      NFT1=NFT+1
      IF (NFT.LT.1 .OR. NFT1.GT.N) RETURN1
      P=F-FLOAT(NF)
      PZ =1.-P
      Y=PZ*YT(NFT)+P*YT(NFT1)
      RETURN
      END
