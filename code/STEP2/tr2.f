      SUBROUTINE TREND2(xc,A,LEN,Xmid,BAD,MIN1,MIN2, SL1,SL2,Ymid,RMS,
     *                                               SL,Y0,RMS0)
C**** finds a fit using regression analysis by a line
C**** with a break in slope at Xmid. Returned are the 2 slopes
C**** SL1,SL2 provided we have at least MIN1,MIN2 data.
C**** Linear regression data are also computed (for emergencies)
      REAL*8 xc(*),A(*)
      REAL*8 sx(2),sxx(2),sxa(2),sa,saa,denom,xnum1,xnum2
      INTEGER kount(2)

      sl1=bad
      sl2=bad
      Ymid=bad
      sa=0.
      saa=0.
      do k=1,2
        kount(k)=0
        sx(k)=0.
        sxx(k)=0.
        sxa(k)=0.
      end do

      do 100 n=1,len
      if(a(n).eq.BAD) go to 100
      x=xc(n)-Xmid
      sa=sa+a(n)
      saa=saa+a(n)**2
      k=1
      if(x.gt.0.) k=2
      kount(k)=kount(k)+1
      sx(k)=sx(k)+x
      sxx(k)=sxx(k)+x**2
      sxa(k)=sxa(k)+x*a(n)
  100 continue

      ntot=kount(1)+kount(2)
      denom=ntot*sxx(1)*sxx(2)-sxx(1)*sx(2)**2-sxx(2)*sx(1)**2
      xnum1=sx(1)*(sx(2)*sxa(2)-sxx(2)*sa)+sxa(1)*(ntot*sxx(2)-sx(2)**2)
      xnum2=sx(2)*(sx(1)*sxa(1)-sxx(1)*sa)+sxa(2)*(ntot*sxx(1)-sx(1)**2)

      if(kount(1).lt.MIN1.or.kount(2).lt.MIN2) return
      sl1=xnum1/denom
      sl2=xnum2/denom
      Ymid=(sa-sl1*sx(1)-sl2*sx(2))/ntot
      RMS=ntot*Ymid**2+saa-2*Ymid*(sa-sl1*sx(1)-sl2*sx(2))+
     *    sl1*sl1*sxx(1)+sl2*sl2*sxx(2)-2*sl1*sxa(1)-2*sl2*sxa(2)

C**** linear regression
      sx(1)=sx(1)+sx(2)
      sxx(1)=sxx(1)+sxx(2)
      sxa(1)=sxa(1)+sxa(2)
      sl=(ntot*sxa(1)-sa*sx(1))/(ntot*sxx(1)-sx(1)**2)
      Y0=(sa-sl*sx(1))/ntot
      RMS0=ntot*Y0**2+saa+sl*sl*sxx(1)-2*Y0*(sa-sl*sx(1))-2*sl*sxa(1)

      return
      end
