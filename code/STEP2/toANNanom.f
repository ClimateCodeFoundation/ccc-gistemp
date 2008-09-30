      parameter (iylast=3000,i4o=(iylast-1701+1),i4=i4o*12)
      integer info(9),idata(0:i4),iann(i4o),llih(4)
      integer infoo(9),ianno(i4o),lliho(4)
      character title*80,name*36,nameo*36,titleo*80/
     *    'ANNUAL MEAN TEMPERATURE ANOMALIES (.01 C)'/

      open(2,file='work/fort.2',form='unformatted')
      open(3,file='work/fort.3',form='unformatted')
      read(2,end=300) info,title
      write(*,*) title,info
      ibad=info(7)
      if(info(1).eq.ibad) stop 'no stations'
      do i=2,8
      infoo(i)=info(i)
      end do
      infoo(3)=5   !  ann.means (6=mon.means)
      infoo(4)=info(4)/12   ! length of time series
      infoo(5)=infoo(4)+info(5)-info(4) ! length of records

      m1=info(1)
      m2=info(9)
   10 do i=0,i4
      idata(i)=ibad
      end do
Cw           write(*,*) 'reading 2',1+m2-m1
      call sread(idata(m1),1+m2-m1,lliho,nameo,m1n,m2n)

      iy1=1+(m1-1)/12
      mon1=12*(iy1-1)   !  = december of year iy1-1
      nyrs=(m2+11-mon1)/12
      if(nyrs-1+iy1.gt.i4o) nyrs=i4o+1-iy1
      if(nyrs.le.0) stop 'no station records - impossible'
Cw    write(*,*) 'm1,m2,mon1,nyrs,iy1',m1,m2,mon1,nyrs,iy1
      call annav(idata(mon1),nyrs, ianno, iy1, i1,i2, ibad)
Cw             write(*,*) 'after annav',i1,i2
      m1=m1n
      m2=m2n
      if(i1.eq.ibad) go to 10 ! stop 11  ! go to 10
      infoo(1)=i1
      infoo(9)=i2
      write(3) infoo,titleo

   20 do i=0,i4
      idata(i)=ibad
      end do
      call sread(idata(m1),1+m2-m1,llih,name,m1n,m2n)

      iy1=1+(m1-1)/12
      mon1=12*(iy1-1)   !  = december of year iy1-1
      nyrs=(m2+11-mon1)/12
      if(nyrs-1+iy1.gt.i4o) nyrs=i4o+1-iy1
      if(nyrs.le.0) stop 'no station records - impossible'
      call annav(idata(mon1),nyrs, iann, iy1, i1n,i2n, ibad)

      if(i1n.ne.ibad) then
          call swrite (ianno(i1),1-i1+i2,lliho,nameo,i1n,i2n)
      else
          m1=m1n
          m2=m2n
          go to 20
      end if

      do 40 m=i1n,i2n
   40 ianno(m)=iann(m)
      do n=1,4
      lliho(n)=llih(n)
      end do
      nameo=name

      i1=i1n
      i2=i2n
      m1=m1n
      m2=m2n
      if(m1.ne.ibad) go to 20

  100 call swrite (ianno(i1),1-i1+i2,lliho,nameo,ibad,ibad)
      stop 0
  300 stop 1
      end

      subroutine sread (idata,ndim,llih,name,m1,m2)
      character*36 name
      dimension idata(ndim),llih(4)
      read(2) idata,llih,name,m1,m2
      return
      end

      subroutine swrite (idata,ndim,llih,name,m1,m2)
      character*36 name
      dimension idata(ndim),llih(4)
      write(3) idata,llih,name,m1,m2
C        do i=1,ndim
C        write(*,*) i,idata(i)
C        end do
C        if(ndim.gt.0) stop 11
      return
      end

      subroutine annav (mon,nyrs, iann, iy1, iy1n,iy2n, ibad)
      dimension mon(12,nyrs),iann(nyrs),seas(4),av(12)
      do m=1,12
      ny=0
      av(m)=0.
      do n=1,nyrs
      if(mon(m,n).ne.ibad) then
         ny=ny+1
         av(m)=av(m)+mon(m,n)
      end if
      end do
      if(ny.eq.0) stop 'station too short - impossible'
      av(m)=av(m)/ny
      end do
Cw                write(*,*) av
      iy1n=ibad
      iy2n=ibad
      nok=0
      do n=1,nyrs
C**** find 4 seasonal means ; mon(1,.)=dec ... mon(12,.)=nov-data
      do is=1,4
      seasis=0.
      nms=0
      do m=1+(is-1)*3,3+(is-1)*3
      if(mon(m,n).ne.ibad) then
        nms=nms+1
        seasis=seasis+(mon(m,n)-av(m))
      end if
      end do
      seas(is)=ibad
      if(nms.gt.1) seas(is)=seasis/nms
      end do
C**** find annual mean anomalies from seasonal means
      sann=0.
      nss=0
      do is=1,4
      if(seas(is).ne.ibad) then
        sann=sann+seas(is)
        nss=nss+1
      end if
      end do
      iann(iy1+n-1)=ibad
      if(nss.gt.2) then
        iy2n=iy1+n-1
        iann(iy2n)=nint(10.*sann/nss)
        if(nok.eq.0) iy1n=iy2n
        nok=nok+1
      end if
Cw              write(*,*) iy1+n-1, iann(iy1+n-1)
      end do
        !    iy1n=ibad ! to stop run
      return
      end
