      integer,allocatable,dimension(:) :: idata
      integer info(9)
      character title*80,name*36,fileo*12

      open(1,file='work/Ts.bin',form='unformatted')
      read(1,end=100) info,title

      do n=51,56
      write(fileo,'(a11,i1)') 'work/Ts.bin',n-50
      open(n,file=fileo,form='unformatted')
      write(n) info,title
      end do
      allocate (idata(info(4)))
   10 read(1,end=100) idata,Lat,Lon,ID,iht,name,m1,m2
!!    write(99,*) Lat,Lon,ID
!!    if(ID.eq.603550002) write(98,'(12I6)') idata
      mu=lat+899
      if(mu.lt.0) mu=0
      mu=56-mu/300
      write(mu) idata,Lat,Lon,ID,iht,name,m1,m2
      go to 10

  100 stop 0
      end
