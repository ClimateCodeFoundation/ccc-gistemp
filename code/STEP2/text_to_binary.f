      parameter (iyear1=1880, lim=20, multm=2, ibad=9999)

      integer info(9), monthly(12), iok(12)
      integer lat, lon, id(0:multm), ht,hto, mmax(0:multm)
      character title*80, name*36,nameo*36,line*64,li*102,sid*12,sidi*11
      integer,allocatable,dimension(:,:) :: idata

      call getarg(1,name)
      read(name,*) lastyr
      MTOT=12*(lastyr-iyear1+1)
      write(*,*) 'last year with data:',lastyr
      allocate( idata(mtot,0:multm) )

      title='GHCN V2 Temperatures (.1 C)'
      info(1) = 1
      info(2) = 1
      info(3) = 6
      info(4) = MTOT
      info(5) = MTOT+15
      info(6) = iyear1
      info(7) = 9999
      info(8) = -9999
      info(9) = MTOT

      open(1,file='work/Ts.txt',form = 'formatted')
      open(2,file='work/Ts.bin',form='unformatted')
      open(3,file='input/v2.inv',form='formatted')
      open(88,file='log/short.station.list',form='formatted')
      open(99,file='log/station.log',form='formatted')

      write(2) info, title

      ict = 0
      mult=0
      id1o=-99

      read(1, '(a64)') line   ; sidi='           '
10    continue
      ict = ict + 1
      if (1000*(ict/1000).eq.ict) print*,ict,'processed so far'
      read(line(2:64), '(i4,i5,a12,i4,a36)') lat, lon, sid, ht, name
!     write(0,*) sid
      do while (sid(1:11) .ne. sidi )
        read(3,'(a)',end=500) li
        sidi=li(1:11)
!     write(0,*) sidi
      end do
      name(31:31)=li(102:102) ! US-brightness index 1/2/3=dark/dim/brite
      name(32:32)=li(68:68)   ! population index (R/S/U=rural/other)
      name(33:33)=li(101:101) ! GHCN-brightness index A/B/C=dark/dim/brt
      name(34:36)=li(1:3)     ! country code (425=US)

C     print*, line(1:64)
      read(sid(4:12), '(i9)') idfull
      read(sid(4:11), '(i8)') id1

      if(id1.eq.id1o) then
        mult=mult+1
        if(mult.gt.multm) then
          write(*,*) 'increase multm from',multm
          stop 13
        end if
      else
        if(id1o.lt.0) go to 15
        do m=0,mult
        if(mmax(m).ge.lim) write(2) (idata(i,m),i=1,MTOT),lato,lono,
     *                               id(m),hto,nameo,1,MTOT
        if(mmax(m).lt.lim) write(88,*) id(m),' ',nameo,' dropped'
        end do
   15   mult=0
        id1o=id1
      end if
      id(mult)=idfull
      do i = 1, MTOT
      idata(i,mult) = ibad
      end do

      do m = 1,12
      iok(m) = 0
      end do
      mmax(mult)=0
      monmin=MTOT+1
      monmax=0

20    continue

      ieof = 1
      read(1, '(a64)', end=30) line

      if (line(1:1).eq.' ') then
      ieof = 0
      goto 30
      end if

      read(line, '(i4, 12i5)') iyr, (monthly(m), m = 1, 12)
C     print*, line

      ix = (iyr - iyear1)*12
      do m = 1, 12
      idatum = monthly(m)
      if (idatum.ne.9999) then
        idata(ix + m,mult) = idatum
        iok(m)=iok(m)+1
        if(mmax(mult).lt.iok(m)) mmax(mult)=iok(m)
        if(monmin.gt.ix+m) monmin=ix+m
        monmax=ix+m
      end if
      end do

      goto 20

30    continue

      write(99,*) mult, id(mult), mmax(mult),lim, monmin,monmax

      lato=lat
      lono=lon
      hto=ht
      nameo=name

      if (ieof.eq.0) goto 10

      do m=0,mult
      if(mmax(m).ge.lim) write(2) (idata(i,m),i=1,MTOT),lato,lono,
     *                             id(m),hto,nameo,1,MTOT
      if(mmax(m).lt.lim) write(88,*) id(m),' ',nameo,' dropped'
      end do
      close(1)
      close(2)
      write(*,*) 'number of station ids:',ict
      stop  0

  500 write(*,*) 'should not happen - inventory too short'
      write(*,*) sid
      stop 'should not happen'

      end
