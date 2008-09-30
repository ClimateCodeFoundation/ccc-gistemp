      parameter (iylast=3000,i4=(iylast-1701+1)*12)
      integer info(9),idata(i4),idato(i4)
      character header*105,title*80,name*36,nameo*36,cc*3,fn*12

      open(1,file='work/fort.1',form='formatted')
      read(1,'(a80)') header
      read(1,'(a3,i9,2f9.3,i5,9x,f9.3,27x,i5,1x,i4,i5,1x,i4,i5)')
     *   cc,IDc,sl1,sl2,knee,sl0,iy1,iy2,iy1e,iy2e,iflag
      do 200 iin=12,17
      io=iin+10
      fn='work/fort.xx'
      write(fn(11:12),'(i2)') iin
      open(iin,file=fn,form='unformatted')
      write(fn(11:12),'(i2)') io
      open(io,file=fn,form='unformatted')
      read(iin,end=300) info,title
      ibad=info(7)
      if(info(1).eq.ibad) go to 200
      m1o=info(1)
      m2o=info(9)
   10 read(iin,end=300) (idato(i),i=m1o,m2o),
     *   Lato,Lono,IDo,ihto,nameo,m1,m2
      if(IDo.ne.IDc) then
         if(nameo(31:32).eq.' R'.or.nameo(31:31).eq.'1') then
Calt     if(nameo(33:33).eq.'A') then
              write(io) info,title
         else            ! skip the station
              info(1)=m1
              info(9)=m2
              m1o=m1
              m2o=m2
              write(*,*) 'station',IDo,' ',nameo,' skipped'
              go to 10
         end if
      else
        call adj(info,idato,sl1,sl2,knee,sl0,iy1e,iy2e,iy1,iy2,iflag,
     *       m1o,m2o)
        write(*,*) 'station',IDo,' ',nameo,' adjusted'
        info(1)=m1o
        info(9)=m2o
        write(io) info,title
        read(1,'(a3,i9,2f9.3,i5,9x,f9.3,27x,i5,1x,i4,i5,1x,i4,2i5)',
     *     end=100) cc,IDc,sl1,sl2,knee,sl0,iy1,iy2,iy1e,iy2e,iflag
      end if

   30 read(iin,end=100) (idata(i),i=m1,m2),Lat,Lon,ID,iht,name,m1n,m2n
      if(ID.ne.IDc) then
          if(name(31:32).eq.' R'.or.name(31:31).eq.'1') then
Calt      if(name(33:33).eq.'A') then
            write(io)(idato(i),i=m1o,m2o),Lato,Lono,IDo,ihto,nameo,m1,m2
            write(*,*) 'station',IDo,' ',nameo,' saved',io-21
          else          ! skip the station
            m1=m1n
            m2=m2n
            write(*,*) 'station',ID,' ',name,' skipped'
            go to 30
          end if
      else
        call adj(info,idata,sl1,sl2,knee,sl0,iy1e,iy2e,iy1,iy2,iflag,
     *       m1,m2)
        write(*,*) 'station',ID,' ',name,' adjusted'
        write(io)(idato(i),i=m1o,m2o),Lato,Lono,IDo,ihto,nameo,m1,m2
        write(*,*) 'station',IDo,' ',nameo,' saved',io-21
        IDc=-9999
        read(1,'(a3,i9,2f9.3,i5,9x,f9.3,27x,i5,1x,i4,i5,1x,i4,2i5)',
     *     end=35)  cc,IDc,sl1,sl2,knee,sl0,iy1,iy2,iy1e,iy2e,iflag
      end if

   35 do 40 m=m1,m2
   40 idato(m)=idata(m)
      Lato=Lat
      Lono=Lon
      IDo=ID
      ihto=iht
      nameo=name
      m1o=m1
      m2o=m2
      m1=m1n
      m2=m2n
      go to 30

  100 write(io) (idato(m),m=m1o,m2o),Lato,Lono,IDo,ihto,nameo,ibad,ibad
      write(*,*) 'station',IDo,' ',nameo,' saved',io-21
  200 continue
      stop 0
  300 stop 1
      end

      subroutine adj(info,idata,sl1,sl2,knee,sl0,iy1,iy2,iy1a,iy2a,
     *       iflag,m1,m2)
      dimension info(*),idata(*)
      if(iflag.ne.0.and.iflag.ne.100) then
C****    Use linear approximation
         sl1=sl0
         sl2=sl0
      end if
      miss=info(7)
      m1o=m1
      m2o=m2
      m1=-100
      m0=12*(iy1-info(6))   !  Dec of year iy1
      do iy=iy1,iy2
        sl=sl1
        if(iy.gt.knee) sl=sl2
        iya=iy
        if(iy.lt.iy1a) iya=iy1a
        if(iy.gt.iy2a) iya=iy2a
        iadj=nint( (iya-knee)*sl-(iy2a-knee)*sl2 )
        do m=m0,m0+11
           if(m.ge.m1o.and.m.le.m2o.and.idata(m).ne.miss) then
               if(m1.lt.0) m1=m
               idata(m)=idata(m)+iadj
               m2=m
           end if
        end do
        m0=m0+12
      end do

      return
      end
