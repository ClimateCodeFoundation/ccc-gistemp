c   To read the  oi SST version 2
c    and the masks.
c
c
      parameter (im=360,jm=180,ioff=im/2)
      parameter (iyrbeg=1880,iyreo=3000, nymax=15)
      parameter (monx=12*(iyreo-iyrbeg+1))
c
c   i is longitude index (values 1 to im)
c   j is the latitude index (values 1 to jm)
c
      real*4 sst(im,jm,12*nymax)
c
c   sst is the input SST real*2 array
c    (values are degrees C)
c
      integer info(8)
c
c   ls is the land/sea array not needed - built into climatology
c
      real to(monx)
      real clim(im,jm,12)
c
c   rsst is the reconstructed SST mean in degrees C
c   Read the EOF reconstruction area mask (unit 4)
c
      character*80 filei/
     *  'input/oiv2mon.yyyymm'/,line,title
      ny=index(filei,'yyyy')

      open(1,file='input/SBBX.HadR2',form='unformatted')
      read(1) info,title
      bad=info(7)
      mnow=info(1)
      info(1)=1  ! output file NOT trimmed/reorganized
      if(info(6).ne.iyrbeg) then
        write(*,*) 'iyrbeg-new',iyrbeg,' iyrbeg-old',info(6)
        stop 'iyrbeg inconsistent'
      end if
      write(*,*) 'SBBX.HadR2  opened'
      open(24,file='input/SBBX_LtSN.LnWE.dat',form='unformatted')

      open(14,file='input/oisstv2_mod4.clim',
     *        form='unformatted')
      read(14) line,clim
      close (14)
c
c   Read in the SST data for recent years (unit 11)
c
      call getarg(1,line)
      read(line,*) iyr1
      if(iyr1.lt.iyrbeg) then
        write(*,*) 'iyr1=',iyr1,' cannot be less than',iyrbeg
        stop 'arg 1 off'
      end if
      if(iyr1.gt.iyreo) then
        write(*,*) 'iyr1=',iyr1,' cannot be less than',iyreo
        stop 'arg 1 off'
      end if
      call getarg(2,line)
      read(line,*) mo1
      if(mo1.lt.1.or.mo1.gt.12) then
        write(*,*) 'mo1=',mo1,' must be between 1 and 12'
        stop 'arg 2 off'
      end if
      call getarg(3,line)
      read(line,*) mo2
      if(mo2.lt.mo1.or.mo2.gt.12) then
        write(*,*) 'mo2=',mo2,' must be between',mo1,' and 12'
        stop 'arg 3 off'
      end if
      info(4)=12*(iyr1-iyrbeg+1) ! wipe out data after new data
      info(5)=info(4)+4     ! old non-trimmed format

      do 100 iyr=iyr1,iyr1
      do 90 mon=mo1,mo2
      write(filei(ny:ny+5),'(I4.4,I2.2)') iyr,mon
      write(*,*) 'trying to read ',filei(1:70)
      open(11,form='unformatted',file=filei,err=110)
      read(11) iyr0,imo
      iyre=iyr
      moe=mon
      if(iyr.ne.iyr0) write(*,*) 'years not ok',iyr,iyr0
      if(mon.ne.imo) write(*,*) 'months not ok',mon,imo
c
c   iyr is year (value 1950 to ....)
c   imo is month (value 1 to 12)
c
      read(11) ((sst(i,j,12*(iyr-iyr1)+mon),i=1,im),j=1,jm)
   90 close (11)
  100 continue
  110 open(2,file='work/SBBX.HadR2.upd',form='unformatted')
      title(41:80)=' Had: 1880-11/1981, oi2: 12/1981-mm/yyyy'
      write(title(74:80),'(I2,''/'',I4)') moe,iyre
      write(2) info,title
c****
c**** Interpolate to Sergei's subbox grid
c****
      moff=(iyr1-iyrbeg)*12
      do 500 n=1,8000
      do 200 m=1,monx
  200 to(m)=bad
      call sread(1,to,mnow,lts,ltn,lnw,lne,mnext)
      mnow=mnext
      js=(18001+(lts+9000)*jm)/18000
      jn=(17999+(ltn+9000)*jm)/18000
      iw=(36001+(lnw+18000)*im)/36000+ioff
      ie=(35999+(lne+18000)*im)/36000+ioff
      if(ie.gt.im) then
         iw=iw-im
         ie=ie-im
      end if
      if(iw.gt.im) stop 'iw>im'
      if(ie.gt.im) stop 'ie>im'
      if(iw.lt.1) stop 'iw<1'
      if(ie.lt.1) stop 'ie<1'


      do 220 m=mo1,mo2
      month=mod(m-1,12)+1
      kt=0
      ts=0.
      do 210 j=js,jn
      do 210 i=iw,ie
      if(sst(i,j,m).le.-1.77.or.clim(i,j,month).eq.bad) go to 210
      kt=kt+1
      ts=ts+(sst(i,j,m)-clim(i,j,month))
c     write(*,*) 'new data:',m,' at',i,j,n
  210 continue
        to(m+moff)=bad
        if(kt.gt.0) to(m+moff)=ts/kt
  220 continue
      do 230 m=(iyre-iyrbeg)*12+moe+1,info(4)
c     write(*,*) 'filling with bad:',m
  230 to(m)=bad

      write(2) (to(i),i=1,info(4)),lts,ltn,lnw,lne
  500 continue
      stop
      end

      subroutine sread(in,a,na,lts,ltn,lnw,lne,mnext)
      real a(na)
      integer info(7)
      read(in) mnext,lts,ltn,lnw,lne,x,x,x,a
      return
      end
