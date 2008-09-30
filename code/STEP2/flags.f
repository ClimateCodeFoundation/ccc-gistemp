      parameter (lshort=7,slplim=1.,slpx=.5)
      character*105 title,cc*3
      sumch=0.
      sumchp=0.
      nsta=0
      nstap=0
      nyrs=0
      nyrsp=0
      nshort=0
      nsl1=0
      nsl2=0
      ndsl=0
      nok=0
      nokx=0
      nswitch=0
      nsw0=0
      open(78,file='work/fort.78',form='formatted')
      open( 2,file='work/fort.2',form='formatted')
   10 read(78,'(a)',end=100) title
      read(title,'(a3,i9,2f9.3,i5,5f9.3,i5,1x,i4,i5,1x,i4)',end=100)
     *    cc,id,sl1,sl2,knee,yk,sl,ylin,rms,rms0,iy1,iy2,iy1e,iy2e
      nsta=nsta+1
      sumch=sumch+sl*(iy2-iy1+1)
      nyrs=nyrs+(iy2-iy1+1)
      if(sl.lt.0.) then
        nstap=nstap+1
        nyrsp=nyrsp+(iy2-iy1+1)
        sumchp=sumchp+sl*(iy2-iy1+1)
      end if
C***  classify : iflag: +1 for short legs etc
      iflag=0
      if(knee.lt.iy1+lshort.or.knee.gt.iy2-lshort) iflag=iflag+1
      if(knee.lt.iy1+lshort.or.knee.gt.iy2-lshort) nshort=nshort+1
      if(abs(sl1).gt.slplim) iflag=iflag+20
      if(abs(sl1).gt.slplim) nsl1=nsl1+1
      if(abs(sl2).gt.slplim) iflag=iflag+10
      if(abs(sl2).gt.slplim) nsl2=nsl2+1
      if(abs(sl2-sl1).gt.slplim) iflag=iflag+100
      if(abs(sl2-sl1).gt.slplim) ndsl=ndsl+1
      if(abs(sl2-sl1).gt.slpx) iflag=iflag+100
      if(iflag.eq.0) nok=nok+1
      if(iflag.eq.100) nokx=nokx+1
      if(sl1*sl2.lt.0..and.abs(sl1).gt..2.and.
     *  abs(sl2).gt..2) iflag=iflag+1000
      if(iflag.ge.1000) nswitch=nswitch+1
      if(iflag.eq.1000) nsw0=nsw0+1
      write(title(101:105),'(i5)') iflag
      write(2,'(a)') title
      go to 10

  100 write (*,*) 'all',nsta,-sumch/nsta,-10.*sumch/nyrs
      write (*,*) 'urb warm',nstap,-sumchp/nstap,-10.*sumchp/nyrsp
      write(*,*) '# short,sl1,sl2,dsl,ok',nshort,nsl1,nsl2,ndsl,nok,nokx
      write (*,*) 'switches: all , else ok ',nswitch,nsw0
      stop
      end
