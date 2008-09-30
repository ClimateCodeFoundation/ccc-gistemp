C*********************************************************************
C ***               unit#
C *** Input files:  31-36 ANN.dTs.GHCN.CL.1 ... ANN.dTs.GHCN.CL.6
C ***
C *** Output file:  78    list of ID's of Urban stations with
C ***                     homogenization info  (text file)
C ***                     Header line is added in subsequent step
C*********************************************************************
C****
C**** This program combines for each urban station the rural stations
C**** within R=1000km and writes out parameters for broken line
C**** approximations to the difference of urban and combined rural
C**** annual anomaly time series.
C****
C**** Input files:    units 31,32,...,36
C****  Record 1: I1,INFO(2),...,INFO(8),I1L,TITLE   header record
C****  Record 2: IDATA(I1-->I1L),LT,LN,ID,HT,NAME,I2,I2L  station 1
C****  Record 3: IDATA(I2-->I2L),LT,LN,ID,HT,NAME,I3,I3L  station 2
C****            etc.       NAME(31:31)=brightnessIndex 1=dark->3=bright
C****            etc.       NAME(32:32)=pop.flag R/S/U rur/sm.town/urban
C****            etc.       NAME(34:36)=country code
C****
C****  IDATA(1) refers to year IYRBEG=INFO(6)
C****  IYRM=INFO(4) is the max. length of an input time series,
C****
C****  INFO(1),...,INFO(8)  are 4-byte integers,
C****  TITLE is an 80-byte character string,
C****      INFO 2 = KQ (quantity flag, see below)
C****           3 = MAVG (time avg flag: 1 - 4 DJF - SON, 5 ANN,
C****                     6 MONTHLY, 7 SEAS, 8 - 19 JAN - DEC  )
C****           4 = IYRM  (length of each time record)
C****           5 = IYRM+4 (size of data record length)
C****           6 = IYRBEG (first year of full time record)
C****           7 = flag for missing data
C****           8 = flag for precipitation trace

C**** Combining of rural stations
C**** ===========================
C**** Stations within Rngbr km of the urban center U contribute
C**** to the mean at U with weight  1.- d/Rngbr  (d = distance
C**** between rural and urban station in  km). To remove the station
C**** bias, station data are shifted before combining them with the
C**** current mean. The shift is such that the means over the time
C**** period they have in common remains unchanged. If that common
C**** period is less than 20(NCRIT) years, the station is disregarded.
C**** To decrease that chance, stations are combined successively in
C**** order of the length of their time record.

C**** The homogeneity adjustment parameters
C**** =====================================
C**** To minimize the impact of the natural local variability, only
C**** that part of the combined rural record is actually used that is
C**** supported by at least 3 stations, i.e. heads and tails of the
C**** record that are based on only 1 or 2 stations are dropped. The
C**** difference between that truncated combination and the non-rural
C**** record is found and the best linear fit and best fit by a broken
C**** line (with a variable "knee") to that difference series are found.
C**** The parameters defining those 2 approximations are tabulated.
C****
C**** Note: No attempt is made to find the longterm trends for urban
C****       and rural combination separately; using the difference only
C****       minimizes the impact of short term regional events that
C****       affect both rural and urban stations, hence cancel out.

C?*** Input parameters (# of input files, time period)
      PARAMETER (INM=6,KM0=1,IYRM0=KM0*(3000-1880+1),ITRIM=1)
C?*** Earth radius, lower overlap limit, min.coverage of approx.range
      PARAMETER (REARTH=6375.,NCRIT=20,NRURM=3,XCRIT=2./3.)
C?*** Work array sizes
      PARAMETER (NSTAM=8000,MSIZE=700000)
      CHARACTER*80 TITLEI,TITLEO,NAME*36,fn*12
C**** Input arrays
      DIMENSION IDATA(MSIZE),INFO(8+ITRIM),ITRL(14+ITRIM)
      REAL RDATA(MSIZE)
      EQUIVALENCE (ITRL(5),NAME),(IDATA,RDATA)
C**** Station dependent arrays
      DIMENSION I1SU(NSTAM),ILSU(NSTAM),MFSU(NSTAM),IDU(NSTAM),
     *          I1SR(NSTAM),ILSR(NSTAM),MFSR(NSTAM),IDR(NSTAM),
     *   LEN(NSTAM),WTI(NSTAM),IORD(NSTAM),ISOFI(NSTAM) ! rural
     *   ,nuseid(nstam),lenis(nstam)   ! for diagn. purposes only
      real*8 r(nstam)
      REAL*4 CSLONU(NSTAM),CSLATU(NSTAM),SNLONU(NSTAM),SNLATU(NSTAM)
      REAL*4 CSLONR(NSTAM),CSLATR(NSTAM),SNLONR(NSTAM),SNLATR(NSTAM)
      CHARACTER*3 CC(NSTAM)
C**** Output-related arrays (dim of output array = IYRM < IYRM0)
      DIMENSION AVG(IYRM0),DNEW(IYRM0),WT(IYRM0),URB(IYRM0),IWT(IYRM0)
      COMMON/LIMIT/XBAD,NLAP
      COMMON/FITCOM/W(900),X(900),F(900),CF(900),DF(900),ZFP(20),ZR(20)
     1             ,FPAR(20),DELTAP,DFSTOP,X0,TMEAN,RMEAN,TFMEAN,RMSFIT
     2             ,YR(900),TS(900),TSFIT(900),RMSP(20),KPFIT(20)
     3             ,KPAR,NXY,NCP,NWCYCL,NITERS,NFIBNO,LISTIT,IDW
      REAL*8 W,X,F,CF,DF,ZFP,ZR,FPAR,DELTAP,DFSTOP,X0,TMEAN,RMEAN,TM3
      REAL*8 TFMEAN,RMSFIT,YR,TS,TSFIT,RMSP
      PI180=DATAN(1.D0)/45.D0
      RngbrF=1000.
      IF(IARGC().GT.0) THEN
        CALL GETARG(1,TITLEI)
        READ(TITLEI,*) IRngbr
        RngbrF=IRngbr
      END IF
      RngbrH=RngbrF/2
      RBYRCF=REARTH/RngbrF
      RBYRCH=REARTH/RngbrH
      CSCRIF=COS(RngbrF/REARTH)
      CSCRIH=COS(RngbrH/REARTH)
      NLAP=NCRIT
      IF(IARGC().GT.1) THEN
        CALL GETARG(2,TITLEI)
        READ(TITLEI,*) NLAP
      END IF

C**** Open the input files
      do in=31,30+inm
        fn='work/fort.xx'
        write(fn(11:12),'(i2)') in
        open(in,file=fn,form='unformatted')
      end do

C**** Output text file
      open(78,file='work/fort.78',form='formatted') ! table of adj parameters

C**** Diagnostic output files in addition to log on standard output
      open(66,file='work/fort.66',form='formatted') ! combination info
      open(77,file='work/fort.77',form='formatted') ! station usage stats
      open(79,file='work/fort.79',form='formatted') ! isolated urban stations

C**** Read and interpret the header record
      READ (31) INFO,TITLEI
      KQ=INFO(2)
      IF(KQ.NE.1) THEN
         WRITE(6,*) ' PROGRAM NOT READY FOR QUANTITY ',KQ
         STOP 'INAPPROPRIATE QUANTITY'
      END IF
C     KM = the number of time frames per year
      KM=1
      IF(INFO(3).EQ.6) KM=12
      IF(INFO(3).EQ.7) KM=4
      IF(KM.NE.KM0) STOP 'ERROR: CHANGE KM0'
      ML=INFO(4)
      iyrm=INFO(4)
      NYRSIN=INFO(4)
      IF(IYRM0.LT.INFO(4)) THEN
         WRITE(6,'('' INCREASE IYRM0 TO AT LEAST'',I5)') INFO(4)
         STOP 'ERROR: IYRM0 NOT OK'
      END IF
      IYOFF=INFO(6)-1
      MBAD=INFO(7)
      LAST=INFO(7)
      XBAD=MBAD
C**** TRACE=INFO(8)  ! not needed here

C**** Collect the station data needed
      ISU=0
      ISR=0
      I1Snow=1
      DO 90 IN=1,INM
      REWIND 30+IN
      READ (30+IN) INFO
      MF=INFO(1)
      IF(ITRIM.GT.0) ML=INFO(8+ITRIM)
   50 IF(MF.GE.LAST) GO TO 90
      IF(I1Snow+ML-MF.GT.MSIZE) STOP 'ERROR: MSIZE TOO SMALL'
      CALL SREAD (30+IN,ML+1-MF,IDATA(I1Snow),ITRL,14+ITRIM)
      MFCUR=MF
      LENGTH=ML-MF+1
      MF=ITRL(14)
      IF(ITRIM.GT.0) ML=ITRL(14+ITRIM)
      XLAT=.1*ITRL(1)
      XLON=.1*ITRL(2)
      IF(NAME(31:32).EQ.' R'.or.NAME(31:31).EQ.'1') THEN
         ISR=ISR+1                          ! count rural stations
         IF(ISR.GT.NSTAM) STOP 'ERROR: NSTAM TOO SMALL'
         MFSR(ISR)=MFCUR
         I1SR(ISR)=I1Snow
         ILSR(ISR)=I1Snow+LENGTH-1
         CSLATR(ISR)=COS(XLAT*PI180)
         SNLATR(ISR)=SIN(XLAT*PI180)
         CSLONR(ISR)=COS(XLON*PI180)
         SNLONR(ISR)=SIN(XLON*PI180)
         IDR(ISR)=ITRL(3)
      ELSE
         ISU=ISU+1                          ! count non-rural stations
         IF(ISU.GT.NSTAM) STOP 'ERROR: NSTAM TOO SMALL'
         MFSU(ISU)=MFCUR
         I1SU(ISU)=I1Snow
         ILSU(ISU)=I1Snow+LENGTH-1
         CSLATU(ISU)=COS(XLAT*PI180)
         SNLATU(ISU)=SIN(XLAT*PI180)
         CSLONU(ISU)=COS(XLON*PI180)
         SNLONU(ISU)=SIN(XLON*PI180)
         IDU(ISU)=ITRL(3)
         CC(ISU)=NAME(34:36)
      END IF
      I1Snow=I1Snow+LENGTH
      GO TO 50
   90 CONTINUE
      NSTAu=ISU  ! total number of bright/urban or dim/sm.town stations
      NSTAr=ISR  ! total number of dark/rural stations
      LDTOT=I1Snow-1   !  total length of IDATA used
      write(*,*) 'number of rural/urban stations',NSTAr,NSTAu

C**** Convert data to real numbers (ann avgs were multiplied by 10)
      DO 115 N=1,LDTOT
      IF(IDATA(N).EQ.MBAD) THEN
         RDATA(N)=XBAD
      ELSE
         RDATA(N)=.1*IDATA(N) !  .01C to .1C
      END IF
  115 CONTINUE

C**** Order the NSTAr rural stations according to length of time record
      DO 121 IS=1,NSTAr
C 121 LEN(IS)=ILSR(IS)-I1SR(IS)+1 !  (counting gaps also)
      LEN(IS)=0
      nuseid(is)=0
      DO 120 M=I1SR(IS),ILSR(IS)
      IF(RDATA(M).EQ.XBAD) GO TO 120
      LEN(IS)=LEN(IS)+1
  120 CONTINUE
      write(*,*) 'rural station:',is,' id:',idr(is),' #ok',len(is)
  121 continue
      do m=1,nstar
        r(m)=len(m)+(m-1d0)/dfloat(nstar) ! make "lengths" unique
      end do
      CALL SORT (IORD,NSTAr,R)
      do m=1,nstar
        len(m)=r(m)
      end do
      do 122 IS=1,NSTAr
  122 write(*,*) 'rural station:',is,' id:',idr(iord(is)),' #ok',len(is)

C**** Combine time series for rural stations around each urban station
      DO 200 NURB=1,NSTAu
C**   Loop over all rural stations in memory; try 1/2 the radius first
      CSCRIT=CSCRIH
      RBYRC=RBYRCH
      Rngbr=RngbrH
  125 IS0=0            ! counter for rural neighbors
      DO 130 N=1,NSTAr
      ISR=IORD(N)
      IYU1=MFSU(NURB)+IYOFF-1   !  subtract 1 for a possible partial yr
      IYU2=IYU1+ILSU(NURB)-I1SU(NURB)+2  ! add 1 for partial year
C**   Find cosine of distance between center and station
      CSDBYR=SNLATR(ISR)*SNLATU(NURB)+CSLATR(ISR)*CSLATU(NURB)*
     *  (CSLONR(ISR)*CSLONU(NURB)+SNLONR(ISR)*SNLONU(NURB))
      X0=1950.
      IF(CSDBYR.LE.CSCRIT) GO TO 130
      DBYRC=0.
C**   The arc is replaced by the smaller chord (to avoid using ACOS)
      IF(CSDBYR.LT.1.) DBYRC=RBYRC*SQRT(2.*(1.-CSDBYR))
      IS0=IS0+1
      WTI(IS0)=1.-DBYRC
      ISOFI(IS0)=ISR
      lenis(IS0)=len(n)
  130 CONTINUE
C****
C**** Combine the station data
C****
      DO 150 M=1,IYRM
      WT(M)=0.
      IWT(M)=0
      URB(M)=XBAD
  150 AVG(M)=XBAD
      IF(IS0.EQ.0) THEN
         IF(CSCRIT.EQ.CSCRIH) THEN
            write(*,*) 'Trying full radius',RngbrF
            CSCRIT=CSCRIF
            RBYRC=RBYRCF
            Rngbr=RngbrF
            GO TO 125
         END IF
         WRITE(79,'('' NO RURAL NEIGHBORS FOR '',I9)') IDU(NURB)
         GO TO 200
      END IF
      ioff=MFSU(NURB)-I1SU(NURB)
      DO 151 M=I1SU(NURB),ILSU(NURB)
c     if(M+IOFF.gt.IYRM0) stop 231
  151 URB(M+IOFF)=RDATA(M)
      write(*,'(a,i9,a,i4,a,2i5,f9.0)') 'urb stnID:',idu(nurb),' # rur:'
     *  ,is0,' ranges:',MFSU(NURB)+IYOFF,ILSU(NURB)+ioff+IYOFF,Rngbr
C**** Start with the station with the longest time record
      IS=ISOFI(1)
      IOFF=MFSR(IS)-I1SR(IS)
      nuseid(is)=nuseid(is)+1
C**** First update of full data and weight arrays
      DO 160 M=I1SR(IS),ILSR(IS)
c     if(M+IOFF.gt.IYRM0) stop 244
      AVG(M+IOFF)=RDATA(M)
      IF(RDATA(M).LT.XBAD) WT(M+IOFF)=WTI(1)
      IF(RDATA(M).LT.XBAD) IWT(M+IOFF)=1
  160 CONTINUE
      write(*,'(a,i5,a,i4,i6,i10)') 'longest rur range:',
     *   MFSR(IS)+IYOFF,'-',ILSR(IS)+ioff+IYOFF,LENis(1),idr(is)
C**** Add in the remaining stations
      DO 190 I=2,IS0
      IS=ISOFI(I)
      IOFF=MFSR(IS)-I1SR(IS)
      write(*,'(a,i5,a,i5,a,i4,i6,i10)')'add stn',i,' range:',
     *   MFSR(IS)+IYOFF,'-',ILSR(IS)+ioff+IYOFF,LENis(i),idr(is)
C**** Extend the new data into a full series
      DO 170 M=1,IYRM
  170 DNEW(M)=XBAD
      DO 180 M=I1SR(IS),ILSR(IS)
c     if(M+IOFF.gt.IYRM0) stop 259
  180 DNEW(M+IOFF)=RDATA(M)
      NF1=MFSR(IS)
      NL1=ILSR(IS)+IOFF
C**** Shift new data, then combine them with current mean
      CALL CMBINE (AVG(1),WT,IWT, DNEW,NF1,NL1,WTI(I),
     *  IDR(IS),NSM,NCOM)
      write(*,*) 'data added: ',nsm,' overlap:',ncom,' years'
      IF(NSM.EQ.0) GO TO 190
      nuseid(is)=nuseid(is)+1
  190 CONTINUE
C**** Subtract urban station and call a curve fitting program
      IY1=1
  191 NXY=0
      N3=0
      N3F=0
      N3L=0
      NXX=0
      TMEAN=0
      if(IY1.eq.1)
     *   write(66,'(a,i9)')'year dTs-urban dTs-rural StnID=',idu(nurb)
      DO 195 IY=IY1,IYRM
      IF(AVG(IY).NE.XBAD.OR.URB(IY).NE.XBAD) NXX=NXX+1
      IF(AVG(IY).EQ.XBAD.OR.URB(IY).EQ.XBAD) GO TO 194
      if(IWT(IY).ge.NRURM) then
        N3L=IY
        N3=N3+1
        if(N3F.eq.0) N3F=IY
      end if
      if(N3.le.0) go to 195
      NXY=NXY+1
      TS(NXY)=AVG(IY)-URB(IY)
      F(NXY)=TS(NXY)
      TMEAN=TMEAN+F(NXY)
      YR(NXY)=IY+IYOFF
      X(NXY)=YR(NXY)-X0
      W(NXY)=1.
      if(IWT(IY).ge.NRURM) THEN
         NXY3=NXY
         TM3=TMEAN
      END IF
  194 if(nxx.gt.0.and.IY1.eq.1) write(66,'(i4,2f10.2)')
     *      IY+IYOFF,URB(IY),AVG(IY)
  195 CONTINUE
      IF(N3.LT.NCRIT) THEN
        IF(RBYRC.NE.RBYRCF) THEN
          write(*,*) 'trying full radius',RngbrF
          RBYRC=RBYRCF
          CSCRIT=CSCRIF
         Rngbr=RngbrF
          GO TO 125
        END IF
        WRITE(79,'(a3,i9.9,a13,i5,a15,i5,a50,a5)') CC(NURB),IDU(NURB),
     *   '  good years:',N3,'   total years:',N3L-N3F+1,
     *   ' too little rural-neighbors-overlap - drop station',' 9999'
        GO TO 200
      ELSE IF(FLOAT(N3).LT.XCRIT*(N3L-N3F+1.-.1)) THEN
        ! Not enough good years for the given range (<66%)
        ! The -0.1 is to prevent equality with potential uncertainty
        ! due to hardware rounding or compiler behaviour.
        ! Nick Barnes, Ravenbrook Limited, 2008-09-10
        ! Try to save cases in which the gaps are in the early part:
        IY1=N3L-(N3-1)/XCRIT
        if(IY1 < N3F+1) IY1=N3F+1                  ! avoid infinite loop
        WRITE(79,'(a3,i9.9,a17,i5,a1,i4)')
     *    CC(NURB),IDU(NURB),' drop early years',1+IYOFF,'-',IY1-1+IYOFF
        GO TO 191
      ELSE
        TMEAN=TM3/NXY3
        NXY=NXY3
        CALL GETFIT(15) ! 15=unit for debug output - currently not used
C**** Find extended range
        IYXTND=NINT(N3/XCRIT)-(N3L-N3F+1)
        write(*,*) 'possible range increase',IYXTND,N3,N3L-N3F+1
        n1x=N3F+IYOFF
        n2x=N3L+IYOFF
        IF(IYXTND.lt.0) stop 'impossible'
        if(IYXTND.gt.0) then
            LXEND=IYU2-(N3L+IYOFF)
            IF(IYXTND.le.LXEND) then
               n2x=n2x+LXEND
            else
               n1x=n1x-(IYXTND-LXEND)
               if(n1x.lt.IYU1) n1x=IYU1
               n2x=IYU2
            end if
        end if

C**** Write out a table entry for the table of adjustment parameters
        write(78,'(a3,i9.9,2f9.3,i5,5f9.3,I5,a1,I4,i5,a1,i4)') CC(nurb),
     *   idu(nurb),(fpar(i),i=1,2),nint(fpar(3)+X0),(fpar(i),i=4,6),
     *   (rmsp(i),i=1,2),N3F+IYOFF,'-',N3L+IYOFF,N1X,'-',N2X
      END IF
  200 CONTINUE
      nuse=0
      do n=1,NSTAr
      if(nuseid(n).gt.0) then
          write(77,*) 'used station ',idr(n),nuseid(n),' times'
          nuse=nuse+1
      end if
      end do
      write(77,*) nuse,' rural stations were used'
      STOP
      END

      SUBROUTINE CMBINE (AVG,WT,IWT, DNEW,NF1,NL1,WT1, ID,NSM,NCOM)
C****
C**** Bias of new data is removed by subtracting the difference
C**** over the common domain. Then the new data are averaged in.
C****
      COMMON/LIMIT/XBAD,NLAP
      DIMENSION AVG(*),DNEW(*),WT(*),IWT(*)
C**** Loop over years
      NSM=0
C**** Find means over common domain to compute bias
      SUMN=0
      SUM=0
      NCOM=0
      DO 10 N=NF1,NL1
      IF(AVG(N).GE.XBAD.OR.DNEW(N).GE.XBAD) GO TO 10
      NCOM=NCOM+1
      SUM=SUM+AVG(N)
      SUMN=SUMN+DNEW(N)
   10 CONTINUE
      IF(NCOM.LT.NLAP) RETURN
      BIAS=(SUM-SUMN)/FLOAT(NCOM)
C**** Update period of valid data, averages and weights
      DO 20 N=NF1,NL1
      IF(DNEW(N).GE.XBAD) GO TO 20
      WTNEW=WT(N)+WT1
      AVG(N)=(WT(N)*AVG(N)+WT1*(DNEW(N)+BIAS))/WTNEW
      WT(N)=WTNEW
      IWT(N)=IWT(N)+1
      NSM=NSM+1
   20 CONTINUE
   50 CONTINUE
      RETURN
      END

      SUBROUTINE SORT (INDX,NDIM,RLNGTH)
C**** Sorts INDEX and LNGTH such that rLNGTH becomes decreasing
      integer INDX(NDIM)
      real*8 RLNGTH(NDIM)
      DO 10 N=1,NDIM
   10 INDX(N)=N
      DO 30 N=1,NDIM-1
C**** Find maximum of LNGTH(N),...LNGTH(NDIM)
      NLMAX=N
      DO 20 NN=N+1,NDIM
      IF(RLNGTH(NN).GT.RLNGTH(NLMAX)) NLMAX=NN
   20 CONTINUE
C**** Switch positions N and NLMAX
      RMAX=RLNGTH(NLMAX)
      RLNGTH(NLMAX)=RLNGTH(N)
      RLNGTH(N)=RMAX
      IMAX=INDX(NLMAX)
      INDX(NLMAX)=INDX(N)
   30 INDX(N)=IMAX
      RETURN
      END

      SUBROUTINE SREAD (IN,LEN,IDATA,ITRL,NTRL)
C**** Speed read routine for input records
      DIMENSION IDATA(LEN),ITRL(NTRL)
      READ (IN) IDATA,ITRL
      RETURN
      END
