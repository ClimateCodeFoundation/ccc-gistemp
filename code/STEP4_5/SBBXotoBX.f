C*********************************************************************
C *** PROGRAM READS NCAR FILES
C *** Input argument: Rland (0-1200km) radius of influence of station
C *** Input file:   10,11    subbox.data (land, ocean based)
C ***
C *** Output file:  12       box.data
C ***
C *** This program combines 8000 subbox to 80 box data
C*********************************************************************
C****
C**** This program interpolates the given station data or their
C**** ANOMALIES with respect to 1951-1980 to a prescribed grid.
C****
C**** Input file:    units 10 and 11
C****
C****  10: Record 1: M1,INFO(2),...,INFO(8),TITLE    header record
C****      Record 2: M2,LTS,LTN,LNW,LNE,n1,n2,d,AVG(1 --> M1)  box 1
C****      Record 3: M3,LTS,LTN,LNW,LNE,n1,n2,d,AVG(1 --> M2)  box 2
C****            etc.
C****
C****  11: Record 1: MO1,INFOO(2),...,INFOO(8),TITLEO   header record
C****      Record 2: MO2,LTS,LTN,LNW,LNE,0,n,0,AVGO(1 --> MO1)  box 1
C****      Record 3: MO3,LTS,LTN,LNW,LNE,0,n,0,AVGO(1 --> MO2)  box 2
C****            etc.
C****
C**** Output file :     unit  12 (regional means)
C****
C****  12: Record 1: INFO(1),...,INFO(8),TITLER     header record
C****      Record 2: AVGR(1-->M0),WTR(1-->M0),NG,LTS,LTN,LNW,LNE box 1
C****      Record 3: AVGR(1-->M0),WTR(1-->M0),NG,LTS,LTN,LNW,LNE box 2
C****            etc.
C****
C****  AVG(1-->M0) is a full time series, starting at January
C****  of year INFO(6), M0=INFO(4) is a multiple of 12.
C****  AVGO(1) corresponds to Jan of year INFOO(6), M0O=INFOO(4).
C****  WTR(M) is the area of the part of the region that contained
C****  valid data for month M (in square meters).
C****
C****  NG is the total number of non-missing data on that record.
C****  LTS,LTN is the latitude of the southern,northern edge of the
C****  (sub)box, LNW,LNE the longitude of the western,eastern edge
C****  in hundredths of degrees (all 4-byte integers).
C****
C****  INFO(1),...,INFO(8)  are 4-byte integers,
C****  TITLE,... are 80-byte character strings,
C****  AVG,AVGO,AVGR,WTR are 4-byte reals.
C****      INFO 1 = 1 (indicates that time series are not trimmed)
C****           2 = KQ (quantity flag, see below)
C****           3 = MAVG (time avg flag: 1 - 4 DJF - SON, 5 ANN,
C****                     6 MONTHLY, 7 SEAS, 8 - 19 JAN - DEC  )
C****           4 = MONM  (length of each time record)
C****           5 = MONM+4 (size of data record length)
C****           6 = IYRBEG (first year of each time record)
C****           7 = flag for missing data
C****           8 = flag for precipitation trace
C****  In the output file missing data are flagged by
C****  the real number  XBAD = FLOAT( INFO(7) )
C****
C**** The spatial averaging is done as follows:
C**** To remove the subbox bias, subbox data are shifted before
C**** combining them with the with the current mean (starting with
C**** subbox with the longest time series. The shift is such that
C**** the means over the time period they have in common remains
C**** unchanged (individually for each month).
C**** A final shift causes the 1961-1990 mean to become zero for
C**** each month. (If no values are available in 1961-1990,
C**** the average over the whole period is set to zero).
C****
C?*** Input parameters (crit.radius in km, limits for time period)
      PARAMETER (RCRIT=1200.,KQM=20)
C?*** Out put parameters (output time period, base period)
      PARAMETER (IYBASE=1961,LYBASE=1990)
C?*** Grid dimensions
      PARAMETER (NRM=80,NCM=100)
C?*** Limits
      PARAMETER (REARTH=6375.)
      CHARACTER*80 TITLE,TITLEO
C**** Input arrays
      DIMENSION INFO(8),INFOO(8)
C**** Output-related arrays (dim of output array = MONM < MONM0)
      DIMENSION WTM(12), XTRL(7),ITRL(7),ITRLO(7),
     *  IORDR(NCM+NCM),WGTC(NCM+NCM),DIST(NCM)
      real, allocatable :: avg(:,:),avgr(:),wtr(:)
C?*** Output grid dependent arrays (regions and centers)
      REAL*4 PI180,
     *  CSLONC(NCM,NRM),CSLATC(NCM,NRM),SNLONC(NCM,NRM),SNLATC(NCM,NRM)
      COMMON/GRIDC/PI180,XS(NRM),XN(NRM),XE(NRM),XW(NRM),
     *  CSLONC,CSLATC,SNLONC,SNLATC,AREA(NCM,NRM),LATLON(4,NCM+1,NRM)
      COMMON/LIMIT/XBAD,NOVRLP,BIAS(12)
      EQUIVALENCE (XTRL,ITRL)

      open(10,file='work/SBBX1880.Ts.GHCN.CL.PA.1200',
     *        form='unformatted')
      open(11,file='input/SBBX.HadR2', form='unformatted')
      open(12,file='result/BX.Ts.ho2.GHCN.CL.PA.1200',
     *        form='unformatted')

      PI180=DATAN(1.D0)/45.D0
      TOMETR=4.*180.*PI180 * REARTH**2/8000.
      NOVRLP=20
      if(iargc().ne.2) stop 'supply Rlnd(-1->1200),Rintrp(0->1)'
      call getarg(1,title)
      read(title,*) Rland
      if(Rland.ge.RCRIT) Rland=9999.
      if(Rland.lt.0.) Rland=-9999.
      call getarg(2,title)
      read(title,*) Rintrp
      if(Rintrp.lt.0.) Rintrp=0.
      if(Rintrp.gt.1.) Rintrp=1.
C****
C**** Read and use the header record of an input file
C****
      READ (10) INFO,TITLE
      READ (11) INFOO,TITLEO
      MNOW=INFO(1)
      MNOWO=INFOO(1)
      KQ=INFO(2)
C     KM = the number of time frames per year
      KM=1
      IF(INFO(3).EQ.6) KM=12
      IF(INFO(3).EQ.7) KM=4
      NYRSIN=INFO(4)/KM
      IYRBEG=INFO(6)
      IYRBGO=INFOO(6)
      IYRBGC=MIN(IYRBEG,IYRBGO)
      I1TIN=1+12*(IYRBEG-IYRBGC)
      I1TINO=1+12*(IYRBGO-IYRBGC)
      MONM=INFO(4)
      MONMO=INFOO(4)
      MONMC=MAX(MONM+I1TIN-1,MONMO+I1TINO-1)
      allocate ( AVG(monmc,NCM+NCM),AVGR(monmc),WTR(monmc) )

      NFB=1+IYBASE-IYRBGC
      NLB=1+LYBASE-IYRBGC
      MBAD=INFO(7)
      LAST=INFO(7)
      XBAD=MBAD
      TRACE=INFO(8)
      WRITE(6,'('' INFO SUBBOX'',8I10)') INFO
      INFO(4)=MONM         ! better MONMC
      INFO(5)=2*MONM+5     ! better 2*MONMC+5
      INFO(6)=IYRBGC
      WRITE(6,'('' INFO BOX   '',8I10)') INFO
      WRITE (12) INFO,TITLE
C**** Find grid-dependent arrays (GRIDC-COMMON)
      CALL GRIDEA (NRM,NCM)
C**** Convert areas to units of square-meters
      DO 20 NR=1,NRM
      DO 20 NC=1,NCM
   20 AREA(NC,NR)=AREA(NC,NR)*TOMETR
C****
C**** Loop over NRM large regions
C****
      DO 300 NR=1,NRM
C**** Collect the subbox data needed for region NR
      DO 200 NC=1,NCM
      DO 100 M=1,MONMC
      AVG(M,NC)=XBAD
  100 AVG(M,NC+NCM)=XBAD
      CALL SREAD(10,MNOW,AVG(I1TIN,NC),XTRL,MNEW)
      MNOW=MNEW
      DIST(NC)=XTRL(7)
      if(DIST(NC).LE.Rland) READ(11) MNEWO
      if(DIST(NC).GT.Rland)
     *   CALL SREAD(11,MNOWO,AVG(I1TINO,NC+NCM),ITRLO,MNEWO)
      MNOWO=MNEWO
      WGTC(NC)=0
      WGTC(NC+NCM)=0
      DO 150 M=1,MONMC
      IF(AVG(M,NC).NE.XBAD) WGTC(NC)=WGTC(NC)+1
      IF(AVG(M,NC+NCM).NE.XBAD) WGTC(NC+NCM)=WGTC(NC+NCM)+1
  150 CONTINUE
      wocn=max(0.,(DIST(NC)-Rland)/(RCRIT-Rland))
      if(Rland.eq.Xbad) wocn=0.
      if(WGTC(nc+ncm).lt.12*NOVRLP) wocn=0.
      if(wocn.gt.Rintrp) wocn=1.
      WGTC(NC)=0.
      WGTC(NC+NCM)=0.
      DO 160 M=1,MONMC
      IF(AVG(M,NC).NE.XBAD) WGTC(NC)=WGTC(NC)+(1-wocn)
      IF(AVG(M,NC+NCM).NE.XBAD) WGTC(NC+NCM)=WGTC(NC+NCM)+wocn
  160 CONTINUE
  200 CONTINUE
C****
C**** Find the regionally averaged time series
C****
      CALL SORT (IORDR,NCM+NCM,WGTC)
      NC=IORDR(1)
      ncr=nc
      if(nc.gt.ncm) ncr=nc-ncm
      wocn=max(0.,(DIST(NCr)-Rland)/(RCRIT-Rland))
      if(Rland.eq.Xbad) wocn=0.
      if(nc.eq.ncr.and.WGTC(nc+ncm).lt.12*NOVRLP) wocn=0.
      if(wocn.gt.Rintrp) wocn=1.
      wnc=wocn
      if(nc.le.ncm) wnc=1.-wocn
      DO 205 K=1,KM
      WTM(K)=AREA(NC,NR)*wnc
  205 BIAS(K)=0.
      DO 210 M=1,MONMC
      WTR(M)=0.
      IF(AVG(M,NC).LT.XBAD) WTR(M)=AREA(NCr,NR)*wnc
  210 AVGR(M)=AVG(M,NC)
C**** Add in the series of the remaining centers in region NR
      DO 220 N=2,NCM+NCM
      NC=IORDR(N)
      IF(WGTC(N).LT.12*NOVRLP) GO TO 230
      ncr=nc
      if(nc.gt.ncm) ncr=nc-ncm
      wocn=max(0.,(DIST(NCr)-Rland)/(RCRIT-Rland))
      if(Rland.eq.Xbad) wocn=0.
      if(nc.eq.ncr.and.WGTC(nc+ncm).lt.12*NOVRLP) wocn=0.
      if(wocn.gt.Rintrp) wocn=1.
      wnc=wocn
      if(nc.le.ncm) wnc=1.-wocn
      WT1=AREA(NCr,NR)*wnc
      CALL CMBINE (AVGR,WTR,AVG(1,NC),1,MONMC/KM,WT1,WTM,KM,
     *  NC,NSM)
  220 CONTINUE
  230 CONTINUE
C**** Set BIAS=time average over the base period if IYBASE > 0
      IF(NFB.GT.0) CALL TAVG (AVGR,KM,NYRSIN, NFB,NLB, NR,0, 0.)
      NGOOD=0
      M=0
      DO 240 IY=1,MONMC/KM
      DO 240 K=1,KM
      M=M+1
      IF(AVGR(M).EQ.XBAD) GO TO 240
      AVGR(M)=AVGR(M)-BIAS(K)
      NGOOD=NGOOD+1
  240 CONTINUE
      WRITE(6,'(1X,12I4,5X,12I4)')(NINT(10.*AVGR(M)),M=1,MONMC)
      WRITE(6,'(1X,12I4,5X,12I4)')(NINT(WTR(M)/TOMETR),M=1,MONMC)
  290 CALL SWRITE(12,MONM,AVGR,WTR,NGOOD,LATLON(1,NCM+1,NR)) ! MONMC !!!
  300 WRITE(6,'('' REGION,Months w/data'',2I9)') NR,NGOOD
      STOP
      END

      SUBROUTINE GRIDEA (NRM1,NCM1)
C****
C**** This output grid dependent routine sets up the latitude and
C**** longitude limits for the 1200 km hull of the large regions
C**** and computes the relevant grid quantities.
C****
C**** Current order of boxes:   1-4   north , west->east ...
C****                            .      to
C****                          76-80  south , west->east
C**** Order of subboxes:        1-10  south , west->east ...
C****                            .      to
C****                         91-100  north , west->east
C****
      PARAMETER (NRM=80,ICM=10,JCM=10,  NCM=ICM*JCM)
      REAL*4 PI180,SNLATJ(JCM),CSLATJ(JCM),LT100S(JCM),LT100N(JCM),
     *  CSLONC(NCM,NRM),CSLATC(NCM,NRM),SNLONC(NCM,NRM),SNLATC(NCM,NRM)
      COMMON/GRIDC/PI180,XS(NRM),XN(NRM),XE(NRM),XW(NRM),
     *  CSLONC,CSLATC,SNLONC,SNLATC,AREA(NCM,NRM),LATLON(4,NCM+1,NRM)
C****
C**** Sergej's equal area grid
C****
C**** Grid constants for latitude zones
CCCCC REAL BANDS(9)/90.,64.2,44.4,23.6,0.,-23.6,-44.4,-64.2,-90./
C**** We need only the SINEs of the (northern) band edges
      REAL SNNEDG(9)/1.,.9,.7,.4,0.,-.4,-.7,-.9,-1./
C     Input data sets - 1: 90N-60N  2: 60N-30N ... 6: 60S-90S 7:short
      INTEGER INFRST(8)/1,1,2,2,3,4,5,5/,NUMJ(8)/4,8,12,16,16,12,8,4/
      INTEGER INLAST(8)/2,2,3,4,5,5,6,6/
C**** Check grid dimensions
      IF(NRM1.NE.NRM.OR.ICM*JCM.NE.NCM1) THEN
         WRITE(6,'(6I9)') NRM1,NRM, ICM*JCM,NCM1
         STOP 'ERROR: GRID SIZES INCONSISTENT'
      END IF
C**** Loop over all BOXES (large regions)
      NR=0
      DO 50 J=1,8
C**** Find the sin of the latitudes of the centers in band J
      SNN=SNNEDG(J)
      SNS=SNNEDG(J+1)
      DSLATJ=(SNN-SNS)/JCM
      DO 20 JC=1,JCM
      LT100S(JC)=NINT(100./PI180*ASIN(SNS+(JC-1)*DSLATJ))
      LT100N(JC)=NINT(100./PI180*ASIN(SNS+JC*DSLATJ))
      SNLATJ(JC)=SNS+(JC-.5)*DSLATJ
   20 CSLATJ(JC)=SQRT(1.-SNLATJ(JC)**2)
      DO 50 I=1,NUMJ(J)
      NR=NR+1
      DLON=360./NUMJ(J)
      XEAST=-180.+I*DLON
      XWEST=XEAST-DLON
C**** Find sine and cosine of the center latitudes and longitudes
      DLONC=DLON/ICM
      DO 40 IC=1,ICM
      LN100W=NINT(100.*(XWEST+(IC-1)*DLONC))
      LN100E=NINT(100.*(XWEST+IC*DLONC))
      RLONI=PI180*(XWEST+(IC-.5)*DLONC)
      SNLONI=SIN(RLONI)
      CSLONI=COS(RLONI)
      DO 40 JC=1,JCM
CSRG  ICJC=JC+JCM*(IC-1)
      ICJC=IC+ICM*(JC-1)
      LATLON(1,ICJC,NR)=LT100S(JC)
      LATLON(2,ICJC,NR)=LT100N(JC)
      LATLON(3,ICJC,NR)=LN100W
      LATLON(4,ICJC,NR)=LN100E
      AREA(ICJC,NR)=1.
      SNLONC(ICJC,NR)=SNLONI
      CSLONC(ICJC,NR)=CSLONI
      SNLATC(ICJC,NR)=SNLATJ(JC)
   40 CSLATC(ICJC,NR)=CSLATJ(JC)
      LATLON(1,NCM+1,NR)=LT100S(1)
      LATLON(2,NCM+1,NR)=LT100N(JCM)
      LATLON(3,NCM+1,NR)=NINT(100.*XWEST)
      LATLON(4,NCM+1,NR)=NINT(100.*(XWEST+DLON))
   50 CONTINUE
      RETURN
      END

      SUBROUTINE CMBINE (AVG,WT, DNEW,NF1,NL1,WT1,WTM, KM, ID,NSM)
C****
C**** Bias of new data is removed by subtracting the difference
C**** over the common domain. Then the new data are averaged in.
C****
      COMMON/LIMIT/XBAD,NOVRLP,BIAS(12)
      DIMENSION AVG(KM,*),DNEW(KM,*),WT(KM,*),WTM(*),MISSNG(12)
C**** Loop over months or seasons if appropriate
      NSM=0
      MISSED=KM
      DO 50 K=1,KM
      MISSNG(K)=1
C**** Find means over common domain to compute bias
      SUMN=0
      SUM=0
      NCOM=0
      DO 10 N=NF1,NL1
      IF(AVG(K,N).GE.XBAD.OR.DNEW(K,N).GE.XBAD) GO TO 10
      NCOM=NCOM+1
      SUM=SUM+AVG(K,N)
      SUMN=SUMN+DNEW(K,N)
   10 CONTINUE
      IF(NCOM.LT.NOVRLP) GO TO 50
      BIASK=(SUM-SUMN)/FLOAT(NCOM)
C**** Find mean bias
      WTMNEW=WTM(K)+WT1
      BIAS(K)=(WTM(K)*BIAS(K)+WT1*BIASK)/WTMNEW
      WTM(K)=WTMNEW
C**** Update period of valid data, averages and weights
      DO 20 N=NF1,NL1
      IF(DNEW(K,N).GE.XBAD) GO TO 20
      WTNEW=WT(K,N)+WT1
      AVG(K,N)=(WT(K,N)*AVG(K,N)+WT1*(DNEW(K,N)+BIASK))/WTNEW
      WT(K,N)=WTNEW
      NSM=NSM+1
   20 CONTINUE
      MISSED=MISSED-1
      MISSNG(K)=0
   50 CONTINUE
      IF(MISSED.GT.0) WRITE(6,90) ID,WT1,MISSNG
   90 FORMAT(' UNUSED DATA - ID/SUBBOX,WT',I8,F5.2,12I2)
      RETURN
      END

      SUBROUTINE TAVG (DATA,KM,NYRS, NFB,NLB, NR,NC, DEFLT)
C****
C**** TAVG computes the time averages (separately for each calendar
C**** month if KM=12) over the base period (year NFB to NLB) and
C**** saves them in BIAS. In case of no data, the average is set to
C**** DEFLT if NR=0 or computed over the whole period if NR>0.
C****
      COMMON/LIMIT/XBAD,NOVRLP,BIAS(12)
      DIMENSION DATA(KM,*),LEN(12)
      MISSED=KM
      DO 50 K=1,KM
      BIAS(K)=DEFLT
      SUM=0.
      M=0
      DO 10 N=NFB,NLB
      IF(DATA(K,N).GE.XBAD) GO TO 10
      M=M+1
      SUM=SUM+DATA(K,N)
   10 CONTINUE
      LEN(K)=M
      IF(M.EQ.0) GO TO 50
      BIAS(K)=SUM/FLOAT(M)
      MISSED=MISSED-1
   50 CONTINUE
      IF(NR*MISSED.EQ.0) RETURN
C**** If base period is data free, use bias with respect to whole series
      DO 100 K=1,KM
      IF(LEN(K).GT.0) GO TO 100
      WRITE(6,'(''0NO DATA IN BASE PERIOD - MONTH,NR,NC'',3I9)') K,NR,NC
      SUM=0.
      M=0
      DO 60 N=1,NYRS
      IF(DATA(K,N).GE.XBAD) GO TO 60
      M=M+1
      SUM=SUM+DATA(K,N)
   60 CONTINUE
      IF(M.EQ.0) GO TO 100
      BIAS(K)=SUM/FLOAT(M)
  100 CONTINUE
      RETURN
      END

      SUBROUTINE SORT (INDEX,NDIM,WEIGHT)
C**** Sorts INDEX and WEIGHT such that WEIGHT becomes decreasing
      DIMENSION INDEX(NDIM),WEIGHT(NDIM)
      DO 10 N=1,NDIM
   10 INDEX(N)=N
      DO 30 N=1,NDIM-1
C**** Find maximum of WEIGHT(N),...WEIGHT(NDIM)
      NWMAX=N
      DO 20 NN=N+1,NDIM
      IF(WEIGHT(NN).GT.WEIGHT(NWMAX)) NWMAX=NN
   20 CONTINUE
C**** Switch positions N and NWMAX
      WMAX=WEIGHT(NWMAX)
      WEIGHT(NWMAX)=WEIGHT(N)
      WEIGHT(N)=WMAX
      IMAX=INDEX(NWMAX)
      INDEX(NWMAX)=INDEX(N)
   30 INDEX(N)=IMAX
      RETURN
      END

      SUBROUTINE SREAD (IN,LEN,RDATA,ITRL,MNEXT)
C**** Speed read routine for input records
      DIMENSION RDATA(LEN),ITRL(7)
      READ (IN) MNEXT,ITRL,RDATA
      RETURN
      END

      SUBROUTINE SWRITE (IOUT,NDIM,ARRAY1,ARRAY2,NG,LATLON)
      DIMENSION ARRAY1(NDIM),ARRAY2(NDIM),LATLON(4)
      WRITE(IOUT) ARRAY1,ARRAY2,NG,LATLON
      RETURN
      END

