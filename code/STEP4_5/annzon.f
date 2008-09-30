C*********************************************************************
C *** program reads ZONAL monthly means and recomputes REGIONAL means
C *** as well as annual means.
C *** Input file:  10    zonal.means (ZON1977.T1200)
C ***
C *** Output files: 11    annual.zonal.means (ANNZON1977.T1200)
C ***               12    zonal.means (ZON1977.T1200) changed
C*********************************************************************
C****
C****
C**** Displays data and their annual means
C****
C**** Input file:    unit 10   (=output of job NCARSURF ZONAVG)
C****
C****  10: Record 1: INFO(1),...,INFO(8),TITLE       header record
C****      Record 2: DATA(1-->MONM),WT(1-->MONM),TITLE1     belt 1
C****      Record 3: DATA(1-->MONM),WT(1-->MONM),TITLE2     belt 2
C****            etc.
C****  DATA(1-->MONM) is a full time series, starting at January
C****  of year IYRBEG and ending at December of year IYREND.
C****  WT is proportional to the area containing valid data.
C****  TITLE1,... describe the latitude belts.
C****
C****  INFO(1),...,INFO(8)  are 4-byte integers,
C****  All TITLEs are 80-byte character strings,
C****  all other entries in records 2,3,... are 4-byte reals.
C****      INFO 1 and 5 are irrelevant
C****           2 = KQ (quantity flag, see below)
C****           3 = MAVG (time avg flag: 1 - 4 DJF - SON, 5 ANN,
C****                     6 MONTHLY, 7 SEAS, 8 - 19 JAN - DEC  )
C****           4 = MONM  (length of each time record)
C****           6 = IYRBEG (first year of each time record)
C****           7 = flag for missing data
C****           8 = flag for precipitation trace
C****  Flag for missing data:  XBAD = FLOAT( INFO(7) )
C****
C?*** Input parameters (time period - use LINFO to find it)
      PARAMETER (IY1TAB=1880)
C?*** Grid dimensions and limit parameters
C?*** Alternate global (IAVGGH=1) and hemispheric means (IAVGGH=2)
      PARAMETER (IAVGGH=2, JBM=8,NZS=3, JZM=JBM+NZS+3, JZP=14, MONMIN=6)
      CHARACTER*80 TITLE,TITL2,TITLEZ(JZM),BLANK/' '/,TITLEO/
     * 'ANNUALLY AVERAGED (7 OR MORE MONTHS) TEMPERATURE ANOMALIES (C)'/
      CHARACTER*10 TIT(3)/'    GLOBAL','N.HEMISPH.','S.HEMISPH.'/
C**** Input arrays
      DIMENSION INFO(8),INFOO(8)
      real*4, allocatable :: DATA(:,:,:),WT(:,:,:) ! (KM,IYRSX,JZM)
      INTEGER IORD(JZP)/14,12,13, 9,10,11, 1,2,3,4,5,6,7,8/,IOUT(18)
C**** Output-related arrays (dim of output array = MONM < MONM0)
      INTEGER JZG(4,2)/9,10,10,11, 9,4,5,11/
      REAL, dimension(4) :: WTSP=(/3.,2.,2.,3./)
      REAL, allocatable :: ANN(:,:),ANNW(:,:) !  (IYRSX,JZM)
C****
C**** Read and use the header record of an input file
C****
      open(10,file='work/ZON.Ts.ho2.GHCN.CL.PA.1200.step1',
     *     form='unformatted')
      open(96, file='result/ZonAnn.Ts.ho2.GHCN.CL.PA.txt',
     *     form='formatted')
      open(97, file='result/GLB.Ts.ho2.GHCN.CL.PA.txt',
     *     form='formatted')
      open(98, file='result/NH.Ts.ho2.GHCN.CL.PA.txt',
     *     form='formatted')
      open(99, file='result/SH.Ts.ho2.GHCN.CL.PA.txt',
     *     form='formatted')
      open(12, file='work/ZON.Ts.ho2.GHCN.CL.PA.1200',
     *     form='unformatted')
      open(11, file='work/ANNZON.Ts.ho2.GHCN.CL.PA.1200',
     *     form='unformatted')
      READ (10) INFO,TITLE,TITL2
      KQ=INFO(2)
C     KM = the number of time frames per year
      KM=1
      IF(INFO(3).EQ.6) KM=12
      IF(INFO(3).EQ.7) KM=4
      IYRBEG=INFO(6)
      MONM=INFO(4)
      IYRS=MONM/KM
      allocate ( data(km,iyrs,jzm),wt(km,iyrs,jzm) )
      allocate ( ann (   iyrs,jzm),annw( iyrs,jzm) )
      IYREND=INFO(4)/KM+IYRBEG-1
      MBAD=INFO(7)
      XBAD=MBAD
C**** Create and write out the header record of the output file
      WRITE(6,'(1X,A80)') TITLE
      WRITE(96,*) 'Annual Temperature Anomalies (.01 C) - ',TITLE(29:80)
      WRITE(97,'(A80)') TITLE
      WRITE(98,'(A80)') TITLE
      WRITE(99,'(A80)') TITLE
      WRITE(6,'(8I10)') INFO
      DO 10 I=1,8
   10 INFOO(I)=INFO(I)
      INFOO(3)=5
      INFOO(4)=INFO(4)/12
      WRITE(TITLEO(20:20),'(I1)') MONMIN
      WRITE(6,'(''1'',A80)') TITLEO
      WRITE(6,'(8I10)') INFOO
C****
C**** Collect the JZM zonal means (basic+special+hemi+global)
C****
      DO 30 JZ=1,JZM
      CALL SREAD(10,MONM,DATA(1,1,JZ),WT(1,1,JZ),TITLEZ(JZ))
   30 CONTINUE
C****
C**** Find the annual means
C****
      DO 100 JZ=1,JZM
      DO 100 IY=1,IYRS
      ANN(IY,JZ)=XBAD
      ANNW(IY,JZ)=0.
      ANNIY=0.
      ANNWIY=0.
      MON=0
      DO 50 M=1,KM
      IF(DATA(M,IY,JZ).EQ.XBAD) GO TO 50
      MON=MON+1
      ANNIY=ANNIY+DATA(M,IY,JZ)
      ANNWIY=ANNWIY+WT(M,IY,JZ)
   50 CONTINUE
      IF(MON.GE.MONMIN) ANN(IY,JZ)=ANNIY/MON
      ANNW(IY,JZ)=ANNWIY/12.
  100 CONTINUE
C****
C**** Alternate global mean (from North.lats,Equ.reg,South.lats)
C****
      IF(IAVGGH.EQ.0) GO TO 180
      IGLB=IAVGGH
      DO 120 IY=1,IYRS
      GLOB=0.
      ANN(IY,JZM)=XBAD
      DO 110 J=1,4
      IF(ANN(IY,JZG(J,IGLB)).EQ.XBAD) GO TO 120
      GLOB=GLOB+ANN(IY,JZG(J,IGLB))*WTSP(J)
  110 CONTINUE
      ANN(IY,JZM)=.1*GLOB
  120 CONTINUE
      DO 130 IY=1,IYRS
      DO 130 M=1,12
      DATA(M,IY,JZM)=XBAD
      GLOB=0.
      DO 125 J=1,4
      IF(DATA(M,IY,JZG(J,IGLB)).EQ.XBAD) GO TO 130
      GLOB=GLOB+DATA(M,IY,JZG(J,IGLB))*WTSP(J)
  125 CONTINUE
      DATA(M,IY,JZM)=.1*GLOB
  130 CONTINUE
C****
C**** Alternate hemispheric means (similar to alt. global mean)
C****
      IF(IAVGGH.EQ.1) GO TO 180
      DO 150 IHEM=1,2
      DO 140 IY=1,IYRS
      ANN(IY,IHEM+11)=XBAD
      IF(ANN(IY,IHEM+3).NE.XBAD.AND.ANN(IY,2*IHEM+7).NE.XBAD)
     *   ANN(IY,IHEM+11)=.4*ANN(IY,IHEM+3)+.6*ANN(IY,2*IHEM+7)
  140 CONTINUE
      DO 150 IY=1,IYRS
      DO 150 M=1,12
      DATA(M,IY,IHEM+11)=XBAD
      IF(DATA(M,IY,IHEM+3).NE.XBAD.AND.DATA(M,IY,2*IHEM+7).NE.XBAD)
     *  DATA(M,IY,IHEM+11)=.4*DATA(M,IY,IHEM+3)+.6*DATA(M,IY,2*IHEM+7)
  150 CONTINUE
C****
C**** Display the annual means
C****
  180 DO 190 JZ=1,JZP
  190 WRITE(6,'(I6,2X,A80)') JZ,TITLEZ(IORD(JZ))
      IYRSP=IYRS
      if(data(12,iyrs,14).gt.8000.) IYRSP=IYRS-1 ! skip incomplete year
      DO 200 IY=IY1TAB-IYRBEG+1,IYRSP
      IF( (IY+IYRBEG.gt.IY1TAB+5.and.MOD(IY+IYRBEG-2,20).EQ.0)
     * .or. IY.eq.IY1TAB-IYRBEG+1 ) THEN
        WRITE(96,*)
        WRITE(96,'(''                           24N   24S   90S   '',
     *        ''  64N   44N   24N   EQU   24S   44S   64S   90S'')')
        WRITE(96,'(''Year  Glob  NHem  SHem    -90N  -24N  -24S   '',
     *    '' -90N  -64N  -44N  -24N  -EQU  -24S  -44S  -64S Year'')')
      ENDIF
      IYR=IYRBEG+IY-1
      WRITE(96,'(I4,3(1X,I5),2X,3(1X,I5),2X,8(1X,I5),I5)') IYR,
     *  (NINT(100.*ANN(IY,IORD(JZ))),JZ=1,JZP),IYR
  200 CONTINUE
      WRITE(96,'(''Year  Glob  NHem  SHem     24N   24S   90S   '',
     *  ''  64N   44N   24N   EQU   24S   44S   64S   90S Year'')')
      WRITE(96,'(''                          -90N  -24N  -24S   '',
     *      '' -90N  -64N  -44N  -24N  -EQU  -24S  -44S  -64S'')')
      WRITE(96,*)
C****
C**** Display the monthly global and hemispheric means
C****
      DO 208 J=1,3
      WRITE(96+J,'(A10,'' Temperature Anomalies'',
     * '' in .01 C     base period: 1951-1980'')') TIT(J)
      DO 207 IY=IY1TAB-IYRBEG+1,IYRS
      IF( (IY+IYRBEG.gt.IY1TAB+5.and.MOD(IY+IYRBEG-2,20).EQ.0)
     *   .or. IY.eq.IY1TAB-IYRBEG+1  ) THEN
        WRITE(96+J,*)
        WRITE(96+J,'(''Year   Jan  Feb  Mar  Apr  May  Jun  Jul  Aug'',
     *''  Sep  Oct  Nov  Dec    J-D D-N    DJF  MAM  JJA  SON  Year'')')
      END IF
      awin=9999.
      if(iy.gt.1)
     *awin=data(12,IY-1,IORD(J))+data(1,IY,IORD(J))+data(2,IY,IORD(J))
      aspr=data(3,IY,IORD(J))+data(4,IY,IORD(J))+data(5,IY,IORD(J))
      asmr=data(6,IY,IORD(J))+data(7,IY,IORD(J))+data(8,IY,IORD(J))
      afl=data(9,IY,IORD(J))+data(10,IY,IORD(J))+data(11,IY,IORD(J))
      IOUT(15)=100*MBAD
      if(awin.lt.8000.) IOUT(15)=NINT(100.*awin/3.)
      IOUT(16)=100*MBAD
      if(aspr.lt.8000.) IOUT(16)=NINT(100.*aspr/3.)
      IOUT(17)=100*MBAD
      if(asmr.lt.8000.) IOUT(17)=NINT(100.*asmr/3.)
      IOUT(18)=100*MBAD
      if(afl.lt.8000.)  IOUT(18)=NINT(100.*afl/3.)
      IOUT(14)=100*MBAD
      ann2=awin+aspr+asmr+afl
      if(ann2.lt.8000.) IOUT(14)=NINT(100.*ann2/12.)
      IOUT(13)=100*MBAD
      ann1=ANN(IY,IORD(J))
      if(iy.eq.iyrs .and. data(12,IY,IORD(J)).gt.8000.) ann1=9999. 
      if(ann1.lt.8000.) IOUT(13)=NINT(100.*ANN(IY,IORD(J)))
      DO 202 M=1,12
  202 IOUT(M)=NINT(100.*DATA(M,IY,IORD(J)))
      IYR=IYRBEG-1+IY
  207 WRITE(96+J,'(I4,1X,12I5,2X,I5,I4,2X,4I5,I6)') IYR,IOUT,IYR
      WRITE(96+J,'(''Year   Jan  Feb  Mar  Apr  May  Jun  Jul  Aug'',
     *''  Sep  Oct  Nov  Dec    J-D D-N    DJF  MAM  JJA  SON  Year'')')
  208 CONTINUE
C****
C**** Save annual means on disk
C****
      WRITE(11) INFOO,TITLEO,TITLE(29:80),BLANK(1:28)
      DO 210 JZ=1,JZM
      CALL SWRITE(11,MONM/12,ANN(1,JZ),ANNW(1,JZ),TITLEZ(JZ))
  210 CONTINUE
C****
C**** Resave monthly means on disk
C****
      REWIND 10
      WRITE(12) INFO,TITLE,TITL2
      DO 310 JZ=1,JZM
      CALL SWRITE(12,MONM,DATA(1,1,JZ),WT(1,1,JZ),TITLEZ(JZ))
  310 CONTINUE
      STOP
      END

      SUBROUTINE SREAD (IN,NDIM,AR,WT,TITLE)
      CHARACTER*80 TITLE
      DIMENSION AR(NDIM),WT(NDIM)
      READ (IN) AR,WT,TITLE
C     WRITE(6,*) 'in ',TITLE,ndim,(AR(n),n=ndim-11,ndim)
      RETURN
      END

      SUBROUTINE SWRITE (IOUT,NDIM,AR,WT,TITLE)
      CHARACTER*80 TITLE
      DIMENSION AR(NDIM),WT(NDIM)
      WRITE (IOUT) AR,WT,TITLE
C     WRITE(6,*) 'out',TITLE,ndim,(AR(n),n=ndim-11,ndim)
      RETURN
      END
