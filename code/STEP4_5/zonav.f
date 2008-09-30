C*********************************************************************
C *** PROGRAM READS BXdata
C *** Input files:  11    box.data (BX1977.T1200)
C ***
C *** Output files: 10    zonal.means (ZON1977.T1200)
C ***
C*********************************************************************
C****
C**** This program combines the given gridded data (anomalies)
C**** to produce AVERAGES over various LATITUDE BELTS.
C****
C**** Input file:    unit 11   (=output of job NCARSURF GRIDDING)
C****
C****  11: Record 1: INFOI(1),...,INFOI(8),TITLEI   header record
C****      Record 2: AR(1-->MONM0),WTR(1-->MONM0),NG  grid point 1
C****      Record 3: AR(1-->MONM0),WTR(1-->MONM0),NG  grid point 2
C****            etc.
C****
C**** Output files:     unit 10
C****
C****  10: Record 1: INFO(1),...,INFO(8),TITLEO,TXT    header record
C****      Record 2: DATA(1-->MONM),WT(1-->MONM),TITLE1     belt 1
C****      Record 3: DATA(1-->MONM),WT(1-->MONM),TITLE2     belt 2
C****            etc.
C****  DATA(1-->MONM) is a full time series, starting at January
C****  of year IYRBEG and ending at December of year IYREND.
C****  WT is proportional to the area containing valid data.
C****  AR(1) refers to Jan of year IYRBG0 which may
C****  be less than IYRBEG, MONM0 is the length of an input time
C****  series, and WTR(M) is the area of the part of the region
C****  that contained valid data for month M.
C****  NG is the total number of non-missing data on that record.
C****  TITLE1,... describe the latitude belts.
C****
C****  INFO(1),...,INFO(8)  are 4-byte integers,
C****  All TITLEs and TXT are 80-byte character strings,
C****  all other entries in records 2,3,... are 4-byte reals.
C****      INFO 1 and 5 are irrelevant
C****           2 = KQ (quantity flag, see below)
C****           3 = MAVG (time avg flag: 1 - 4 DJF - SON, 5 ANN,
C****                     6 MONTHLY, 7 SEAS, 8 - 19 JAN - DEC  )
C****           4 = MONM  (length of each time record)
C****           6 = IYRBEG (first year of each time record)
C****           7 = flag for missing data
C****           8 = flag for precipitation trace
C****  INFO(I)=INFOI(I) except perhaps for I=4 and I=6.
C****  In the output file missing data are flagged by
C****  the real number  XBAD = FLOAT( INFO(7) )
C****
C**** JBM zonal means are computed first, combining successively
C**** the appropriate regional data (AR with weight WTR). To remove
C**** the regional bias, the data of a new region are shifted
C**** so that the mean over the common period remains unchanged
C**** after its addition. If that common period is less than
C**** 20(NCRIT) years, the region is disregarded. To avoid that
C**** case as much as possible, regions are combined in order of
C**** the length of their time record. A final shift causes the
C**** 1951-1980 mean to become zero (for each month).
C****
C**** All other means (incl. hemispheric and global means) are
C**** computed from these zonal means using the same technique.
C**** NOTE: the weight of a zone may be smaller than its area
C**** since data-less parts are disregarded; this also causes the
C**** global mean to be different from the mean of the hemispheric
C**** means.
C****
C?*** Input parameters (# of input files, time period)
C?*** Out put parameters (output time period, base period)
      PARAMETER (IYBASE=1951,LYBASE=1980)
C?*** Grid dimensions and limit parameters
      PARAMETER (NRM=80,JBM=8,IBMM=16,NZS=3 ,NCRIT=20)
      CHARACTER*80 TITLEI,TITLEO,TITLEZ(JBM+NZS+3)/
     * JBM*'see subr GRIDx',NZS*'see subr GRIDx',
     * 'NORTHERN HEMISPHERE','SOUTHERN HEMISPHERE','GLOBAL MEANS'/
C?*** Special title for Jim's plots
      CHARACTER*40 TXT(2)/
     *  ' zones:  90->64.2->44.4->23.6->0->-23.6-',
     *  '>-44.4->-64.2->-90                      '/
C**** Input arrays
      DIMENSION INFOI(8),LENR(IBMM)
      real*4, allocatable :: ar(:,:),wtr(:,:)
C**** Output-related arrays (dim of output array = MONM < MONM0)
      DIMENSION INFO(8),IORD(JBM+IBMM),LENZ(JBM)
      real*4, allocatable :: AVG(:,:),WT(:,:),AVGG(:),WTG(:)
      COMMON/GRIDC/IBM(JBM),KZONE(JBM,NZS+3)
      COMMON/LIMIT/XBAD,NOVRLP,BIAS(12)
      NOVRLP=NCRIT
      open(10,file='work/ZON.Ts.ho2.GHCN.CL.PA.1200.step1',
     *        form='unformatted')

C****
C**** Read and use the header record of an input file
C****
      open(11,file='result/BX.Ts.ho2.GHCN.CL.PA.1200',
     *        form='unformatted')

      READ (11) INFOI,TITLEI
      KQ=INFOI(2)
C     KM = the number of time frames per year
      KM=1
      IF(INFOI(3).EQ.6) KM=12
      IF(INFOI(3).EQ.7) KM=4
      ML=INFOI(4)
      NYRSIN=INFOI(4)/KM
      IYRBEG=INFOI(6)
      MONM=ML ; IYREND = MONM/KM + IYRBEG - 1
      allocate ( ar(ml,IBMM),wtr(ml,IBMM), AVG(ml,JBM),WT(ml,JBM) )
      allocate ( AVGG(ml),WTG(ml) )
      MFOUT=1+(IYRBEG-INFOI(6))*KM
      NFB=1+IYBASE-INFOI(6)
      NLB=1+LYBASE-INFOI(6)
      MBAD=INFOI(7)
      LAST=INFOI(7)
      XBAD=MBAD
      TRACE=INFOI(8)
C**** Create and write out the header record of the output file
      DO 10 I=1,8
   10 INFO(I)=INFOI(I)
      INFO(4)=MONM
      INFO(6)=IYRBEG
      TITLEO=TITLEI
      WRITE(6,'(1X,A80)') TITLEI,TITLEO
      WRITE(6,'(8I10)') INFOI
      WRITE(6,'(8I10)') INFO
      WRITE (10) INFO,TITLEO,TXT
C**** Find grid-dependent arrays (GRIDC-COMMON)
      CALL GRIDEA (NRM,TITLEZ,IBMM,NZS)
C****
C**** Find the means over the JBM basic latitude belts
C****
      DO 300 JB=1,JBM
C**** Collect the data for the IBM(JB) regions in the belt JB
      DO 30 IB=1,IBM(JB)
      CALL SREAD(11,ML   ,AR(1,IB),WTR(1,IB),LENR(IB))
   30 CONTINUE
C**** Check whether there are any data
      LNTOT=0
      DO 40 N=1,IBM(JB)
      LNTOT=LNTOT+LENR(N)
   40 CONTINUE
      IF(LNTOT.EQ.0) THEN
         WRITE(6,*) ' **** NO DATA FOR ZONE ',JB,TITLEZ(JB)
         GO TO 300
      END IF
C**** Order the IBM regions according to length of time record
      CALL SORT (IORD,IBM(JB),LENR)
      NR=IORD(1)
C**** Start with longest time series
      DO 210 M=1,MONM
      WT(M,JB)=WTR(M,NR)
  210 AVG(M,JB)=AR(M,NR)
C**** Add in the series of the remaining regions in belt JB
      DO 220 N=2,IBM(JB)
      NR=IORD(N)
      IF(LENR(N).EQ.0) GO TO 230
      CALL CMBINE (AVG(1,JB),WT(1,JB), AR(1,NR), 1,NYRSIN,
     *             WTR(1,NR),WTM, KM, JB)
  220 CONTINUE
  230 CONTINUE
C**** Set BIAS=time average over the base period if IYBASE > 0
      CALL TAVG(AVG(1,JB),KM,NYRSIN, NFB,NLB, JB, 0.)
      LENZ(JB)=0
      M=0
      DO 240 IY=1,NYRSIN
      DO 240 K=1,KM
      M=M+1
      IF(AVG(M,JB).EQ.XBAD) GO TO 240
      AVG(M,JB)=AVG(M,JB)-BIAS(K)
      LENZ(JB)=LENZ(JB)+1
  240 CONTINUE
      write(*,*) 'zonal mean',JB,TITLEZ(JB)
      WRITE(6,'(1X,12I4,5X,12I4)')(NINT(10.*AVG(M,JB)),M=1,MONM)
      write(*,*) 'weights'
      WRITE(6,'(1X,12I4,5X,12I4)')(NINT(WT(M,JB)/60000.+.499),M=1,MONM)
      CALL SWRITE(10,MONM,AVG(MFOUT,JB),WT(MFOUT,JB),TITLEZ(JB))
      WRITE(6,*) AVG(MFOUT,JB),WT(MFOUT,JB),TITLEZ(JB),' SAVED ON DISK'
  300 CONTINUE
C****
C**** Find the means for the remaining NZS+3 bands
C****
      CALL SORT (IORD,JBM,LENZ)
      DO 400 JZ=1,NZS+3
C**** Start with the longest time series
      IF(LENZ(1).EQ.0) THEN
         WRITE(6,*) ' **** NO DATA FOR ZONE ',JBM+JZ,TITLEZ(JBM+JZ)
         GO TO 400
      END IF
      DO 310 J1=1,JBM
      IF(KZONE(IORD(J1),JZ).GT.0) GO TO 320
  310 CONTINUE
  320 JB=IORD(J1)
      DO 330 M=1,MONM
      WTG(M)=WT(M,JB)
  330 AVGG(M)=AVG(M,JB)
C**** Add in the remaining latitude belts JB with KZONE(JB,JZ)=1
      DO 340 J=J1+1,JBM
      IF(KZONE(IORD(J),JZ).EQ.0) GO TO 340
      JB=IORD(J)
      IF(LENZ(J).EQ.0) GO TO 360
      CALL CMBINE (AVGG,WTG,AVG(1,JB), 1,NYRSIN, WT(1,JB),WTM,KM,JBM+JZ)
  340 CONTINUE
C**** Set BIAS=time average over the base period if IYBASE > 0
      IF(NFB.GT.0) CALL TAVG(AVGG,KM,NYRSIN, NFB,NLB, JBM+JZ, 0.)
      M=0
      DO 350 IY=1,NYRSIN
      DO 350 K=1,KM
      M=M+1
      IF(AVGG(M).NE.XBAD) AVGG(M)=AVGG(M)-BIAS(K)
  350 CONTINUE
  360 CONTINUE
      write(*,*) 'zone',JBM+JZ
      WRITE(6,'(1X,12I4,5X,12I4)')(NINT(10.*AVGG(M)),M=1,MONM)
      write(*,*) 'weights',JBM+JZ
      WRITE(6,'(1X,12I4,5X,12I4)')(NINT(WTG(M)/60000.+.49999),M=1,MONM)
      CALL SWRITE(10,MONM,AVGG(MFOUT),WTG(MFOUT),TITLEZ(JBM+JZ))
      WRITE(6,*) AVGG(MFOUT),WTG(MFOUT),TITLEZ(JBM+JZ),' SAVED ON DISK'
  400 CONTINUE
      STOP
      END

      SUBROUTINE GRIDEA (NRM0,TITLEZ, IBMM0,NZS0)
C****
C**** This output grid dependent routine sets the parameters
C**** IBM(J) (= number of regions in latitude belt J)  and
C**** KZONE(J,N) (=1 if belt J belongs to domain N) and the TITLES
C****
C**** Current order of boxes:   1-4   north , west->east ...
C****                            .      to
C****                          77-80  south , west->east
C**** Order of subboxes:        1-10  south , west->east ...
C****                            .      to
C****                         91-100  north , west->east
C****   Sergei's ordering       1-10  west , south->north
C****                            .      to
C****                         91-100  east , south->north
C****
      PARAMETER (NRM=80,ICM=10,JCM=10, NZS=3, NCM=ICM*JCM)
      CHARACTER*80 TITLEZ(*),TITLES(1+NZS)/
     *    '  LATITUDE BELT FROM  XX.X ? TO  XX.X ?',
     *    '  NORTHERN LATITUDES: 23.6 N TO  90.0 N',
     *    '       LOW LATITUDES: 23.6 S TO  23.6 N',
     *    '  SOUTHERN LATITUDES: 90.0 S TO  23.6 S'/
      COMMON/GRIDC/IBM(8),KZONE(8,NZS+3)
C****
C**** Sergej's equal area grid
C****
C**** Grid constants for latitude zones
      REAL BANDS(9)/90.,64.2,44.4,23.6,0.,-23.6,-44.4,-64.2,-90./
      INTEGER NUMJ(8)/4,8,12,16,16,12,8,4/
C**** Check grid dimensions
      IF(NRM0.NE.NRM.OR.NZS.NE.NZS0) THEN
         WRITE(6,'(4I9)') NRM0,NRM, NZS0,NZS
         STOP 'ERROR: INCONSISTENCY FOR NRM NCM NZS'
      END IF
C**** Loop over all ZONES
      DO 50 J=1,8
      IBM(J)=NUMJ(J)
      IF(IBM(J).LE.IBMM0) GO TO 20
      WRITE(6,*) 'IBMM SHOULD BE AT LEAST EQUAL TO',IBM(J),' J=',J
      STOP 'IBMM TOO SMALL'
C**** Define the large bands as combination of the basic zones
   20 KZONE(J,1)=0
      IF(J.LE.3) KZONE(J,1)=1
      KZONE(J,2)=0
      IF(J.EQ.4.OR.J.EQ.5) KZONE(J,2)=1
      KZONE(J,3)=0
      IF(J.GE.6) KZONE(J,3)=1
      KZONE(J,NZS+1)=0
      IF(J.LE.4) KZONE(J,NZS+1)=1
      KZONE(J,NZS+2)=1-KZONE(J,NZS+1)
      KZONE(J,NZS+3)=1
C**** Set the titles TITLEZ for the basic latitude bands
      TITLEZ(J)=TITLES(1)
      WRITE(TITLEZ(J)(23:26),'(F4.1)') ABS(BANDS(J+1))
      WRITE(TITLEZ(J)(28:28),'(A1)') 'N'
      IF(BANDS(J+1).LT.0.) WRITE(TITLEZ(J)(28:28),'(A1)') 'S'
      IF(BANDS(J+1).EQ.0.) WRITE(TITLEZ(J)(22:28),'(A7)') 'EQUATOR'
      WRITE(TITLEZ(J)(34:37),'(F4.1)') ABS(BANDS(J))
      WRITE(TITLEZ(J)(39:39),'(A1)') 'N'
      IF(BANDS(J).LT.0.) WRITE(TITLEZ(J)(39:39),'(A1)') 'S'
      IF(BANDS(J).EQ.0.) WRITE(TITLEZ(J)(33:39),'(A7)') 'EQUATOR'
   50 CONTINUE
C**** Set the titles for the special zones
      DO 60 J=1,NZS
   60 TITLEZ(8+J)=TITLES(1+J)
      RETURN
      END

      SUBROUTINE CMBINE (AVG,WT, DNEW,NF1,NL1,WT1,WTM, KM, ID)
C****
C**** Bias of new data is removed by subtracting the difference
C**** over the common domain. Then the new data are averaged in.
C****
      COMMON/LIMIT/XBAD,NOVRLP,BIAS(12)
      DIMENSION AVG(KM,*),DNEW(KM,*),WT(KM,*),WT1(KM,*),MISSNG(12)
C**** Loop over months or seasons if appropriate
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
C?*** Find mean bias
C?    WTMNEW=WTM+WT1
C?    BIAS(K)=(WTM*BIAS(K)+WT1*BIASK)/WTMNEW
C?    WTM=WTMNEW
C**** Update period of valid data, averages and weights
      DO 20 N=NF1,NL1
      IF(DNEW(K,N).GE.XBAD) GO TO 20
      WTNEW=WT(K,N)+WT1(K,N)
      AVG(K,N)=(WT(K,N)*AVG(K,N)+WT1(K,N)*(DNEW(K,N)+BIASK))/WTNEW
      WT(K,N)=WTNEW
   20 CONTINUE
      MISSED=MISSED-1
      MISSNG(K)=0
   50 CONTINUE
      IF(MISSED.GT.0) WRITE(6,90) ID,MISSNG
   90 FORMAT(' UNUSED DATA - ID/SUBBOX',I8,12I2)
      RETURN
      END

      SUBROUTINE TAVG (AVG,KM,NYRS, NFB,NLB, JB, DEFLT)
C****
C**** Data are shifted by KM constants such that the KM averages
C**** over the base period (year NFB to NLB) are zero
C****
      COMMON/LIMIT/XBAD,NOVRLP,BIAS(12)
      DIMENSION AVG(KM,*),LEN(12)
      MISSED=KM
      DO 50 K=1,KM
      BIAS(K)=DEFLT
      SUM=0.
      M=0
      DO 10 N=NFB,NLB
      IF(AVG(K,N).GE.XBAD) GO TO 10
      M=M+1
      SUM=SUM+AVG(K,N)
   10 CONTINUE
      LEN(K)=M
      IF(M.EQ.0) GO TO 50
      BIAS(K)=SUM/FLOAT(M)
      MISSED=MISSED-1
   50 CONTINUE
      IF(JB*MISSED.EQ.0) RETURN
C**** If base period is data free, use bias with respect to whole series
      DO 100 K=1,KM
      IF(LEN(K).GT.0) GO TO 100
      WRITE(6,'(''0NO DATA IN BASE PERIOD - MONTH,JB'',3I9)') K,JB
      SUM=0.
      M=0
      DO 60 N=1,NYRS
      IF(AVG(K,N).GE.XBAD) GO TO 60
      M=M+1
      SUM=SUM+AVG(K,N)
   60 CONTINUE
      IF(M.EQ.0) GO TO 100
      BIAS(K)=SUM/FLOAT(M)
  100 CONTINUE
      RETURN
      END

      SUBROUTINE SORT (INDEX,NDIM,LNGTH)
C**** Sorts INDEX and LNGTH such that LNGTH becomes decreasing
      DIMENSION INDEX(NDIM),LNGTH(NDIM)
      DO 10 N=1,NDIM
   10 INDEX(N)=N
      DO 30 N=1,NDIM-1
C**** Find maximum of LNGTH(N),...LNGTH(NDIM)
      NLMAX=N
      DO 20 NN=N+1,NDIM
      IF(LNGTH(NN).GT.LNGTH(NLMAX)) NLMAX=NN
   20 CONTINUE
C**** Switch positions N and NLMAX
      LMAX=LNGTH(NLMAX)
      LNGTH(NLMAX)=LNGTH(N)
      LNGTH(N)=LMAX
      IMAX=INDEX(NLMAX)
      INDEX(NLMAX)=INDEX(N)
   30 INDEX(N)=IMAX
      RETURN
      END

      SUBROUTINE SREAD (IN,LEN,DATA,WT,NG)
C**** Speed read routine for input records
      DIMENSION DATA(LEN),WT(LEN)
      READ (IN) DATA,WT,NG
      RETURN
      END

      SUBROUTINE SWRITE (IOUT,NDIM,AR,WT,TITLE)
      DIMENSION AR(NDIM),WT(NDIM),TITLE(20)
      WRITE(IOUT) AR,WT,TITLE
      RETURN
      END
