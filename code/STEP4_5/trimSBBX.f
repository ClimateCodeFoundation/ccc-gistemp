C**** This program trims SBBX files by replacing a totally missing
C**** time series by its first element. The number of elements of the
C**** next time series is added at the BEGINNING of the previous record.
C**** Input file:    unit  10
C****
C****  Record 1: 1,INFOI(2),...,INFOI(8),TITLE             header record
C****  Record 2: AVG(1-->NM),LTS,LTN,LNW,LNE(,NSt,NstMn,Dmin)    sbbx 1
C****  Record 3: AVG(1-->NM),LTS,LTN,LNW,LNE(,NSt,NstMn,Dmin)    sbbx 2
C****            etc.                        not for ocn data
C****
C**** Output file:   unit  11     rearranged !!!
C****
C****  Record 1: NM1,INFO(2),...,INFO(8),TITLE             header record
C****  Record 2: NM2,LTS,LTN,LNW,LNE,NSt,NstMn,Dmin,AVG(1-->NM1)  sbbx 1
C****  Record 3: NM3,LTS,LTN,LNW,LNE,NSt,NstMn,Dmin,AVG(1-->NM2)  sbbx 2
C****            etc.
C****  AVG(1-->NM=INFO(4)) is a full time series, starting at January
C****  of year IYRBEGONFO(6) and ending at December of year IYREND.
C****  NSt   = # of stations contributing to the sub box (0 for ocnfile)
C****  NstMn = # of stations months contributing to the sub box (oc:#ok)
C****  Dmin  = distance of center from nearest contributing station (km)
C****  LTS,LTN is the latitude of the southern,northern edge,
C****  LNW,LNE the longitude of the western,eastern edge of the sub box
C****          in hundredths of degrees (all 4-byte integers).
C****
C****  INFO(1),...,INFO(8)  are 4-byte integers,
C****  TITLE is an 80-byte character string,
C****  AVG are 4-byte reals.
C****      INFO 1 = NM or 1 dep. on whether sbbx 1 has data or not
C****           2 = KQ (quantity flag, see below)
C****           3 = MAVG (time avg flag: 1 - 4 DJF - SON, 5 ANN,
C****                     6 MONTHLY, 7 SEAS, 8 - 19 JAN - DEC  )
C****           4 = NM   (length of each time record)
C****           5 = NM+8 (size of data record length)
C****           6 = IYRBEG (first year of each time record)
C****           7 = flag for missing data
C****           8 = flag for precipitation trace
C****  INFOI(I)=INFO(I) for I=2,3,4,6,7,8   INFOI(5)=NM+7
C****  Missing data are flagged by the real number XBAD=FLOAT(INFO(7))
C****
C****  NOTE: Trimmed Land- and Ocean SBBXfiles have the SAME structure
C****
C?*** Input parameters (# of input files, time period)
      PARAMETER (KM0=12,MONM0=KM0*(3000-1850+1),NSBBX=8000)
      CHARACTER*80 TITLE
      DIMENSION AVG(MONM0),AVGO(MONM0)
      INTEGER INFO(8),INFOO(8),LAT(7),LATO(7)
C****
C**** Read and use the header record of an input file
C****
      open(10,file='work/fort.10',form='unformatted')
      open(11,file='work/fort.11',form='unformatted')
      READ (10) INFO,TITLE
C     KM = the number of time frames per year
      KM=1
      IF(INFO(3).EQ.6) KM=12
      IF(INFO(3).EQ.7) KM=4
      IF(KM.NE.KM0) STOP 'ERROR: CHANGE KM0'
      ML=INFO(4)
      IF(MONM0.LT.ML) THEN
         WRITE(6,'('' SET MONM0 at least TO '',I5)') ML
         STOP 'ERROR: MONM0 NOT OK'
      END IF
      MBAD=INFO(7)
      LAST=INFO(7)
      XBAD=MBAD
C**** Create and write out the header record of the output file
      DO 10 I=2,8
   10 INFOO(I)=INFO(I)
      INFOO(5)=ML+8
      NLAT=INFO(5)-INFO(4)
      LAT(6)=1
      LATO(6)=1
      LATO(5)=0
      LAT(5)=0
      LATO(7)=0
      LAT(7)=0
C****
C**** Loop over NSBBX Sub Boxes
C****
      CALL SREAD(10,ML,AVGO,LATO,NLAT)
      INFOO(1)=1   !  or ML if first SBBX has data
      IF(NLAT.EQ.4) CALL FINDL6(AVGO,ML,XBAD,LATO(6))
      IF(LATO(6).GT.0) INFOO(1)=ML
      WRITE(11) INFOO,TITLE
      DO 300 N=2,NSBBX
      CALL SREAD(10,ML,AVG,LAT,NLAT)
      MLN=1
      IF(NLAT.EQ.4) CALL FINDL6(AVG,ML,XBAD,LAT(6))
      IF(LAT(6).GT.0) MLN=ML
      IF(LATO(6).EQ.0) THEN
         WRITE(11) MLN,LATO,XBAD
      ELSE
         CALL SWRITE(11,ML,AVGO,LATO,MLN)
      END IF
      DO M=1,ML
      AVGO(M)=AVG(M)
      END DO
      DO M=1,7
      LATO(M)=LAT(M)
      END DO
  300 CONTINUE
      MLN=0
      IF(LATO(6).EQ.0) THEN
         WRITE(11) MLN,LATO,XBAD
      ELSE
         CALL SWRITE(11,ML,AVGO,LATO,MLN)
      END IF
      STOP
      END

      SUBROUTINE SWRITE (IOUT,NDIM,ARRAY,LATLON,ML)
      DIMENSION ARRAY(NDIM),LATLON(7)
      WRITE(IOUT) ML,LATLON,ARRAY
      RETURN
      END

      SUBROUTINE SREAD (IIN,NDIM,ARRAY,LATLON,NLAT)
      DIMENSION ARRAY(NDIM),LATLON(NLAT)
      READ (IIN) ARRAY,LATLON
      RETURN
      END

      SUBROUTINE FINDL6 (AV,ML,XBAD,IOK)
      REAL AV(ML)
      IOK=0
      do M=1,ML
      IF(AV(M).NE.XBAD) IOK=IOK+1
      end do
      return
      end

