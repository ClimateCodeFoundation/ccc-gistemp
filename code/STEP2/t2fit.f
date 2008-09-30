      subroutine getfit(nw)
      COMMON/FITCOM/W(900),X(900),F(900),CF(900),DF(900),ZFP(20),ZR(20)
     1             ,FPAR(20),DELTAP,DFSTOP,X0,TMEAN,RMEAN,TFMEAN,RMSFIT
     2             ,YR(900),TS(900),TSFIT(900),RMSP(20),KPFIT(20)
     3             ,KPAR,NXY,NCP,NWCYCL,NITERS,NFIBNO,LISTIT,IDW
      REAL*8 W,X,F,CF,DF,ZFP,ZR,FPAR,DELTAP,DFSTOP,X0,TMEAN,RMEAN,TM3
      REAL*8 TFMEAN,RMSFIT,YR,TS,TSFIT,RMSP

      nhalf=nxy/2

      RMSmin=1.e20

      do n=6,nxy-5
        Xknee=x(n)
        call TREND2(x,F,nxy,Xknee,9999.,2,2,               ! input
     *    sl1,sl2,Yknee,RMS,sl,Y0,RMS0)                   ! output

        if(RMS.lt.RMSmin) then
           RMSmin=RMS
           xmin=Xknee+X0
           fpar(1)=sl1
           fpar(2)=sl2
           fpar(3)=Xknee
           fpar(4)=Yknee
           fpar(5)=sl
           fpar(6)=Y0
           rmsp(1)=RMS/nxy
           rmsp(2)=RMS0/nxy
        end if

      end do

c     write(nw,*) xmin,RMSmin/nxy

      return
      end
