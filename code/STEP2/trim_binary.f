      integer,allocatable,dimension(:) :: idata,idato
      integer info(9)
      character title*80,name*36,nameo*36,filei*12,fileo*17

      do n=1,6
        write(filei,'(a11,i1)') 'work/Ts.bin',n
        open(2,file=filei,form='unformatted')
        write(fileo,'(a12,a5)') filei,'.trim'
        open(3,file=fileo,form='unformatted')

        read(2,end=300) info,title
        i4=info(4)
        allocate (idata(i4),idato(i4))
        read(2,end=200) idato,Lato,Lono,IDo,ihto,nameo,mx,my

        m1o=info(7)
        m2o=info(7)
        do m=1,i4
          if(abs(idato(m)).gt.8000) then
            idato(m)=info(7)
          else
            m2o=m
            if(m1o.eq.info(7)) m1o=m
          end if
        end do
        info(1)=m1o
        info(9)=m2o
        write(3)info,title

   10   m1=info(7)
        m2=info(7)
        read(2,end=100) idata,Lat,Lon,ID,iht,name,mx,my
        do m=1,i4
          if(abs(idata(m)).gt.8000) then
             idata(m)=info(7)
          else
             m2=m
             if(m1.eq.info(7)) m1=m
          end if
        end do
        write(3) (idato(m),m=m1o,m2o),Lato,Lono,IDo,ihto,nameo,m1,m2
        do 20 m=1,i4
   20   idato(m)=idata(m)
        Lato=Lat
        Lono=Lon
        IDo=ID
        ihto=iht
        nameo=name
        m1o=m1
        m2o=m2
        go to 10

  100   write(3) (idato(m),m=m1o,m2o),Lato,Lono,IDo,ihto,nameo,m1,m2
        go to 250
  200   info(1)=info(7)
        info(9)=info(7)
        write(3)info,title
  250   deallocate(idata,idato)
        close (2)
        close (3)
      end do
      stop 0
  300 stop 1
      end
