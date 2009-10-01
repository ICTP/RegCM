      implicit none
      include 'icbc.param'
      REAL    SST(JX,IY)
      integer NDAY,NMO,NYEAR
      integer number,n
c
      number=8        ! number need be changed for your own SST
c
      open(10,file='SST.RCM',form='unformatted')
      open(20,file='newSST.RCM',form='unformatted')
      do n=1,number
        read(10) NDAY,NMO,NYEAR, SST
        print*,'NDAY=',NDAY,' NMO=',NMO,' NYEAR=',NYEAR

c     write your own code for changing SST if necessary

        write(20) NDAY,NMO,NYEAR, SST
      enddo
c
      stop
      end
