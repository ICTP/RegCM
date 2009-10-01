      implicit none
      integer im,ny,ncut,nvar
      parameter (im=158,ny=107,ncut=5)
      real*4  a(im,ny)
      integer i,j,n
      open(10,file='SRF.1991mon',form='unformatted'
     &       ,recl=im*ny*4,access='direct')
      open(20,file='SRF.dat',form='unformatted'
     &       ,recl=(im-ncut*2)*(ny-ncut*2)*4,access='direct')
      nvar=27
      do n=1,nvar*12
         read(10,rec=n) a
         write(20,rec=n) ((a(i,j),i=ncut+1,im-ncut),j=ncut+1,ny-ncut)
      enddo
      close(10)
      close(20)
      open(10,file='ATM.1991mon',form='unformatted'
     &       ,recl=im*ny*4,access='direct')
      open(20,file='ATM.dat',form='unformatted'
     &       ,recl=(im-ncut*2)*(ny-ncut*2)*4,access='direct')
      nvar=113
      do n=1,nvar*12
         read(10,rec=n) a
         write(20,rec=n) ((a(i,j),i=ncut+1,im-ncut),j=ncut+1,ny-ncut)
      enddo
      close(10)
      close(20)
      open(10,file='RAD.1991mon',form='unformatted'
     &       ,recl=im*ny*4,access='direct')
      open(20,file='RAD.dat',form='unformatted'
     &       ,recl=(im-ncut*2)*(ny-ncut*2)*4,access='direct')
      nvar=81
      do n=1,nvar*12
         read(10,rec=n) a
         write(20,rec=n) ((a(i,j),i=ncut+1,im-ncut),j=ncut+1,ny-ncut)
      enddo
      close(10)
      close(20)
      open(10,file='OUT_HEAD',form='unformatted'
     &       ,recl=im*ny*4,access='direct')
      open(20,file='HEAD.dat',form='unformatted'
     &       ,recl=(im-ncut*2)*(ny-ncut*2)*4,access='direct')
      do n=2,11
         read(10,rec=n) a
         write(20,rec=n) ((a(i,j),i=ncut+1,im-ncut),j=ncut+1,ny-ncut)
      enddo
      close(10)
      close(20)
      stop
      end
