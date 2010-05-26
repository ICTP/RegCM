!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! This code regrids most of RegCM Input/Output (DOMAIN, ICBC, OUT_HEAD, 
!  ATM, RAD, and SRF) files from projected grid to a regular lat-lon grid.
!  Written by Xunqiang Bi, ESP/ICTP
!
      program regrid
      implicit none
      integer iy,jx,kz,nsg,ntr,ibyte
      logical regDOMAIN,regICBC,regOUT_HEAD,regATM,regSRF,regRAD
      integer idate0,idate1,idate2
      character*128 Path_Input,Path_Output
      character*20 DomainName

      real*4  xminlon,xmaxlon,xminlat,xmaxlat,ddeg
      COMMON /WINDOW/ xminlon,xmaxlon,xminlat,xmaxlat,ddeg

      namelist /shareparam/ iy,jx,kz,nsg,ntr,ibyte,Path_Input &
                           ,DomainName,Path_Output
      namelist /dateparam/ idate0,idate1,idate2
      namelist /regridparam/ & 
                regDOMAIN,regICBC,regOUT_HEAD,regATM,regSRF,regRAD  &
               ,xminlon,xmaxlon,xminlat,xmaxlat,ddeg

      xminlon  = -9999.
      xmaxlon  = -9999.
      xminlat  = -9999.
      xmaxlat  = -9999.      
      ddeg     = 0.5

      read(*,shareparam) 
      read(*,dateparam) 
      read(*,regridparam) 

      if(regDOMAIN) call regrid_DOMAIN(iy,jx,kz,ibyte &
                                      ,Path_Input,DomainName)

      if(regICBC) call regrid_ICBC(iy,jx,kz,ibyte,idate0,idate1,idate2 &
                                  ,Path_Input,DomainName)

      if(regOUT_HEAD) call regrid_OUT_HEAD(iy,jx,kz,ibyte,Path_Output)

      if(regATM) call regrid_ATM(iy,jx,kz,ibyte,idate0,idate1,idate2 &
                                ,Path_Output)

      if(regSRF) call regrid_SRF(iy,jx,kz,ibyte,idate0,idate1,idate2 &
                                ,Path_Output)

      if(regRAD) call regrid_RAD(iy,jx,kz,ibyte,idate0,idate1,idate2 &
                                ,Path_Output)

      stop
      end

      subroutine regrid_DOMAIN(iy,jx,kz,ibyte,Path_Input,DomainName)
      implicit none
      integer iy,jx,kz,ibyte
      character*128 Path_Input
      character*20 DomainName
      integer iiy,jjx,kkz
      real*4  dsinm,clat,clon,plat,plon,GRDFAC,ptop
      character*6 iproj
      real*4, allocatable ::  sigma(:)
      integer igrads,ibigend
      real*4  truelatL,truelatH

      integer i,j,k,l,m,n

      real*4  size_2d(4),ddeg
      COMMON /WINDOW/ size_2d,ddeg

      real*4, allocatable ::  o(:,:,:),xlat(:,:),xlon(:,:)
      real*4  xmaxlat,xminlat,xmaxlon,xminlon

      integer mlat,nlat,mlon,nlon
      logical there

      integer MUR,NUR,MUL,NUL,MDR,NDR,MDL,NDL
      real*4  DISTa,DISTb,DISTc,DISTd,AAA,DIST
      real*4  PI
      real*4, allocatable :: out(:,:)
      real*4, allocatable :: olat(:),olon(:)

      integer, allocatable :: I1UR(:,:), J1UR(:,:)
      integer, allocatable :: I1UL(:,:), J1UL(:,:)
      integer, allocatable :: I1DR(:,:), J1DR(:,:)
      integer, allocatable :: I1DL(:,:), J1DL(:,:)
      real*4, allocatable :: D1XT(:,:), D1Xa(:,:), D1Xb(:,:)
      real*4, allocatable ::            D1Xc(:,:), D1Xd(:,:)

      PI = atan(1.)*4.

      allocate(sigma(kz+1))
      allocate(o(jx,iy,9))
      allocate(xlat(jx,iy))
      allocate(xlon(jx,iy))

      inquire(file=trim(Path_Input)//trim(DomainName)//'.INFO' &
             ,exist=there)
      if(.not.there) then
         write(*,*) trim(Path_Input)//trim(DomainName)//'.INFO' &
                  ,' is not avaiable'
         stop
      endif
      open(10,file=trim(Path_Input)//trim(DomainName)//'.INFO' &
          ,form='unformatted',recl=jx*iy*ibyte,access='direct')
      read(10,rec=1) iiy,jjx,kkz,dsinm,clat,clon,plat,plon,GRDFAC  &
                    ,iproj,(sigma(k),k=1,kz+1),ptop,igrads,ibigend &
                    ,truelatL,truelatH
      if(iiy.ne.iy.or.jjx.ne.jx.or.kkz.ne.kz) then
         write(*,*) 'iy,jx,kz in parameter = ',iy,jx,kz
         write(*,*) 'iy,jx,kz in DOMAIN.INFO ',iiy,jjx,kkz
         write(*,*) 'They are not consistent'
         stop
      endif
      read(10,rec=5) xlat         ! xlat
      if(size_2d(3).gt.-9990.) then
         xminlat = size_2d(3)
         xmaxlat = size_2d(4)
      else
         xmaxlat = -1.e10
         xminlat =  1.e10
         do i=2,iy-1
         do j=2,jx-1
            if(xlat(j,i).gt.xmaxlat) xmaxlat=xlat(j,i)
            if(xlat(j,i).lt.xminlat) xminlat=xlat(j,i)
         enddo
         enddo
      endif
      read(10,rec=6) xlon         ! xlon
      if(size_2d(1).gt.-9990.) then
         xminlon=size_2d(1)
         xmaxlon=size_2d(2)
      else
         xmaxlon = -1.e10
         xminlon =  1.e10
         do i=2,iy-1
         do j=2,jx-1
            if(xlon(j,i).gt.xmaxlon) xmaxlon=xlon(j,i)
            if(xlon(j,i).lt.xminlon) xminlon=xlon(j,i)
         enddo
         enddo
      endif
      
      mlat = xmaxlat/ddeg - 1
      nlat = xminlat/ddeg + 1
      mlon = xmaxlon/ddeg - 1
      nlon = xminlon/ddeg + 2

      write(*,*) xminlat,mlat-nlat+1,(nlat-1)*ddeg
      write(*,*) xminlon,mlon-nlon+1,(nlon-1)*ddeg

      open(20,file=trim(Path_Input)//trim(DomainName)//'_LL.dat' &
         ,form='unformatted',recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte &
         ,access='direct')

      allocate(out(nlon:mlon,nlat:mlat))
      allocate(olat(nlat:mlat))
      allocate(olon(nlon:mlon))
      do n=nlat,mlat
         olat(n) = real(n-1)*ddeg
      enddo
      do m=nlon,mlon
         olon(m) = real(m-1)*ddeg
      enddo

      read(10,rec=2)  ((o(j,i,1),j=1,jx),i=1,iy)   ! ht
      read(10,rec=5)  ((o(j,i,2),j=1,jx),i=1,iy)   ! xlat
      read(10,rec=6)  ((o(j,i,3),j=1,jx),i=1,iy)   ! xlon
      read(10,rec=7)  ((o(j,i,4),j=1,jx),i=1,iy)   ! dlat
      read(10,rec=8)  ((o(j,i,5),j=1,jx),i=1,iy)   ! dlon
      read(10,rec=9)  ((o(j,i,6),j=1,jx),i=1,iy)   ! xmap
      read(10,rec=10) ((o(j,i,7),j=1,jx),i=1,iy)   ! dmap
      read(10,rec=11) ((o(j,i,8),j=1,jx),i=1,iy)   ! coriol
      read(10,rec=13) ((o(j,i,9),j=1,jx),i=1,iy)   ! mask

      if(iproj.eq.'NORMER') then
         write(*,*)'regrid: NORMER projection, DOMAIN.INFO'
         do l=1,9
         do n=nlat,mlat
         do m=nlon,mlon
            out(m,n) = -9999.

            do i=2,iy-2
            do j=2,jx-2
               if((xlat(j,i).lt.olat(n).and.      &
                   xlon(j,i).lt.olon(m))          &
             .and.(xlat(j+1,i+1).ge.olat(n).and.  &
                   xlon(j+1,i+1).ge.olon(m))) then
                   out(m,n) = &
       ( o(j,i,l)*(xlon(j+1,i+1)-olon(m))*(xlat(j+1,i+1)-olat(n)) &
       + o(j+1,i,l)*(olon(m)-xlon(j,i+1))*(xlat(j+1,i+1)-olat(n)) &
       + o(j,i+1,l)*(xlon(j+1,i)-olon(m))*(olat(n)-xlat(j+1,i)) &
       + o(j+1,i+1,l)*(olon(m)-xlon(j,i))*(olat(n)-xlat(j,i)) ) &
       / ( (xlon(j+1,i+1)-xlon(j,i))*(xlat(j+1,i+1)-xlat(j,i)) )
               endif
            enddo
            enddo
         enddo
         enddo
         write(20,rec=l) ((out(m,n),m=nlon,mlon),n=nlat,mlat)
         enddo
         close(20)

      else if(iproj.eq.'LAMCON'.or.iproj.eq.'ROTMER') then
         if(iproj.eq.'LAMCON') then
            write(*,*)'regrid: LAMCON projection, DOMAIN.INFO'
         else if(iproj.eq.'ROTMER') then
            write(*,*)'regrid: ROTMER projection, DOMAIN.INFO'
         endif

         allocate(I1UR(nlon:mlon,nlat:mlat))
         allocate(J1UR(nlon:mlon,nlat:mlat))
         allocate(I1UL(nlon:mlon,nlat:mlat))
         allocate(J1UL(nlon:mlon,nlat:mlat))
         allocate(I1DR(nlon:mlon,nlat:mlat))
         allocate(J1DR(nlon:mlon,nlat:mlat))
         allocate(I1DL(nlon:mlon,nlat:mlat))
         allocate(J1DL(nlon:mlon,nlat:mlat))
         allocate(D1XT(nlon:mlon,nlat:mlat))
         allocate(D1Xa(nlon:mlon,nlat:mlat))
         allocate(D1Xb(nlon:mlon,nlat:mlat))
         allocate(D1Xc(nlon:mlon,nlat:mlat))
         allocate(D1Xd(nlon:mlon,nlat:mlat))

         do n=nlat,mlat
         do m=nlon,mlon

            MUR=1000
            NUR=1000
            MUL=1000
            NUL=1000
            MDR=1000
            NDR=1000
            MDL=1000
            NDL=1000
            DISTa=1.E8
            DISTb=1.E8
            DISTc=1.E8
            DISTd=1.E8

            do i=2,iy-1
            do j=2,jx-1
         IF((xlon(j,i).ge.olon(m).and.xlon(j,i)-olon(m).lt.10.) .and. &
            (xlat(j,i).ge.olat(n).and.xlat(j,i)-olat(n).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTa.gt.AAA) then
                DISTa = AAA
                MUR = j
                NUR = i
             endif
         ENDIF
         IF((xlon(j,i).lt.olon(m).and.olon(m)-xlon(j,i).lt.10.) .and. &
            (xlat(j,i).ge.olat(n).and.xlat(j,i)-olat(n).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTb.gt.AAA) then
                DISTb = AAA
                MUL = j
                NUL = i
             endif
         ENDIF
         IF((xlon(j,i).ge.olon(m).and.xlon(j,i)-olon(m).lt.10.) .and. &
            (xlat(j,i).lt.olat(n).and.olat(n)-xlat(j,i).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTc.gt.AAA) then
                DISTc = AAA
                MDR = j
                NDR = i
             endif
         ENDIF
         IF((xlon(j,i).lt.olon(m).and.olon(m)-xlon(j,i).lt.10.) .and. &
            (xlat(j,i).lt.olat(n).and.olat(n)-xlat(j,i).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTd.gt.AAA) then
                DISTd = AAA
                MDL = j
                NDL = i
             endif
         ENDIF
            enddo
            enddo
            
            DIST=amin1(DISTa,DISTb,DISTc,DISTd)
            I1UR(m,n) = MUR
            J1UR(m,n) = NUR
            I1UL(m,n) = MUL
            J1UL(m,n) = NUL
            I1DR(m,n) = MDR
            J1DR(m,n) = NDR
            I1DL(m,n) = MDL
            J1DL(m,n) = NDL
            D1XT(m,n) = DIST
            D1Xa(m,n) = DISTa
            D1Xb(m,n) = DISTb
            D1Xc(m,n) = DISTc
            D1Xd(m,n) = DISTd
         enddo
         enddo

         do l=1,9
         do n=nlat,mlat
         do m=nlon,mlon
            out(m,n) = -9999.

            if(I1UR(m,n).lt.999.and.J1UR(m,n).lt.999.and. &
               I1UL(m,n).lt.999.and.J1UL(m,n).lt.999.and. &
               I1DR(m,n).lt.999.and.J1DR(m,n).lt.999.and. &
               I1DL(m,n).lt.999.and.J1DL(m,n).lt.999) then
               if(D1XT(m,n).gt.0.0001) then
      out(m,n) = ( o(I1UR(m,n),J1UR(m,n),l)/D1Xa(m,n) &
                  +o(I1UL(m,n),J1UL(m,n),l)/D1Xb(m,n) &
                  +o(I1DR(m,n),J1DR(m,n),l)/D1Xc(m,n) &
                  +o(I1DL(m,n),J1DL(m,n),l)/D1Xd(m,n) )  &
          /( 1./D1Xa(m,n)+1./D1Xb(m,n)+1./D1Xc(m,n)+1./D1Xd(m,n) )
               else
                  if(D1Xa(m,n).eq.D1XT(m,n)) then
                     out(m,n) = o(I1UR(m,n),J1UR(m,n),l)
                  else if(D1Xb(m,n).eq.D1XT(m,n)) then
                     out(m,n) = o(I1UL(m,n),J1UL(m,n),l)
                  else if(D1Xc(m,n).eq.D1XT(m,n)) then
                     out(m,n) = o(I1DR(m,n),J1DR(m,n),l)
                  else if(D1Xd(m,n).eq.D1XT(m,n)) then
                     out(m,n) = o(I1DL(m,n),J1DL(m,n),l)
                  endif
               endif
            endif
         enddo
         enddo
         write(20,rec=l) ((out(m,n),m=nlon,mlon),n=nlat,mlat)
         enddo
         close(20)
            
         deallocate(I1UR)
         deallocate(J1UR)
         deallocate(I1UL)
         deallocate(J1UL)
         deallocate(I1DR)
         deallocate(J1DR)
         deallocate(I1DL)
         deallocate(J1DL)
         deallocate(D1XT)
         deallocate(D1Xa)
         deallocate(D1Xb)
         deallocate(D1Xc)
         deallocate(D1Xd)
      endif
      deallocate(out)
      deallocate(olat)
      deallocate(olon)

      deallocate(sigma)
      deallocate(o)
      deallocate(xlat)
      deallocate(xlon)

      if(igrads.eq.1) then
         inquire(file=trim(Path_Input)//trim(DomainName)//'_LL.ctl' &
                ,exist=there)
         if(there) then
         open(31,file=trim(Path_Input)//trim(DomainName)//'_LL.ctl' &
                ,form='formatted',status='replace')
         else
         open(31,file=trim(Path_Input)//trim(DomainName)//'_LL.ctl' &
                ,form='formatted',status='new')
         endif
         write(31,10) '^'//trim(DomainName)//'_LL.dat'
  10  format('dset ',A18)
         write(31,20)
  20  format('title RegCM domain information')
         if(ibigend.eq.1) then
            write(31,30)
  30  format('options big_endian')
         else
            write(31,40)
  40  format('options little_endian')
         endif
         write(31,50)
  50  format('undef -9999.')
         write(31,200) mlon-nlon+1,(nlon-1)*ddeg,ddeg
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
         write(31,210) mlat-nlat+1,(nlat-1)*ddeg,ddeg
 210  format('ydef ',I3,' linear ',f9.4,' ',f9.4)
         write(31,300) 1,1000.
 300  format('zdef ',I1,' levels ',f7.2)
         write(31,400) 1
 400  format('tdef ',I1,' linear 00z01Jan2001 1mo')
         write(31,500) 9
 500  format('vars ',I2)
         write(31,600) 'ht      ','surface elevation          '
         write(31,600) 'xlat    ','latitude  of cross points  '
         write(31,600) 'xlon    ','longitude of cross points  '
         write(31,600) 'dlat    ','latitude  of dot points    '
         write(31,600) 'dlon    ','longitude of dot points    '
         write(31,600) 'xmap    ','map factors of cross points'
         write(31,600) 'dmap    ','map factors of dot points  '
         write(31,600) 'coriol  ','coriol force               '
         write(31,600) 'mask    ','land/sea mask              '
 600  format(A8,'0 99 ',A26)           
         write(31,700)
 700  format('endvars')
         close(31)
      endif
      return
      end

      subroutine regrid_OUT_HEAD(iy,jx,kz,ibyte,Path_Output)
      implicit none
      integer iy,jx,kz,ibyte
      character*128 Path_Output
      integer iiy,jjx,kkz
      integer mdate0,ibltyp,icup,ipptls,iboudy
      real*4  dxsp,ptsp,clat,clon,plat,plon
      real*4  dto,dtb,dtr,dtc
      integer iotyp
      character*6 iproj
      real*4, allocatable ::  sigma(:)

      integer i,j,k,l,m,n
      integer igrads,ibigend

      real*4  size_2d(4),ddeg
      COMMON /WINDOW/ size_2d,ddeg

      real*4, allocatable ::  o(:,:,:),xlat(:,:),xlon(:,:)
      real*4  xmaxlat,xminlat,xmaxlon,xminlon

      integer mlat,nlat,mlon,nlon
      logical there

      integer MUR,NUR,MUL,NUL,MDR,NDR,MDL,NDL
      real*4  DISTa,DISTb,DISTc,DISTd,AAA,DIST
      real*4  PI
      real*4, allocatable :: out(:,:)
      real*4, allocatable :: olat(:),olon(:)

      integer, allocatable :: I1UR(:,:), J1UR(:,:)
      integer, allocatable :: I1UL(:,:), J1UL(:,:)
      integer, allocatable :: I1DR(:,:), J1DR(:,:)
      integer, allocatable :: I1DL(:,:), J1DL(:,:)
      real*4, allocatable :: D1XT(:,:), D1Xa(:,:), D1Xb(:,:)
      real*4, allocatable ::            D1Xc(:,:), D1Xd(:,:)

      igrads = 1
      ibigend= 1

      PI = atan(1.)*4.

      allocate(sigma(kz+1))
      allocate(o(jx-2,iy-2,7))
      allocate(xlat(jx-2,iy-2))
      allocate(xlon(jx-2,iy-2))

      inquire(file=trim(Path_Output)//'OUT_HEAD',exist=there)
      if(.not.there) then
         write(*,*) trim(Path_Output)//'OUT_HEAD',' is not avaiable'
         stop
      endif
      open(10,file=trim(Path_Output)//'OUT_HEAD',form='unformatted' &
             ,recl=(jx-2)*(iy-2)*ibyte,access='direct')
      read(10,rec=1) mdate0,ibltyp,icup,ipptls,iboudy  &
                    ,iiy,jjx,kkz,(sigma(k),k=1,kz+1)   &
                    ,dxsp,ptsp,clat,clon,plat,plon          &
                    ,iproj,dto,dtb,dtr,dtc,iotyp
      if(iiy.ne.iy.or.jjx.ne.jx.or.kkz.ne.kz) then
         write(*,*) 'iy,jx,kz in parameter = ',iy,jx,kz
         write(*,*) 'iy,jx,kz in OUT_HEAD ',iiy,jjx,kkz
         write(*,*) 'They are not consistent'
         stop
      endif
      read(10,rec=6) xlat         ! xlat
      if(size_2d(3).gt.-9990.) then
         xminlat = size_2d(3)
         xmaxlat = size_2d(4)
      else
         xmaxlat = -1.e10
         xminlat =  1.e10
         do i=1,iy-2
         do j=1,jx-2
            if(xlat(j,i).gt.xmaxlat) xmaxlat=xlat(j,i)
            if(xlat(j,i).lt.xminlat) xminlat=xlat(j,i)
         enddo
         enddo
      endif
      read(10,rec=7) xlon         ! xlon
      if(size_2d(1).gt.-9990.) then
         xminlon=size_2d(1)
         xmaxlon=size_2d(2)
      else
         xmaxlon = -1.e10
         xminlon =  1.e10
         do i=1,iy-2
         do j=1,jx-2
            if(xlon(j,i).gt.xmaxlon) xmaxlon=xlon(j,i)
            if(xlon(j,i).lt.xminlon) xminlon=xlon(j,i)
         enddo
         enddo
      endif
      
      mlat = xmaxlat/ddeg - 1
      nlat = xminlat/ddeg + 1
      mlon = xmaxlon/ddeg - 1
      nlon = xminlon/ddeg + 2

      write(*,*) xminlat,mlat-nlat+1,(nlat-1)*ddeg
      write(*,*) xminlon,mlon-nlon+1,(nlon-1)*ddeg

      open(20,file=trim(Path_Output)//'OUTHEAD_LL.dat' &
             ,form='unformatted' &
             ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')

      allocate(out(nlon:mlon,nlat:mlat))
      allocate(olat(nlat:mlat))
      allocate(olon(nlon:mlon))
      do n=nlat,mlat
         olat(n) = real(n-1)*ddeg
      enddo
      do m=nlon,mlon
         olon(m) = real(m-1)*ddeg
      enddo

      read(10,rec=2)  ((o(j,i,1),j=1,jx-2),i=1,iy-2)   ! ht
      read(10,rec=6)  ((o(j,i,2),j=1,jx-2),i=1,iy-2)   ! xlat
      read(10,rec=7)  ((o(j,i,3),j=1,jx-2),i=1,iy-2)   ! xlon
      read(10,rec=8)  ((o(j,i,4),j=1,jx-2),i=1,iy-2)   ! xmap
      read(10,rec=9)  ((o(j,i,5),j=1,jx-2),i=1,iy-2)   ! dmap
      read(10,rec=10) ((o(j,i,6),j=1,jx-2),i=1,iy-2)   ! coriol
      read(10,rec=11) ((o(j,i,7),j=1,jx-2),i=1,iy-2)   ! mask

      if(iproj.eq.'NORMER') then
         write(*,*)'regrid: NORMER projection, OUT_HEAD'
         do l=1,7
         do n=nlat,mlat
         do m=nlon,mlon
            out(m,n) = -9999.

            do i=1,iy-3
            do j=1,jx-3
               if((xlat(j,i).lt.olat(n).and.      &
                   xlon(j,i).lt.olon(m))          &
             .and.(xlat(j+1,i+1).ge.olat(n).and.  &
                   xlon(j+1,i+1).ge.olon(m))) then
                   out(m,n) = &
       ( o(j,i,l)*(xlon(j+1,i+1)-olon(m))*(xlat(j+1,i+1)-olat(n)) &
       + o(j+1,i,l)*(olon(m)-xlon(j,i+1))*(xlat(j+1,i+1)-olat(n)) &
       + o(j,i+1,l)*(xlon(j+1,i)-olon(m))*(olat(n)-xlat(j+1,i)) &
       + o(j+1,i+1,l)*(olon(m)-xlon(j,i))*(olat(n)-xlat(j,i)) ) &
       / ( (xlon(j+1,i+1)-xlon(j,i))*(xlat(j+1,i+1)-xlat(j,i)) )
               endif
            enddo
            enddo
         enddo
         enddo
         write(20,rec=l) ((out(m,n),m=nlon,mlon),n=nlat,mlat)
         enddo
         close(20)

      else if(iproj.eq.'LAMCON'.or.iproj.eq.'ROTMER') then
         if(iproj.eq.'LAMCON') then
            write(*,*)'regrid: LAMCON projection, OUT_HEAD'
         else if(iproj.eq.'ROTMER') then
            write(*,*)'regrid: ROTMER projection, OUT_HEAD'
         endif

         allocate(I1UR(nlon:mlon,nlat:mlat))
         allocate(J1UR(nlon:mlon,nlat:mlat))
         allocate(I1UL(nlon:mlon,nlat:mlat))
         allocate(J1UL(nlon:mlon,nlat:mlat))
         allocate(I1DR(nlon:mlon,nlat:mlat))
         allocate(J1DR(nlon:mlon,nlat:mlat))
         allocate(I1DL(nlon:mlon,nlat:mlat))
         allocate(J1DL(nlon:mlon,nlat:mlat))
         allocate(D1XT(nlon:mlon,nlat:mlat))
         allocate(D1Xa(nlon:mlon,nlat:mlat))
         allocate(D1Xb(nlon:mlon,nlat:mlat))
         allocate(D1Xc(nlon:mlon,nlat:mlat))
         allocate(D1Xd(nlon:mlon,nlat:mlat))

         do n=nlat,mlat
         do m=nlon,mlon

            MUR=1000
            NUR=1000
            MUL=1000
            NUL=1000
            MDR=1000
            NDR=1000
            MDL=1000
            NDL=1000
            DISTa=1.E8
            DISTb=1.E8
            DISTc=1.E8
            DISTd=1.E8

            do i=1,iy-2
            do j=1,jx-2
         IF((xlon(j,i).ge.olon(m).and.xlon(j,i)-olon(m).lt.10.) .and. &
            (xlat(j,i).ge.olat(n).and.xlat(j,i)-olat(n).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTa.gt.AAA) then
                DISTa = AAA
                MUR = j
                NUR = i
             endif
         ENDIF
         IF((xlon(j,i).lt.olon(m).and.olon(m)-xlon(j,i).lt.10.) .and. &
            (xlat(j,i).ge.olat(n).and.xlat(j,i)-olat(n).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTb.gt.AAA) then
                DISTb = AAA
                MUL = j
                NUL = i
             endif
         ENDIF
         IF((xlon(j,i).ge.olon(m).and.xlon(j,i)-olon(m).lt.10.) .and. &
            (xlat(j,i).lt.olat(n).and.olat(n)-xlat(j,i).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTc.gt.AAA) then
                DISTc = AAA
                MDR = j
                NDR = i
             endif
         ENDIF
         IF((xlon(j,i).lt.olon(m).and.olon(m)-xlon(j,i).lt.10.) .and. &
            (xlat(j,i).lt.olat(n).and.olat(n)-xlat(j,i).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTd.gt.AAA) then
                DISTd = AAA
                MDL = j
                NDL = i
             endif
         ENDIF
            enddo
            enddo
            
            DIST=amin1(DISTa,DISTb,DISTc,DISTd)
            I1UR(m,n) = MUR
            J1UR(m,n) = NUR
            I1UL(m,n) = MUL
            J1UL(m,n) = NUL
            I1DR(m,n) = MDR
            J1DR(m,n) = NDR
            I1DL(m,n) = MDL
            J1DL(m,n) = NDL
            D1XT(m,n) = DIST
            D1Xa(m,n) = DISTa
            D1Xb(m,n) = DISTb
            D1Xc(m,n) = DISTc
            D1Xd(m,n) = DISTd
         enddo
         enddo

         do l=1,7
         do n=nlat,mlat
         do m=nlon,mlon
            out(m,n) = -9999.

            if(I1UR(m,n).lt.999.and.J1UR(m,n).lt.999.and. &
               I1UL(m,n).lt.999.and.J1UL(m,n).lt.999.and. &
               I1DR(m,n).lt.999.and.J1DR(m,n).lt.999.and. &
               I1DL(m,n).lt.999.and.J1DL(m,n).lt.999) then
               if(D1XT(m,n).gt.0.0001) then
      out(m,n) = ( o(I1UR(m,n),J1UR(m,n),l)/D1Xa(m,n) &
                  +o(I1UL(m,n),J1UL(m,n),l)/D1Xb(m,n) &
                  +o(I1DR(m,n),J1DR(m,n),l)/D1Xc(m,n) &
                  +o(I1DL(m,n),J1DL(m,n),l)/D1Xd(m,n) )  &
          /( 1./D1Xa(m,n)+1./D1Xb(m,n)+1./D1Xc(m,n)+1./D1Xd(m,n) )
               else
                  if(D1Xa(m,n).eq.D1XT(m,n)) then
                     out(m,n) = o(I1UR(m,n),J1UR(m,n),l)
                  else if(D1Xb(m,n).eq.D1XT(m,n)) then
                     out(m,n) = o(I1UL(m,n),J1UL(m,n),l)
                  else if(D1Xc(m,n).eq.D1XT(m,n)) then
                     out(m,n) = o(I1DR(m,n),J1DR(m,n),l)
                  else if(D1Xd(m,n).eq.D1XT(m,n)) then
                     out(m,n) = o(I1DL(m,n),J1DL(m,n),l)
                  endif
               endif
            endif
         enddo
         enddo
         write(20,rec=l) ((out(m,n),m=nlon,mlon),n=nlat,mlat)
         enddo
         close(20)
            
         deallocate(I1UR)
         deallocate(J1UR)
         deallocate(I1UL)
         deallocate(J1UL)
         deallocate(I1DR)
         deallocate(J1DR)
         deallocate(I1DL)
         deallocate(J1DL)
         deallocate(D1XT)
         deallocate(D1Xa)
         deallocate(D1Xb)
         deallocate(D1Xc)
         deallocate(D1Xd)
      endif
      deallocate(out)
      deallocate(olat)
      deallocate(olon)

      deallocate(sigma)
      deallocate(o)
      deallocate(xlat)
      deallocate(xlon)

      if(igrads.eq.1) then
         inquire(file=trim(Path_Output)//'OUTHEAD_LL.ctl' &
                ,exist=there)
         if(there) then
         open(31,file=trim(Path_Output)//'OUTHEAD_LL.ctl' &
                ,form='formatted',status='replace')
         else
         open(31,file=trim(Path_Output)//'OUTHEAD_LL.ctl' &
                ,form='formatted',status='new')
         endif
         write(31,10)
  10  format('dset ^OUTHEAD_LL.dat')
         write(31,20)
  20  format('title RegCM domain information')
         if(ibigend.eq.1) then
            write(31,30)
  30  format('options big_endian')
         else
            write(31,40)
  40  format('options little_endian')
         endif
         write(31,50)
  50  format('undef -9999.')
         write(31,200) mlon-nlon+1,(nlon-1)*ddeg,ddeg
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
         write(31,210) mlat-nlat+1,(nlat-1)*ddeg,ddeg
 210  format('ydef ',I3,' linear ',f9.4,' ',f9.4)
         write(31,300) 1,1000.
 300  format('zdef ',I1,' levels ',f7.2)
         write(31,400) 1
 400  format('tdef ',I1,' linear 00z01Jan2001 1mo')
         write(31,500) 7
 500  format('vars ',I2)
         write(31,600) 'ht      ','surface elevation          '
         write(31,600) 'xlat    ','latitude  of cross points  '
         write(31,600) 'xlon    ','longitude of cross points  '
         write(31,600) 'xmap    ','map factors of cross points'
         write(31,600) 'dmap    ','map factors of dot points  '
         write(31,600) 'coriol  ','coriol force               '
         write(31,600) 'mask    ','land/sea mask              '
 600  format(A8,'0 99 ',A26)           
         write(31,700)
 700  format('endvars')
         close(31)
      endif
      return
      end

      subroutine regrid_ICBC(iy,jx,kz,ibyte,idate0,idate1,idate2 &
                            ,Path_Input,DomainName)
      implicit none
      integer iy,jx,kz,ibyte,idate0,idate1,idate2
      character*128 Path_Input
      character*20 DomainName
      integer iiy,jjx,kkz
      real*4  dsinm,clat,clon,plat,plon,GRDFAC,ptop
      character*6 iproj
      real*4, allocatable ::  sigma(:)
      integer igrads,ibigend
      real*4  truelatL,truelatH
      character*4 :: chy(1941:2100)
      data chy/'1941','1942','1943','1944','1945','1946','1947','1948',&
        '1949','1950','1951','1952','1953','1954','1955','1956','1957',&
        '1958','1959','1960','1961','1962','1963','1964','1965','1966',&
        '1967','1968','1969','1970','1971','1972','1973','1974','1975',&
        '1976','1977','1978','1979','1980','1981','1982','1983','1984',&
        '1985','1986','1987','1988','1989','1990','1991','1992','1993',&
        '1994','1995','1996','1997','1998','1999','2000','2001','2002',&
        '2003','2004','2005','2006','2007','2008','2009','2010','2011',&
        '2012','2013','2014','2015','2016','2017','2018','2019','2020',&
        '2021','2022','2023','2024','2025','2026','2027','2028','2029',&
        '2030','2031','2032','2033','2034','2035','2036','2037','2038',&
        '2039','2040','2041','2042','2043','2044','2045','2046','2047',&
        '2048','2049','2050','2051','2052','2053','2054','2055','2056',&
        '2057','2058','2059','2060','2061','2062','2063','2064','2065',&
        '2066','2067','2068','2069','2070','2071','2072','2073','2074',&
        '2075','2076','2077','2078','2079','2080','2081','2082','2083',&
        '2084','2085','2086','2087','2088','2089','2090','2091','2092',&
        '2093','2094','2095','2096','2097','2098','2099','2100'/
      character*2 chm(12)
      data chm/'01','02','03','04','05','06','07','08','09','10', &
               '11','12'/
      character*3 chmc(12)
      data chmc/'jan','feb','mar','apr','may','jun'  &
               ,'jul','aug','sep','oct','nov','dec'/
      character*14 filein
      character*18 fileout
      integer ntype,nfile,nyear,month,n_slice,mrec,nrec,nnn

      integer i,j,k,l,m,n

      real*4  size_2d(4),ddeg
      COMMON /WINDOW/ size_2d,ddeg

      real*4, allocatable ::  o(:,:,:),xlat(:,:),xlon(:,:)
      real*4, allocatable ::           dlat(:,:),dlon(:,:)
      real*4  xmaxlat,xminlat,xmaxlon,xminlon

      integer mlat,nlat,mlon,nlon
      integer n_month
      logical there

      integer MUR,NUR,MUL,NUL,MDR,NDR,MDL,NDL
      real*4  DISTa,DISTb,DISTc,DISTd,AAA,DIST
      integer IDATE
      real*4  PI,PIR180
      real*4  POLLAM,POLPHI,POLCPHI,POLSPHI,ZPHI,ZRLA,ZRLAP
      real*4  ZARG1,ZARG2,ZNORM,SINDEL,COSDEL,US,VS,X
      real*4  SIGN0,GRIDFC
      real*4, allocatable :: out(:,:)
      real*4, allocatable :: olat(:),olon(:)

      integer, allocatable :: I1UR(:,:), J1UR(:,:)
      integer, allocatable :: I1UL(:,:), J1UL(:,:)
      integer, allocatable :: I1DR(:,:), J1DR(:,:)
      integer, allocatable :: I1DL(:,:), J1DL(:,:)
      real*4, allocatable :: D1XT(:,:), D1Xa(:,:), D1Xb(:,:)
      real*4, allocatable ::            D1Xc(:,:), D1Xd(:,:)

      integer, allocatable :: I2UR(:,:), J2UR(:,:)
      integer, allocatable :: I2UL(:,:), J2UL(:,:)
      integer, allocatable :: I2DR(:,:), J2DR(:,:)
      integer, allocatable :: I2DL(:,:), J2DL(:,:)
      real*4, allocatable :: D2XT(:,:), D2Xa(:,:), D2Xb(:,:)
      real*4, allocatable ::            D2Xc(:,:), D2Xd(:,:)

      PI = atan(1.)*4.
      PIR180 = atan(1.)/45.

      allocate(sigma(kz+1))
      allocate(o(jx,iy,kz*4+2))
      allocate(xlat(jx,iy))
      allocate(xlon(jx,iy))
      allocate(dlat(jx,iy))
      allocate(dlon(jx,iy))

      inquire(file=trim(Path_Input)//trim(DomainName)//'.INFO' &
             ,exist=there)
      if(.not.there) then
              write(*,*) trim(Path_Input)//trim(DomainName)//'.INFO' &
                        ,' is not avaiable'
         stop
      endif
      open(10,file=trim(Path_Input)//trim(DomainName)//'.INFO'  &
             ,form='unformatted',recl=jx*iy*ibyte,access='direct')
      read(10,rec=1) iiy,jjx,kkz,dsinm,clat,clon,plat,plon,GRDFAC  &
                    ,iproj,(sigma(k),k=1,kz+1),ptop,igrads,ibigend &
                    ,truelatL,truelatH
      if(iiy.ne.iy.or.jjx.ne.jx.or.kkz.ne.kz) then
         write(*,*) 'iy,jx,kz in parameter = ',iy,jx,kz
         write(*,*) 'iy,jx,kz in DOMAIN.INFO ',iiy,jjx,kkz
         write(*,*) 'They are not consistent'
         stop
      endif

      read(10,rec=5) xlat         ! xlat
      if(size_2d(3).gt.-9990.) then
         xminlat = size_2d(3)
         xmaxlat = size_2d(4)
      else
         xmaxlat = -1.e10
         xminlat =  1.e10
         do i=2,iy-2
         do j=2,jx-2
            if(xlat(j,i).gt.xmaxlat) xmaxlat=xlat(j,i)
            if(xlat(j,i).lt.xminlat) xminlat=xlat(j,i)
         enddo
         enddo
      endif
      read(10,rec=6) xlon         ! xlon
      if(size_2d(1).gt.-9990.) then
         xminlon=size_2d(1)
         xmaxlon=size_2d(2)
      else
         xmaxlon = -1.e10
         xminlon =  1.e10
         do i=2,iy-2
         do j=2,jx-2
            if(xlon(j,i).gt.xmaxlon) xmaxlon=xlon(j,i)
            if(xlon(j,i).lt.xminlon) xminlon=xlon(j,i)
         enddo
         enddo
      endif
      read(10,rec=7) dlat         ! dlat
      read(10,rec=8) dlon         ! dlon
      close(10)
      
      mlat = xmaxlat/ddeg - 1
      nlat = xminlat/ddeg + 1
      mlon = xmaxlon/ddeg - 1
      nlon = xminlon/ddeg + 2

      write(*,*) xminlat,mlat-nlat+1,(nlat-1)*ddeg
      write(*,*) xminlon,mlon-nlon+1,(nlon-1)*ddeg

      allocate(out(nlon:mlon,nlat:mlat))
      allocate(olat(nlat:mlat))
      allocate(olon(nlon:mlon))
      do n=nlat,mlat
         olat(n) = real(n-1)*ddeg
      enddo
      do m=nlon,mlon
         olon(m) = real(m-1)*ddeg
      enddo

      if(idate1.gt.1940010100) then        ! original file
         n_month = (idate2/1000000-idate1/1000000)*12 &
                 + (mod(idate2/10000,100)-mod(idate1/10000,100))   
         if(mod(idate2,10000).gt.0100) n_month = n_month+1
         ntype = 0
      else if(idate1.gt.19400101) then     ! daily mean file
         n_month = (idate2/10000-idate1/10000)*12 &
                 + (mod(idate2/100,100)-mod(idate1/100,100))   
         if(mod(idate2,100).gt.01) n_month = n_month+1
         ntype = 1
      else if(idate1.gt.194001) then       ! monthly mean file
         n_month = (idate2/100-idate1/100)*12 &
                 + (mod(idate2,100)-mod(idate1,100))+1 
         ntype = 2
      endif

      if(iproj.eq.'NORMER') then
         do nfile=1,n_month
            if(ntype.eq.0) then
               nyear = idate1/1000000
               month = idate1/10000-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(nfile.eq.1.or.month.eq.1) then
                  write(*,*)'regrid: NORMER ICBC Orig.',nyear,month
               else
                  write(*,*)'                         ',nyear,month
               endif

               if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
                  month.eq.7.or.month.eq.8.or.month.eq.10.or. &
                  month.eq.12) then
                  n_slice = 31*24/6+1
               else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                       month.eq.11) then
                  n_slice = 30*24/6+1
               else
                  n_slice = 28*24/6+1
                  if(mod(nyear,4).eq.0) n_slice = 29*24/6+1
                  if(mod(nyear,100).eq.0) n_slice = 28*24/6+1
                  if(mod(nyear,400).eq.0) n_slice = 29*24/6+1
               endif
               if(nfile.eq.n_month.and.mod(idate2/100-1,100).ne.0) &
               n_slice=min(n_slice,nint(mod(idate2/100-1,100)*24./6))+1
               filein = 'ICBC'//chy(nyear)//chm(month)//'0100'
               fileout= 'ICBC_LL.'//chy(nyear)//chm(month)//'0100'

               inquire(file= &
                    trim(Path_Input)//trim(DomainName)//'_'//filein &
                      ,exist=there)
               if(.not.there) then
                  write(*,*) &
                    trim(Path_Input)//trim(DomainName)//'_'//filein &
                        ,' is not avaiable'
                  stop
               endif
               open(10,file= &
                    trim(Path_Input)//trim(DomainName)//'_'//filein  &
                   ,form='unformatted',recl=iy*jx*ibyte,access='direct')
               mrec = 0
               open(20,file= &
                    trim(Path_Input)//trim(DomainName)//'_'//fileout &
                   ,form='unformatted' &
               ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
               nrec = 0
            else if(ntype.eq.1) then
               nyear = idate1/10000
               month = idate1/100-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(nfile.eq.1.or.month.eq.1) then
                  write(*,*)'regrid: NORMER ICBC Daily',nyear,month
               else
                  write(*,*)'                         ',nyear,month
               endif

               if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
                  month.eq.7.or.month.eq.8.or.month.eq.10.or. &
                  month.eq.12) then
                  n_slice = 31
               else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                       month.eq.11) then
                  n_slice = 30
               else
                  n_slice = 28
                  if(mod(nyear,4).eq.0) n_slice = 29
                  if(mod(nyear,100).eq.0) n_slice = 28
                  if(mod(nyear,400).eq.0) n_slice = 29
               endif
               if(nfile.eq.n_month.and.mod(idate2-1,100).ne.0) &
               n_slice=min(n_slice,mod(idate2-1,100)+1)
               filein = 'ICBC'//chy(nyear)//chm(month)//'01'
               fileout= 'ICBC_LL.'//chy(nyear)//chm(month)//'01'

               inquire(file= &
                 trim(Path_Input)//trim(DomainName)//'_'//filein(1:12) &
                      ,exist=there)
               if(.not.there) then
                  write(*,*) &
                 trim(Path_Input)//trim(DomainName)//'_'//filein(1:12) &
                        ,' is not avaiable'
                  stop
               endif
               open(10,file= &
                 trim(Path_Input)//trim(DomainName)//'_'//filein(1:12) &
                   ,form='unformatted' &
                   ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
               mrec = 0
               open(20,file= &
                trim(Path_Input)//trim(DomainName)//'_'//fileout(1:16) &
                ,form='unformatted' &
                ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
               nrec = 0
            else if(ntype.eq.2) then
               nyear = idate1/100
               month = idate1-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(month.eq.1.or.nfile.eq.1) then
                  write(*,*)'regrid: NORMER ICBC Month',nyear,month
               else
                  write(*,*)'                         ',nyear,month
               endif

               n_slice = 1

               if(nfile.eq.1) then
                  filein = 'ICBC'
                  fileout= 'ICBC_LL'

                  inquire(file= &
          trim(Path_Input)//trim(DomainName)//'_'//filein(1:4)//'.mon' &
                         ,exist=there)
                  if(.not.there) then
                  write(*,*) &
          trim(Path_Input)//trim(DomainName)//'_'//filein(1:7)//'.mon' &
                        ,' is not avaiable'
                     stop
                  endif
                  open(10,file= &
          trim(Path_Input)//trim(DomainName)//'_'//filein(1:7)//'.mon' &
                      ,form='unformatted' &
                      ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
                  nrec=0
                  open(20,file= &
         trim(Path_Input)//trim(DomainName)//'_'//fileout(1:7)//'.mon' &
                ,form='unformatted' &
                ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
                  mrec=0
               endif
            endif

            do nnn=1,n_slice
               if(ntype.eq.0) then
                  mrec = mrec+1
                  read(10,rec=mrec) IDATE
!                 write(*,*) 'IDATE = ',IDATE
                  do l=1,kz*4+2
                     mrec = mrec+1
                     read(10,rec=mrec) ((o(j,i,l),j=1,jx),i=1,iy)
                  enddo
               else
                  do l=1,kz*4+2
                     mrec = mrec+1
                     read(10,rec=mrec) ((o(j,i,l),j=2,jx-1),i=2,iy-1)
                  enddo
               endif
               do l=1,kz*2
                  do n=nlat,mlat
                  do m=nlon,mlon
                     out(m,n) = -9999.

                     do i=2,iy-2
                     do j=2,jx-2
                        if((dlat(j,i).lt.olat(n).and.      &
                            dlon(j,i).lt.olon(m))          &
                      .and.(dlat(j+1,i+1).ge.olat(n).and.  &
                            dlon(j+1,i+1).ge.olon(m))) then
                            out(m,n) = &
       ( o(j,i,l)*(dlon(j+1,i+1)-olon(m))*(dlat(j+1,i+1)-olat(n)) &
       + o(j+1,i,l)*(olon(m)-dlon(j,i+1))*(dlat(j+1,i+1)-olat(n)) &
       + o(j,i+1,l)*(dlon(j+1,i)-olon(m))*(olat(n)-dlat(j+1,i)) &
       + o(j+1,i+1,l)*(olon(m)-dlon(j,i))*(olat(n)-dlat(j,i)) ) &
       / ( (dlon(j+1,i+1)-dlon(j,i))*(dlat(j+1,i+1)-dlat(j,i)) )
                        endif
                     enddo
                     enddo
                  enddo
                  enddo
                  nrec = nrec+1
                  write(20,rec=nrec)((out(m,n),m=nlon,mlon),n=nlat,mlat)
               enddo
               do l=kz*2+1,kz*4+2
                  do n=nlat,mlat
                  do m=nlon,mlon
                     out(m,n) = -9999.

                     do i=2,iy-2
                     do j=2,jx-2
                        if((xlat(j,i).lt.olat(n).and.      &
                            xlon(j,i).lt.olon(m))          &
                      .and.(xlat(j+1,i+1).ge.olat(n).and.  &
                            xlon(j+1,i+1).ge.olon(m))) then
                            out(m,n) = &
       ( o(j,i,l)*(xlon(j+1,i+1)-olon(m))*(xlat(j+1,i+1)-olat(n)) &
       + o(j+1,i,l)*(olon(m)-xlon(j,i+1))*(xlat(j+1,i+1)-olat(n)) &
       + o(j,i+1,l)*(xlon(j+1,i)-olon(m))*(olat(n)-xlat(j+1,i)) &
       + o(j+1,i+1,l)*(olon(m)-xlon(j,i))*(olat(n)-xlat(j,i)) ) &
       / ( (xlon(j+1,i+1)-xlon(j,i))*(xlat(j+1,i+1)-xlat(j,i)) )
                        endif
                     enddo
                     enddo
                  enddo
                  enddo
                  nrec = nrec+1
                  write(20,rec=nrec)((out(m,n),m=nlon,mlon),n=nlat,mlat)
               enddo
            enddo
            if(.not.ntype.eq.2) close(20)
            if(.not.ntype.eq.2) close(10)
            if(igrads.eq.1.and.(ntype.eq.0.or.ntype.eq.1.or. &
               (ntype.eq.2.and.nfile.eq.1))) then
               if(ntype.eq.0) then
                  inquire(file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:18)//'.ctl' &
                         ,exist=there)
                  if(there) then
                     open(31,file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:18)//'.ctl' &
                      ,form='formatted',status='replace')
                  else
                     open(31,file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:18)//'.ctl' &
                      ,form='formatted',status='new')
                  endif
                  write(31,10) '^'//trim(DomainName)//'_'//fileout(1:18)
               else if(ntype.eq.1) then
                  inquire(file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:16)//'.ctl' &
                         ,exist=there)
                  if(there) then
                     open(31,file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:16)//'.ctl' &
                      ,form='formatted',status='replace')
                  else
                     open(31,file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:16)//'.ctl' &
                      ,form='formatted',status='new')
                  endif
                  write(31,11) '^'//trim(DomainName)//'_'//fileout(1:16)
               else if(ntype.eq.2.and.nfile.eq.1) then
                  inquire(file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:7)//'.ctl' &
                         ,exist=there)
                  if(there) then
                     open(31,file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:7)//'.ctl' &
                      ,form='formatted',status='replace')
                  else
                     open(31,file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:7)//'.ctl' &
                      ,form='formatted',status='new')
                  endif
          write(31,12) '^'//trim(DomainName)//'_'//fileout(1:7)//'.mon'
               endif
               write(31,20)
               if(ibigend.eq.1) then
                  write(31,30)
               else
                  write(31,40)
               endif
               write(31,50)
               write(31,200) mlon-nlon+1,(nlon-1)*ddeg,ddeg
               write(31,210) mlat-nlat+1,(nlat-1)*ddeg,ddeg
               write(31,300) kz, &
       ((1013.25-ptop*10.)*(sigma(k)+sigma(k+1))*0.5+ptop*10.,k=kz,1,-1)
               if(ntype.eq.0) then
                  write(31,400) n_slice,0,'01',chmc(month),nyear,6
               else if(ntype.eq.1) then
                  write(31,401) n_slice,'01',chmc(month),nyear
               else if(ntype.eq.2) then
                  write(31,402) n_month,'16',chmc(month),nyear
               endif
               write(31,500) 6
               write(31,650) 'u       ',kz,'westerly wind (m/s)        '
               write(31,650) 'v       ',kz,'southerly wind (m/s)       '
               write(31,650) 't       ',kz,'air temperature (degree)   '
               write(31,650) 'q       ',kz,'water vapor mixing ratio   '
               write(31,600) 'ps      ',   'surface pressure (hPa)     '
               write(31,600) 'tas     ',   'surface air temperature, K '
               write(31,700)
               close(31)
            endif
         enddo

      else if(iproj.eq.'LAMCON'.or.iproj.eq.'ROTMER') then

         allocate(I1UR(nlon:mlon,nlat:mlat))
         allocate(J1UR(nlon:mlon,nlat:mlat))
         allocate(I1UL(nlon:mlon,nlat:mlat))
         allocate(J1UL(nlon:mlon,nlat:mlat))
         allocate(I1DR(nlon:mlon,nlat:mlat))
         allocate(J1DR(nlon:mlon,nlat:mlat))
         allocate(I1DL(nlon:mlon,nlat:mlat))
         allocate(J1DL(nlon:mlon,nlat:mlat))
         allocate(D1XT(nlon:mlon,nlat:mlat))
         allocate(D1Xa(nlon:mlon,nlat:mlat))
         allocate(D1Xb(nlon:mlon,nlat:mlat))
         allocate(D1Xc(nlon:mlon,nlat:mlat))
         allocate(D1Xd(nlon:mlon,nlat:mlat))

         do n=nlat,mlat
         do m=nlon,mlon

            MUR=1000
            NUR=1000
            MUL=1000
            NUL=1000
            MDR=1000
            NDR=1000
            MDL=1000
            NDL=1000
            DISTa=1.E8
            DISTb=1.E8
            DISTc=1.E8
            DISTd=1.E8

            do i=2,iy-1
            do j=2,jx-1
         IF((xlon(j,i).ge.olon(m).and.xlon(j,i)-olon(m).lt.10.) .and. &
            (xlat(j,i).ge.olat(n).and.xlat(j,i)-olat(n).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTa.gt.AAA) then
                DISTa = AAA
                MUR = j
                NUR = i
             endif
         ENDIF
         IF((xlon(j,i).lt.olon(m).and.olon(m)-xlon(j,i).lt.10.) .and. &
            (xlat(j,i).ge.olat(n).and.xlat(j,i)-olat(n).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTb.gt.AAA) then
                DISTb = AAA
                MUL = j
                NUL = i
             endif
         ENDIF
         IF((xlon(j,i).ge.olon(m).and.xlon(j,i)-olon(m).lt.10.) .and. &
            (xlat(j,i).lt.olat(n).and.olat(n)-xlat(j,i).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTc.gt.AAA) then
                DISTc = AAA
                MDR = j
                NDR = i
             endif
         ENDIF
         IF((xlon(j,i).lt.olon(m).and.olon(m)-xlon(j,i).lt.10.) .and. &
            (xlat(j,i).lt.olat(n).and.olat(n)-xlat(j,i).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTd.gt.AAA) then
                DISTd = AAA
                MDL = j
                NDL = i
             endif
         ENDIF
            enddo
            enddo
            
            DIST=amin1(DISTa,DISTb,DISTc,DISTd)
            I1UR(m,n) = MUR
            J1UR(m,n) = NUR
            I1UL(m,n) = MUL
            J1UL(m,n) = NUL
            I1DR(m,n) = MDR
            J1DR(m,n) = NDR
            I1DL(m,n) = MDL
            J1DL(m,n) = NDL
            D1XT(m,n) = DIST
            D1Xa(m,n) = DISTa
            D1Xb(m,n) = DISTb
            D1Xc(m,n) = DISTc
            D1Xd(m,n) = DISTd
         enddo
         enddo

         allocate(I2UR(nlon:mlon,nlat:mlat))
         allocate(J2UR(nlon:mlon,nlat:mlat))
         allocate(I2UL(nlon:mlon,nlat:mlat))
         allocate(J2UL(nlon:mlon,nlat:mlat))
         allocate(I2DR(nlon:mlon,nlat:mlat))
         allocate(J2DR(nlon:mlon,nlat:mlat))
         allocate(I2DL(nlon:mlon,nlat:mlat))
         allocate(J2DL(nlon:mlon,nlat:mlat))
         allocate(D2XT(nlon:mlon,nlat:mlat))
         allocate(D2Xa(nlon:mlon,nlat:mlat))
         allocate(D2Xb(nlon:mlon,nlat:mlat))
         allocate(D2Xc(nlon:mlon,nlat:mlat))
         allocate(D2Xd(nlon:mlon,nlat:mlat))

         do n=nlat,mlat
         do m=nlon,mlon

            MUR=1000
            NUR=1000
            MUL=1000
            NUL=1000
            MDR=1000
            NDR=1000
            MDL=1000
            NDL=1000
            DISTa=1.E8
            DISTb=1.E8
            DISTc=1.E8
            DISTd=1.E8

            do i=2,iy-1
            do j=2,jx-1
         IF((dlon(j,i).ge.olon(m).and.dlon(j,i)-olon(m).lt.10.) .and. &
            (dlat(j,i).ge.olat(n).and.dlat(j,i)-olat(n).lt.10.)) then
             AAA = ((dlon(j,i)-olon(m)) &
                 *cos((dlat(j,i)+olat(n))/360.*pi))**2 &
                 +(dlat(j,i)-olat(n))**2
             if(DISTa.gt.AAA) then
                DISTa = AAA
                MUR = j
                NUR = i
             endif
         ENDIF
         IF((dlon(j,i).lt.olon(m).and.olon(m)-dlon(j,i).lt.10.) .and. &
            (dlat(j,i).ge.olat(n).and.dlat(j,i)-olat(n).lt.10.)) then
             AAA = ((dlon(j,i)-olon(m)) &
                 *cos((dlat(j,i)+olat(n))/360.*pi))**2 &
                 +(dlat(j,i)-olat(n))**2
             if(DISTb.gt.AAA) then
                DISTb = AAA
                MUL = j
                NUL = i
             endif
         ENDIF
         IF((dlon(j,i).ge.olon(m).and.dlon(j,i)-olon(m).lt.10.) .and. &
            (dlat(j,i).lt.olat(n).and.olat(n)-dlat(j,i).lt.10.)) then
             AAA = ((dlon(j,i)-olon(m)) &
                 *cos((dlat(j,i)+olat(n))/360.*pi))**2 &
                 +(dlat(j,i)-olat(n))**2
             if(DISTc.gt.AAA) then
                DISTc = AAA
                MDR = j
                NDR = i
             endif
         ENDIF
         IF((dlon(j,i).lt.olon(m).and.olon(m)-dlon(j,i).lt.10.) .and. &
            (dlat(j,i).lt.olat(n).and.olat(n)-dlat(j,i).lt.10.)) then
             AAA = ((dlon(j,i)-olon(m)) &
                 *cos((dlat(j,i)+olat(n))/360.*pi))**2 &
                 +(dlat(j,i)-olat(n))**2
             if(DISTd.gt.AAA) then
                DISTd = AAA
                MDL = j
                NDL = i
             endif
         ENDIF
            enddo
            enddo
            
            DIST=amin1(DISTa,DISTb,DISTc,DISTd)
            I2UR(m,n) = MUR
            J2UR(m,n) = NUR
            I2UL(m,n) = MUL
            J2UL(m,n) = NUL
            I2DR(m,n) = MDR
            J2DR(m,n) = NDR
            I2DL(m,n) = MDL
            J2DL(m,n) = NDL
            D2XT(m,n) = DIST
            D2Xa(m,n) = DISTa
            D2Xb(m,n) = DISTb
            D2Xc(m,n) = DISTc
            D2Xd(m,n) = DISTd
         enddo
         enddo

         do nfile=1,n_month
            if(ntype.eq.0) then
               nyear = idate1/1000000
               month = idate1/10000-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(iproj.eq.'LAMCON') then
                  if(nfile.eq.1.or.month.eq.1) then
                     write(*,*)'regrid: LAMCON ICBC Orig.',nyear,month
                  else
                     write(*,*)'                         ',nyear,month
                  endif
               else if(iproj.eq.'ROTMER') then
                  if(nfile.eq.1.or.month.eq.1) then
                     write(*,*)'regrid: ROTMER ICBC Orig.',nyear,month
                  else
                     write(*,*)'                         ',nyear,month
                  endif
               endif

               if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
                  month.eq.7.or.month.eq.8.or.month.eq.10.or. &
                  month.eq.12) then
                  n_slice = 31*24/6+1
               else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                       month.eq.11) then
                  n_slice = 30*24/6+1
               else
                  n_slice = 28*24/6+1
                  if(mod(nyear,4).eq.0) n_slice = 29*24/6+1
                  if(mod(nyear,100).eq.0) n_slice = 28*24/6+1
                  if(mod(nyear,400).eq.0) n_slice = 29*24/6+1
               endif
               if(nfile.eq.n_month.and.mod(idate2/100-1,100).ne.0) &
               n_slice=min(n_slice,nint(mod(idate2/100-1,100)*24./6))+1
               filein = 'ICBC'//chy(nyear)//chm(month)//'0100'
               fileout= 'ICBC_LL.'//chy(nyear)//chm(month)//'0100'

               inquire(file= &
                    trim(Path_Input)//trim(DomainName)//'_'//filein &
                      ,exist=there)
               if(.not.there) then
                  write(*,*) &
                    trim(Path_Input)//trim(DomainName)//'_'//filein &
                        ,' is not avaiable'
                  stop
               endif
               open(10,file= &
                    trim(Path_Input)//trim(DomainName)//'_'//filein  &
                   ,form='unformatted',recl=iy*jx*ibyte,access='direct')
               mrec = 0
               open(20,file= &
                    trim(Path_Input)//trim(DomainName)//'_'//fileout &
                   ,form='unformatted' &
               ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
               nrec = 0
            else if(ntype.eq.1) then
               nyear = idate1/10000
               month = idate1/100-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(iproj.eq.'LAMCON') then
                  if(nfile.eq.1.or.month.eq.1) then
                     write(*,*)'regrid: LAMCON ICBC Daily',nyear,month
                  else
                     write(*,*)'                         ',nyear,month
                  endif
               else if(iproj.eq.'ROTMER') then
                  if(nfile.eq.1.or.month.eq.1) then
                     write(*,*)'regrid: ROTMER ICBC Daily',nyear,month
                  else
                     write(*,*)'                         ',nyear,month
                  endif
               endif

               if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
                  month.eq.7.or.month.eq.8.or.month.eq.10.or. &
                  month.eq.12) then
                  n_slice = 31
               else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                       month.eq.11) then
                  n_slice = 30
               else
                  n_slice = 28
                  if(mod(nyear,4).eq.0) n_slice = 29
                  if(mod(nyear,100).eq.0) n_slice = 28
                  if(mod(nyear,400).eq.0) n_slice = 29
               endif
               if(nfile.eq.n_month.and.mod(idate2-1,100).ne.0) &
               n_slice=min(n_slice,mod(idate2-1,100)+1)
               filein = 'ICBC'//chy(nyear)//chm(month)//'01'
               fileout= 'ICBC_LL.'//chy(nyear)//chm(month)//'01'

               inquire(file= &
                 trim(Path_Input)//trim(DomainName)//'_'//filein(1:12) &
                      ,exist=there)
               if(.not.there) then
                  write(*,*) &
                 trim(Path_Input)//trim(DomainName)//'_'//filein(1:12) &
                        ,' is not avaiable'
                  stop
               endif
               open(10,file= &
                trim(Path_Input)//trim(DomainName)//'_'//filein(1:12) &
                   ,form='unformatted' &
                   ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
               mrec = 0
               open(20,file= &
                trim(Path_Input)//trim(DomainName)//'_'//fileout(1:16) &
                   ,form='unformatted' &
               ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
               nrec = 0
            else if(ntype.eq.2) then
               nyear = idate1/100
               month = idate1-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(iproj.eq.'LAMCON') then
                  if(month.eq.1.or.nfile.eq.1) then
                     write(*,*)'regrid: LAMCON ICBC Month',nyear,month
                  else
                     write(*,*)'                         ',nyear,month
                  endif
               else if(iproj.eq.'ROTMER') then
                  if(month.eq.1.or.nfile.eq.1) then
                     write(*,*)'regrid: ROTMER ICBC Month',nyear,month
                  else
                     write(*,*)'                         ',nyear,month
                  endif
               endif

               n_slice = 1

               if(nfile.eq.1) then
                  filein = 'ICBC'
                  fileout= 'ICBC_LL'

                  inquire(file= &
          trim(Path_Input)//trim(DomainName)//'_'//filein(1:4)//'.mon' &
                         ,exist=there)
                  if(.not.there) then
                  write(*,*) &
          trim(Path_Input)//trim(DomainName)//'_'//filein(1:4)//'.mon' &
                        ,' is not avaiable'
                     stop
                  endif
                  open(10,file= &
          trim(Path_Input)//trim(DomainName)//'_'//filein(1:4)//'.mon' &
                      ,form='unformatted' &
                      ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
                  mrec = 0
                open(20,file= &
         trim(Path_Input)//trim(DomainName)//'_'//fileout(1:7)//'.mon' &
                    ,form='unformatted' &
                ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
                  nrec = 0
               endif
            endif

            do nnn=1,n_slice
               if(ntype.eq.0) then
                  mrec = mrec+1
                  read(10,rec=mrec) IDATE
!                 write(*,*) 'IDATE = ',IDATE
                  do l=1,kz*4+2
                     mrec = mrec+1
                     read(10,rec=mrec) ((o(j,i,l),j=1,jx),i=1,iy)
                  enddo
               else
                  do l=1,kz*4+2
                     mrec = mrec+1
                     read(10,rec=mrec) ((o(j,i,l),j=2,jx-1),i=2,iy-1)
                  enddo
               endif
               IF(iproj.eq.'ROTMER') THEN
                  IF(PLAT.GT.0.) THEN
                     POLLAM = PLON + 180.
                     POLPHI = 90. - PLAT
                  ELSE
                     POLLAM = PLON
                     POLPHI = 90. + PLAT
                  ENDIF
                  IF(POLLAM.GT.180.) POLLAM = POLLAM - 360.
                  POLCPHI = cos(PIR180*POLPHI)
                  POLSPHI = sin(PIR180*POLPHI)
                  do i=2,iy-1
                  do j=2,jx-1
                     ZPHI = DLAT(j,i)*PIR180
                     ZRLA = DLON(j,i)*PIR180
                     IF(DLAT(j,i).gt.89.999999) ZRLA = 0.0
                     ZRLAP = POLLAM*PIR180 - ZRLA
                     ZARG1  = POLCPHI*sin(ZRLAP)
                     ZARG2  = POLSPHI*cos(ZPHI)  &
                            - POLCPHI*sin(ZPHI)*cos(ZRLAP)
                     ZNORM  = 1.0/sqrt(ZARG1**2+ZARG2**2)
                     SINDEL = ZARG1*ZNORM
                     COSDEL = ZARG2*ZNORM
                     do k=1,kz
                        US =  o(j,i,k)*COSDEL + o(j,i,k+kz)*SINDEL
                        VS = -o(j,i,k)*SINDEL + o(j,i,k+kz)*COSDEL
                        o(j,i,k) = US
                        o(j,i,k+kz) = VS
                     enddo
                  enddo
                  enddo
               ELSE IF(iproj.eq.'LAMCON') THEN
                  IF(CLAT.lt.0.) THEN
                     SIGN0= -1.           ! SOUTH HEMESPHERE
                  ELSE
                     SIGN0=  1.           ! NORTH HEMESPHERE
                  ENDIF
                  IF(abs(truelatL-truelatH).gt.0.1) THEN
                     GRIDFC=(alog10(cos(truelatL*PIR180)) &
                            -alog10(cos(truelatH*PIR180)))  &
               /(alog10(tan(45.0-SIGN0*truelatL/2.0*PIR180))  &
                -alog10(tan(45.0-SIGN0*truelatH/2.0*PIR180)))
                  ELSE
                     GRIDFC=SIGN0*sin(truelatL*PIR180)
                  ENDIF
                  do i=2,iy-1
                  do j=2,jx-1
                     IF((CLON.ge.0.0.and.DLON(j,i).ge.0.).or. &
                        (CLON.lt.0.0.and.DLON(j,i).lt.0.)) THEN
                        X=(CLON-DLON(j,i))*PIR180*GRIDFC
                     ELSE
                        IF(CLON.ge.0.0) THEN
                           IF(abs(CLON-(DLON(j,i)+360.)).lt.   &
                              abs(CLON- DLON(j,i)) ) THEN
                              X=(CLON-(DLON(j,i)+360.))*PIR180*GRIDFC
                           ELSE
                              X=(CLON-DLON(j,i))*PIR180*GRIDFC
                           ENDIF
                        ELSE
                           IF(abs(CLON-(DLON(j,i)-360.)).LT. &
                              abs(CLON- DLON(j,i)) ) THEN
                              X=(CLON-(DLON(j,i)-360.))*PIR180*GRIDFC
                           ELSE
                              X=(CLON-DLON(j,i))*PIR180*GRIDFC
                           ENDIF
                        ENDIF
                     ENDIF
                     SINDEL=sin(X)
                     COSDEL=cos(X)
                     IF(CLAT.ge.0.) THEN
                        do k=1,kz
                           US= o(j,i,k)*COSDEL - o(j,i,k+kz)*SINDEL
                           VS=-o(j,i,k)*SINDEL + o(j,i,k+kz)*COSDEL
                           o(j,i,k) = US
                           o(j,i,k+kz) = VS
                        enddo
                     ELSE
                        do k=1,kz
                           US= o(j,i,k)*COSDEL + o(j,i,k+kz)*SINDEL
                           VS=-o(j,i,k)*SINDEL + o(j,i,k+kz)*COSDEL
                           o(j,i,k) = US
                           o(j,i,k+kz) = VS
                        enddo
                     ENDIF
                  enddo
                  enddo
               ENDIF
               do l=1,kz*2
                  do n=nlat,mlat
                  do m=nlon,mlon
                     out(m,n) = -9999.

            if(I2UR(m,n).lt.999.and.J2UR(m,n).lt.999.and. &
               I2UL(m,n).lt.999.and.J2UL(m,n).lt.999.and. &
               I2DR(m,n).lt.999.and.J2DR(m,n).lt.999.and. &
               I2DL(m,n).lt.999.and.J2DL(m,n).lt.999) then
               if(D2XT(m,n).gt.0.0001) then
      out(m,n) = ( o(I2UR(m,n),J2UR(m,n),l)/D2Xa(m,n) &
                  +o(I2UL(m,n),J2UL(m,n),l)/D2Xb(m,n) &
                  +o(I2DR(m,n),J2DR(m,n),l)/D2Xc(m,n) &
                  +o(I2DL(m,n),J2DL(m,n),l)/D2Xd(m,n) )  &
          /( 1./D2Xa(m,n)+1./D2Xb(m,n)+1./D2Xc(m,n)+1./D2Xd(m,n) )
               else
                  if(D2Xa(m,n).eq.D2XT(m,n)) then
                     out(m,n) = o(I2UR(m,n),J2UR(m,n),l)
                  else if(D2Xb(m,n).eq.D2XT(m,n)) then
                     out(m,n) = o(I2UL(m,n),J2UL(m,n),l)
                  else if(D2Xc(m,n).eq.D2XT(m,n)) then
                     out(m,n) = o(I2DR(m,n),J2DR(m,n),l)
                  else if(D2Xd(m,n).eq.D2XT(m,n)) then
                     out(m,n) = o(I2DL(m,n),J2DL(m,n),l)
                  endif
               endif
            endif
                  enddo
                  enddo
                  nrec=nrec+1
                  write(20,rec=nrec)((out(m,n),m=nlon,mlon),n=nlat,mlat)
               enddo
               do l=kz*2+1,kz*4+2
                  do n=nlat,mlat
                  do m=nlon,mlon
                     out(m,n) = -9999.

            if(I1UR(m,n).lt.999.and.J1UR(m,n).lt.999.and. &
               I1UL(m,n).lt.999.and.J1UL(m,n).lt.999.and. &
               I1DR(m,n).lt.999.and.J1DR(m,n).lt.999.and. &
               I1DL(m,n).lt.999.and.J1DL(m,n).lt.999) then
               if(D1XT(m,n).gt.0.0001) then
      out(m,n) = ( o(I1UR(m,n),J1UR(m,n),l)/D1Xa(m,n) &
                  +o(I1UL(m,n),J1UL(m,n),l)/D1Xb(m,n) &
                  +o(I1DR(m,n),J1DR(m,n),l)/D1Xc(m,n) &
                  +o(I1DL(m,n),J1DL(m,n),l)/D1Xd(m,n) )  &
          /( 1./D1Xa(m,n)+1./D1Xb(m,n)+1./D1Xc(m,n)+1./D1Xd(m,n) )
               else
                  if(D1Xa(m,n).eq.D1XT(m,n)) then
                     out(m,n) = o(I1UR(m,n),J1UR(m,n),l)
                  else if(D1Xb(m,n).eq.D1XT(m,n)) then
                     out(m,n) = o(I1UL(m,n),J1UL(m,n),l)
                  else if(D1Xc(m,n).eq.D1XT(m,n)) then
                     out(m,n) = o(I1DR(m,n),J1DR(m,n),l)
                  else if(D1Xd(m,n).eq.D1XT(m,n)) then
                     out(m,n) = o(I1DL(m,n),J1DL(m,n),l)
                  endif
               endif
            endif
                  enddo
                  enddo
                  nrec=nrec+1
                  write(20,rec=nrec)((out(m,n),m=nlon,mlon),n=nlat,mlat)
               enddo
            enddo
            if(.not.ntype.eq.2) close(20)
            if(.not.ntype.eq.2) close(10)
            if(igrads.eq.1.and.(ntype.eq.0.or.ntype.eq.1.or. &
               (ntype.eq.2.and.nfile.eq.1))) then
               if(ntype.eq.0) then
                  inquire(file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:18)//'.ctl' &
                         ,exist=there)
                  if(there) then
                     open(31,file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:18)//'.ctl' &
                      ,form='formatted',status='replace')
                  else
                     open(31,file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:18)//'.ctl' &
                      ,form='formatted',status='new')
                  endif
                  write(31,10) '^'//trim(DomainName)//'_'//fileout(1:18)
               else if(ntype.eq.1) then
                  inquire(file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:16)//'.ctl' &
                         ,exist=there)
                  if(there) then
                     open(31,file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:16)//'.ctl' &
                      ,form='formatted',status='replace')
                  else
                     open(31,file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:16)//'.ctl' &
                      ,form='formatted',status='new')
                  endif
                  write(31,11) '^'//trim(DomainName)//'_'//fileout(1:16)
               else if(ntype.eq.2) then
                  inquire(file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:7)//'.ctl' &
                         ,exist=there)
                  if(there) then
                     open(31,file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:7)//'.ctl' &
                      ,form='formatted',status='replace')
                  else
                     open(31,file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:7)//'.ctl' &
                      ,form='formatted',status='new')
                  endif
          write(31,12) '^'//trim(DomainName)//'_'//fileout(1:7)//'.mon'
               endif
               write(31,20)
               if(ibigend.eq.1) then
                  write(31,30)
               else
                  write(31,40)
               endif
               write(31,50)
               write(31,200) mlon-nlon+1,(nlon-1)*ddeg,ddeg
               write(31,210) mlat-nlat+1,(nlat-1)*ddeg,ddeg
               write(31,300) kz, &
       ((1013.25-ptop*10.)*(sigma(k)+sigma(k+1))*0.5+ptop*10.,k=kz,1,-1)
               if(ntype.eq.0) then
                  write(31,400) n_slice,0,'01',chmc(month),nyear,6
               else if(ntype.eq.1) then
                  write(31,401) n_slice,'01',chmc(month),nyear
               else if(ntype.eq.2) then
                  write(31,402) n_month,'16',chmc(month),nyear
               endif
               write(31,500) 6
               write(31,650) 'u       ',kz,'westerly wind (m/s)        '
               write(31,650) 'v       ',kz,'southerly wind (m/s)       '
               write(31,650) 't       ',kz,'air temperature (degree, K)'
               write(31,650) 'q       ',kz,'water vapor mixing ratio   '
               write(31,600) 'ps      ',   'surface pressure (hPa)     '
               write(31,600) 'tas     ',   'surface air temperature, K '
               write(31,700)
               close(31)
            endif
         enddo
            
         deallocate(I1UR)
         deallocate(J1UR)
         deallocate(I1UL)
         deallocate(J1UL)
         deallocate(I1DR)
         deallocate(J1DR)
         deallocate(I1DL)
         deallocate(J1DL)
         deallocate(D1XT)
         deallocate(D1Xa)
         deallocate(D1Xb)
         deallocate(D1Xc)
         deallocate(D1Xd)
            
         deallocate(I2UR)
         deallocate(J2UR)
         deallocate(I2UL)
         deallocate(J2UL)
         deallocate(I2DR)
         deallocate(J2DR)
         deallocate(I2DL)
         deallocate(J2DL)
         deallocate(D2XT)
         deallocate(D2Xa)
         deallocate(D2Xb)
         deallocate(D2Xc)
         deallocate(D2Xd)
      endif
      deallocate(out)
      deallocate(olat)
      deallocate(olon)

      deallocate(sigma)
      deallocate(o)
      deallocate(xlat)
      deallocate(xlon)
      deallocate(dlat)
      deallocate(dlon)

  10  format('dset ',A38)
  11  format('dset ',A36)
  12  format('dset ',A32)
  20  format('title RegCM domain information')
  30  format('options big_endian')
  40  format('options little_endian')
  50  format('undef -9999.')
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
 210  format('ydef ',I3,' linear ',f9.4,' ',f9.4)
 300  format('zdef ',I2,' levels ',30f7.2)
 400  format('tdef ',I4,' linear ',I2,'z',A2,A3,I4,' ',I2,'hr')
 401  format('tdef ',I4,' linear ',A2,A3,I4,' ','1dy')
 402  format('tdef ',I4,' linear ',A2,A3,I4,' ','1mo')
 500  format('vars ',I2)
 600  format(A8,'0 99 ',A26)           
 650  format(A8,I2,' 0 ',A26)
 700  format('endvars')
      return
      end

      subroutine regrid_ATM(iy,jx,kz,ibyte,idate0,idate1,idate2 &
                           ,Path_Output)
      implicit none
      integer iy,jx,kz,ibyte,idate0,idate1,idate2
      character*128 Path_Output
      real*4  truelatL,truelatH
      integer iiy,jjx,kkz
      integer mdate0,ibltyp,icup,ipptls,iboudy
      real*4  dxsp,ptsp,clat,clon,plat,plon
      real*4  dto,dtb,dtr,dtc
      integer iotyp
      character*6 iproj
      real*4, allocatable ::  sigma(:)
      character*4 :: chy(1941:2100)
      data chy/'1941','1942','1943','1944','1945','1946','1947','1948',&
        '1949','1950','1951','1952','1953','1954','1955','1956','1957',&
        '1958','1959','1960','1961','1962','1963','1964','1965','1966',&
        '1967','1968','1969','1970','1971','1972','1973','1974','1975',&
        '1976','1977','1978','1979','1980','1981','1982','1983','1984',&
        '1985','1986','1987','1988','1989','1990','1991','1992','1993',&
        '1994','1995','1996','1997','1998','1999','2000','2001','2002',&
        '2003','2004','2005','2006','2007','2008','2009','2010','2011',&
        '2012','2013','2014','2015','2016','2017','2018','2019','2020',&
        '2021','2022','2023','2024','2025','2026','2027','2028','2029',&
        '2030','2031','2032','2033','2034','2035','2036','2037','2038',&
        '2039','2040','2041','2042','2043','2044','2045','2046','2047',&
        '2048','2049','2050','2051','2052','2053','2054','2055','2056',&
        '2057','2058','2059','2060','2061','2062','2063','2064','2065',&
        '2066','2067','2068','2069','2070','2071','2072','2073','2074',&
        '2075','2076','2077','2078','2079','2080','2081','2082','2083',&
        '2084','2085','2086','2087','2088','2089','2090','2091','2092',&
        '2093','2094','2095','2096','2097','2098','2099','2100'/
      character*2 chm(12)
      data chm/'01','02','03','04','05','06','07','08','09','10', &
               '11','12'/
      character*3 chmc(12)
      data chmc/'jan','feb','mar','apr','may','jun'  &
               ,'jul','aug','sep','oct','nov','dec'/
      character*14 filein
      character*17 fileout
      integer ntype,nfile,nyear,month,n_slice,mrec,nrec,nnn

      integer i,j,k,l,m,n
      integer igrads,ibigend

      real*4  size_2d(4),ddeg
      COMMON /WINDOW/ size_2d,ddeg

      real*4, allocatable ::  o(:,:,:),xlat(:,:),xlon(:,:)
      real*4  xmaxlat,xminlat,xmaxlon,xminlon

      integer mlat,nlat,mlon,nlon
      integer n_month
      logical there

      integer MUR,NUR,MUL,NUL,MDR,NDR,MDL,NDL
      real*4  DISTa,DISTb,DISTc,DISTd,AAA,DIST
      real*4  PI,PIR180
      real*4  POLLAM,POLPHI,POLCPHI,POLSPHI,ZPHI,ZRLA,ZRLAP
      real*4  ZARG1,ZARG2,ZNORM,SINDEL,COSDEL,US,VS,X
      real*4  SIGN0,GRIDFC
      real*4, allocatable :: out(:,:)
      real*4, allocatable :: olat(:),olon(:)

      integer, allocatable :: I1UR(:,:), J1UR(:,:)
      integer, allocatable :: I1UL(:,:), J1UL(:,:)
      integer, allocatable :: I1DR(:,:), J1DR(:,:)
      integer, allocatable :: I1DL(:,:), J1DL(:,:)
      real*4, allocatable :: D1XT(:,:), D1Xa(:,:), D1Xb(:,:)
      real*4, allocatable ::            D1Xc(:,:), D1Xd(:,:)

      igrads = 1
      ibigend= 1

      PI = atan(1.)*4.
      PIR180 = atan(1.)/45.

      allocate(sigma(kz+1))
      allocate(o(jx-2,iy-2,kz*6+5))
      allocate(xlat(jx-2,iy-2))
      allocate(xlon(jx-2,iy-2))

      inquire(file=trim(Path_Output)//'OUT_HEAD',exist=there)
      if(.not.there) then
         write(*,*) trim(Path_Output)//'OUT_HEAD',' is not avaiable'
         stop
      endif
      open(10,file=trim(Path_Output)//'OUT_HEAD',form='unformatted' &
             ,recl=(jx-2)*(iy-2)*ibyte,access='direct')
      read(10,rec=1) mdate0,ibltyp,icup,ipptls,iboudy  &
                    ,iiy,jjx,kkz,(sigma(k),k=1,kz+1)   &
                    ,dxsp,ptsp,clat,clon,plat,plon          &
                    ,iproj,dto,dtb,dtr,dtc,iotyp,truelatL,truelatH
      if(iiy.ne.iy.or.jjx.ne.jx.or.kkz.ne.kz) then
         write(*,*) 'iy,jx,kz in parameter = ',iy,jx,kz
         write(*,*) 'iy,jx,kz in OUT_HEAD ',iiy,jjx,kkz
         write(*,*) 'They are not consistent'
         stop
      endif
      read(10,rec=6) xlat         ! xlat
      if(size_2d(3).gt.-9990.) then
         xminlat = size_2d(3)
         xmaxlat = size_2d(4)
      else
         xmaxlat = -1.e10
         xminlat =  1.e10
         do i=1,iy-2
         do j=1,jx-2
            if(xlat(j,i).gt.xmaxlat) xmaxlat=xlat(j,i)
            if(xlat(j,i).lt.xminlat) xminlat=xlat(j,i)
         enddo
         enddo
      endif
      read(10,rec=7) xlon         ! xlon
      if(size_2d(1).gt.-9990.) then
         xminlon=size_2d(1)
         xmaxlon=size_2d(2)
      else
         xmaxlon = -1.e10
         xminlon =  1.e10
         do i=1,iy-2
         do j=1,jx-2
            if(xlon(j,i).gt.xmaxlon) xmaxlon=xlon(j,i)
            if(xlon(j,i).lt.xminlon) xminlon=xlon(j,i)
         enddo
         enddo
      endif
      close(10)
      
      mlat = xmaxlat/ddeg - 1
      nlat = xminlat/ddeg + 1
      mlon = xmaxlon/ddeg - 1
      nlon = xminlon/ddeg + 2

      write(*,*) xminlat,mlat-nlat+1,(nlat-1)*ddeg
      write(*,*) xminlon,mlon-nlon+1,(nlon-1)*ddeg

      allocate(out(nlon:mlon,nlat:mlat))
      allocate(olat(nlat:mlat))
      allocate(olon(nlon:mlon))
      do n=nlat,mlat
         olat(n) = real(n-1)*ddeg
      enddo
      do m=nlon,mlon
         olon(m) = real(m-1)*ddeg
      enddo

      if(idate1.gt.1940010100) then        ! original file
         n_month = (idate2/1000000-idate1/1000000)*12 &
                 + (mod(idate2/10000,100)-mod(idate1/10000,100))   
         if(mod(idate2,10000).gt.0100) n_month = n_month+1
         ntype = 0
      else if(idate1.gt.19400101) then     ! daily mean file
         n_month = (idate2/10000-idate1/10000)*12 &
                 + (mod(idate2/100,100)-mod(idate1/100,100))   
         if(mod(idate2,100).gt.01) n_month = n_month+1
         ntype = 1
      else if(idate1.gt.194001) then       ! monthly mean file
         n_month = (idate2/100-idate1/100)*12 &
                 + (mod(idate2,100)-mod(idate1,100))+1
         ntype = 2
      endif

      if(iproj.eq.'NORMER') then
         do nfile=1,n_month
            if(ntype.eq.0) then
               nyear = idate1/1000000
               month = idate1/10000-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(nfile.eq.1.or.month.eq.1) then
                  write(*,*)'regrid: NORMER ATM Orig.',nyear,month
               else
                  write(*,*)'                        ',nyear,month
               endif

               if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
                  month.eq.7.or.month.eq.8.or.month.eq.10.or. &
                  month.eq.12) then
                  n_slice = 31*24/dto
               else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                       month.eq.11) then
                  n_slice = 30*24/dto
               else
                  n_slice = 28*24/dto
                  if(mod(nyear,4).eq.0) n_slice = 29*24/dto
                  if(mod(nyear,100).eq.0) n_slice = 28*24/dto
                  if(mod(nyear,400).eq.0) n_slice = 29*24/dto
               endif
               if(nfile.eq.n_month.and.mod(idate2/100-1,100).ne.0) &
               n_slice=min(n_slice,nint(mod(idate2/100-1,100)*24./dto))
               if(nfile.eq.1.and.idate0.eq.idate1) n_slice=n_slice+1
               filein = 'ATM.'//chy(nyear)//chm(month)//'0100'
               fileout= 'ATM_LL.'//chy(nyear)//chm(month)//'0100'
               inquire(file=trim(Path_Output)//filein,exist=there)
               if(.not.there) then
                  write(*,*) trim(Path_Output)//filein,' is not avaiable'
                  stop
               endif
             open(10,file=trim(Path_Output)//filein,form='unformatted' &
                      ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
               mrec = 0
            open(20,file=trim(Path_Output)//fileout,form='unformatted' &
               ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
               nrec = 0
            else if(ntype.eq.1) then
               nyear = idate1/10000
               month = idate1/100-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(nfile.eq.1.or.month.eq.1) then
                  write(*,*)'regrid: NORMER ATM Daily',nyear,month
               else
                  write(*,*)'                        ',nyear,month
               endif
               if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
                  month.eq.7.or.month.eq.8.or.month.eq.10.or. &
                  month.eq.12) then
                  n_slice = 31
               else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                       month.eq.11) then
                  n_slice = 30
               else
                  n_slice = 28
                  if(mod(nyear,4).eq.0) n_slice = 29
                  if(mod(nyear,100).eq.0) n_slice = 28
                  if(mod(nyear,400).eq.0) n_slice = 29
               endif
               if(nfile.eq.n_month.and.mod(idate2-1,100).ne.0) &
               n_slice=min(n_slice,mod(idate2-1,100)+1)
               filein = 'ATM.'//chy(nyear)//chm(month)//'01'
               fileout= 'ATM_LL.'//chy(nyear)//chm(month)//'01'
               inquire(file=trim(Path_Output)//filein,exist=there)
               if(.not.there) then
                  write(*,*) trim(Path_Output)//filein,' is not avaiable'
                  stop
               endif
             open(10,file=trim(Path_Output)//filein,form='unformatted' &
                      ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
               mrec = 0
            open(20,file=trim(Path_Output)//fileout,form='unformatted' &
               ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
               nrec = 0

            else if(ntype.eq.2) then
               nyear = idate1/100
               month = idate1-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(month.eq.1.or.nfile.eq.1) then
                  write(*,*)'regrid: NORMER ATM Month',nyear,month
               else
                  write(*,*)'                        ',nyear,month
               endif

               n_slice = 1

               if(nfile.eq.1) then
                  filein = 'ATM.mon'
                  fileout= 'ATM_LL.mon'
                  inquire(file=trim(Path_Output)//filein(1:7),exist=there)
                  if(.not.there) then
                     write(*,*) trim(Path_Output)//filein(1:7),' is not avaiable'
                     stop
                  endif
             open(10,file=trim(Path_Output)//filein(1:7),form='unformatted' &
                      ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
                  nrec = 0
            open(20,file=trim(Path_Output)//fileout(1:10),form='unformatted' &
                ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
                  mrec = 0
               endif
            endif

            do nnn=1,n_slice
               do l=1,kz*6+5
                  mrec = mrec+1
                  read(10,rec=mrec) ((o(j,i,l),j=1,jx-2),i=1,iy-2)
               enddo
               do l=1,kz*6+5
                  do n=nlat,mlat
                  do m=nlon,mlon
                     out(m,n) = -9999.

                     do i=1,iy-3
                     do j=1,jx-3
                        if((xlat(j,i).lt.olat(n).and.      &
                            xlon(j,i).lt.olon(m))          &
                      .and.(xlat(j+1,i+1).ge.olat(n).and.  &
                            xlon(j+1,i+1).ge.olon(m))) then
                            out(m,n) = &
       ( o(j,i,l)*(xlon(j+1,i+1)-olon(m))*(xlat(j+1,i+1)-olat(n)) &
       + o(j+1,i,l)*(olon(m)-xlon(j,i+1))*(xlat(j+1,i+1)-olat(n)) &
       + o(j,i+1,l)*(xlon(j+1,i)-olon(m))*(olat(n)-xlat(j+1,i)) &
       + o(j+1,i+1,l)*(olon(m)-xlon(j,i))*(olat(n)-xlat(j,i)) ) &
       / ( (xlon(j+1,i+1)-xlon(j,i))*(xlat(j+1,i+1)-xlat(j,i)) )
                        endif
                     enddo
                     enddo
                  enddo
                  enddo
                  nrec = nrec+1
                  write(20,rec=nrec)((out(m,n),m=nlon,mlon),n=nlat,mlat)
               enddo
            enddo
            if(.not.ntype.eq.2) close(20)
            if(.not.ntype.eq.2) close(10)
            if(igrads.eq.1.and.(ntype.eq.0.or.ntype.eq.1.or. &
               (ntype.eq.2.and.nfile.eq.1))) then
               if(ntype.eq.0) then
                 inquire(file=trim(Path_Output)//fileout(1:17)//'.ctl' &
                        ,exist=there)
                  if(there) then
                     open(31,file= &
                          trim(Path_Output)//fileout(1:17)//'.ctl' &
                         ,form='formatted',status='replace')
                  else
                     open(31,file= &
                          trim(Path_Output)//fileout(1:17)//'.ctl' &
                         ,form='formatted',status='new')
                  endif
                  write(31,10) fileout(1:17)
               else if(ntype.eq.1) then
                 inquire(file=trim(Path_Output)//fileout(1:15)//'.ctl' &
                        ,exist=there)
                  if(there) then
                     open(31,file= &
                          trim(Path_Output)//fileout(1:15)//'.ctl' &
                         ,form='formatted',status='replace')
                  else
                     open(31,file= &
                          trim(Path_Output)//fileout(1:15)//'.ctl' &
                         ,form='formatted',status='new')
                  endif
                  write(31,11) fileout(1:15)
               else if(ntype.eq.2) then
                 inquire(file=trim(Path_Output)//fileout(1:10)//'.ctl' &
                        ,exist=there)
                  if(there) then
                     open(31,file= &
                          trim(Path_Output)//fileout(1:10)//'.ctl' &
                         ,form='formatted',status='replace')
                  else
                     open(31,file= &
                          trim(Path_Output)//fileout(1:10)//'.ctl' &
                         ,form='formatted',status='new')
                  endif
                  write(31,12) fileout(1:10)
               endif
               write(31,20)
               if(ibigend.eq.1) then
                  write(31,30)
               else
                  write(31,40)
               endif
               write(31,50)
               write(31,200) mlon-nlon+1,(nlon-1)*ddeg,ddeg
               write(31,210) mlat-nlat+1,(nlat-1)*ddeg,ddeg
               write(31,300) kz, &
          ((1013.25-ptsp*10.)*(sigma(k)+sigma(k+1))*0.5+ptsp*10.,k=1,kz)
               if(ntype.eq.0) then
               if(nfile.eq.1.and.idate0.eq.idate1) then
               write(31,400) n_slice,0,'01',chmc(month),nyear,nint(dto)
               else
       write(31,400) n_slice,nint(dto),'01',chmc(month),nyear,nint(dto)
               endif
               else if(ntype.eq.1) then
                  write(31,401) n_slice,'01',chmc(month),nyear
               else if(ntype.eq.2) then
                  write(31,402) n_month,'16',chmc(month),nyear
               endif
               write(31,500) 11
               write(31,650) 'u       ',kz,'westerly wind (m/s)        '
               write(31,650) 'v       ',kz,'southerly wind (m/s)       '
               write(31,650) 'w       ',kz,'omega (hPa/s)   p-velocity '
               write(31,650) 't       ',kz,'air temperature (degree)   '
               write(31,650) 'qv      ',kz,'water vapor mixing ratio   '
               write(31,650) 'qc      ',kz,'cloud water mixing ratio   '
               write(31,600) 'ps      ',   'surface pressure (hPa)     '
               write(31,600) 'tpr     ',   'total precipitation(mm/day)'
               write(31,600) 'tgb     ',   'lower groud temp. in BATS  '
               write(31,600) 'swt     ',   'total soil water in mm H2O '
               write(31,600) 'rno     ',   'accumulated infiltration   '
               write(31,700)
               close(31)
            endif
         enddo

      else if(iproj.eq.'LAMCON'.or.iproj.eq.'ROTMER') then

         allocate(I1UR(nlon:mlon,nlat:mlat))
         allocate(J1UR(nlon:mlon,nlat:mlat))
         allocate(I1UL(nlon:mlon,nlat:mlat))
         allocate(J1UL(nlon:mlon,nlat:mlat))
         allocate(I1DR(nlon:mlon,nlat:mlat))
         allocate(J1DR(nlon:mlon,nlat:mlat))
         allocate(I1DL(nlon:mlon,nlat:mlat))
         allocate(J1DL(nlon:mlon,nlat:mlat))
         allocate(D1XT(nlon:mlon,nlat:mlat))
         allocate(D1Xa(nlon:mlon,nlat:mlat))
         allocate(D1Xb(nlon:mlon,nlat:mlat))
         allocate(D1Xc(nlon:mlon,nlat:mlat))
         allocate(D1Xd(nlon:mlon,nlat:mlat))

         do n=nlat,mlat
         do m=nlon,mlon

            MUR=1000
            NUR=1000
            MUL=1000
            NUL=1000
            MDR=1000
            NDR=1000
            MDL=1000
            NDL=1000
            DISTa=1.E8
            DISTb=1.E8
            DISTc=1.E8
            DISTd=1.E8

            do i=1,iy-2
            do j=1,jx-2
         IF((xlon(j,i).ge.olon(m).and.xlon(j,i)-olon(m).lt.10.) .and. &
            (xlat(j,i).ge.olat(n).and.xlat(j,i)-olat(n).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTa.gt.AAA) then
                DISTa = AAA
                MUR = j
                NUR = i
             endif
         ENDIF
         IF((xlon(j,i).lt.olon(m).and.olon(m)-xlon(j,i).lt.10.) .and. &
            (xlat(j,i).ge.olat(n).and.xlat(j,i)-olat(n).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTb.gt.AAA) then
                DISTb = AAA
                MUL = j
                NUL = i
             endif
         ENDIF
         IF((xlon(j,i).ge.olon(m).and.xlon(j,i)-olon(m).lt.10.) .and. &
            (xlat(j,i).lt.olat(n).and.olat(n)-xlat(j,i).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTc.gt.AAA) then
                DISTc = AAA
                MDR = j
                NDR = i
             endif
         ENDIF
         IF((xlon(j,i).lt.olon(m).and.olon(m)-xlon(j,i).lt.10.) .and. &
            (xlat(j,i).lt.olat(n).and.olat(n)-xlat(j,i).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTd.gt.AAA) then
                DISTd = AAA
                MDL = j
                NDL = i
             endif
         ENDIF
            enddo
            enddo
            
            DIST=amin1(DISTa,DISTb,DISTc,DISTd)
            I1UR(m,n) = MUR
            J1UR(m,n) = NUR
            I1UL(m,n) = MUL
            J1UL(m,n) = NUL
            I1DR(m,n) = MDR
            J1DR(m,n) = NDR
            I1DL(m,n) = MDL
            J1DL(m,n) = NDL
            D1XT(m,n) = DIST
            D1Xa(m,n) = DISTa
            D1Xb(m,n) = DISTb
            D1Xc(m,n) = DISTc
            D1Xd(m,n) = DISTd
         enddo
         enddo

         do nfile=1,n_month
            if(ntype.eq.0) then
               nyear = idate1/1000000
               month = idate1/10000-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(iproj.eq.'LAMCON') then
                  if(nfile.eq.1.or.month.eq.1) then
                     write(*,*)'regrid: LAMCON ATM Orig.',nyear,month
                  else
                     write(*,*)'                        ',nyear,month
                  endif
               else if(iproj.eq.'ROTMER') then
                  if(nfile.eq.1.or.month.eq.1) then
                     write(*,*)'regrid: ROTMER ATM Orig.',nyear,month
                  else
                     write(*,*)'                        ',nyear,month
                  endif
               endif

               if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
                  month.eq.7.or.month.eq.8.or.month.eq.10.or. &
                  month.eq.12) then
                  n_slice = 31*24/dto
               else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                       month.eq.11) then
                  n_slice = 30*24/dto
               else
                  n_slice = 28*24/dto
                  if(mod(nyear,4).eq.0) n_slice = 29*24/dto
                  if(mod(nyear,100).eq.0) n_slice = 28*24/dto
                  if(mod(nyear,400).eq.0) n_slice = 29*24/dto
               endif
               if(nfile.eq.n_month.and.mod(idate2/100-1,100).ne.0) &
               n_slice=min(n_slice,nint(mod(idate2/100-1,100)*24./dto))
               if(nfile.eq.1.and.idate0.eq.idate1) n_slice=n_slice+1
               filein = 'ATM.'//chy(nyear)//chm(month)//'0100'
               fileout= 'ATM_LL.'//chy(nyear)//chm(month)//'0100'
               inquire(file=trim(Path_Output)//filein,exist=there)
               if(.not.there) then
                  write(*,*) trim(Path_Output)//filein,' is not avaiable'
                  stop
               endif
             open(10,file=trim(Path_Output)//filein,form='unformatted' &
                      ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
               mrec = 0
            open(20,file=trim(Path_Output)//fileout,form='unformatted' &
               ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
               nrec = 0
            else if(ntype.eq.1) then
               nyear = idate1/10000
               month = idate1/100-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(iproj.eq.'LAMCON') then
                  if(nfile.eq.1.or.month.eq.1) then
                     write(*,*)'regrid: LAMCON ATM Daily',nyear,month
                  else
                     write(*,*)'                        ',nyear,month
                  endif
               else if(iproj.eq.'ROTMER') then
                  if(nfile.eq.1.or.month.eq.1) then
                     write(*,*)'regrid: ROTMER ATM Daily',nyear,month
                  else
                     write(*,*)'                        ',nyear,month
                  endif
               endif

               if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
                  month.eq.7.or.month.eq.8.or.month.eq.10.or. &
                  month.eq.12) then
                  n_slice = 31
               else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                       month.eq.11) then
                  n_slice = 30
               else
                  n_slice = 28
                  if(mod(nyear,4).eq.0) n_slice = 29
                  if(mod(nyear,100).eq.0) n_slice = 28
                  if(mod(nyear,400).eq.0) n_slice = 29
               endif
               if(nfile.eq.n_month.and.mod(idate2-1,100).ne.0) &
               n_slice=min(n_slice,mod(idate2-1,100)+1)
               filein = 'ATM.'//chy(nyear)//chm(month)//'01'
               fileout= 'ATM_LL.'//chy(nyear)//chm(month)//'01'
               inquire(file=trim(Path_Output)//filein,exist=there)
               if(.not.there) then
                  write(*,*) trim(Path_Output)//filein,' is not avaiable'
                  stop
               endif
             open(10,file=trim(Path_Output)//filein,form='unformatted' &
                      ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
               mrec = 0
            open(20,file=trim(Path_Output)//fileout,form='unformatted' &
               ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
               nrec = 0
            else if(ntype.eq.2) then
               nyear = idate1/100
               month = idate1-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(iproj.eq.'LAMCON') then
                  if(month.eq.1.or.nfile.eq.1) then
                     write(*,*)'regrid: LAMCON ATM Month',nyear,month
                  else
                     write(*,*)'                        ',nyear,month
                  endif
               else if(iproj.eq.'ROTMER') then
                  if(month.eq.1.or.nfile.eq.1) then
                     write(*,*)'regrid: ROTMER ATM Month',nyear,month
                  else
                     write(*,*)'                        ',nyear,month
                  endif
               endif

               n_slice = 1

               if(nfile.eq.1) then
                  filein = 'ATM.mon'
                  fileout= 'ATM_LL.mon'
                  inquire(file=trim(Path_Output)//filein(1:7)  &
                         ,exist=there)
                  if(.not.there) then
                     write(*,*) trim(Path_Output)//filein(1:7)  &
                              ,' is not avaiable'
                     stop
                  endif
                  open(10,file=trim(Path_Output)//filein(1:7)  &
                         ,form='unformatted' &
                         ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
                  nrec = 0
                  open(20,file=trim(Path_Output)//fileout(1:10)  &
                         ,form='unformatted' &
                ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
                  mrec = 0
               endif
            endif

            do nnn=1,n_slice
               do l=1,kz*6+5
                  mrec = mrec+1
                  read(10,rec=mrec) ((o(j,i,l),j=1,jx-2),i=1,iy-2)
               enddo
               
               IF(iproj.eq.'ROTMER') THEN
                  IF(PLAT.GT.0.) THEN
                     POLLAM = PLON + 180.
                     POLPHI = 90. - PLAT
                  ELSE
                     POLLAM = PLON
                     POLPHI = 90. + PLAT
                  ENDIF
                  IF(POLLAM.GT.180.) POLLAM = POLLAM - 360.
                  POLCPHI = cos(PIR180*POLPHI)
                  POLSPHI = sin(PIR180*POLPHI)
                  do i=1,iy-2
                  do j=1,jx-2
                     ZPHI = XLAT(j,i)*PIR180
                     ZRLA = XLON(j,i)*PIR180
                     IF(XLAT(j,i).gt.89.999999) ZRLA = 0.0
                     ZRLAP = POLLAM*PIR180 - ZRLA
                     ZARG1  = POLCPHI*sin(ZRLAP)
                     ZARG2  = POLSPHI*cos(ZPHI)  &
                            - POLCPHI*sin(ZPHI)*cos(ZRLAP)
                     ZNORM  = 1.0/sqrt(ZARG1**2+ZARG2**2)
                     SINDEL = ZARG1*ZNORM
                     COSDEL = ZARG2*ZNORM
                     do k=1,kz
                        US =  o(j,i,k)*COSDEL + o(j,i,k+kz)*SINDEL
                        VS = -o(j,i,k)*SINDEL + o(j,i,k+kz)*COSDEL
                        o(j,i,k) = US
                        o(j,i,k+kz) = VS
                     enddo
                  enddo
                  enddo
               ELSE IF(iproj.eq.'LAMCON') THEN
                  IF(CLAT.lt.0.) THEN
                     SIGN0= -1.           ! SOUTH HEMESPHERE
                  ELSE
                     SIGN0=  1.           ! NORTH HEMESPHERE
                  ENDIF
                  IF(abs(truelatL-truelatH).gt.0.1) THEN
                     GRIDFC=(alog10(cos(truelatL*PIR180)) &
                            -alog10(cos(truelatH*PIR180)))  &
               /(alog10(tan(45.0-SIGN0*truelatL/2.0*PIR180))  &
                -alog10(tan(45.0-SIGN0*truelatH/2.0*PIR180)))
                  ELSE
                     GRIDFC=SIGN0*sin(truelatL*PIR180)
                  ENDIF
                  do i=1,iy-2
                  do j=1,jx-2
                     IF((CLON.ge.0.0.and.XLON(j,i).ge.0.).or. &
                        (CLON.lt.0.0.and.XLON(j,i).lt.0.)) THEN
                        X=(CLON-XLON(j,i))*PIR180*GRIDFC
                     ELSE
                        IF(CLON.ge.0.0) THEN
                           IF(abs(CLON-(XLON(j,i)+360.)).lt.   &
                              abs(CLON- XLON(j,i)) ) THEN
                              X=(CLON-(XLON(j,i)+360.))*PIR180*GRIDFC
                           ELSE
                              X=(CLON-XLON(j,i))*PIR180*GRIDFC
                           ENDIF
                        ELSE
                           IF(abs(CLON-(XLON(j,i)-360.)).LT. &
                              abs(CLON- XLON(j,i)) ) THEN
                              X=(CLON-(XLON(j,i)-360.))*PIR180*GRIDFC
                           ELSE
                              X=(CLON-XLON(j,i))*PIR180*GRIDFC
                           ENDIF
                        ENDIF
                     ENDIF
                     SINDEL=sin(X)
                     COSDEL=cos(X)
                     IF(CLAT.ge.0.) THEN
                        do k=1,kz
                           US= o(j,i,k)*COSDEL - o(j,i,k+kz)*SINDEL
                           VS=-o(j,i,k)*SINDEL + o(j,i,k+kz)*COSDEL
                           o(j,i,k) = US
                           o(j,i,k+kz) = VS
                        enddo
                     ELSE
                        do k=1,kz
                           US= o(j,i,k)*COSDEL + o(j,i,k+kz)*SINDEL
                           VS=-o(j,i,k)*SINDEL + o(j,i,k+kz)*COSDEL
                           o(j,i,k) = US
                           o(j,i,k+kz) = VS
                        enddo
                     ENDIF
                  enddo
                  enddo
               ENDIF

               do l=1,kz*6+5
                  do n=nlat,mlat
                  do m=nlon,mlon
                     out(m,n) = -9999.

            if(I1UR(m,n).lt.999.and.J1UR(m,n).lt.999.and. &
               I1UL(m,n).lt.999.and.J1UL(m,n).lt.999.and. &
               I1DR(m,n).lt.999.and.J1DR(m,n).lt.999.and. &
               I1DL(m,n).lt.999.and.J1DL(m,n).lt.999) then
               if(D1XT(m,n).gt.0.0001) then
      out(m,n) = ( o(I1UR(m,n),J1UR(m,n),l)/D1Xa(m,n) &
                  +o(I1UL(m,n),J1UL(m,n),l)/D1Xb(m,n) &
                  +o(I1DR(m,n),J1DR(m,n),l)/D1Xc(m,n) &
                  +o(I1DL(m,n),J1DL(m,n),l)/D1Xd(m,n) )  &
          /( 1./D1Xa(m,n)+1./D1Xb(m,n)+1./D1Xc(m,n)+1./D1Xd(m,n) )
               else
                  if(D1Xa(m,n).eq.D1XT(m,n)) then
                     out(m,n) = o(I1UR(m,n),J1UR(m,n),l)
                  else if(D1Xb(m,n).eq.D1XT(m,n)) then
                     out(m,n) = o(I1UL(m,n),J1UL(m,n),l)
                  else if(D1Xc(m,n).eq.D1XT(m,n)) then
                     out(m,n) = o(I1DR(m,n),J1DR(m,n),l)
                  else if(D1Xd(m,n).eq.D1XT(m,n)) then
                     out(m,n) = o(I1DL(m,n),J1DL(m,n),l)
                  endif
               endif
            endif
                  enddo
                  enddo
                  nrec=nrec+1
                  write(20,rec=nrec)((out(m,n),m=nlon,mlon),n=nlat,mlat)
               enddo
            enddo
            if(.not.ntype.eq.2) close(20)
            if(.not.ntype.eq.2) close(10)
            if(igrads.eq.1.and.(ntype.eq.0.or.ntype.eq.1.or. &
               (ntype.eq.2.and.nfile.eq.1))) then
               if(ntype.eq.0) then
                 inquire(file=trim(Path_Output)//fileout(1:17)//'.ctl' &
                        ,exist=there)
                  if(there) then
                     open(31,file= &
                          trim(Path_Output)//fileout(1:17)//'.ctl' &
                         ,form='formatted',status='replace')
                  else
                     open(31,file= &
                          trim(Path_Output)//fileout(1:17)//'.ctl' &
                         ,form='formatted',status='new')
                  endif
                  write(31,10) fileout(1:17)
               else if(ntype.eq.1) then
                 inquire(file=trim(Path_Output)//fileout(1:15)//'.ctl' &
                        ,exist=there)
                  if(there) then
                     open(31,file= &
                          trim(Path_Output)//fileout(1:15)//'.ctl' &
                         ,form='formatted',status='replace')
                  else
                     open(31,file= &
                          trim(Path_Output)//fileout(1:15)//'.ctl' &
                         ,form='formatted',status='new')
                  endif
                  write(31,11) fileout(1:15)
               else if(ntype.eq.2) then
                 inquire(file=trim(Path_Output)//fileout(1:10)//'.ctl' &
                        ,exist=there)
                  if(there) then
                     open(31,file= &
                          trim(Path_Output)//fileout(1:10)//'.ctl' &
                         ,form='formatted',status='replace')
                  else
                     open(31,file= &
                          trim(Path_Output)//fileout(1:10)//'.ctl' &
                         ,form='formatted',status='new')
                  endif
                  write(31,12) fileout(1:10)
               endif
               write(31,20)
               if(ibigend.eq.1) then
                  write(31,30)
               else
                  write(31,40)
               endif
               write(31,50)
               write(31,200) mlon-nlon+1,(nlon-1)*ddeg,ddeg
               write(31,210) mlat-nlat+1,(nlat-1)*ddeg,ddeg
               write(31,300) kz, &
          ((1013.25-ptsp*10.)*(sigma(k)+sigma(k+1))*0.5+ptsp*10.,k=1,kz)
               if(ntype.eq.0) then
               if(nfile.eq.1.and.idate0.eq.idate1) then
               write(31,400) n_slice,0,'01',chmc(month),nyear,nint(dto)
               else
       write(31,400) n_slice,nint(dto),'01',chmc(month),nyear,nint(dto)
               endif
               else if(ntype.eq.1) then
                  write(31,401) n_slice,'01',chmc(month),nyear
               else if(ntype.eq.2) then
                  write(31,402) n_month,'16',chmc(month),nyear
               endif
               write(31,500) 11
               write(31,650) 'u       ',kz,'westerly wind (m/s)        '
               write(31,650) 'v       ',kz,'southerly wind (m/s)       '
               write(31,650) 'w       ',kz,'omega (hPa/s)   p-velocity '
               write(31,650) 't       ',kz,'air temperature (degree)   '
               write(31,650) 'qv      ',kz,'water vapor mixing ratio   '
               write(31,650) 'qc      ',kz,'cloud water mixing ratio   '
               write(31,600) 'ps      ',   'surface pressure (hPa)     '
               write(31,600) 'tpr     ',   'total precipitation(mm/day)'
               write(31,600) 'tgb     ',   'lower groud temp. in BATS  '
               write(31,600) 'swt     ',   'total soil water in mm H2O '
               write(31,600) 'rno     ',   'accumulated infiltration   '
               write(31,700)
               close(31)
            endif
         enddo
            
         deallocate(I1UR)
         deallocate(J1UR)
         deallocate(I1UL)
         deallocate(J1UL)
         deallocate(I1DR)
         deallocate(J1DR)
         deallocate(I1DL)
         deallocate(J1DL)
         deallocate(D1XT)
         deallocate(D1Xa)
         deallocate(D1Xb)
         deallocate(D1Xc)
         deallocate(D1Xd)
      endif
      deallocate(out)
      deallocate(olat)
      deallocate(olon)

      deallocate(sigma)
      deallocate(o)
      deallocate(xlat)
      deallocate(xlon)

  10  format('dset ^',A17)
  11  format('dset ^',A15)
  12  format('dset ^',A10)
  20  format('title RegCM domain information')
  30  format('options big_endian')
  40  format('options little_endian')
  50  format('undef -9999.')
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
 210  format('ydef ',I3,' linear ',f9.4,' ',f9.4)
 300  format('zdef ',I2,' levels ',30f7.2)
 400  format('tdef ',I4,' linear ',I2,'z',A2,A3,I4,' ',I2,'hr')
 401  format('tdef ',I4,' linear ',A2,A3,I4,' ','1dy')
 402  format('tdef ',I4,' linear ',A2,A3,I4,' ','1mo')
 500  format('vars ',I2)
 600  format(A8,'0 99 ',A26)           
 650  format(A8,I2,' 0 ',A26)
 700  format('endvars')
      return
      end

      subroutine regrid_SRF(iy,jx,kz,ibyte,idate0,idate1,idate2 &
                           ,Path_Output)
      implicit none
      integer iy,jx,kz,ibyte,idate0,idate1,idate2
      character*128 Path_Output
      real*4  truelatL,truelatH
      integer iiy,jjx,kkz
      integer mdate0,ibltyp,icup,ipptls,iboudy
      real*4  dxsp,ptsp,clat,clon,plat,plon
      real*4  dto,dtb,dtr,dtc
      integer iotyp
      character*6 iproj
      real*4, allocatable ::  sigma(:)
      character*4 :: chy(1941:2100)
      data chy/'1941','1942','1943','1944','1945','1946','1947','1948',&
        '1949','1950','1951','1952','1953','1954','1955','1956','1957',&
        '1958','1959','1960','1961','1962','1963','1964','1965','1966',&
        '1967','1968','1969','1970','1971','1972','1973','1974','1975',&
        '1976','1977','1978','1979','1980','1981','1982','1983','1984',&
        '1985','1986','1987','1988','1989','1990','1991','1992','1993',&
        '1994','1995','1996','1997','1998','1999','2000','2001','2002',&
        '2003','2004','2005','2006','2007','2008','2009','2010','2011',&
        '2012','2013','2014','2015','2016','2017','2018','2019','2020',&
        '2021','2022','2023','2024','2025','2026','2027','2028','2029',&
        '2030','2031','2032','2033','2034','2035','2036','2037','2038',&
        '2039','2040','2041','2042','2043','2044','2045','2046','2047',&
        '2048','2049','2050','2051','2052','2053','2054','2055','2056',&
        '2057','2058','2059','2060','2061','2062','2063','2064','2065',&
        '2066','2067','2068','2069','2070','2071','2072','2073','2074',&
        '2075','2076','2077','2078','2079','2080','2081','2082','2083',&
        '2084','2085','2086','2087','2088','2089','2090','2091','2092',&
        '2093','2094','2095','2096','2097','2098','2099','2100'/
      character*2 chm(12)
      data chm/'01','02','03','04','05','06','07','08','09','10', &
               '11','12'/
      character*3 chmc(12)
      data chmc/'jan','feb','mar','apr','may','jun'  &
               ,'jul','aug','sep','oct','nov','dec'/
      character*14 filein
      character*17 fileout
      integer ntype,nfile,nyear,month,n_slice,mrec,nrec,nnn

      integer i,j,k,l,m,n
      integer igrads,ibigend

      real*4  size_2d(4),ddeg
      COMMON /WINDOW/ size_2d,ddeg

      real*4, allocatable ::  o(:,:,:),xlat(:,:),xlon(:,:)
      real*4  xmaxlat,xminlat,xmaxlon,xminlon

      integer mlat,nlat,mlon,nlon
      integer n_month
      logical there

      integer MUR,NUR,MUL,NUL,MDR,NDR,MDL,NDL
      real*4  DISTa,DISTb,DISTc,DISTd,AAA,DIST
      real*4  PI,PIR180
      real*4  POLLAM,POLPHI,POLCPHI,POLSPHI,ZPHI,ZRLA,ZRLAP
      real*4  ZARG1,ZARG2,ZNORM,SINDEL,COSDEL,US,VS,X
      real*4  SIGN0,GRIDFC
      real*4  sum0,weight
      real*4, allocatable :: out(:,:)
      real*4, allocatable :: olat(:),olon(:)

      integer, allocatable :: I1UR(:,:), J1UR(:,:)
      integer, allocatable :: I1UL(:,:), J1UL(:,:)
      integer, allocatable :: I1DR(:,:), J1DR(:,:)
      integer, allocatable :: I1DL(:,:), J1DL(:,:)
      real*4, allocatable :: D1XT(:,:), D1Xa(:,:), D1Xb(:,:)
      real*4, allocatable ::            D1Xc(:,:), D1Xd(:,:)

      igrads = 1
      ibigend= 1

      PI = atan(1.)*4.
      PIR180 = atan(1.)/45.

      allocate(sigma(kz+1))
      allocate(o(jx-2,iy-2,27))
      allocate(xlat(jx-2,iy-2))
      allocate(xlon(jx-2,iy-2))

      inquire(file=trim(Path_Output)//'OUT_HEAD',exist=there)
      if(.not.there) then
         write(*,*) trim(Path_Output)//'OUT_HEAD',' is not avaiable'
         stop
      endif
      open(10,file=trim(Path_Output)//'OUT_HEAD',form='unformatted' &
             ,recl=(jx-2)*(iy-2)*ibyte,access='direct')
      read(10,rec=1) mdate0,ibltyp,icup,ipptls,iboudy  &
                    ,iiy,jjx,kkz,(sigma(k),k=1,kz+1)   &
                    ,dxsp,ptsp,clat,clon,plat,plon          &
                    ,iproj,dto,dtb,dtr,dtc,iotyp,truelatL,truelatH
      if(iiy.ne.iy.or.jjx.ne.jx.or.kkz.ne.kz) then
         write(*,*) 'iy,jx,kz in parameter = ',iy,jx,kz
         write(*,*) 'iy,jx,kz in OUT_HEAD ',iiy,jjx,kkz
         write(*,*) 'They are not consistent'
         stop
      endif
      read(10,rec=6) xlat         ! xlat
      if(size_2d(3).gt.-9990.) then
         xminlat = size_2d(3)
         xmaxlat = size_2d(4)
      else
         xmaxlat = -1.e10
         xminlat =  1.e10
         do i=1,iy-2
         do j=1,jx-2
            if(xlat(j,i).gt.xmaxlat) xmaxlat=xlat(j,i)
            if(xlat(j,i).lt.xminlat) xminlat=xlat(j,i)
         enddo
         enddo
      endif
      read(10,rec=7) xlon         ! xlon
      if(size_2d(1).gt.-9990.) then
         xminlon=size_2d(1)
         xmaxlon=size_2d(2)
      else
         xmaxlon = -1.e10
         xminlon =  1.e10
         do i=1,iy-2
         do j=1,jx-2
            if(xlon(j,i).gt.xmaxlon) xmaxlon=xlon(j,i)
            if(xlon(j,i).lt.xminlon) xminlon=xlon(j,i)
         enddo
         enddo
      endif
      close(10)
      
      mlat = xmaxlat/ddeg - 1
      nlat = xminlat/ddeg + 1
      mlon = xmaxlon/ddeg - 1
      nlon = xminlon/ddeg + 2

      write(*,*) xminlat,mlat-nlat+1,(nlat-1)*ddeg
      write(*,*) xminlon,mlon-nlon+1,(nlon-1)*ddeg

      allocate(out(nlon:mlon,nlat:mlat))
      allocate(olat(nlat:mlat))
      allocate(olon(nlon:mlon))
      do n=nlat,mlat
         olat(n) = real(n-1)*ddeg
      enddo
      do m=nlon,mlon
         olon(m) = real(m-1)*ddeg
      enddo

      if(idate1.gt.1940010100) then        ! original file
         n_month = (idate2/1000000-idate1/1000000)*12 &
                 + (mod(idate2/10000,100)-mod(idate1/10000,100))   
         if(mod(idate2,10000).gt.0100) n_month = n_month+1
         ntype = 0
      else if(idate1.gt.19400101) then     ! daily mean file
         n_month = (idate2/10000-idate1/10000)*12 &
                 + (mod(idate2/100,100)-mod(idate1/100,100))   
         if(mod(idate2,100).gt.01) n_month = n_month+1
         ntype = 1
      else if(idate1.gt.194001) then       ! monthly mean file
         n_month = (idate2/100-idate1/100)*12 &
                 + (mod(idate2,100)-mod(idate1,100))+1
         ntype = 2
      endif

      if(iproj.eq.'NORMER') then
         do nfile=1,n_month
            if(ntype.eq.0) then
               nyear = idate1/1000000
               month = idate1/10000-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(nfile.eq.1.or.month.eq.1) then
                  write(*,*)'regrid: NORMER SRF Orig.',nyear,month
               else
                  write(*,*)'                        ',nyear,month
               endif

               if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
                  month.eq.7.or.month.eq.8.or.month.eq.10.or. &
                  month.eq.12) then
                  n_slice = 31*24/dtb
               else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                       month.eq.11) then
                  n_slice = 30*24/dtb
               else
                  n_slice = 28*24/dtb
                  if(mod(nyear,4).eq.0) n_slice = 29*24/dtb
                  if(mod(nyear,100).eq.0) n_slice = 28*24/dtb
                  if(mod(nyear,400).eq.0) n_slice = 29*24/dtb
               endif
               if(nfile.eq.n_month.and.mod(idate2/100-1,100).ne.0) &
               n_slice=min(n_slice,nint(mod(idate2/100-1,100)*24./dtb))
               if(nfile.eq.1.and.idate0.eq.idate1) n_slice=n_slice+1
               filein = 'SRF.'//chy(nyear)//chm(month)//'0100'
               fileout= 'SRF_LL.'//chy(nyear)//chm(month)//'0100'
               inquire(file=trim(Path_Output)//filein,exist=there)
               if(.not.there) then
                  write(*,*) trim(Path_Output)//filein,' is not avaiable'
                  stop
               endif
             open(10,file=trim(Path_Output)//filein,form='unformatted' &
                      ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
               mrec = 0
            open(20,file=trim(Path_Output)//fileout,form='unformatted' &
               ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
               nrec = 0
            else if(ntype.eq.1) then
               nyear = idate1/10000
               month = idate1/100-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(nfile.eq.1.or.month.eq.1) then
                  write(*,*)'regrid: NORMER SRF Daily',nyear,month
               else
                  write(*,*)'                        ',nyear,month
               endif
               if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
                  month.eq.7.or.month.eq.8.or.month.eq.10.or. &
                  month.eq.12) then
                  n_slice = 31
               else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                       month.eq.11) then
                  n_slice = 30
               else
                  n_slice = 28
                  if(mod(nyear,4).eq.0) n_slice = 29
                  if(mod(nyear,100).eq.0) n_slice = 28
                  if(mod(nyear,400).eq.0) n_slice = 29
               endif
               if(nfile.eq.n_month.and.mod(idate2-1,100).ne.0) &
               n_slice=min(n_slice,mod(idate2-1,100)+1)
               filein = 'SRF.'//chy(nyear)//chm(month)//'01'
               fileout= 'SRF_LL.'//chy(nyear)//chm(month)//'01'
               inquire(file=trim(Path_Output)//filein,exist=there)
               if(.not.there) then
                  write(*,*) trim(Path_Output)//filein,' is not avaiable'
                  stop
               endif
             open(10,file=trim(Path_Output)//filein,form='unformatted' &
                      ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
               mrec = 0
            open(20,file=trim(Path_Output)//fileout,form='unformatted' &
               ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
               nrec = 0
            else if(ntype.eq.2) then
               nyear = idate1/100
               month = idate1-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(month.eq.1.or.nfile.eq.1) then
                  write(*,*)'regrid: NORMER SRF Month',nyear,month
               else
                  write(*,*)'                        ',nyear,month
               endif

               n_slice = 1

               if(nfile.eq.1) then
                  filein = 'SRF.mon'
                  fileout= 'SRF_LL.mon'
                inquire(file=trim(Path_Output)//filein(1:7),exist=there)
                  if(.not.there) then
                     write(*,*) trim(Path_Output)//filein(1:7)  &
                               ,' is not avaiable'
                     stop
                  endif
                  open(10,file=trim(Path_Output)//filein(1:7)  &
                         ,form='unformatted' &
                         ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
                  nrec = 0
                  open(20,file=trim(Path_Output)//fileout(1:10)  &
                         ,form='unformatted' &
                ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
                     mrec = 0
               endif
            endif

            do nnn=1,n_slice
               do l=1,27
                  mrec = mrec+1
                  read(10,rec=mrec) ((o(j,i,l),j=1,jx-2),i=1,iy-2)
               enddo
               do l=1,27
                  do n=nlat,mlat
                  do m=nlon,mlon
                     out(m,n) = -9999.

                     do i=1,iy-3
                     do j=1,jx-3
                        if((xlat(j,i).lt.olat(n).and.      &
                            xlon(j,i).lt.olon(m))          &
                      .and.(xlat(j+1,i+1).ge.olat(n).and.  &
                            xlon(j+1,i+1).ge.olon(m))) then
                      sum0 = 0.0
                      weight = 0.0
                      if(o(j,i,l).gt.-9990.) then
            sum0 = sum0       &
            + o(j,i,l)*(xlon(j+1,i+1)-olon(m))*(xlat(j+1,i+1)-olat(n))
            weight = weight   &
            + (xlon(j+1,i+1)-olon(m))*(xlat(j+1,i+1)-olat(n))
                      endif
                      if(o(j+1,i,l).gt.-9990.) then
            sum0 = sum0       &
            + o(j+1,i,l)*(olon(m)-xlon(j,i+1))*(xlat(j+1,i+1)-olat(n))
            weight = weight   &
            + (olon(m)-xlon(j,i+1))*(xlat(j+1,i+1)-olat(n))
                      endif
                      if(o(j,i+1,l).gt.-9990.) then
            sum0 = sum0       &
            + o(j,i+1,l)*(xlon(j+1,i)-olon(m))*(olat(n)-xlat(j+1,i))
            weight = weight   &
            + (xlon(j+1,i)-olon(m))*(olat(n)-xlat(j+1,i))
                      endif
                      if(o(j,i,l).gt.-9990.) then
            sum0 = sum0       &
            + o(j+1,i+1,l)*(olon(m)-xlon(j,i))*(olat(n)-xlat(j,i))
            weight = weight   &
            + (olon(m)-xlon(j,i))*(olat(n)-xlat(j,i))
                      endif
                           if(weight.gt.1.e-10) then
                              out(m,n) = sum0/weight
                           endif
                        endif
                     enddo
                     enddo
                  enddo
                  enddo
                  nrec = nrec+1
                  write(20,rec=nrec)((out(m,n),m=nlon,mlon),n=nlat,mlat)
               enddo
            enddo
            if(.not.ntype.eq.2) close(20)
            if(.not.ntype.eq.2) close(10)
            if(igrads.eq.1.and.(ntype.eq.0.or.ntype.eq.1.or. &
               (ntype.eq.2.and.nfile.eq.1))) then
               if(ntype.eq.0) then
                 inquire(file=trim(Path_Output)//fileout(1:17)//'.ctl' &
                        ,exist=there)
                  if(there) then
                     open(31,file= &
                          trim(Path_Output)//fileout(1:17)//'.ctl' &
                         ,form='formatted',status='replace')
                  else
                     open(31,file= &
                          trim(Path_Output)//fileout(1:17)//'.ctl' &
                         ,form='formatted',status='new')
                  endif
                  write(31,10) fileout(1:17)
               else if(ntype.eq.1) then
                 inquire(file=trim(Path_Output)//fileout(1:15)//'.ctl' &
                        ,exist=there)
                  if(there) then
                     open(31,file= &
                          trim(Path_Output)//fileout(1:15)//'.ctl' &
                         ,form='formatted',status='replace')
                  else
                     open(31,file= &
                          trim(Path_Output)//fileout(1:15)//'.ctl' &
                         ,form='formatted',status='new')
                  endif
                  write(31,11) fileout(1:15)
               else if(ntype.eq.2) then
                 inquire(file=trim(Path_Output)//fileout(1:10)//'.ctl' &
                        ,exist=there)
                  if(there) then
                     open(31,file= &
                          trim(Path_Output)//fileout(1:10)//'.ctl' &
                         ,form='formatted',status='replace')
                  else
                     open(31,file= &
                          trim(Path_Output)//fileout(1:10)//'.ctl' &
                         ,form='formatted',status='new')
                  endif
                  write(31,12) fileout(1:10)
               endif
               write(31,20)
               if(ibigend.eq.1) then
                  write(31,30)
               else
                  write(31,40)
               endif
               write(31,50)
               write(31,200) mlon-nlon+1,(nlon-1)*ddeg,ddeg
               write(31,210) mlat-nlat+1,(nlat-1)*ddeg,ddeg
               write(31,300) &
            (1013.25-ptsp*10.)*(sigma(kz)+sigma(kz+1))*0.5+ptsp*10.
               if(ntype.eq.0) then
               if(nfile.eq.1.and.idate0.eq.idate1) then
               write(31,400) n_slice,0,'01',chmc(month),nyear,nint(dtb)
               else
       write(31,400) n_slice,nint(dtb),'01',chmc(month),nyear,nint(dtb)
               endif
               else if(ntype.eq.1) then
                  write(31,401) n_slice,'01',chmc(month),nyear
               else if(ntype.eq.2) then
                  write(31,402) n_month,'16',chmc(month),nyear
               endif
               write(31,500) 27
        write(31,600) 'u10m    ','westerly  wind at 10m (m/s)          '
        write(31,600) 'v10m    ','southerly wind at 10m (m/s)          '
        write(31,600) 'uvdrag  ','surface drag stress                  '
        write(31,600) 'tg      ','ground temperature (degree)          '
        write(31,600) 'tlef    ','temperature of foliage               '
        write(31,600) 't2m     ','air temperature at 2m (K)            '
        write(31,600) 'q2m     ','water vapor mixing ratio at 2m(kg/kg)'
        write(31,600) 'ssw     ','upper layer soil water               '
        write(31,600) 'rsw     ','root zone soil water                 '
        write(31,600) 'tpr     ','total precipitation (mm/day)         '
        write(31,600) 'evp     ','evapotranspiration (mm/day)          '
        write(31,600) 'runoff  ','surface runoff (mm/day)              '
        write(31,600) 'scv     ','snow amount (mm, water equivalent)   '
        write(31,600) 'sena    ','sensible heat flux (W/m2)            '
        write(31,600) 'flw     ','net infrared energy flux (W/m2)      '
        write(31,600) 'fsw     ','net absorbed solar energy flux (W/m2)'
        write(31,600) 'flwd    ','downward infrared energy flux (W/m2) '
        write(31,600) 'sina    ','incident solar energy flux (W/m2)    '
        write(31,600) 'prcv    ','convective precipitation (mm/day)    '
        write(31,600) 'ps      ','surface pressure (hPa)               '
        write(31,600) 'zpbl    ','PBL layer height                     '
        write(31,600) 'tgmax   ','maximum ground temperature (K)       '
        write(31,600) 'tgmin   ','minimum ground temperature (K)       '
        write(31,600) 't2max   ','maximum 2m air temperature (K)       '
        write(31,600) 't2min   ','minimum 2m air temperature (K)       '
        write(31,600) 'w10max  ','maximum 10m wind speed (m/s)         '
        write(31,600) 'ps_min  ','minimum surface pressure (hPa)       '
               write(31,700)
               close(31)
            endif
         enddo

      else if(iproj.eq.'LAMCON'.or.iproj.eq.'ROTMER') then

         allocate(I1UR(nlon:mlon,nlat:mlat))
         allocate(J1UR(nlon:mlon,nlat:mlat))
         allocate(I1UL(nlon:mlon,nlat:mlat))
         allocate(J1UL(nlon:mlon,nlat:mlat))
         allocate(I1DR(nlon:mlon,nlat:mlat))
         allocate(J1DR(nlon:mlon,nlat:mlat))
         allocate(I1DL(nlon:mlon,nlat:mlat))
         allocate(J1DL(nlon:mlon,nlat:mlat))
         allocate(D1XT(nlon:mlon,nlat:mlat))
         allocate(D1Xa(nlon:mlon,nlat:mlat))
         allocate(D1Xb(nlon:mlon,nlat:mlat))
         allocate(D1Xc(nlon:mlon,nlat:mlat))
         allocate(D1Xd(nlon:mlon,nlat:mlat))

         do n=nlat,mlat
         do m=nlon,mlon

            MUR=1000
            NUR=1000
            MUL=1000
            NUL=1000
            MDR=1000
            NDR=1000
            MDL=1000
            NDL=1000
            DISTa=1.E8
            DISTb=1.E8
            DISTc=1.E8
            DISTd=1.E8

            do i=1,iy-2
            do j=1,jx-2
         IF((xlon(j,i).ge.olon(m).and.xlon(j,i)-olon(m).lt.10.) .and. &
            (xlat(j,i).ge.olat(n).and.xlat(j,i)-olat(n).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTa.gt.AAA) then
                DISTa = AAA
                MUR = j
                NUR = i
             endif
         ENDIF
         IF((xlon(j,i).lt.olon(m).and.olon(m)-xlon(j,i).lt.10.) .and. &
            (xlat(j,i).ge.olat(n).and.xlat(j,i)-olat(n).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTb.gt.AAA) then
                DISTb = AAA
                MUL = j
                NUL = i
             endif
         ENDIF
         IF((xlon(j,i).ge.olon(m).and.xlon(j,i)-olon(m).lt.10.) .and. &
            (xlat(j,i).lt.olat(n).and.olat(n)-xlat(j,i).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTc.gt.AAA) then
                DISTc = AAA
                MDR = j
                NDR = i
             endif
         ENDIF
         IF((xlon(j,i).lt.olon(m).and.olon(m)-xlon(j,i).lt.10.) .and. &
            (xlat(j,i).lt.olat(n).and.olat(n)-xlat(j,i).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTd.gt.AAA) then
                DISTd = AAA
                MDL = j
                NDL = i
             endif
         ENDIF
            enddo
            enddo
            
            DIST=amin1(DISTa,DISTb,DISTc,DISTd)
            I1UR(m,n) = MUR
            J1UR(m,n) = NUR
            I1UL(m,n) = MUL
            J1UL(m,n) = NUL
            I1DR(m,n) = MDR
            J1DR(m,n) = NDR
            I1DL(m,n) = MDL
            J1DL(m,n) = NDL
            D1XT(m,n) = DIST
            D1Xa(m,n) = DISTa
            D1Xb(m,n) = DISTb
            D1Xc(m,n) = DISTc
            D1Xd(m,n) = DISTd
         enddo
         enddo

         do nfile=1,n_month
            if(ntype.eq.0) then
               nyear = idate1/1000000
               month = idate1/10000-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(iproj.eq.'LAMCON') then
                  if(nfile.eq.1.or.month.eq.1) then
                     write(*,*)'regrid: LAMCON SRF Orig.',nyear,month
                  else
                     write(*,*)'                        ',nyear,month
                  endif
               else if(iproj.eq.'ROTMER') then
                  if(nfile.eq.1.or.month.eq.1) then
                     write(*,*)'regrid: ROTMER SRF Orig.',nyear,month
                  else
                     write(*,*)'                        ',nyear,month
                  endif
               endif

               if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
                  month.eq.7.or.month.eq.8.or.month.eq.10.or. &
                  month.eq.12) then
                  n_slice = 31*24/dtb
               else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                       month.eq.11) then
                  n_slice = 30*24/dtb
               else
                  n_slice = 28*24/dtb
                  if(mod(nyear,4).eq.0) n_slice = 29*24/dtb
                  if(mod(nyear,100).eq.0) n_slice = 28*24/dtb
                  if(mod(nyear,400).eq.0) n_slice = 29*24/dtb
               endif
               if(nfile.eq.n_month.and.mod(idate2/100-1,100).ne.0) &
               n_slice=min(n_slice,nint(mod(idate2/100-1,100)*24./dtb))
               if(nfile.eq.1.and.idate0.eq.idate1) n_slice=n_slice+1
               filein = 'SRF.'//chy(nyear)//chm(month)//'0100'
               fileout= 'SRF_LL.'//chy(nyear)//chm(month)//'0100'
               inquire(file=trim(Path_Output)//filein,exist=there)
               if(.not.there) then
                  write(*,*) trim(Path_Output)//filein,' is not avaiable'
                  stop
               endif
             open(10,file=trim(Path_Output)//filein,form='unformatted' &
                      ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
               mrec = 0
            open(20,file=trim(Path_Output)//fileout,form='unformatted' &
               ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
               nrec = 0
            else if(ntype.eq.1) then
               nyear = idate1/10000
               month = idate1/100-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(iproj.eq.'LAMCON') then
                  if(nfile.eq.1.or.month.eq.1) then
                     write(*,*)'regrid: LAMCON SRF Daily',nyear,month
                  else
                     write(*,*)'                        ',nyear,month
                  endif
               else if(iproj.eq.'ROTMER') then
                  if(nfile.eq.1.or.month.eq.1) then
                     write(*,*)'regrid: ROTMER SRF Daily',nyear,month
                  else
                     write(*,*)'                        ',nyear,month
                  endif
               endif

               if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
                  month.eq.7.or.month.eq.8.or.month.eq.10.or. &
                  month.eq.12) then
                  n_slice = 31
               else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                       month.eq.11) then
                  n_slice = 30
               else
                  n_slice = 28
                  if(mod(nyear,4).eq.0) n_slice = 29
                  if(mod(nyear,100).eq.0) n_slice = 28
                  if(mod(nyear,400).eq.0) n_slice = 29
               endif
               if(nfile.eq.n_month.and.mod(idate2-1,100).ne.0) &
               n_slice=min(n_slice,mod(idate2-1,100)+1)
               filein = 'SRF.'//chy(nyear)//chm(month)//'01'
               fileout= 'SRF_LL.'//chy(nyear)//chm(month)//'01'
               inquire(file=trim(Path_Output)//filein,exist=there)
               if(.not.there) then
                  write(*,*) trim(Path_Output)//filein,' is not avaiable'
                  stop
               endif
             open(10,file=trim(Path_Output)//filein,form='unformatted' &
                      ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
               mrec = 0
            open(20,file=trim(Path_Output)//fileout,form='unformatted' &
               ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
               nrec = 0
            else if(ntype.eq.2) then
               nyear = idate1/100
               month = idate1-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(iproj.eq.'LAMCON') then
                  if(month.eq.1.or.nfile.eq.1) then
                     write(*,*)'regrid: LAMCON SRF Month',nyear,month
                  else
                     write(*,*)'                        ',nyear,month
                  endif
               else if(iproj.eq.'ROTMER') then
                  if(month.eq.1.or.nfile.eq.1) then
                     write(*,*)'regrid: ROTMER SRF Month',nyear,month
                  else
                     write(*,*)'                        ',nyear,month
                  endif
               endif

               n_slice = 1

               if(nfile.eq.1) then
                  filein = 'SRF.mon'
                  fileout= 'SRF_LL.mon'
                inquire(file=trim(Path_Output)//filein(1:7),exist=there)
                  if(.not.there) then
                     write(*,*) trim(Path_Output)//filein(1:7)  &
                               ,' is not avaiable'
                     stop
                  endif
                  open(10,file=trim(Path_Output)//filein(1:7)  &
                         ,form='unformatted' &
                         ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
                  nrec = 0
                  open(20,file=trim(Path_Output)//fileout(1:10)  &
                         ,form='unformatted' &
                ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
                  mrec = 0
               endif
            endif

            do nnn=1,n_slice
               do l=1,27
                  mrec = mrec+1
                  read(10,rec=mrec) ((o(j,i,l),j=1,jx-2),i=1,iy-2)
               enddo
               IF(iproj.eq.'ROTMER') THEN
                  IF(PLAT.GT.0.) THEN
                     POLLAM = PLON + 180.
                     POLPHI = 90. - PLAT
                  ELSE
                     POLLAM = PLON
                     POLPHI = 90. + PLAT
                  ENDIF
                  IF(POLLAM.GT.180.) POLLAM = POLLAM - 360.
                  POLCPHI = cos(PIR180*POLPHI)
                  POLSPHI = sin(PIR180*POLPHI)
                  do i=1,iy-2
                  do j=1,jx-2
                     ZPHI = XLAT(j,i)*PIR180
                     ZRLA = XLON(j,i)*PIR180
                     IF(XLAT(j,i).gt.89.999999) ZRLA = 0.0
                     ZRLAP = POLLAM*PIR180 - ZRLA
                     ZARG1  = POLCPHI*sin(ZRLAP)
                     ZARG2  = POLSPHI*cos(ZPHI)  &
                            - POLCPHI*sin(ZPHI)*cos(ZRLAP)
                     ZNORM  = 1.0/sqrt(ZARG1**2+ZARG2**2)
                     SINDEL = ZARG1*ZNORM
                     COSDEL = ZARG2*ZNORM
                     US =  o(j,i,1)*COSDEL + o(j,i,2)*SINDEL
                     VS = -o(j,i,1)*SINDEL + o(j,i,2)*COSDEL
                     o(j,i,1) = US
                     o(j,i,2) = VS
                  enddo
                  enddo
               ELSE IF(iproj.eq.'LAMCON') THEN
                  IF(CLAT.lt.0.) THEN
                     SIGN0= -1.           ! SOUTH HEMESPHERE
                  ELSE
                     SIGN0=  1.           ! NORTH HEMESPHERE
                  ENDIF
                  IF(abs(truelatL-truelatH).gt.0.1) THEN
                     GRIDFC=(alog10(cos(truelatL*PIR180)) &
                            -alog10(cos(truelatH*PIR180)))  &
               /(alog10(tan(45.0-SIGN0*truelatL/2.0*PIR180))  &
                -alog10(tan(45.0-SIGN0*truelatH/2.0*PIR180)))
                  ELSE
                     GRIDFC=SIGN0*sin(truelatL*PIR180)
                  ENDIF
                  do i=1,iy-2
                  do j=1,jx-2
                     IF((CLON.ge.0.0.and.XLON(j,i).ge.0.).or. &
                        (CLON.lt.0.0.and.XLON(j,i).lt.0.)) THEN
                        X=(CLON-XLON(j,i))*PIR180*GRIDFC
                     ELSE
                        IF(CLON.ge.0.0) THEN
                           IF(abs(CLON-(XLON(j,i)+360.)).lt.   &
                              abs(CLON- XLON(j,i)) ) THEN
                              X=(CLON-(XLON(j,i)+360.))*PIR180*GRIDFC
                           ELSE
                              X=(CLON-XLON(j,i))*PIR180*GRIDFC
                           ENDIF
                        ELSE
                           IF(abs(CLON-(XLON(j,i)-360.)).LT. &
                              abs(CLON- XLON(j,i)) ) THEN
                              X=(CLON-(XLON(j,i)-360.))*PIR180*GRIDFC
                           ELSE
                              X=(CLON-XLON(j,i))*PIR180*GRIDFC
                           ENDIF
                        ENDIF
                     ENDIF
                     SINDEL=sin(X)
                     COSDEL=cos(X)
                     IF(CLAT.ge.0.) THEN
                        US= o(j,i,1)*COSDEL - o(j,i,2)*SINDEL
                        VS=-o(j,i,1)*SINDEL + o(j,i,2)*COSDEL
                        o(j,i,1) = US
                        o(j,i,2) = VS
                     ELSE
                        US= o(j,i,1)*COSDEL + o(j,i,2)*SINDEL
                        VS=-o(j,i,1)*SINDEL + o(j,i,2)*COSDEL
                        o(j,i,1) = US
                        o(j,i,2) = VS
                     ENDIF
                  enddo
                  enddo
               ENDIF
               do l=1,27
                  do n=nlat,mlat
                  do m=nlon,mlon
                     out(m,n) = -9999.

            if(I1UR(m,n).lt.999.and.J1UR(m,n).lt.999.and. &
               I1UL(m,n).lt.999.and.J1UL(m,n).lt.999.and. &
               I1DR(m,n).lt.999.and.J1DR(m,n).lt.999.and. &
               I1DL(m,n).lt.999.and.J1DL(m,n).lt.999) then
               if(D1XT(m,n).gt.0.0001) then
                  sum0=0.0
                  weight = 0.0
                  if(o(I1UR(m,n),J1UR(m,n),l).gt.-9990.) then
                     sum0 = sum0 + o(I1UR(m,n),J1UR(m,n),l)/D1Xa(m,n)
                     weight = weight + 1./D1Xa(m,n)
                  endif
                  if(o(I1UL(m,n),J1UL(m,n),l).gt.-9990.) then
                     sum0 = sum0 + o(I1UL(m,n),J1UL(m,n),l)/D1Xb(m,n)
                     weight = weight + 1./D1Xb(m,n)
                  endif
                  if(o(I1DR(m,n),J1DR(m,n),l).gt.-9990.) then
                     sum0 = sum0 + o(I1DR(m,n),J1DR(m,n),l)/D1Xc(m,n)
                     weight = weight + 1./D1Xc(m,n)
                  endif
                  if(o(I1DL(m,n),J1DL(m,n),l).gt.-9990.) then
                     sum0 = sum0 + o(I1DL(m,n),J1DL(m,n),l)/D1Xd(m,n)
                     weight = weight + 1./D1Xd(m,n)
                  endif
                  if(weight.gt.1.e-10) then
                     out(m,n) = sum0/weight
                  endif
               else
                  if(D1Xa(m,n).eq.D1XT(m,n).and. &
                     o(I1UR(m,n),J1UR(m,n),l).gt.-9990.) then
                     out(m,n) = o(I1UR(m,n),J1UR(m,n),l)
                  else if(D1Xb(m,n).eq.D1XT(m,n).and. &
                          o(I1UL(m,n),J1UL(m,n),l).gt.-9990.) then
                     out(m,n) = o(I1UL(m,n),J1UL(m,n),l)
                  else if(D1Xc(m,n).eq.D1XT(m,n).and. &
                          o(I1DR(m,n),J1DR(m,n),l).gt.-9990.) then
                     out(m,n) = o(I1DR(m,n),J1DR(m,n),l)
                  else if(D1Xd(m,n).eq.D1XT(m,n).and. &
                          o(I1DL(m,n),J1DL(m,n),l).gt.-9990.) then
                     out(m,n) = o(I1DL(m,n),J1DL(m,n),l)
                  endif
               endif
            endif
                  enddo
                  enddo
                  nrec=nrec+1
                  write(20,rec=nrec)((out(m,n),m=nlon,mlon),n=nlat,mlat)
               enddo
            enddo
            if(.not.ntype.eq.2) close(20)
            if(.not.ntype.eq.2) close(10)
            if(igrads.eq.1.and.(ntype.eq.0.or.ntype.eq.1.or. &
               (ntype.eq.2.and.nfile.eq.1))) then
               if(ntype.eq.0) then
                 inquire(file=trim(Path_Output)//fileout(1:17)//'.ctl' &
                        ,exist=there)
                  if(there) then
                     open(31,file= &
                          trim(Path_Output)//fileout(1:17)//'.ctl' &
                         ,form='formatted',status='replace')
                  else
                     open(31,file= &
                          trim(Path_Output)//fileout(1:17)//'.ctl' &
                         ,form='formatted',status='new')
                  endif
                  write(31,10) fileout(1:17)
               else if(ntype.eq.1) then
                 inquire(file=trim(Path_Output)//fileout(1:15)//'.ctl' &
                        ,exist=there)
                  if(there) then
                     open(31,file= &
                          trim(Path_Output)//fileout(1:15)//'.ctl' &
                         ,form='formatted',status='replace')
                  else
                     open(31,file= &
                          trim(Path_Output)//fileout(1:15)//'.ctl' &
                         ,form='formatted',status='new')
                  endif
                  write(31,11) fileout(1:15)
               else if(ntype.eq.2) then
                 inquire(file=trim(Path_Output)//fileout(1:10)//'.ctl' &
                        ,exist=there)
                  if(there) then
                     open(31,file= &
                          trim(Path_Output)//fileout(1:10)//'.ctl' &
                         ,form='formatted',status='replace')
                  else
                     open(31,file= &
                          trim(Path_Output)//fileout(1:10)//'.ctl' &
                         ,form='formatted',status='new')
                  endif
                  write(31,12) fileout(1:10)
               endif
               write(31,20)
               if(ibigend.eq.1) then
                  write(31,30)
               else
                  write(31,40)
               endif
               write(31,50)
               write(31,200) mlon-nlon+1,(nlon-1)*ddeg,ddeg
               write(31,210) mlat-nlat+1,(nlat-1)*ddeg,ddeg
               write(31,300) &
            (1013.25-ptsp*10.)*(sigma(kz)+sigma(kz+1))*0.5+ptsp*10.
               if(ntype.eq.0) then
               if(nfile.eq.1.and.idate0.eq.idate1) then
               write(31,400) n_slice,0,'01',chmc(month),nyear,nint(dtb)
               else
       write(31,400) n_slice,nint(dtb),'01',chmc(month),nyear,nint(dtb)
               endif
               else if(ntype.eq.1) then
                  write(31,401) n_slice,'01',chmc(month),nyear
               else if(ntype.eq.2) then
                  write(31,402) n_month,'16',chmc(month),nyear
               endif
               write(31,500) 27
        write(31,600) 'u10m    ','westerly  wind at 10m (m/s)          '
        write(31,600) 'v10m    ','southerly wind at 10m (m/s)          '
        write(31,600) 'uvdrag  ','surface drag stress                  '
        write(31,600) 'tg      ','ground temperature (degree)          '
        write(31,600) 'tlef    ','temperature of foliage               '
        write(31,600) 't2m     ','air temperature at 2m (K)            '
        write(31,600) 'q2m     ','water vapor mixing ratio at 2m(kg/kg)'
        write(31,600) 'ssw     ','upper layer soil water               '
        write(31,600) 'rsw     ','root zone soil water                 '
        write(31,600) 'tpr     ','total precipitation (mm/day)         '
        write(31,600) 'evp     ','evapotranspiration (mm/day)          '
        write(31,600) 'runoff  ','surface runoff (mm/day)              '
        write(31,600) 'scv     ','snow amount (mm, water equivalent)   '
        write(31,600) 'sena    ','sensible heat flux (W/m2)            '
        write(31,600) 'flw     ','net infrared energy flux (W/m2)      '
        write(31,600) 'fsw     ','net absorbed solar energy flux (W/m2)'
        write(31,600) 'flwd    ','downward infrared energy flux (W/m2) '
        write(31,600) 'sina    ','incident solar energy flux (W/m2)    '
        write(31,600) 'prcv    ','convective precipitation (mm/day)    '
        write(31,600) 'ps      ','surface pressure (hPa)               '
        write(31,600) 'zpbl    ','PBL layer height                     '
        write(31,600) 'tgmax   ','maximum ground temperature (K)       '
        write(31,600) 'tgmin   ','minimum ground temperature (K)       '
        write(31,600) 't2max   ','maximum 2m air temperature (K)       '
        write(31,600) 't2min   ','minimum 2m air temperature (K)       '
        write(31,600) 'w10max  ','maximum 10m wind speed (m/s)         '
        write(31,600) 'ps_min  ','minimum surface pressure (hPa)       '
               write(31,700)
               close(31)
            endif
         enddo
            
         deallocate(I1UR)
         deallocate(J1UR)
         deallocate(I1UL)
         deallocate(J1UL)
         deallocate(I1DR)
         deallocate(J1DR)
         deallocate(I1DL)
         deallocate(J1DL)
         deallocate(D1XT)
         deallocate(D1Xa)
         deallocate(D1Xb)
         deallocate(D1Xc)
         deallocate(D1Xd)
      endif
      deallocate(out)
      deallocate(olat)
      deallocate(olon)

      deallocate(sigma)
      deallocate(o)
      deallocate(xlat)
      deallocate(xlon)

  10  format('dset ^',A17)
  11  format('dset ^',A15)
  12  format('dset ^',A10)
  20  format('title RegCM domain information')
  30  format('options big_endian')
  40  format('options little_endian')
  50  format('undef -9999.')
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
 210  format('ydef ',I3,' linear ',f9.4,' ',f9.4)
 300  format('zdef 1',' levels ',f7.2)
 400  format('tdef ',I4,' linear ',I2,'z',A2,A3,I4,' ',I2,'hr')
 401  format('tdef ',I4,' linear ',A2,A3,I4,' ','1dy')
 402  format('tdef ',I4,' linear ',A2,A3,I4,' ','1mo')
 500  format('vars ',I2)
 600  format(A8,'0 99 ',A36)           
 700  format('endvars')
      return
      end

      subroutine regrid_RAD(iy,jx,kz,ibyte,idate0,idate1,idate2 &
                           ,Path_Output)
      implicit none
      integer iy,jx,kz,ibyte,idate0,idate1,idate2
      character*128 Path_Output
      integer iiy,jjx,kkz
      integer mdate0,ibltyp,icup,ipptls,iboudy
      real*4  dxsp,ptsp,clat,clon,plat,plon
      real*4  dto,dtb,dtr,dtc
      integer iotyp
      character*6 iproj
      real*4, allocatable ::  sigma(:)
      character*4 :: chy(1941:2100)
      data chy/'1941','1942','1943','1944','1945','1946','1947','1948',&
        '1949','1950','1951','1952','1953','1954','1955','1956','1957',&
        '1958','1959','1960','1961','1962','1963','1964','1965','1966',&
        '1967','1968','1969','1970','1971','1972','1973','1974','1975',&
        '1976','1977','1978','1979','1980','1981','1982','1983','1984',&
        '1985','1986','1987','1988','1989','1990','1991','1992','1993',&
        '1994','1995','1996','1997','1998','1999','2000','2001','2002',&
        '2003','2004','2005','2006','2007','2008','2009','2010','2011',&
        '2012','2013','2014','2015','2016','2017','2018','2019','2020',&
        '2021','2022','2023','2024','2025','2026','2027','2028','2029',&
        '2030','2031','2032','2033','2034','2035','2036','2037','2038',&
        '2039','2040','2041','2042','2043','2044','2045','2046','2047',&
        '2048','2049','2050','2051','2052','2053','2054','2055','2056',&
        '2057','2058','2059','2060','2061','2062','2063','2064','2065',&
        '2066','2067','2068','2069','2070','2071','2072','2073','2074',&
        '2075','2076','2077','2078','2079','2080','2081','2082','2083',&
        '2084','2085','2086','2087','2088','2089','2090','2091','2092',&
        '2093','2094','2095','2096','2097','2098','2099','2100'/
      character*2 chm(12)
      data chm/'01','02','03','04','05','06','07','08','09','10', &
               '11','12'/
      character*3 chmc(12)
      data chmc/'jan','feb','mar','apr','may','jun'  &
               ,'jul','aug','sep','oct','nov','dec'/
      character*14 filein
      character*17 fileout
      integer ntype,nfile,nyear,month,n_slice,mrec,nrec,nnn

      integer i,j,k,l,m,n
      integer igrads,ibigend

      real*4  size_2d(4),ddeg
      COMMON /WINDOW/ size_2d,ddeg

      real*4, allocatable ::  o(:,:,:),xlat(:,:),xlon(:,:)
      real*4  xmaxlat,xminlat,xmaxlon,xminlon

      integer mlat,nlat,mlon,nlon
      integer n_month
      logical there

      integer MUR,NUR,MUL,NUL,MDR,NDR,MDL,NDL
      real*4  DISTa,DISTb,DISTc,DISTd,AAA,DIST
      real*4  PI
      real*4, allocatable :: out(:,:)
      real*4, allocatable :: olat(:),olon(:)

      integer, allocatable :: I1UR(:,:), J1UR(:,:)
      integer, allocatable :: I1UL(:,:), J1UL(:,:)
      integer, allocatable :: I1DR(:,:), J1DR(:,:)
      integer, allocatable :: I1DL(:,:), J1DL(:,:)
      real*4, allocatable :: D1XT(:,:), D1Xa(:,:), D1Xb(:,:)
      real*4, allocatable ::            D1Xc(:,:), D1Xd(:,:)

      igrads = 1
      ibigend= 1

      PI = atan(1.)*4.

      allocate(sigma(kz+1))
      allocate(o(jx-2,iy-2,kz*4+9))
      allocate(xlat(jx-2,iy-2))
      allocate(xlon(jx-2,iy-2))

      inquire(file=trim(Path_Output)//'OUT_HEAD',exist=there)
      if(.not.there) then
         write(*,*) trim(Path_Output)//'OUT_HEAD',' is not avaiable'
         stop
      endif
      open(10,file=trim(Path_Output)//'OUT_HEAD',form='unformatted' &
             ,recl=(jx-2)*(iy-2)*ibyte,access='direct')
      read(10,rec=1) mdate0,ibltyp,icup,ipptls,iboudy  &
                    ,iiy,jjx,kkz,(sigma(k),k=1,kz+1)   &
                    ,dxsp,ptsp,clat,clon,plat,plon          &
                    ,iproj,dto,dtb,dtr,dtc,iotyp
      if(iiy.ne.iy.or.jjx.ne.jx.or.kkz.ne.kz) then
         write(*,*) 'iy,jx,kz in parameter = ',iy,jx,kz
         write(*,*) 'iy,jx,kz in OUT_HEAD ',iiy,jjx,kkz
         write(*,*) 'They are not consistent'
         stop
      endif
      read(10,rec=6) xlat         ! xlat
      if(size_2d(3).gt.-9990.) then
         xminlat = size_2d(3)
         xmaxlat = size_2d(4)
      else
         xmaxlat = -1.e10
         xminlat =  1.e10
         do i=1,iy-2
         do j=1,jx-2
            if(xlat(j,i).gt.xmaxlat) xmaxlat=xlat(j,i)
            if(xlat(j,i).lt.xminlat) xminlat=xlat(j,i)
         enddo
         enddo
      endif
      read(10,rec=7) xlon         ! xlon
      if(size_2d(1).gt.-9990.) then
         xminlon=size_2d(1)
         xmaxlon=size_2d(2)
      else
         xmaxlon = -1.e10
         xminlon =  1.e10
         do i=1,iy-2
         do j=1,jx-2
            if(xlon(j,i).gt.xmaxlon) xmaxlon=xlon(j,i)
            if(xlon(j,i).lt.xminlon) xminlon=xlon(j,i)
         enddo
         enddo
      endif
      close(10)
      
      mlat = xmaxlat/ddeg - 1
      nlat = xminlat/ddeg + 1
      mlon = xmaxlon/ddeg - 1
      nlon = xminlon/ddeg + 2

      write(*,*) xminlat,mlat-nlat+1,(nlat-1)*ddeg
      write(*,*) xminlon,mlon-nlon+1,(nlon-1)*ddeg

      allocate(out(nlon:mlon,nlat:mlat))
      allocate(olat(nlat:mlat))
      allocate(olon(nlon:mlon))
      do n=nlat,mlat
         olat(n) = real(n-1)*ddeg
      enddo
      do m=nlon,mlon
         olon(m) = real(m-1)*ddeg
      enddo

      if(idate1.gt.1940010100) then        ! original file
         n_month = (idate2/1000000-idate1/1000000)*12 &
                 + (mod(idate2/10000,100)-mod(idate1/10000,100))   
         if(mod(idate2,10000).gt.0100) n_month = n_month+1
         ntype = 0
      else if(idate1.gt.19400101) then     ! daily mean file
         n_month = (idate2/10000-idate1/10000)*12 &
                 + (mod(idate2/100,100)-mod(idate1/100,100))   
         if(mod(idate2,100).gt.01) n_month = n_month+1
         ntype = 1
      else if(idate1.gt.194001) then       ! monthly mean file
         n_month = (idate2/100-idate1/100)*12 &
                 + (mod(idate2,100)-mod(idate1,100))+1
         ntype = 2
      endif

      if(iproj.eq.'NORMER') then
         do nfile=1,n_month
            if(ntype.eq.0) then
               nyear = idate1/1000000
               month = idate1/10000-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(nfile.eq.1.or.month.eq.1) then
                  write(*,*)'regrid: NORMER RAD Orig.',nyear,month
               else
                  write(*,*)'                        ',nyear,month
               endif

               if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
                  month.eq.7.or.month.eq.8.or.month.eq.10.or. &
                  month.eq.12) then
                  n_slice = 31*24/dtr
               else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                       month.eq.11) then
                  n_slice = 30*24/dtr
               else
                  n_slice = 28*24/dtr
                  if(mod(nyear,4).eq.0) n_slice = 29*24/dtr
                  if(mod(nyear,100).eq.0) n_slice = 28*24/dtr
                  if(mod(nyear,400).eq.0) n_slice = 29*24/dtr
               endif
               if(nfile.eq.n_month.and.mod(idate2/100-1,100).ne.0) &
               n_slice=min(n_slice,nint(mod(idate2/100-1,100)*24./dtr))
               if(nfile.eq.1.and.idate0.eq.idate1) n_slice=n_slice+1
               filein = 'RAD.'//chy(nyear)//chm(month)//'0100'
               fileout= 'RAD_LL.'//chy(nyear)//chm(month)//'0100'
               inquire(file=trim(Path_Output)//filein,exist=there)
               if(.not.there) then
                  write(*,*) trim(Path_Output)//filein,' is not avaiable'
                  stop
               endif
             open(10,file=trim(Path_Output)//filein,form='unformatted' &
                      ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
               mrec = 0
            open(20,file=trim(Path_Output)//fileout,form='unformatted' &
               ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
               nrec = 0
            else if(ntype.eq.1) then
               nyear = idate1/10000
               month = idate1/100-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(nfile.eq.1.or.month.eq.1) then
                  write(*,*)'regrid: NORMER RAD Daily',nyear,month
               else
                  write(*,*)'                        ',nyear,month
               endif

               if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
                  month.eq.7.or.month.eq.8.or.month.eq.10.or. &
                  month.eq.12) then
                  n_slice = 31
               else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                       month.eq.11) then
                  n_slice = 30
               else
                  n_slice = 28
                  if(mod(nyear,4).eq.0) n_slice = 29
                  if(mod(nyear,100).eq.0) n_slice = 28
                  if(mod(nyear,400).eq.0) n_slice = 29
               endif
               if(nfile.eq.n_month.and.mod(idate2-1,100).ne.0) &
               n_slice=min(n_slice,mod(idate2-1,100)+1)
               filein = 'RAD.'//chy(nyear)//chm(month)//'01'
               fileout= 'RAD_LL.'//chy(nyear)//chm(month)//'01'
               inquire(file=trim(Path_Output)//filein,exist=there)
               if(.not.there) then
                  write(*,*) trim(Path_Output)//filein,' is not avaiable'
                  stop
               endif
             open(10,file=trim(Path_Output)//filein,form='unformatted' &
                      ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
               mrec = 0
            open(20,file=trim(Path_Output)//fileout,form='unformatted' &
               ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
               nrec = 0
            else if(ntype.eq.2) then
               nyear = idate1/100
               month = idate1-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(month.eq.1.or.nfile.eq.1) then
                  write(*,*)'regrid: NORMER RAD Month',nyear,month
               else
                  write(*,*)'                        ',nyear,month
               endif

               n_slice = 1

               if(nfile.eq.1) then
                  filein = 'RAD.mon'
                  fileout= 'RAD_LL.mon'
                  inquire(file=trim(Path_Output)//filein(1:7)  &
                         ,exist=there)
                  if(.not.there) then
                     write(*,*) trim(Path_Output)//filein(1:7)  &
                               ,' is not avaiable'
                     stop
                  endif
                  open(10,file=trim(Path_Output)//filein(1:7)  &
                         ,form='unformatted' &
                         ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
                  nrec = 0
                  open(20,file=trim(Path_Output)//fileout(1:10)  &
                         ,form='unformatted' &
                ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
                  mrec = 0
               endif
            endif

            do nnn=1,n_slice
               do l=1,kz*4+9
                  mrec = mrec+1
                  read(10,rec=mrec) ((o(j,i,l),j=1,jx-2),i=1,iy-2)
               enddo
               do l=1,kz*4+9
                  do n=nlat,mlat
                  do m=nlon,mlon
                     out(m,n) = -9999.

                     do i=1,iy-3
                     do j=1,jx-3
                        if((xlat(j,i).lt.olat(n).and.      &
                            xlon(j,i).lt.olon(m))          &
                      .and.(xlat(j+1,i+1).ge.olat(n).and.  &
                            xlon(j+1,i+1).ge.olon(m))) then
                            out(m,n) = &
       ( o(j,i,l)*(xlon(j+1,i+1)-olon(m))*(xlat(j+1,i+1)-olat(n)) &
       + o(j+1,i,l)*(olon(m)-xlon(j,i+1))*(xlat(j+1,i+1)-olat(n)) &
       + o(j,i+1,l)*(xlon(j+1,i)-olon(m))*(olat(n)-xlat(j+1,i)) &
       + o(j+1,i+1,l)*(olon(m)-xlon(j,i))*(olat(n)-xlat(j,i)) ) &
       / ( (xlon(j+1,i+1)-xlon(j,i))*(xlat(j+1,i+1)-xlat(j,i)) )
                        endif
                     enddo
                     enddo
                  enddo
                  enddo
                  nrec = nrec+1
                  write(20,rec=nrec)((out(m,n),m=nlon,mlon),n=nlat,mlat)
               enddo
            enddo
            if(.not.ntype.eq.2) close(20)
            if(.not.ntype.eq.2) close(10)
            if(igrads.eq.1.and.(ntype.eq.0.or.ntype.eq.1.or. &
               (ntype.eq.2.and.nfile.eq.1))) then
               if(ntype.eq.0) then
                 inquire(file=trim(Path_Output)//fileout(1:17)//'.ctl' &
                        ,exist=there)
                  if(there) then
                     open(31,file= &
                          trim(Path_Output)//fileout(1:17)//'.ctl' &
                         ,form='formatted',status='replace')
                  else
                     open(31,file= &
                          trim(Path_Output)//fileout(1:17)//'.ctl' &
                         ,form='formatted',status='new')
                  endif
                  write(31,10) fileout(1:17)
               else if(ntype.eq.1) then
                 inquire(file=trim(Path_Output)//fileout(1:15)//'.ctl' &
                        ,exist=there)
                  if(there) then
                     open(31,file= &
                          trim(Path_Output)//fileout(1:15)//'.ctl' &
                         ,form='formatted',status='replace')
                  else
                     open(31,file= &
                          trim(Path_Output)//fileout(1:15)//'.ctl' &
                         ,form='formatted',status='new')
                  endif
                  write(31,11) fileout(1:15)
               else if(ntype.eq.2) then
                 inquire(file=trim(Path_Output)//fileout(1:10)//'.ctl' &
                        ,exist=there)
                  if(there) then
                     open(31,file= &
                          trim(Path_Output)//fileout(1:10)//'.ctl' &
                         ,form='formatted',status='replace')
                  else
                     open(31,file= &
                          trim(Path_Output)//fileout(1:10)//'.ctl' &
                         ,form='formatted',status='new')
                  endif
                  write(31,12) fileout(1:10)
               endif
               write(31,20)
               if(ibigend.eq.1) then
                  write(31,30)
               else
                  write(31,40)
               endif
               write(31,50)
               write(31,200) mlon-nlon+1,(nlon-1)*ddeg,ddeg
               write(31,210) mlat-nlat+1,(nlat-1)*ddeg,ddeg
               write(31,300) kz, &
          ((1013.25-ptsp*10.)*(sigma(k)+sigma(k+1))*0.5+ptsp*10.,k=1,kz)
               if(ntype.eq.0) then
               if(nfile.eq.1.and.idate0.eq.idate1) then
               write(31,400) n_slice,0,'01',chmc(month),nyear,nint(dtr)
               else
       write(31,400) n_slice,nint(dtr),'01',chmc(month),nyear,nint(dtr)
               endif
               else if(ntype.eq.1) then
                  write(31,401) n_slice,'01',chmc(month),nyear
               else if(ntype.eq.2) then
                  write(31,402) n_month,'16',chmc(month),nyear
               endif
               write(31,500) 13
        write(31,650)'cld   ',kz,'cloud fractional cover               '
        write(31,650)'clwp  ',kz,'cloud liquid water path              '
        write(31,650)'qrs   ',kz,'solar heating rate                   '
        write(31,650)'qrl   ',kz,'longwave cooling rate                '
        write(31,600)'frsa  ',   'surface absorbed solar flux          '
        write(31,600)'frla  ',   'longwave cooling of surface          '
        write(31,600)'clrst ',   'clearsky total column abs solar flux '
        write(31,600)'clrss ',   'clearsky surface absorbed solar flux '
        write(31,600)'clrlt ',   'clearsky net upward LW flux at TOA   '
        write(31,600)'clrls ',   'clearsky LW cooling at surface (W/m2)'
        write(31,600)'solin ',   'instantaneous incident solar (W/m2)  '
        write(31,600)'sabtp ',   'total column absorbed solar flux W/m2'
        write(31,600)'firtp ',   'net upward LW flux at TOA (W/m2)     '
               write(31,700)
               close(31)
            endif
         enddo

      else if(iproj.eq.'LAMCON'.or.iproj.eq.'ROTMER') then

         allocate(I1UR(nlon:mlon,nlat:mlat))
         allocate(J1UR(nlon:mlon,nlat:mlat))
         allocate(I1UL(nlon:mlon,nlat:mlat))
         allocate(J1UL(nlon:mlon,nlat:mlat))
         allocate(I1DR(nlon:mlon,nlat:mlat))
         allocate(J1DR(nlon:mlon,nlat:mlat))
         allocate(I1DL(nlon:mlon,nlat:mlat))
         allocate(J1DL(nlon:mlon,nlat:mlat))
         allocate(D1XT(nlon:mlon,nlat:mlat))
         allocate(D1Xa(nlon:mlon,nlat:mlat))
         allocate(D1Xb(nlon:mlon,nlat:mlat))
         allocate(D1Xc(nlon:mlon,nlat:mlat))
         allocate(D1Xd(nlon:mlon,nlat:mlat))

         do n=nlat,mlat
         do m=nlon,mlon

            MUR=1000
            NUR=1000
            MUL=1000
            NUL=1000
            MDR=1000
            NDR=1000
            MDL=1000
            NDL=1000
            DISTa=1.E8
            DISTb=1.E8
            DISTc=1.E8
            DISTd=1.E8

            do i=1,iy-2
            do j=1,jx-2
         IF((xlon(j,i).ge.olon(m).and.xlon(j,i)-olon(m).lt.10.) .and. &
            (xlat(j,i).ge.olat(n).and.xlat(j,i)-olat(n).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTa.gt.AAA) then
                DISTa = AAA
                MUR = j
                NUR = i
             endif
         ENDIF
         IF((xlon(j,i).lt.olon(m).and.olon(m)-xlon(j,i).lt.10.) .and. &
            (xlat(j,i).ge.olat(n).and.xlat(j,i)-olat(n).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTb.gt.AAA) then
                DISTb = AAA
                MUL = j
                NUL = i
             endif
         ENDIF
         IF((xlon(j,i).ge.olon(m).and.xlon(j,i)-olon(m).lt.10.) .and. &
            (xlat(j,i).lt.olat(n).and.olat(n)-xlat(j,i).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTc.gt.AAA) then
                DISTc = AAA
                MDR = j
                NDR = i
             endif
         ENDIF
         IF((xlon(j,i).lt.olon(m).and.olon(m)-xlon(j,i).lt.10.) .and. &
            (xlat(j,i).lt.olat(n).and.olat(n)-xlat(j,i).lt.10.)) then
             AAA = ((xlon(j,i)-olon(m)) &
                 *cos((xlat(j,i)+olat(n))/360.*pi))**2 &
                 +(xlat(j,i)-olat(n))**2
             if(DISTd.gt.AAA) then
                DISTd = AAA
                MDL = j
                NDL = i
             endif
         ENDIF
            enddo
            enddo
            
            DIST=amin1(DISTa,DISTb,DISTc,DISTd)
            I1UR(m,n) = MUR
            J1UR(m,n) = NUR
            I1UL(m,n) = MUL
            J1UL(m,n) = NUL
            I1DR(m,n) = MDR
            J1DR(m,n) = NDR
            I1DL(m,n) = MDL
            J1DL(m,n) = NDL
            D1XT(m,n) = DIST
            D1Xa(m,n) = DISTa
            D1Xb(m,n) = DISTb
            D1Xc(m,n) = DISTc
            D1Xd(m,n) = DISTd
         enddo
         enddo

         do nfile=1,n_month
            if(ntype.eq.0) then
               nyear = idate1/1000000
               month = idate1/10000-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(iproj.eq.'LAMCON') then
                  if(nfile.eq.1.or.month.eq.1) then
                     write(*,*)'regrid: LAMCON RAD Orig.',nyear,month
                  else
                     write(*,*)'                        ',nyear,month
                  endif
               else if(iproj.eq.'ROTMER') then
                  if(nfile.eq.1.or.month.eq.1) then
                     write(*,*)'regrid: ROTMER RAD Orig.',nyear,month
                  else
                     write(*,*)'                        ',nyear,month
                  endif
               endif

               if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
                  month.eq.7.or.month.eq.8.or.month.eq.10.or. &
                  month.eq.12) then
                  n_slice = 31*24/dtr
               else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                       month.eq.11) then
                  n_slice = 30*24/dtr
               else
                  n_slice = 28*24/dtr
                  if(mod(nyear,4).eq.0) n_slice = 29*24/dtr
                  if(mod(nyear,100).eq.0) n_slice = 28*24/dtr
                  if(mod(nyear,400).eq.0) n_slice = 29*24/dtr
               endif
               if(nfile.eq.n_month.and.mod(idate2/100-1,100).ne.0) &
               n_slice=min(n_slice,nint(mod(idate2/100-1,100)*24./dtr))
               if(nfile.eq.1.and.idate0.eq.idate1) n_slice=n_slice+1
               filein = 'RAD.'//chy(nyear)//chm(month)//'0100'
               fileout= 'RAD_LL.'//chy(nyear)//chm(month)//'0100'
               inquire(file=trim(Path_Output)//filein,exist=there)
               if(.not.there) then
                  write(*,*) trim(Path_Output)//filein,' is not avaiable'
                  stop
               endif
             open(10,file=trim(Path_Output)//filein,form='unformatted' &
                      ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
               mrec = 0
            open(20,file=trim(Path_Output)//fileout,form='unformatted' &
               ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
               nrec = 0
            else if(ntype.eq.1) then
               nyear = idate1/10000
               month = idate1/100-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(iproj.eq.'LAMCON') then
                  if(nfile.eq.1.or.month.eq.1) then
                     write(*,*)'regrid: LAMCON RAD Daily',nyear,month
                  else
                     write(*,*)'                        ',nyear,month
                  endif
               else if(iproj.eq.'ROTMER') then
                  if(nfile.eq.1.or.month.eq.1) then
                     write(*,*)'regrid: ROTMER RAD Daily',nyear,month
                  else
                     write(*,*)'                        ',nyear,month
                  endif
               endif

               if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
                  month.eq.7.or.month.eq.8.or.month.eq.10.or. &
                  month.eq.12) then
                  n_slice = 31
               else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                       month.eq.11) then
                  n_slice = 30
               else
                  n_slice = 28
                  if(mod(nyear,4).eq.0) n_slice = 29
                  if(mod(nyear,100).eq.0) n_slice = 28
                  if(mod(nyear,400).eq.0) n_slice = 29
               endif
               if(nfile.eq.n_month.and.mod(idate2-1,100).ne.0) &
               n_slice=min(n_slice,mod(idate2-1,100)+1)
               filein = 'RAD.'//chy(nyear)//chm(month)//'01'
               fileout= 'RAD_LL.'//chy(nyear)//chm(month)//'01'
               inquire(file=trim(Path_Output)//filein,exist=there)
               if(.not.there) then
                  write(*,*) trim(Path_Output)//filein,' is not avaiable'
                  stop
               endif
             open(10,file=trim(Path_Output)//filein,form='unformatted' &
                      ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
               mrec = 0
            open(20,file=trim(Path_Output)//fileout,form='unformatted' &
               ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
               nrec = 0
            else if(ntype.eq.2) then
               nyear = idate1/100
               month = idate1-nyear*100 +nfile-1
               nyear = nyear + (month-1)/12
               month = mod(month,12)
               if(month.eq.0) month = 12

               if(iproj.eq.'LAMCON') then
                  if(month.eq.1.or.nfile.eq.1) then
                     write(*,*)'regrid: LAMCON RAD Month',nyear,month
                  else
                     write(*,*)'                        ',nyear,month
                  endif
               else if(iproj.eq.'ROTMER') then
                  if(month.eq.1.or.nfile.eq.1) then
                     write(*,*)'regrid: ROTMER RAD Month',nyear,month
                  else
                     write(*,*)'                        ',nyear,month
                  endif
               endif

               n_slice = 1

               if(nfile.eq.1) then
                  filein = 'RAD.mon'
                  fileout= 'RAD_LL.mon'
                  inquire(file=trim(Path_Output)//filein(1:7)  &
                         ,exist=there)
                  if(.not.there) then
                     write(*,*) trim(Path_Output)//filein(1:7)  &
                               ,' is not avaiable'
                     stop
                  endif
                  open(10,file=trim(Path_Output)//filein(1:7)  &
                         ,form='unformatted' &
                         ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
                  nrec = 0
                  open(20,file=trim(Path_Output)//fileout(1:10)  &
                         ,form='unformatted' &
                ,recl=(mlon-nlon+1)*(mlat-nlat+1)*ibyte,access='direct')
                  mrec = 0
               endif
            endif

            do nnn=1,n_slice
               do l=1,kz*4+9
                  mrec = mrec+1
                  read(10,rec=mrec) ((o(j,i,l),j=1,jx-2),i=1,iy-2)
               enddo
               do l=1,kz*4+9
                  do n=nlat,mlat
                  do m=nlon,mlon
                     out(m,n) = -9999.

            if(I1UR(m,n).lt.999.and.J1UR(m,n).lt.999.and. &
               I1UL(m,n).lt.999.and.J1UL(m,n).lt.999.and. &
               I1DR(m,n).lt.999.and.J1DR(m,n).lt.999.and. &
               I1DL(m,n).lt.999.and.J1DL(m,n).lt.999) then
               if(D1XT(m,n).gt.0.0001) then
      out(m,n) = ( o(I1UR(m,n),J1UR(m,n),l)/D1Xa(m,n) &
                  +o(I1UL(m,n),J1UL(m,n),l)/D1Xb(m,n) &
                  +o(I1DR(m,n),J1DR(m,n),l)/D1Xc(m,n) &
                  +o(I1DL(m,n),J1DL(m,n),l)/D1Xd(m,n) )  &
          /( 1./D1Xa(m,n)+1./D1Xb(m,n)+1./D1Xc(m,n)+1./D1Xd(m,n) )
               else
                  if(D1Xa(m,n).eq.D1XT(m,n)) then
                     out(m,n) = o(I1UR(m,n),J1UR(m,n),l)
                  else if(D1Xb(m,n).eq.D1XT(m,n)) then
                     out(m,n) = o(I1UL(m,n),J1UL(m,n),l)
                  else if(D1Xc(m,n).eq.D1XT(m,n)) then
                     out(m,n) = o(I1DR(m,n),J1DR(m,n),l)
                  else if(D1Xd(m,n).eq.D1XT(m,n)) then
                     out(m,n) = o(I1DL(m,n),J1DL(m,n),l)
                  endif
               endif
            endif
                  enddo
                  enddo
                  nrec=nrec+1
                  write(20,rec=nrec)((out(m,n),m=nlon,mlon),n=nlat,mlat)
               enddo
            enddo
            if(.not.ntype.eq.2) close(20)
            if(.not.ntype.eq.2) close(10)
            if(igrads.eq.1.and.(ntype.eq.0.or.ntype.eq.1.or. &
               (ntype.eq.2.and.nfile.eq.1))) then
               if(ntype.eq.0) then
                 inquire(file=trim(Path_Output)//fileout(1:17)//'.ctl' &
                        ,exist=there)
                  if(there) then
                     open(31,file= &
                          trim(Path_Output)//fileout(1:17)//'.ctl' &
                         ,form='formatted',status='replace')
                  else
                     open(31,file= &
                          trim(Path_Output)//fileout(1:17)//'.ctl' &
                         ,form='formatted',status='new')
                  endif
                  write(31,10) fileout(1:17)
               else if(ntype.eq.1) then
                 inquire(file=trim(Path_Output)//fileout(1:15)//'.ctl' &
                        ,exist=there)
                  if(there) then
                     open(31,file= &
                          trim(Path_Output)//fileout(1:15)//'.ctl' &
                         ,form='formatted',status='replace')
                  else
                     open(31,file= &
                          trim(Path_Output)//fileout(1:15)//'.ctl' &
                         ,form='formatted',status='new')
                  endif
                  write(31,11) fileout(1:15)
               else if(ntype.eq.2) then
                 inquire(file=trim(Path_Output)//fileout(1:10)//'.ctl' &
                        ,exist=there)
                  if(there) then
                     open(31,file= &
                          trim(Path_Output)//fileout(1:10)//'.ctl' &
                         ,form='formatted',status='replace')
                  else
                     open(31,file= &
                          trim(Path_Output)//fileout(1:10)//'.ctl' &
                         ,form='formatted',status='new')
                  endif
                  write(31,12) fileout(1:10)
               endif
               write(31,20)
               if(ibigend.eq.1) then
                  write(31,30)
               else
                  write(31,40)
               endif
               write(31,50)
               write(31,200) mlon-nlon+1,(nlon-1)*ddeg,ddeg
               write(31,210) mlat-nlat+1,(nlat-1)*ddeg,ddeg
               write(31,300) kz, &
          ((1013.25-ptsp*10.)*(sigma(k)+sigma(k+1))*0.5+ptsp*10.,k=1,kz)
               if(ntype.eq.0) then
               if(nfile.eq.1.and.idate0.eq.idate1) then
               write(31,400) n_slice,0,'01',chmc(month),nyear,nint(dtr)
               else
       write(31,400) n_slice,nint(dtr),'01',chmc(month),nyear,nint(dtr)
               endif
               else if(ntype.eq.1) then
                  write(31,401) n_slice,'01',chmc(month),nyear
               else if(ntype.eq.2) then
                  write(31,402) n_month,'16',chmc(month),nyear
               endif
               write(31,500) 13
        write(31,650)'cld   ',kz,'cloud fractional cover               '
        write(31,650)'clwp  ',kz,'cloud liquid water path              '
        write(31,650)'qrs   ',kz,'solar heating rate                   '
        write(31,650)'qrl   ',kz,'longwave cooling rate                '
        write(31,600)'frsa  ',   'surface absorbed solar flux          '
        write(31,600)'frla  ',   'longwave cooling of surface          '
        write(31,600)'clrst ',   'clearsky total column abs solar flux '
        write(31,600)'clrss ',   'clearsky surface absorbed solar flux '
        write(31,600)'clrlt ',   'clearsky net upward LW flux at TOA   '
        write(31,600)'clrls ',   'clearsky LW cooling at surface (W/m2)'
        write(31,600)'solin ',   'instantaneous incident solar (W/m2)  '
        write(31,600)'sabtp ',   'total column absorbed solar flux W/m2'
        write(31,600)'firtp ',   'net upward LW flux at TOA (W/m2)     '
               write(31,700)
               close(31)
            endif
         enddo
            
         deallocate(I1UR)
         deallocate(J1UR)
         deallocate(I1UL)
         deallocate(J1UL)
         deallocate(I1DR)
         deallocate(J1DR)
         deallocate(I1DL)
         deallocate(J1DL)
         deallocate(D1XT)
         deallocate(D1Xa)
         deallocate(D1Xb)
         deallocate(D1Xc)
         deallocate(D1Xd)
      endif
      deallocate(out)
      deallocate(olat)
      deallocate(olon)

      deallocate(sigma)
      deallocate(o)
      deallocate(xlat)
      deallocate(xlon)

  10  format('dset ^',A17)
  11  format('dset ^',A15)
  12  format('dset ^',A10)
  20  format('title RegCM domain information')
  30  format('options big_endian')
  40  format('options little_endian')
  50  format('undef -9999.')
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
 210  format('ydef ',I3,' linear ',f9.4,' ',f9.4)
 300  format('zdef ',I2,' levels ',30f8.2)
 400  format('tdef ',I4,' linear ',I2,'z',A2,A3,I4,' ',I2,'hr')
 401  format('tdef ',I4,' linear ',A2,A3,I4,' ','1dy')
 402  format('tdef ',I4,' linear ',A2,A3,I4,' ','1mo')
 500  format('vars ',I2)
 600  format(A6,'0 99 ',A36)           
 650  format(A6,I2,' 0 ',A36)
 700  format('endvars')
      return
      end
