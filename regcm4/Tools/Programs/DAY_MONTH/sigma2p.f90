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

! This code converts the 3D fields of T, Q, U and V from model sigma levels,
!  to several assigned presure levels for ICBC and ATM files. The geopotential 
!  height and relative humidity fields are outputed at presure levels as well.  
!  Written by Xunqiang Bi, ESP/ICTP
!
      program sigma2p
      implicit none
      integer iy,jx,kz,np,nsg,ntr,ibyte
      integer idate0,idate1,idate2
      real*4, save ::  plev(11)
      logical s2p_ICBC,s2p_ATM
      character*128 Path_Input,Path_Output
      character*20 DomainName

      namelist /shareparam/ iy,jx,kz,np,nsg,ntr,ibyte,Path_Input &
                           ,DomainName,Path_Output
      namelist /dateparam/ idate0,idate1,idate2
      namelist /sigma2p_param/ plev,s2p_ICBC,s2p_ATM

      read(*,shareparam) 
      if(np.ne.11) then
         write(*,*) 'Number of pressure levels does not equal to 11'
         write(*,*) 'Please reset np, or reset plev'
      endif

      read(*,dateparam) 
      read(*,sigma2p_param) 

      if(s2p_ICBC) call sigma2p_ICBC(iy,jx,kz,np,plev,ibyte &
                        ,idate0,idate1,idate2,Path_Input,DomainName)

      if(s2p_ATM) call sigma2p_ATM(iy,jx,kz,np,plev,ibyte &
                                  ,idate0,idate1,idate2,Path_Output)

      stop
      end

      subroutine sigma2p_ICBC(iy,jx,kz,np,plev,ibyte &
                 ,idate0,idate1,idate2,Path_Input,DomainName)
      implicit none
      integer iy,jx,kz,np,ibyte,idate0,idate1,idate2
      real*4  plev(np)
      character*128 Path_Input
      character*20 DomainName
      integer iiy,jjx,kkz
      real*4  dsinm,clat,clon,plat,plon,GRDFAC
      character*6 iproj
      real*4, allocatable,save ::  sigma(:),sig(:)
      integer igrads,ibigend
      real*4  truelatL,truelatH
      character*4 :: chy
      character*2 cday(31)
      data cday/'01','02','03','04','05','06','07','08','09','10', &
                '11','12','13','14','15','16','17','18','19','20', &
                '21','22','23','24','25','26','27','28','29','30','31'/
      integer nday
      character*2 chm(12)
      data chm/'01','02','03','04','05','06','07','08','09','10', &
               '11','12'/
      character*3 chmc(12)
      data chmc/'jan','feb','mar','apr','may','jun'  &
               ,'jul','aug','sep','oct','nov','dec'/
      character*14 filein
      character*16 fileout
      integer ntype,nfile,nyear,month,n_slice,mrec,nrec,nnn

      integer i,j,k,l,m,n

      real*4, allocatable,save ::  xlat(:,:),xlon(:,:)

      real*4, allocatable ::  fin(:,:)

      real*4, allocatable,save :: u(:,:,:),v(:,:,:),t(:,:,:),q(:,:,:)
      real*4, allocatable,save :: ps(:,:),tgb(:,:),ht(:,:),h(:,:,:)
      real*4, allocatable,save :: slp(:,:)

      real*4, allocatable,save :: up(:,:,:),vp(:,:,:),tp(:,:,:)
      real*4, allocatable,save :: qp(:,:,:),hp(:,:,:)

      real*4, save :: ptop
      real*4, parameter, save :: rgas = 287.04
      real*4, parameter, save :: grav = 9.80616
      real*4, parameter, save :: bltop = 0.96
      real*4, parameter, save :: tlapse = -6.5E-3

      integer n_month
      logical there
      real*4  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real*4  centerj,centeri
      integer ny,nx

      integer IDATE

      allocate(sigma(kz+1))
      allocate(sig(kz))
      allocate(fin(jx,iy))

      allocate(xlat(jx-2,iy-2))
      allocate(xlon(jx-2,iy-2))
      allocate(u(jx-2,iy-2,kz))
      allocate(v(jx-2,iy-2,kz))
      allocate(t(jx-2,iy-2,kz))
      allocate(q(jx-2,iy-2,kz))
      allocate(h(jx-2,iy-2,kz))
      allocate(ps(jx-2,iy-2))
      allocate(tgb(jx-2,iy-2))
      allocate(ht(jx-2,iy-2))

      allocate(slp(jx-2,iy-2))
      allocate(up(jx-2,iy-2,np))
      allocate(vp(jx-2,iy-2,np))
      allocate(tp(jx-2,iy-2,np))
      allocate(qp(jx-2,iy-2,np))
      allocate(hp(jx-2,iy-2,np))

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

      ptop = ptop*10.
      do k=1,kz
         sig(k) = (sigma(k)+sigma(k+1))*0.5
      enddo

      read(10,rec=2) fin          ! ht
      do i=1,iy-2
      do j=1,jx-1
         ht(j,i) = fin(j+1,i+1)
      enddo
      enddo
      read(10,rec=5) fin          ! xlat
      do i=1,iy-2
      do j=1,jx-1
         xlat(j,i) = fin(j+1,i+1)
      enddo
      enddo
      read(10,rec=6) fin          ! xlon
      do i=1,iy-2
      do j=1,jx-1
         xlon(j,i) = fin(j+1,i+1)
      enddo
      enddo
      close(10)
      
      if(idate1.ge.1000010100) then        ! original file
         n_month = (idate2/1000000-idate1/1000000)*12 &
                 + (mod(idate2/10000,100)-mod(idate1/10000,100))   
         if(mod(idate2,10000).gt.0100) n_month = n_month+1
         ntype = 0
      else if(idate1.ge.10000101) then     ! daily mean file
         n_month = (idate2/10000-idate1/10000)*12 &
                 + (mod(idate2/100,100)-mod(idate1/100,100))   
         if(mod(idate2,100).gt.01) n_month = n_month+1
         ntype = 1
      else if(idate1.ge.100001) then       ! monthly mean file
         n_month = (idate2/100-idate1/100)*12 &
                 + (mod(idate2,100)-mod(idate1,100))+1 
         ntype = 2
      endif

      do nfile=1,n_month
         if(ntype.eq.0) then
            nyear = idate1/1000000
            month = idate1/10000-nyear*100 +nfile-1
            nyear = nyear + (month-1)/12
            month = mod(month,12)
            if(month.eq.0) month = 12

            if(nfile.eq.1.or.month.eq.1) then
               write(*,*)'sigma2p: ICBC Orig.',nyear,month
            else
               write(*,*)'                   ',nyear,month
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
            write(chy,199) nyear
            filein = 'ICBC'//chy//chm(month)//'0100'
            fileout= 'ICBC_P'//chy//chm(month)//'0100'

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
           ,form='unformatted',recl=(iy-2)*(jx-2)*ibyte,access='direct')
            nrec = 0
         else if(ntype.eq.1) then
            nyear = idate1/10000
            month = idate1/100-nyear*100 +nfile-1
            nyear = nyear + (month-1)/12
            month = mod(month,12)
            if(month.eq.0) month = 12

            if(nfile.eq.1.or.month.eq.1) then
               write(*,*)'sigma2p: ICBC Daily',nyear,month
            else
               write(*,*)'                   ',nyear,month
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
            write(chy,199) nyear
            filein = 'ICBC'//chy//chm(month)//'01'
            fileout= 'ICBC_P'//chy//chm(month)//'01'

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
           ,form='unformatted',recl=(iy-2)*(jx-2)*ibyte,access='direct')
            mrec = 0
            open(20,file= &
                trim(Path_Input)//trim(DomainName)//'_'//fileout(1:14) &
           ,form='unformatted',recl=(iy-2)*(jx-2)*ibyte,access='direct')
            nrec = 0
         else if(ntype.eq.2) then
            nyear = idate1/100
            month = idate1-nyear*100 +nfile-1
            nyear = nyear + (month-1)/12
            month = mod(month,12)
            if(month.eq.0) month = 12

            if(nfile.eq.1.or.month.eq.1) then
               write(*,*)'sigma2p: ICBC Month',nyear,month
            else
               write(*,*)'                   ',nyear,month
            endif

            n_slice = 1

            if(nfile.eq.1) then
               filein = 'ICBC'
               fileout= 'ICBC_P'

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
         trim(Path_Input)//trim(DomainName)//'_'//fileout(1:6)//'.mon' &
                    ,form='unformatted' &
                    ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
               nrec = 0
            endif
         endif

         do nnn=1,n_slice
            if(ntype.eq.0) then
               mrec = mrec+1
               read(10,rec=mrec) IDATE
!              write(*,*) 'IDATE = ',IDATE
               do l=kz,1,-1
                  mrec = mrec+1
                  read(10,rec=mrec) fin
                  do i=1,iy-2
                  do j=1,jx-2
                     u(j,i,l) = fin(j+1,i+1)
                  enddo
                  enddo
               enddo
               do l=kz,1,-1
                  mrec = mrec+1
                  read(10,rec=mrec) fin
                  do i=1,iy-2
                  do j=1,jx-2
                     v(j,i,l) = fin(j+1,i+1)
                  enddo
                  enddo
               enddo
               do l=kz,1,-1
                  mrec = mrec+1
                  read(10,rec=mrec) fin
                  do i=1,iy-2
                  do j=1,jx-2
                     t(j,i,l) = fin(j+1,i+1)
                  enddo
                  enddo
               enddo
               do l=kz,1,-1
                  mrec = mrec+1
                  read(10,rec=mrec) fin
                  do i=1,iy-2
                  do j=1,jx-2
                     q(j,i,l) = fin(j+1,i+1)
                  enddo
                  enddo
               enddo
               mrec = mrec+1
               read(10,rec=mrec) fin
               do i=1,iy-2
               do j=1,jx-2
                  ps(j,i) = fin(j+1,i+1)*10.
               enddo
               enddo
               mrec = mrec+1
               read(10,rec=mrec) fin
               do i=1,iy-2
               do j=1,jx-2
                  tgb(j,i) = fin(j+1,i+1)
                  tgb(j,i) = t(j,i,kz)
               enddo
               enddo
            else
               do l=kz,1,-1
                  mrec = mrec+1
                  read(10,rec=mrec) ((u(j,i,l),j=1,jx-2),i=1,iy-2)
               enddo
               do l=kz,1,-1
                  mrec = mrec+1
                  read(10,rec=mrec) ((v(j,i,l),j=1,jx-2),i=1,iy-2)
               enddo
               do l=kz,1,-1
                  mrec = mrec+1
                  read(10,rec=mrec) ((t(j,i,l),j=1,jx-2),i=1,iy-2)
               enddo
               do l=kz,1,-1
                  mrec = mrec+1
                  read(10,rec=mrec) ((q(j,i,l),j=1,jx-2),i=1,iy-2)
               enddo
               mrec = mrec+1
               read(10,rec=mrec) ((ps(j,i),j=1,jx-2),i=1,iy-2)
               do i=1,iy-2
               do j=1,jx-2
                  ps(j,i) = ps(j,i)*10.
               enddo
               enddo
               mrec = mrec+1
               read(10,rec=mrec) ((tgb(j,i),j=1,jx-2),i=1,iy-2)
               do i=1,iy-2
               do j=1,jx-2
                  tgb(j,i) = t(j,i,kz)
               enddo
               enddo
            endif

          ! to calculate Heights on sigma surfaces.
            call htsig(t,h,ps,ht,sig,jx-2,iy-2,kz,PTOP,RGAS,GRAV)

          ! to calculate Sea-Level Pressure using
          !  1. ERRICO's solution described in height
          !  2. a simple formulae
          !  3. MM5 method
            call slpres(h,t,ps,ht,tgb,slp,sig,jx-2,iy-2,kz  &
                       ,RGAS,GRAV,BLTOP,TLAPSE)

          ! to interpolate H,U,V,T,Q and QC
          !  1. For Heights
            call height(hp,h,t,ps,ht,sig,jx-2,iy-2,kz,plev,np &
                       ,PTOP,RGAS,GRAV,BLTOP,TLAPSE)
            do k=1,np
               nrec=nrec+1
               write(20,rec=nrec)((hp(j,i,k),j=1,jx-2),i=1,iy-2)
            enddo
          !  2. For Zonal and Meridional Winds
            call intlin(up,u,ps,sig,jx-2,iy-2,kz,plev,np,PTOP)
            call intlin(vp,v,ps,sig,jx-2,iy-2,kz,plev,np,PTOP)
          !  3. For Temperatures
            call intlog(tp,t,ps,sig,jx-2,iy-2,kz,plev,np &
                       ,PTOP,RGAS,GRAV,BLTOP,TLAPSE)
            do k=1,np
               nrec=nrec+1
               write(20,rec=nrec)((tp(j,i,k),j=1,jx-2),i=1,iy-2)
            enddo
            do k=1,np
               nrec=nrec+1
               write(20,rec=nrec)((up(j,i,k),j=1,jx-2),i=1,iy-2)
            enddo
            do k=1,np
               nrec=nrec+1
               write(20,rec=nrec)((vp(j,i,k),j=1,jx-2),i=1,iy-2)
            enddo
          !  4. For Moisture qv
            call humid1(t,q,ps,sig,jx-2,iy-2,kz,ptop)
            call intlin(qp,q,ps,sig,jx-2,iy-2,kz,plev,np,PTOP)
            do k=1,np
               nrec=nrec+1
               write(20,rec=nrec)((qp(j,i,k)*100.,j=1,jx-2),i=1,iy-2)
            enddo

            call humid2(tp,qp,plev,jx-2,iy-2,np)
            do k=1,np
               nrec=nrec+1
               write(20,rec=nrec)((qp(j,i,k),j=1,jx-2),i=1,iy-2)
            enddo
            nrec=nrec+1
            write(20,rec=nrec) ps
            nrec=nrec+1
            write(20,rec=nrec) slp
         enddo
         if(.not.ntype.eq.2) close(20)
         if(.not.ntype.eq.2) close(10)
         if(igrads.eq.1.and.(ntype.eq.0.or.ntype.eq.1.or. &
            (ntype.eq.2.and.nfile.eq.1))) then
            if(ntype.eq.0) then
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
               write(31,10) '^'//trim(DomainName)//'_'//fileout(1:16)
            else if(ntype.eq.1) then
               inquire(file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:14)//'.ctl' &
                      ,exist=there)
               if(there) then
                  open(31,file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:14)//'.ctl' &
                      ,form='formatted',status='replace')
               else
                  open(31,file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:14)//'.ctl' &
                      ,form='formatted',status='new')
               endif
               write(31,11) '^'//trim(DomainName)//'_'//fileout(1:14)
            else if(ntype.eq.2) then
               inquire(file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:6)//'.ctl' &
                         ,exist=there)
               if(there) then
                  open(31,file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:6)//'.ctl' &
                      ,form='formatted',status='replace')
               else
                  open(31,file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout(1:6)//'.ctl' &
                      ,form='formatted',status='new')
               endif
          write(31,12) '^'//trim(DomainName)//'_'//fileout(1:6)//'.mon'
            endif
            write(31,20)
            if(ibigend.eq.1) then
               write(31,30)
            else
               write(31,40)
            endif
            write(31,50)
            if(iproj.eq.'LAMCON'.or.iproj.eq.'ROTMER') then
               alatmin= 999999.
               alatmax=-999999.
               do j=1,jx-2
                  if(xlat(j,1).lt.alatmin) alatmin=xlat(j,1)
                  if(xlat(j,iy-2).gt.alatmax) alatmax=xlat(j,iy-2)
               enddo
               alonmin= 999999.
               alonmax=-999999.
               do i=1,iy-2
               do j=1,jx-2
                  if(clon.ge.0.0) then
                     if(xlon(j,i).ge.0.0) then
                        alonmin = amin1(alonmin,xlon(j,i))
                        alonmax = amax1(alonmax,xlon(j,i))
                     else
                        if(abs(clon-xlon(j,i)).lt. &
                           abs(clon-(xlon(j,i)+360.))) then
                           alonmin = amin1(alonmin,xlon(j,i))
                           alonmax = amax1(alonmax,xlon(j,i))
                        else
                           alonmin = amin1(alonmin,xlon(j,i)+360.)
                           alonmax = amax1(alonmax,xlon(j,i)+360.)
                        endif
                     endif
                  else
                     if(xlon(j,i).lt.0.0) then
                        alonmin = amin1(alonmin,xlon(j,i))
                        alonmax = amax1(alonmax,xlon(j,i))
                     else
                        if(abs(clon-xlon(j,i)).lt. &
                           abs(clon-(xlon(j,i)-360.))) then
                           alonmin = amin1(alonmin,xlon(j,i))
                           alonmax = amax1(alonmax,xlon(j,i))
                        else
                           alonmin = amin1(alonmin,xlon(j,i)-360.)
                           alonmax = amax1(alonmax,xlon(j,i)-360.)
                        endif
                     endif
                  endif
               enddo
               enddo
               rlatinc=dsinm*0.001/111./2.
               rloninc=dsinm*0.001/111./2.
               ny=2+nint(abs(alatmax-alatmin)/rlatinc)
               nx=1+nint(abs((alonmax-alonmin)/rloninc))
               centerj=(jx-2)/2.
               centeri=(iy-2)/2.
            endif
            if(iproj.eq.'LAMCON') then        ! Lambert projection
               write(31,100) jx-2,iy-2,clat,clon,centerj,centeri, &
                             truelatL,truelatH,clon,dsinm,dsinm
 100  format('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
               write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
               write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
            elseif(iproj.eq.'POLSTR') then    !
            elseif(iproj.eq.'NORMER') then
               write(31,200)  jx-2,xlon(1,1),xlon(2,1)-xlon(1,1)
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
               write(31,210) iy-2
 210  format('ydef ',I3,' levels')
               write(31,220) (xlat(1,i),i=1,iy-2)
 220  format(10f7.2)
            elseif(iproj.eq.'ROTMER') then
            if(nfile.eq.1.or.month.eq.1) then
               write(*,*) 'Note that rotated Mercartor (ROTMER)' &
                        ,' projections are not supported by GrADS.'
               write(*,*) '  Although not exact, the eta.u projection' &
                        ,' in GrADS is somewhat similar.'
            write(*,*) ' FERRET, however, does support this projection.'
            endif
               write(31,230) jx-2,iy-2,plon,plat,dsinm/111000. &
                                           ,dsinm/111000.*.95238
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
               write(31,110) nx+2,alonmin-rloninc,rloninc
               write(31,120) ny+2,alatmin-rlatinc,rlatinc
            else
               write(*,*) 'Are you sure your map projection is right ?'
               stop
            endif
            write(31,300) np, (plev(k),k=1,np)
            if(ntype.eq.0) then
               if(nfile.eq.1) then
                  nday = mod(idate1,10000)/100
                  write(31,400) n_slice,0,cday(nday),chmc(month),nyear,6
               else
                  write(31,400) n_slice,0,'01',chmc(month),nyear,6
               endif
            else if(ntype.eq.1) then
               write(31,401) n_slice,'01',chmc(month),nyear
            else if(ntype.eq.2) then
               write(31,402) n_month,'16',chmc(month),nyear
            endif
            write(31,500) 8
            write(31,650) 'h       ',np,'geopotential height (m)   '
            write(31,650) 't       ',np,'air temperature (deg, K)  '
            if(iproj.eq.'LAMCON') then
               write(31,651) 'u       ',np,'westerly wind (m/s)       '
               write(31,652) 'v       ',np,'southerly wind (m/s)      '
            else
               write(31,650) 'u       ',np,'westerly wind (m/s)       '
               write(31,650) 'v       ',np,'southerly wind (m/s)      '
            endif
            write(31,650) 'rh      ',np,'relative moisture (%)     '
            write(31,650) 'q       ',np,'specific moisture (kg/kg) '
            write(31,600) 'ps      ',   'surface pressure (hPa)    '
            write(31,600) 'slp     ',   'sea level pressure (hPa)  '
            write(31,700)
            close(31)
         endif
      enddo
            
      deallocate(sigma)
      deallocate(sig)
      deallocate(fin)
      deallocate(xlat)
      deallocate(xlon)
      deallocate(u)
      deallocate(v)
      deallocate(t)
      deallocate(q)
      deallocate(h)
      deallocate(ps)
      deallocate(tgb)
      deallocate(ht)
      deallocate(slp)
      deallocate(up)
      deallocate(vp)
      deallocate(tp)
      deallocate(qp)
      deallocate(hp)
            
 199  format(I4)
  10  format('dset ',A38)
  11  format('dset ',A36)
  12  format('dset ',A32)
  20  format('title RegCM domain information')
  30  format('options big_endian')
  40  format('options little_endian')
  50  format('undef -9999.')
 300  format('zdef ',I2,' levels ',30f7.2)
 400  format('tdef ',I4,' linear ',I2,'z',A2,A3,I4,' ',I2,'hr')
 401  format('tdef ',I4,' linear ',A2,A3,I4,' ','1dy')
 402  format('tdef ',I4,' linear ',A2,A3,I4,' ','1mo')
 500  format('vars ',I2)
 600  format(A8,'0 99 ',A26)           
 650  format(A8,I2,' 0 ',A26)
 651  format(A8,I2,' 33,100 ',A36)
 652  format(A8,I2,' 34,100 ',A36)
 700  format('endvars')
      return
      end

      subroutine sigma2p_ATM(iy,jx,kz,np,plev,ibyte &
                 ,idate0,idate1,idate2,Path_Output)
      implicit none
      integer iy,jx,kz,np,ibyte,idate0,idate1,idate2
      real*4  plev(np)
      character*128 Path_Output
      integer iiy,jjx,kkz
      integer mdate0,ibltyp,icup,ipptls,iboudy
      real*4  truelatL,truelatH
      real*4  dxsp,clat,clon,plat,plon
      real*4  dto,dtb,dtr,dtc
      integer iotyp
      character*6 iproj
      real*4, allocatable,save ::  sigma(:),sig(:)
      character*4 :: chy
      character*2 cday(31)
      data cday/'01','02','03','04','05','06','07','08','09','10', &
                '11','12','13','14','15','16','17','18','19','20', &
                '21','22','23','24','25','26','27','28','29','30','31'/
      integer nday
      character*2 chm(12)
      data chm/'01','02','03','04','05','06','07','08','09','10', &
               '11','12'/
      character*3 chmc(12)
      data chmc/'jan','feb','mar','apr','may','jun'  &
               ,'jul','aug','sep','oct','nov','dec'/
      character*14 filein
      character*16 fileout
      integer ntype,nfile,nyear,month,n_slice,mrec,nrec,nnn

      integer i,j,k,l,m,n

      real*4, allocatable,save ::  xlat(:,:),xlon(:,:)

      real*4, allocatable,save :: u(:,:,:),v(:,:,:),t(:,:,:),q(:,:,:)
      real*4, allocatable,save :: w(:,:,:),c(:,:,:),h(:,:,:)
      real*4, allocatable,save :: ps(:,:),tgb(:,:),ht(:,:)
      real*4, allocatable,save :: slp(:,:),tpr(:,:),swt(:,:),rno(:,:)

      real*4, allocatable,save :: up(:,:,:),vp(:,:,:),tp(:,:,:)
      real*4, allocatable,save :: qp(:,:,:),hp(:,:,:),wp(:,:,:)
      real*4, allocatable,save :: cp(:,:,:)

      real*4, save :: ptop
      real*4, parameter, save :: rgas = 287.04
      real*4, parameter, save :: grav = 9.80616
      real*4, parameter, save :: bltop = 0.96
      real*4, parameter, save :: tlapse = -6.5E-3

      integer n_month
      logical there
      real*4  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real*4  centerj,centeri
      integer ny,nx

      integer IDATE

      allocate(sigma(kz+1))
      allocate(sig(kz))
      allocate(xlat(jx-2,iy-2))
      allocate(xlon(jx-2,iy-2))

      allocate(u(jx-2,iy-2,kz))
      allocate(v(jx-2,iy-2,kz))
      allocate(t(jx-2,iy-2,kz))
      allocate(q(jx-2,iy-2,kz))
      allocate(w(jx-2,iy-2,kz))
      allocate(c(jx-2,iy-2,kz))
      allocate(h(jx-2,iy-2,kz))

      allocate(ps(jx-2,iy-2))
      allocate(tgb(jx-2,iy-2))
      allocate(ht(jx-2,iy-2))
      allocate(slp(jx-2,iy-2))
      allocate(tpr(jx-2,iy-2))   ! 
      allocate(swt(jx-2,iy-2))
      allocate(rno(jx-2,iy-2))

      allocate(up(jx-2,iy-2,np))
      allocate(vp(jx-2,iy-2,np))
      allocate(tp(jx-2,iy-2,np))
      allocate(qp(jx-2,iy-2,np))
      allocate(wp(jx-2,iy-2,np))
      allocate(cp(jx-2,iy-2,np))
      allocate(hp(jx-2,iy-2,np))

      inquire(file=trim(Path_Output)//'OUT_HEAD',exist=there)
      if(.not.there) then
         write(*,*) trim(Path_Output)//'OUT_HEAD',' is not avaiable'
         stop
      endif
      open(10,file=trim(Path_Output)//'OUT_HEAD' &
             ,recl=(jx-2)*(iy-2)*ibyte,access='direct')
      read(10,rec=1) mdate0,ibltyp,icup,ipptls,iboudy  &
                    ,iiy,jjx,kkz,(sigma(k),k=kz+1,1,-1)   &
                    ,dxsp,ptop,clat,clon,plat,plon     &
                    ,iproj,dto,dtb,dtr,dtc,iotyp,truelatL,truelatH
      if(iiy.ne.iy.or.jjx.ne.jx.or.kkz.ne.kz) then
         write(*,*) 'iy,jx,kz in parameter = ',iy,jx,kz
         write(*,*) 'iy,jx,kz in OUT_HEAD ',iiy,jjx,kkz
         write(*,*) 'They are not consistent'
         stop
      endif

      ptop = ptop*10.
      do k=1,kz
         sig(k) = (sigma(k)+sigma(k+1))*0.5
      enddo

      read(10,rec=2) ht          ! ht
      read(10,rec=6) xlat        ! xlat
      read(10,rec=7) xlon        ! xlon
      close(10)
      
      if(idate1.ge.1000010100) then        ! original file
         n_month = (idate2/1000000-idate1/1000000)*12 &
                 + (mod(idate2/10000,100)-mod(idate1/10000,100))   
         if(mod(idate2,10000).gt.0100) n_month = n_month+1
         ntype = 0
      else if(idate1.ge.10000101) then     ! daily mean file
         n_month = (idate2/10000-idate1/10000)*12 &
                 + (mod(idate2/100,100)-mod(idate1/100,100))   
         if(mod(idate2,100).gt.01) n_month = n_month+1
         ntype = 1
      else if(idate1.ge.100001) then       ! monthly mean file
         n_month = (idate2/100-idate1/100)*12 &
                 + (mod(idate2,100)-mod(idate1,100))+1 
         ntype = 2
      endif

      do nfile=1,n_month
         if(ntype.eq.0) then
            nyear = idate1/1000000
            month = idate1/10000-nyear*100 +nfile-1
            nyear = nyear + (month-1)/12
            month = mod(month,12)
            if(month.eq.0) month = 12

            if(nfile.eq.1.or.month.eq.1) then
               write(*,*)'sigma2p: ATM Orig.',nyear,month
            else
               write(*,*)'                  ',nyear,month
            endif

            if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
               month.eq.7.or.month.eq.8.or.month.eq.10.or. &
               month.eq.12) then
               n_slice = nint(31*24/dto)
            else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                    month.eq.11) then
               n_slice = nint(30*24/dto)
            else
               n_slice = nint(28*24/dto)
               if(mod(nyear,4).eq.0) n_slice = nint(29*24/dto)
               if(mod(nyear,100).eq.0) n_slice = nint(28*24/dto)
               if(mod(nyear,400).eq.0) n_slice = nint(29*24/dto)
            endif
            if(nfile.eq.n_month.and.mod(idate2/100-1,100).ne.0) &
            n_slice=min(n_slice,nint(mod(idate2/100-1,100)*24./dto))
            if(nfile.eq.1.and.idate0.eq.idate1) n_slice=n_slice+1
            write(chy,199) nyear
            filein = 'ATM.'//chy//chm(month)//'0100'
            fileout= 'ATM_P.'//chy//chm(month)//'0100'

            inquire(file=trim(Path_Output)//filein,exist=there)
            if(.not.there) then
               write(*,*) trim(Path_Output)//filein &
                     ,' is not avaiable'
               stop
            endif
            if(iotyp.eq.1) then
               open(10,file=trim(Path_Output)//filein  &
           ,form='unformatted',recl=(iy-2)*(jx-2)*ibyte,access='direct')
               mrec = 0
            else
              open(10,file=trim(Path_Output)//filein,form='unformatted')
            endif
            open(20,file=trim(Path_Output)//fileout &
           ,form='unformatted',recl=(iy-2)*(jx-2)*ibyte,access='direct')
            nrec = 0
         else if(ntype.eq.1) then
            nyear = idate1/10000
            month = idate1/100-nyear*100 +nfile-1
            nyear = nyear + (month-1)/12
            month = mod(month,12)
            if(month.eq.0) month = 12

            if(nfile.eq.1.or.month.eq.1) then
               write(*,*)'sigma2p: ATM Daily',nyear,month
            else
               write(*,*)'                  ',nyear,month
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
            write(chy,199) nyear
            filein = 'ATM.'//chy//chm(month)//'01'
            fileout= 'ATM_P.'//chy//chm(month)//'01'

            inquire(file=trim(Path_Output)//filein(1:12),exist=there)
            if(.not.there) then
               write(*,*) trim(Path_Output)//filein(1:12) &
                        ,' is not avaiable'
               stop
            endif
            open(10,file=trim(Path_Output)//filein(1:12) &
           ,form='unformatted',recl=(iy-2)*(jx-2)*ibyte,access='direct')
            mrec = 0
            open(20,file=trim(Path_Output)//'_'//fileout(1:14) &
           ,form='unformatted',recl=(iy-2)*(jx-2)*ibyte,access='direct')
            nrec = 0
         else if(ntype.eq.2) then
            nyear = idate1/100
            month = idate1-nyear*100 +nfile-1
            nyear = nyear + (month-1)/12
            month = mod(month,12)
            if(month.eq.0) month = 12

            if(nfile.eq.1.or.month.eq.1) then
               write(*,*)'sigma2p: ATM Month',nyear,month
            else
               write(*,*)'                  ',nyear,month
            endif

            n_slice = 1

            if(nfile.eq.1) then
               filein = 'ATM.mon'
               fileout= 'ATM_P.mon'

               inquire(file=trim(Path_Output)//filein(1:7),exist=there)
               if(.not.there) then
                  write(*,*) trim(Path_Output)//filein(1:7) &
                        ,' is not avaiable'
                  stop
               endif
               open(10,file=trim(Path_Output)//filein(1:7) &
                      ,form='unformatted' &
                      ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
               mrec = 0
               open(20,file=trim(Path_Output)//fileout(1:9) &
                    ,form='unformatted' &
                    ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
               nrec = 0
            endif
         endif

         do nnn=1,n_slice
            if(ntype.eq.0) then
               if(iotyp.eq.1) then
                  do l=kz,1,-1
                     mrec = mrec+1
                     read(10,rec=mrec) ((u(j,i,l),j=1,jx-2),i=1,iy-2)
                  enddo
                  do l=kz,1,-1
                     mrec = mrec+1
                     read(10,rec=mrec) ((v(j,i,l),j=1,jx-2),i=1,iy-2)
                  enddo
                  do l=kz,1,-1
                     mrec = mrec+1
                     read(10,rec=mrec) ((w(j,i,l),j=1,jx-2),i=1,iy-2)
                  enddo
                  do l=kz,1,-1
                     mrec = mrec+1
                     read(10,rec=mrec) ((t(j,i,l),j=1,jx-2),i=1,iy-2)
                  enddo
                  do l=kz,1,-1
                     mrec = mrec+1
                     read(10,rec=mrec) ((q(j,i,l),j=1,jx-2),i=1,iy-2)
                  enddo
                  do l=kz,1,-1
                     mrec = mrec+1
                     read(10,rec=mrec) ((c(j,i,l),j=1,jx-2),i=1,iy-2)
                  enddo
                  mrec = mrec+1
                  read(10,rec=mrec) ps
                  mrec = mrec+1
                  read(10,rec=mrec) tpr
                  mrec = mrec+1
                  read(10,rec=mrec) tgb
                  mrec = mrec+1
                  read(10,rec=mrec) swt
                  mrec = mrec+1
                  read(10,rec=mrec) rno
               else if(iotyp.eq.2) then
                  read(10) IDATE
!                 write(*,*) 'IDATE = ',IDATE
                  do l=kz,1,-1
                     read(10) ((u(j,i,l),j=1,jx-2),i=1,iy-2)
                  enddo
                  do l=kz,1,-1
                     read(10) ((v(j,i,l),j=1,jx-2),i=1,iy-2)
                  enddo
                  do l=kz,1,-1
                     read(10) ((w(j,i,l),j=1,jx-2),i=1,iy-2)
                  enddo
                  do l=kz,1,-1
                     read(10) ((t(j,i,l),j=1,jx-2),i=1,iy-2)
                  enddo
                  do l=kz,1,-1
                     read(10) ((q(j,i,l),j=1,jx-2),i=1,iy-2)
                  enddo
                  do l=kz,1,-1
                     read(10) ((c(j,i,l),j=1,jx-2),i=1,iy-2)
                  enddo
                  read(10) ps
                  read(10) tpr
                  read(10) tgb
                  read(10) swt
                  read(10) rno
               endif
            else
               do l=kz,1,-1
                  mrec = mrec+1
                  read(10,rec=mrec) ((u(j,i,l),j=1,jx-2),i=1,iy-2)
               enddo
               do l=kz,1,-1
                  mrec = mrec+1
                  read(10,rec=mrec) ((v(j,i,l),j=1,jx-2),i=1,iy-2)
               enddo
               do l=kz,1,-1
                  mrec = mrec+1
                  read(10,rec=mrec) ((w(j,i,l),j=1,jx-2),i=1,iy-2)
               enddo
               do l=kz,1,-1
                  mrec = mrec+1
                  read(10,rec=mrec) ((t(j,i,l),j=1,jx-2),i=1,iy-2)
               enddo
               do l=kz,1,-1
                  mrec = mrec+1
                  read(10,rec=mrec) ((q(j,i,l),j=1,jx-2),i=1,iy-2)
               enddo
               do l=kz,1,-1
                  mrec = mrec+1
                  read(10,rec=mrec) ((c(j,i,l),j=1,jx-2),i=1,iy-2)
               enddo
               mrec = mrec+1
               read(10,rec=mrec) ps
               mrec = mrec+1
               read(10,rec=mrec) tpr
               mrec = mrec+1
               read(10,rec=mrec) tgb
               mrec = mrec+1
               read(10,rec=mrec) swt
               mrec = mrec+1
               read(10,rec=mrec) rno
            endif

          ! to calculate Heights on sigma surfaces.
            call htsig(t,h,ps,ht,sig,jx-2,iy-2,kz,PTOP,RGAS,GRAV)

          ! to calculate Sea-Level Pressure using
          !  1. ERRICO's solution described in height
          !  2. a simple formulae
          !  3. MM5 method
            call slpres(h,t,ps,ht,tgb,slp,sig,jx-2,iy-2,kz  &
                       ,RGAS,GRAV,BLTOP,TLAPSE)

          ! to interpolate H,U,V,T,Q and QC
          !  1. For Heights
            call height(hp,h,t,ps,ht,sig,jx-2,iy-2,kz,plev,np &
                       ,PTOP,RGAS,GRAV,BLTOP,TLAPSE)
            do k=1,np
               nrec=nrec+1
               write(20,rec=nrec)((hp(j,i,k),j=1,jx-2),i=1,iy-2)
            enddo
          !  2. For Zonal and Meridional Winds
            call intlin(up,u,ps,sig,jx-2,iy-2,kz,plev,np,PTOP)
            call intlin(vp,v,ps,sig,jx-2,iy-2,kz,plev,np,PTOP)
            call intlin(wp,w,ps,sig,jx-2,iy-2,kz,plev,np,PTOP)
          !  3. For Temperatures
            call intlog(tp,t,ps,sig,jx-2,iy-2,kz,plev,np &
                       ,PTOP,RGAS,GRAV,BLTOP,TLAPSE)
            do k=1,np
               nrec=nrec+1
               write(20,rec=nrec)((tp(j,i,k),j=1,jx-2),i=1,iy-2)
            enddo
            do k=1,np
               nrec=nrec+1
               write(20,rec=nrec)((up(j,i,k),j=1,jx-2),i=1,iy-2)
            enddo
            do k=1,np
               nrec=nrec+1
               write(20,rec=nrec)((vp(j,i,k),j=1,jx-2),i=1,iy-2)
            enddo
            do k=1,np
               nrec=nrec+1
               write(20,rec=nrec)((wp(j,i,k),j=1,jx-2),i=1,iy-2)
            enddo
          !  4. For Moisture qv
            call humid1(t,q,ps,sig,jx-2,iy-2,kz,ptop)
            call intlin(qp,q,ps,sig,jx-2,iy-2,kz,plev,np,PTOP)
            do k=1,np
               nrec=nrec+1
               write(20,rec=nrec)((qp(j,i,k)*100.,j=1,jx-2),i=1,iy-2)
            enddo

            call humid2(tp,qp,plev,jx-2,iy-2,np)
            do k=1,np
               nrec=nrec+1
               write(20,rec=nrec)((qp(j,i,k),j=1,jx-2),i=1,iy-2)
            enddo
            call intlin(cp,c,ps,sig,jx-2,iy-2,kz,plev,np,PTOP)
            do k=1,np
               nrec=nrec+1
               write(20,rec=nrec)((cp(j,i,k),j=1,jx-2),i=1,iy-2)
            enddo
            nrec=nrec+1
            write(20,rec=nrec) ps
            nrec=nrec+1
            write(20,rec=nrec) slp
            nrec=nrec+1
            write(20,rec=nrec) tpr
            nrec=nrec+1
            write(20,rec=nrec) tgb
            nrec=nrec+1
            write(20,rec=nrec) swt
            nrec=nrec+1
            write(20,rec=nrec) rno
         enddo
         if(.not.ntype.eq.2) close(20)
         if(.not.ntype.eq.2) close(10)
         if(ntype.eq.0.or.ntype.eq.1.or.(ntype.eq.2.and.nfile.eq.1))then
            if(ntype.eq.0) then
               inquire(file=trim(Path_Output)//fileout(1:16)//'.ctl' &
                      ,exist=there)
               if(there) then
                open(31,file=trim(Path_Output)//fileout(1:16)//'.ctl' &
                      ,form='formatted',status='replace')
               else
                open(31,file=trim(Path_Output)//fileout(1:16)//'.ctl' &
                      ,form='formatted',status='new')
               endif
               write(31,10) '^'//fileout(1:16)
            else if(ntype.eq.1) then
               inquire(file=trim(Path_Output)//fileout(1:14)//'.ctl' &
                      ,exist=there)
               if(there) then
                open(31,file=trim(Path_Output)//fileout(1:14)//'.ctl' &
                      ,form='formatted',status='replace')
               else
                open(31,file=trim(Path_Output)//fileout(1:14)//'.ctl' &
                      ,form='formatted',status='new')
               endif
               write(31,11) '^'//fileout(1:14)
            else if(ntype.eq.2) then
               inquire(file=trim(Path_Output)//fileout(1:6)//'.ctl' &
                         ,exist=there)
               if(there) then
                  open(31,file=trim(Path_Output)//fileout(1:6)//'.ctl' &
                      ,form='formatted',status='replace')
               else
                  open(31,file=trim(Path_Output)//fileout(1:6)//'.ctl' &
                      ,form='formatted',status='new')
               endif
               write(31,12) '^'//fileout(1:6)//'.mon'
            endif
            write(31,20)
            write(31,30)
            write(31,50)
            if(iproj.eq.'LAMCON'.or.iproj.eq.'ROTMER') then
               alatmin= 999999.
               alatmax=-999999.
               do j=1,jx-2
                  if(xlat(j,1).lt.alatmin) alatmin=xlat(j,1)
                  if(xlat(j,iy-2).gt.alatmax) alatmax=xlat(j,iy-2)
               enddo
               alonmin= 999999.
               alonmax=-999999.
               do i=1,iy-2
               do j=1,jx-2
                  if(clon.ge.0.0) then
                     if(xlon(j,i).ge.0.0) then
                        alonmin = amin1(alonmin,xlon(j,i))
                        alonmax = amax1(alonmax,xlon(j,i))
                     else
                        if(abs(clon-xlon(j,i)).lt. &
                           abs(clon-(xlon(j,i)+360.))) then
                           alonmin = amin1(alonmin,xlon(j,i))
                           alonmax = amax1(alonmax,xlon(j,i))
                        else
                           alonmin = amin1(alonmin,xlon(j,i)+360.)
                           alonmax = amax1(alonmax,xlon(j,i)+360.)
                        endif
                     endif
                  else
                     if(xlon(j,i).lt.0.0) then
                        alonmin = amin1(alonmin,xlon(j,i))
                        alonmax = amax1(alonmax,xlon(j,i))
                     else
                        if(abs(clon-xlon(j,i)).lt. &
                           abs(clon-(xlon(j,i)-360.))) then
                           alonmin = amin1(alonmin,xlon(j,i))
                           alonmax = amax1(alonmax,xlon(j,i))
                        else
                           alonmin = amin1(alonmin,xlon(j,i)-360.)
                           alonmax = amax1(alonmax,xlon(j,i)-360.)
                        endif
                     endif
                  endif
               enddo
               enddo
               rlatinc=dxsp/111./2.
               rloninc=dxsp/111./2.
               ny=2+nint(abs(alatmax-alatmin)/rlatinc)
               nx=1+nint(abs((alonmax-alonmin)/rloninc))
               centerj=(jx-2)/2.
               centeri=(iy-2)/2.
            endif
            if(iproj.eq.'LAMCON') then        ! Lambert projection
               write(31,100) jx-2,iy-2,clat,clon,centerj,centeri, &
                      truelatL,truelatH,clon,dxsp*1000.,dxsp*1000.
 100  format('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
               write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
               write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
            elseif(iproj.eq.'POLSTR') then    !
            elseif(iproj.eq.'NORMER') then
               write(31,200)  jx-2,xlon(1,1),xlon(2,1)-xlon(1,1)
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
               write(31,210) iy-2
 210  format('ydef ',I3,' levels')
               write(31,220) (xlat(1,i),i=1,iy-2)
 220  format(10f7.2)
            elseif(iproj.eq.'ROTMER') then
            if(nfile.eq.1.or.month.eq.1) then
               write(*,*) 'Note that rotated Mercartor (ROTMER)' &
                        ,' projections are not supported by GrADS.'
               write(*,*) '  Although not exact, the eta.u projection' &
                        ,' in GrADS is somewhat similar.'
            write(*,*) ' FERRET, however, does support this projection.'
            endif
               write(31,230) jx-2,iy-2,plon,plat,dxsp/111. &
                                           ,dxsp/111.*.95238
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
               write(31,110) nx+2,alonmin-rloninc,rloninc
               write(31,120) ny+2,alatmin-rlatinc,rlatinc
            else
               write(*,*) 'Are you sure your map projection is right ?'
               stop
            endif
            write(31,300) np, (plev(k),k=1,np)
            if(ntype.eq.0) then
               if(nfile.eq.1) then
                  nday = mod(idate1,10000)/100
                  if(idate0.eq.idate1) then
                     write(31,400) &
                     n_slice,0,'01',chmc(month),nyear,nint(dto)
                  else
                     write(31,400) n_slice, &
               nint(dto),cday(nday),chmc(month),nyear,nint(dto)
                  endif
               else
                  write(31,400) n_slice &
                    ,nint(dto),'01',chmc(month),nyear,nint(dto)
               endif
            else if(ntype.eq.1) then
               write(31,401) n_slice,'01',chmc(month),nyear
            else if(ntype.eq.2) then
               write(31,402) n_month,'16',chmc(month),nyear
            endif
            write(31,500) 14
            write(31,650) 'h       ',np,'geopotential height (m)   '
            write(31,650) 't       ',np,'air temperature (deg, K)  '
            if(iproj.eq.'LAMCON') then
               write(31,651) 'u       ',np,'westerly wind (m/s)       '
               write(31,652) 'v       ',np,'southerly wind (m/s)      '
            else
               write(31,650) 'u       ',np,'westerly wind (m/s)       '
               write(31,650) 'v       ',np,'southerly wind (m/s)      '
            endif
            write(31,650) 'omega   ',np,'vertical velocity of p    '
            write(31,650) 'rh      ',np,'relative moisture (%)     '
            write(31,650) 'q       ',np,'specific moisture (kg/kg) '
            write(31,650) 'qc      ',np,'cloud water       (kg/kg) '
            write(31,600) 'ps      ',   'surface pressure (hPa)    '
            write(31,600) 'slp     ',   'sea level pressure (hPa)  '
            write(31,600) 'tpr     ',   'precipitation, (mm/day)   '
            write(31,600) 'tgb     ',   'lower soil temperture, K  '
            write(31,600) 'swt     ',   'total soil water in mm H2O'
            write(31,600) 'rno     ',   'accumulated infiltration  '
            write(31,700)
            close(31)
         endif
      enddo
            
      deallocate(sigma)
      deallocate(sig)
      deallocate(xlat)
      deallocate(xlon)
      deallocate(u)
      deallocate(v)
      deallocate(t)
      deallocate(q)
      deallocate(w)
      deallocate(c)
      deallocate(h)
      deallocate(ps)
      deallocate(tgb)
      deallocate(ht)
      deallocate(slp)
      deallocate(tpr)
      deallocate(swt)
      deallocate(rno)
      deallocate(up)
      deallocate(vp)
      deallocate(tp)
      deallocate(qp)
      deallocate(wp)
      deallocate(cp)
      deallocate(hp)
            
 199  format(I4)
  10  format('dset ',A38)
  11  format('dset ',A36)
  12  format('dset ',A32)
  20  format('title RegCM domain information')
  30  format('options big_endian')
  40  format('options little_endian')
  50  format('undef -9999.')
 300  format('zdef ',I2,' levels ',30f7.2)
 400  format('tdef ',I4,' linear ',I2,'z',A2,A3,I4,' ',I2,'hr')
 401  format('tdef ',I4,' linear ',A2,A3,I4,' ','1dy')
 402  format('tdef ',I4,' linear ',A2,A3,I4,' ','1mo')
 500  format('vars ',I2)
 600  format(A8,'0 99 ',A26)           
 650  format(A8,I2,' 0 ',A26)
 651  format(A8,I2,' 33,100 ',A36)
 652  format(A8,I2,' 34,100 ',A36)
 700  format('endvars')
      return
      end

      SUBROUTINE HUMID2(T,Q,PRESLV,IM,JM,KP)
      implicit none
      INTEGER IM,JM,KP
      REAL    TR,QMIN
      PARAMETER (TR=1./273.16)
      PARAMETER (QMIN=0.0)   ! MINIMUM VALUE OF SPECIFIC HUMIDITY
      REAL    T(IM,JM,KP),Q(IM,JM,KP)
      REAL    PRESLV(KP)
      INTEGER I,J,K
      REAL    HL,SATVP,QS
!
!  THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
!  DATA ON SIGMA LEVELS
!
      DO K=1,KP
        DO J=1,JM
          DO I=1,IM
            HL=597.3-.566*(T(I,J,K)-273.16)
            SATVP=6.11*EXP(9.045*HL*(TR-1./T(I,J,K)))
            QS=.622*SATVP/(PRESLV(K)-SATVP)           ! preslv (hPa)
            IF (Q(I,J,K).LT.QMIN) Q(I,J,K)=QMIN       ! SPECIFIED MINIMUM
            Q(I,J,K)=Q(I,J,K)*QS
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE HUMID1(T,Q,PS,SIGMA,IM,JM,KM,PTOP)
      implicit none
      INTEGER IM,JM,KM
      REAL    TR,QMIN
      PARAMETER (TR=1./273.16)
      PARAMETER (QMIN=0.0)   ! MINIMUM VALUE OF SPECIFIC HUMIDITY
      REAL    T(IM,JM,KM),Q(IM,JM,KM)
      REAL    PS(IM,JM)
      REAL    SIGMA(KM)
      REAL    PTOP
      INTEGER I,J,K
      REAL    HL,SATVP,QS,P
!
!  THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
!  DATA ON SIGMA LEVELS
!
      DO K=1,KM
        DO J=1,JM
          DO I=1,IM
            P=SIGMA(K)*(PS(I,J)-PTOP)+PTOP
            HL=597.3-.566*(T(I,J,K)-273.16)           ! LATENT HEAT OF EVAP.
            SATVP=6.11*EXP(9.045*HL*(TR-1./T(I,J,K))) ! SATURATION VAP PRESS.
            QS=.622*SATVP/(P-SATVP)                   ! SAT. MIXING RATIO
            Q(I,J,K)=Q(I,J,K)/QS
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
      SUBROUTINE INTLOG(FP,F,PSTAR,SIG,IM,JM,KM,P,KP &
                       ,PTOP,RGAS,GRAV,BLTOP,TLAPSE)
      implicit none
      INTEGER IM,JM,KM,KP
      REAL    FP(IM,JM,KP),F(IM,JM,KM)
      REAL    PSTAR(IM,JM)
      REAL    SIG(KM),P(KP)
      REAL    PTOP,RGAS,GRAV,BLTOP,TLAPSE
      INTEGER I,J,K,N
      INTEGER K1,K1P,KBC
      REAL    SIGP,WP,W1
!
!  INTLOG IS FOR VERTICAL INTERPOLATION OF T.  THE INTERPOLATION IS
!        LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
!        THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!        WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
!        CONSIDERED TO HAVE A LAPSE RATE OF TLAPSE (K/M), AND THE
!        THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
!        TWO EXTREME TEMPERATURES IN THE LAYER.

!
!** FIND FIRST SIGMA LEVEL ABOVE BOUNDARY LAYER (LESS THAN SIG=BLTOP)
      DO K=1,KM
        IF(SIG(K).LT.BLTOP) KBC = K
      ENDDO
      DO J=1,JM
        DO I=1,IM
          DO N=1,KP
            SIGP = (P(N)-PTOP) / (PSTAR(I,J)-PTOP)
            K1=0
            DO K=1,KM
              IF (SIGP.GT.SIG(K)) K1=K
            ENDDO
            IF(SIGP.LE.SIG(1)) THEN
              FP(I,J,N) = F(I,J,1)
            ELSE IF((SIGP.GT.SIG(1)).AND.(SIGP.LT.SIG(KM))) THEN
              K1P = K1 + 1
              WP  = LOG(SIGP/SIG(K1)) / LOG(SIG(K1P)/SIG(K1))
              W1  = 1. - WP
              FP(I,J,N)= W1*F(I,J,K1) + WP*F(I,J,K1P)
            ELSE IF((SIGP.GE.SIG(KM)).AND.(SIGP.LE.1.))THEN
              FP(I,J,N)= F(I,J,KM)
            ELSE IF(SIGP.GT.1.) THEN
              FP(I,J,N) = F(I,J,KBC) &
              * EXP(-RGAS*TLAPSE*LOG(SIGP/SIG(KBC))/GRAV)
!                  ***** FROM R. ERRICO, SEE ROUTINE HEIGHT *****
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      
      RETURN
      END
!
      SUBROUTINE INTLIN(FP,F,PSTAR,SIG,IM,JM,KM,P,KP,PTOP)
      implicit none
      INTEGER IM,JM,KM,KP
      REAL    FP(IM,JM,KP),F(IM,JM,KM)
      REAL    PSTAR(IM,JM)
      REAL    SIG(KM),P(KP)
      REAL    PTOP
      INTEGER I,J,K,N
      INTEGER K1,K1P
      REAL    SIGP,WP,W1
!
!  INTLIN IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE HUMIDITY.
!        THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION IS
!        NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.

      DO J=1,JM
        DO I=1,IM
          DO N=1,KP
            SIGP = (P(N)-PTOP) / (PSTAR(I,J)-PTOP)
            K1=0
            DO K=1,KM
              IF (SIGP.GT.SIG(K)) K1=K
            ENDDO
            IF(SIGP.LE.SIG(1)) THEN
              FP(I,J,N) = F(I,J,1)
            ELSE IF((SIGP.GT.SIG(1)).AND.(SIGP.LT.SIG(KM))) THEN
              K1P = K1 + 1
              WP  = (SIGP-SIG(K1))/(SIG(K1P)-SIG(K1))
              W1  = 1.-WP
              FP(I,J,N)  = W1*F(I,J,K1)+WP*F(I,J,K1P)
            ELSE IF(SIGP.GE.SIG(KM)) THEN
              FP(I,J,N)  = F(I,J,KM)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
      SUBROUTINE HEIGHT(HP,H,T,PSTAR,HT,SIG,IM,JM,KM,P,KP &
                       ,PTOP,RGAS,GRAV,BLTOP,TLAPSE)

!  HEIGHT DETERMINES THE HEIGHT OF PRESSURE LEVELS.
!     ON INPUT:
!        H AND T ARE HEIGHT AND TEMPERATURE ON SIGMA, RESPECTIVELY.
!        PSTAR = SURFACE PRESSURE - MODEL TOP PRESSURE.
!        SIG = SIGMA LEVELS.
!        P = PRESSURE LEVELS DESIRED.
!     ON OUTPUT:
!        ALL FIELDS EXCEPT H ARE UNCHANGED.
!        H HAS HEIGHT FIELDS AT KP PRESSURE LEVELS.
!
!  FOR UPWARD EXTRAPOLATION, T IS CONSIDERED TO HAVE 0 VERITCAL DERIV.
!  FOR DOWNWARD EXTRAPOLATION, T HAS LAPSE RATE OF TLAPSE (K/KM)
!     AND EXTRAPOLATION IS DONE FROM THE LOWEST SIGMA LEVEL ABOVE
!     THE BOUNDARY LAYER (TOP ARBITRARILY TAKEN AT SIGMA = BLTOP).
!     EQUATION USED IS EXACT SOLUTION TO HYDROSTATIC RELATION,
!     GOTTEN FROM R. ERRICO (ALSO USED IN SLPRES ROUTINE):
!      Z = Z0 - (T0/TLAPSE) * (1.-EXP(-R*TLAPSE*LN(P/P0)/G))
!
      implicit none
      INTEGER IM,JM,KM,KP
      REAL    T(IM,JM,KM),H(IM,JM,KM),HP(IM,JM,KP)
      REAL    PSTAR(IM,JM),HT(IM,JM)
      REAL    SIG(KM),P(KP)
      REAL    PTOP,RGAS,GRAV,BLTOP,TLAPSE
      REAL    PSIG(100)
      INTEGER I,J,K,KBC,N,KT,KB
      REAL    PSFC,TEMP,WT,WB
!
      DO K=1,KM
         IF(SIG(K).LT.BLTOP) THEN
           KBC=K
         ENDIF
      ENDDO
!     PRINT *,'FIRST SIGMA LEVEL ABOVE BNDY LAYER:', SIG(KBC)
!
      DO J=1,JM
        DO I=1,IM
          DO K=1,KM
            PSIG(K) = SIG(K) * (PSTAR(I,J)-PTOP) + PTOP
          ENDDO
          PSFC = PSTAR(I,J)
          DO N = 1,KP
            KT = 1
            DO K=1,KM
              IF (PSIG(K).LT.P(N)) KT=K
            ENDDO
            KB = KT + 1
            IF(P(N).LE.PSIG(1)) THEN
              TEMP = T(I,J,1)
              HP(I,J,N) =H(I,J,1)+RGAS*TEMP*LOG(PSIG(1)/P(N))/GRAV
            ELSE IF((P(N).GT.PSIG(1)).AND.(P(N).LT.PSIG(KM))) THEN
              WT = LOG(PSIG(KB)/P(N)) / LOG(PSIG(KB)/PSIG(KT))
              WB = LOG(P(N)/PSIG(KT)) / LOG(PSIG(KB)/PSIG(KT))
              TEMP = WT * T(I,J,KT) + WB * T(I,J,KB)
              TEMP = ( TEMP + T(I,J,KB) ) / 2.
              HP(I,J,N) =H(I,J,KB)+RGAS*TEMP*LOG(PSIG(KB)/P(N))/GRAV
            ELSE IF((P(N).GE.PSIG(KM)).AND.(P(N).LE.PSFC)) THEN
              TEMP = T(I,J,KM)
              HP(I,J,N) =HT(I,J)+RGAS*TEMP*LOG(PSFC/P(N))/GRAV
            ELSE IF(P(N).GT.PSFC) THEN
              TEMP = T(I,J,KBC) - TLAPSE * (H(I,J,KBC)-HT(I,J))
              HP(I,J,N) =HT(I,J)-(TEMP/TLAPSE)  &
                    * ( 1.-EXP(-RGAS*TLAPSE*LOG(P(N)/PSFC)/GRAV))
!
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
      SUBROUTINE SLPRES(H,T,PSTAR,HT,TG,SLP,SIG,IM,JM,KM  &
                       ,RGAS,GRAV,BLTOP,TLAPSE)
      implicit none
      INTEGER IM,JM,KM
      REAL    T(IM,JM,KM),H(IM,JM,KM)
      REAL    PSTAR(IM,JM),HT(IM,JM),TG(IM,JM)
      REAL    SLP(IM,JM)
      REAL    SIG(KM)
      REAL    RGAS,GRAV,BLTOP,TLAPSE
      INTEGER KBC,I,J,K
      REAL    TSFC
!
      DO K=1,KM
         IF(SIG(K).LT.BLTOP) THEN
           KBC=K
         ENDIF
      ENDDO
      DO J=1,JM
         DO I=1,IM
            TSFC = T(I,J,KBC)-TLAPSE*(H(I,J,KBC)-HT(I,J))
            SLP(I,J) = PSTAR(I,J)  &
            * EXP( -GRAV/(RGAS*TLAPSE)*LOG(1.-HT(I,J)*TLAPSE/TSFC))
         ENDDO
      ENDDO

!     DO J=1,JM
!        DO I=1,IM
!           SLP2(I,J) = PSTAR(I,J)
!    &      * EXP( GRAV*HT(I,J)/(RGAS*0.5*(TG(I,J)+288.15)))
!        ENDDO
!     ENDDO

      RETURN
      END
!
      SUBROUTINE HTSIG(T,H,PSTAR,HT,SIG,IM,JM,KM,PTOP,RGAS,GRAV)
      implicit none
      INTEGER IM,JM,KM
      REAL    T(IM,JM,KM),H(IM,JM,KM)
      REAL    PSTAR(IM,JM),HT(IM,JM)
      REAL    SIG(KM)
      REAL    PTOP,RGAS,GRAV
      INTEGER I,J,K
      REAL    TBAR
!
      DO J=1,JM
      DO I=1,IM
         H(I,J,KM) = HT(I,J) + RGAS/GRAV*T(I,J,KM) &
                   * LOG(PSTAR(I,J)/((PSTAR(I,J)-PTOP)*SIG(KM)+PTOP))
      ENDDO
      ENDDO
      DO K=KM-1,1,-1
      DO J=1,JM
      DO I=1,IM
         TBAR = 0.5*( T(I,J,K)+T(I,J,K+1) )
         H(I,J,K) = H(I,J,K+1) +RGAS/GRAV*TBAR  &
                  * LOG(((PSTAR(I,J)-PTOP)*SIG(K+1)+PTOP)  &
                       /((PSTAR(I,J)-PTOP)*SIG(K)+PTOP))
      ENDDO
      ENDDO
      ENDDO
      RETURN
      END
