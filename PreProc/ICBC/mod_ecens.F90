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

module mod_ecens

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_stdio
  use mod_memutil
  use mod_grid
  use mod_write
  use mod_interp
  use mod_vertint
  use mod_hgt
  use mod_humid
  use mod_mksst
  use mod_uvrot
  use mod_vectutil
  use mod_message
  use mod_nchelper
  use netcdf

  implicit none

  private

  integer(ik4) :: mlev , nlat , nlon , ntime
  integer(ik4) , parameter :: nplev = 18

  real(rkx) , pointer , dimension(:) :: vlat
  real(rkx) , pointer , dimension(:) :: rlat
  real(rkx) , pointer , dimension(:) :: vlon
  real(rkx) , pointer , dimension(:) :: ak , bk
  real(rkx) , pointer , dimension(:,:,:) :: work

  real(rkx) , pointer , dimension(:,:,:) :: bb
  real(rkx) , pointer , dimension(:,:,:) :: b2
  real(rkx) , pointer , dimension(:,:,:) :: d2
  real(rkx) , pointer , dimension(:,:,:) :: b3
  real(rkx) , pointer , dimension(:,:,:) :: d3
  real(rkx) , pointer , dimension(:,:,:) :: pp3d , z1

  real(rkx) , pointer , dimension(:,:) :: zs2
  real(rkx) , pointer , dimension(:,:) :: ps2
  real(rkx) , pointer , dimension(:,:,:) :: q2 , t2 , u2 , v2
  real(rkx) , pointer , dimension(:,:,:) :: tp , qp , hp
  real(rkx) , pointer , dimension(:,:,:) :: up , vp
  real(rkx) , pointer , dimension(:,:,:) :: t3 , q3 , h3
  real(rkx) , pointer , dimension(:,:,:) :: u3 , v3

  real(rkx) , dimension(nplev) :: pplev , sigmar
  real(rkx) :: pss

  integer(ik4) :: mdlver , ensnum
  integer(ik4) , parameter :: ibctime = 6
  integer(ik4) , dimension(5) :: inet5, ivar5
  integer(ik4) , dimension(4) :: icount , istart
  type(rcm_time_and_date) , save :: ilastdate
  type(rcm_time_and_date) , pointer , dimension(:) :: enstime
  real(rkx) , pointer , dimension(:) :: xtime
  character(len=256) :: pathaddname , fname
  character(len=4) :: ensbase = 'f8qe'
  type(rcm_time_and_date) , save :: fmon
  integer(ik4) :: ifmon

  public :: getecens , headerecens

  contains

  subroutine getecens(idate)
    implicit none
    type(rcm_time_and_date) :: idate
    integer(ik4) :: i , j , k

    call ecens_6hour(idate)

    write (*,*) 'READ IN fields at DATE:' , tochar(idate)

    do k = 1 , mlev
      do j = 1 , nlat
        do i = 1 , nlon
          pp3d(i,j,k) = ak(k) + ps2(i,j)*bk(k)
        end do
      end do
    end do

    ps2  = ps2 * 0.01   ! Get to millibars
    pp3d = pp3d * 0.01  ! Get to millibars
    !
    ! to calculate Heights on sigma surfaces.
    call htsig(t2,z1,pp3d,ps2,zs2,nlon,nlat,mlev)
    !
    ! to interpolate H,U,V,T,Q and QC
    ! 1. For Heights
    call height(hp,z1,t2,ps2,pp3d,zs2,nlon,nlat,mlev,pplev,nplev)
    ! 2. For Zonal and Meridional Winds
    call intlin(up,u2,ps2,pp3d,nlon,nlat,mlev,pplev,nplev)
    call intlin(vp,v2,ps2,pp3d,nlon,nlat,mlev,pplev,nplev)
    ! 3. For Temperatures
    call intlog(tp,t2,ps2,pp3d,nlon,nlat,mlev,pplev,nplev)
    ! 4. For Moisture qva & qca
    call mxr2rh(t2,q2,pp3d,nlon,nlat,mlev,-9999.0_rkx)
    call intlin(qp,q2,ps2,pp3d,nlon,nlat,mlev,pplev,nplev)
    !
    ! Horizontal interpolation of both the scalar and vector fields
    !
    call bilinx2(b3,b2,xlon,xlat,vlon,vlat,nlon,nlat,jx,iy,nplev*3)
    call bilinx2(d3,d2,dlon,dlat,vlon,vlat,nlon,nlat,jx,iy,nplev*2)
    !
    ! Rotate U-V fields after horizontal interpolation
    !
    call uvrot4(u3,v3,dlon,dlat,clon,clat,xcone,jx,iy,nplev,plon,plat,iproj)
    !
    ! Vertical interpolation
    !
    call top2btm(t3,jx,iy,nplev)
    call top2btm(q3,jx,iy,nplev)
    call top2btm(h3,jx,iy,nplev)
    call top2btm(u3,jx,iy,nplev)
    call top2btm(v3,jx,iy,nplev)
    !
    ! New calculation of P* on RegCM topography.
    !
    call intgtb(pa,za,tlayer,topogm,t3,h3,pss,sigmar,jx,iy,nplev)
    call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
    call crs2dot(pd4,ps4,jx,iy,i_band)
    !
    ! Determine surface temps on RegCM topography.
    ! interpolation from pressure levels
    !
    call intv3(ts4,t3,ps4,pss,sigmar,ptop,jx,iy,nplev)
    call readsst(ts4,idate)
    ! Interpolate U, V, T, and Q.
    call intv1(u4,u3,pd4,sigmah,pss,sigmar,ptop,jx,iy,kz,nplev)
    call intv1(v4,v3,pd4,sigmah,pss,sigmar,ptop,jx,iy,kz,nplev)
    call intv2(t4,t3,ps4,sigmah,pss,sigmar,ptop,jx,iy,kz,nplev)
    call intv1(q4,q3,ps4,sigmah,pss,sigmar,ptop,jx,iy,kz,nplev)
    call rh2mxr(t4,q4,ps4,ptop,sigmah,jx,iy,kz)
  end subroutine getecens
  !
  !-------------------------------------------------------------
  !
  subroutine ecens_6hour(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate

    type(rcm_time_and_date) :: fmon
    type(rcm_time_interval) :: tdif
    integer(ik4) :: i , inet , it , j , kkrec , istat , ivar , jdim
    logical :: lfirst
    character(len=6) , dimension(5) :: varname
    character(len=6) , dimension(5) :: vfname
    character(len=64) :: cunit , ccal

    data varname /'var152','var130','var133','var131','var132'/
    data vfname  /'lnsp','air','qhum','uwnd','vwnd'/
    data lfirst  /.true./

    ! Open file and get files id and variable id(s)

    fmon = monfirst(idate)

    ! ECMWF Fix : they do not have calendar attribute to time variable
    ccal = 'gregorian'

    if ( lfirst .or. (idate > ilastdate) ) then
      do kkrec = 1 , 5
        fmon = monfirst(idate)
        ifmon = toint10(fmon)/100
        write (fname,'(a,i1,a,i0.8,a,a,a,a,a,i1,a,i0.8,a,i0.2,a)') &
           'ENS', mdlver, pthsep, ifmon , pthsep, trim(ensbase), '.', &
           trim(vfname(kkrec)), '.', ensnum, '.', ifmon, '.', 0, '.nc'
        pathaddname = trim(inpglob)//pthsep//trim(fname)
        istat=nf90_open(pathaddname,nf90_nowrite,inet5(kkrec))
        call checkncerr(istat,__FILE__,__LINE__, &
                        'Error open '//trim(pathaddname))
        istat = nf90_inq_varid(inet5(kkrec),varname(kkrec),ivar5(kkrec))
        call checkncerr(istat,__FILE__,__LINE__, &
                        'Error find var '//trim(varname(kkrec)))
        if ( kkrec == 1 ) then
          istat = nf90_inq_dimid(inet5(kkrec),'time',jdim)
          call checkncerr(istat,__FILE__,__LINE__, &
                          'Error find dim time')
          istat = nf90_inquire_dimension(inet5(kkrec),jdim,len=ntime)
          call checkncerr(istat,__FILE__,__LINE__, &
                          'Error read dim time')
          call getmem1d(xtime,1,ntime,'mod_ecens:xtime')
          call getmem1d(enstime,1,ntime,'mod_ecens:enstime')
          istat = nf90_inq_varid(inet5(kkrec),'time',ivar)
          call checkncerr(istat,__FILE__,__LINE__, &
                          'Error find var time')
          istat = nf90_get_att(inet5(kkrec),ivar,'units',cunit)
          call checkncerr(istat,__FILE__,__LINE__, &
                          'Error read var units')
          istat = nf90_get_var(inet5(kkrec),ivar,xtime)
          call checkncerr(istat,__FILE__,__LINE__, &
                          'Error read var time')
          do i = 1 , ntime
            enstime(i) = timeval2date(xtime(i),cunit,ccal)
          end do
          ilastdate = enstime(ntime)
        end if
        write (*,*) inet5(kkrec), trim(pathaddname)
      enddo
      lfirst = .false.
    endif

    tdif = idate - enstime(1)
    it = nint(tohours(tdif))/ibctime+1

    do kkrec = 1 , 5

      inet = inet5(kkrec)
      ivar = ivar5(kkrec)

      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = it

      icount(1) = nlon
      icount(2) = nlat
      icount(3) = mlev
      icount(4) = 1

      if ( kkrec == 1 ) then
        icount(3) = 1
        istat = nf90_get_var(inet,ivar,work(:,:,1),istart,icount)
        call checkncerr(istat,__FILE__,__LINE__, &
                        'Err read var '//vfname(kkrec))
      else
        istat = nf90_get_var(inet,ivar,work,istart,icount)
        call checkncerr(istat,__FILE__,__LINE__, &
                        'Err read var '//vfname(kkrec))
      endif

      if ( kkrec == 1 ) then
        do j = 1 , nlat
          ps2(:,j) = exp(work(:,nlat-j+1,1))
        end do
      else if ( kkrec==2 ) then
        do j = 1 , nlat
          t2(:,j,:) = work(:,nlat-j+1,:)
        end do
      else if ( kkrec==3 ) then
        do j = 1 , nlat
          q2(:,j,:) = work(:,nlat-j+1,:)
        end do
      else if ( kkrec==4 ) then
        do j = 1 , nlat
          u2(:,j,:) = work(:,nlat-j+1,:)
        end do
      else
        do j = 1 , nlat
          v2(:,j,:) = work(:,nlat-j+1,:)
        end do
      end if
    end do
  end subroutine ecens_6hour

  subroutine headerecens
    implicit none
    integer(ik4) :: j , k
    integer(ik4) :: istat , inet , jdim , ivar

    read (dattyp(4:4),'(i1)') mdlver
    read (dattyp(5:5),'(i1)') ensnum
    fmon = monfirst(globidate1)
    ifmon = toint10(fmon)/100
    write (fname,'(a,i1,a,i0.8,a,a,a,a,a,i1,a,i0.8,a,i0.2,a)') &
    'ENS', mdlver, pthsep, ifmon , pthsep, trim(ensbase), '.', 'hgt', &
                  '.', ensnum, '.', ifmon, '.', 0, '.nc'
    pathaddname = trim(inpglob)//pthsep//trim(fname)

    istat = nf90_open(pathaddname,nf90_nowrite,inet)
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error open '//trim(pathaddname))
    istat = nf90_inq_dimid(inet,'lon',jdim)
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error find dim lon')
    istat = nf90_inquire_dimension(inet,jdim,len=nlon)
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error read dim lon')
    istat = nf90_inq_dimid(inet,'lat',jdim)
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error find dim lat')
    istat = nf90_inquire_dimension(inet,jdim,len=nlat)
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error read dim lat')
    istat = nf90_inq_dimid(inet,'mlev',jdim)
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error find dim mlev')
    istat = nf90_inquire_dimension(inet,jdim,len=mlev)
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error read dim mlev')

    call getmem1d(vlat,1,nlat,'mod_ecens:vlat')
    call getmem1d(rlat,1,nlat,'mod_ecens:rlat')
    call getmem1d(vlon,1,nlon,'mod_ecens:vlon')
    call getmem1d(ak,1,mlev,'mod_ecens:ak')
    call getmem1d(bk,1,mlev,'mod_ecens:bk')
    call getmem3d(work,1,nlon,1,nlat,1,mlev,'mod_ecens:work')

    istat = nf90_inq_varid(inet,'lon',ivar)
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error find var lon')
    istat = nf90_get_var(inet,ivar,vlon)
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error read var lon')
    istat = nf90_inq_varid(inet,'lat',ivar)
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error find var lat')
    istat = nf90_get_var(inet,ivar,rlat)
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error read var lat')
    do j = 1 , nlat
      vlat(j) = rlat(nlat-j+1)
    end do
    istat = nf90_inq_varid(inet,'hyam',ivar)
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error find var hyam')
    istat = nf90_get_var(inet,ivar,ak)
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error read var hyam')
    istat = nf90_inq_varid(inet,'hybm',ivar)
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error find var hyam')
    istat = nf90_get_var(inet,ivar,bk)
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error read var hybm')

    pplev(1) = 30.
    pplev(2) = 50.
    pplev(3) = 70.
    pplev(4) = 100.
    pplev(5) = 150.
    pplev(6) = 200.
    pplev(7) = 250.
    pplev(8) = 300.
    pplev(9) = 350.
    pplev(10) = 420.
    pplev(11) = 500.
    pplev(12) = 600.
    pplev(13) = 700.
    pplev(14) = 780.
    pplev(15) = 850.
    pplev(16) = 920.
    pplev(17) = 960.
    pplev(18) = 1000.

    do k = 1 , nplev
      sigmar(k) = pplev(k)/pplev(nplev)
    end do
    pss = pplev(nplev) / d_10 ! centibars

    call getmem3d(bb,1,nlon,1,nlat,1,mlev*4+2,'mod_ecens:bb')
    call getmem3d(pp3d,1,nlon,1,nlat,1,mlev,'mod_ecens:pp3d')
    call getmem3d(z1,1,nlon,1,nlat,1,mlev,'mod_ecens:z1')
    call getmem3d(b2,1,nlon,1,nlat,1,nplev*3,'mod_ecens:b2')
    call getmem3d(d2,1,nlon,1,nlat,1,nplev*2,'mod_ecens:d2')
    call getmem3d(b3,1,jx,1,iy,1,nplev*3,'mod_ecens:b3')
    call getmem3d(d3,1,jx,1,iy,1,nplev*2,'mod_ecens:d3')

    ! Set up pointers

    t2 => bb(:,:,1:mlev)
    q2 => bb(:,:,mlev+1:2*mlev)
    u2 => bb(:,:,2*mlev+1:3*mlev)
    v2 => bb(:,:,3*mlev+1:4*mlev)
    zs2 => bb(:,:,4*mlev+1)
    ps2 => bb(:,:,4*mlev+2)

    tp => b2(:,:,1:nplev)
    qp => b2(:,:,nplev+1:2*nplev)
    hp => b2(:,:,2*nplev+1:3*nplev)
    up => d2(:,:,1:nplev)
    vp => d2(:,:,nplev+1:2*nplev)
    t3 => b3(:,:,1:nplev)
    q3 => b3(:,:,nplev+1:2*nplev)
    h3 => b3(:,:,2*nplev+1:3*nplev)
    u3 => d3(:,:,1:nplev)
    v3 => d3(:,:,nplev+1:2*nplev)

    istat = nf90_inq_varid(inet,'var129',ivar)
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error find var var129')

    istart(:) = 1
    icount(1) = nlon
    icount(2) = nlat
    icount(3) = 1
    icount(4) = 1
    istat = nf90_get_var(inet,ivar,work,istart,icount)
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error read var var129')
    do j = 1 , nlat
      zs2(:,j) = work(:,nlat-j+1,1)
    end do
    where ( zs2 > 0.0 )
      zs2 = zs2 / 9.80616
    elsewhere
      zs2  = 0.0
    end where

    istat = nf90_close(inet)
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error close file')
  end subroutine headerecens

end module mod_ecens

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
