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

      module mod_nest

      use mod_dynparam
      use mod_constants
      use mod_date
      use mod_grid
      use mod_write
      use mod_interp
      use mod_vertint
      use mod_hgt
      use mod_humid
      use mod_mksst
      use mod_uvrot
      use mod_vectutil

      private

      integer , parameter :: np = 15

      integer :: nrec

      real(4) , allocatable , target , dimension(:,:,:) :: b3
      real(4) , allocatable , target , dimension(:,:,:) :: d3
      real(4) , allocatable , dimension(:,:) :: xb3pd
      real(4) , allocatable , dimension(:,:,:) :: z1

      real(4) , allocatable , target , dimension(:,:,:) :: b2
      real(4) , allocatable , target , dimension(:,:,:) :: d2

      real(4) , allocatable , dimension(:,:,:) :: c , q , t
      real(4) , allocatable , dimension(:,:,:) :: u , v
      real(4) , allocatable , dimension(:,:) :: ps
      real(4) , allocatable , dimension(:,:) :: ht_in
      real(4) , allocatable , dimension(:,:) :: xlat_in , xlon_in

      real(4) , pointer , dimension(:,:,:) :: c3 , h3 , q3 , t3
      real(4) , pointer , dimension(:,:,:) :: u3 , v3
 
      real(4) , pointer , dimension(:,:,:) :: cp , hp , qp , tp
      real(4) , pointer , dimension(:,:,:) :: up , vp

      real(4) , dimension(np) :: plev , sigmar
      real(4) , allocatable , dimension(:) :: sig

      character(6) :: iproj_in

      integer :: iy_in , jx_in , kz_in , iotyp_in , idate0

      real(4) :: clat_in , clon_in , plat_in , plon_in , ptop_in

      public :: get_nest , headnest

      character(14) :: fillin
      character(256) :: inpfile
!
      integer :: ncid
      integer , dimension(:) , allocatable :: itimes
      character(64) :: timeunits
!
      contains

      subroutine get_nest(idate)
      use netcdf
      implicit none
!
      integer , intent(in) :: idate
!
      real(8) , dimension(:) , allocatable :: xtimes
      integer :: i , istatus , ivarid , idimid , irec
      integer , dimension(4) :: istart , icount
!
      if (.not. allocated(b2)) then
        write (*,*) 'Called get_nest before headnest !'
        stop
      end if
!
      if ( idate > itimes(nrec) ) then
        istatus = nf90_close(ncid)
        call check_ok(istatus, 'Error close')
        write (fillin,'(a,i10)') 'ATM.', idate
        inpfile = trim(inpglob)//pthsep//'RegCM'//pthsep//fillin//'.nc'
        istatus = nf90_open(inpfile, nf90_nowrite, ncid)
        call check_ok(istatus, 'Error opening '//trim(inpfile))
        istatus = nf90_inq_dimid(ncid, 'time', idimid)
        call check_ok(istatus,'Dimension time missing')
        istatus = nf90_inquire_dimension(ncid, idimid, len=nrec)
        call check_ok(istatus,'Dimension time read error')
        istatus = nf90_inq_varid(ncid, 'time', ivarid)
        call check_ok(istatus,'variable time missing')
        deallocate(itimes)
        allocate(itimes(nrec))
        allocate(xtimes(nrec))
        istatus = nf90_get_var(ncid, ivarid, xtimes)
        call check_ok(istatus,'variable time read error')
        do i = 1 , nrec
          itimes(i) = timeval2idate(xtimes(i), timeunits)
        end do
        deallocate(xtimes)
      end if

      irec = -1
      do irec = 1 , nrec
        if (idate == itimes(irec)) exit
      end do
      if (irec < 0) then
        write (6,*) 'Error : time ', idate, ' not in file'
        stop
      end if

      istart(4) = irec
      istart(3) = 1
      istart(2) = 1
      istart(1) = 1
      icount(4) = 1
      icount(3) = kz_in
      icount(2) = iy_in
      icount(1) = jx_in
      istatus = nf90_inq_varid(ncid, 'u', ivarid)
      call check_ok(istatus,'variable u missing')
      istatus = nf90_get_var(ncid, ivarid, u, istart, icount)
      call check_ok(istatus,'variable u read error')
      istatus = nf90_inq_varid(ncid, 'v', ivarid)
      call check_ok(istatus,'variable v missing')
      istatus = nf90_get_var(ncid, ivarid, v, istart, icount)
      call check_ok(istatus,'variable v read error')
      istatus = nf90_inq_varid(ncid, 't', ivarid)
      call check_ok(istatus,'variable t missing')
      istatus = nf90_get_var(ncid, ivarid, t, istart, icount)
      call check_ok(istatus,'variable t read error')
      istatus = nf90_inq_varid(ncid, 'qv', ivarid)
      call check_ok(istatus,'variable qv missing')
      istatus = nf90_get_var(ncid, ivarid, q, istart, icount)
      call check_ok(istatus,'variable qv read error')
      istatus = nf90_inq_varid(ncid, 'qc', ivarid)
      call check_ok(istatus,'variable qc missing')
      istatus = nf90_get_var(ncid, ivarid, c, istart, icount)
      call check_ok(istatus,'variable qc read error')
      istart(3) = irec
      istart(2) = 1
      istart(1) = 1
      icount(3) = 1
      icount(2) = iy_in
      icount(1) = jx_in
      istatus = nf90_inq_varid(ncid, 'ps', ivarid)
      call check_ok(istatus,'variable ps missing')
      istatus = nf90_get_var(ncid, ivarid, ps, istart(1:3), icount(1:3))
      call check_ok(istatus,'variable ps read error')

      write (*,*) 'READ IN fields at DATE:' , idate , ' from ' , fillin

!     to calculate Heights on sigma surfaces.
      call htsig_o(t,z1,ps,ht_in,sig,ptop_in,jx_in,iy_in,kz_in)
!
!     to interpolate H,U,V,T,Q and QC
!     1. For Heights
      call height_o(hp,z1,t,ps,ht_in,sig,ptop_in,jx_in,iy_in,kz_in,    &
                  & plev,np)
!     2. For Zonal and Meridional Winds
      call intlin_o(up,u,ps,sig,ptop_in,jx_in,iy_in,kz_in,plev,np)
      call intlin_o(vp,v,ps,sig,ptop_in,jx_in,iy_in,kz_in,plev,np)
!     3. For Temperatures
      call intlog_o(tp,t,ps,sig,ptop_in,jx_in,iy_in,kz_in,plev,np)
!     4. For Moisture qva & qca
      call humid1_o(t,q,ps,sig,ptop_in,jx_in,iy_in,kz_in)
      call intlin_o(qp,q,ps,sig,ptop_in,jx_in,iy_in,kz_in,plev,np)
      call intlog_o(cp,c,ps,sig,ptop_in,jx_in,iy_in,kz_in,plev,np)
      call uvrot4nx(up,vp,xlon_in,xlat_in,clon_in,clat_in,grdfac,       &
             &      jx_in,iy_in,np,plon_in,plat_in,iproj_in)
!
!     HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
      call cressmcr(b3,b2,xlon,xlat,xlon_in,xlat_in,jx,iy,              &
                  & jx_in,iy_in,np,4)
      call cressmdt(d3,d2,dlon,dlat,xlon_in,xlat_in,jx,iy,              &
                  & jx_in,iy_in,np,2)
!
!     ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
      call uvrot4(u3,v3,dlon,dlat,clon,clat,grdfac,jx,iy,np,plon,plat,  &
                & iproj)
!
!     X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!     X X
!     V E R T I C A L   I N T E R P O L A T I O N
!
!     X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!     X X
!HH:  CHANGE THE VERTICAL ORDER.
      call top2btm(t3,jx,iy,np)
      call top2btm(q3,jx,iy,np)
      call top2btm(c3,jx,iy,np)
      call top2btm(h3,jx,iy,np)
      call top2btm(u3,jx,iy,np)
      call top2btm(v3,jx,iy,np)
!HH:OVER
!
!     ******           NEW CALCULATION OF P* ON RegCM TOPOGRAPHY.
      call intgtb(pa,za,tlayer,topogm,t3,h3,sigmar,jx,iy,np)
 
      call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
      if(i_band.eq.1) then
         call p1p2_band(xb3pd,ps4,jx,iy)
      else
         call p1p2(xb3pd,ps4,jx,iy)
      endif
!
!     F0    DETERMINE SURFACE TEMPS ON RegCM TOPOGRAPHY.
!     INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
      call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,np)
 
      call readsst(ts4,idate)

!     F2     DETERMINE P* AND HEIGHT.
!
!     F3     INTERPOLATE U, V, T, AND Q.
      call intv1(u4,u3,xb3pd,sigma2,sigmar,ptop,jx,iy,kz,np)
      call intv1(v4,v3,xb3pd,sigma2,sigmar,ptop,jx,iy,kz,np)
!
      call intv2(t4,t3,ps4,sigma2,sigmar,ptop,jx,iy,kz,np)
 
      call intv1(q4,q3,ps4,sigma2,sigmar,ptop,jx,iy,kz,np)
      call humid2fv(t4,q4,ps4,ptop,sigma2,jx,iy,kz)
      call intv1(c4,c3,ps4,sigma2,sigmar,ptop,jx,iy,kz,np)
!
!     F4     DETERMINE H
      call hydrost(h4,t4,topogm,ps4,ptop,sigmaf,sigma2,dsigma,jx,iy,kz)
!
!     G      WRITE AN INITIAL FILE FOR THE RegCM
      call writef(idate)
!
      end subroutine get_nest
!
!
!
      subroutine headnest
      use netcdf
      implicit none
!
! Local variables
!
      real(4) :: xsign
      integer :: i , k , istatus , idimid , ivarid
      logical :: there
      real(8) , dimension(:) , allocatable :: xtimes
      real(4) , dimension(2) :: trlat
!
! Setup interp
!
      imxmn = 0
      lcross = 0
      ldot = 0
!
      plev(1) = 50.
      plev(2) = 70.
      plev(3) = 100.
      plev(4) = 150.
      plev(5) = 200.
      plev(6) = 250.
      plev(7) = 300.
      plev(8) = 400.
      plev(9) = 500.
      plev(10) = 600.
      plev(11) = 700.
      plev(12) = 775.
      plev(13) = 850.
      plev(14) = 925.
      plev(15) = 1000.
 
      do k = 1 , np
        sigmar(k) = plev(k)*0.001
      end do
 
      write (fillin,'(a,i10)') 'ATM.', imonfirst(globidate1)
      inpfile = trim(inpglob)//pthsep//'RegCM'//pthsep//fillin//'.nc'
      inquire (file=inpfile,exist=there)
      if ( .not.there ) then
        write (*,*) trim(inpfile), ' is not available'
        write (*,*) 'please copy (or link)' , trim(inpfile)
        stop
      end if
      istatus = nf90_open(inpfile, nf90_nowrite, ncid)
      call check_ok(istatus, 'Error opening '//trim(inpfile))

      istatus = nf90_inq_dimid(ncid, 'iy', idimid)
      call check_ok(istatus,'Dimension iy missing')
      istatus = nf90_inquire_dimension(ncid, idimid, len=iy_in)
      call check_ok(istatus,'Dimension iy read error')
      istatus = nf90_inq_dimid(ncid, 'jx', idimid)
      call check_ok(istatus,'Dimension jx missing')
      istatus = nf90_inquire_dimension(ncid, idimid, len=jx_in)
      call check_ok(istatus,'Dimension jx read error')
      istatus = nf90_inq_dimid(ncid, 'kz', idimid)
      call check_ok(istatus,'Dimension kz missing')
      istatus = nf90_inquire_dimension(ncid, idimid, len=kz_in)
      call check_ok(istatus,'Dimension kz read error')
      istatus = nf90_inq_dimid(ncid, 'time', idimid)
      call check_ok(istatus,'Dimension time missing')
      istatus = nf90_inquire_dimension(ncid, idimid, len=nrec)
      call check_ok(istatus,'Dimension time read error')
      istatus = nf90_inq_varid(ncid, 'time', ivarid)
      call check_ok(istatus,'variable time missing')
      istatus = nf90_get_att(ncid, ivarid, 'units', timeunits)
      call check_ok(istatus,'variable time units missing')
      allocate(itimes(nrec))
      allocate(xtimes(nrec))
      istatus = nf90_get_var(ncid, ivarid, xtimes)
      call check_ok(istatus,'variable time read error')
      do i = 1 , nrec
        itimes(i) = timeval2idate(xtimes(i), timeunits)
      end do
      deallocate(xtimes)

!     Reserve space for I/O

      allocate(sig(kz_in), stat=istatus)
      if (istatus /= 0) stop 'Allocation Error in headnest: sig'
      allocate(b2(jx_in,iy_in,np*4), stat=istatus)
      if (istatus /= 0) stop 'Allocation Error in headnest: b2'
      allocate(d2(jx_in,iy_in,np*2), stat=istatus)
      if (istatus /= 0) stop 'Allocation Error in headnest: d2'
      allocate(c(jx_in,iy_in,kz_in), stat=istatus)
      if (istatus /= 0) stop 'Allocation Error in headnest: c'
      allocate(q(jx_in,iy_in,kz_in), stat=istatus)
      if (istatus /= 0) stop 'Allocation Error in headnest: q'
      allocate(t(jx_in,iy_in,kz_in), stat=istatus)
      if (istatus /= 0) stop 'Allocation Error in headnest: t'
      allocate(u(jx_in,iy_in,kz_in), stat=istatus)
      if (istatus /= 0) stop 'Allocation Error in headnest: u'
      allocate(v(jx_in,iy_in,kz_in), stat=istatus)
      if (istatus /= 0) stop 'Allocation Error in headnest: v'
      allocate(ps(jx_in,iy_in), stat=istatus)
      if (istatus /= 0) stop 'Allocation Error in headnest: ps'
      allocate(xlat_in(jx_in,iy_in), stat=istatus)
      if (istatus /= 0) stop 'Allocation Error in headnest: xlat_in'
      allocate(xlon_in(jx_in,iy_in), stat=istatus)
      if (istatus /= 0) stop 'Allocation Error in headnest: xlon_in'
      allocate(ht_in(jx_in,iy_in), stat=istatus)
      if (istatus /= 0) stop 'Allocation Error in headnest: ht_in'

      istatus = nf90_inq_varid(ncid, 'sigma', ivarid) 
      call check_ok(istatus,'variable sigma error')
      istatus = nf90_get_var(ncid, ivarid, sig)
      call check_ok(istatus,'variable sigma read error')
      istatus = nf90_inq_varid(ncid, 'xlat', ivarid) 
      call check_ok(istatus,'variable xlat error')
      istatus = nf90_get_var(ncid, ivarid, xlat_in)
      call check_ok(istatus,'variable xlat read error')
      istatus = nf90_inq_varid(ncid, 'xlon', ivarid) 
      call check_ok(istatus,'variable xlon error')
      istatus = nf90_get_var(ncid, ivarid, xlon_in)
      call check_ok(istatus,'variable xlon read error')
      istatus = nf90_inq_varid(ncid, 'topo', ivarid) 
      call check_ok(istatus,'variable topo error')
      istatus = nf90_get_var(ncid, ivarid, ht_in)
      call check_ok(istatus,'variable topo read error')
      istatus = nf90_inq_varid(ncid, 'ptop', ivarid) 
      call check_ok(istatus,'variable ptop error')
      istatus = nf90_get_var(ncid, ivarid, ptop_in)
      call check_ok(istatus,'variable ptop read error')

      istatus = nf90_get_att(ncid, nf90_global, &
                   &    'projection', iproj_in)
      call check_ok(istatus,'attribure iproj read error')
      istatus = nf90_get_att(ncid, nf90_global, &
                   &    'latitude_of_projection_origin', clat_in)
      call check_ok(istatus,'attribure clat read error')
      istatus = nf90_get_att(ncid, nf90_global, &
                   &    'longitude_of_projection_origin', clon_in)
      call check_ok(istatus,'attribure clat read error')

      if ( iproj_in=='LAMCON' ) then
        istatus = nf90_get_att(ncid, nf90_global, &
                   &    'standard_parallel', trlat)
        call check_ok(istatus,'attribure truelat read error')
        if ( clat_in<0. ) then
          xsign = -1.       ! SOUTH HEMESPHERE
        else
          xsign = 1.        ! NORTH HEMESPHERE
        end if
        if ( abs(trlat(1)-trlat(2))>1.E-1 ) then
          grdfac = (log10(cos(trlat(1)*degrad))                         &
                   & -log10(cos(trlat(2)*degrad)))                      &
                   & /(log10(tan((45.0-xsign*trlat(1)/2.0)*degrad))     &
                   & -log10(tan((45.0-xsign*trlat(2)/2.0)*degrad)))
        else
          grdfac = xsign*sin(trlat(1)*degrad)
        end if
      else if ( iproj_in=='POLSTR' ) then
        grdfac = 1.0
      else if ( iproj_in=='NORMER' ) then
        grdfac = 0.0
      else
        istatus = nf90_get_att(ncid, nf90_global, &
                   &    'grid_north_pole_latitude', plat_in)
        call check_ok(istatus,'attribure plat read error')
        istatus = nf90_get_att(ncid, nf90_global, &
                   &    'grid_north_pole_longitude', plon_in)
        call check_ok(istatus,'attribure plon read error')
        grdfac = 0.0
      end if
 
      if (allocated(b3)) deallocate(b3)
      if (allocated(d3)) deallocate(d3)
      if (allocated(xb3pd)) deallocate(xb3pd)
      if (allocated(z1)) deallocate(z1)
      allocate(b3(iy_in,jx_in,np*4))
      allocate(d3(iy_in,jx_in,np*2))
      allocate(xb3pd(iy_in,jx_in))
      allocate(z1(iy_in,jx_in,kz_in))

!     Set up pointers
 
      tp => b2(:,:,1:np)
      qp => b2(:,:,np+1:2*np)
      cp => b2(:,:,2*np+1:3*np)
      hp => b2(:,:,3*np+1:4*np)
      up => d2(:,:,1:np)
      vp => d2(:,:,np+1:2*np)
      t3 => b3(:,:,1:np)
      q3 => b3(:,:,np+1:2*np)
      c3 => b3(:,:,2*np+1:3*np)
      h3 => b3(:,:,3*np+1:4*np)
      u3 => d3(:,:,1:np)
      v3 => d3(:,:,np+1:2*np)

      end subroutine headnest

      end module mod_nest
