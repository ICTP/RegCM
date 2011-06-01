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
      real(4) , allocatable , dimension(:,:,:) :: z1

      real(4) , allocatable , target , dimension(:,:,:) :: b2
      real(4) , allocatable , target , dimension(:,:,:) :: d2

      real(4) , allocatable , dimension(:,:,:) :: q , t
      real(4) , allocatable , dimension(:,:,:) :: u , v
      real(4) , allocatable , dimension(:,:) :: ps
      real(4) , allocatable , dimension(:,:) :: ht_in
      real(4) , allocatable , dimension(:,:) :: xlat_in , xlon_in

      real(4) , pointer , dimension(:,:,:) :: h3 , q3 , t3
      real(4) , pointer , dimension(:,:,:) :: u3 , v3
 
      real(4) , pointer , dimension(:,:,:) :: hp , qp , tp
      real(4) , pointer , dimension(:,:,:) :: up , vp

      real(4) , dimension(np) :: plev , sigmar
      real(4) , allocatable , dimension(:) :: sig

      character(6) :: iproj_in

      integer :: iy_in , jx_in , kz_in , iotyp_in , idate0

      real(8) :: clat_in , clon_in , plat_in , plon_in , ptop_in , &
                 grdfac_in

      public :: get_nest , headnest

      character(14) :: fillin
      character(256) :: inpfile
!
      integer :: ncinp
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
        istatus = nf90_close(ncinp)
        call check_ok(istatus, 'Error close')
        write (fillin,'(a,i10)') 'ATM.', imonfirst(idate)
        inpfile = trim(inpglob)//pthsep//'RegCM'//pthsep//fillin//'.nc'
        istatus = nf90_open(inpfile, nf90_nowrite, ncinp)
        call check_ok(istatus, 'Error opening '//trim(inpfile))
        istatus = nf90_inq_dimid(ncinp, 'time', idimid)
        call check_ok(istatus,'Dimension time missing')
        istatus = nf90_inquire_dimension(ncinp, idimid, len=nrec)
        call check_ok(istatus,'Dimension time read error')
        istatus = nf90_inq_varid(ncinp, 'time', ivarid)
        call check_ok(istatus,'variable time missing')
        deallocate(itimes)
        allocate(itimes(nrec))
        allocate(xtimes(nrec))
        istatus = nf90_get_var(ncinp, ivarid, xtimes)
        call check_ok(istatus,'variable time read error')
        do i = 1 , nrec
          itimes(i) = timeval2idate(xtimes(i), timeunits)
        end do
        deallocate(xtimes)
      end if

      irec = -1
      do i = 1 , nrec
        if (idate == itimes(i)) then
          irec = i
          exit
        end if
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
      istatus = nf90_inq_varid(ncinp, 'u', ivarid)
      call check_ok(istatus,'variable u missing')
      istatus = nf90_get_var(ncinp, ivarid, u, istart, icount)
      call check_ok(istatus,'variable u read error')
      istatus = nf90_inq_varid(ncinp, 'v', ivarid)
      call check_ok(istatus,'variable v missing')
      istatus = nf90_get_var(ncinp, ivarid, v, istart, icount)
      call check_ok(istatus,'variable v read error')
      istatus = nf90_inq_varid(ncinp, 't', ivarid)
      call check_ok(istatus,'variable t missing')
      istatus = nf90_get_var(ncinp, ivarid, t, istart, icount)
      call check_ok(istatus,'variable t read error')
      istatus = nf90_inq_varid(ncinp, 'qv', ivarid)
      call check_ok(istatus,'variable qv missing')
      istatus = nf90_get_var(ncinp, ivarid, q, istart, icount)
      call check_ok(istatus,'variable qv read error')
      istart(3) = irec
      istart(2) = 1
      istart(1) = 1
      icount(3) = 1
      icount(2) = iy_in
      icount(1) = jx_in
      istatus = nf90_inq_varid(ncinp, 'ps', ivarid)
      call check_ok(istatus,'variable ps missing')
      istatus = nf90_get_var(ncinp, ivarid, ps, istart(1:3), icount(1:3))
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
      call uvrot4nx(up,vp,xlon_in,xlat_in,clon_in,clat_in,grdfac_in,    &
             &      jx_in,iy_in,np,plon_in,plat_in,iproj_in)
!
!     HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
      call cressmcr(b3,b2,xlon,xlat,xlon_in,xlat_in,jx,iy,              &
                  & jx_in,iy_in,np,3)
      call cressmdt(d3,d2,dlon,dlat,xlon_in,xlat_in,jx,iy,              &
                  & jx_in,iy_in,np,2)
!
!     ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
      call uvrot4(u3,v3,dlon,dlat,clon,clat,grdfac,jx,iy,np, &
                  plon,plat,iproj)
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
      call top2btm(h3,jx,iy,np)
      call top2btm(u3,jx,iy,np)
      call top2btm(v3,jx,iy,np)
!HH:OVER
!
!     NEW CALCULATION OF P* ON RegCM TOPOGRAPHY.
      call intgtb(pa,za,tlayer,topogm,t3,h3,sigmar,jx,iy,np)
 
      call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
      if(i_band.eq.1) then
         call p1p2_band(b3pd,ps4,jx,iy)
      else
         call p1p2(b3pd,ps4,jx,iy)
      endif
!
!     F0    DETERMINE SURFACE TEMPS ON RegCM TOPOGRAPHY.
!     INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
      call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,np)
 
      call readsst(ts4,idate)

!     F2     DETERMINE P* AND HEIGHT.
!
!     F3     INTERPOLATE U, V, T, AND Q.
      call intv1(u4,u3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,np)
      call intv1(v4,v3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,np)
!
      call intv2(t4,t3,ps4,sigma2,sigmar,ptop,jx,iy,kz,np)
 
      call intv1(q4,q3,ps4,sigma2,sigmar,ptop,jx,iy,kz,np)
      call humid2(t4,q4,ps4,ptop,sigma2,jx,iy,kz)
!
!     F4     DETERMINE H
      call hydrost(h4,t4,topogm,ps4,ptop,sigmaf,sigma2,dsigma,jx,iy,kz)
!
!     G      WRITE AN INITIAL FILE FOR THE RegCM
      call writef(idate)
!
      end subroutine get_nest
!
      subroutine headnest
      use netcdf
      implicit none
!
      real(8) :: xsign
      integer :: i , k , istatus , idimid , ivarid
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
      istatus = nf90_open(inpfile, nf90_nowrite, ncinp)
      call check_ok(istatus, 'Error opening '//trim(inpfile))

      istatus = nf90_inq_dimid(ncinp, 'iy', idimid)
      call check_ok(istatus,'Dimension iy missing')
      istatus = nf90_inquire_dimension(ncinp, idimid, len=iy_in)
      call check_ok(istatus,'Dimension iy read error')
      istatus = nf90_inq_dimid(ncinp, 'jx', idimid)
      call check_ok(istatus,'Dimension jx missing')
      istatus = nf90_inquire_dimension(ncinp, idimid, len=jx_in)
      call check_ok(istatus,'Dimension jx read error')
      istatus = nf90_inq_dimid(ncinp, 'kz', idimid)
      call check_ok(istatus,'Dimension kz missing')
      istatus = nf90_inquire_dimension(ncinp, idimid, len=kz_in)
      call check_ok(istatus,'Dimension kz read error')
      istatus = nf90_inq_dimid(ncinp, 'time', idimid)
      call check_ok(istatus,'Dimension time missing')
      istatus = nf90_inquire_dimension(ncinp, idimid, len=nrec)
      call check_ok(istatus,'Dimension time read error')
      istatus = nf90_inq_varid(ncinp, 'time', ivarid)
      call check_ok(istatus,'variable time missing')
      istatus = nf90_get_att(ncinp, ivarid, 'units', timeunits)
      call check_ok(istatus,'variable time units missing')
      allocate(itimes(nrec))
      allocate(xtimes(nrec))
      istatus = nf90_get_var(ncinp, ivarid, xtimes)
      call check_ok(istatus,'variable time read error')
      do i = 1 , nrec
        itimes(i) = timeval2idate(xtimes(i), timeunits)
      end do
      deallocate(xtimes)

!     Reserve space for I/O

      allocate(sig(kz_in), stat=istatus)
      if (istatus /= 0) stop 'Allocation Error in headnest: sig'
      allocate(b2(jx_in,iy_in,np*3), stat=istatus)
      if (istatus /= 0) stop 'Allocation Error in headnest: b2'
      allocate(d2(jx_in,iy_in,np*2), stat=istatus)
      if (istatus /= 0) stop 'Allocation Error in headnest: d2'
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
      allocate(b3(iy,jx,np*3),stat=istatus)
      if (istatus /= 0) stop 'Allocation Error in headnest: b3'
      allocate(d3(iy,jx,np*2),stat=istatus)
      if (istatus /= 0) stop 'Allocation Error in headnest: d3'
      allocate(z1(iy_in,jx_in,kz_in),stat=istatus)
      if (istatus /= 0) stop 'Allocation Error in headnest: z1'

      istatus = nf90_inq_varid(ncinp, 'sigma', ivarid) 
      call check_ok(istatus,'variable sigma error')
      istatus = nf90_get_var(ncinp, ivarid, sig)
      call check_ok(istatus,'variable sigma read error')
      istatus = nf90_inq_varid(ncinp, 'xlat', ivarid) 
      call check_ok(istatus,'variable xlat error')
      istatus = nf90_get_var(ncinp, ivarid, xlat_in)
      call check_ok(istatus,'variable xlat read error')
      istatus = nf90_inq_varid(ncinp, 'xlon', ivarid) 
      call check_ok(istatus,'variable xlon error')
      istatus = nf90_get_var(ncinp, ivarid, xlon_in)
      call check_ok(istatus,'variable xlon read error')
      istatus = nf90_inq_varid(ncinp, 'topo', ivarid) 
      call check_ok(istatus,'variable topo error')
      istatus = nf90_get_var(ncinp, ivarid, ht_in)
      call check_ok(istatus,'variable topo read error')
      istatus = nf90_inq_varid(ncinp, 'ptop', ivarid) 
      call check_ok(istatus,'variable ptop error')
      istatus = nf90_get_var(ncinp, ivarid, ptop_in)
      call check_ok(istatus,'variable ptop read error')

      istatus = nf90_get_att(ncinp, nf90_global, &
                   &    'projection', iproj_in)
      call check_ok(istatus,'attribure iproj read error')
      istatus = nf90_get_att(ncinp, nf90_global, &
                   &    'latitude_of_projection_origin', clat_in)
      call check_ok(istatus,'attribure clat read error')
      istatus = nf90_get_att(ncinp, nf90_global, &
                   &    'longitude_of_projection_origin', clon_in)
      call check_ok(istatus,'attribure clat read error')

      if ( iproj_in=='LAMCON' ) then
        istatus = nf90_get_att(ncinp, nf90_global, &
                   &    'standard_parallel', trlat)
        call check_ok(istatus,'attribure truelat read error')
        if ( clat_in<0. ) then
          xsign = -1.0D0       ! SOUTH HEMESPHERE
        else
          xsign =  1.0D0       ! NORTH HEMESPHERE
        end if
        if ( abs(trlat(1)-trlat(2))>1.E-1 ) then
          grdfac_in = (log10(cos(trlat(1)*degrad))                      &
                   & -log10(cos(trlat(2)*degrad)))                      &
                   & /(log10(tan((45.0-xsign*trlat(1)/2.0)*degrad))     &
                   & -log10(tan((45.0-xsign*trlat(2)/2.0)*degrad)))
        else
          grdfac_in = xsign*dsin(dble(trlat(1))*degrad)
        end if
      else if ( iproj_in=='POLSTR' ) then
        grdfac_in = 1.0D0
      else if ( iproj_in=='NORMER' ) then
        grdfac_in = 0.0D0
      else
        istatus = nf90_get_att(ncinp, nf90_global, &
                   &    'grid_north_pole_latitude', plat_in)
        call check_ok(istatus,'attribure plat read error')
        istatus = nf90_get_att(ncinp, nf90_global, &
                   &    'grid_north_pole_longitude', plon_in)
        call check_ok(istatus,'attribure plon read error')
        grdfac_in = 0.0D0
      end if
 
!     Set up pointers
 
      tp => b2(:,:,1:np)
      qp => b2(:,:,np+1:2*np)
      hp => b2(:,:,2*np+1:3*np)
      up => d2(:,:,1:np)
      vp => d2(:,:,np+1:2*np)
      t3 => b3(:,:,1:np)
      q3 => b3(:,:,np+1:2*np)
      h3 => b3(:,:,2*np+1:3*np)
      u3 => d3(:,:,1:np)
      v3 => d3(:,:,np+1:2*np)

      end subroutine headnest
!
      subroutine check_ok(ierr,message)
        use netcdf
        implicit none
        integer , intent(in) :: ierr
        character(*) :: message
        if (ierr /= nf90_noerr) then
          write (6,*) message
          write (6,*) nf90_strerror(ierr)
          stop
        end if
      end subroutine check_ok
!
      end module mod_nest
