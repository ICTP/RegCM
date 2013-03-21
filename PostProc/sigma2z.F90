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

subroutine myabort
  implicit none
  call abort
end subroutine myabort

program sigma2p

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_message
  use mod_vertint
  use mod_derived
  use mod_nchelper
  use netcdf

  implicit none

  integer(ik4) , parameter :: nz = 14
  real(rk4) , dimension(nz) :: zlevs

  character(256) :: prgname , ncsfile , ncpfile
  character(128) :: attname , dimname , varname
  integer(ik4) :: numarg , istatus , ncid , ncout

  integer(ik4) , allocatable , dimension(:) :: dimids , dimlen
  real(rk4) , allocatable , dimension(:) :: sigma
  real(rk4) , allocatable , dimension(:,:,:) :: xvar , tempvar
  real(rk4) , allocatable , dimension(:,:,:) :: zvar , tazvar , qazvar , hazvar
  real(rk4) , allocatable , dimension(:,:) :: ps , topo , mslp
  real(rk4) , allocatable , dimension(:) :: avar
  character , allocatable , dimension(:) :: tvar
  real(rk4) , allocatable , dimension(:) :: azvar
  real(rk8) , allocatable , dimension(:) :: times
  real(rk8) , allocatable , dimension(:) :: sigfix
  logical , allocatable , dimension(:) :: lkvarflag , ltvarflag , lchnameflag
  integer(ik4) , allocatable , dimension(:) :: varsize
  integer(ik4) , allocatable , dimension(:) :: intscheme
  integer(ik4) , allocatable , dimension(:) :: nvdims
  integer(ik4) , allocatable , dimension(:,:) :: dimsize
  integer(ik4) , allocatable , dimension(:) :: istart , icount
  integer(ik4) :: ndims , nvars , natts , udimid , nvatts
  integer(ik4) :: ivarid , idimid , xtype
  integer(ik4) :: jxdimid , iydimid , kzdimid , itdimid , itvarid , ikvarid
  integer(ik4) :: ipsvarid , ishvarid
  integer(ik4) :: jx , iy , kz , nt
  real(rk8) :: ptop
  integer(ik4) , dimension(4) :: tdimids
  integer(ik4) , dimension(3) :: psdimids
  integer(ik4) :: i , j , k , it , iv , iid1 , iid2 , ii , i3d , p3d , ich
  integer(ik4) :: tvarid , qvarid , irhvar , ihgvar , imslzvar , ircm_map
  logical :: has_t , has_q
  logical :: make_rh , make_hgt
  integer(ik4) :: n3d , iz3d
#if defined ( __PGI ) || defined ( IBM ) || defined ( __OPENCC__ )
  integer(ik4) , external :: iargc
#endif

  data has_t /.false./
  data has_q /.false./
  data make_rh /.false./
  data make_hgt /.false./

  data zlevs /20.0,50.0,80.0,100.0,150.0,200.0,500.0,750.0,1000.0, &
              1500.0,2000.0,5000.0,7000.0,10000.0/

  call getarg(0, prgname)
  numarg = iargc( )
  if (numarg < 1) then
    write (6,*) 'Not enough arguments.'
    write (6,*) ' '
    write (6,*) 'Usage : ', trim(prgname), ' Rcmfile.nc'
    write (6,*) ' '
    stop
  end if

  call getarg(1, ncsfile)
  iid1 = scan(ncsfile, '/', .true.)
  iid2 = scan(ncsfile, '.', .true.)
  ncpfile = trim(ncsfile(iid1+1:iid2-1))//'_hgt.nc'

  istatus = nf90_open(ncsfile, nf90_nowrite, ncid)
  call checkncerr(istatus,__FILE__,__LINE__,'Error Opening Input file '//trim(ncsfile))

  jxdimid = -1
  iydimid = -1
  kzdimid = -1
  itdimid = -1
  itvarid = -1
  ikvarid = -1
  iy = -1
  jx = -1
#ifdef NETCDF4_HDF5
  istatus = nf90_create(ncpfile, &
            ior(ior(nf90_clobber,nf90_hdf5),nf90_classic_model), &
            ncout)
#else
  istatus = nf90_create(ncpfile, nf90_clobber, ncout)
#endif
  call checkncerr(istatus,__FILE__,__LINE__,'Error Opening Output file '//trim(ncpfile))

  istatus = nf90_inquire(ncid,ndims,nvars,natts,udimid)
  call checkncerr(istatus,__FILE__,__LINE__,'Error Reading Netcdf file '//trim(ncsfile))

  do i = 1 , natts
    istatus = nf90_inq_attname(ncid, nf90_global, i, attname)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read global attribute')
    istatus = nf90_copy_att(ncid, nf90_global, attname, ncout, nf90_global)
    call checkncerr(istatus,__FILE__,__LINE__,'Error copy attribute '//trim(attname))
  end do

  allocate(dimlen(ndims), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'dimlen')
 
  kz = 0
  nt = 0
  do i = 1 , ndims
    istatus = nf90_inquire_dimension(ncid, i, dimname, dimlen(i))
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading dimension info')
    if (dimname == 'iy' .or. dimname == 'y') then
      iy = dimlen(i)
      iydimid = i
    else if (dimname == 'jx' .or. dimname == 'x') then
      jx = dimlen(i)
      jxdimid = i
    else if (dimname == 'kz' .or. dimname == 'z' .or. dimname == 'lev') then
      kz = dimlen(i)
      kzdimid = i
    end if
    if (dimname == 'time') then
      itdimid = i
      nt = dimlen(i)
      istatus = nf90_def_dim(ncout, dimname, nf90_unlimited, idimid)
    else if (dimname == 'kz' .or. dimname == 'z' .or. dimname == 'lev') then
      istatus = nf90_def_dim(ncout, 'zlev', nz, idimid)
    else
      istatus = nf90_def_dim(ncout, dimname, dimlen(i), idimid)
    end if
    call checkncerr(istatus,__FILE__,__LINE__,'Error define dimension '//trim(dimname))
  end do

  i3d = jx*iy*kz
  p3d = jx*iy*nz

  if (kz == 0 .or. nt == 0) then
    write (6,*) 'Nothing to do: no vertical dimension or no time'
    istatus = nf90_close(ncout)
    istatus = nf90_close(ncid)
    stop
  end if

  allocate(dimids(ndims), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'dimids')
  allocate(lkvarflag(nvars), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'lkvarflag')
  allocate(ltvarflag(nvars), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'ltvarflag')
  allocate(lchnameflag(nvars), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'lchnameflag')
  allocate(varsize(nvars), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'varsize')
  allocate(intscheme(nvars), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'intscheme')
  allocate(nvdims(nvars), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'nvdims')
  allocate(dimsize(ndims,nvars), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'dimsize')
  allocate(istart(ndims), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'istart')
  allocate(icount(ndims), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'icount')
  allocate(ps(jx,iy), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'ps')
  allocate(xvar(jx,iy,kz), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'xvar')
  allocate(tempvar(jx,iy,kz), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'tempvar')
  allocate(zvar(jx,iy,nz), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'zvar')
  allocate(times(nt), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'times')
  allocate(sigma(kz), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'sigma')

  itvarid = 0
  do i = 1 , nvars
    lkvarflag(i) = .false.
    ltvarflag(i) = .false.
    lchnameflag(i) = .false.
    varsize(i) = 1
    intscheme(i) = 1
    dimsize(:,i) = 1
    istatus = nf90_inquire_variable(ncid,i,varname,xtype,nvdims(i), &
                                    dimids,nvatts)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire variable')
    if (varname == 'time') then
      itvarid = i
    else if (varname == 'sigma' .or. varname == 'lev') then
      varname = 'zlev'
      ikvarid = i
    else if (varname == 'ptop') then
      istatus = nf90_get_var(ncid,i,ptop)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read variable ptop')
    else if (varname == 'rcm_map') then
      ircm_map = i
    else if (varname == 'ps') then
      ipsvarid = i
      psdimids = dimids(1:3)
    else if (varname == 'topo') then
      ishvarid = i
    else if (varname == 'ta' .or. varname == 't') then
      has_t = .true.
      intscheme(i) = 2
      tvarid = i
      tdimids = dimids(1:4)
    else if (varname == 'qas' .or. varname == 'qv') then
      has_q = .true.
      qvarid = i
    else if (varname == 'chtrname') then
      lchnameflag(i) = .true.
    end if
    iv = 1
    do j = 1 , nvdims(i)
      if (dimids(j) == kzdimid) then
        lkvarflag(i) = .true.
      else if (dimids(j) == itdimid) then
        ltvarflag(i) = .true.
        iv = iv + 1
        cycle
      end if
      dimsize(iv,i) = dimlen(dimids(j))
      iv = iv + 1
      varsize(i) = varsize(i) * dimlen(dimids(j))
    end do
    istatus = nf90_def_var(ncout, varname, xtype, dimids(1:nvdims(i)), ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error define variable '//trim(varname))
#ifdef NETCDF4_HDF5
    if (nvdims(i) > 2) then
      istatus = nf90_def_var_deflate(ncout, ivarid, 1, 1, 9)
      call checkncerr(istatus,__FILE__,__LINE__,'Error set deflate for '//trim(varname))
    end if
#endif
    if (varname == 'zlev') then
      istatus = nf90_put_att(ncout, ivarid, 'standard_name', 'height')
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding hgt standard name')
      istatus = nf90_put_att(ncout, ivarid, 'long_name', &
                             'elevation')
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding hgt long name')
      istatus = nf90_put_att(ncout, ivarid, 'units', 'm')
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding hgt units')
      istatus = nf90_put_att(ncout, ivarid, 'axis', 'Z')
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding hgt axis')
      istatus = nf90_put_att(ncout, ivarid, 'positive', 'up')
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding hgt positive')
      cycle
    end if
    do j = 1 , nvatts
      istatus = nf90_inq_attname(ncid, i, j, attname)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read att for '//trim(varname))
      istatus = nf90_copy_att(ncid, i, attname, ncout, ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error copy att '//trim(attname)//' for '//trim(varname))
    end do
  end do

  if ( has_t ) then
    make_hgt = .true.
    istart(1) = 1
    istart(2) = 1
    icount(1) = jx
    icount(2) = iy
    allocate(topo(jx,iy), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'topo')
    allocate(mslp(jx,iy), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'mslp')
    allocate(hazvar(jx,iy,nz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'hazvar')
    istatus = nf90_get_var(ncid, ishvarid, topo, istart(1:2), icount(1:2))
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading topo.')
    istatus = nf90_def_var(ncout, 'hgt', nf90_float, tdimids, ihgvar)
    call checkncerr(istatus,__FILE__,__LINE__,'Error define variable rh')
#ifdef NETCDF4_HDF5
    istatus = nf90_def_var_deflate(ncout, ihgvar, 1, 1, 9)
    call checkncerr(istatus,__FILE__,__LINE__,'Error set deflate for hgt')
#endif
    istatus = nf90_put_att(ncout, ihgvar, 'standard_name', 'height')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding standard name')
    istatus = nf90_put_att(ncout, ihgvar, 'long_name', 'Vertical distance above surface')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding long name')
    istatus = nf90_put_att(ncout, ihgvar, 'units', 'm')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding units')
    istatus = nf90_put_att(ncout, ihgvar, '_FillValue', smissval)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding missval')
    istatus = nf90_def_var(ncout, 'mslp', nf90_float, psdimids, imslzvar)
    call checkncerr(istatus,__FILE__,__LINE__,'Error define variable rh')
#ifdef NETCDF4_HDF5
    istatus = nf90_def_var_deflate(ncout, imslzvar, 1, 1, 9)
    call checkncerr(istatus,__FILE__,__LINE__,'Error set deflate for hgt')
#endif
    istatus = nf90_put_att(ncout, imslzvar, 'standard_name', &
                     'air_pressure_at_sea_level')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding standard name')
    istatus = nf90_put_att(ncout, imslzvar, 'long_name', &
                     'Sea Level pressure')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding long name')
    istatus = nf90_put_att(ncout, imslzvar, 'units', 'hPa')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding units')
  end if
  if ( has_t .and. has_q ) then
    make_rh = .true.
    allocate(tazvar(jx,iy,nz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'tazvar')
    allocate(qazvar(jx,iy,nz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'qazvar')
    istatus = nf90_def_var(ncout, 'rh', nf90_float, tdimids, irhvar)
    call checkncerr(istatus,__FILE__,__LINE__,'Error define variable rh')
#ifdef NETCDF4_HDF5
    istatus = nf90_def_var_deflate(ncout, irhvar, 1, 1, 9)
    call checkncerr(istatus,__FILE__,__LINE__,'Error set deflate for rh')
#endif
    istatus = nf90_put_att(ncout, irhvar, 'standard_name', 'relative_humidity')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding standard name')
    istatus = nf90_put_att(ncout, irhvar, 'long_name', 'Relative Humidity')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding long name')
    istatus = nf90_put_att(ncout, irhvar, 'units', '%')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding units')
  end if

  istatus = nf90_inq_varid(ncid, "sigma", ivarid)
  if ( istatus /= nf90_noerr ) then
    istatus = nf90_inq_varid(ncid, "lev", ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading variable sigma.')
  end if
  istatus = nf90_get_var(ncid, ivarid, sigma)
  call checkncerr(istatus,__FILE__,__LINE__,'Error reading variable sigma.')
  if ( sigma(1) < dlowval ) then
    ! Fix for a buggy RegCM 4.3.x revision
    allocate(sigfix(kz+1))
    sigfix(1:kz) = sigma
    sigfix(kz+1) = d_one
    do k = 1 , kz
      sigma(k) = 0.5*real(sigfix(k)+sigfix(k+1))
    end do
    deallocate(sigfix)
  end if

  istatus = nf90_inq_varid(ncid, "time", ivarid)
  call checkncerr(istatus,__FILE__,__LINE__,'Error reading variable time.')
  istatus = nf90_get_var(ncid, ivarid, times)
  call checkncerr(istatus,__FILE__,__LINE__,'Error reading variable time.')

  istatus = nf90_enddef(ncout)
  call checkncerr(istatus,__FILE__,__LINE__,'Error preparing for write on output')

! Write time independent variables

  do i = 1 , nvars
    if (i == itvarid) cycle
    if (i == ircm_map) cycle
    if (i == ikvarid) then
      istatus = nf90_put_var(ncout, i, zlevs)
      call checkncerr(istatus,__FILE__,__LINE__,'Error writing variable zlev')
      cycle
    end if
    if (.not. ltvarflag(i)) then
      if ( .not. lchnameflag(i) ) then
        allocate(avar(varsize(i)), stat=istatus)
        call checkalloc(istatus,__FILE__,__LINE__,'avar')
        iv = nvdims(i)
        istart(:) = 1
        icount(1:iv) = dimsize(1:iv,i)
        istatus = nf90_get_var(ncid, i, avar, istart(1:iv), icount(1:iv))
        call checkncerr(istatus,__FILE__,__LINE__,'Error reading variable.')
        istatus = nf90_put_var(ncout, i, avar, istart(1:iv), icount(1:iv))
        call checkncerr(istatus,__FILE__,__LINE__,'Error writing variable.')
        deallocate(avar)
      else
        allocate(tvar(varsize(i)), stat=istatus)
        call checkalloc(istatus,__FILE__,__LINE__,'tvar')
        iv = nvdims(i)
        istart(:) = 1
        icount(1:iv) = dimsize(1:iv,i)
        istatus = nf90_get_var(ncid, i, tvar, istart(1:iv), icount(1:iv))
        call checkncerr(istatus,__FILE__,__LINE__,'Error reading variable.')
        istatus = nf90_put_var(ncout, i, tvar, istart(1:iv), icount(1:iv))
        call checkncerr(istatus,__FILE__,__LINE__,'Error writing variable.')
        deallocate(tvar)
      end if
    end if
  end do

! Write time dependent variables

  do it = 1 , nt
    istart(1) = it
    icount(1) = 1
    istatus = nf90_put_var(ncout, itvarid, times(it:it), &
                           istart(1:1), icount(1:1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error writing time.')
    istart(1) = 1
    istart(2) = 1
    istart(3) = it
    icount(1) = jx
    icount(2) = iy
    icount(3) = 1
    istatus = nf90_get_var(ncid, ipsvarid, ps, istart(1:3), icount(1:3))
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading ps.')
    istatus = nf90_put_var(ncout, ipsvarid, ps, istart(1:3), icount(1:3))
    call checkncerr(istatus,__FILE__,__LINE__,'Error writing ps.')
    istart(1) = 1
    istart(2) = 1
    istart(3) = 1
    istart(4) = it
    icount(1) = jx
    icount(2) = iy
    icount(3) = kz
    icount(4) = 1
    istatus = nf90_get_var(ncid, tvarid, tempvar, istart(1:4), icount(1:4))
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading temp.')
    do i = 1 , nvars
      if (.not. ltvarflag(i)) cycle
      if (i == itvarid) cycle
      if (i == ipsvarid) cycle

      if (lkvarflag(i)) then

!       Do interpolation

        iv = nvdims(i)
        if ( iv == 4 ) then
          istart(iv) = it
          icount(iv) = 1
          istart(1:iv-1) = 1
          icount(1:iv-1) = dimsize(1:iv-1,i)
          allocate(avar(varsize(i)),stat=istatus)
          call checkalloc(istatus,__FILE__,__LINE__,'avar')
          istatus = nf90_get_var(ncid, i, avar, istart(1:iv), icount(1:iv))
          call checkncerr(istatus,__FILE__,__LINE__,'Error reading var to interpolate.')

          n3d = varsize(i) / i3d
          iz3d = p3d*n3d
          allocate(azvar(iz3d),stat=istatus)
          call checkalloc(istatus,__FILE__,__LINE__,'azvar')
          do ii = 1 , n3d
            xvar = reshape(avar((ii-1)*i3d+1:ii*i3d),(/jx,iy,kz/))
            call intlin_o(zvar,xvar,tempvar,ps,sigma,ptop,jx,iy,kz,zlevs,nz)
            azvar((ii-1)*iz3d+1:ii*iz3d) = reshape(zvar,(/iz3d/))
          end do
          if ( i == tvarid ) then
            tazvar = zvar
          else if ( i == qvarid ) then
            qazvar = zvar
          end if

          icount(iv-1) = nz
          istatus = nf90_put_var(ncout, i, azvar, istart(1:iv), icount(1:iv))
          call checkncerr(istatus,__FILE__,__LINE__,'Error writing interp variable.')

          deallocate(avar)
          deallocate(azvar)
        else if ( iv == 5 ) then
          istart(iv) = it
          icount(iv) = 1
          istart(1:iv-1) = 1
          icount(1:iv-1) = dimsize(1:iv-1,i)
          allocate(avar(varsize(i)),stat=istatus)
          call checkalloc(istatus,__FILE__,__LINE__,'avar')
          istatus = nf90_get_var(ncid, i, avar, istart(1:iv), icount(1:iv))
          call checkncerr(istatus,__FILE__,__LINE__,'Error reading var to interpolate.')

          n3d = varsize(i) / dimsize(iv-1,i) / i3d
          iz3d = p3d*n3d
          do ich = 1 , dimsize(iv-1,i)
            istart(iv) = it
            icount(iv) = 1
            istart(iv-1) = ich
            icount(iv-1) = 1
            istart(1:iv-2) = 1
            icount(1:iv-2) = dimsize(1:iv-2,i)
            icount(iv-2) = nz
            allocate(azvar(iz3d),stat=istatus)
            call checkalloc(istatus,__FILE__,__LINE__,'azvar')
            do ii = 1 , n3d
              xvar = reshape(avar((ii-1)*i3d+(ich-1)*i3d+1:(ii+ich-1)*i3d), &
                             (/jx,iy,kz/))
              call intlin_o(zvar,xvar,tempvar,ps,sigma,ptop,jx,iy,kz,zlevs,nz)
              azvar((ii-1)*iz3d+1:ii*iz3d) = reshape(zvar,(/iz3d/))
            end do
            istatus = nf90_put_var(ncout, i, azvar, istart(1:iv), icount(1:iv))
            call checkncerr(istatus,__FILE__,__LINE__,'Error writing interp variable.')
            deallocate(azvar)
          end do
          deallocate(avar)
        end if

      else

!       No interpolation

        iv = nvdims(i)
        istart(iv) = it
        icount(iv) = 1
        istart(1:iv-1) = 1
        icount(1:iv-1) = dimsize(1:iv-1,i)
        allocate(avar(varsize(i)),stat=istatus)
        call checkalloc(istatus,__FILE__,__LINE__,'avar')
        istatus = nf90_get_var(ncid, i, avar, istart(1:iv), icount(1:iv))
        call checkncerr(istatus,__FILE__,__LINE__,'Error reading variable to pass')
        istatus = nf90_put_var(ncout, i, avar, istart(1:iv), icount(1:iv))
        call checkncerr(istatus,__FILE__,__LINE__,'Error writing variable to pass')
        deallocate(avar)

      end if

    end do

    if ( make_rh ) then
      call calc_rh(tazvar,qazvar,zlevs,jx,iy,nz)
      iv = 4
      istart(iv) = it
      icount(iv) = 1
      istart(1:iv-1) = 1
      icount(1:iv-1) = dimsize(1:iv-1,tvarid)
      icount(iv-1) = nz
      istatus = nf90_put_var(ncout, irhvar, qazvar, istart(1:iv), icount(1:iv))
      call checkncerr(istatus,__FILE__,__LINE__,'Error writing rh variable.')
    end if
    if ( make_hgt ) then
      call calc_hgt(hazvar,tazvar,ps,topo,zlevs,jx,iy,nz)
      iv = 4
      istart(iv) = it
      icount(iv) = 1
      istart(1:iv-1) = 1
      icount(1:iv-1) = dimsize(1:iv-1,tvarid)
      icount(iv-1) = nz
      istatus = nf90_put_var(ncout, ihgvar, hazvar, istart(1:iv), icount(1:iv))
      call checkncerr(istatus,__FILE__,__LINE__,'Error writing hgt variable.')
      call calc_mslpres(tazvar,ps,topo,mslp,zlevs,jx,iy,nz)
      iv = 3
      istart(iv) = it
      icount(iv) = 1
      istart(1:iv-1) = 1
      icount(1:iv-1) = dimsize(1:iv-1,tvarid)
      istatus = nf90_put_var(ncout, imslzvar, mslp, istart(1:iv), icount(1:iv))
      call checkncerr(istatus,__FILE__,__LINE__,'Error writing mslp variable.')
    end if

  end do

  deallocate(dimlen)
  deallocate(dimids)
  deallocate(lkvarflag)
  deallocate(ltvarflag)
  deallocate(varsize)
  deallocate(nvdims)
  deallocate(dimsize)
  deallocate(sigma)
  deallocate(times)
  deallocate(ps)
  deallocate(xvar)
  deallocate(zvar)
  if ( make_rh ) then
    deallocate(qazvar)
    deallocate(tazvar)
  end if
  if ( make_hgt ) then
    deallocate(topo)
  end if

  istatus = nf90_close(ncid)
  call checkncerr(istatus,__FILE__,__LINE__,'Error close input file '//trim(ncsfile))

  istatus = nf90_close(ncout)
  call checkncerr(istatus,__FILE__,__LINE__,'Error close output file '//trim(ncpfile))
!
end program sigma2p
