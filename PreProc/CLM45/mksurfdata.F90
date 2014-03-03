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

program mksurfdata
 
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_message
  use mod_grid
  use mod_date
  use mod_memutil
  use mod_stdio
  use mod_domain
  use mod_nchelper
  use mod_mkglacier
  use mod_mkwetland
  use mod_mkurban
  use mod_mkpft
  use netcdf

  implicit none

  integer , parameter :: maxd3 = 17
  integer , parameter :: maxd4 = 12
  integer , parameter :: npft = 17
  integer , parameter :: nsoil = 10

  real(rk4) , parameter :: vmisdat = -9999.0
  logical , parameter :: bvoc = .false.

  integer(ik4) :: istatus , ncid , ndim , nvar
  integer(ik4) , dimension(6) :: idims , ivdims
  integer(ik4) :: ivartime , iglcvar , iwetvar , ilakevar , iurbanvar , ipftvar
  type(rcm_time_and_date) :: irefdate , imondate
  type(rcm_time_interval) :: tdif
  real(rk4) , pointer , dimension(:) :: yiy
  real(rk4) , pointer , dimension(:) :: xjx
  real(rk4) :: hptop
  real(rk8) , dimension(1) :: xdate
  integer(ik4) , dimension(2) :: ihvar
  integer(ik4) , dimension(2) :: illvar
  integer(ik4) , dimension(2) :: izvar
  integer(ik4) , dimension(1) :: istart1 , icount1
  real(rk4) :: spft , diff
  real(rk8) :: operat
  integer(ik4) :: ierr
  integer(ik4) :: i , j , np , it , ipnt
  character(len=256) :: namelistfile , prgname
  character(len=256) :: terfile , outfile
  character(len=64) :: csdate
  real(rk4) , dimension(:,:) , pointer :: pctspec
  real(rk4) , pointer , dimension(:,:,:,:) :: var4d
  logical , dimension(npft) :: pft_gt0

  call get_command_argument(0,value=prgname)
  call get_command_argument(1,value=namelistfile)
  call initparam(namelistfile, ierr)
  if ( ierr/=0 ) then
    write(stderr,*) 'Parameter initialization not completed'
    write(stderr,*) 'Usage : '
    write(stderr,*) '          ', trim(prgname), ' regcm.in'
    write(stderr,*) ' '
    call die('clm2rcm','Check argument and namelist syntax',1)
  end if

  call memory_init

  call init_domain
  !
  ! Get latitudes, longitudes and mask from DOMAIN file
  !
  if ( nsg > 1 ) then
    write (terfile,'(a,i0.3,a)') &
         trim(dirter)//pthsep//trim(domname)//'_DOMAIN',nsg,'.nc'
  else
    terfile = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
  end if
  call openfile_withname(terfile,ncid)
  call read_domain(ncid,sigx,xlat,xlon,ht=topo,mask=xmask)
  call closefile(ncid)

  ! Open Output in NetCDF format

  outfile = trim(dirglob)//pthsep//trim(domname)//'_CLM45_surface.nc'

  call createfile_withname(outfile,ncid)
  call add_common_global_params(ncid,'clm2rcm',.false.)
  ndim = 1
  call define_basic_dimensions(ncid,jx,iy,kzp1,ndim,idims)
  call add_dimension(ncid,'time',nf90_unlimited,ndim,idims)
  call add_dimension(ncid,'lsmpft',npft,ndim,idims)
  call add_dimension(ncid,'nlevsoi',nsoil,ndim,idims)
  call define_horizontal_coord(ncid,jx,iy,xjx,yiy,idims,ihvar)
  call define_vertical_coord(ncid,idims,izvar)

  nvar = 1
  call define_cross_geolocation_coord(ncid,idims,nvar,illvar)

  istatus = nf90_def_var(ncid, 'time', nf90_double, idims(4:4), ivartime)
  call checkncerr(istatus,__FILE__,__LINE__,'Error add variable time')
  irefdate = yrfirst(globidate1)
  irefdate = monmiddle(irefdate)
  csdate = tochar(irefdate)
  istatus = nf90_put_att(ncid, ivartime, 'units', 'hours since '//csdate)
  call checkncerr(istatus,__FILE__,__LINE__,'Error add time units')
  istatus = nf90_put_att(ncid, ivartime, 'calendar', calendar)
  call checkncerr(istatus,__FILE__,__LINE__,'Error add time calendar')

  ! Variables

  istatus = nf90_def_var(ncid, 'PCT_GLACIER', nf90_float, idims(1:2), iglcvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var glacier')
  istatus = nf90_put_att(ncid, iglcvar, 'long_name','percent glacier')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add glacier long_name')
  istatus = nf90_put_att(ncid, iglcvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add glacier units')
  istatus = nf90_put_att(ncid, iglcvar, '_FillValue',vmisdat)
  call checkncerr(istatus,__FILE__,__LINE__,'Error add glacier _FillValue')
  istatus = nf90_put_att(ncid, iglcvar, 'coordinates','xlat xlon')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add glacier coordinates')

  istatus = nf90_def_var(ncid, 'PCT_WETLAND', nf90_float, idims(1:2), iwetvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var wetland')
  istatus = nf90_put_att(ncid, iwetvar, 'long_name','percent wetland')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add wetland long_name')
  istatus = nf90_put_att(ncid, iwetvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add wetland units')
  istatus = nf90_put_att(ncid, iwetvar, '_FillValue',vmisdat)
  call checkncerr(istatus,__FILE__,__LINE__,'Error add wetland _FillValue')
  istatus = nf90_put_att(ncid, iwetvar, 'coordinates','xlat xlon')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add wetland coordinates')

  istatus = nf90_def_var(ncid, 'PCT_LAKE', nf90_float, idims(1:2), ilakevar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var lake')
  istatus = nf90_put_att(ncid, ilakevar, 'long_name','percent lake')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add lake long_name')
  istatus = nf90_put_att(ncid, ilakevar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add lake units')
  istatus = nf90_put_att(ncid, ilakevar, '_FillValue',vmisdat)
  call checkncerr(istatus,__FILE__,__LINE__,'Error add lake _FillValue')
  istatus = nf90_put_att(ncid, ilakevar, 'coordinates','xlat xlon')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add lake coordinates')

  istatus = nf90_def_var(ncid, 'PCT_URBAN', nf90_float, idims(1:2), iurbanvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var urban')
  istatus = nf90_put_att(ncid, iurbanvar, 'long_name','percent urban')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add urban long_name')
  istatus = nf90_put_att(ncid, iurbanvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add urban units')
  istatus = nf90_put_att(ncid, iurbanvar, '_FillValue',vmisdat)
  call checkncerr(istatus,__FILE__,__LINE__,'Error add urban _FillValue')
  istatus = nf90_put_att(ncid, iurbanvar, 'coordinates','xlat xlon')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add urban coordinates')

  ivdims(1:2) = idims(1:2)
  ivdims(3) = idims(5)
  istatus = nf90_def_var(ncid, 'PCT_PFT', nf90_float, ivdims(1:3), ipftvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var pft')
  istatus = nf90_put_att(ncid, ipftvar, 'long_name','percent pft')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add pft long_name')
  istatus = nf90_put_att(ncid, ipftvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add pft units')
  istatus = nf90_put_att(ncid, ipftvar, '_FillValue',vmisdat)
  call checkncerr(istatus,__FILE__,__LINE__,'Error add pft _FillValue')
  istatus = nf90_put_att(ncid, ipftvar, 'coordinates','xlat xlon')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add pft coordinates')

  istatus = nf90_enddef(ncid)
  call checkncerr(istatus,__FILE__,__LINE__,'Error exit define mode')

  hptop = real(ptop*10.0D0)
  call write_vertical_coord(ncid,sigx,hptop,izvar)
  call write_horizontal_coord(ncid,xjx,yiy,ihvar)
  ipnt = 1
  call write_var2d_static(ncid,'xlat',xlat,ipnt,illvar)
  call write_var2d_static(ncid,'xlon',xlon,ipnt,illvar)

  imondate = irefdate
  do it = 1 , 12
    istart1(1) = it
    icount1(1) = 1
    tdif = imondate-irefdate 
    xdate(1) = tohours(tdif)
    istatus = nf90_put_var(ncid, ivartime, xdate, istart1, icount1)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable time write')
    imondate = nextmon(imondate)
    imondate = monmiddle(imondate)
  end do

  call getmem2d(pctspec,1,jx,1,iy,'mksurfdata: pctspec')
  call getmem4d(var4d,1,jx,1,iy,1,maxd3,1,maxd4,'mksurfdata: var4d')
  pctspec(:,:) = 0.0
 
  call mkglacier('mksrf_glacier.nc',var4d(:,:,1,1))
  where ( xmask < 0.5 )
    var4d(:,:,1,1) = vmisdat
  end where
  where (var4d(:,:,1,1) >= 0.0)
    pctspec(:,:) = pctspec(:,:) + var4d(:,:,1,1)
  end where
  istatus = nf90_put_var(ncid, iglcvar, var4d(:,:,1,1))
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write glacier')

  call mkwetland('mksrf_lanwat.nc',var4d(:,:,1,1),var4d(:,:,2,1))
  where ( xmask < 0.5 )
    var4d(:,:,1,1) = vmisdat
    var4d(:,:,2,1) = vmisdat
  end where
  where (var4d(:,:,1,1) >= 0.0)
    pctspec(:,:) = pctspec(:,:) + var4d(:,:,1,1)
  end where
  where (var4d(:,:,2,1) >= 0.0)
    pctspec(:,:) = pctspec(:,:) + var4d(:,:,2,1)
  end where
  istatus = nf90_put_var(ncid, iwetvar, var4d(:,:,1,1))
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write wetland')
  istatus = nf90_put_var(ncid, ilakevar, var4d(:,:,2,1))
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write lake')

  call mkurban('mksrf_urban.nc',var4d(:,:,1,1))
  where ( xmask < 0.5 )
    var4d(:,:,1,1) = vmisdat
  end where
  where (var4d(:,:,1,1) >= 0.0)
    pctspec(:,:) = pctspec(:,:) + var4d(:,:,1,1)
  end where
  istatus = nf90_put_var(ncid, iurbanvar, var4d(:,:,1,1))
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write urban')

  call mkpft('mksrf_pft.nc',var4d(:,:,1:npft,1))
  do np = 1 , npft
    where ( xmask < 0.5 )
      var4d(:,:,np,1) = vmisdat
    end where
  end do
  ! Here adjustment !
  do i = 1 , iy
    do j = 1 , jx
      if ( xmask(j,i) > 0.5 ) then
        if ( pctspec(j,i) > 99.99999 ) then
          var4d(j,i,:,1) = vmisdat
        else if ( pctspec(j,i) < 0.00001 ) then
          var4d(j,i,:,1) = 0.0
          var4d(j,i,1,1) = 100.0
        else
          spft = 0.0
          if ( pctspec(j,i) > 0.00001 ) then
            do np = 1 , npft
              if ( var4d(j,i,np,1) > 0.0 ) then
                spft = spft + var4d(j,i,np,1) * 100.0/(100.0-pctspec(j,i))
              end if
            end do
          else 
            do np = 1 , npft
              if ( var4d(j,i,np,1) > 0.0 ) then
                spft = spft + var4d(j,i,np,1)
              end if
            end do
          end if
          diff = spft - 100.0
          if ( abs(diff) > 1.0E-5 ) then
            ! Normalize it !
            if ( abs(diff) > 0.0001 ) then
              ! Not worth doing if too small a diff...
              do np = 1 , npft
                if ( var4d(j,i,np,1) > 0.0 ) then
                  operat = dble(diff)*(dble(var4d(j,i,np,1))/dble(spft))
                  var4d(j,i,np,1) = var4d(j,i,np,1)-real(operat)
                  if ( var4d(j,i,np,1) < 0.0 ) then
                    var4d(j,i,np,1) = 0.0
                  end if
                end if
              end do
              ! Re-compute diff
              spft = 0.0
              if ( pctspec(j,i) > 0.00001 ) then
                do np = 1 , npft
                  if ( var4d(j,i,np,1) > 0.0 ) then
                    spft = spft + var4d(j,i,np,1) * 100.0/(100.0-pctspec(j,i))
                  end if
                end do
              else
                do np = 1 , npft
                  if ( var4d(j,i,np,1) > 0.0 ) then
                    spft = spft + var4d(j,i,np,1)
                  end if
                end do
              end if
              diff = spft - 100.0
            end if
            if ( abs(diff) > 1.0E-5 ) then
              pft_gt0 = (var4d(j,i,:,1) > diff/real(npft))
              where (pft_gt0)
                var4d(j,i,:,1) = var4d(j,i,:,1) - diff/real(count(pft_gt0))
              end where
            end if
          end if
        end if
      end if
    end do
  end do
  istatus = nf90_put_var(ncid, ipftvar, var4d(:,:,1:npft,1))
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write pft')
  istatus = nf90_close(ncid)
  call checkncerr(istatus,__FILE__,__LINE__,  &
    'Error close file '//trim(outfile))

  call memory_destroy

  write(stdout,*) 'Successfully completed CLM preprocessing.'

  contains

  recursive subroutine sortpatch(vals,svals,ird,lsub)
    implicit none
    real(rk4) , dimension(:) , intent(in) :: vals
    real(rk4) , dimension(:) , intent(inout) :: svals
    integer(ik4) , dimension(:) , intent(inout) :: ird
    logical , optional :: lsub
    integer(ik4) :: i , iswap
    real(rk4) :: rswap
    if ( .not. present(lsub) ) then
      do i = 1 , size(vals)
        ird(i) = i
        svals(i) = vals(i)
      end do
    end if
    do i = 1 , size(vals)-1
      if ( svals(i) < svals(i+1) ) then
        rswap = svals(i+1)
        iswap = ird(i+1)
        svals(i+1) = svals(i)
        ird(i+1) = ird(i)
        svals(i) = rswap
        ird(i) = iswap
        call sortpatch(vals,svals,ird,.true.)
      end if
    end do
  end subroutine sortpatch
 
end program mksurfdata
