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
  use mod_constants , only : raddeg
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
  use mod_mklaisai
  use mod_mkfmax
  use mod_mksoilcol
  use mod_mksoitex
  use mod_mkgdp
  use mod_mkpeatf
  use mod_mkabm
  use mod_mkvocef
  use mod_mkorganic
  use netcdf

  implicit none

  integer , parameter :: npft = 17
  integer , parameter :: nmon = 12
  integer , parameter :: nsoil = 10

  integer :: ngcells

  integer , parameter :: maxd3 = max(npft,nsoil)
  integer , parameter :: maxd4 = nmon

  real(rk4) , parameter :: vmisdat = -9999.0

  integer(ik4) :: istatus , ncid , ndim , nvar
  integer(ik4) , dimension(7) :: idims , ivdims
  integer(ik4) :: ivartime , iglcvar , iwetvar , ilakevar , iurbanvar
  integer(ik4) :: ipftvar , ilaivar , isaivar , ivgtopvar , ivgbotvar
  integer(ik4) :: ifmaxvar , isoilcolvar , isandvar , iclayvar
  integer(ik4) :: islopevar , istdvar , igdpvar , ipeatfvar , iabmvar
  integer(ik4) :: ief_btrvar , ief_crpvar , ief_fdtvar , ief_fetvar
  integer(ik4) :: ief_grsvar , ief_shrvar , iorganicvar
  integer(ik4) :: ilndvar
  type(rcm_time_and_date) :: irefdate , imondate
  type(rcm_time_interval) :: tdif
  real(rk4) , pointer , dimension(:) :: yiy
  real(rk4) , pointer , dimension(:) :: xjx
  real(rk4) :: hptop
  real(rk8) , dimension(1) :: xdate
  integer(ik4) , dimension(3) :: istart , icount
  integer(ik4) , dimension(2) :: ihvar
  integer(ik4) , dimension(2) :: illvar
  integer(ik4) , dimension(2) :: izvar
  integer(ik4) , dimension(1) :: istart1 , icount1
  real(rk4) :: spft , diff , mean
  real(rk8) :: operat
  integer(ik4) :: ierr
  integer(ik4) :: i , j , ip , il , np , nm , it , ipnt
  character(len=256) :: namelistfile , prgname
  character(len=256) :: terfile , outfile
  character(len=64) :: csdate
  real(rk4) , dimension(:,:) , pointer :: pctspec
  real(rk4) , pointer , dimension(:,:,:,:,:) :: var5d
  real(rk4) , pointer , dimension(:) :: gcvar
  integer , pointer , dimension(:) :: landpoint
  logical , dimension(npft) :: pft_gt0
  logical , dimension(:,:) , pointer :: lmask

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
  ngcells = count(xmask(2:jx-2,2:iy-2) > 0.5)
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
  call add_dimension(ncid,'gridcell',ngcells,ndim,idims)
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

  istatus = nf90_def_var(ncid, 'gridcell', nf90_int, idims(7), ilndvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var gridcell')
  istatus = nf90_put_att(ncid, ilndvar, 'compress','xlat xlon')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add gridcell compress')

  istatus = nf90_def_var(ncid, 'SLOPE', nf90_float, idims(7), islopevar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var slope')
  istatus = nf90_put_att(ncid, islopevar, 'long_name','Elevation slope')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add slope long_name')
  istatus = nf90_put_att(ncid, islopevar, 'units','degree')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add slope units')
  istatus = nf90_def_var(ncid, 'STD_ELEV', nf90_float, idims(7), istdvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var stddev')
  istatus = nf90_put_att(ncid, istdvar, 'long_name','Elevation std dev')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add stddev long_name')
  istatus = nf90_put_att(ncid, istdvar, 'units','m')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add stddev units')

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

  ivdims(1) = idims(7)
  ivdims(2) = idims(5)
  ivdims(3) = idims(4)
  istatus = nf90_def_var(ncid, 'MONTHLY_LAI', nf90_float, ivdims(1:3), ilaivar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var monthly_lai')
  istatus = nf90_put_att(ncid, ilaivar, 'long_name','monthly leaf area index')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_lai long_name')
  istatus = nf90_put_att(ncid, ilaivar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_lai units')
  istatus = nf90_def_var(ncid, 'MONTHLY_SAI', nf90_float, ivdims(1:3), isaivar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var monthly_sai')
  istatus = nf90_put_att(ncid, isaivar, 'long_name','monthly stem area index')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_sai long_name')
  istatus = nf90_put_att(ncid, isaivar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_sai units')
  istatus = nf90_def_var(ncid, 'MONTHLY_HEIGHT_TOP', nf90_float, &
          ivdims(1:3), ivgtopvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var monthly_top')
  istatus = nf90_put_att(ncid, ivgtopvar, 'long_name','monthly height top')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_top long_name')
  istatus = nf90_put_att(ncid, ivgtopvar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_top units')
  istatus = nf90_def_var(ncid, 'MONTHLY_HEIGHT_BOT', nf90_float, &
          ivdims(1:3), ivgbotvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var monthly_bot')
  istatus = nf90_put_att(ncid, ivgbotvar, 'long_name','monthly height bottom')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_bot long_name')
  istatus = nf90_put_att(ncid, ivgbotvar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_bot units')

  istatus = nf90_def_var(ncid, 'FMAX', nf90_float, idims(7), ifmaxvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var fmax')
  istatus = nf90_put_att(ncid, ifmaxvar, 'long_name', &
          'maximum fractional saturated area at 1/8 degree')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add fmax long_name')
  istatus = nf90_put_att(ncid, ifmaxvar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add fmax units')

  istatus = nf90_def_var(ncid, 'SOIL_COLOR', nf90_float, idims(7),isoilcolvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var soilcol')
  istatus = nf90_put_att(ncid, isoilcolvar, 'long_name', &
          'maximum fractional saturated area at 1/8 degree')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add soilcol long_name')
  istatus = nf90_put_att(ncid, isoilcolvar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add soilcol units')

  istatus = nf90_def_var(ncid, 'gdp', nf90_float, idims(7),igdpvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var gdp')
  istatus = nf90_put_att(ncid, igdpvar, 'long_name', 'real GDP')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add gdp long_name')
  istatus = nf90_put_att(ncid, igdpvar, 'units','K 1995US$/capita')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add gdp units')

  istatus = nf90_def_var(ncid, 'peatf', nf90_float, idims(7),ipeatfvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var peatf')
  istatus = nf90_put_att(ncid, ipeatfvar, 'long_name', &
          'global fraction of peatland')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add peatf long_name')
  istatus = nf90_put_att(ncid, ipeatfvar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add peatf units')

  istatus = nf90_def_var(ncid, 'abm', nf90_float, idims(7),iabmvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var abm')
  istatus = nf90_put_att(ncid, iabmvar, 'long_name', &
          'peak month for agri fire')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add abm long_name')
  istatus = nf90_put_att(ncid, iabmvar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add abm units')

  istatus = nf90_def_var(ncid, 'EF1_BTR', nf90_float, idims(7),ief_btrvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ef_btr')
  istatus = nf90_put_att(ncid, ief_btrvar, 'long_name', &
          'broadleaf tree emission factor')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_btr long_name')
  istatus = nf90_put_att(ncid,ief_btrvar,'units','micrograms isoprene m-2 h-1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_btr units')
  istatus = nf90_def_var(ncid, 'EF1_CRP', nf90_float, idims(7),ief_crpvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ef_crp')
  istatus = nf90_put_att(ncid, ief_crpvar, 'long_name', &
          'crop emission factor')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_crp long_name')
  istatus = nf90_put_att(ncid,ief_crpvar,'units','micrograms isoprene m-2 h-1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_crp units')
  istatus = nf90_def_var(ncid, 'EF1_FDT', nf90_float, idims(7),ief_fdtvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ef_fdt')
  istatus = nf90_put_att(ncid, ief_fdtvar, 'long_name', &
          'Fineleaf deciduous tree emission factor')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_fdt long_name')
  istatus = nf90_put_att(ncid,ief_fdtvar,'units','micrograms isoprene m-2 h-1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_fdt units')
  istatus = nf90_def_var(ncid, 'EF1_FET', nf90_float, idims(7),ief_fetvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ef_fet')
  istatus = nf90_put_att(ncid, ief_fetvar, 'long_name', &
          'Fineleaf evergreen tree emission factor')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_fet long_name')
  istatus = nf90_put_att(ncid,ief_fetvar,'units','micrograms isoprene m-2 h-1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_fet units')
  istatus = nf90_def_var(ncid, 'EF1_GRS', nf90_float, idims(7),ief_grsvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ef_grs')
  istatus = nf90_put_att(ncid, ief_grsvar, 'long_name', &
          'grass, non-vascular plants and other ground cover emission factor')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_grs long_name')
  istatus = nf90_put_att(ncid,ief_grsvar,'units','micrograms isoprene m-2 h-1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_grs units')
  istatus = nf90_def_var(ncid, 'EF1_SHR', nf90_float, idims(7),ief_shrvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ef_shr')
  istatus = nf90_put_att(ncid, ief_shrvar, 'long_name', &
          'shrub emission factor')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_shr long_name')
  istatus = nf90_put_att(ncid,ief_shrvar,'units','micrograms isoprene m-2 h-1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_shr units')

  ivdims(1) = idims(7)
  ivdims(2) = idims(6)
  istatus = nf90_def_var(ncid, 'PCT_SAND', nf90_float, ivdims(1:2), isandvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var sand')
  istatus = nf90_put_att(ncid, isandvar, 'long_name','percent sand')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add sand long_name')
  istatus = nf90_put_att(ncid, isandvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add sand units')
  istatus = nf90_def_var(ncid, 'PCT_CLAY', nf90_float, ivdims(1:2), iclayvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var clay')
  istatus = nf90_put_att(ncid, iclayvar, 'long_name','percent clay')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add clay long_name')
  istatus = nf90_put_att(ncid, iclayvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add clay units')

  istatus = nf90_def_var(ncid, 'ORGANIC', nf90_float, ivdims(1:2), iorganicvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var organic')
  istatus = nf90_put_att(ncid, iorganicvar, 'long_name', &
          'organic soil density at soil levels')
  istatus = nf90_put_att(ncid, iorganicvar, 'comment', &
          'assumed carbon content 0.58 gC per gOM')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add organic long_name')
  istatus = nf90_put_att(ncid, iorganicvar, 'units','kg m-3')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add organic units')

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
  call getmem2d(lmask,2,jx-2,2,iy-2,'mksurfdata: lmask')
  call getmem1d(gcvar,1,ngcells,'mksurfdata: gcvar')
  call getmem1d(landpoint,1,ngcells,'mksurfdata: gcvar')
  call getmem5d(var5d,1,jx,1,iy,1,maxd3,1,maxd4,1,4,'mksurfdata: var5d')
  pctspec(:,:) = 0.0
  lmask = (xmask(2:jx-2,2:iy-2) > 0.5)
  ip = 1
  do i = 2 , iy-2
    do j = 2 , jx-2
      if ( xmask(j,i) > 0.5 ) then
        landpoint(ip) = (i-1)*jx+j
        ip = ip + 1
      end if
    end do
  end do
  istatus = nf90_put_var(ncid, ilndvar, landpoint)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write gridcell')
 
  var5d(:,:,1,1,1) = 0.0
  ! Calculate slope and std
  do i = 2 , iy-2
    do j = 2 , jx-2
      var5d(j,i,1,1,1) = &
          atan((sum(topo(j-1:j+1,i-1:i+1)-topo(j,i))/8.0)/real(ds*1000.0D0))
      mean = sum(topo(j-1:j+1,i-1:i+1))/9.0
      var5d(j,i,2,1,1) = sqrt(sum((topo(j-1:j+1,i-1:i+1)-mean)**2)/8.0)
    end do
  end do
  where ( xmask < 0.5 )
    var5d(:,:,1,1,1) = vmisdat
    var5d(:,:,2,1,1) = vmisdat
  end where
  gcvar = pack(var5d(2:jx-2,2:iy-2,1,1,1),lmask)*real(raddeg)
  istatus = nf90_put_var(ncid, islopevar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write slope')
  gcvar = pack(var5d(2:jx-2,2:iy-2,2,1,1),lmask)
  istatus = nf90_put_var(ncid, istdvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write stddev')

  call mkglacier('mksrf_glacier.nc',var5d(:,:,1,1,1))
  where ( xmask < 0.5 )
    var5d(:,:,1,1,1) = vmisdat
  end where
  where (var5d(:,:,1,1,1) >= 0.0)
    pctspec(:,:) = pctspec(:,:) + var5d(:,:,1,1,1)
  end where
  istatus = nf90_put_var(ncid, iglcvar, var5d(:,:,1,1,1))
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write glacier')

  call mkwetland('mksrf_lanwat.nc',var5d(:,:,1,1,1),var5d(:,:,2,1,1))
  where ( xmask < 0.5 )
    var5d(:,:,1,1,1) = vmisdat
    var5d(:,:,2,1,1) = vmisdat
  end where
  where (var5d(:,:,1,1,1) >= 0.0)
    pctspec(:,:) = pctspec(:,:) + var5d(:,:,1,1,1)
  end where
  where (var5d(:,:,2,1,1) >= 0.0)
    pctspec(:,:) = pctspec(:,:) + var5d(:,:,2,1,1)
  end where
  istatus = nf90_put_var(ncid, iwetvar, var5d(:,:,1,1,1))
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write wetland')
  istatus = nf90_put_var(ncid, ilakevar, var5d(:,:,2,1,1))
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write lake')

  call mkurban('mksrf_urban.nc',var5d(:,:,1,1,1))
  where ( xmask < 0.5 )
    var5d(:,:,1,1,1) = vmisdat
  end where
  where (var5d(:,:,1,1,1) >= 0.0)
    pctspec(:,:) = pctspec(:,:) + var5d(:,:,1,1,1)
  end where
  istatus = nf90_put_var(ncid, iurbanvar, var5d(:,:,1,1,1))
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write urban')

  call mkpft('mksrf_pft.nc',var5d(:,:,1:npft,1,1))
  do np = 1 , npft
    where ( xmask < 0.5 )
      var5d(:,:,np,1,1) = vmisdat
    end where
  end do
  ! Here adjustment !
  do i = 1 , iy
    jloop: do j = 1 , jx
      if ( xmask(j,i) > 0.5 ) then
        if ( pctspec(j,i) > 99.99999 ) then
          var5d(j,i,:,1,1) = vmisdat
        else
          spft = 0.0
          if ( pctspec(j,i) > 0.00001 ) then
            do np = 1 , npft
              if ( var5d(j,i,np,1,1) > 0.0 ) then
                spft = spft + var5d(j,i,np,1,1) * 100.0/(100.0-pctspec(j,i))
              end if
            end do
          else 
            do np = 1 , npft
              if ( var5d(j,i,np,1,1) > 0.0 ) then
                spft = spft + var5d(j,i,np,1,1)
              end if
            end do
          end if
          if ( spft < 0.00001 ) then
            var5d(j,i,1,1,1) = 100.0-pctspec(j,i)
            var5d(j,i,2:npft,1,1) = 0.0
            cycle jloop
          end if
          diff = spft - 100.0
          if ( abs(diff) > 1.0E-5 ) then
            ! Normalize it !
            if ( abs(diff) > 0.0001 ) then
              ! Not worth doing if too small a diff...
              do np = 1 , npft
                if ( var5d(j,i,np,1,1) > 0.0 ) then
                  operat = dble(diff)*(dble(var5d(j,i,np,1,1))/dble(spft))
                  var5d(j,i,np,1,1) = var5d(j,i,np,1,1)-real(operat)
                  if ( var5d(j,i,np,1,1) < 0.0 ) then
                    var5d(j,i,np,1,1) = 0.0
                  end if
                end if
              end do
              ! Re-compute diff
              spft = 0.0
              if ( pctspec(j,i) > 0.00001 ) then
                do np = 1 , npft
                  if ( var5d(j,i,np,1,1) > 0.0 ) then
                    spft = spft + var5d(j,i,np,1,1) * 100.0/(100.0-pctspec(j,i))
                  end if
                end do
              else
                do np = 1 , npft
                  if ( var5d(j,i,np,1,1) > 0.0 ) then
                    spft = spft + var5d(j,i,np,1,1)
                  end if
                end do
              end if
              diff = spft - 100.0
            end if
            if ( abs(diff) > 1.0E-5 ) then
              pft_gt0 = (var5d(j,i,:,1,1) > diff/real(npft))
              where (pft_gt0)
                var5d(j,i,:,1,1) = var5d(j,i,:,1,1) - diff/real(count(pft_gt0))
              end where
            end if
          end if
        end if
      end if
    end do jloop
  end do
  istatus = nf90_put_var(ncid, ipftvar, var5d(:,:,1:npft,1,1))
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write pft')

  call mklaisai('mksrf_lai.nc',var5d(:,:,1:npft,1:nmon,1), &
                               var5d(:,:,1:npft,1:nmon,2), &
                               var5d(:,:,1:npft,1:nmon,3), &
                               var5d(:,:,1:npft,1:nmon,4))
  where (var5d(:,:,1:npft,1:nmon,1:4) < 0.0 )
    var5d(:,:,1:npft,1:nmon,1:4) = 0.0
  end where
  do nm = 1 , nmon
    do np = 1 , npft
      where ( xmask < 0.5 )
        var5d(:,:,np,nm,1) = vmisdat
        var5d(:,:,np,nm,2) = vmisdat
        var5d(:,:,np,nm,3) = vmisdat
        var5d(:,:,np,nm,4) = vmisdat
      end where
    end do
  end do
  do it = 1 , nmon
    do ip = 1 , npft
      istart(1) = 1
      icount(1) = ngcells
      istart(2) = ip
      icount(2) = 1
      istart(3) = it
      icount(3) = 1
      gcvar = pack(var5d(2:jx-2,2:iy-2,ip,it,1),lmask)
      istatus = nf90_put_var(ncid, isaivar, gcvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, 'Error write monthly_lai')
      gcvar = pack(var5d(2:jx-2,2:iy-2,ip,it,2),lmask)
      istatus = nf90_put_var(ncid, ilaivar, gcvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, 'Error write monthly_sai')
      gcvar = pack(var5d(2:jx-2,2:iy-2,ip,it,3),lmask)
      istatus = nf90_put_var(ncid, ivgtopvar, gcvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, 'Error write monthly_top')
      gcvar = pack(var5d(2:jx-2,2:iy-2,ip,it,4),lmask)
      istatus = nf90_put_var(ncid, ivgbotvar, gcvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, 'Error write monthly_bot')
    end do
  end do

  call mkfmax('mksrf_fmax.nc',var5d(:,:,1,1,1))
  where ( xmask < 0.5 )
    var5d(:,:,1,1,1) = vmisdat
  end where
  gcvar = pack(var5d(2:jx-2,2:iy-2,1,1,1),lmask)
  where ( gcvar < 0.0 ) gcvar = 0.05
  istatus = nf90_put_var(ncid, ifmaxvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write fmax')

  call mksoilcol('mksrf_soicol.nc',var5d(:,:,1,1,1))
  where ( xmask < 0.5 )
    var5d(:,:,1,1,1) = vmisdat
  end where
  gcvar = pack(var5d(2:jx-2,2:iy-2,1,1,1),lmask)
  istatus = nf90_put_var(ncid, isoilcolvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write soil color')

  call mksoitex('mksrf_soitex.nc',var5d(:,:,1:nsoil,1,1),var5d(:,:,1:nsoil,2,1))
  do il = 1 , nsoil
    where ( xmask < 0.5 )
      var5d(:,:,il,1,1) = vmisdat
      var5d(:,:,il,2,1) = vmisdat
    end where
  end do
  do il = 1 , nsoil
    istart(1) = 1
    icount(1) = ngcells
    istart(2) = il
    icount(2) = 1
    gcvar = pack(var5d(2:jx-2,2:iy-2,il,1,1),lmask)
    istatus = nf90_put_var(ncid, isandvar, gcvar,istart(1:2),icount(1:2))
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write sand pct')
    gcvar = pack(var5d(2:jx-2,2:iy-2,il,2,1),lmask)
    istatus = nf90_put_var(ncid, iclayvar, gcvar,istart(1:2),icount(1:2))
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write clay pct')
  end do

  call mkgdp('mksrf_gdp.nc',var5d(:,:,1,1,1))
  where ( xmask < 0.5 )
    var5d(:,:,1,1,1) = vmisdat
  end where
  gcvar = pack(var5d(2:jx-2,2:iy-2,1,1,1),lmask)
  istatus = nf90_put_var(ncid, igdpvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write gdp')

  call mkpeatf('mksrf_peatf.nc',var5d(:,:,1,1,1))
  where ( xmask < 0.5 )
    var5d(:,:,1,1,1) = vmisdat
  end where
  gcvar = pack(var5d(2:jx-2,2:iy-2,1,1,1),lmask)
  istatus = nf90_put_var(ncid, ipeatfvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write peatf')

  call mkabm('mksrf_abm.nc',var5d(:,:,1,1,1))
  where ( xmask < 0.5 )
    var5d(:,:,1,1,1) = vmisdat
  end where
  gcvar = pack(var5d(2:jx-2,2:iy-2,1,1,1),lmask)
  istatus = nf90_put_var(ncid, iabmvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write abm')

  call mkvocef('mksrf_vocef.nc',var5d(:,:,1:6,1,1))
  where ( xmask < 0.5 )
    var5d(:,:,1,1,1) = vmisdat
    var5d(:,:,2,1,1) = vmisdat
    var5d(:,:,3,1,1) = vmisdat
    var5d(:,:,4,1,1) = vmisdat
    var5d(:,:,5,1,1) = vmisdat
    var5d(:,:,6,1,1) = vmisdat
  end where
  gcvar = pack(var5d(2:jx-2,2:iy-2,1,1,1),lmask)
  istatus = nf90_put_var(ncid, ief_btrvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ef_btr')
  gcvar = pack(var5d(2:jx-2,2:iy-2,2,1,1),lmask)
  istatus = nf90_put_var(ncid, ief_crpvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ef_crp')
  gcvar = pack(var5d(2:jx-2,2:iy-2,3,1,1),lmask)
  istatus = nf90_put_var(ncid, ief_fdtvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ef_fdt')
  gcvar = pack(var5d(2:jx-2,2:iy-2,4,1,1),lmask)
  istatus = nf90_put_var(ncid, ief_fetvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ef_fet')
  gcvar = pack(var5d(2:jx-2,2:iy-2,5,1,1),lmask)
  istatus = nf90_put_var(ncid, ief_grsvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ef_grs')
  gcvar = pack(var5d(2:jx-2,2:iy-2,6,1,1),lmask)
  istatus = nf90_put_var(ncid, ief_shrvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ef_shr')

  call mkorganic('mksrf_organic.nc',var5d(:,:,1:nsoil,1,1))
  do il = 1 , nsoil
    where ( xmask < 0.5 )
      var5d(:,:,il,1,1) = vmisdat
    end where
  end do
  do il = 1 , nsoil
    istart(1) = 1
    icount(1) = ngcells
    istart(2) = il
    icount(2) = 1
    gcvar = pack(var5d(2:jx-2,2:iy-2,il,1,1),lmask)
    istatus = nf90_put_var(ncid, iorganicvar, gcvar,istart(1:2),icount(1:2))
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write organic')
  end do

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
