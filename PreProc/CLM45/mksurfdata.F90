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
  use mod_mklake
  use netcdf

  implicit none

  integer , parameter :: npft = 17
  integer , parameter :: nmon = 12
  integer , parameter :: nsoil = 10

  integer :: ngcells

  integer , parameter :: maxd3 = max(npft,nsoil)
  integer , parameter :: maxd4 = nmon

  real(rk8) , parameter :: vmisdat = -9999.0D0

  integer(ik4) :: istatus , ncid , ndim , nvar
  integer(ik4) , dimension(7) :: idims , ivdims
  integer(ik4) :: ivartime , iglcvar , iwetvar , ilakevar , iurbanvar
  integer(ik4) :: ipftvar , ilaivar , isaivar , ivgtopvar , ivgbotvar
  integer(ik4) :: ifmaxvar , isoilcolvar , isandvar , iclayvar
  integer(ik4) :: islopevar , istdvar , igdpvar , ipeatfvar , iabmvar
  integer(ik4) :: ief_btrvar , ief_crpvar , ief_fdtvar , ief_fetvar
  integer(ik4) :: ief_grsvar , ief_shrvar , iorganicvar , idepthvar
  integer(ik4) :: ilndvar , iscvar , ilatvar , ilonvar , itopovar
  integer(ik4) :: ijxvar , iiyvar
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
  integer(ik4) , dimension(1) :: istart1 , icount1 , mxsoil_color , iloc
  real(rk8) :: spft , mean
  real(rk8) :: operat , diff
  integer(ik4) :: ierr
  integer(ik4) :: i , j , n , ip , il , np , nm , it , ipnt
  integer(ik4) :: jgstart , jgstop
  character(len=256) :: namelistfile , prgname
  character(len=256) :: terfile , outfile
  character(len=64) :: csdate
  real(rk8) , dimension(:,:) , pointer :: pctspec , pctslake
  real(rk8) , pointer , dimension(:,:,:,:,:) :: var5d
  real(rk8) , pointer , dimension(:) :: gcvar
  integer(ik4) , pointer , dimension(:) :: iiy , ijx
  integer , pointer , dimension(:) :: landpoint
  logical , dimension(npft) :: pft_gt0

  call get_command_argument(0,value=prgname)
  call get_command_argument(1,value=namelistfile)
  call initparam(namelistfile, ierr)
  if ( ierr/=0 ) then
    write(stderr,*) 'Parameter initialization not completed'
    write(stderr,*) 'Usage : '
    write(stderr,*) '          ', trim(prgname), ' regcm.in'
    write(stderr,*) ' '
    call die(__FILE__,'Check argument and namelist syntax',__LINE__)
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
  rxlat = real(xlat)
  rxlon = real(xlon)
  rsigx = real(sigx)
  call setup_pack()
  if ( i_band == 1 ) then
    jgstart = 1
    jgstop = jxsg
  else
    jgstart = 2
    jgstop = jxsg-2
  end if
  ngcells = count(xmask(jgstart:jgstop,2:iysg-2) > 0.5D0)
  call closefile(ncid)

  ! Open Output in NetCDF format

  outfile = trim(dirglob)//pthsep//trim(domname)//'_CLM45_surface.nc'

  call createfile_withname(outfile,ncid)
  call add_common_global_params(ncid,'mksurfdata',.false.)
  ndim = 1
  call define_basic_dimensions(ncid,jxsg,iysg,kzp1,ndim,idims)
  call add_dimension(ncid,'time',nf90_unlimited,ndim,idims)
  call add_dimension(ncid,'lsmpft',npft,ndim,idims)
  call add_dimension(ncid,'nlevsoi',nsoil,ndim,idims)
  call add_dimension(ncid,'gridcell',ngcells,ndim,idims)
  call define_horizontal_coord(ncid,jxsg,iysg,xjx,yiy,idims,ihvar)
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
  istatus = nf90_def_var(ncid, 'mxsoil_color', nf90_int, iscvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var mxsoil_color')

  istatus = nf90_def_var(ncid, 'gridcell', nf90_int, idims(7), ilndvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var gridcell')
  istatus = nf90_put_att(ncid, ilndvar, 'compress','xlat xlon')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add gridcell compress')

  istatus = nf90_def_var(ncid, 'xclon', nf90_double, idims(7), ilonvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var xclon')
  istatus = nf90_put_att(ncid, ilndvar, 'compress','xlat xlon')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add xclon')

  istatus = nf90_def_var(ncid, 'yclat', nf90_double, idims(7), ilatvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var yclat')
  istatus = nf90_put_att(ncid, ilndvar, 'compress','xlat xlon')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add yclat')

  istatus = nf90_def_var(ncid, 'topo', nf90_double, idims(7), itopovar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var topo')
  istatus = nf90_put_att(ncid, ilndvar, 'compress','xlat xlon')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add topo')

  istatus = nf90_def_var(ncid, 'jx2d', nf90_int, idims(7), ijxvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var jx2d')
  istatus = nf90_put_att(ncid, ilndvar, 'compress','xlat xlon')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add jx2d')

  istatus = nf90_def_var(ncid, 'iy2d', nf90_int, idims(7), iiyvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var iy2d')
  istatus = nf90_put_att(ncid, ilndvar, 'compress','xlat xlon')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add iy2d')

  istatus = nf90_def_var(ncid, 'SLOPE', nf90_double, idims(7), islopevar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var slope')
  istatus = nf90_put_att(ncid, islopevar, 'long_name','Elevation slope')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add slope long_name')
  istatus = nf90_put_att(ncid, islopevar, 'units','degree')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add slope units')
  istatus = nf90_def_var(ncid, 'STD_ELEV', nf90_double, idims(7), istdvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var stddev')
  istatus = nf90_put_att(ncid, istdvar, 'long_name','Elevation std dev')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add stddev long_name')
  istatus = nf90_put_att(ncid, istdvar, 'units','m')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add stddev units')

  istatus = nf90_def_var(ncid, 'PCT_GLACIER', nf90_double, idims(7), iglcvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var glacier')
  istatus = nf90_put_att(ncid, iglcvar, 'long_name','percent glacier')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add glacier long_name')
  istatus = nf90_put_att(ncid, iglcvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add glacier units')

  istatus = nf90_def_var(ncid, 'PCT_WETLAND', nf90_double, idims(7), iwetvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var wetland')
  istatus = nf90_put_att(ncid, iwetvar, 'long_name','percent wetland')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add wetland long_name')
  istatus = nf90_put_att(ncid, iwetvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add wetland units')

  istatus = nf90_def_var(ncid, 'PCT_LAKE', nf90_double, idims(7), ilakevar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var lake')
  istatus = nf90_put_att(ncid, ilakevar, 'long_name','percent lake')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add lake long_name')
  istatus = nf90_put_att(ncid, ilakevar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add lake units')

  istatus = nf90_def_var(ncid, 'PCT_URBAN', nf90_double, idims(7), iurbanvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var urban')
  istatus = nf90_put_att(ncid, iurbanvar, 'long_name','percent urban')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add urban long_name')
  istatus = nf90_put_att(ncid, iurbanvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add urban units')

  ivdims(1) = idims(7)
  ivdims(2) = idims(5)
  istatus = nf90_def_var(ncid, 'PCT_PFT', nf90_double, ivdims(1:2), ipftvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var pft')
  istatus = nf90_put_att(ncid, ipftvar, 'long_name','percent pft')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add pft long_name')
  istatus = nf90_put_att(ncid, ipftvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add pft units')
  istatus = nf90_put_att(ncid, ipftvar, '_FillValue',vmisdat)
  call checkncerr(istatus,__FILE__,__LINE__,'Error add pft Fill Value')

  ivdims(1) = idims(7)
  ivdims(2) = idims(5)
  ivdims(3) = idims(4)
  istatus = nf90_def_var(ncid, 'MONTHLY_LAI', nf90_double, ivdims(1:3), ilaivar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var monthly_lai')
  istatus = nf90_put_att(ncid, ilaivar, 'long_name','monthly leaf area index')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_lai long_name')
  istatus = nf90_put_att(ncid, ilaivar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_lai units')
  istatus = nf90_def_var(ncid, 'MONTHLY_SAI', nf90_double, ivdims(1:3), isaivar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var monthly_sai')
  istatus = nf90_put_att(ncid, isaivar, 'long_name','monthly stem area index')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_sai long_name')
  istatus = nf90_put_att(ncid, isaivar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_sai units')
  istatus = nf90_def_var(ncid, 'MONTHLY_HEIGHT_TOP', nf90_double, &
          ivdims(1:3), ivgtopvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var monthly_top')
  istatus = nf90_put_att(ncid, ivgtopvar, 'long_name','monthly height top')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_top long_name')
  istatus = nf90_put_att(ncid, ivgtopvar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_top units')
  istatus = nf90_def_var(ncid, 'MONTHLY_HEIGHT_BOT', nf90_double, &
          ivdims(1:3), ivgbotvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var monthly_bot')
  istatus = nf90_put_att(ncid, ivgbotvar, 'long_name','monthly height bottom')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_bot long_name')
  istatus = nf90_put_att(ncid, ivgbotvar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_bot units')

  istatus = nf90_def_var(ncid, 'FMAX', nf90_double, idims(7), ifmaxvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var fmax')
  istatus = nf90_put_att(ncid, ifmaxvar, 'long_name', &
          'maximum fractional saturated area at 1/8 degree')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add fmax long_name')
  istatus = nf90_put_att(ncid, ifmaxvar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add fmax units')

  istatus = nf90_def_var(ncid, 'SOIL_COLOR', nf90_double, idims(7),isoilcolvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var soilcol')
  istatus = nf90_put_att(ncid, isoilcolvar, 'long_name', &
          'maximum fractional saturated area at 1/8 degree')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add soilcol long_name')
  istatus = nf90_put_att(ncid, isoilcolvar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add soilcol units')

  istatus = nf90_def_var(ncid, 'gdp', nf90_double, idims(7),igdpvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var gdp')
  istatus = nf90_put_att(ncid, igdpvar, 'long_name', 'real GDP')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add gdp long_name')
  istatus = nf90_put_att(ncid, igdpvar, 'units','K 1995US$/capita')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add gdp units')

  istatus = nf90_def_var(ncid, 'peatf', nf90_double, idims(7),ipeatfvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var peatf')
  istatus = nf90_put_att(ncid, ipeatfvar, 'long_name', &
          'global fraction of peatland')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add peatf long_name')
  istatus = nf90_put_att(ncid, ipeatfvar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add peatf units')

  istatus = nf90_def_var(ncid, 'abm', nf90_double, idims(7),iabmvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var abm')
  istatus = nf90_put_att(ncid, iabmvar, 'long_name', &
          'peak month for agri fire')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add abm long_name')
  istatus = nf90_put_att(ncid, iabmvar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add abm units')

  istatus = nf90_def_var(ncid, 'EF1_BTR', nf90_double, idims(7),ief_btrvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ef_btr')
  istatus = nf90_put_att(ncid, ief_btrvar, 'long_name', &
          'broadleaf tree emission factor')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_btr long_name')
  istatus = nf90_put_att(ncid,ief_btrvar,'units','micrograms isoprene m-2 h-1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_btr units')
  istatus = nf90_def_var(ncid, 'EF1_CRP', nf90_double, idims(7),ief_crpvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ef_crp')
  istatus = nf90_put_att(ncid, ief_crpvar, 'long_name', &
          'crop emission factor')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_crp long_name')
  istatus = nf90_put_att(ncid,ief_crpvar,'units','micrograms isoprene m-2 h-1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_crp units')
  istatus = nf90_def_var(ncid, 'EF1_FDT', nf90_double, idims(7),ief_fdtvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ef_fdt')
  istatus = nf90_put_att(ncid, ief_fdtvar, 'long_name', &
          'Fineleaf deciduous tree emission factor')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_fdt long_name')
  istatus = nf90_put_att(ncid,ief_fdtvar,'units','micrograms isoprene m-2 h-1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_fdt units')
  istatus = nf90_def_var(ncid, 'EF1_FET', nf90_double, idims(7),ief_fetvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ef_fet')
  istatus = nf90_put_att(ncid, ief_fetvar, 'long_name', &
          'Fineleaf evergreen tree emission factor')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_fet long_name')
  istatus = nf90_put_att(ncid,ief_fetvar,'units','micrograms isoprene m-2 h-1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_fet units')
  istatus = nf90_def_var(ncid, 'EF1_GRS', nf90_double, idims(7),ief_grsvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ef_grs')
  istatus = nf90_put_att(ncid, ief_grsvar, 'long_name', &
          'grass, non-vascular plants and other ground cover emission factor')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_grs long_name')
  istatus = nf90_put_att(ncid,ief_grsvar,'units','micrograms isoprene m-2 h-1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_grs units')
  istatus = nf90_def_var(ncid, 'EF1_SHR', nf90_double, idims(7),ief_shrvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ef_shr')
  istatus = nf90_put_att(ncid, ief_shrvar, 'long_name', &
          'shrub emission factor')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_shr long_name')
  istatus = nf90_put_att(ncid,ief_shrvar,'units','micrograms isoprene m-2 h-1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_shr units')

  ivdims(1) = idims(7)
  ivdims(2) = idims(6)
  istatus = nf90_def_var(ncid, 'PCT_SAND', nf90_double, ivdims(1:2), isandvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var sand')
  istatus = nf90_put_att(ncid, isandvar, 'long_name','percent sand')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add sand long_name')
  istatus = nf90_put_att(ncid, isandvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add sand units')
  istatus = nf90_def_var(ncid, 'PCT_CLAY', nf90_double, ivdims(1:2), iclayvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var clay')
  istatus = nf90_put_att(ncid, iclayvar, 'long_name','percent clay')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add clay long_name')
  istatus = nf90_put_att(ncid, iclayvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add clay units')

  istatus = nf90_def_var(ncid, 'ORGANIC', nf90_double, ivdims(1:2), iorganicvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var organic')
  istatus = nf90_put_att(ncid, iorganicvar, 'long_name', &
          'organic soil density at soil levels')
  istatus = nf90_put_att(ncid, iorganicvar, 'comment', &
          'assumed carbon content 0.58 gC per gOM')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add organic long_name')
  istatus = nf90_put_att(ncid, iorganicvar, 'units','kg m-3')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add organic units')

  istatus = nf90_def_var(ncid, 'LAKEDEPTH',nf90_double,idims(7),idepthvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var lake')
  istatus = nf90_put_att(ncid, idepthvar, 'long_name', &
          'Lake Depth (default = 10m)')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add lake long_name')
  istatus = nf90_put_att(ncid, idepthvar, 'units','m')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add lake units')

  istatus = nf90_enddef(ncid)
  call checkncerr(istatus,__FILE__,__LINE__,'Error exit define mode')

  hptop = real(ptop*10.0D0)
  call write_vertical_coord(ncid,rsigx,hptop,izvar)
  call write_horizontal_coord(ncid,xjx,yiy,ihvar)
  ipnt = 1
  call write_var2d_static(ncid,'xlat',rxlat,ipnt,illvar)
  call write_var2d_static(ncid,'xlon',rxlon,ipnt,illvar)

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

  mxsoil_color(1) = 20
  istatus = nf90_put_var(ncid, iscvar, mxsoil_color)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write mxsoil_color')

  call getmem2d(pctspec,1,jxsg,1,iysg,'mksurfdata: pctspec')
  call getmem2d(pctslake,1,jxsg,1,iysg,'mksurfdata: pctslake')
  call getmem1d(gcvar,1,ngcells,'mksurfdata: gcvar')
  call getmem1d(iiy,1,ngcells,'mksurfdata: iiy')
  call getmem1d(ijx,1,ngcells,'mksurfdata: ijx')
  call getmem1d(landpoint,1,ngcells,'mksurfdata: landpoint')
  call getmem5d(var5d,1,jxsg,1,iysg,1,maxd3,1,maxd4,1,4,'mksurfdata: var5d')
  pctspec(:,:) = 0.0D0
  pctslake(:,:) = 0.0D0
  ip = 1
  do i = 2 , iysg-2
    do j = jgstart , jgstop
      if ( xmask(j,i) > 0.5D0 ) then
        landpoint(ip) = (i-1)*jxsg+j
        ip = ip + 1
      end if
    end do
  end do
  istatus = nf90_put_var(ncid, ilndvar, landpoint)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write gridcell')
  call mypack(xlat,gcvar)
  istatus = nf90_put_var(ncid, ilatvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write yclat')
  call mypack(xlon,gcvar)
  istatus = nf90_put_var(ncid, ilonvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write xclon')
  call mypack(topo,gcvar)
  istatus = nf90_put_var(ncid, itopovar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write topo')
  ip = 1
  do i = 2 , iysg-2
    do j = jgstart , jgstop
      if ( xmask(j,i) > 0.5D0 ) then
        iiy(ip) = i
        ijx(ip) = j
        ip = ip + 1
      end if
    end do
  end do
  istatus = nf90_put_var(ncid, iiyvar, iiy)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write iiy')
  istatus = nf90_put_var(ncid, ijxvar, ijx)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ijx')

  var5d(:,:,1,1,1) = 0.0D0
  ! Calculate slope and std
  do i = 2 , iysg-2
    do j = jgstart , jgstop
      var5d(j,i,1,1,1) = &
          atan((sum(topo(j-1:j+1,i-1:i+1)-topo(j,i))/8.0D0)/(ds*1000.0D0))
      mean = sum(topo(j-1:j+1,i-1:i+1))/9.0D0
      var5d(j,i,2,1,1) = sqrt(sum((topo(j-1:j+1,i-1:i+1)-mean)**2)/8.0D0)
    end do
  end do
  where ( xmask < 0.5D0 )
    var5d(:,:,1,1,1) = vmisdat
    var5d(:,:,2,1,1) = vmisdat
  end where
  call mypack(var5d(:,:,1,1,1),gcvar)
  gcvar = gcvar*real(raddeg)
  istatus = nf90_put_var(ncid, islopevar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write slope')
  call mypack(var5d(:,:,2,1,1),gcvar)
  istatus = nf90_put_var(ncid, istdvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write stddev')

  call mkglacier('mksrf_glacier.nc',var5d(:,:,1,1,1))
  call mkwetland('mksrf_lanwat.nc',var5d(:,:,2,1,1),var5d(:,:,3,1,1))
  call mkurban('mksrf_urban.nc',var5d(:,:,4,1,1))
  var5d(:,:,1:4,1,1) = max(var5d(:,:,1:4,1,1),0.0D0)
  where ( xmask < 0.5D0 )
    var5d(:,:,1,1,1) = vmisdat
    var5d(:,:,2,1,1) = vmisdat
    var5d(:,:,3,1,1) = vmisdat
    var5d(:,:,4,1,1) = vmisdat
  end where

  where ( var5d(:,:,1,1,1) > 100.0D0 ) var5d(:,:,1,1,1) = 100.0D0
  where ( var5d(:,:,2,1,1) > 100.0D0 ) var5d(:,:,2,1,1) = 100.0D0
  where ( var5d(:,:,3,1,1) > 100.0D0 ) var5d(:,:,3,1,1) = 100.0D0
  where ( var5d(:,:,4,1,1) > 100.0D0 ) var5d(:,:,4,1,1) = 100.0D0

  do i = 2 , iysg-2
    do j = jgstart , jgstop
      if ( xmask(j,i) > 0.5D0 ) then
        pctspec(j,i) = var5d(j,i,1,1,1) + var5d(j,i,2,1,1) + &
                       var5d(j,i,3,1,1) + var5d(j,i,4,1,1)
        if ( pctspec(j,i) < 0.0D0 ) then
          write(stderr,*) 'Negative pctspec ',pctspec(j,i),' at j,i ', j , i
          call die(__FILE__,'PCTSPEC error',__LINE__)
        end if
        if ( pctspec(j,i) > 99.9D0 ) then
          diff = 100.0D0 - pctspec(j,i)
          iloc = maxloc(var5d(j,i,:,1,1))
          var5d(j,i,iloc(1),1,1) = var5d(j,i,iloc(1),1,1) + diff
          pctspec(j,i) = var5d(j,i,1,1,1) + var5d(j,i,2,1,1) + &
                         var5d(j,i,3,1,1) + var5d(j,i,4,1,1)
          if ( pctspec(j,i) /= 100.0D0 ) then
            diff = 100.0D0 - pctspec(j,i)
            iloc = maxloc(var5d(j,i,:,1,1))
            var5d(j,i,iloc(1),1,1) = var5d(j,i,iloc(1),1,1) + diff
            pctspec(j,i) = var5d(j,i,1,1,1) + var5d(j,i,2,1,1) + &
                           var5d(j,i,3,1,1) + var5d(j,i,4,1,1)
            if ( pctspec(j,i) /= 100.0D0 ) then
              write(stderr,*) 'Cannot normalize pctspec at j,i ', j , i
              call die(__FILE__,'PCTSPEC normalization error',__LINE__)
            end if
          end if
        end if
      end if
    end do
  end do
  call mypack(var5d(:,:,1,1,1),gcvar)
  istatus = nf90_put_var(ncid, iglcvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write glacier')
  call mypack(var5d(:,:,2,1,1),gcvar)
  istatus = nf90_put_var(ncid, iwetvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write wetland')
  call mypack(var5d(:,:,3,1,1),gcvar)
  istatus = nf90_put_var(ncid, ilakevar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write lake')
  call mypack(var5d(:,:,4,1,1),gcvar)
  istatus = nf90_put_var(ncid, iurbanvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write urban')
  where ( var5d(:,:,3,1,1) > 0.0D0 )
    pctslake = var5d(:,:,3,1,1)
  end where

  call mkpft('mksrf_pft.nc',var5d(:,:,1:npft,1,1))
  var5d(:,:,1:npft,1,1) = max(var5d(:,:,1:npft,1,1),0.0D0)
  do np = 1 , npft
    where ( xmask < 0.5D0 )
      var5d(:,:,np,1,1) = vmisdat
    end where
  end do
  ! Here adjustment !
  do i = 2 , iysg-2
    jloop: do j = jgstart , jgstop
      if ( xmask(j,i) > 0.5D0 ) then
        if ( pctspec(j,i) > 99.99999D0 ) then
          var5d(j,i,:,1,1) = vmisdat
        else
          spft = 0.0D0
          if ( pctspec(j,i) > 0.00001D0 ) then
            do np = 1 , npft
              if ( var5d(j,i,np,1,1) > 0.0D0 ) then
                spft = spft + var5d(j,i,np,1,1) * 100.0D0/(100.0D0-pctspec(j,i))
              end if
            end do
          else 
            do np = 1 , npft
              if ( var5d(j,i,np,1,1) > 0.0D0 ) then
                spft = spft + var5d(j,i,np,1,1)
              end if
            end do
          end if
          if ( spft < 0.00001D0 ) then
            ! Substitute with something around it
            call bestaround(var5d(:,:,:,1,1),i,j)
            spft = 0.0D0
            if ( pctspec(j,i) > 0.00001D0 ) then
              do np = 1 , npft
                if ( var5d(j,i,np,1,1) > 0.0D0 ) then
                  spft = spft + var5d(j,i,np,1,1) * &
                          100.0D0/(100.0D0-pctspec(j,i))
                end if
              end do
            else 
              do np = 1 , npft
                if ( var5d(j,i,np,1,1) > 0.0D0 ) then
                  spft = spft + var5d(j,i,np,1,1)
                end if
              end do
            end if
            if ( spft < 0.00001D0 ) then
              call die(__FILE__,'No points around !',__LINE__)
            end if
          end if
          diff = spft - 100.0D0
          if ( abs(diff) > 1.0D-5 ) then
            ! Normalize it !
            if ( abs(diff) > 0.0001D0 ) then
              ! Not worth doing if too small a diff...
              do np = 1 , npft
                if ( var5d(j,i,np,1,1) > 0.0D0 ) then
                  operat = diff*(var5d(j,i,np,1,1)/spft)
                  var5d(j,i,np,1,1) = var5d(j,i,np,1,1) - operat
                  if ( var5d(j,i,np,1,1) < 0.0D0 ) then
                    var5d(j,i,np,1,1) = 0.0D0
                  end if
                end if
              end do
              ! Re-compute diff
              spft = 0.0D0
              if ( pctspec(j,i) > 0.00001D0 ) then
                do np = 1 , npft
                  if ( var5d(j,i,np,1,1) > 0.0D0 ) then
                    spft = spft + var5d(j,i,np,1,1) * &
                            100.0D0/(100.0D0-pctspec(j,i))
                  end if
                end do
              else
                do np = 1 , npft
                  if ( var5d(j,i,np,1,1) > 0.0D0 ) then
                    spft = spft + var5d(j,i,np,1,1)
                  end if
                end do
              end if
              diff = spft - 100.0D0
            end if
            if ( abs(diff) > 1.0D-5 ) then
              pft_gt0 = (var5d(j,i,:,1,1) > diff/dble(npft))
              do n = 1 , npft
                if ( pft_gt0(n) .and. count(pft_gt0) > 0 ) then
                  var5d(j,i,n,1,1) = (var5d(j,i,n,1,1) - &
                     diff/dble(count(pft_gt0)))
                end if
              end do
            end if
          end if
        end if
      end if
    end do jloop
  end do
  do ip = 1 , npft
    istart(1) = 1
    icount(1) = ngcells
    istart(2) = ip
    icount(2) = 1
    call mypack(var5d(:,:,ip,1,1),gcvar)
    istatus = nf90_put_var(ncid, ipftvar, gcvar, istart(1:2), icount(1:2))
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write pft')
  end do

  call mklaisai('mksrf_lai.nc',var5d(:,:,1:npft,1:nmon,1), &
                               var5d(:,:,1:npft,1:nmon,2), &
                               var5d(:,:,1:npft,1:nmon,3), &
                               var5d(:,:,1:npft,1:nmon,4))
  where (var5d(:,:,1:npft,1:nmon,1:4) < 0.0D0 )
    var5d(:,:,1:npft,1:nmon,1:4) = 0.0D0
  end where
  do nm = 1 , nmon
    do np = 1 , npft
      where ( xmask < 0.5D0 )
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
      call mypack(var5d(:,:,ip,it,1),gcvar)
      istatus = nf90_put_var(ncid, isaivar, gcvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, 'Error write monthly_lai')
      call mypack(var5d(:,:,ip,it,2),gcvar)
      istatus = nf90_put_var(ncid, ilaivar, gcvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, 'Error write monthly_sai')
      call mypack(var5d(:,:,ip,it,3),gcvar)
      istatus = nf90_put_var(ncid, ivgtopvar, gcvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, 'Error write monthly_top')
      call mypack(var5d(:,:,ip,it,4),gcvar)
      istatus = nf90_put_var(ncid, ivgbotvar, gcvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, 'Error write monthly_bot')
    end do
  end do

  call mkfmax('mksrf_fmax.nc',var5d(:,:,1,1,1))
  where ( xmask < 0.5D0 )
    var5d(:,:,1,1,1) = vmisdat
  end where
  call mypack(var5d(:,:,1,1,1),gcvar)
  where ( gcvar < 0.0D0 ) gcvar = 0.05D0
  istatus = nf90_put_var(ncid, ifmaxvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write fmax')

  call mksoilcol('mksrf_soicol.nc',var5d(:,:,1,1,1))
  where ( xmask < 0.5D0 )
    var5d(:,:,1,1,1) = vmisdat
  end where
  call mypack(var5d(:,:,1,1,1),gcvar)
  if ( any(gcvar < 0.0D0) ) call fillvar(gcvar)
  istatus = nf90_put_var(ncid, isoilcolvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write soil color')

  call mksoitex('mksrf_soitex.nc',var5d(:,:,1:nsoil,1,1),var5d(:,:,1:nsoil,2,1))
  where ( var5d(:,:,:,1,1) < 0.0D0 )
    var5d(:,:,:,1,1) = 50.0D0
    var5d(:,:,:,2,1) = 50.0D0
  end where
  do il = 1 , nsoil
    where ( xmask < 0.5D0 )
      var5d(:,:,il,1,1) = vmisdat
      var5d(:,:,il,2,1) = vmisdat
    end where
  end do
  do il = 1 , nsoil
    istart(1) = 1
    icount(1) = ngcells
    istart(2) = il
    icount(2) = 1
    call mypack(var5d(:,:,il,1,1),gcvar)
    istatus = nf90_put_var(ncid, isandvar, gcvar,istart(1:2),icount(1:2))
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write sand pct')
    call mypack(var5d(:,:,il,2,1),gcvar)
    istatus = nf90_put_var(ncid, iclayvar, gcvar,istart(1:2),icount(1:2))
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write clay pct')
  end do

  call mkgdp('mksrf_gdp.nc',var5d(:,:,1,1,1))
  where ( xmask < 0.5D0 )
    var5d(:,:,1,1,1) = vmisdat
  end where
  call mypack(var5d(:,:,1,1,1),gcvar)
  if ( any(gcvar < 0.0D0) ) call fillvar(gcvar)
  istatus = nf90_put_var(ncid, igdpvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write gdp')

  call mkpeatf('mksrf_peatf.nc',var5d(:,:,1,1,1))
  where ( var5d(:,:,1,1,1) < 0.0D0 )
   var5d(:,:,1,1,1) = 0.0D0
  end where
  where ( xmask < 0.5D0 )
    var5d(:,:,1,1,1) = vmisdat
  end where
  call mypack(var5d(:,:,1,1,1),gcvar)
  if ( any(gcvar < 0.0D0) ) call fillvar(gcvar)
  istatus = nf90_put_var(ncid, ipeatfvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write peatf')

  call mkabm('mksrf_abm.nc',var5d(:,:,1,1,1))
  where ( xmask < 0.5D0 )
    var5d(:,:,1,1,1) = vmisdat
  end where
  call mypack(var5d(:,:,1,1,1),gcvar)
  if ( any(gcvar < 0.0D0) ) call fillvar(gcvar)
  istatus = nf90_put_var(ncid, iabmvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write abm')

  call mkvocef('mksrf_vocef.nc',var5d(:,:,1:6,1,1))
  where ( var5d(:,:,1,1,1) < 0.0D0 )
    var5d(:,:,1,1,1) = 0.0D0
    var5d(:,:,2,1,1) = 0.0D0
    var5d(:,:,3,1,1) = 0.0D0
    var5d(:,:,4,1,1) = 0.0D0
    var5d(:,:,5,1,1) = 0.0D0
    var5d(:,:,6,1,1) = 0.0D0
  end where
  where ( xmask < 0.5D0 )
    var5d(:,:,1,1,1) = vmisdat
    var5d(:,:,2,1,1) = vmisdat
    var5d(:,:,3,1,1) = vmisdat
    var5d(:,:,4,1,1) = vmisdat
    var5d(:,:,5,1,1) = vmisdat
    var5d(:,:,6,1,1) = vmisdat
  end where
  call mypack(var5d(:,:,1,1,1),gcvar)
  istatus = nf90_put_var(ncid, ief_btrvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ef_btr')
  call mypack(var5d(:,:,2,1,1),gcvar)
  istatus = nf90_put_var(ncid, ief_crpvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ef_crp')
  call mypack(var5d(:,:,3,1,1),gcvar)
  istatus = nf90_put_var(ncid, ief_fdtvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ef_fdt')
  call mypack(var5d(:,:,4,1,1),gcvar)
  istatus = nf90_put_var(ncid, ief_fetvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ef_fet')
  call mypack(var5d(:,:,5,1,1),gcvar)
  istatus = nf90_put_var(ncid, ief_grsvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ef_grs')
  call mypack(var5d(:,:,6,1,1),gcvar)
  istatus = nf90_put_var(ncid, ief_shrvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ef_shr')

  call mkorganic('mksrf_organic.nc',var5d(:,:,1:nsoil,1,1))
  where ( var5d(:,:,1:nsoil,1,1) < 0.0D0 )
    var5d(:,:,1:nsoil,1,1) = 1.0D0
  end where
  do il = 1 , nsoil
    where ( xmask < 0.5D0 )
      var5d(:,:,il,1,1) = vmisdat
    end where
  end do
  do il = 1 , nsoil
    istart(1) = 1
    icount(1) = ngcells
    istart(2) = il
    icount(2) = 1
    call mypack(var5d(:,:,il,1,1),gcvar)
    istatus = nf90_put_var(ncid, iorganicvar, gcvar,istart(1:2),icount(1:2))
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write organic')
  end do

  call mklake('mksrf_lake.nc',var5d(:,:,1,1,1))
  where ( var5d(:,:,1,1,1) < 10.0D0 )
    var5d(:,:,1,1,1) = 10.0D0
  end where
  call mypack(var5d(:,:,1,1,1),gcvar)
  istatus = nf90_put_var(ncid, idepthvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write lake')

  istatus = nf90_close(ncid)
  call checkncerr(istatus,__FILE__,__LINE__,  &
    'Error close file '//trim(outfile))

  call memory_destroy

  write(stdout,*) 'Successfully completed CLM preprocessing.'

  contains

  subroutine fillvar(xvar)
    implicit none
    real(rk8) , dimension(:) , intent(inout) :: xvar
    integer(ik4) :: nvar , i , ip , left , right
    real(rk8) :: leftval , rightval
    nvar = size(xvar)
    do i = 1 , nvar
      if ( xvar(i) < 0.0D0 ) then
        leftval = -1.0D0
        rightval = -1.0D0
        ip = 1
        do
          if ( leftval < 0.0D0 ) left = i - ip
          if ( left < 1 ) left = np - left
          if ( rightval < 0.0D0 ) right = i + ip
          if ( right > nvar ) right = right - nvar
          if ( xvar(left) > 0.0D0 ) leftval = xvar(left)
          if ( xvar(right) > 0.0D0 ) rightval = xvar(right)
          if ( leftval > 0.0D0 .and. rightval > 0.0D0 ) then
            xvar(i) = (leftval+rightval)/2.0D0
          end if
          if ( ip == nvar / 2 ) then
            call die(__FILE__,'Not finding anything around !',__LINE__)
            exit
          end if
          if ( xvar(i) > 0.0D0 ) exit
          ip = ip + 1
        end do
      end if
    end do
  end subroutine fillvar

  subroutine bestaround(pft,i,j)
    implicit none
    real(rk8) , dimension(:,:,:) , intent(inout) :: pft
    integer(ik4) , intent(in) :: i , j
    real(rk8) , dimension (npft,(2*minval(shape(pft(:,:,1)))/2+1)**2) :: vals
    integer(ik4) :: ii , jj , js , is , ip , n , il , maxil
    il = 1
    maxil = minval(shape(pft(:,:,1)))/2
    do
      ip = 0
      vals(:,:) = 0.0D0
      do ii = i - il , i + il
        do jj = j - il , j + il
          is = ii
          js = jj
          if ( js < 1 ) js = 1-js
          if ( js > jxsg ) js = 2*jxsg - js
          if ( is < 1 ) is = 1-is
          if ( is > iysg ) is = 2*iysg - is
          do n = 1 , npft
            if ( pft(js,is,n) > 0.00001D0 ) then
              ip = ip + 1
              vals(n,ip) = vals(n,ip) + pft(js,is,n)
            end if
          end do
        end do
      end do
      if ( ip > 0 ) then
        do n = 1 , npft
          pft(j,i,n) = sum(vals(n,1:ip))/real(ip)
        end do
        exit
      else
        il = il + 1
        if ( il == maxil ) then
          write(stderr,*) 'At point lat = ',xlat(j,i)
          write(stderr,*) '         lon = ',xlon(j,i)
          call die(__FILE__,'Not finding anything around !',__LINE__)
          exit
        end if
      end if
    end do
  end subroutine bestaround

  recursive subroutine sortpatch(vals,svals,ird,lsub)
    implicit none
    real(rk8) , dimension(:) , intent(in) :: vals
    real(rk8) , dimension(:) , intent(inout) :: svals
    integer(ik4) , dimension(:) , intent(inout) :: ird
    logical , optional :: lsub
    integer(ik4) :: i , iswap
    real(rk8) :: rswap
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
