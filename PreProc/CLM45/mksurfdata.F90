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
 
#ifndef CLM45
  use mod_stdio
  write(stdout,*) 'RegCM built without CLM45 support.'
  write(stdout,*) 'Please recompile it using --enable-clm45 flag'
#else

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
#ifdef CN
  use mod_mklightning
  use mod_mkpopd
#ifdef DYNPFT
  use mod_mkharvest
  use mod_mkdynpft
#endif
#ifdef LCH4
  use mod_mklch4
#endif
#endif
#ifdef VICHYDRO
  use mod_mkvic
#endif

  implicit none

  integer(ik4) :: npft
  integer(ik4) , parameter :: nlurb = 5
  integer(ik4) , parameter :: nsoil = 10

  integer(ik4) , parameter :: nmon = 12
  integer(ik4) , parameter :: nrad = 2
  integer(ik4) , parameter :: nsol = 2
  integer(ik4) , parameter :: numurbl = 3

#ifdef CN
  integer(ik4) , parameter :: noleap_yday_3h = 365*8
  integer(ik4) , parameter :: nyears = 2100-1850+1
#ifdef DYNPFT
  integer(ik4) :: y1 , y2 , mon , day , hour
  character(len=4) :: cy
  character(len=32) :: p1 , p2
#endif
#endif

  integer(ik4) :: ngcells

  real(rk8) , parameter :: vmisdat = -9999.0D0

  integer(ik4) :: istatus , ncid , ndim , nvar
  integer(ik4) , dimension(12) :: idims , ivdims
  integer(ik4) :: ivartime , iglcvar , iwetvar , ilakevar , iurbanvar
  integer(ik4) :: ipftvar , ilaivar , isaivar , ivgtopvar , ivgbotvar
  integer(ik4) :: ifmaxvar , isoilcolvar , isandvar , iclayvar
  integer(ik4) :: islopevar , istdvar , igdpvar , ipeatfvar , iabmvar
  integer(ik4) :: ief_btrvar , ief_crpvar , ief_fdtvar , ief_fetvar
  integer(ik4) :: ief_grsvar , ief_shrvar , iorganicvar , idepthvar
  integer(ik4) :: ilndvar , iscvar , ilatvar , ilonvar , itopovar
  integer(ik4) , dimension(npu2d) :: iurb2d
  integer(ik4) , dimension(npu3d) :: iurb3d
  integer(ik4) , dimension(npu4d*2) :: iurb4d
#ifdef CN
  integer(ik4) :: ilightning , ipopden
#ifdef DYNPFT
  integer(ik4) :: iharvvh1 , iharvvh2 , iharvsh1 , iharvsh2 , iharvsh3 , igraz
  character(len=256) :: dynfile
#endif
#ifdef LCH4
  integer(ik4) :: if0 , ip3 , izwt0
#endif
#endif
#ifdef VICHYDRO
  integer(ik4) :: ibinfl , ids , idsmax ,iws
#endif
  integer(ik4) :: ijxvar , iiyvar
  type(rcm_time_and_date) :: irefdate , imondate
  type(rcm_time_interval) :: tdif
  real(rk4) , pointer , dimension(:) :: yiy
  real(rk4) , pointer , dimension(:) :: xjx
  real(rk4) :: hptop
  real(rk8) :: smallnum
  real(rk8) , dimension(1) :: xdate
  integer(ik4) , dimension(3) :: istart , icount
  integer(ik4) , dimension(2) :: ihvar
  integer(ik4) , dimension(2) :: illvar
  integer(ik4) , dimension(2) :: izvar
  integer(ik4) , dimension(1) :: istart1 , icount1 , mxsoil_color , iloc
  real(rk8) :: spft , mean
  real(rk8) :: operat , diff
  integer(ik4) :: ierr
  integer(ik4) :: i , j , n , ip , il , ir , iu , np , nm , it , ipnt , iurbmax
  integer(ik4) :: jgstart , jgstop , igstart , igstop
  character(len=256) :: namelistfile , prgname
  character(len=256) :: terfile , outfile
  character(len=64) :: csdate , pftfile , laifile
  real(rk8) , dimension(:,:) , pointer :: pctspec , pctslake
  real(rk8) , pointer , dimension(:,:) :: var2d
  real(rk8) , pointer , dimension(:,:,:) :: var3d
  real(rk8) , pointer , dimension(:,:,:,:) :: var4d
  real(rk8) , pointer , dimension(:,:,:,:,:) :: var5d
  real(rk8) , pointer , dimension(:,:,:,:,:,:) :: var6d
  real(rk8) , pointer , dimension(:) :: gcvar
  integer(ik4) , pointer , dimension(:) :: iiy , ijx
  integer(ik4) , pointer , dimension(:) :: landpoint
  logical , pointer , dimension(:) :: pft_gt0
  logical :: subgrid

  smallnum = dble(epsilon(1.0))/4.0D0

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

  if ( .not. enable_more_crop_pft ) then
#ifdef CROP
#ifdef DYNPFT
    call die(__FILE__, &
            'Dynamic pft and CROP are not permitted together',__LINE__)
#endif
    npft = 25
    pftfile = 'mksrf_24pft.nc'
    laifile = 'mksrf_24lai.nc'
#else
    npft = 17
#ifdef DYNPFT
    call split_idate(globidate1,y1,mon,day,hour)
    p1 = 'dynamic'
    p2 = '.'
    write(cy,'(i0.4)') y1
    if ( y1 > 2005 ) then
      select case (dattyp(4:5))
        case ('RF')
          continue
        case ('26')
          p2 = 'SCENARIO'//pthsep//'RCP2.6'
        case ('45')
          p2 = 'SCENARIO'//pthsep//'RCP4.5'
        case ('60')
          p2 = 'SCENARIO'//pthsep//'RCP6.0'
        case ('85')
          p2 = 'SCENARIO'//pthsep//'RCP8.5'
        case default
          call die(__FILE__,'Dynamic landuse only supported for CMIP5',__LINE__)
      end select
    end if
    pftfile = trim(p1)//pthsep//trim(p2)//pthsep//'mksrf_landuse_'//cy//'.nc'
#else
    pftfile = 'mksrf_pft.nc'
#endif
    laifile = 'mksrf_lai.nc'
#endif
  else
#ifdef DYNPFT
    call die(__FILE__, &
            'Dynamic pft and enable_more_crop_pft are not permitted together', &
            __LINE__)
#endif
    write(stdout,*) 'This will fail if you do not modify the model code!'
    write(stdout,*) 'Edit Main/clmlib/clm4.5/mod_clm_varpar.F90'
    write(stdout,*) 'and search for ADD_MORE_CROP_PFT comment!'
    npft = 21
    pftfile = 'mksrf_pft_crop.nc'
    laifile = 'mksrf_lai_crop.nc'
  end if
  call getmem1d(pft_gt0,1,npft,'pft_gt0')

  call init_domain
  !
  ! Get latitudes, longitudes and mask from DOMAIN file
  !
  if ( nsg > 1 ) then
    write (terfile,'(a,i0.3,a)') &
         trim(dirter)//pthsep//trim(domname)//'_DOMAIN',nsg,'.nc'
    subgrid = .true.
  else
    terfile = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
    subgrid = .false.
  end if
  call openfile_withname(terfile,ncid)
  call read_domain(ncid,sigx,xlat,xlon,ht=topo,mask=xmask,lsubgrid=subgrid)
  rxlat = real(xlat)
  rxlon = real(xlon)
  rsigx = real(sigx)
  if ( i_band == 1 ) then
    jgstart = 1
    jgstop = jxsg
    igstart = nsg+1
    igstop = iysg-2*nsg
    call setup_pack(1,jx,2,iy-2)
  else
    jgstart = nsg+1
    jgstop = jxsg-2*nsg
    igstart = nsg+1
    igstop = iysg-2*nsg
    call setup_pack(2,jx-2,2,iy-2)
  end if
  ngcells = count(xmask(jgstart:jgstop,igstart:igstop) > 0.5D0)
  call closefile(ncid)

  ! Open Output in NetCDF format

  outfile = trim(dirglob)//pthsep//trim(domname)//'_CLM45_surface.nc'

  call createfile_withname(outfile,ncid)
  call add_common_global_params(ncid,'mksurfdata',.false.)
  ndim = 1
  call define_basic_dimensions(ncid,jxsg,iysg,kzp1,ndim,idims)
  call add_dimension(ncid,'nmon',nmon,ndim,idims)
  call add_dimension(ncid,'lsmpft',npft,ndim,idims)
  call add_dimension(ncid,'nlevsoi',nsoil,ndim,idims)
  call add_dimension(ncid,'gridcell',ngcells,ndim,idims)
  call add_dimension(ncid,'numurbl',numurbl,ndim,idims)
  call add_dimension(ncid,'nlevurb',nlurb,ndim,idims)
  call add_dimension(ncid,'numrad',nrad,ndim,idims)
#ifdef CN
  call add_dimension(ncid,'noleap_3h',noleap_yday_3h,ndim,idims)
  call add_dimension(ncid,'year',nyears,ndim,idims)
#endif
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

  ivdims(1) = idims(7)
  ivdims(2) = idims(8)
  istatus = nf90_def_var(ncid, 'PCT_URBAN', nf90_double, ivdims(1:2), iurbanvar)
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

  ivdims(1) = idims(7)
  ivdims(2) = idims(8)
  do i = 1 , npu2d
    istatus = nf90_def_var(ncid, parm2d(i), nf90_double, ivdims(1:2), iurb2d(i))
    call checkncerr(istatus,__FILE__,__LINE__, 'Error add var')
    istatus = nf90_put_att(ncid, iurb2d(i), 'long_name', lngn2d(i))
    call checkncerr(istatus,__FILE__,__LINE__,'Error add long_name')
    istatus = nf90_put_att(ncid, iurb2d(i), 'units', unit2d(i) )
    call checkncerr(istatus,__FILE__,__LINE__,'Error add units')
    istatus = nf90_put_att(ncid, iurb2d(i), '_FillValue',vmisdat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error add _FillValue')
  end do

  ivdims(1) = idims(7)
  ivdims(2) = idims(8)
  ivdims(3) = idims(9)
  do i = 1 , npu3d
    istatus = nf90_def_var(ncid, parm3d(i), nf90_double, ivdims(1:3), iurb3d(i))
    call checkncerr(istatus,__FILE__,__LINE__, 'Error add var')
    istatus = nf90_put_att(ncid, iurb3d(i), 'long_name', lngn3d(i))
    call checkncerr(istatus,__FILE__,__LINE__,'Error add long_name')
    istatus = nf90_put_att(ncid, iurb3d(i), 'units', unit3d(i) )
    call checkncerr(istatus,__FILE__,__LINE__,'Error add units')
    istatus = nf90_put_att(ncid, iurb3d(i), '_FillValue',vmisdat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error add _FillValue')
  end do

  ivdims(1) = idims(7)
  ivdims(2) = idims(8)
  ivdims(3) = idims(10)
  do i = 1 , npu4d
    istatus = nf90_def_var(ncid, trim(parm4d(i))//'_DIR', &
      nf90_double, ivdims(1:3), iurb4d(2*i-1))
    call checkncerr(istatus,__FILE__,__LINE__, 'Error add var')
    istatus = nf90_put_att(ncid, iurb4d(2*i-1), 'long_name', lngn4d(i))
    call checkncerr(istatus,__FILE__,__LINE__,'Error add long_name')
    istatus = nf90_put_att(ncid, iurb4d(2*i-1), 'units', unit4d(i) )
    call checkncerr(istatus,__FILE__,__LINE__,'Error add units')
    istatus = nf90_put_att(ncid, iurb4d(2*i-1), '_FillValue',vmisdat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error add _FillValue')
    istatus = nf90_def_var(ncid, trim(parm4d(i))//'_DIF', &
      nf90_double, ivdims(1:3), iurb4d(2*i))
    call checkncerr(istatus,__FILE__,__LINE__, 'Error add var')
    istatus = nf90_put_att(ncid, iurb4d(2*i), 'long_name', lngn4d(i))
    call checkncerr(istatus,__FILE__,__LINE__,'Error add long_name')
    istatus = nf90_put_att(ncid, iurb4d(2*i), 'units', unit4d(i) )
    call checkncerr(istatus,__FILE__,__LINE__,'Error add units')
    istatus = nf90_put_att(ncid, iurb4d(2*i), '_FillValue',vmisdat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error add _FillValue')
  end do

#ifdef CN
  ivdims(1) = idims(7)
  ivdims(2) = idims(11)
  istatus = nf90_def_var(ncid, 'LNFM',nf90_double,ivdims(1:2),ilightning)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var lnfm')
  istatus = nf90_put_att(ncid, ilightning, 'long_name', &
          'Lightning frequency')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add lnfm long_name')
  istatus = nf90_put_att(ncid, ilightning, 'units','counts/km^2/hr')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add lnfm units')
  istatus = nf90_put_att(ncid, ilightning, '_FillValue',vmisdat)
  call checkncerr(istatus,__FILE__,__LINE__,'Error add _FillValue')
  ivdims(1) = idims(7)
  ivdims(2) = idims(12)
  istatus = nf90_def_var(ncid, 'HDM',nf90_double,ivdims(1:2),ipopden)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var hdm')
  istatus = nf90_put_att(ncid, ipopden, 'long_name', &
          'human population density')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add hdm long_name')
  istatus = nf90_put_att(ncid, ipopden, 'units','counts/km^2')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add hdm units')
#ifdef LCH4
  istatus = nf90_def_var(ncid, 'F0',nf90_double,idims(7),if0)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var F0')
  istatus = nf90_def_var(ncid, 'P3',nf90_double,idims(7),ip3)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var P3')
  istatus = nf90_def_var(ncid, 'ZWT0',nf90_double,idims(7),izwt0)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ZWT0')
#endif
#endif

#ifdef VICHYDRO
  istatus = nf90_def_var(ncid, 'binfl',nf90_double,idims(7),ibinfl)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var binfl')
  istatus = nf90_put_att(ncid, ibinfl, 'long_name', &
          'VIC b parameter for the Variable Infiltration Capacity Curve')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add binfl long_name')
  istatus = nf90_put_att(ncid, ibinfl, 'units','unitless')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add binfl units')
  istatus = nf90_def_var(ncid, 'Ds',nf90_double,idims(7),ids)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ds')
  istatus = nf90_put_att(ncid, ids, 'long_name', &
          'VIC Ds parameter for the ARNO curve')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ds long_name')
  istatus = nf90_put_att(ncid, ids, 'units','unitless')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ds units')
  istatus = nf90_def_var(ncid, 'Dsmax',nf90_double,idims(7),idsmax)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var dsmax')
  istatus = nf90_put_att(ncid, idsmax, 'long_name', &
          'VIC Dsmax parameter for the ARNO curve')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add dsmax long_name')
  istatus = nf90_put_att(ncid, idsmax, 'units','unitless')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add dsmax units')
  istatus = nf90_def_var(ncid, 'Ws',nf90_double,idims(7),iws)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ws')
  istatus = nf90_put_att(ncid, iws, 'long_name', &
          'VIC Ws parameter for the ARNO curve')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ws long_name')
  istatus = nf90_put_att(ncid, iws, 'units','unitless')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ws units')
#endif

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
  pctspec(:,:) = 0.0D0
  pctslake(:,:) = 0.0D0
  ip = 1
  do i = igstart , igstop
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
  do i = igstart , igstop
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

  allocate(var3d(jxsg,iysg,2))
  var3d(:,:,1) = 0.0D0
  ! Calculate slope and std
  do i = igstart , igstop
    do j = jgstart , jgstop
      var3d(j,i,1) = &
          atan((sum(topo(j-1:j+1,i-1:i+1)-topo(j,i))/8.0D0)/(ds*1000.0D0))
      mean = sum(topo(j-1:j+1,i-1:i+1))/9.0D0
      var3d(j,i,2) = sqrt(sum((topo(j-1:j+1,i-1:i+1)-mean)**2)/8.0D0)
    end do
  end do
  where ( xmask < 0.5D0 )
    var3d(:,:,1) = vmisdat
    var3d(:,:,2) = vmisdat
  end where
  call mypack(var3d(:,:,1),gcvar)
  gcvar = gcvar*real(raddeg)
  istatus = nf90_put_var(ncid, islopevar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write slope')
  call mypack(var3d(:,:,2),gcvar)
  istatus = nf90_put_var(ncid, istdvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write stddev')
  deallocate(var3d)

  iurbmax = numurbl+3
  allocate(var3d(jxsg,iysg,iurbmax))
  call mkglacier('mksrf_glacier.nc',var3d(:,:,1))
  call mkwetland('mksrf_lanwat.nc',var3d(:,:,2),var3d(:,:,3))
  call mkurban_base('mksrf_urban.nc',var3d(:,:,4:iurbmax))
  if ( .not. enable_urban_landunit ) then
    write (stderr,*) 'Disable URBAN Areas in CLM4.5 Model !'
    var3d(:,:,4:iurbmax) = 0.0D0
  end if
  var3d = max(var3d,0.0D0)
  do i = 1 , iurbmax
    where ( xmask < 0.5D0 )
      var3d(:,:,i) = vmisdat
    end where
  end do

  where ( var3d(:,:,1) > 100.0D0 ) var3d(:,:,1) = 100.0D0
  where ( var3d(:,:,2) > 100.0D0 ) var3d(:,:,2) = 100.0D0
  where ( var3d(:,:,3) > 100.0D0 ) var3d(:,:,3) = 100.0D0
  where ( var3d(:,:,4:iurbmax) > 100.0D0 ) var3d(:,:,4:iurbmax) = 100.0D0
  var3d(:,:,:) = dint(var3d(:,:,:))

  do i = igstart , igstop
    do j = jgstart , jgstop
      if ( xmask(j,i) > 0.5D0 ) then
        pctspec(j,i) = sum(var3d(j,i,:))
        if ( pctspec(j,i) < 0.0D0 ) then
          write(stderr,*) 'Negative pctspec ',pctspec(j,i),' at j,i ', j , i
          call die(__FILE__,'PCTSPEC error',__LINE__)
        end if
        if ( pctspec(j,i) > (200.0D0 - smallnum) ) then
          var3d(j,i,:) = var3d(j,i,:) / (pctspec(j,i)/100.0D0)
          pctspec(j,i) = sum(var3d(j,i,:))
        end if
        if ( pctspec(j,i) > (100.0D0 - smallnum) ) then
          diff = 100.0D0 - pctspec(j,i)
          iloc = maxloc(var3d(j,i,:))
          var3d(j,i,iloc(1)) = var3d(j,i,iloc(1)) + diff
          pctspec(j,i) = sum(var3d(j,i,:))
          if ( pctspec(j,i) /= 100.0D0 ) then
            diff = 100.0D0 - pctspec(j,i)
            iloc = maxloc(var3d(j,i,:))
            var3d(j,i,iloc(1)) = var3d(j,i,iloc(1)) + diff
            pctspec(j,i) = sum(var3d(j,i,:))
            if ( pctspec(j,i) /= 100.0D0 ) then
              write(stderr,*) 'Cannot normalize pctspec at j,i ', &
                      j , i , pctspec(j,i)
              call die(__FILE__,'PCTSPEC normalization error',__LINE__)
            end if
          end if
        end if
        if ( pctspec(j,i) < smallnum ) then
          pctspec(j,i) = 0.0D0
          var3d(j,i,:) = 0.0D0
        end if
      end if
    end do
  end do
  call mypack(var3d(:,:,1),gcvar)
  istatus = nf90_put_var(ncid, iglcvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write glacier')
  call mypack(var3d(:,:,2),gcvar)
  istatus = nf90_put_var(ncid, iwetvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write wetland')
  call mypack(var3d(:,:,3),gcvar)
  istatus = nf90_put_var(ncid, ilakevar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write lake')
  do i = 1 , numurbl
    istart(1) = 1
    icount(1) = ngcells
    istart(2) = i
    icount(2) = 1
    call mypack(var3d(:,:,3+i),gcvar)
    istatus = nf90_put_var(ncid, iurbanvar, gcvar, istart(1:2), icount(1:2))
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write urban')
  end do
  where ( var3d(:,:,3) > 0.0D0 )
    pctslake = var3d(:,:,3)
  end where
  deallocate(var3d)

  allocate(var3d(jxsg,iysg,npft))
  call mkpft(pftfile,var3d(:,:,:))
  var3d(:,:,:) = dint(max(var3d(:,:,:),0.0D0))
  do np = 1 , npft
    where ( xmask < 0.5D0 )
      var3d(:,:,np) = vmisdat
    end where
  end do

  ! Here adjustment !
  do i = igstart , igstop
    do j = jgstart , jgstop
      if ( xmask(j,i) > 0.5D0 ) then
        if ( pctspec(j,i) > (100.0D0 - smallnum) ) then
          var3d(j,i,:) = 0.0D0
        else
          spft = 0.0D0
          do np = 1 , npft
            if ( var3d(j,i,np) > 0.0D0 ) then
              spft = spft + (var3d(j,i,np)*100.0D0)/(100.0D0-pctspec(j,i))
            end if
          end do
          if ( spft < smallnum ) then
            ! Substitute with something around it
            call bestaround_pft(var3d,i,j)
            spft = 0.0D0
            do np = 1 , npft
              if ( var3d(j,i,np) > 0.0D0 ) then
                spft = spft + (var3d(j,i,np) * 100.0D0)/(100.0D0-pctspec(j,i))
              end if
            end do
            if ( spft < smallnum ) then
              call die(__FILE__,'No points around !',__LINE__)
            end if
          end if
          diff = spft - 100.0D0
          if ( abs(diff) > smallnum ) then
            ! Normalize it !
            do np = 1 , npft
              if ( var3d(j,i,np) > 0.0D0 ) then
                operat = diff*(var3d(j,i,np)/spft)
                var3d(j,i,np) = var3d(j,i,np) - operat
                if ( var3d(j,i,np) < 0.0D0 ) then
                  var3d(j,i,np) = 0.0D0
                end if
              end if
            end do
            ! Re-compute diff
            spft = 0.0D0
            do np = 1 , npft
              if ( var3d(j,i,np) > 0.0D0 ) then
                spft = spft + (var3d(j,i,np) * 100.0D0)/(100.0D0-pctspec(j,i))
              end if
            end do
            diff = spft - 100.0D0
            if ( abs(diff) > smallnum ) then
              pft_gt0 = (var3d(j,i,:) > diff/dble(npft))
              do n = 1 , npft
                if ( pft_gt0(n) .and. count(pft_gt0) > 0 ) then
                  var3d(j,i,n) = (var3d(j,i,n) - diff/dble(count(pft_gt0)))
                end if
              end do
              spft = 0.0D0
              do np = 1 , npft
                if ( var3d(j,i,np) > 0.0D0 ) then
                  spft = spft + (var3d(j,i,np) * 100.0D0)/(100.0D0-pctspec(j,i))
                end if
              end do
              diff = spft - 100.0D0
              if ( abs(diff) > smallnum ) then
                do np = 1 , npft
                  if ( var3d(j,i,np) > 0.0D0 ) then
                    operat = diff*(var3d(j,i,np)/spft)
                    var3d(j,i,np) = var3d(j,i,np) - operat
                    if ( var3d(j,i,np) < 0.0D0 ) then
                      var3d(j,i,np) = 0.0D0
                    end if
                  end if
                end do
                spft = 0.0D0
                do np = 1 , npft
                  if ( var3d(j,i,np) > 0.0D0 ) then
                    spft = spft + &
                            (var3d(j,i,np) * 100.0D0)/(100.0D0-pctspec(j,i))
                  end if
                end do
                diff = spft - 100.0D0
                if ( abs(diff) > smallnum ) then
                  write (0,*) abs(diff)
                  call die(__FILE__,'Cannot normalize...',__LINE__)
                end if
              end if
            end if
          end if
        end if
      end if
    end do
  end do
  do ip = 1 , npft
    istart(1) = 1
    icount(1) = ngcells
    istart(2) = ip
    icount(2) = 1
    call mypack(var3d(:,:,ip),gcvar)
    istatus = nf90_put_var(ncid, ipftvar, gcvar, istart(1:2), icount(1:2))
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write pft')
  end do
  deallocate(var3d)

  allocate(var5d(jxsg,iysg,npft,nmon,4))
  call mklaisai(laifile,var5d(:,:,:,:,1), var5d(:,:,:,:,2), &
                        var5d(:,:,:,:,3), var5d(:,:,:,:,4))
  where (var5d(:,:,:,:,:) < 0.0D0 )
    var5d(:,:,:,:,:) = 0.0D0
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
  deallocate(var5d)

  allocate(var2d(jxsg,iysg))
  call mkfmax('mksrf_fmax.nc',var2d)
  where ( xmask < 0.5D0 )
    var2d = vmisdat
  end where
  call mypack(var2d,gcvar)
  where ( gcvar < 0.0D0 ) gcvar = 0.05D0
  istatus = nf90_put_var(ncid, ifmaxvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write fmax')
  deallocate(var2d)

  allocate(var2d(jxsg,iysg))
  call mksoilcol('mksrf_soicol.nc',var2d)
  where ( xmask < 0.5D0 )
    var2d = vmisdat
  end where
  call mypack(var2d,gcvar)
  if ( any(gcvar < 0.0D0) ) call fillvar(gcvar)
  istatus = nf90_put_var(ncid, isoilcolvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write soil color')
  deallocate(var2d)

  allocate(var4d(jxsg,iysg,nsoil,2))
  call mksoitex('mksrf_soitex.nc',var4d(:,:,:,1),var4d(:,:,:,2))
  where ( var4d(:,:,:,1) < 0.0D0 )
    var4d(:,:,:,1) = 50.0D0
    var4d(:,:,:,2) = 50.0D0
  end where
  do il = 1 , nsoil
    where ( xmask < 0.5D0 )
      var4d(:,:,il,1) = vmisdat
      var4d(:,:,il,2) = vmisdat
    end where
  end do
  do il = 1 , nsoil
    istart(1) = 1
    icount(1) = ngcells
    istart(2) = il
    icount(2) = 1
    call mypack(var4d(:,:,il,1),gcvar)
    istatus = nf90_put_var(ncid, isandvar, gcvar,istart(1:2),icount(1:2))
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write sand pct')
    call mypack(var4d(:,:,il,2),gcvar)
    istatus = nf90_put_var(ncid, iclayvar, gcvar,istart(1:2),icount(1:2))
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write clay pct')
  end do
  deallocate(var4d)

  allocate(var2d(jxsg,iysg))
  call mkgdp('mksrf_gdp.nc',var2d)
  where ( xmask < 0.5D0 )
    var2d = vmisdat
  end where
  call mypack(var2d,gcvar)
  if ( any(gcvar < 0.0D0) ) call fillvar(gcvar)
  istatus = nf90_put_var(ncid, igdpvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write gdp')
  deallocate(var2d)

  allocate(var2d(jxsg,iysg))
  call mkpeatf('mksrf_peatf.nc',var2d)
  where ( var2d < 0.0D0 )
   var2d = 0.0D0
  end where
  where ( xmask < 0.5D0 )
    var2d = vmisdat
  end where
  call mypack(var2d,gcvar)
  if ( any(gcvar < 0.0D0) ) call fillvar(gcvar)
  istatus = nf90_put_var(ncid, ipeatfvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write peatf')
  deallocate(var2d)

  allocate(var2d(jxsg,iysg))
  call mkabm('mksrf_abm.nc',var2d)
  where ( xmask < 0.5D0 )
    var2d = vmisdat
  end where
  call mypack(var2d,gcvar)
  if ( any(gcvar < 0.0D0) ) call fillvar(gcvar)
  istatus = nf90_put_var(ncid, iabmvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write abm')
  deallocate(var2d)

  allocate(var3d(jxsg,iysg,6))
  call mkvocef('mksrf_vocef.nc',var3d)
  where ( var3d(:,:,1) < 0.0D0 )
    var3d(:,:,1) = 0.0D0
    var3d(:,:,2) = 0.0D0
    var3d(:,:,3) = 0.0D0
    var3d(:,:,4) = 0.0D0
    var3d(:,:,5) = 0.0D0
    var3d(:,:,6) = 0.0D0
  end where
  where ( xmask < 0.5D0 )
    var3d(:,:,1) = vmisdat
    var3d(:,:,2) = vmisdat
    var3d(:,:,3) = vmisdat
    var3d(:,:,4) = vmisdat
    var3d(:,:,5) = vmisdat
    var3d(:,:,6) = vmisdat
  end where
  call mypack(var3d(:,:,1),gcvar)
  istatus = nf90_put_var(ncid, ief_btrvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ef_btr')
  call mypack(var3d(:,:,2),gcvar)
  istatus = nf90_put_var(ncid, ief_crpvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ef_crp')
  call mypack(var3d(:,:,3),gcvar)
  istatus = nf90_put_var(ncid, ief_fdtvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ef_fdt')
  call mypack(var3d(:,:,4),gcvar)
  istatus = nf90_put_var(ncid, ief_fetvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ef_fet')
  call mypack(var3d(:,:,5),gcvar)
  istatus = nf90_put_var(ncid, ief_grsvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ef_grs')
  call mypack(var3d(:,:,6),gcvar)
  istatus = nf90_put_var(ncid, ief_shrvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ef_shr')
  deallocate(var3d)

  allocate(var3d(jxsg,iysg,nsoil))
  call mkorganic('mksrf_organic.nc',var3d)
  where ( var3d < 0.0D0 )
    var3d = 1.0D0
  end where
  do il = 1 , nsoil
    where ( xmask < 0.5D0 )
      var3d(:,:,il) = vmisdat
    end where
  end do
  do il = 1 , nsoil
    istart(1) = 1
    icount(1) = ngcells
    istart(2) = il
    icount(2) = 1
    call mypack(var3d(:,:,il),gcvar)
    istatus = nf90_put_var(ncid, iorganicvar, gcvar,istart(1:2),icount(1:2))
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write organic')
  end do
  deallocate(var3d)

  allocate(var2d(jxsg,iysg))
  call mklake('mksrf_lake.nc',var2d)
  where ( var2d < 10.0D0 )
    var2d = 10.0D0
  end where
  call mypack(var2d,gcvar)
  istatus = nf90_put_var(ncid, idepthvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write lake')
  deallocate(var2d)

  allocate(var4d(jxsg,iysg,numurbl,npu2d))
  allocate(var5d(jxsg,iysg,nlurb,numurbl,npu3d))
  allocate(var6d(jxsg,iysg,nsol,nrad,numurbl,npu4d))
  var4d = vmisdat
  var5d = vmisdat
  var6d = vmisdat
  call mkurban_param('mksrf_urban.nc',var4d,var5d,var6d)
  istart(1) = 1
  icount(1) = ngcells
  do iu = 1 , numurbl
    istart(2) = iu
    icount(2) = 1
    do i = 1 , npu2d
      call mypack(var4d(:,:,iu,ip2d(parm2d(i))),gcvar)
      istatus = nf90_put_var(ncid, iurb2d(i), gcvar, istart(1:2), icount(1:2))
      call checkncerr(istatus,__FILE__,__LINE__, 'Error write '//parm2d(i))
    end do
  end do
  do iu = 1 , numurbl
    istart(2) = iu
    icount(2) = 1
    do il = 1 , nlurb
      istart(3) = il
      icount(3) = 1
      do i = 1 , npu3d
        call mypack(var5d(:,:,il,iu,ip3d(parm3d(i))),gcvar)
        where ( gcvar < 0.0D0 ) gcvar = vmisdat
        istatus = nf90_put_var(ncid, iurb3d(i), gcvar, istart(1:3), icount(1:3))
        call checkncerr(istatus,__FILE__,__LINE__, 'Error write '//parm3d(i))
      end do
    end do
  end do
  do iu = 1 , numurbl
    istart(2) = iu
    icount(2) = 1
    do ir = 1 , nrad
      istart(3) = ir
      icount(3) = 1
      ! Here we assume (NO other information in input file), that DIRECT
      ! albedo is solar 1, and DIF is solar 2. Hope this is correct. As
      ! of now they are equal in the input file, so no fear. MAY change,
      ! so this comment is here as a reminder. In case, just switch the
      ! iurb4d indexes in the two nf90_put_var calls.
      do i = 1 , npu4d
        call mypack(var6d(:,:,1,ir,iu,ip4d(parm4d(i))),gcvar)
        where ( gcvar < 0.0D0 ) gcvar = vmisdat
        istatus = nf90_put_var(ncid, iurb4d(2*i-1), &
                gcvar, istart(1:3), icount(1:3))
        call checkncerr(istatus,__FILE__,__LINE__, 'Error write '//parm4d(i))
        call mypack(var6d(:,:,2,ir,iu,ip4d(parm4d(i))),gcvar)
        where ( gcvar < 0.0D0 ) gcvar = vmisdat
        istatus = nf90_put_var(ncid, iurb4d(2*i), &
                gcvar, istart(1:3), icount(1:3))
        call checkncerr(istatus,__FILE__,__LINE__, 'Error write '//parm4d(i))
      end do
    end do
  end do
  deallocate(var4d)
  deallocate(var5d)
  deallocate(var6d)

#ifdef CN
  allocate(var2d(jxsg,iysg))
  istart(1) = 1
  icount(1) = ngcells
  do it = 1 , noleap_yday_3h
    call mklightning('mksrf_lightning.nc',var2d,it)
    where ( xmask < 0.5D0 )
      var2d = vmisdat
    end where
    call mypack(var2d,gcvar)
    istart(2) = it
    icount(2) = 1
    istatus = nf90_put_var(ncid, ilightning, gcvar, istart(1:2), icount(1:2))
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write lnfm')
  end do
  do it = 1 , nyears
    call mkpopd('mksrf_popd.nc',var2d,it)
    where ( xmask < 0.5D0 )
      var2d = vmisdat
    end where
    call mypack(var2d,gcvar)
    istart(2) = it
    icount(2) = 1
    istatus = nf90_put_var(ncid, ipopden, gcvar, istart(1:2), icount(1:2))
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write hdm')
  end do
  deallocate(var2d)
#ifdef LCH4
  allocate(var3d(jxsg,iysg,3))
  call mklch4('mksrf_ch4inversion.nc',var3d)
  do i = 1 , 3
    where ( xmask < 0.5D0 )
      var3d(:,:,i) = vmisdat
    end where
  end do
  call mypack(var3d(:,:,1),gcvar)
  where ( gcvar < 0.0 )
    gcvar = 0.01D0
  end where
  istatus = nf90_put_var(ncid, if0, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write F0')
  call mypack(var3d(:,:,2),gcvar)
  where ( gcvar < 0.0 )
    gcvar = 10.0D0
  end where
  istatus = nf90_put_var(ncid, ip3, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write P3')
  call mypack(var3d(:,:,3),gcvar)
  where ( gcvar < 0.0 )
    gcvar = 0.01D0
  end where
  istatus = nf90_put_var(ncid, izwt0, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ZWt0')
  deallocate(var3d)
#endif
#endif

#ifdef VICHYDRO
  allocate(var3d(jxsg,iysg,4))
  call mkvic('mksrf_vic.nc',var3d)
  do i = 1 , 4
    where ( xmask < 0.5D0 )
      var3d(:,:,i) = vmisdat
    end where
  end do
  call mypack(var3d(:,:,1),gcvar)
  where ( gcvar < 0.0 )
    gcvar = 0.1D0
  end where
  istatus = nf90_put_var(ncid, ibinfl, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write binfl')
  call mypack(var3d(:,:,2),gcvar)
  where ( gcvar < 0.0 )
    gcvar = 0.1D0
  end where
  istatus = nf90_put_var(ncid, ids, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ds')
  call mypack(var3d(:,:,3),gcvar)
  where ( gcvar < 0.0 )
    gcvar = 10.0D0
  end where
  istatus = nf90_put_var(ncid, idsmax, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write dsmax')
  call mypack(var3d(:,:,4),gcvar)
  where ( gcvar < 0.0 )
    gcvar = 0.75D0
  end where
  istatus = nf90_put_var(ncid, iws, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ws')
  deallocate(var3d)
#endif

  istatus = nf90_close(ncid)
  call checkncerr(istatus,__FILE__,__LINE__,  &
    'Error close file '//trim(outfile))

#ifdef CN
#ifdef DYNPFT
  call split_idate(globidate1,y1,mon,day,hour)
  call split_idate(globidate2,y2,mon,day,hour)
  do it = y1 - 1 , y2 + 1
    write (cy,'(i0.4)') it
    dynfile = trim(dirglob)//pthsep//trim(domname)//'_CLM45_surface_'//cy//'.nc'
    call createfile_withname(dynfile,ncid)
    call add_common_global_params(ncid,'mksurfdata',.false.)
    ndim = 1
    idims(:) = -1
    call define_basic_dimensions(ncid,jxsg,iysg,kzp1,ndim,idims)
    call add_dimension(ncid,'lsmpft',npft,ndim,idims)
    call add_dimension(ncid,'gridcell',ngcells,ndim,idims)
    call define_horizontal_coord(ncid,jxsg,iysg,xjx,yiy,idims,ihvar)
    call define_vertical_coord(ncid,idims,izvar)

    nvar = 1
    call define_cross_geolocation_coord(ncid,idims,nvar,illvar)

    ! Variables
    istatus = nf90_def_var(ncid, 'gridcell', nf90_int, idims(5), ilndvar)
    call checkncerr(istatus,__FILE__,__LINE__,  'Error add var gridcell')
    istatus = nf90_put_att(ncid, ilndvar, 'compress','xlat xlon')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add gridcell compress')

    istatus = nf90_def_var(ncid, 'xclon', nf90_double, idims(5), ilonvar)
    call checkncerr(istatus,__FILE__,__LINE__,  'Error add var xclon')
    istatus = nf90_put_att(ncid, ilndvar, 'compress','xlat xlon')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add xclon')

    istatus = nf90_def_var(ncid, 'yclat', nf90_double, idims(5), ilatvar)
    call checkncerr(istatus,__FILE__,__LINE__,  'Error add var yclat')
    istatus = nf90_put_att(ncid, ilndvar, 'compress','xlat xlon')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add yclat')

    istatus = nf90_def_var(ncid, 'jx2d', nf90_int, idims(5), ijxvar)
    call checkncerr(istatus,__FILE__,__LINE__,  'Error add var jx2d')
    istatus = nf90_put_att(ncid, ilndvar, 'compress','xlat xlon')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add jx2d')

    istatus = nf90_def_var(ncid, 'iy2d', nf90_int, idims(5), iiyvar)
    call checkncerr(istatus,__FILE__,__LINE__,  'Error add var iy2d')
    istatus = nf90_put_att(ncid, ilndvar, 'compress','xlat xlon')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add iy2d')

    ivdims(1) = idims(5)
    ivdims(2) = idims(4)
    istatus = nf90_def_var(ncid, 'PCT_PFT', nf90_double, ivdims(1:2), ipftvar)
    call checkncerr(istatus,__FILE__,__LINE__,  'Error add var pft')
    istatus = nf90_put_att(ncid, ipftvar, 'long_name','percent pft')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add pft long_name')
    istatus = nf90_put_att(ncid, ipftvar, 'units','%')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add pft units')
    istatus = nf90_put_att(ncid, ipftvar, '_FillValue',vmisdat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error add pft Fill Value')
    istatus = nf90_def_var(ncid, 'HARVEST_VH1', nf90_double, idims(5), iharvvh1)
    call checkncerr(istatus,__FILE__,__LINE__,  'Error add var HARVEST_VH1')
    istatus = nf90_put_att(ncid, iharvvh1, 'long_name', &
            'harvest from primary forest')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add HARVEST_VH1 long_name')
    istatus = nf90_def_var(ncid, 'HARVEST_VH2', nf90_double, idims(5), iharvvh2)
    call checkncerr(istatus,__FILE__,__LINE__,  'Error add var HARVEST_VH2')
    istatus = nf90_put_att(ncid, iharvvh2, 'long_name', &
            'harvest from primary non-forest')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add HARVEST_VH2 long_name')
    istatus = nf90_def_var(ncid, 'HARVEST_SH1', nf90_double, idims(5), iharvsh1)
    call checkncerr(istatus,__FILE__,__LINE__,  'Error add var HARVEST_SH1')
    istatus = nf90_put_att(ncid, iharvsh1, 'long_name', &
            'harvest from secondary mature-forest')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add HARVEST_SH1 long_name')
    istatus = nf90_def_var(ncid, 'HARVEST_SH2', nf90_double, idims(5), iharvsh2)
    call checkncerr(istatus,__FILE__,__LINE__,  'Error add var HARVEST_SH2')
    istatus = nf90_put_att(ncid, iharvsh2, 'long_name', &
            'harvest from secondary young-forest')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add HARVEST_SH2 long_name')
    istatus = nf90_def_var(ncid, 'HARVEST_SH3', nf90_double, idims(5), iharvsh3)
    call checkncerr(istatus,__FILE__,__LINE__,  'Error add var HARVEST_SH3')
    istatus = nf90_put_att(ncid, iharvsh3, 'long_name', &
            'harvest from secondary non-forest')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add HARVEST_SH3 long_name')
    istatus = nf90_def_var(ncid, 'GRAZING', nf90_double, idims(5), igraz)
    call checkncerr(istatus,__FILE__,__LINE__,  'Error add var GRAZING')
    istatus = nf90_put_att(ncid, igraz, 'long_name', &
            'grazing of herbacous pfts')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add GRAZING long_name')

    istatus = nf90_enddef(ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error exit define mode')

    hptop = real(ptop*10.0D0)
    call write_vertical_coord(ncid,rsigx,hptop,izvar)
    call write_horizontal_coord(ncid,xjx,yiy,ihvar)
    ipnt = 1
    call write_var2d_static(ncid,'xlat',rxlat,ipnt,illvar)
    call write_var2d_static(ncid,'xlon',rxlon,ipnt,illvar)

    istatus = nf90_put_var(ncid, ilndvar, landpoint)
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write gridcell')
    call mypack(xlat,gcvar)
    istatus = nf90_put_var(ncid, ilatvar, gcvar)
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write yclat')
    call mypack(xlon,gcvar)
    istatus = nf90_put_var(ncid, ilonvar, gcvar)
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write xclon')
    istatus = nf90_put_var(ncid, iiyvar, iiy)
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write iiy')
    istatus = nf90_put_var(ncid, ijxvar, ijx)
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write ijx')

    allocate(var3d(jxsg,iysg,npft))
    call mkdynpft(var3d(:,:,:),it)
    var3d(:,:,:) = dint(max(var3d(:,:,:),0.0D0))
    do np = 1 , npft
      where ( xmask < 0.5D0 )
        var3d(:,:,np) = vmisdat
      end where
    end do

    ! Here adjustment !
    do i = igstart , igstop
      do j = jgstart , jgstop
        if ( xmask(j,i) > 0.5D0 ) then
          if ( pctspec(j,i) > (100.0D0 - smallnum) ) then
            var3d(j,i,:) = 0.0D0
          else
            spft = 0.0D0
            do np = 1 , npft
              if ( var3d(j,i,np) > 0.0D0 ) then
                spft = spft + (var3d(j,i,np)*100.0D0)/(100.0D0-pctspec(j,i))
              end if
            end do
            if ( spft < smallnum ) then
              ! Substitute with something around it
              call bestaround_pft(var3d,i,j)
              spft = 0.0D0
              do np = 1 , npft
                if ( var3d(j,i,np) > 0.0D0 ) then
                  spft = spft + (var3d(j,i,np) * 100.0D0)/(100.0D0-pctspec(j,i))
                end if
              end do
              if ( spft < smallnum ) then
                call die(__FILE__,'No points around !',__LINE__)
              end if
            end if
            diff = spft - 100.0D0
            if ( abs(diff) > smallnum ) then
              ! Normalize it !
              do np = 1 , npft
                if ( var3d(j,i,np) > 0.0D0 ) then
                  operat = diff*(var3d(j,i,np)/spft)
                  var3d(j,i,np) = var3d(j,i,np) - operat
                  if ( var3d(j,i,np) < 0.0D0 ) then
                    var3d(j,i,np) = 0.0D0
                  end if
                end if
              end do
              ! Re-compute diff
              spft = 0.0D0
              do np = 1 , npft
                if ( var3d(j,i,np) > 0.0D0 ) then
                  spft = spft + (var3d(j,i,np) * 100.0D0)/(100.0D0-pctspec(j,i))
                end if
              end do
              diff = spft - 100.0D0
              if ( abs(diff) > smallnum ) then
                pft_gt0 = (var3d(j,i,:) > diff/dble(npft))
                do n = 1 , npft
                  if ( pft_gt0(n) .and. count(pft_gt0) > 0 ) then
                    var3d(j,i,n) = (var3d(j,i,n) - diff/dble(count(pft_gt0)))
                  end if
                end do
                spft = 0.0D0
                do np = 1 , npft
                  if ( var3d(j,i,np) > 0.0D0 ) then
                    spft = spft + &
                            (var3d(j,i,np) * 100.0D0)/(100.0D0-pctspec(j,i))
                  end if
                end do
                diff = spft - 100.0D0
                if ( abs(diff) > smallnum ) then
                  do np = 1 , npft
                    if ( var3d(j,i,np) > 0.0D0 ) then
                      operat = diff*(var3d(j,i,np)/spft)
                      var3d(j,i,np) = var3d(j,i,np) - operat
                      if ( var3d(j,i,np) < 0.0D0 ) then
                        var3d(j,i,np) = 0.0D0
                      end if
                    end if
                  end do
                  spft = 0.0D0
                  do np = 1 , npft
                    if ( var3d(j,i,np) > 0.0D0 ) then
                      spft = spft + &
                              (var3d(j,i,np) * 100.0D0)/(100.0D0-pctspec(j,i))
                    end if
                  end do
                  diff = spft - 100.0D0
                  if ( abs(diff) > smallnum ) then
                    write (0,*) abs(diff)
                    call die(__FILE__,'Cannot normalize...',__LINE__)
                  end if
                end if
              end if
            end if
          end if
        end if
      end do
    end do
    do ip = 1 , npft
      istart(1) = 1
      icount(1) = ngcells
      istart(2) = ip
      icount(2) = 1
      call mypack(var3d(:,:,ip),gcvar)
      istatus = nf90_put_var(ncid, ipftvar, gcvar, istart(1:2), icount(1:2))
      call checkncerr(istatus,__FILE__,__LINE__, 'Error write pft')
    end do
    deallocate(var3d)

    allocate(var3d(jxsg,iysg,6))
    call mkharvest(var3d,it)
    do i = 1 , 6
      where ( xmask < 0.5D0 )
        var3d(:,:,i) = vmisdat
      end where
    end do
    call mypack(var3d(:,:,1),gcvar)
    istatus = nf90_put_var(ncid, iharvvh1, gcvar)
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write HARVEST_VH1')
    call mypack(var3d(:,:,2),gcvar)
    istatus = nf90_put_var(ncid, iharvvh2, gcvar)
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write HARVEST_VH2')
    call mypack(var3d(:,:,3),gcvar)
    istatus = nf90_put_var(ncid, iharvsh1, gcvar)
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write HARVEST_SH1')
    call mypack(var3d(:,:,4),gcvar)
    istatus = nf90_put_var(ncid, iharvsh2, gcvar)
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write HARVEST_SH2')
    call mypack(var3d(:,:,5),gcvar)
    istatus = nf90_put_var(ncid, iharvsh3, gcvar)
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write HARVEST_SH3')
    call mypack(var3d(:,:,6),gcvar)
    istatus = nf90_put_var(ncid, igraz, gcvar)
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write GRAZING')

    istatus = nf90_close(ncid)
    call checkncerr(istatus,__FILE__,__LINE__,  &
      'Error close file '//trim(outfile))
  end do
#endif
#endif

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

  subroutine bestaround_pft(pft,i,j)
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
            if ( pft(js,is,n) > smallnum ) then
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
  end subroutine bestaround_pft

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

#endif
end program mksurfdata
