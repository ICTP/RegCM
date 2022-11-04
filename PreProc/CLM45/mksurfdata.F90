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
  stop ' Execution terminated because of runtime error'
end subroutine myabort

program mksurfdata

#ifndef CLM45
  use mod_stdio
  write(stdout,*) 'RegCM built without CLM45 support.'
  write(stdout,*) 'Please recompile it using --enable-clm45 flag'
#else

#ifdef CROP
#ifndef CN
  ERROR : CN MUST BE DEFINED ON ACTIVATING CROP
#endif
#endif

#ifdef LCH4
#ifndef CN
  ERROR : CN MUST BE DEFINED ON ACTIVATING LCH4
#endif
#endif

  use mod_intkinds
  use mod_constants , only : raddeg , dlowval
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
  use mod_zita
  use netcdf
#ifdef CN
  use mod_mklightning
  use mod_mkpopd
  use mod_mkq10soil
  use mod_mkndep
#endif
#ifdef DYNPFT
#ifdef CN
  use mod_mkharvest
#endif
  use mod_mkdynpft
#endif
#ifdef CN
#ifdef LCH4
  use mod_mkch4topm
  use mod_mksoilph
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
#endif
#ifdef DYNPFT
  integer(ik4) , parameter :: noleap_yday_3h = 365*8
  integer(ik4) , parameter :: nyears = 2100-1850+1
  integer(ik4) :: y1 , y2 , mon , day , hour
  character(len=4) :: cy
  character(len=32) :: p1 , p2
#endif
#ifdef LUCASPFT
  integer(ik4) :: y1 , y2 , mon , day , hour
  character(len=32) :: p1 , p2
#endif

  integer(ik4) :: ngcells

  real(rkx) , parameter :: vmisdat = -9999.0_rkx

  integer(ik4) :: istatus , ncid , ndim , nvar
  integer(ik4) , dimension(12) :: idims , ivdims
  integer(ik4) :: ivartime , iglcvar , iwetvar , ilakevar , iurbanvar
  integer(ik4) :: islope2d , istddev2d
  integer(ik4) :: ipft2dvar , iglc2dvar , iwet2dvar , ilake2dvar , iurban2dvar
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
  integer(ik4) :: q10 , ndep
#ifdef LCH4
  integer(ik4) :: k , q , v , maxf
  integer(ik4) :: isoilphvar
#endif
#endif
#ifdef DYNPFT
#ifdef CN
  integer(ik4) :: iharvvh1 , iharvvh2 , iharvsh1 , iharvsh2 , iharvsh3 , igraz
#endif
  character(len=256) :: dynfile
#endif
#ifdef VICHYDRO
  integer(ik4) :: ibinfl , ids , idsmax ,iws
#endif
  integer(ik4) :: ijxvar , iiyvar
  type(rcm_time_and_date) :: irefdate , imondate
  type(rcm_time_interval) :: tdif
  real(rk4) , pointer , dimension(:) :: yiy => null( )
  real(rk4) , pointer , dimension(:) :: xjx => null( )
  real(rk4) :: hptop
  real(rk8) , dimension(1) :: xdate
  integer(ik4) , dimension(3) :: istart , icount
  integer(ik4) , dimension(2) :: ihvar
  integer(ik4) , dimension(2) :: illvar
  integer(ik4) , dimension(3) :: izvar
  integer(ik4) , dimension(1) :: istart1 , icount1 , mxsoil_color , iloc
  real(rkx) :: spft , mean , diff , argf
  integer(ik4) :: ierr
  integer(ik4) :: i , j , ip , il , ir , iu , it , ipnt , iurbmax
  integer(ik4) :: jgstart , jgstop , igstart , igstop
  character(len=256) :: namelistfile , prgname
  character(len=256) :: terfile , outfile
  character(len=64) :: csdate , pftfile , laifile
  real(rkx) , dimension(:,:) , pointer :: pctspec => null( )
  real(rkx) , dimension(:,:) , pointer :: pctslake => null( )
  real(rkx) , dimension(:,:) , pointer :: pctbare => null( )
  real(rkx) , pointer , dimension(:,:) :: var2d => null( )
  integer(ik4) , pointer , dimension(:,:) :: ivar2d => null( )
  real(rkx) , pointer , dimension(:,:,:) :: var3d => null( )
  real(rkx) , pointer , dimension(:,:,:) :: var3dp => null( )
  real(rkx) , pointer , dimension(:,:,:,:) :: var4d => null( )
  real(rkx) , pointer , dimension(:,:,:,:,:) :: var5d => null( )
  real(rkx) , pointer , dimension(:,:,:,:,:,:) :: var6d => null( )
  real(rkx) , pointer , dimension(:) :: gcvar => null( )
  integer(ik4) , pointer , dimension(:) :: igcvar => null( )
  integer(ik4) , pointer , dimension(:) :: iiy => null( )
  integer(ik4) , pointer , dimension(:) :: ijx => null( )
  integer(ik4) , pointer , dimension(:) :: landpoint => null( )
  logical , pointer , dimension(:) :: pft_gt0 => null( )
  logical :: subgrid
  integer(ik4) :: hostnm
  integer(ik4) :: ihost , idir
  integer(ik4) :: getcwd
  integer(ik4) , dimension(8) :: tval
  character (len=32) :: cdata='?'
  character (len=5) :: czone='?'
  character (len=32) :: hostname='?'
  character (len=32) :: user='?'
  character (len=128) :: directory='?'
  character (len=*) , parameter :: f99001 = &
          '(2x," GIT Revision: ",a," compiled at: data : ",a,"  time: ",a,/)'

  write (stdout,  &
     "(/,2x,'This is mksurfdata part of RegCM package version 4')")
  write (stdout,f99001)  GIT_VER, __DATE__ , __TIME__

#ifdef IBM
  hostname='ibm platform '
  user= 'Unknown'
#else
  ihost = hostnm(hostname)
  call getlog(user)
#endif
  call date_and_time(zone=czone,values=tval)
  idir = getcwd(directory)

  write(cdata,'(i0.4,"-",i0.2,"-",i0.2," ",i0.2,":",i0.2,":",i0.2,a)') &
     tval(1), tval(2), tval(3), tval(5), tval(6), tval(7), czone
  write(stdout,*) ": this run start at  : ",trim(cdata)
  write(stdout,*) ": it is submitted by : ",trim(user)
  write(stdout,*) ": it is running on   : ",trim(hostname)
  write(stdout,*) ": in directory       : ",trim(directory)
  write(stdout,*) "                     "

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
        case ('85', '15')
          p2 = 'SCENARIO'//pthsep//'RCP8.5'
        case default
          if ( dattyp /= "EIN15" .and. &
               dattyp(1:4) /= "NNRP" .and. &
               dattyp /= "JRA55" ) then
            call die(__FILE__, &
              'Dynamic landuse only supported for CMIP5',__LINE__)
          end if
      end select
    end if
    pftfile = trim(p1)//pthsep//trim(p2)//pthsep//'mksrf_landuse_'//cy//'.nc'
#else
#ifdef LUCAS_PFT
    call split_idate(globidate1,y1,mon,day,hour)
    if ( y1 < 2015 ) then
      pftfile = 'alternative'//pthsep// &
        'LUCAS_LUC_v09_ESACCI_LUH2_historical_1950_2015.nc'
    else
      pftfile = 'alternative'//pthsep// &
        'LUCAS_LUC_v09_ESACCI_LUH2_rcp26_2015_2100.nc'
    end if
#else
    pftfile = 'mksrf_pft.nc'
#endif
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
  ngcells = count(xmask(jgstart:jgstop,igstart:igstop) > 0.5_rkx)
  call closefile(ncid)
  !
  ! Open Output in NetCDF format
  !
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

  istatus = nf90_def_var(ncid, 'time', regcm_vartype, idims(4:4), ivartime)
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

  istatus = nf90_def_var(ncid, 'xclon', regcm_vartype, idims(7), ilonvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var xclon')
  istatus = nf90_put_att(ncid, ilndvar, 'compress','xlat xlon')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add xclon')

  istatus = nf90_def_var(ncid, 'yclat', regcm_vartype, idims(7), ilatvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var yclat')
  istatus = nf90_put_att(ncid, ilndvar, 'compress','xlat xlon')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add yclat')

  istatus = nf90_def_var(ncid, 'topo', regcm_vartype, idims(7), itopovar)
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

  istatus = nf90_def_var(ncid, 'SLOPE', regcm_vartype, idims(7), islopevar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var slope')
  istatus = nf90_put_att(ncid, islopevar, 'long_name','Elevation slope')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add slope long_name')
  istatus = nf90_put_att(ncid, islopevar, 'units','degree')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add slope units')
  istatus = nf90_def_var(ncid, 'STD_ELEV', regcm_vartype, idims(7), istdvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var stddev')
  istatus = nf90_put_att(ncid, istdvar, 'long_name','Elevation std dev')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add stddev long_name')
  istatus = nf90_put_att(ncid, istdvar, 'units','m')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add stddev units')

  istatus = nf90_def_var(ncid, 'PCT_GLACIER', regcm_vartype, idims(7), iglcvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var glacier')
  istatus = nf90_put_att(ncid, iglcvar, 'long_name','percent glacier')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add glacier long_name')
  istatus = nf90_put_att(ncid, iglcvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add glacier units')

  istatus = nf90_def_var(ncid, 'PCT_WETLAND', regcm_vartype, idims(7), iwetvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var wetland')
  istatus = nf90_put_att(ncid, iwetvar, 'long_name','percent wetland')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add wetland long_name')
  istatus = nf90_put_att(ncid, iwetvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add wetland units')

  istatus = nf90_def_var(ncid, 'PCT_LAKE', regcm_vartype, idims(7), ilakevar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var lake')
  istatus = nf90_put_att(ncid, ilakevar, 'long_name','percent lake')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add lake long_name')
  istatus = nf90_put_att(ncid, ilakevar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add lake units')

  ivdims(1) = idims(7)
  ivdims(2) = idims(8)
  istatus = nf90_def_var(ncid, 'PCT_URBAN', regcm_vartype, &
                         ivdims(1:2), iurbanvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var urban')
  istatus = nf90_put_att(ncid, iurbanvar, 'long_name','percent urban')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add urban long_name')
  istatus = nf90_put_att(ncid, iurbanvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add urban units')

  ivdims(1) = idims(7)
  ivdims(2) = idims(5)
  istatus = nf90_def_var(ncid, 'PCT_PFT', regcm_vartype, ivdims(1:2), ipftvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var pft')
  istatus = nf90_put_att(ncid, ipftvar, 'long_name','percent pft')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add pft long_name')
  istatus = nf90_put_att(ncid, ipftvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add pft units')
  istatus = nf90_put_att(ncid, ipftvar, '_FillValue',vmisdat)
  call checkncerr(istatus,__FILE__,__LINE__,'Error add pft Fill Value')

  ivdims(1) = idims(1)
  ivdims(2) = idims(2)
  ivdims(3) = idims(5)
  istatus = nf90_def_var(ncid, 'pft_2d', regcm_vartype, ivdims(1:3), ipft2dvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var pft')
  istatus = nf90_put_att(ncid, ipft2dvar, 'long_name','percent pft')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add pft long_name')
  istatus = nf90_put_att(ncid, ipft2dvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add pft units')
  istatus = nf90_put_att(ncid, ipft2dvar, '_FillValue',vmisdat)
  call checkncerr(istatus,__FILE__,__LINE__,'Error add pft Fill Value')

  istatus = nf90_def_var(ncid, 'slope', regcm_vartype, ivdims(1:2), islope2d)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var mxc_2d')
  istatus = nf90_put_att(ncid, islope2d, '_FillValue',vmisdat)
  call checkncerr(istatus,__FILE__,__LINE__,'Error add mxc_2d Fill Value')

  istatus = nf90_def_var(ncid, 'stddv', regcm_vartype, ivdims(1:2), istddev2d)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var mxc_2d')
  istatus = nf90_put_att(ncid, istddev2d, '_FillValue',vmisdat)
  call checkncerr(istatus,__FILE__,__LINE__,'Error add mxc_2d Fill Value')

  istatus = nf90_def_var(ncid, 'glc_2d', regcm_vartype, ivdims(1:2), iglc2dvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var glacier')
  istatus = nf90_put_att(ncid, iglc2dvar, 'long_name','percent glacier')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add glacier long_name')
  istatus = nf90_put_att(ncid, iglc2dvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add glacier units')
  istatus = nf90_put_att(ncid, iglc2dvar, '_FillValue',vmisdat)
  call checkncerr(istatus,__FILE__,__LINE__,'Error add glacier Fill Value')

  istatus = nf90_def_var(ncid, 'wet_2d', regcm_vartype, ivdims(1:2), iwet2dvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var wetland')
  istatus = nf90_put_att(ncid, iwet2dvar, 'long_name','percent wetland')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add wetland long_name')
  istatus = nf90_put_att(ncid, iwet2dvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add wetland units')
  istatus = nf90_put_att(ncid, iwet2dvar, '_FillValue',vmisdat)
  call checkncerr(istatus,__FILE__,__LINE__,'Error add wetland Fill Value')

  istatus = nf90_def_var(ncid, 'lak_2d', regcm_vartype, ivdims(1:2), ilake2dvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var lake')
  istatus = nf90_put_att(ncid, ilake2dvar, 'long_name','percent lake')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add lake long_name')
  istatus = nf90_put_att(ncid, ilake2dvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add lake units')
  istatus = nf90_put_att(ncid, ilake2dvar, '_FillValue',vmisdat)
  call checkncerr(istatus,__FILE__,__LINE__,'Error add lake Fill Value')

  ivdims(1) = idims(1)
  ivdims(2) = idims(2)
  ivdims(3) = idims(8)
  istatus = nf90_def_var(ncid, 'urb_2d', regcm_vartype, &
                         ivdims(1:3), iurban2dvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var urban')
  istatus = nf90_put_att(ncid, iurban2dvar, 'long_name','percent urban')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add urban long_name')
  istatus = nf90_put_att(ncid, iurban2dvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add urban units')
  istatus = nf90_put_att(ncid, iurban2dvar, '_FillValue',vmisdat)
  call checkncerr(istatus,__FILE__,__LINE__,'Error add urban Fill Value')

  ivdims(1) = idims(7)
  ivdims(2) = idims(5)
  ivdims(3) = idims(4)
  istatus = nf90_def_var(ncid, 'MONTHLY_LAI', regcm_vartype, &
                         ivdims(1:3), ilaivar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var monthly_lai')
  istatus = nf90_put_att(ncid, ilaivar, 'long_name','monthly leaf area index')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_lai long_name')
  istatus = nf90_put_att(ncid, ilaivar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_lai units')
  istatus = nf90_def_var(ncid, 'MONTHLY_SAI', regcm_vartype, &
                         ivdims(1:3), isaivar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var monthly_sai')
  istatus = nf90_put_att(ncid, isaivar, 'long_name','monthly stem area index')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_sai long_name')
  istatus = nf90_put_att(ncid, isaivar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_sai units')
  istatus = nf90_def_var(ncid, 'MONTHLY_HEIGHT_TOP', regcm_vartype, &
          ivdims(1:3), ivgtopvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var monthly_top')
  istatus = nf90_put_att(ncid, ivgtopvar, 'long_name','monthly height top')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_top long_name')
  istatus = nf90_put_att(ncid, ivgtopvar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_top units')
  istatus = nf90_def_var(ncid, 'MONTHLY_HEIGHT_BOT', regcm_vartype, &
          ivdims(1:3), ivgbotvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var monthly_bot')
  istatus = nf90_put_att(ncid, ivgbotvar, 'long_name','monthly height bottom')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_bot long_name')
  istatus = nf90_put_att(ncid, ivgbotvar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add monthly_bot units')

  istatus = nf90_def_var(ncid, 'FMAX', regcm_vartype, idims(7), ifmaxvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var fmax')
  istatus = nf90_put_att(ncid, ifmaxvar, 'long_name', &
          'maximum fractional saturated area')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add fmax long_name')
  istatus = nf90_put_att(ncid, ifmaxvar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add fmax units')

  istatus = nf90_def_var(ncid, 'SOIL_COLOR', nf90_int, &
                         idims(7),isoilcolvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var soilcol')
  istatus = nf90_put_att(ncid, isoilcolvar, 'long_name', 'Soil color')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add soilcol long_name')
  istatus = nf90_put_att(ncid, isoilcolvar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add soilcol units')

  istatus = nf90_def_var(ncid, 'gdp', regcm_vartype, idims(7),igdpvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var gdp')
  istatus = nf90_put_att(ncid, igdpvar, 'long_name', 'real GDP')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add gdp long_name')
  istatus = nf90_put_att(ncid, igdpvar, 'units','K 1995US$/capita')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add gdp units')

  istatus = nf90_def_var(ncid, 'peatf', regcm_vartype, idims(7),ipeatfvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var peatf')
  istatus = nf90_put_att(ncid, ipeatfvar, 'long_name', &
          'global fraction of peatland')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add peatf long_name')
  istatus = nf90_put_att(ncid, ipeatfvar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add peatf units')

  istatus = nf90_def_var(ncid, 'abm', nf90_int, idims(7),iabmvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var abm')
  istatus = nf90_put_att(ncid, iabmvar, 'long_name', &
          'peak month for agri fire')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add abm long_name')
  istatus = nf90_put_att(ncid, iabmvar, 'units','1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add abm units')

  istatus = nf90_def_var(ncid, 'EF1_BTR', regcm_vartype, idims(7),ief_btrvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ef_btr')
  istatus = nf90_put_att(ncid, ief_btrvar, 'long_name', &
          'broadleaf tree emission factor')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_btr long_name')
  istatus = nf90_put_att(ncid,ief_btrvar,'units','micrograms isoprene m-2 h-1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_btr units')
  istatus = nf90_def_var(ncid, 'EF1_CRP', regcm_vartype, idims(7),ief_crpvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ef_crp')
  istatus = nf90_put_att(ncid, ief_crpvar, 'long_name', &
          'crop emission factor')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_crp long_name')
  istatus = nf90_put_att(ncid,ief_crpvar,'units','micrograms isoprene m-2 h-1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_crp units')
  istatus = nf90_def_var(ncid, 'EF1_FDT', regcm_vartype, idims(7),ief_fdtvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ef_fdt')
  istatus = nf90_put_att(ncid, ief_fdtvar, 'long_name', &
          'Fineleaf deciduous tree emission factor')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_fdt long_name')
  istatus = nf90_put_att(ncid,ief_fdtvar,'units','micrograms isoprene m-2 h-1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_fdt units')
  istatus = nf90_def_var(ncid, 'EF1_FET', regcm_vartype, idims(7),ief_fetvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ef_fet')
  istatus = nf90_put_att(ncid, ief_fetvar, 'long_name', &
          'Fineleaf evergreen tree emission factor')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_fet long_name')
  istatus = nf90_put_att(ncid,ief_fetvar,'units','micrograms isoprene m-2 h-1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_fet units')
  istatus = nf90_def_var(ncid, 'EF1_GRS', regcm_vartype, idims(7),ief_grsvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ef_grs')
  istatus = nf90_put_att(ncid, ief_grsvar, 'long_name', &
          'grass, non-vascular plants and other ground cover emission factor')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_grs long_name')
  istatus = nf90_put_att(ncid,ief_grsvar,'units','micrograms isoprene m-2 h-1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_grs units')
  istatus = nf90_def_var(ncid, 'EF1_SHR', regcm_vartype, idims(7),ief_shrvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ef_shr')
  istatus = nf90_put_att(ncid, ief_shrvar, 'long_name', &
          'shrub emission factor')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_shr long_name')
  istatus = nf90_put_att(ncid,ief_shrvar,'units','micrograms isoprene m-2 h-1')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ef_shr units')

  ivdims(1) = idims(7)
  ivdims(2) = idims(6)
  istatus = nf90_def_var(ncid, 'PCT_SAND', regcm_vartype, ivdims(1:2), isandvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var sand')
  istatus = nf90_put_att(ncid, isandvar, 'long_name','percent sand')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add sand long_name')
  istatus = nf90_put_att(ncid, isandvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add sand units')
  istatus = nf90_def_var(ncid, 'PCT_CLAY', regcm_vartype, ivdims(1:2), iclayvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var clay')
  istatus = nf90_put_att(ncid, iclayvar, 'long_name','percent clay')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add clay long_name')
  istatus = nf90_put_att(ncid, iclayvar, 'units','%')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add clay units')

  istatus = nf90_def_var(ncid, 'ORGANIC', regcm_vartype, &
                         ivdims(1:2), iorganicvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var organic')
  istatus = nf90_put_att(ncid, iorganicvar, 'long_name', &
          'organic soil density at soil levels')
  istatus = nf90_put_att(ncid, iorganicvar, 'comment', &
          'assumed carbon content 0.58 gC per gOM')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add organic long_name')
  istatus = nf90_put_att(ncid, iorganicvar, 'units','kg m-3')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add organic units')

  istatus = nf90_def_var(ncid, 'LAKEDEPTH',regcm_vartype,idims(7),idepthvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var lake')
  istatus = nf90_put_att(ncid, idepthvar, 'long_name', &
          'Lake Depth (default = 10m)')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add lake long_name')
  istatus = nf90_put_att(ncid, idepthvar, 'units','m')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add lake units')

  ivdims(1) = idims(7)
  ivdims(2) = idims(8)
  do i = 1 , npu2d
    if ( parm2d(i) == 'NLEV_IMPROAD' ) then
      istatus = nf90_def_var(ncid, parm2d(i), nf90_int, &
                             ivdims(1:2), iurb2d(i))
    else
      istatus = nf90_def_var(ncid, parm2d(i), regcm_vartype, &
                             ivdims(1:2), iurb2d(i))
    end if
    call checkncerr(istatus,__FILE__,__LINE__, 'Error add var')
    istatus = nf90_put_att(ncid, iurb2d(i), 'long_name', lngn2d(i))
    call checkncerr(istatus,__FILE__,__LINE__,'Error add long_name')
    istatus = nf90_put_att(ncid, iurb2d(i), 'units', unit2d(i) )
    call checkncerr(istatus,__FILE__,__LINE__,'Error add units')
    if ( parm2d(i) /= 'NLEV_IMPROAD' ) then
      istatus = nf90_put_att(ncid, iurb2d(i), '_FillValue',vmisdat)
      call checkncerr(istatus,__FILE__,__LINE__,'Error add _FillValue')
    end if
  end do

  ivdims(1) = idims(7)
  ivdims(2) = idims(8)
  ivdims(3) = idims(9)
  do i = 1 , npu3d
    istatus = nf90_def_var(ncid, parm3d(i), regcm_vartype, &
                           ivdims(1:3), iurb3d(i))
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
      regcm_vartype, ivdims(1:3), iurb4d(2*i-1))
    call checkncerr(istatus,__FILE__,__LINE__, 'Error add var')
    istatus = nf90_put_att(ncid, iurb4d(2*i-1), 'long_name', lngn4d(i))
    call checkncerr(istatus,__FILE__,__LINE__,'Error add long_name')
    istatus = nf90_put_att(ncid, iurb4d(2*i-1), 'units', unit4d(i) )
    call checkncerr(istatus,__FILE__,__LINE__,'Error add units')
    istatus = nf90_put_att(ncid, iurb4d(2*i-1), '_FillValue',vmisdat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error add _FillValue')
    istatus = nf90_def_var(ncid, trim(parm4d(i))//'_DIF', &
      regcm_vartype, ivdims(1:3), iurb4d(2*i))
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
  istatus = nf90_def_var(ncid, 'LNFM',regcm_vartype,ivdims(1:2),ilightning)
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
  istatus = nf90_def_var(ncid, 'HDM',regcm_vartype,ivdims(1:2),ipopden)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var hdm')
  istatus = nf90_put_att(ncid, ipopden, 'long_name', &
          'human population density')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add hdm long_name')
  istatus = nf90_put_att(ncid, ipopden, 'units','counts/km^2')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add hdm units')

#if defined(CN)
  ! samy : reading global gridded Q10 soil respiration parameter
  istatus = nf90_def_var(ncid, 'Q10',regcm_vartype, idims(7),q10)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var Q10')
  istatus = nf90_put_att(ncid, q10, 'long_name', &
          'Global Q10 soil respiration')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add Q10 long_name')
  istatus = nf90_put_att(ncid, q10, 'units','unitless')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add Q10 units')
  istatus = nf90_def_var(ncid, 'NDEP',regcm_vartype, idims(7),ndep)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ndep')
  istatus = nf90_put_att(ncid, ndep, 'long_name', &
            'Global Nitrogen deposition')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add NDEP long_name')
  istatus = nf90_put_att(ncid, ndep, 'units','g(N)/m2/yr')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add NDEP units')
#endif

#ifdef LCH4
  istatus = nf90_def_var(ncid, 'K_PAR',regcm_vartype,idims(7),k)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var K_PAR')
  istatus = nf90_def_var(ncid, 'XM_PAR',regcm_vartype,idims(7),q)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var XM_PAR')
  istatus = nf90_def_var(ncid, 'V_PAR',regcm_vartype,idims(7),v)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var V_PAR')
  istatus = nf90_def_var(ncid, 'MAXF',regcm_vartype,idims(7),maxf)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var MAXF')
  ! samy : for soil ph for CH4 emission
  istatus = nf90_def_var(ncid, 'PH', regcm_vartype, idims(7),isoilphvar)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var PH')
  istatus = nf90_put_att(ncid, isoilphvar, 'long_name', &
          'Global soil pH')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add PH long_name')
  istatus = nf90_put_att(ncid, isoilphvar, 'units','unitless')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add PH units')
#endif
#endif

#ifdef VICHYDRO
  istatus = nf90_def_var(ncid, 'binfl',regcm_vartype,idims(7),ibinfl)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var binfl')
  istatus = nf90_put_att(ncid, ibinfl, 'long_name', &
          'VIC b parameter for the Variable Infiltration Capacity Curve')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add binfl long_name')
  istatus = nf90_put_att(ncid, ibinfl, 'units','unitless')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add binfl units')
  istatus = nf90_def_var(ncid, 'Ds',regcm_vartype,idims(7),ids)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ds')
  istatus = nf90_put_att(ncid, ids, 'long_name', &
          'VIC Ds parameter for the ARNO curve')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ds long_name')
  istatus = nf90_put_att(ncid, ids, 'units','unitless')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ds units')
  istatus = nf90_def_var(ncid, 'Dsmax',regcm_vartype,idims(7),idsmax)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var dsmax')
  istatus = nf90_put_att(ncid, idsmax, 'long_name', &
          'VIC Dsmax parameter for the ARNO curve')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add dsmax long_name')
  istatus = nf90_put_att(ncid, idsmax, 'units','unitless')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add dsmax units')
  istatus = nf90_def_var(ncid, 'Ws',regcm_vartype,idims(7),iws)
  call checkncerr(istatus,__FILE__,__LINE__,  'Error add var ws')
  istatus = nf90_put_att(ncid, iws, 'long_name', &
          'VIC Ws parameter for the ARNO curve')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ws long_name')
  istatus = nf90_put_att(ncid, iws, 'units','unitless')
  call checkncerr(istatus,__FILE__,__LINE__,'Error add ws units')
#endif

  istatus = nf90_enddef(ncid)
  call checkncerr(istatus,__FILE__,__LINE__,'Error exit define mode')

  if ( idynamic < 3 ) then
    hptop = real(ptop*10.0_rkx)
    call write_vertical_coord(ncid,rsigx,hptop,izvar)
  else
    call model_zitah(zita)
    ax = real(md_ak(zita),rk4)
    bx = real(md_bk(zita),rk4)
    call write_vertical_coord_zita(ncid,rsigx,ax,bx,izvar)
  end if
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
  call getmem2d(pctbare,1,jxsg,1,iysg,'mksurfdata: pctbare')
  call getmem1d(gcvar,1,ngcells,'mksurfdata: gcvar')
  call getmem1d(igcvar,1,ngcells,'mksurfdata: igcvar')
  call getmem1d(iiy,1,ngcells,'mksurfdata: iiy')
  call getmem1d(ijx,1,ngcells,'mksurfdata: ijx')
  call getmem1d(landpoint,1,ngcells,'mksurfdata: landpoint')
  pctspec(:,:) = 0.0_rkx
  pctslake(:,:) = 0.0_rkx
  pctbare(:,:) = 0.0_rkx
  ip = 1
  do i = igstart , igstop
    do j = jgstart , jgstop
      if ( xmask(j,i) > 0.5_rkx ) then
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
      if ( xmask(j,i) > 0.5_rkx ) then
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
  var3d(:,:,1) = 0.0_rkx
  ! Calculate slope and std
  do i = igstart , igstop
    do j = jgstart , jgstop
      argf = sum(topo(j-1:j+1,i-1:i+1)-topo(j,i))/8.0_rkx
      if ( abs(argf) > dlowval ) then
        var3d(j,i,1) = raddeg * &
          atan((sum(topo(j-1:j+1,i-1:i+1)-topo(j,i))/8.0_rkx)/(ds*1000.0_rkx))
      else
        var3d(j,i,1) = 0.0_rkx
      end if
      mean = sum(topo(j-1:j+1,i-1:i+1))/9.0_rkx
      var3d(j,i,2) = sqrt(sum((topo(j-1:j+1,i-1:i+1)-mean)**2)/8.0_rkx)
    end do
  end do
  where ( xmask < 0.5_rkx )
    var3d(:,:,1) = vmisdat
    var3d(:,:,2) = vmisdat
  end where
  call mypack(var3d(:,:,1),gcvar)
  gcvar = gcvar*real(raddeg,rkx)
  istatus = nf90_put_var(ncid, islopevar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write slope')
  call mypack(var3d(:,:,2),gcvar)
  istatus = nf90_put_var(ncid, istdvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write stddev')
  istatus = nf90_put_var(ncid, islope2d, var3d(:,:,1))
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write slope')
  istatus = nf90_put_var(ncid, istddev2d, var3d(:,:,2))
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write sttdev')
  deallocate(var3d)

  write(stdout,*) 'Created topographic informations...'

  iurbmax = numurbl+3
  allocate(var3d(jxsg,iysg,iurbmax))
  call mkglacier('mksrf_glacier.nc',xmask,var3d(:,:,1))
  call mkwetland('mksrf_lanwat.nc',xmask,var3d(:,:,2),var3d(:,:,3))
  call mkurban_base('mksrf_urban.nc',xmask,var3d(:,:,4:iurbmax))
  var3d = nint(var3d)
  if ( .not. enable_urban_landunit ) then
    write (stderr,*) 'Disable URBAN Areas in CLM4.5 Model !'
    var3d(:,:,4:iurbmax) = 0.0_rkx
  end if

  !
  ! Normalize the sum of pctspec to range 0-100
  !
  do i = igstart , igstop
    do j = jgstart , jgstop
      if ( xmask(j,i) > 0.5_rkx ) then
        pctspec(j,i) = sum(var3d(j,i,:))
        ! If 2 or more special classes at 100 % on same point
        if ( pctspec(j,i) >= 200.0_rkx ) then
          var3d(j,i,:) = var3d(j,i,:) / (pctspec(j,i)/100.0_rkx)
          pctspec(j,i) = sum(var3d(j,i,:))
        end if
        ! Reduce to biggest if the sum is not less equal 100
        if ( pctspec(j,i) > 100.0_rkx ) then
          iloc = maxloc(var3d(j,i,:))
          var3d(j,i,:) = 0.0_rkx
          var3d(j,i,iloc(1)) = 100.0_rkx
          pctspec(j,i) = 100.0_rkx
        end if
      end if
    end do
  end do

  if ( any(pctspec < 0.0_rkx) .or. any(pctspec > 100.0_rkx) ) then
    write(stderr,*) 'Error in special categories.'
    call die(__FILE__,'Cannot continue!',__LINE__)
  end if

  allocate(var3dp(jxsg,iysg,npft))
  call mkpft(pftfile,xmask,var3dp(:,:,:))
  var3dp = nint(var3dp)

  ! Here adjustment !
  do i = igstart , igstop
    do j = jgstart , jgstop
      if ( xmask(j,i) > 0.5_rkx ) then
        if ( pctspec(j,i) > 99.9_rkx ) then
          var3dp(j,i,:) = 0.0_rkx
        else
          spft = sum(var3dp(j,i,:))
          if ( spft < -1.0_rkx ) then
            call die(__FILE__,'No points around !',__LINE__)
          end if
          diff = spft - (100.0_rkx-pctspec(j,i))
          if ( abs(diff) > 0.5_rkx ) then
            do while ( diff >= 1.0_rkx )
              iloc = maxloc(var3dp(j,i,:))
              var3dp(j,i,iloc(1)) = nint(var3dp(j,i,iloc(1)) - 1._rkx)
              spft = sum(var3dp(j,i,:))
              diff = spft - (100.0_rkx-pctspec(j,i))
            end do
            do while ( diff <= -1.0_rkx )
              iloc = maxloc(var3dp(j,i,:))
              var3dp(j,i,iloc(1)) = nint(var3dp(j,i,iloc(1)) + 1._rkx)
              spft = sum(var3dp(j,i,:))
              diff = spft - (100.0_rkx-pctspec(j,i))
            end do
          end if
        end if
        if ( var3dp(j,i,1) > 50.0_rkx ) then
          pctbare(j,i) = var3dp(j,i,1)
        end if
      end if
    end do
  end do
  !var3dp(:,:,15) = var3dp(:,:,15) + var3dp(:,:,5)
  !var3dp(:,:,5) = 0.0_rkx
  !var3dp(:,:,15) = var3dp(:,:,15) + var3dp(:,:,7)
  !var3dp(:,:,7) = 0.0_rkx

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
  where ( var3d(:,:,3) > 0.0_rkx )
    pctslake = var3d(:,:,3)
  end where
  istatus = nf90_put_var(ncid, iglc2dvar, var3d(:,:,1))
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write glc2d')
  istatus = nf90_put_var(ncid, iwet2dvar, var3d(:,:,2))
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write wet2d')
  istatus = nf90_put_var(ncid, ilake2dvar, var3d(:,:,3))
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write lak2d')
  istatus = nf90_put_var(ncid, iurban2dvar, var3d(:,:,4:))
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write urb2d')
  deallocate(var3d)

  write(stdout,*) 'Created special categories informations...'

  do ip = 1 , npft
    istart(1) = 1
    icount(1) = ngcells
    istart(2) = ip
    icount(2) = 1
    call mypack(var3dp(:,:,ip),gcvar)
    istatus = nf90_put_var(ncid, ipftvar, gcvar, istart(1:2), icount(1:2))
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write pft')
  end do
  istatus = nf90_put_var(ncid, ipft2dvar, var3dp)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write pft2d')
  deallocate(var3dp)

  write(stdout,*) 'Created pft informations...'

  allocate(var5d(jxsg,iysg,npft,nmon,4))
  call mklaisai(laifile,xmask,var5d(:,:,:,:,1), var5d(:,:,:,:,2), &
                var5d(:,:,:,:,3), var5d(:,:,:,:,4))

  do it = 1 , nmon
    do ip = 1 , npft
      istart(1) = 1
      icount(1) = ngcells
      istart(2) = ip
      icount(2) = 1
      istart(3) = it
      icount(3) = 1
      call mypack(var5d(:,:,ip,it,1),gcvar)
      istatus = nf90_put_var(ncid, ilaivar, gcvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, 'Error write monthly_lai')
      call mypack(var5d(:,:,ip,it,2),gcvar)
      istatus = nf90_put_var(ncid, isaivar, gcvar,istart,icount)
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

  write(stdout,*) 'Created lai/sai informations...'

  allocate(var2d(jxsg,iysg))
  call mkfmax('mksrf_fmax.nc',xmask,var2d)
  call mypack(var2d,gcvar)
  istatus = nf90_put_var(ncid, ifmaxvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write fmax')
  deallocate(var2d)

  allocate(ivar2d(jxsg,iysg))
  call mksoilcol('mksrf_soicol.nc',xmask,ivar2d)
  call mypack(ivar2d,igcvar)
  istatus = nf90_put_var(ncid, isoilcolvar, igcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write soil color')
  deallocate(ivar2d)

  allocate(var4d(jxsg,iysg,nsoil,2))
  call mksoitex('mksrf_soitex.nc',xmask,var4d(:,:,:,1),var4d(:,:,:,2))
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

  write(stdout,*) 'Created soil texture/color informations...'

  allocate(var2d(jxsg,iysg))
  call mkgdp('mksrf_gdp.nc',xmask,var2d)
  call mypack(var2d,gcvar)
  istatus = nf90_put_var(ncid, igdpvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write gdp')
  deallocate(var2d)

  allocate(var2d(jxsg,iysg))
  call mkpeatf('mksrf_peatf.nc',xmask,var2d)
  call mypack(var2d,gcvar)
  istatus = nf90_put_var(ncid, ipeatfvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write peatf')
  deallocate(var2d)

  allocate(ivar2d(jxsg,iysg))
  call mkabm('mksrf_abm.nc',xmask,ivar2d)
  call mypack(ivar2d,igcvar)
  istatus = nf90_put_var(ncid, iabmvar, igcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write abm')
  deallocate(ivar2d)

  write(stdout,*) 'Created human population informations...'

  allocate(var3d(jxsg,iysg,6))
  call mkvocef('mksrf_vocef.nc',xmask,var3d)
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
  call mkorganic('mksrf_organic.nc',xmask,var3d)
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

  write(stdout,*) 'Created VOC/horganic informations...'

  allocate(var2d(jxsg,iysg))
  call mklake('mksrf_lake.nc',xmask,var2d)
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
  call mkurban_param('mksrf_urban.nc',xmask,var4d,var5d,var6d)
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
        istatus = nf90_put_var(ncid, iurb4d(2*i-1), &
                gcvar, istart(1:3), icount(1:3))
        call checkncerr(istatus,__FILE__,__LINE__, 'Error write '//parm4d(i))
        call mypack(var6d(:,:,2,ir,iu,ip4d(parm4d(i))),gcvar)
        istatus = nf90_put_var(ncid, iurb4d(2*i), &
                gcvar, istart(1:3), icount(1:3))
        call checkncerr(istatus,__FILE__,__LINE__, 'Error write '//parm4d(i))
      end do
    end do
  end do
  deallocate(var4d)
  deallocate(var5d)
  deallocate(var6d)

  write(stdout,*) 'Created LAKE/URBAN informations...'

#ifdef CN
  allocate(var2d(jxsg,iysg))
  istart(1) = 1
  icount(1) = ngcells

  call mklightning_init('mksrf_lightning.nc')
  do it = 1 , noleap_yday_3h
    call mklightning(var2d,xmask,it)
    call mypack(var2d,gcvar)
    istart(2) = it
    icount(2) = 1
    istatus = nf90_put_var(ncid, ilightning, gcvar, istart(1:2), icount(1:2))
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write lnfm')
  end do
  call mklightning_close

  call mkpopd_init('mksrf_popd.nc')
  do it = 1 , nyears
    call mkpopd(var2d,xmask,it)
    call mypack(var2d,gcvar)
    istart(2) = it
    icount(2) = 1
    istatus = nf90_put_var(ncid, ipopden, gcvar, istart(1:2), icount(1:2))
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write hdm')
  end do
  call mkpopd_close
  deallocate(var2d)

  write(stdout,*) 'Created FIRE informations...'

#if defined(CN)
  ! q10 soil respiration
  allocate(var2d(jxsg,iysg))
  call mkq10soil('mksrf_q10soil.nc',xmask,var2d)
  call mypack(var2d,gcvar)
  istatus = nf90_put_var(ncid, q10, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write Q10')
  write(stdout,*) 'Created SOIL RESPIRATION informations...'
  ! ndep nitrogen total deposition
  call mkndep('mksrf_ndep.nc',xmask,var2d)
  call mypack(var2d,gcvar)
  istatus = nf90_put_var(ncid, ndep, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write NDEP')
  write(stdout,*) 'Created NITROGEN Deposition informations...'
  deallocate(var2d)
#endif

#ifdef LCH4
  allocate(var3d(jxsg,iysg,4))
  call mkch4topm('mksrf_ch4topm.nc',xmask,var3d)
  call mypack(var3d(:,:,1),gcvar)
  istatus = nf90_put_var(ncid, k, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write K_PAR')
  call mypack(var3d(:,:,2),gcvar)
  istatus = nf90_put_var(ncid, q, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write XM_PAR')
  call mypack(var3d(:,:,3),gcvar)
  istatus = nf90_put_var(ncid, v, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write V_PAR')
  call mypack(var3d(:,:,4),gcvar)
  istatus = nf90_put_var(ncid, maxf, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write MAXF')
  deallocate(var3d)

  ! soil ph
  allocate(var2d(jxsg,iysg))
  call mksoilph('mksrf_soilph.nc',xmask,var2d)
  call mypack(var2d,gcvar)
  istatus = nf90_put_var(ncid, isoilphvar, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write PH')
  deallocate(var2d)

  write(stdout,*) 'Created CH4 informations...'
#endif
#endif

#ifdef VICHYDRO
  allocate(var3d(jxsg,iysg,4))
  call mkvic('mksrf_vic.nc',xmask,var3d)
  call mypack(var3d(:,:,1),gcvar)
  istatus = nf90_put_var(ncid, ibinfl, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write binfl')
  call mypack(var3d(:,:,2),gcvar)
  istatus = nf90_put_var(ncid, ids, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ds')
  call mypack(var3d(:,:,3),gcvar)
  istatus = nf90_put_var(ncid, idsmax, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write dsmax')
  call mypack(var3d(:,:,4),gcvar)
  istatus = nf90_put_var(ncid, iws, gcvar)
  call checkncerr(istatus,__FILE__,__LINE__, 'Error write ws')
  deallocate(var3d)

  write(stdout,*) 'Created VICHYDRO informations...'
#endif

  istatus = nf90_close(ncid)
  call checkncerr(istatus,__FILE__,__LINE__,  &
    'Error close file '//trim(outfile))

#ifdef DYNPFT
  call split_idate(globidate1,y1,mon,day,hour)
  call split_idate(globidate2,y2,mon,day,hour)
  do it = y1 - 1 , y2 + 1
    write(stdout,*) 'Creating year ', it, ' informations...'
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

    ivdims(1) = idims(1)
    ivdims(2) = idims(2)
    ivdims(3) = idims(4)
    istatus = nf90_def_var(ncid, 'pft_2d', regcm_vartype, &
                           ivdims(1:3), ipft2dvar)
    call checkncerr(istatus,__FILE__,__LINE__,  'Error add var pft')
    istatus = nf90_put_att(ncid, ipft2dvar, 'long_name','percent pft')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add pft long_name')
    istatus = nf90_put_att(ncid, ipft2dvar, 'units','%')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add pft units')
    istatus = nf90_put_att(ncid, ipft2dvar, '_FillValue',vmisdat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error add pft Fill Value')

    ! Variables
    istatus = nf90_def_var(ncid, 'gridcell', nf90_int, idims(5), ilndvar)
    call checkncerr(istatus,__FILE__,__LINE__,  'Error add var gridcell')
    istatus = nf90_put_att(ncid, ilndvar, 'compress','xlat xlon')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add gridcell compress')

    istatus = nf90_def_var(ncid, 'xclon', regcm_vartype, idims(5), ilonvar)
    call checkncerr(istatus,__FILE__,__LINE__,  'Error add var xclon')
    istatus = nf90_put_att(ncid, ilndvar, 'compress','xlat xlon')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add xclon')

    istatus = nf90_def_var(ncid, 'yclat', regcm_vartype, idims(5), ilatvar)
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
    istatus = nf90_def_var(ncid, 'PCT_PFT', regcm_vartype, ivdims(1:2), ipftvar)
    call checkncerr(istatus,__FILE__,__LINE__,  'Error add var pft')
    istatus = nf90_put_att(ncid, ipftvar, 'long_name','percent pft')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add pft long_name')
    istatus = nf90_put_att(ncid, ipftvar, 'units','%')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add pft units')
    istatus = nf90_put_att(ncid, ipftvar, '_FillValue',vmisdat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error add pft Fill Value')
#ifdef CN
    istatus = nf90_def_var(ncid, 'HARVEST_VH1', regcm_vartype, &
                           idims(5), iharvvh1)
    call checkncerr(istatus,__FILE__,__LINE__,  'Error add var HARVEST_VH1')
    istatus = nf90_put_att(ncid, iharvvh1, 'long_name', &
            'harvest from primary forest')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add HARVEST_VH1 long_name')
    istatus = nf90_def_var(ncid, 'HARVEST_VH2', regcm_vartype, &
                           idims(5), iharvvh2)
    call checkncerr(istatus,__FILE__,__LINE__,  'Error add var HARVEST_VH2')
    istatus = nf90_put_att(ncid, iharvvh2, 'long_name', &
            'harvest from primary non-forest')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add HARVEST_VH2 long_name')
    istatus = nf90_def_var(ncid, 'HARVEST_SH1', regcm_vartype, &
                           idims(5), iharvsh1)
    call checkncerr(istatus,__FILE__,__LINE__,  'Error add var HARVEST_SH1')
    istatus = nf90_put_att(ncid, iharvsh1, 'long_name', &
            'harvest from secondary mature-forest')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add HARVEST_SH1 long_name')
    istatus = nf90_def_var(ncid, 'HARVEST_SH2', regcm_vartype, &
                           idims(5), iharvsh2)
    call checkncerr(istatus,__FILE__,__LINE__,  'Error add var HARVEST_SH2')
    istatus = nf90_put_att(ncid, iharvsh2, 'long_name', &
            'harvest from secondary young-forest')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add HARVEST_SH2 long_name')
    istatus = nf90_def_var(ncid, 'HARVEST_SH3', regcm_vartype, &
                           idims(5), iharvsh3)
    call checkncerr(istatus,__FILE__,__LINE__,  'Error add var HARVEST_SH3')
    istatus = nf90_put_att(ncid, iharvsh3, 'long_name', &
            'harvest from secondary non-forest')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add HARVEST_SH3 long_name')
    istatus = nf90_def_var(ncid, 'GRAZING', regcm_vartype, idims(5), igraz)
    call checkncerr(istatus,__FILE__,__LINE__,  'Error add var GRAZING')
    istatus = nf90_put_att(ncid, igraz, 'long_name', &
            'grazing of herbacous pfts')
    call checkncerr(istatus,__FILE__,__LINE__,'Error add GRAZING long_name')
#endif

    istatus = nf90_enddef(ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error exit define mode')

    if ( idynamic < 3 ) then
      hptop = real(ptop*10.0_rkx)
      call write_vertical_coord(ncid,rsigx,hptop,izvar)
    else
      call model_zitah(zita)
      ax = md_ak(zita)
      bx = md_bk(zita)
      call write_vertical_coord_sigma(ncid,rsigx,ax,bx,izvar)
    end if
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
    call mkdynpft(xmask,var3d(:,:,:),it)

    do i = igstart , igstop
      do j = jgstart , jgstop
        if ( xmask(j,i) > 0.5_rkx ) then
          if ( pctspec(j,i) > 99.9_rkx ) then
            var3d(j,i,:) = 0.0_rkx
          else
            spft = sum(var3d(j,i,:))
            if ( spft < -1.0_rkx ) then
              call die(__FILE__,'No points around !',__LINE__)
            end if
            diff = spft - (100.0_rkx-pctspec(j,i))
            if ( abs(diff) > 0.5_rkx ) then
              do while ( diff >= 1.0_rkx )
                iloc = maxloc(var3d(j,i,:))
                var3d(j,i,iloc(1)) = nint(var3d(j,i,iloc(1)) - 1._rkx)
                spft = sum(var3d(j,i,:))
                diff = spft - (100.0_rkx-pctspec(j,i))
              end do
              do while ( diff <= -1.0_rkx )
                iloc = maxloc(var3d(j,i,:))
                var3d(j,i,iloc(1)) = nint(var3d(j,i,iloc(1)) + 1._rkx)
                spft = sum(var3d(j,i,:))
                diff = spft - (100.0_rkx-pctspec(j,i))
              end do
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
    istatus = nf90_put_var(ncid, ipft2dvar, var3d)
    call checkncerr(istatus,__FILE__,__LINE__, 'Error write pft2d')
    deallocate(var3d)

#ifdef CN
    allocate(var3d(jxsg,iysg,6))
    call mkharvest(xmask,var3d,it)
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
#endif

    istatus = nf90_close(ncid)
    call checkncerr(istatus,__FILE__,__LINE__,  &
      'Error close file '//trim(outfile))
  end do
#endif

  call memory_destroy

  write(stdout,*) 'Successfully completed CLM preprocessing.'

  contains

  recursive subroutine sortpatch(vals,svals,ird,lsub)
    implicit none
    real(rkx) , dimension(:) , intent(in) :: vals
    real(rkx) , dimension(:) , intent(inout) :: svals
    integer(ik4) , dimension(:) , intent(inout) :: ird
    logical , optional :: lsub
    integer(ik4) :: i , iswap
    real(rkx) :: rswap
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
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
