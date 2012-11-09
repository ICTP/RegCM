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

program clm2rcm
 
  use mod_intkinds
  use mod_realkinds
  use mod_nclib
  use mod_dynparam
  use mod_message
  use mod_grid
  use mod_param_clm
  use mod_date
  use mod_clm3grid
  use mod_memutil
  use mod_stdio
  use mod_domain
  use mod_nchelper
  use netcdf

  implicit none
!
  real(rk4) , parameter :: vmisdat=-9999.0
  integer(ik4) , parameter :: ndim = 3
  logical , parameter :: bvoc = .false.
!
  integer(ik4) :: istatus , ncid , incin , idum
  integer(ik4) , dimension(4) :: idims
  integer(ik4) , dimension(4) :: ivdims
  integer(ik4) :: ldim , ivar
  type(rcm_time_and_date) :: irefdate , imondate
  type(rcm_time_interval) :: tdif
  real(rk4) , pointer , dimension(:) :: yiy
  real(rk4) , pointer , dimension(:) :: xjx
  real(rk4) :: hptop , xmiss
  real(rk8) , dimension(1) :: xdate
  integer(ik4) , dimension(2) :: ihvar
  integer(ik4) , dimension(2) :: illvar
  integer(ik4) , dimension(2) :: izvar
  integer(ik4) , dimension(4) :: icount , istart
  integer(ik4) , dimension(1) :: istart1 , icount1
  integer(ik4) , dimension(3) :: iadim
  character(len=64) , dimension(nfld) :: lnam
  character(len=64) :: cdum
  character(len=64) , dimension(nfld) :: units
  real(rk4) , dimension(3) :: varmax , varmin
  real(rk8) :: xhr
  real(rk4) :: offset , xscale , xlatmin , xlatmax , xlonmin , xlonmax
  real(rk4) :: pxerr , pmax
  real(rk4) , pointer , dimension(:) :: glat , glon , zlat ,      &
                                       zlev , zlon
  real(rk4) , pointer , dimension(:,:) :: mpu
  real(rk4) , pointer , dimension(:,:,:) :: regxyz
  real(rk4) , pointer , dimension(:,:,:,:) :: regyxzt , zoom
  real(rk4) , pointer , dimension(:,:,:) :: save_regyz
  real(rk4) , pointer , dimension(:,:) :: new_forest , old_forest
  logical :: hassmask = .false.
  real(rk4) , pointer , dimension(:,:) :: landmask , sandclay
  integer(ik4) :: ipathdiv , ierr
  integer(ik4) :: i , iz , it , j , k , l , kmax , ipnt
  integer(ik4) :: jotyp , idin , idout , ifield , ifld , imap
  character(len=256) :: namelistfile , prgname
  character(len=256) :: inpfile , terfile , checkfile
  character(len=256) :: outfil_nc
  character(len=64) :: csdate , cldim
  integer(ik4) , dimension(8) :: ilevs
  integer(ik4) , parameter :: iforest = 5
  integer(ik4) , parameter :: igrass  = 14
!
  data ilevs /-1,-1,-1,-1,-1,-1,-1,-1/
  data xmiss /-9999.0/
!
!     Read input global namelist
!
  call getarg(0, prgname)
  call getarg(1, namelistfile)
  call initparam(namelistfile, ierr)
  if ( ierr/=0 ) then
    write(stderr,*) 'Parameter initialization not completed'
    write(stderr,*) 'Usage : '
    write(stderr,*) '          ', trim(prgname), ' regcm.in'
    write(stderr,*) ' '
    call die('clm2rcm','Check argument and namelist syntax',1)
  end if
!
  call memory_init

  if ( nsg/=1 ) then
    write (stderr,*) 'CLM does not work with subgridding enable.'
    write (stderr,*) 'Please set nsg=1 in regcm.in'
    call die('clm2rcm','Check argument and namelist syntax',1)
  end if

  call init_domain
!
!     ** Get latitudes and longitudes from DOMAIN file
!
  terfile = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
  call openfile_withname(terfile,incin)
  call read_domain(incin,sigx,xlat,xlon)
  istatus = nf90_close(incin)
  call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error closing file '//trim(terfile))
!     ** Set output variables
  jotyp = 2
  xscale = 1.
  offset = 0.
 
!     ** Open Output checkfile in NetCDF format

  checkfile = trim(dirglob)//pthsep//trim(domname)//'_CLM3.nc'

  call createfile_withname(checkfile,ncid)
  call add_common_global_params(ncid,'clm2rcm',.false.)
  ipnt = 1
  call define_basic_dimensions(ncid,jx,iy,kzp1,ipnt,idims)
  call add_dimension(ncid,'time',nf90_unlimited,ipnt,idims)

  call define_horizontal_coord(ncid,jx,iy,xjx,yiy,idims,ihvar)
  call define_vertical_coord(ncid,idims,izvar)

  ipnt = 1
  call define_cross_geolocation_coord(ncid,idims,ipnt,illvar)

  istatus = nf90_def_var(ncid, 'time', nf90_double, idims(4:4), ivar)
  call checkncerr(istatus,__FILE__,__LINE__,'Error add variable time')
  irefdate = yrfirst(globidate1)
  irefdate = monmiddle(irefdate)
  csdate = tochar(irefdate)
  istatus = nf90_put_att(ncid, ivar, 'units', 'hours since '//csdate)
  call checkncerr(istatus,__FILE__,__LINE__,'Error add time units')
  istatus = nf90_put_att(ncid, ivar, 'calendar', calendar)
  call checkncerr(istatus,__FILE__,__LINE__,'Error add time calendar')
!
  istatus = nf90_enddef(ncid)
  call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error End Definitions NetCDF output')
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
    istatus = nf90_put_var(ncid, ivar, xdate, istart1, icount1)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable time write')
    imondate = nextmon(imondate)
    imondate = monmiddle(imondate)
  end do

!     ** determine which files to create (emission factor map or not)
  call comp(ifield,bvoc)
 
!     ** Loop over the fields defined in clm.param
  do ifld = 1 , ifield
 
!       ** Open input and output files
!       **   Some files have more than one required variable. 
!       Therefore, **   the output file should only be opened once.
    inpfile = trim(inpglob)//infil(ifld)
    if ( ifld==ipft .or. ifld==ilai .or. ifld==ilak .or.   &
         ifld==iglc .or. ifld==iurb .or. ifld==isnd .or.   &
         ifld==icol .or. ifld==ioro .or. ifld==iiso .or.   &
         ifld==ifma .or. ifld==imbo .or. ifld==ibpin .or.  &
         ifld==iapin ) then
!         ************************ CHANGED LINE ABOVE to include iiso
!         ************************
      write(stdout,*) 'OPENING Input NetCDF FILE: ' , trim(inpfile)
      ierr = nf90_open(inpfile,nf90_nowrite,idin)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Cannot open file '//trim(inpfile))
      ipathdiv = scan(inpfile, pthsep, .true.)
      if ( ipathdiv/=0 ) then
        outfil_nc = trim(dirglob)//pthsep//trim(domname)//  &
                  '_RCM'//inpfile(ipathdiv+7:)
      else
        outfil_nc = trim(dirglob)//pthsep//trim(domname)//  &
                 '_RCM'//inpfile(7:)
      endif
!         CALL FEXIST(outfil_nc)
      write(stdout,*) 'OPENING Output NetCDF FILE: ',trim(outfil_nc)
      call rcrecdf(outfil_nc,idout,varmin,varmax,3,ierr)
    end if
 
!   ** Setup RegCM3 grid variables
    call param(jx,iy,nlev(ifld),xlat,xlon,varmin,varmax,        &
               xlat1d,xlon1d,xlonmin,xlonmax,xlatmin,xlatmax,   &
               iadim,ndim)
 
!   ** Setup CLM3 grid variables, including subgrid region
    call getmem1d(glon,1,nlon(ifld),'clm2rcm:glon')
    call getmem1d(glat,1,nlat(ifld),'clm2rcm:glat')
 
    call clm3grid1(nlon(ifld),nlat(ifld),nlev(ifld),ntim(ifld),     &
                   glon1(ifld),glon2(ifld),glat1(ifld),glat2(ifld), &
                   xlonmin,xlonmax,xlatmin,xlatmax,glon,glat,istart,&
                   icount)
 
    if ( ifld==isnd .or. ifld==icly ) then
      istart(4) = 1
      icount(4) = 1
    end if
 
    if ( ifld==ipft ) then
      call getmem3d(save_regyz,1,iy,1,jx,1,nlev(ifld),'clm2rcm:save_regyz')
      call getmem2d(new_forest,1,iy,1,jx,'clm2rcm:new_forest')
      call getmem2d(old_forest,1,iy,1,jx,'clm2rcm:old_forest')
    end if
    call getmem4d(zoom,1,icount(1),1,icount(2),1,icount(3),1,icount(4), &
                  'clm2rcm:zoom')
    call getmem1d(zlon,1,icount(1),'clm2rcm:zlon')
    call getmem1d(zlat,1,icount(2),'clm2rcm:zlat')
    call getmem1d(zlev,1,icount(3),'clm2rcm:zlev')
    call getmem2d(landmask,1,icount(1),1,icount(2),'clm2rcm:landmask')
!
    call clm3grid2(nlon(ifld),nlat(ifld),glon,glat,istart,icount,zlon,zlat,zlev)
!
    write(stdout,*) 'Reading variables from input file'
!
!       ** Read in the variables.
!       In some cases, special reads need to be performed:
!       - Sand/clay fractions need a mapping variable
!       - Lakes, wetlands, soil color, and orography need a
!       180 degree longitiude shift.
!       - Soil color and Orography do not have landmasks (must be made)
    if ( ifld==isnd .or. ifld==icly ) then
      call getmem2d(sandclay,1,ntim(ifld),1,nlev(ifld),'clm2rcm:sandclay')
      call getmem2d(mpu,1,icount(1),1,icount(2),'clm2rcm:mpu')
      call readcdfr4(idin,vnam(ifld),lnam(ifld),units(ifld),1,      &
                     ntim(ifld),1,nlev(ifld),1,1,1,1,sandclay)
      write(stdout,*) 'Read ', trim(lnam(ifld))
      call readcdfr4(idin,vnam_lm,cdum,cdum,istart(1),              &
                     icount(1),istart(2),icount(2),1,1,1,1,landmask)
      call readcdfr4(idin,vnam_st,cdum,cdum,istart(1),              &
                     icount(1),istart(2),icount(2),1,1,1,1,mpu)
      do j = 1 , icount(2)
        do i = 1 , icount(1)
          imap = nint(mpu(i,j))
          do k = 1 , icount(3)
            if ( imap>0 .and. landmask(i,j)>0.5 ) then
              zoom(i,j,k,1) = sandclay(imap,k)
            else
              zoom(i,j,k,1) = vmisdat
            end if
          end do
        end do
      end do
      do k = 1 , icount(3)
        zlev(k) = glev_st(k)
      end do
      ntim(ifld) = 1
 
    else if ( ifld==iiso .or. ifld==ibpin .or. ifld==iapin .or.     &
              ifld==imbo ) then
      call readcdfr4_iso(idin,vnam(ifld),lnam(ifld),units(ifld),    &
                         istart(1),icount(1),istart(2),icount(2),   &
                         istart(3),icount(3),istart(4),icount(4),   &
                         zoom)
      write(stdout,*) 'Read ', trim(lnam(ifld))
    else
      if ( ifld/=icol ) then
        call readcdfr4(idin,vnam_lm,lnam(ifld),units(ifld),         &
               istart(1),icount(1),istart(2),icount(2),1,1,1,1,     &
               landmask)
        write(stdout,*) 'Read ', trim(lnam(ifld))
      end if
      call readcdfr4(idin,vnam(ifld),lnam(ifld),units(ifld),        &
                     istart(1),icount(1),istart(2),icount(2),       &
                     istart(3),icount(3),istart(4),icount(4),zoom)
      write(stdout,*) 'Read ', trim(lnam(ifld))
    end if
 
    if ( ifld==icol .or. ifld==iiso .or. ifld==ibpin .or.           &
         ifld==imbo .or. ifld==iapin ) then
      write(stdout,*) 'Adjusting landmask'
      do j = 1 , icount(2)
        do i = 1 , icount(1)
          if ( zoom(i,j,1,1)>vmin(ifld) ) then
            landmask(i,j) = 1.0
          else
            landmask(i,j) = vmisdat
          end if
        end do
      end do
    end if

    write(stdout,*) 'READ/WRITE: ', trim(vnam(ifld)), ', ', &
               trim(lnam(ifld)), ', ', trim(units(ifld))
 
!       ** Set the non-land values to missing for interpolation purposes

    if ( ifld/=ioro ) call maskme(landmask,zoom,vmisdat,icount(1),  &
                                  icount(2),icount(3),icount(4))
 
!       ** Interpolate data to RegCM3 grid

    call getmem4d(regyxzt,1,jx,1,iy,1,nlev(ifld),1,ntim(ifld),'clm2rcm:regyxzt')
    call getmem3d(regxyz,1,jx,1,iy,1,nlev(ifld),'clm2rcm:regxyz')

    call bilinx4d(zoom,zlon,zlat,icount(1),icount(2),regyxzt,xlon,  &
                  xlat,jx,iy,icount(3),icount(4),vmin(ifld),vmisdat)
 
!   ** Write the interpolated data to NetCDF for CLM and checkfile
    imondate = irefdate
    if ( ifld == ipft ) then
      inquire(file=trim(dirglob)//pthsep//trim(domname)//'_CLM_SAGEMASK', &
             exist=hassmask)
      if ( hassmask ) then
        open(163,file=trim(dirglob)//pthsep//trim(domname)//'_CLM_SAGEMASK', &
             form='unformatted',status='old')
        read(163) save_regyz
      end if
    end if
    do l = 1 , ntim(ifld)
      if ( ifld==ipft ) then
        do i = 1 , iy
          do j = 1 , jx
            if ( regyxzt(j,i,1,l)>-99. ) then
              do k = 1 , nlev(ifld)
                regyxzt(j,i,k,l) = aint(regyxzt(j,i,k,l))
              end do
              pxerr = 100.
              kmax = -1
              pmax = -99.
              do k = 1 , nlev(ifld)
                pxerr = pxerr - regyxzt(j,i,k,l)
                if ( regyxzt(j,i,k,l)>pmax ) then
                  pmax = regyxzt(j,i,k,l)
                  kmax = k
                end if
              end do
              regyxzt(j,i,kmax,l) = regyxzt(j,i,kmax,l) + pxerr
            end if
          end do
        end do
        if ( .not. hassmask ) then
          save_regyz(:,:,:) = regyxzt(:,:,:,1)
          open(163,file=trim(dirglob)//pthsep//trim(domname)//'_CLM_SAGEMASK', &
               form='unformatted',status='new')
          write(163) save_regyz
        else
          new_forest = regyxzt(:,:,iforest,1)
          old_forest = save_regyz(:,:,iforest)
          regyxzt(:,:,:,1) = save_regyz(:,:,:)
          regyxzt(:,:,iforest,1) = new_forest(:,:)
          regyxzt(:,:,igrass,1) = save_regyz(:,:,igrass)+ &
                           (old_forest(:,:)-new_forest(:,:))
        end if
      end if
      do k = 1 , nlev(ifld)
        do j = 1 , jx
          do i = 1 , iy
            if ( ifld==icol ) then
              regyxzt(j,i,k,l) = anint(regyxzt(j,i,k,l))
            end if
            if ( regyxzt(j,i,k,l)>vmin(ifld) ) then
              regxyz(j,i,k) = regyxzt(j,i,k,l)
            else
              regxyz(j,i,k) = 0.0
            end if
          end do
        end do
      end do
      tdif = imondate-irefdate
      xhr = tohours(tdif)
      call writecdf(idout,vnam(ifld),regxyz,jx,iy,nlev(ifld),iadim, &
                    xhr,lnam(ifld),units(ifld),xscale,offset,varmin,&
                    varmax,xlat1d,xlon1d,zlev,0,vmisdat,jotyp)
      imondate = nextmon(imondate)
      imondate = monmiddle(imondate)
    end do

    istatus = nf90_redef(ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error in file '//trim(checkfile))

    ldim = 0
    if (nlev(ifld) > 1) then
      ldim = -1
      do iz = 1 , 8
        if (ilevs(iz) < 0) then
          ldim = -iz
          exit
        end if
        if (nlev(ifld) == ilevs(iz)) then
          ldim = iz
          exit
        end if
      end do
      if (ldim < 0) then
        ldim = -ldim
        ilevs(ldim) = nlev(ifld)
        write (cldim,'(a,i0.3)') 'level_', nlev(ifld)
        istatus = nf90_def_dim(ncid, cldim, nlev(ifld), idum)
        call checkncerr(istatus,__FILE__,__LINE__,  &
          'Error create dim '//trim(cldim))
      end if
    end if

    if (vnam(ifld)(1:8) == 'MONTHLY_') then
      vnam(ifld) = vnam(ifld)(9:)
    end if

    if (ntim(ifld) > 1) then
      if (nlev(ifld) > 1) then
        ivdims(1) = idims(1)
        ivdims(2) = idims(2)
        ivdims(3) = idims(4)+ldim
        ivdims(4) = idims(3)
        istatus = nf90_def_var(ncid, vnam(ifld), nf90_float, ivdims, ivar)
      else
        ivdims(1) = idims(1)
        ivdims(2) = idims(2)
        ivdims(3) = idims(3)
        istatus = nf90_def_var(ncid, vnam(ifld), nf90_float, ivdims(1:3), ivar)
      end if
    else
      if (nlev(ifld) > 1) then
        ivdims(1) = idims(1)
        ivdims(2) = idims(2)
        ivdims(3) = idims(4)+ldim
        istatus = nf90_def_var(ncid, vnam(ifld), nf90_float, ivdims(1:3), ivar)
      else
        ivdims(1) = idims(1)
        ivdims(2) = idims(2)
        istatus = nf90_def_var(ncid, vnam(ifld), nf90_float, ivdims(1:2), ivar)
      end if
    end if
    call checkncerr(istatus,__FILE__,__LINE__,  &
      'Error add var '//trim(vnam(ifld)))
    idum = len_trim(lnam(ifld))
    istatus = nf90_put_att(ncid, ivar, 'long_name', lnam(ifld)(1:idum))
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error add '//trim(vnam(ifld))//' long_name')
    idum = len_trim(units(ifld))
    istatus = nf90_put_att(ncid, ivar, 'units', units(ifld)(1:idum))
    call checkncerr(istatus,__FILE__,__LINE__, &
                  ('Error add '//trim(vnam(ifld))//' units'))
    istatus = nf90_put_att(ncid, ivar, '_FillValue', xmiss)
    call checkncerr(istatus,__FILE__,__LINE__, &
                  ('Error add '//trim(vnam(ifld))//' xmiss'))
    istatus = nf90_put_att(ncid, ivar, 'coordinates', 'xlon xlat')
    call checkncerr(istatus,__FILE__,__LINE__, &
                  ('Error add '//trim(vnam(ifld))//' coords'))

    istatus = nf90_enddef(ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error End redef NetCDF output')

    istatus = nf90_put_var(ncid, ivar, regyxzt)
    call checkncerr(istatus,__FILE__,__LINE__, &
      ('Error '//trim(vnam(ifld))// ' write'))

    if ( associated(save_regyz) ) call relmem3d(save_regyz)

  end do  ! End nfld loop

!     ** Close files

  call clscdf(idin,ierr)
  call clscdf(idout,ierr)

  istatus = nf90_close(ncid)
  call checkncerr(istatus,__FILE__,__LINE__,  &
    'Error close file '//trim(checkfile))

  call memory_destroy

  write(stdout,*) 'Successfully completed CLM preprocessing.'
 
end program clm2rcm
