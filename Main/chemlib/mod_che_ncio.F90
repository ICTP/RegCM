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
!
module mod_che_ncio
!
  use mod_realkinds
  use mod_dynparam
  use mod_mpmessage
  use mod_che_indices
  use mod_che_common
  use netcdf
!
  private
!
  public :: read_texture , read_aerosol , read_emission , recc
  public :: prepare_chem_out , init_mod_che_ncio , writerec_che2
  public :: open_chbc , close_chbc , chbc_search , read_chbc

  public :: chbc_ivar

  integer :: istatus
  integer :: recc

  integer , parameter :: n_chevar = 18
  integer , parameter :: n_optvar = 9
  integer , parameter :: n_chbcvar = 25
  integer :: ichin 
  integer :: ioxcl     

  integer , dimension(:) , allocatable :: ncche     
  integer , dimension(n_chevar) :: ichevar
  integer , dimension(n_optvar) ::ioptvar 

  character(len=8) , dimension(n_chbcvar) :: chbcname
  integer , dimension(n_chbcvar) :: chbc_ivar
  
  type(rcm_time_and_date) , dimension(:) , allocatable :: oxcl_idate
  type(rcm_time_and_date) , dimension(:) , allocatable :: chbc_idate
  type(rcm_time_and_date) , save :: icherefdate
  integer , dimension(9) :: idims 
  integer ::idmin , icherec , ioptrec
  integer :: ibcrec , ibcnrec
  real(dp) :: tpd, cfd
  real(dp) :: rpt

  integer :: o_is
  integer :: o_ie
  integer :: o_js
  integer :: o_je
  integer :: o_ni
  integer :: o_nj
  integer :: o_nz
  logical :: lwrap  
  character(256) :: dname , icbcname
  real(sp) , dimension(:) , pointer :: hsigma
  real(sp) , dimension(:,:) , pointer :: ioxlat
  real(sp) , dimension(:,:) , pointer :: ioxlon
  real(sp) , dimension(:,:) , pointer :: iotopo
  real(sp) , dimension(:,:) , pointer :: iomask
  real(sp) , dimension(:,:,:) , pointer :: dumio

  real(sp) , dimension(:,:) , pointer :: sp2d
  real(sp) , dimension(:,:) , pointer :: iolnds

  data ichin   /-1/
  data ioxcl   /-1/
  data ibcrec  / 1/
  data ibcnrec / 0/

  data chbcname /'O3      ','NO      ','NO2     ','HNO3    ', &
                 'N2O5    ','H2O2    ','CH4     ','CO      ', &
                 'CH2O    ','CH3OH   ','C2H5OH  ','C2H4    ', &
                 'C2H6    ','CH3CHO  ','CH3COCH3','BIGENE  ', &
                 'BIGALK  ','C3H6    ','C3H8    ','ISOP    ', &
                 'TOLUE   ','PAN     ','SO2     ','SO4     ', &
                 'DMS     '/

  contains

    subroutine init_mod_che_ncio
      implicit none
      dname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
      if (lcband) then
        o_is = 2
        o_ie = iy-1
        o_js = 1
        o_je = jx
        o_ni = iy-2
        o_nj = jx
        o_nz = kz
        lwrap = .true.
      else
        o_is = 2
        o_ie = iy-1
        o_js = 2
        o_je = jx-1
        o_ni = iy-2
        o_nj = jx-2
        o_nz = kz
        lwrap = .false.
      end if

      call getmem1d(hsigma,1,o_nz,'ncio:hsigma')
      call getmem2d(ioxlat,1,o_nj,1,o_ni,'ncio:ioxlat')
      call getmem2d(ioxlon,1,o_nj,1,o_ni,'ncio:ioxlon')
      call getmem2d(iotopo,1,o_nj,1,o_ni,'ncio:iotopo')
      call getmem2d(iomask,1,o_nj,1,o_ni,'ncio:iomask')
      call getmem3d(dumio,1,o_nj,1,o_ni,1,o_nz,'ncio:dumio')
      call getmem2d(sp2d,1,jx,1,iy,'ncio:sp2d')
      call getmem2d(iolnds,1,o_nj,1,o_ni,'ncio:iolnds')
    end subroutine init_mod_che_ncio


! IMPORTANT : note the 2 following routine have to be duplicated in the che_ncio_module 
! one day think about merging chem_ncio and ncio module ( Graziano would like that )  
 subroutine open_domain
    use netcdf
    implicit none

    real(dp)                   :: dx
    real(dp) , dimension(kzp1) :: sigma

    integer :: ivarid , idimid
    integer :: iyy , jxx , kzz , k
    character(6) :: proj
    real(sp) :: dsx , iclat , iclon , ptsp
    real(sp) , dimension(kzp1) :: rsdum

    write (aline,*) 'open_domain: READING HEADER FILE:', dname
    call say
    istatus = nf90_open(dname, nf90_nowrite, idmin)
    call check_ok(__FILE__,__LINE__,'Error Opening Domain file '//trim(dname), &
                  'DOMAIN FILE')
!!$    if ( nsg > 1 ) then
!!$      write (aline,*) 'READING HEADER SUBDOMAIN FILE:', sdname
!!$      call say
!!$      istatus = nf90_open(sdname, nf90_nowrite, isdmin)
!!$      call check_ok(__FILE__,__LINE__, &
!!$           'Error Opening SubDomain file '//trim(sdname), 'SUBDOM FILE')
!!$    end if
    istatus = nf90_inq_dimid(idmin, 'iy', idimid)
    call check_ok(__FILE__,__LINE__,'Dimension iy miss', 'DOMAIN FILE')
    istatus = nf90_inquire_dimension(idmin, idimid, len=iyy)
    call check_ok(__FILE__,__LINE__,'Dimension iy read error', 'DOMAIN FILE')
    istatus = nf90_inq_dimid(idmin, 'jx', idimid)
    call check_ok(__FILE__,__LINE__,'Dimension jx miss', 'DOMAIN FILE')
    istatus = nf90_inquire_dimension(idmin, idimid, len=jxx)
    call check_ok(__FILE__,__LINE__,'Dimension jx read error', 'DOMAIN FILE')
    istatus = nf90_inq_dimid(idmin, 'kz', idimid)
    call check_ok(__FILE__,__LINE__,'Dimension kz miss', 'DOMAIN FILE')
    istatus = nf90_inquire_dimension(idmin, idimid, len=kzz)
    call check_ok(__FILE__,__LINE__,'Dimension kz read error', 'DOMAIN FILE')
    istatus = nf90_inq_varid(idmin, 'ptop', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable ptop miss', 'DOMAIN FILE')
    istatus = nf90_get_var(idmin, ivarid, ptsp)
    call check_ok(__FILE__,__LINE__,'Variable ptop read error', 'DOMAIN FILE')
    istatus = nf90_get_att(idmin, nf90_global, 'projection', proj)
    call check_ok(__FILE__,__LINE__,'Attribute projection miss', &
                  'DOMAIN FILE')
    istatus = nf90_get_att(idmin, nf90_global,'grid_size_in_meters', dsx)
    call check_ok(__FILE__,__LINE__,'Attribute gridsize miss','DOMAIN FILE')
    istatus = nf90_get_att(idmin, nf90_global, &
                           'latitude_of_projection_origin', iclat)
    call check_ok(__FILE__,__LINE__,'Attribute clat miss', 'DOMAIN FILE')
    istatus = nf90_get_att(idmin, nf90_global, &
                           'longitude_of_projection_origin', iclon)
    call check_ok(__FILE__,__LINE__,'Attribute clon miss', 'DOMAIN FILE')
!
!         Consistency Check
!
    if ( iyy /= iy .or. jxx /= jx .or. kzz /= kzp1 ) then
      write (6,*) 'Error: dims from regcm.in and DOMAIN file differ.'
      write (aline,*) 'Input namelist : IY=', iy , '  JX=', jx , '  KZ=', kz
      call say
      write (aline,*) 'DOMAIN file    : IY=', iyy ,'  JX=', jxx, '  KZ=', kzz-1
      call say
      call fatal(__FILE__,__LINE__,'DIMENSION MISMATCH')
    end if
    if (dabs(dble(ptsp*d_r10)-dble(ptop)) > 0.001D+00) then
      write (6,*) 'Error: ptop from regcm.in and DOMAIN file differ.'
      write (6,*) 'Input namelist = ', ptop
      write (6,*) 'DOMAIN file    = ', ptsp*d_r10
      call fatal(__FILE__,__LINE__, 'DOMAIN ptop')
    end if
    if (proj /= iproj) then
      write (6,*) 'Error: proj from regcm.in and DOMAIN file differ.'
      write (6,*) 'Input namelist = ', iproj
      write (6,*) 'DOMAIN file    = ', proj
      call fatal(__FILE__,__LINE__, 'DOMAIN proj')
    end if
    if (dabs(dble(dsx*d_r1000)-dble(ds)) > 0.001D+00) then
      write (6,*) 'Error: ds from regcm.in and DOMAIN file differ.'
      write (6,*) 'Input namelist = ', ds
      write (6,*) 'DOMAIN file    = ', dsx*d_r1000
      call fatal(__FILE__,__LINE__, 'DOMAIN ds')
    end if
    if (dabs(dble(iclat)-dble(clat)) > 0.001D+00) then
      write (6,*) 'Error: clat from regcm.in and DOMAIN file differ.'
      write (6,*) 'Input namelist = ', clat
      write (6,*) 'DOMAIN file    = ', iclat
      call fatal(__FILE__,__LINE__, 'DOMAIN clat')
    end if
    if (dabs(dble(iclon)-dble(clon)) > 0.001D+00) then
      write (6,*) 'Error: clon from regcm.in and DOMAIN file differ.'
      write (6,*) 'Input namelist = ', clon
      write (6,*) 'DOMAIN file    = ', iclon
      call fatal(__FILE__,__LINE__, 'DOMAIN clon')
    end if
!
!         Assign values in the top data modules
!
    tpd = houpd/atmfrq
    cfd = houpd/chemfrq
    dx = dble(dsx)
    istatus = nf90_inq_varid(idmin, 'sigma', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable sigma miss', 'DOMAIN FILE')
    istatus = nf90_get_var(idmin, ivarid, rsdum)
    call check_ok(__FILE__,__LINE__,'Variable sigma read error','DOMAIN FILE')
    sigma = dble(rsdum)
    do k = 1 , kz
      hsigma(k) = real((sigma(k)+sigma(k+1))/2.0D0)
    end do

  end subroutine open_domain

  subroutine read_domain
  !subroutine read_domain(ht,lnd,xlat,xlon,xmap,dmap,f)
  ! I shut down passing the argument for now , but maybe could be usefull for chemistry module   

   use netcdf
    implicit none
 !real(dp) , dimension(jx,iy), intent(out)  :: ht
 !....
    real(dp) , dimension(jx,iy)  :: ht
    real(dp) , dimension(jx,iy)  :: lnd
    real(dp) , dimension(jx,iy)  :: xlat
    real(dp) , dimension(jx,iy)  :: xlon
    real(dp) , dimension(jx,iy)  :: xmap
    real(dp) , dimension(jx,iy)  :: dmap
    real(dp) , dimension(jx,iy)  :: f

    integer :: ivarid

    if (idmin < 0) then
      write (6,*) 'Error : Domain file not in open state'
      call fatal(__FILE__,__LINE__, 'DOMAIN FILE')
    end if

    istatus = nf90_inq_varid(idmin, 'topo', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable topo miss', 'DOMAIN FILE')
    istatus = nf90_get_var(idmin, ivarid, sp2d)
    call check_ok(__FILE__,__LINE__,'Variable topo read error', 'DOMAIN FILE')
    ht = dble(sp2d)
    iotopo = sp2d(o_js:o_je,o_is:o_ie)
    istatus = nf90_inq_varid(idmin, 'landuse', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable landuse miss','DOMAIN FILE')
    istatus = nf90_get_var(idmin, ivarid, sp2d)
    call check_ok(__FILE__,__LINE__,'Variable landuse read error','DOMAIN FILE')
    lnd = dble(sp2d)
    iolnds = sp2d(o_js:o_je,o_is:o_ie)
    istatus = nf90_inq_varid(idmin, 'xlat', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable xlat miss', 'DOMAIN FILE')
    istatus = nf90_get_var(idmin, ivarid, sp2d)
    call check_ok(__FILE__,__LINE__,'Variable xlat read error', 'DOMAIN FILE')
    xlat = dble(sp2d)
    ioxlat = sp2d(o_js:o_je,o_is:o_ie)
    istatus = nf90_inq_varid(idmin, 'xlon', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable xlon miss', 'DOMAIN FILE')
    istatus = nf90_get_var(idmin, ivarid, sp2d)
    call check_ok(__FILE__,__LINE__,'Variable xlon read error', 'DOMAIN FILE')
    xlon = dble(sp2d)
    ioxlon = sp2d(o_js:o_je,o_is:o_ie)
    istatus = nf90_inq_varid(idmin, 'xmap', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable xmap miss', 'DOMAIN FILE')
    istatus = nf90_get_var(idmin, ivarid, sp2d)
    call check_ok(__FILE__,__LINE__,'Variable xmap read error','DOMAIN FILE')
    xmap = dble(sp2d)
    istatus = nf90_inq_varid(idmin, 'dmap', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable dmap miss','DOMAIN FILE')
    istatus = nf90_get_var(idmin, ivarid, sp2d)
    call check_ok(__FILE__,__LINE__,'Variable dmap read error', 'DOMAIN FILE')
    dmap = dble(sp2d)
    istatus = nf90_inq_varid(idmin, 'coriol', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable coriol miss', 'DOMAIN FILE')
    istatus = nf90_get_var(idmin, ivarid, sp2d)
    call check_ok(__FILE__,__LINE__,'Variable coriol read error','DOMAIN FILE')
    f = dble(sp2d)
    istatus = nf90_inq_varid(idmin, 'mask', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable mask miss', 'DOMAIN FILE')
    istatus = nf90_get_var(idmin, ivarid, sp2d)
    call check_ok(__FILE__,__LINE__,'Variable mask read error','DOMAIN FILE')
    iomask = sp2d(o_js:o_je,o_is:o_ie)
  end subroutine read_domain

    subroutine close_domain
      implicit none
      if (idmin >= 0) then
        istatus = nf90_close(idmin)
        call check_ok(__FILE__,__LINE__,'Domain file close error','DOMAIN FILE')
        idmin = -1
      end if
    end subroutine close_domain

    subroutine read_texture(nats,texture)
      implicit none
      integer , intent(in) :: nats
      real(dp) , pointer , dimension(:,:,:) , intent(out) :: texture

      integer :: ivarid
      integer :: i , j , n
      integer :: idmin
      integer , dimension(3) :: istart , icount
      character(256) :: dname
      real(sp), dimension(jx,iy) ::  toto

      dname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
      istatus = nf90_open(dname, nf90_nowrite, idmin)
      call check_ok(__FILE__,__LINE__, &
                  'Error Opening Domain file '//trim(dname),'DOMAIN FILE OPEN')
      istatus = nf90_inq_varid(idmin, 'texture_fraction', ivarid)
      call check_ok(__FILE__,__LINE__,'Variable texture_fraction miss', &
                  'DOMAIN FILE')
      istart(2) = 1
      istart(1) = 1
      icount(3) = 1
      icount(2) = iy
      icount(1) = jx
      do n = 1 , nats
        istart(3) = n
        istatus = nf90_get_var(idmin, ivarid, toto, istart, icount)
        call check_ok(__FILE__,__LINE__,'Variable texture_frac read error', &
                      'DOMAIN FILE')
        do i = 1 , iy
          do j = 1 , jx
            texture(j,i,n) = dble(toto(j,i))*0.01D0
            if ( texture(j,i,n) < d_zero ) texture(j,i,n) = d_zero
          end do
        end do
      end do
      istatus = nf90_close(idmin)
      call check_ok(__FILE__,__LINE__,'Domain file close error','DOMAIN FILE')
    end subroutine read_texture

    subroutine read_aerosol(chtrname,chemsrc)
      implicit none
      character(256) :: aername
      character(5) , dimension(ntr) , intent(in) :: chtrname
      real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: chemsrc

      integer :: ncid , ivarid
      real(sp) , dimension(jx,iy) :: toto
      character(5) :: aerctl
      integer , dimension(3) :: istart , icount
      integer :: itr , i , j , m

      aername = trim(dirglob)//pthsep//trim(domname)//'_AERO.nc'
      istatus = nf90_open(aername, nf90_nowrite, ncid)
      call check_ok(__FILE__,__LINE__, &
           'Error Opening Aerosol file '//trim(aername),'AEROSOL FILE OPEN')

      do itr = 1 , ntr
        aerctl = chtrname(itr)
        write (aline, *) itr , aerctl
        call say
        if ( aerctl(1:4) /= 'DUST') then
          if ( aerctl(1:3) == 'SO2' ) then
            if ( aertyp(4:4) == '1' ) then
              istatus = nf90_inq_varid(ncid, 'so2', ivarid)
              call check_ok(__FILE__,__LINE__, &
                            'Variable so2 miss','AEROSOL FILE')
              istatus = nf90_get_var(ncid, ivarid, toto)
              call check_ok(__FILE__,__LINE__, &
                            'Variable so2 read error','AEROSOL FILE')
              do m = 1 , 12
                do i = 1 , iy
                  do j = 1 , jx
                    chemsrc(j,i,m,itr) = dble(toto(j,i))
                  end do
                end do
              end do
            end if
            if ( aertyp(5:5) == '1' ) then
              istatus = nf90_inq_varid(ncid, 'so2_monthly', ivarid)
              call check_ok(__FILE__,__LINE__, &
                            'Variable so2_mon miss','AEROSOL FILE')
              istart(1) = 1
              istart(2) = 1
              icount(1) = jx
              icount(2) = iy
              icount(3) = 1
              do m = 1 , 12
                istart(3) = m
                istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
                call check_ok(__FILE__,__LINE__, &
                              'Variable so2_mon read err','AEROSOL FILE')
                do i = 1 , iy
                  do j = 1 , jx
                    chemsrc(j,i,m,itr) = chemsrc(j,i,m,itr) + dble(toto(j,i))
                  end do
                end do
              end do
            end if
          else if ( aerctl(1:2) == 'BC' ) then
            if ( aertyp(4:4) == '1' ) then
              istatus = nf90_inq_varid(ncid, 'bc', ivarid)
              call check_ok(__FILE__,__LINE__, &
                            'Variable bc miss','AEROSOL FILE')
              istatus = nf90_get_var(ncid, ivarid, toto)
              call check_ok(__FILE__,__LINE__, &
                            'Variable bc read error','AEROSOL FILE')
              do m = 1 , 12
                do i = 1 , iy
                  do j = 1 , jx
                    chemsrc(j,i,m,itr) = dble(toto(j,i))
                  end do
                end do
              end do
            end if
            if ( aertyp(5:5) == '1' ) then
              istatus = nf90_inq_varid(ncid, 'bc_monthly', ivarid)
              call check_ok(__FILE__,__LINE__, &
                            'Variable bc_mon miss','AEROSOL FILE')
              istart(1) = 1
              istart(2) = 1
              icount(1) = jx
              icount(2) = iy
              icount(3) = 1
              do m = 1 , 12
                istart(3) = m
                istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
                call check_ok(__FILE__,__LINE__, &
                              'Variable bc_mon read err','AEROSOL FILE')
                do i = 1 , iy
                  do j = 1 , jx
                    chemsrc(j,i,m,itr) = chemsrc(j,i,m,itr) + dble(toto(j,i))
                  end do
                end do
              end do
            end if
          else if ( aerctl(1:2) == 'OC' ) then
            if ( aertyp(4:4) == '1' ) then
              istatus = nf90_inq_varid(ncid, 'oc', ivarid)
              call check_ok(__FILE__,__LINE__, &
                            'Variable oc miss','AEROSOL FILE')
              istatus = nf90_get_var(ncid, ivarid, toto)
              call check_ok(__FILE__,__LINE__, &
                            'Variable oc read error','AEROSOL FILE')
              do m = 1 , 12
                do i = 1 , iy
                  do j = 1 , jx
                    chemsrc(j,i,m,itr) = dble(toto(j,i))
                  end do
                end do
              end do
            end if
            if ( aertyp(5:5) == '1' ) then
              istatus = nf90_inq_varid(ncid, 'oc_monthly', ivarid)
              call check_ok(__FILE__,__LINE__, &
                            'Variable oc_mon miss','AEROSOL FILE')
              istart(1) = 1
              istart(2) = 1
              icount(1) = jx
              icount(2) = iy
              icount(3) = 1
              do m = 1 , 12
                istart(3) = m
                istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
                call check_ok(__FILE__,__LINE__, &
                              'Variable oc_mon read err','AEROSOL FILE')
                do i = 1 , iy
                  do j = 1 , jx
                    chemsrc(j,i,m,itr) = chemsrc(j,i,m,itr) + dble(toto(j,i))
                  end do
                end do
              end do
            end if
          end if
        end if
      end do

      istatus = nf90_close(ncid)
      call check_ok(__FILE__,__LINE__, &
                    'Error Close Aerosol file '//trim(aername), &
                    'AEROSOL FILE CLOSE')

    end subroutine read_aerosol

    subroutine read_emission(lmonth,echemsrc)
      implicit none
      integer , intent(in) :: lmonth
      real(dp) , pointer , dimension(:,:,:,:) , intent(inout) :: echemsrc
      character(256) :: aername
      integer :: ncid 
      integer , dimension(3) :: istart , icount
      integer :: ivarid
 
! FAB: remember for now, we have 1 emission file containing all monthly
! emission for the whole simulation period
! change that in the future. Also lmonth is not really necessary here,
! but KEEP THIS DIMENSION FOR HIGHER TEMPORAL RESOLUTION INVENTORIES 
! all aggregations / lumping should in the future be done in emission preproc 

      aername = trim(dirglob)//pthsep//trim(domname)//'_CHEMISS.nc'
      print *, 'Opening ch. emission file ', aername

      istatus = nf90_open(aername, nf90_nowrite, ncid)
      call check_ok(__FILE__,__LINE__, &
                    'Error Opening chem emissiom file '//trim(aername), &
                    'CHE EMISS FILE OPEN ERROR')

      !*** intialized in start_chem
      !*** Advice record counter
      recc = recc + 1
      istart(1) = 1
      istart(2) = 1
      istart(3) = recc
      icount(1) = jx
      icount(2) = iy
      icount(3) = 1

!FAB VERY IMPORTANT : THIS READING SECTION SHOULD BE FIXED WITH
!                     HOMOGENEOUS PREPROC

      ! FAB pour esource amma
      if ( igaschem == 1 ) then  
        ! NO emission                  
        if ( ino /= 0 ) then
          call rvar(ncid,istart,icount,ino,lmonth,echemsrc, &
                    'a_NOX',.false.,'b_NOX')
        end if
        ! CO emission
        if ( ico /= 0 ) then
          call rvar(ncid,istart,icount,ico,lmonth,echemsrc, &
                    'a_CO',.false.,'b_CO')
          print*, 'FAB emis testco','ico', maxval(echemsrc)
        end if
        ! HCHO emission                  
        if ( ihcho /= 0 ) then
          ! call rvar(ncid,istart,icount,ihcho,lmonth,echemsrc,'a_HCHO',.false.)
        end if
        ! ACET emission                  
        if ( iacet /= 0 ) then
          call rvar(ncid,istart,icount,iacet,lmonth,echemsrc, &
                    'a_acet',.false.,'b_acet')
        end if
        ! SO2 emission
        if ( iso2 /= 0 ) then
          call rvar(ncid,istart,icount,iso2,lmonth,echemsrc, &
                    'a_so2',.false.)
        end if
        ! CH4
        if ( ich4 /= 0 ) then
          call rvar(ncid,istart,icount,ich4,lmonth,echemsrc, &
                    'a_ch4',.false.,'b_CH4')
        end if
        ! Ethane
        if ( ic2h6 /= 0 ) then
          call rvar(ncid,istart,icount,ic2h6,lmonth,echemsrc, &
                    'a_ETHANE',.false.,'b_ETHANE')
        end if
        ! PAR
        if ( ipar /= 0 ) then
          ! call rvar(ncid,istart,icount,ipar,lmonth,echemsrc, &
          !           'a_c3h8',.false.,'a_butane','bio_c3h8','o_c3h8')
        end if
        ! Ethene
        if ( iethe /= 0 ) then
          call rvar(ncid,istart,icount,iethe,lmonth,echemsrc, &
                    'a_ETHENE',.false.,'b_ETHENE')
        end if
        ! Termenal Alkene
        if ( iolt /= 0 ) then
          call rvar(ncid,istart,icount,iolt,lmonth,echemsrc, &
                    'a_PROPENE',.false.,'b_PROPENE')
        end if
        ! Internal Alkene
        if ( ioli /= 0 ) then
          ! call rvar(ncid,istart,icount,ioli,lmonth,echemsrc,'a_BIGENE',.true.)
        end if
        ! Isoprene
        if ( iisop /= 0 ) then
          call rvar(ncid,istart,icount,iisop,lmonth,echemsrc,'bio_isop',.false.)
        end if
        ! Toluene
        if ( itolue /= 0 ) then
          istatus = nf90_inq_varid(ncid, 'a_XYLENE', ivarid)
          if ( istatus == nf90_noerr ) then
            call rvar(ncid,istart,icount,itolue,lmonth,echemsrc, &
                      'a_TOLUENE',.true.,'b_TOLUENE')
          end if
        end if
        ! Xylene
        if ( ixyl /= 0 ) then
          istatus = nf90_inq_varid(ncid, 'a_TOLUENE', ivarid)
          if ( istatus == nf90_noerr ) then
            call rvar(ncid,istart,icount,ixyl,lmonth,echemsrc, &
                      'a_XYLENE',.true.)
          else
            call rvar(ncid,istart,icount,ixyl,lmonth,echemsrc, &
                      'a_TOLUENE',.true.,'b_XYLENE')
          end if
        end if
        ! Acetaldehyde
        if ( iald2 /= 0 ) then
          ! call rvar(ncid,istart,icount,iald2,lmonth,echemsrc,'a_ALD2',.false.)
        end if
        ! Methanol + Ethanol
        if ( imoh /= 0 ) then
          call rvar(ncid,istart,icount,imoh,lmonth,echemsrc, &
                    'a_ch3oh',.false.,'b_ch3oh','bio_ch3oh')
        end if           
        ! DMS
        if ( idms /= 0 ) then
          ! call rvar(ncid,istart,icount,idms,lmonth,echemsrc,'o_DMS',.false.)
        end if
        ! OC and BC anthropogenic + biomass burning
        if ( ibchb /= 0 ) then
          call rvar(ncid,istart,icount,ibchb,lmonth,echemsrc, &
                    'a_bc',.false.)
        end if
        if ( iochb /= 0 ) then
          call rvar(ncid,istart,icount,iochb,lmonth,echemsrc, &
                    'a_oc',.false.)
        end if
      else 
        ! uniquement les sources BB AMMA
        ! OC and BC anthropogenic + biomass burning
        if ( ibchb /= 0 ) then
          call rvar(ncid,istart,icount,ibchb,lmonth,echemsrc, &
                    'b_BC',.false.)
        end if
        if ( iochb /= 0 ) then
          call rvar(ncid,istart,icount,iochb,lmonth,echemsrc, &
                    'b_OC',.false.)
        end if
        if ( iso2 /= 0 ) then
          call rvar(ncid,istart,icount,iso2,lmonth,echemsrc, &
                    'b_SO2',.false.)
        end if
      end if
      istatus = nf90_close(ncid)
      call check_ok(__FILE__,__LINE__, &
                    'Error Closing Chem emission file '//trim(aername), &
                    'CH EMISS FILE CLOSE ERROR')
    end subroutine read_emission

    subroutine rvar(ncid,istart,icount,ind,lmonth,echemsrc,cna,lh,cnb,cnc,cnd)
      implicit none
      integer , intent(in) :: ncid
      integer , dimension(3) , intent(in) :: istart , icount
      integer , intent(in) :: lmonth
      real(dp) , pointer , dimension(:,:,:,:) , intent(inout) :: echemsrc
      logical , intent(in) :: lh
      character(len=*) , intent(in) :: cna
      character(len=*) , intent(in) , optional :: cnb
      character(len=*) , intent(in) , optional :: cnc
      character(len=*) , intent(in) , optional :: cnd
      integer :: ivarid 
      real(sp) , dimension(jx,iy) :: toto
      integer :: i , j , ind

      istatus = nf90_inq_varid(ncid, cna, ivarid)     
      call check_ok(__FILE__,__LINE__, &
                    'Variable '//cna//' miss','CHEM_EMISS FILE')
      istatus = nf90_get_var(ncid,ivarid,toto)
      call check_ok(__FILE__,__LINE__, &
                    'Variable '//cna//' read err','CHEM_EMISS FILE')
      if ( lh ) then  ! half of lumped Aromatics
        do i = 1 , iy
          do j = 1 , jx
            echemsrc(j,i,lmonth,ind) = d_half*toto(j,i)
          end do
        end do
      else
        do i = 1 , iy
          do j = 1 , jx
            echemsrc(j,i,lmonth,ind) = toto(j,i)
          end do
        end do
      end if
      if ( present(cnb) ) then
        istatus = nf90_inq_varid(ncid, cnb, ivarid)
        call check_ok(__FILE__,__LINE__, &
                      'Variable '//cnb//' miss','CHEM_EMISS FILE')
        istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
        call check_ok(__FILE__,__LINE__, &
                      'Variable '//cnb//' read err','CHEM_EMISS FILE')
        do i = 1 , iy
          do j = 1 , jx
            echemsrc(j,i,lmonth,ind) = toto(j,i) + echemsrc(j,i,lmonth,ind)
          end do
        end do
      end if
      if ( present(cnc) ) then
        istatus = nf90_inq_varid(ncid, cnc, ivarid)
        call check_ok(__FILE__,__LINE__, &
                      'Variable '//cnc//' miss','CHEM_EMISS FILE')
        istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
        call check_ok(__FILE__,__LINE__, &
                      'Variable '//cnc//' read err','CHEM_EMISS FILE')
        do i = 1 , iy
          do j = 1 , jx
            echemsrc(j,i,lmonth,ind) = toto(j,i) + echemsrc(j,i,lmonth,ind)
          end do
        end do
      end if
      if ( present(cnd) ) then
        istatus = nf90_inq_varid(ncid, cnd, ivarid)
        call check_ok(__FILE__,__LINE__, &
                      'Variable '//cnd//' miss','CHEM_EMISS FILE')
        istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
        call check_ok(__FILE__,__LINE__, &
                      'Variable '//cnd//' read err','CHEM_EMISS FILE')
        do i = 1 , iy
          do j = 1 , jx
            echemsrc(j,i,lmonth,ind) = toto(j,i) + echemsrc(j,i,lmonth,ind)
          end do
        end do
      end if
    end subroutine rvar

!------------------------------------------------------------------------------      
!       IROUTINE: prepare_chem_out
!       SUBROUTINE INTERFACE:

    subroutine prepare_chem_out(idate, ifrest)
    
!   !DESCRIPTION:
!   prepare the dimensions variables, write global attributes
!   define the chemistry variables

    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    logical, intent(in) :: ifrest
    integer  :: itr
    character(128)::cdum
    character(8) ::chevarnam
    character(64) :: title
    character(32) :: fbname
    character(16) :: fterr
    character(256) :: ofname , history
    integer , dimension(8) :: tvals
    real(sp) :: hptop
    real(sp) , dimension(2) :: trlat
    real(sp) , dimension(iysg) :: yiy
    real(sp) , dimension(jxsg) :: xjx
    integer :: ncid
    integer , dimension(3) :: izvar
    integer , dimension(2) :: ivvar
    integer , dimension(5) :: illtpvar
    integer :: itvar , i , j
    integer :: ibin , jbin , noutf

    integer , dimension(9) :: tyx
    integer , dimension(9) :: tzyx
          
    character(len=36) :: ctime

    ctime = tochar(idate)

! total number of output  = ntr + 1 for optical properties file
    noutf = ntr
    if (iaerosol == 1 ) noutf = ntr + 1 

    if ( .not. allocated(ncche) ) then
      allocate(ncche(noutf))
    end if
    ibin = 0
    jbin = 0
!   tracer loop , since we are generating one output per
!   tracer + 1 output for OPT
    do itr = 1, noutf
      if ( itr < noutf ) then  
        ncid = ncche(itr)     
        chevarnam =  chtrname(itr)

        if ( chtrname(itr) == 'DUST' )  then
          ibin = ibin+1
          write( chevarnam(5:6),'(I1)') ibin
        end if
        if ( chtrname(itr) == 'SSLT')  then
          jbin = jbin+1
          write( chevarnam(5:6),'(I1)') jbin
        end if
      else if ( itr == noutf) then
        chevarnam = 'OPT'
      end if 
      title = 'ICTP Regional Climatic model V4  '//chevarnam//' output'
      icherefdate = idate
      icherec = 1

!     call close_chem(ncid,itr)

      write (fterr, '(a3,a)') chevarnam, ' FILE'
      write (fbname,'(a,a,i10)') trim(chevarnam), '.', toint10(idate)
      ofname = trim(dirout)//pthsep//trim(domname)// &
                 '_'//trim(fbname)//'.nc'

      write (aline, *) 'Opening new output file ', trim(ofname)
      call say

#ifdef NETCDF4_HDF5
      istatus = nf90_create(ofname, &
                   ior(ior(nf90_clobber,nf90_hdf5),nf90_classic_model),ncid)
#else
      istatus = nf90_create(ofname, nf90_clobber, ncid)
#endif
      ncche(itr) = ncid

      ! Start Global Attributes
      istatus = nf90_put_att(ncid, nf90_global, 'title', title)
      call check_ok(__FILE__,__LINE__,'Error add title', fterr)
      istatus = nf90_put_att(ncid, nf90_global, 'institution', 'ICTP')
      call check_ok(__FILE__,__LINE__,'Error add institution', fterr)
      istatus = nf90_put_att(ncid, nf90_global, 'source', &
                  'RegCM Model '//'SVN_REV'//' simulation output')
      call check_ok(__FILE__,__LINE__,'Error add source', fterr)
      istatus = nf90_put_att(ncid, nf90_global, 'Conventions', 'CF-1.4')
      call check_ok(__FILE__,__LINE__,'Error add Conventions', fterr)
      call date_and_time(values=tvals)
      write (history,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a)') &
             tvals(1) , '-' , tvals(2) , '-' , tvals(3) , ' ' ,       &
             tvals(5) , ':' , tvals(6) , ':' , tvals(7) ,             &
             ' : Created by RegCM model'
      istatus = nf90_put_att(ncid, nf90_global, 'history', history)
      call check_ok(__FILE__,__LINE__,'Error add history', fterr)
      istatus = nf90_put_att(ncid, nf90_global, 'references', &
               'http://eforge.escience-lab.org/gf/project/regcm')
      call check_ok(__FILE__,__LINE__,'Error add references', fterr)
      istatus = nf90_put_att(ncid, nf90_global, 'experiment', domname)
      call check_ok(__FILE__,__LINE__,'Error add experiment', fterr)
      istatus = nf90_put_att(ncid, nf90_global, 'projection', iproj)
      call check_ok(__FILE__,__LINE__,'Error add projection', fterr)
      if (iproj == 'LAMCON') then
        istatus = nf90_put_att(ncid, nf90_global, &
                        'grid_mapping_name', 'lambert_conformal_conic')
        call check_ok(__FILE__,__LINE__,'Error add grid_mapping_name',fterr)
      else if (iproj == 'POLSTR') then
        istatus = nf90_put_att(ncid, nf90_global, &
                        'grid_mapping_name', 'stereographic')
        call check_ok(__FILE__,__LINE__,'Error add grid_mapping_name',fterr)
      else if (iproj == 'NORMER') then
        istatus = nf90_put_att(ncid, nf90_global, &
                        'grid_mapping_name', 'mercator')
        call check_ok(__FILE__,__LINE__,'Error add grid_mapping_name',fterr)
      else if (iproj == 'ROTMER') then
        istatus = nf90_put_att(ncid, nf90_global, &
                        'grid_mapping_name', 'rotated_latitude_longitude')
        call check_ok(__FILE__,__LINE__,'Error add grid_mapping_name',fterr)
      end if
      istatus = nf90_put_att(ncid, nf90_global,'grid_size_in_meters', ds*d_1000)
      call check_ok(__FILE__,__LINE__,'Error add gridsize', fterr)
      istatus = nf90_put_att(ncid, nf90_global, &
                       'latitude_of_projection_origin', clat)
      call check_ok(__FILE__,__LINE__,'Error add clat', fterr)
      istatus = nf90_put_att(ncid, nf90_global,   &
                       'longitude_of_projection_origin', clon)
      call check_ok(__FILE__,__LINE__,'Error add clon', fterr)
      istatus = nf90_put_att(ncid, nf90_global,   &
                       'longitude_of_central_meridian', clon)
      call check_ok(__FILE__,__LINE__,'Error add gmtllon', fterr)
      if (iproj == 'ROTMER') then
        istatus = nf90_put_att(ncid, nf90_global, &
                       'grid_north_pole_latitude', plat)
        call check_ok(__FILE__,__LINE__,'Error add plat', fterr)
        istatus = nf90_put_att(ncid, nf90_global, &
                       'grid_north_pole_longitude', plon)
        call check_ok(__FILE__,__LINE__,'Error add plon', fterr)
      else if (iproj == 'LAMCON') then
        trlat(1) = real(truelatl)
        trlat(2) = real(truelath)
        istatus = nf90_put_att(ncid, nf90_global, 'standard_parallel', trlat)
        call check_ok(__FILE__,__LINE__,'Error add truelat', fterr)
      else if (iproj == 'NORMER') then
        istatus = nf90_put_att(ncid, nf90_global, 'standard_parallel', clat)
        call check_ok(__FILE__,__LINE__,'Error add truelat', fterr)
      else if (iproj == 'POLSTR') then
        trlat(1) = 1.0
        istatus = nf90_put_att(ncid, nf90_global, &
                   'scale_factor_at_projection_origin', trlat(1:1))
        call check_ok(__FILE__,__LINE__,'Error add scfac', fterr)
      end if
!
!     ADD RUN PARAMETERS
!
!FAB  istatus = nf90_put_att(ncid, nf90_global, 'model_IPCC_scenario', scenario)
!            call check_ok(__FILE__,__LINE__,'Error add scenario', fterr)

!!$                        !!!Je skip ce bloc pour l'instant !!!
!!$    call cdumlbcs
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_boundary_conditions' , trim(cdum))
!!$    call check_ok(__FILE__,__LINE__,'Error add lbcs', fterr)
!!$    call cdumcums
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_cumulous_convection_scheme' , trim(cdum))
!!$    call check_ok(__FILE__,__LINE__,'Error add icup', fterr)
!!$    if (icup == 2 .or. icup == 99 .or. icup == 98) then
!!$      call cdumcumcl
!!$      istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_convective_closure_assumption' , trim(cdum))
!!$      call check_ok(__FILE__,__LINE__,'Error add igcc', fterr)
!!$    end if
!!$    call cdumpbl
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_boundary_layer_scheme' , trim(cdum))
!!$    call check_ok(__FILE__,__LINE__,'Error add ibltyp', fterr)
!!$    call cdummoist
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_moist_physics_scheme' , trim(cdum))
!!$    call check_ok(__FILE__,__LINE__,'Error add ipptls', fterr)
!!$    call cdumocnflx
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_ocean_flux_scheme' , trim(cdum))
!!$    call check_ok(__FILE__,__LINE__,'Error add iocnflx', fterr)
!!$    call cdumpgfs
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_pressure_gradient_force_scheme' , trim(cdum))
!!$    call check_ok(__FILE__,__LINE__,'Error add ipgf', fterr)
!!$    call cdumemiss
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_use_emission_factor' , trim(cdum))
!!$    call check_ok(__FILE__,__LINE__,'Error add iemiss', fterr)
!!$    call cdumlakes
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_use_lake_model' , trim(cdum))
!!$    call check_ok(__FILE__,__LINE__,'Error add lakemod', fterr)
!!$    call cdumchems
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_chemistry' , trim(cdum))
!!$    call check_ok(__FILE__,__LINE__,'Error add ichem', fterr)
!!$    call cdumdcsst
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_diurnal_cycle_sst' , trim(cdum))
!!$    call check_ok(__FILE__,__LINE__,'Error add dcsst', fterr)
!!$    call cdumseaice
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_seaice_effect' , trim(cdum))
!!$    call check_ok(__FILE__,__LINE__,'Error add seaice', fterr)
!!$    call cdumdesseas
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_seasonal_desert_albedo_effect' , trim(cdum))
!!$
!!$    call check_ok(__FILE__,__LINE__,'Error add desseas', fterr)

      istatus = nf90_put_att(ncid, nf90_global,  &
               'model_simulation_initial_start' , tochar(globidate1))
      call check_ok(__FILE__,__LINE__,'Error add globidate1', fterr)

!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_simulation_start' , tochar(idate1))
!!$    call check_ok(__FILE__,__LINE__,'Error add idate1', fterr)
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_simulation_expected_end' , tochar(idate2))
!!$    call check_ok(__FILE__,__LINE__,'Error add idate2', fterr)


    if (ifrest) then
      cdum = 'Yes'
    else
      cdum = 'No'
    end if
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_simulation_is_a_restart' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add ifrest', fterr)


!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_timestep_in_seconds' , dt)
!!$    call check_ok(__FILE__,__LINE__,'Error add dt', fterr)
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_timestep_in_minutes_solar_rad_calc' , dtrad)
!!$    call check_ok(__FILE__,__LINE__,'Error add dtrad', fterr)
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_timestep_in_seconds_bats_calc' , dtsrf)
!!$    call check_ok(__FILE__,__LINE__,'Error add dtsrf', fterr)
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_timestep_in_hours_radiation_calc' , dtabem)
!!$    call check_ok(__FILE__,__LINE__,'Error add dtabem', fterr)
!!$    istatus = nf90_put_att(ncid, nf90_global,  &
!!$            'model_timestep_in_hours_boundary_input' , ibdyfrq)
!!$    call check_ok(__FILE__,__LINE__,'Error add dtabem', fterr)


!End of Global Attributes
!
!         ADD DIMENSIONS
!
      istatus = nf90_def_dim(ncid, 'iy', o_ni, idims(2))
      call check_ok(__FILE__,__LINE__,'Error create dim iy', fterr)
      istatus = nf90_def_dim(ncid, 'jx', o_nj, idims(1))
      call check_ok(__FILE__,__LINE__,'Error create dim jx', fterr)
      istatus = nf90_def_dim(ncid, 'time', nf90_unlimited, idims(3))
      call check_ok(__FILE__,__LINE__,'Error create dim time', fterr)
      istatus = nf90_def_dim(ncid, 'kz', kz, idims(4))
      call check_ok(__FILE__,__LINE__,'Error create dim kz', fterr)

!
!         OUT TYPE DEPENDENT DIMENSIONS
!
    istatus = nf90_def_var(ncid, 'sigma', nf90_float, idims(4), izvar(1))
    call check_ok(__FILE__,__LINE__,'Error add var sigma', fterr)
    istatus = nf90_put_att(ncid, izvar(1), 'standard_name', &
                         'atmosphere_sigma_coordinate')
    call check_ok(__FILE__,__LINE__,'Error add sigma standard_name', fterr)
    istatus = nf90_put_att(ncid, izvar(1), 'long_name', &
                         'Sigma at model layers')
    call check_ok(__FILE__,__LINE__,'Error add sigma long_name', fterr)
    istatus = nf90_put_att(ncid, izvar(1), 'units', '1')
    call check_ok(__FILE__,__LINE__,'Error add sigma units', fterr)
    istatus = nf90_put_att(ncid, izvar(1), 'axis', 'Z')
    call check_ok(__FILE__,__LINE__,'Error add sigma axis', fterr)
    istatus = nf90_put_att(ncid, izvar(1), 'positive', 'down')
    call check_ok(__FILE__,__LINE__,'Error add sigma positive', fterr)
    istatus = nf90_put_att(ncid, izvar(1), 'formula_terms',  &
                         'sigma: sigma ps: ps ptop: ptop')
    call check_ok(__FILE__,__LINE__,'Error add sigma formula_terms', fterr)
   
    istatus = nf90_def_var(ncid, 'ptop', nf90_float, varid=izvar(2))
    call check_ok(__FILE__,__LINE__,'Error add var ptop', fterr)
    istatus = nf90_put_att(ncid, izvar(2), 'standard_name', 'air_pressure')
    call check_ok(__FILE__,__LINE__,'Error add ptop standard_name', fterr)
    istatus = nf90_put_att(ncid, izvar(2), 'long_name', 'Pressure at model top')
    call check_ok(__FILE__,__LINE__,'Error add ptop long_name', fterr)
    istatus = nf90_put_att(ncid, izvar(2), 'units', 'hPa')
    call check_ok(__FILE__,__LINE__,'Error add ptop units', fterr)
    istatus = nf90_def_var(ncid, 'iy', nf90_float, idims(2), ivvar(1))
    call check_ok(__FILE__,__LINE__,'Error add var iy', fterr)
    istatus = nf90_put_att(ncid, ivvar(1), 'standard_name', &
                           'projection_y_coordinate')
    call check_ok(__FILE__,__LINE__,'Error add iy standard_name', fterr)
    istatus = nf90_put_att(ncid, ivvar(1), 'long_name', &
                           'y-coordinate in Cartesian system')
    call check_ok(__FILE__,__LINE__,'Error add iy long_name', fterr)
    istatus = nf90_put_att(ncid, ivvar(1), 'units', 'km')
    call check_ok(__FILE__,__LINE__,'Error add iy units', fterr)
    istatus = nf90_def_var(ncid, 'jx', nf90_float, idims(1), ivvar(2))
    call check_ok(__FILE__,__LINE__,'Error add var jx', fterr)
    istatus = nf90_put_att(ncid, ivvar(2), 'standard_name', &
                         'projection_x_coordinate')
    call check_ok(__FILE__,__LINE__,'Error add jx standard_name', fterr)
    istatus = nf90_put_att(ncid, ivvar(2), 'long_name', &
                         'x-coordinate in Cartesian system')
    call check_ok(__FILE__,__LINE__,'Error add jx long_name', fterr)
    istatus = nf90_put_att(ncid, ivvar(2), 'units', 'km')
    call check_ok(__FILE__,__LINE__,'Error add jx units', fterr)
    istatus = nf90_def_var(ncid, 'xlat', nf90_float, idims(1:2), illtpvar(1))
    call check_ok(__FILE__,__LINE__,'Error add var xlat', fterr)
    istatus = nf90_put_att(ncid, illtpvar(1), 'standard_name', 'latitude')
    call check_ok(__FILE__,__LINE__,'Error add xlat standard_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(1), 'long_name', &
                         'Latitude at cross points')
    call check_ok(__FILE__,__LINE__,'Error add xlat long_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(1), 'units', 'degrees_north')
    call check_ok(__FILE__,__LINE__,'Error add xlat units', fterr)
    istatus = nf90_def_var(ncid, 'xlon', nf90_float, idims(1:2), illtpvar(2))
    call check_ok(__FILE__,__LINE__,'Error add var xlon', fterr)
    istatus = nf90_put_att(ncid, illtpvar(2), 'standard_name', 'longitude')
    call check_ok(__FILE__,__LINE__,'Error add xlon standard_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(2), 'long_name', &
                         'Longitude at cross points')
    call check_ok(__FILE__,__LINE__,'Error add xlon long_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(2), 'units', 'degrees_east')
    call check_ok(__FILE__,__LINE__,'Error add xlon units', fterr)
    istatus = nf90_def_var(ncid, 'topo', nf90_float, idims(1:2), illtpvar(3))
    call check_ok(__FILE__,__LINE__,'Error add var topo', fterr)
    istatus = nf90_put_att(ncid, illtpvar(3), 'standard_name', &
                         'surface_altitude')
    call check_ok(__FILE__,__LINE__,'Error add topo standard_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(3), 'long_name',     &
                         'Domain surface elevation')
    call check_ok(__FILE__,__LINE__,'Error add topo long_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(3), 'units', 'm')
    call check_ok(__FILE__,__LINE__,'Error add topo units', fterr)
    istatus = nf90_put_att(ncid, illtpvar(3), 'coordinates', 'xlat xlon')
    call check_ok(__FILE__,__LINE__,'Error add topo coord', fterr)
    istatus = nf90_put_att(ncid, illtpvar(3), 'grid_mapping', 'rcm_map')
    call check_ok(__FILE__,__LINE__,'Error add topo grid_mapping', fterr)
    istatus = nf90_def_var(ncid, 'mask', nf90_float, idims(1:2), illtpvar(4))
    call check_ok(__FILE__,__LINE__,'Error add var mask', fterr)
    istatus = nf90_put_att(ncid, illtpvar(4), 'standard_name', 'landmask')
    call check_ok(__FILE__,__LINE__,'Error add mask standard_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(4), 'long_name',     &
                         'Domain land/ocean mask')
    call check_ok(__FILE__,__LINE__,'Error add mask long_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(4), 'units', '1')
    call check_ok(__FILE__,__LINE__,'Error add mask units', fterr)
    istatus = nf90_put_att(ncid, illtpvar(4), 'coordinates', 'xlat xlon')
    call check_ok(__FILE__,__LINE__,'Error add mask coord', fterr)
    istatus = nf90_put_att(ncid, illtpvar(4), 'grid_mapping', 'rcm_map')
    call check_ok(__FILE__,__LINE__,'Error add mask grid_mapping', fterr)
    istatus = nf90_def_var(ncid, 'time', nf90_double, idims(3:3), itvar)
    call check_ok(__FILE__,__LINE__,'Error add var time', fterr)
    istatus = nf90_put_att(ncid, itvar, 'standard_name', 'time')
    call check_ok(__FILE__,__LINE__,'Error add time standard_name', fterr)
    istatus = nf90_put_att(ncid, itvar, 'long_name', 'time')
    call check_ok(__FILE__,__LINE__,'Error add time long_name', fterr)
    istatus = nf90_put_att(ncid, itvar, 'calendar', calstr(idate%calendar))
    call check_ok(__FILE__,__LINE__,'Error add time calendar', fterr)
    istatus = nf90_put_att(ncid, itvar, 'units', 'hours since '//ctime)
    call check_ok(__FILE__,__LINE__,'Error add time units', fterr)
   
    istatus = nf90_def_var(ncid, 'ps', nf90_float, idims(1:3), illtpvar(5))
    call check_ok(__FILE__,__LINE__,'Error add var ps', fterr)
    istatus = nf90_put_att(ncid, illtpvar(5), 'standard_name', &
                         'surface_air_pressure')
    call check_ok(__FILE__,__LINE__,'Error add ps standard_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(5), 'long_name', 'Surface pressure')
    call check_ok(__FILE__,__LINE__,'Error add ps long_name', fterr)
    istatus = nf90_put_att(ncid, illtpvar(5), 'units', 'hPa')
    call check_ok(__FILE__,__LINE__,'Error add ps units', fterr)
    istatus = nf90_put_att(ncid, illtpvar(5), 'coordinates', 'xlat xlon')
    call check_ok(__FILE__,__LINE__,'Error add ps coord', fterr)
    istatus = nf90_put_att(ncid, illtpvar(5), 'grid_mapping', 'rcm_map')
    call check_ok(__FILE__,__LINE__,'Error add ps grid_mapping', fterr)
    istatus = nf90_put_att(ncid, illtpvar(5), 'cell_methods', 'time: point')
    call check_ok(__FILE__,__LINE__,'Error add ps cell_methods', fterr)

    tyx = (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/)
    tzyx = (/idims(1),idims(2),idims(4),idims(3),-1,-1,-1,-1,-1/)
         

    if (itr < noutf) then      

            ichevar = -1
            ichevar(1) = itvar
            ichevar(2) = illtpvar(5)

            call ch_addvara(ncid,chevarnam,chevarnam, &
                'atmosphere_mixing_ratio_of_tracer', &
                'Tracers mixing ratios','kg kg-1', &
                tzyx,.false.,ichevar(3))

            call ch_addvara(ncid,chevarnam,'wetdep_ls_flx', &
                'wet_deposition_from_large_scale_precip', &
                'Wet deposition LS','mg/m2/d', &
                tyx,.false.,ichevar(4))

            call ch_addvara(ncid,chevarnam,'wetdep_conv_flx', &
                'wet_deposition_from_convective_precip', &
                'Wet deposition CONV','mg/m2/d', &
                tyx,.false.,ichevar(5))

            call ch_addvara(ncid,chevarnam,'drydep_flx', &
                'dry_deposition', &
                'Dry deposition rate','mg/m2/d', &
                tyx,.false.,ichevar(6))

            call ch_addvara(ncid,chevarnam,'emiss_flx', &
                'surface_emission_rate', &
                'Emission rate','mg/m2/d', &
                tyx,.false.,ichevar(7))

             call ch_addvara(ncid,chevarnam,'drydep_vel', &
                'dry deposition velocity', &
                'dr. dep. vel','m.s-1', &
                tyx,.false.,ichevar(8))

     else if (itr==noutf) then 

              
            ioptvar = -1
            ioptvar(1) = itvar
            ioptvar(2) = illtpvar(4)

            call ch_addvara(ncid,chevarnam,'aext8', &
                'aerosol_optical_depth', &
                'aer mix. aod.','1',tzyx,.false.,ioptvar(3))
            call ch_addvara(ncid,chevarnam,'assa8', &
                'aerosol_single_scattering_albedo', &
                'aer mix. sin. scat. alb','1',tzyx,.false.,ioptvar(4))
            call ch_addvara(ncid,chevarnam,'agfu8', &
                'aerosol_asymmetry_parameter', &
                'aer mix. sin. scat. alb','1',tzyx,.false.,ioptvar(5))

            call ch_addvara(ncid,chevarnam,'acstoarf', &
                'toa_instantaneous_shortwave_radiative_forcing', &
                'TOArad SW forcing av.','W m-2', &
                tyx,.false.,ioptvar(6))
            call ch_addvara(ncid,chevarnam,'acstsrrf', &
                'surface_shortwave_radiative_forcing', &
                'SRFrad SW forcing av.','W m-2', &
                tyx,.false.,ioptvar(7))
            call ch_addvara(ncid,chevarnam,'acstalrf', &
                'toa_longwave_radiative_forcing', &
                'TOArad LW forcing av.','W m-2', &
                tyx,.false.,ioptvar(8))
            call ch_addvara(ncid,chevarnam,'acssrlrf', &
                'surface_longwave_radiative_forcing', &
                'SRFrad LW forcing av.','W m-2', &
                tyx,.false.,ioptvar(9))

     end if 



    istatus = nf90_enddef(ncid)
    call check_ok(__FILE__,__LINE__,'Error End Definitions NetCDF output',fterr)

! write variables which are not time dependant in the file

    istatus = nf90_put_var(ncid, izvar(1), hsigma)
    call check_ok(__FILE__,__LINE__,'Error var sigma write', fterr)
    hptop = real(ptop*d_10)
    istatus = nf90_put_var(ncid, izvar(2), hptop)
    call check_ok(__FILE__,__LINE__,'Error var ptop write', fterr)
 
    
      yiy(1) = -real((dble(o_ni-1)/2.0D0)*ds)
      xjx(1) = -real((dble(o_nj-1)/2.0D0)*ds)
      do i = 2 , o_ni
        yiy(i) = yiy(i-1)+real(ds)
      end do
      do j = 2 , o_nj
        xjx(j) = xjx(j-1)+real(ds)
      end do
      istatus = nf90_put_var(ncid, ivvar(1), yiy(1:o_ni))
      call check_ok(__FILE__,__LINE__,'Error var iy write', fterr)
      istatus = nf90_put_var(ncid, ivvar(2), xjx(1:o_nj))
      call check_ok(__FILE__,__LINE__,'Error var jx write', fterr)
      istatus = nf90_put_var(ncid, illtpvar(1), ioxlat)
      call check_ok(__FILE__,__LINE__,'Error var xlat write', fterr)
      istatus = nf90_put_var(ncid, illtpvar(2), ioxlon)
      call check_ok(__FILE__,__LINE__,'Error var xlon write', fterr)
      istatus = nf90_put_var(ncid, illtpvar(3), iotopo)
      call check_ok(__FILE__,__LINE__,'Error var topo write', fterr)
      istatus = nf90_put_var(ncid, illtpvar(4), iomask)
      call check_ok(__FILE__,__LINE__,'Error var mask write', fterr)
  
      ncche(itr) = ncid

! end of tracer loop       
       end do      

        end subroutine prepare_chem_out
!==========================================================================


!=====================================================================
    subroutine ch_addvara(ncid,ctype,vname,vst,vln,vuni,idims,lmiss,ivar)
      implicit none
      integer , intent(in) :: ncid
      character(3) , intent(in) :: ctype
      character(len=*) , intent(in) :: vname
      character(len=*) , intent(in) :: vst , vln , vuni
      integer , dimension(5) , intent(in) :: idims
      logical , intent(in) :: lmiss
      integer , intent(out) :: ivar
      character(64) :: cdum
      real(sp) , parameter :: fillv = +1E+20
      integer :: i , ndims
       
      ndims = 0
      do i = 1 , 5
        if ( idims(i) > 0 ) ndims = ndims+1
      end do

      cdum = vname
      istatus = nf90_def_var(ncid, cdum, nf90_float, &
                             idims(1:ndims), ivar)
      call check_ok(__FILE__,__LINE__,'Error adding variable '//vname, &
                    ctype//' FILE ERROR')
#ifdef NETCDF4_HDF5
      istatus = nf90_def_var_deflate(ncid, ivar, 1, 1, 9)
      call check_ok(__FILE__,__LINE__, &
           'Error setting deflate on variable '//vname,ctype//' FILE ERROR')
#endif
      cdum = vst
      istatus = nf90_put_att(ncid, ivar, 'standard_name', cdum)
      call check_ok(__FILE__,__LINE__, &
            'Error adding '//vname//' standard_name',ctype//' FILE ERROR')
      cdum = vln
      istatus = nf90_put_att(ncid, ivar, 'long_name', cdum)
      call check_ok(__FILE__,__LINE__,'Error adding '//vname//' long_name', &
                     ctype//' FILE ERROR')
      cdum = vuni
      istatus = nf90_put_att(ncid, ivar, 'units', cdum)
      call check_ok(__FILE__,__LINE__,'Error adding '//vname//' units', &
                     ctype//' FILE ERROR')
      istatus = nf90_put_att(ncid, ivar, 'coordinates','xlat xlon')
      call check_ok(__FILE__,__LINE__,'Error adding '//vname//' coordinates', &
                     ctype//' FILE ERROR')
      istatus = nf90_put_att(ncid, ivar, 'grid_mapping','rcm_map')
      call check_ok(__FILE__,__LINE__,'Error adding '//vname//' grid_mapping', &
                    ctype//' FILE ERROR')
      if (lmiss) then
        istatus = nf90_put_att(ncid, ivar, '_FillValue',fillv)
        call check_ok(__FILE__,__LINE__,'Error adding '//vname//' coordinates',&
                      ctype//' FILE ERROR')
      end if
    end subroutine ch_addvara

!============================================================================

    subroutine writerec_che2(nx,ny,nnx,nny,nz,nt,chia,wdlsc,wdcvc, &
                             ddsfc,cemtrac,drydepv,ext,ssa,asp,    &
                             tarf,ssrf,talwrf,srlwrf,ps,idate)
      implicit none
          
      type(rcm_time_and_date) , intent(in) :: idate
      integer , intent(in) :: nx , ny , nnx , nny , nz , nt 
      real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: chia
      real(dp) , pointer , dimension(:,:) , intent(in) :: ps
      real(dp) , pointer , dimension(:,:,:) , intent(in) :: wdlsc , wdcvc , &
                                                   ddsfc , cemtrac , drydepv
      real(dp) , pointer , dimension(:,:,:) , intent(in) :: ext , ssa , asp 
      real(dp) , pointer , dimension(:,:) , intent(in) :: tarf , ssrf , &
                                                          talwrf , srlwrf

      integer :: n , k, noutf
      integer , dimension(5) :: istart , icount
                    
      real(dp) , dimension(1) :: nctime
      type(rcm_time_interval) :: tdif
      character(len=36) :: ctime
      real(dp) :: cfd2

      noutf = ntr

      if (iaerosol == 1 ) noutf = ntr + 1 
 
      do n = 1 , noutf        
        istart(1) = icherec
        icount(1) = 1
        ctime = tochar(idate)
        tdif = idate-icherefdate
        nctime(1) = tohours(tdif)
        istatus = nf90_put_var(ncche(n), ichevar(1), nctime, &
                               istart(1:1), icount(1:1))
        call check_ok(__FILE__,__LINE__, &
                      'Error writing itime '//ctime, 'CHE FILE ERROR')
        dumio(:,:,1) = real((ps(o_js:o_je,o_is:o_ie)+rpt)*10.0D0)
        istart(3) = icherec
        istart(2) = 1
        istart(1) = 1
        icount(3) = 1
        icount(2) = o_ni
        icount(1) = o_nj
        istatus = nf90_put_var(ncche(n), ichevar(2), &
                               dumio(:,:,1), istart(1:3), icount(1:3))
        call check_ok(__FILE__,__LINE__, &
                      'Error writing ps at '//ctime, 'CHE FILE ERROR')

        if ( n < noutf ) then 
          istart(4) = icherec
          istart(3) = 1
          istart(2) = 1
          istart(1) = 1
          icount(5) = 1
          icount(4) = 1
          icount(3) = o_nz
          icount(2) = o_ni
          icount(1) = o_nj

          !*** tracer concentration
          do k = 1 , nz
            dumio(:,:,k) = real(chia(o_js:o_je,o_is:o_ie,nz-k+1,n) / &
                                ps(o_js:o_je,o_is:o_ie))
          end do
          istatus = nf90_put_var(ncche(n), ichevar(3), &
                                 dumio, istart, icount)
          call check_ok(__FILE__,__LINE__, &
               'Error writing '//chtrname(n)//' at '//ctime,'CHE FILE ERROR')
          
          !*** output 2-D gas chem fields
          istart(4) = 1
          istart(3) = icherec
          istart(2) = 1
          istart(1) = 1
          icount(4) = 1
          icount(3) = 1
          icount(2) = o_ni
          icount(1) = o_nj

          ! accumulated quantities between two output steps are converted
          ! to deposition/emission mean rate (mg /m2/per day)  
          cfd = 24.0D0/chemfrq
          cfd2 = dtche / (chemfrq *3600.0D0)

          !*** wet deposition from large-scale precip
          dumio(:,:,1) = real(wdlsc(o_js:o_je,o_is:o_ie,n)*cfd)
          istatus = nf90_put_var(ncche(n), ichevar(4), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__, &
               'Error writing wet dep LS at '//ctime,'CHE FILE ERROR')

          !*** wet deposition from convective precip
          dumio(:,:,1) = real(wdcvc(o_js:o_je,o_is:o_ie,n)*cfd)
          istatus = nf90_put_var(ncche(n), ichevar(5), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__, &
               'Error writing wet dep CONV at '//ctime,'CHE FILE ERROR')

          !*** dry deposition
          dumio(:,:,1) = real(ddsfc(o_js:o_je,o_is:o_ie,n)*cfd)
          istatus = nf90_put_var(ncche(n), ichevar(6), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__, &
               'Error writing dry dep '//ctime, 'CHE FILE ERROR')

          !*** emission rates
          dumio(:,:,1) = real(cemtrac(o_js:o_je,o_is:o_ie,n)*cfd)
          istatus = nf90_put_var(ncche(n), ichevar(7), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__, &
               'Error writing emission rate '//ctime, 'CHE FILE ERROR')

          !*** dry dep vel 
          dumio(:,:,1) = real(drydepv(o_js:o_je,o_is:o_ie,n)*cfd2)
          istatus = nf90_put_var(ncche(n), ichevar(8), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__, &
               'Error writing dr.dep.vel '//ctime, 'CHE FILE ERROR')

          istatus = nf90_sync(ncche(n))
          call check_ok(__FILE__,__LINE__, &
                        'Error sync at '//ctime, 'CHE FILE ERROR')

        else if ( n == noutf ) then 

          istart(4) = icherec
          istart(3) = 1
          istart(2) = 1
          istart(1) = 1
          icount(5) = 1
          icount(4) = 1
          icount(3) = o_nz
          icount(2) = o_ni
          icount(1) = o_nj

          !*** extinction
          do k = 1 , nz
            dumio(:,:,k) = real(ext(o_js:o_je,o_is:o_ie,nz-k+1))
          end do
          istatus = nf90_put_var(ncche(n), ioptvar(3), &
                                 dumio, istart(1:4), icount(1:4))
          call check_ok(__FILE__,__LINE__,'Error writing EXT at '//ctime,&
                         'CHE FILE ERROR')

          !*** SSAE
          do k = 1 , nz
            dumio(:,:,k) = real(ssa(o_js:o_je,o_is:o_ie,nz-k+1))
          end do
          istatus = nf90_put_var(ncche(n), ioptvar(4), &
                                 dumio, istart(1:4), icount(1:4))
          call check_ok(__FILE__,__LINE__,'Error writing SSAE at '//ctime,&
                        'CHE FILE ERROR')

          !*** ASP
          do k = 1 , nz
            dumio(:,:,k) = real(asp(o_js:o_je,o_is:o_ie,nz-k+1))
          end do
          istatus = nf90_put_var(ncche(n), ioptvar(5), &
                                 dumio, istart(1:4), icount(1:4))
          call check_ok(__FILE__,__LINE__,'Error writing ASP at '//ctime,&
                        'CHE FILE ERROR')

          !*** output 2-D Aerosl Radiative forcing : NB these variables
          ! are already accumulated and averagesd in aerout

          istart(3) = icherec
          istart(2) = 1
          istart(1) = 1
          icount(3) = 1
          icount(2) = o_ni
          icount(1) = o_nj

          dumio(:,:,1) = real(tarf(o_js:o_je,o_is:o_ie))
          istatus = nf90_put_var(ncche(n), ioptvar(6) , &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__,'Error writing aertarf at '//ctime, &
                       'OPT FILE ERROR')
          dumio(:,:,1) = real(ssrf(o_js:o_je,o_is:o_ie))
          istatus = nf90_put_var(ncche(n),  ioptvar(7), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__,'Error writing aersrrf at '//ctime, &
                       'OPT FILE ERROR')
          dumio(:,:,1) = real(talwrf(o_js:o_je,o_is:o_ie))
          istatus = nf90_put_var(ncche(n),ioptvar(8), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__, &
                     'Error writing aertalwrf at '//ctime,'OPT FILE ERROR')
          dumio(:,:,1) = real(srlwrf(o_js:o_je,o_is:o_ie))
          istatus = nf90_put_var(ncche(n), ioptvar(9), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__, &
                     'Error writing aersrlwrf at '//ctime,'OPT FILE ERROR')
          istatus = nf90_sync(ncche(n))
          call check_ok(__FILE__,__LINE__, &
                     'Error sync at '//ctime, 'OPT FILE ERROR')
        end if 
      end do     !main species looop
      icherec = icherec + 1
    end subroutine writerec_che2

!===========================================

    integer function chbc_search(idate)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      type(rcm_time_interval) :: tdif
      character(len=32) :: appdat1, appdat2
      if (idate > chbc_idate(ibcnrec) .or. idate < chbc_idate(1)) then
        chbc_search = -1
      else
        tdif = idate - chbc_idate(1)
        ibcrec = (idnint(tohours(tdif))/ibdyfrq)+1 
        if ( ibcrec < 1 .or. ibcrec > ibcnrec ) then
          appdat1 = tochar(idate)
          write (6,*) 'Record is not found in CHBC file for ',appdat1
          appdat1 = tochar(chbc_idate(1))
          appdat2 = tochar(chbc_idate(ibcnrec))
          write (6,*) 'Range is : ', appdat1, '-', appdat2
          call fatal(__FILE__,__LINE__,'CHBC READ')
        end if
        chbc_search = ibcrec
      end if 
    end function chbc_search

!============================================

    subroutine open_chbc(idate)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      character(10) :: ctime
      integer :: idimid , itvar , i , chkdiff
      real(dp) , dimension(:) , allocatable :: icbc_nctime
      character(64) :: icbc_timeunits , icbc_timecal
      integer :: iyy , jxx , kzz

      call close_chbc
      write (ctime, '(i10)') toint10(idate)
      icbcname = trim(dirglob)//pthsep//trim(domname)//'_CHBC.'//ctime//'.nc'
      istatus = nf90_open(icbcname, nf90_nowrite, ichin)
      call check_ok(__FILE__,__LINE__, &
            'Error Opening ICBC file '//trim(icbcname),'CHBC FILE OPEN')
      ibcrec = 1
      ibcnrec = 0
      istatus = nf90_inq_dimid(ichin, 'iy', idimid)
      call check_ok(__FILE__,__LINE__,'Dimension iy miss', 'ICBC FILE')
      istatus = nf90_inquire_dimension(ichin, idimid, len=iyy)
      call check_ok(__FILE__,__LINE__,'Dimension iy read error','ICBC FILE')
      istatus = nf90_inq_dimid(ichin, 'jx', idimid)
      call check_ok(__FILE__,__LINE__,'Dimension jx miss', 'ICBC FILE')
      istatus = nf90_inquire_dimension(ichin, idimid, len=jxx)
      call check_ok(__FILE__,__LINE__,'Dimension jx read error', 'ICBC FILE')
      istatus = nf90_inq_dimid(ichin, 'kz', idimid)
      call check_ok(__FILE__,__LINE__,'Dimension kz miss', 'ICBC FILE')
      istatus = nf90_inquire_dimension(ichin, idimid, len=kzz)
      call check_ok(__FILE__,__LINE__,'Dimension kz read error', 'ICBC FILE')
      if ( iyy /= iy .or. jxx /= jx .or. kzz /= kz ) then
        write (6,*) 'Error: dims from regcm.in and ICBC file differ.'
        write (aline,*) 'Input namelist : IY=', iy , '  JX=', jx , '  KZ=', kz
        call say
        write (aline,*) 'ICBC file      : IY=', iyy, '  JX=', jxx, '  KZ=', kzz
        call say
        call fatal(__FILE__,__LINE__,'DIMENSION MISMATCH')
      end if
      istatus = nf90_inq_dimid(ichin, 'time', idimid)
      call check_ok(__FILE__,__LINE__,'Dimension time miss', 'ICBC FILE')
      istatus = nf90_inquire_dimension(ichin, idimid, len=ibcnrec)
      call check_ok(__FILE__,__LINE__,'Dimension time read error', 'ICBC FILE')
      if ( ibcnrec < 1 ) then
        write (6,*) 'Time var in ICBC has zero dim.'
        call fatal(__FILE__,__LINE__,'ICBC READ')
      end if
      istatus = nf90_inq_varid(ichin, 'time', itvar)
      call check_ok(__FILE__,__LINE__,'variable time miss', 'ICBC FILE')
      istatus = nf90_get_att(ichin, itvar, 'units', icbc_timeunits)
      call check_ok(__FILE__,__LINE__,'variable time units miss','ICBC FILE')
!!$   istatus = nf90_get_att(ichin, itvar, 'calendar', icbc_timecal)
!!$   call check_ok(__FILE__,__LINE__,'variable time calendar miss','ICBC FILE')
      allocate(icbc_nctime(ibcnrec), stat=istatus)
      if ( istatus /= 0 ) then
        write(6,*) 'Memory allocation error in ICBC for time real values'
        call fatal(__FILE__,__LINE__,'ICBC READ')
      end if
      allocate(chbc_idate(ibcnrec), stat=istatus)
      if ( istatus /= 0 ) then
        write(6,*) 'Memory allocation error in ICBC for time array'
        call fatal(__FILE__,__LINE__,'ICBC READ')
      end if
      istatus = nf90_get_var(ichin, itvar, icbc_nctime)
      call check_ok(__FILE__,__LINE__,'variable time read error', 'ICBC FILE')
      do i = 1 , ibcnrec
        chbc_idate(i) = timeval2date(icbc_nctime(i), &
                                     icbc_timeunits,icbc_timecal)
      end do
      if ( ibcnrec > 1 ) then
        chkdiff = idnint(icbc_nctime(2) - icbc_nctime(1))
        if (chkdiff /= ibdyfrq) then
          write (6,*) 'Time var in ICBC inconsistency.'
          write (6,*) 'Expecting ibdyfrq = ', ibdyfrq
          write (6,*) 'Found     ibdyfrq = ', chkdiff
          call fatal(__FILE__,__LINE__,'ICBC READ')
        end if
      end if
      deallocate(icbc_nctime)
      do i = 1 , n_chbcvar
        istatus = nf90_inq_varid(ichin, trim(chbcname(i)), chbc_ivar(i))
        call check_ok(__FILE__,__LINE__, &
             'variable '//trim(chbcname(i))//' missing','CHBC FILE ERROR')
      end do
    end subroutine open_chbc 

    subroutine read_chbc(chebdio)
      implicit none
      real(dp) , dimension (:,:,:,:), intent(out) :: chebdio 
      integer , dimension(4) :: istart , icount
      real(sp) , dimension(jx,iy,kz) :: xread
      integer :: i , j , k, n

      istart(4) = ibcrec
      istart(3) = 1
      istart(2) = 1
      istart(1) = 1
      icount(4) = 1
      icount(3) = kz
      icount(2) = iy
      icount(1) = jx
      do n = 1 , n_chbcvar
        istatus = nf90_get_var(ichin, chbc_ivar(n), xread, &
                  &            istart, icount)
        call check_ok(__FILE__,__LINE__, &
             'variable '//trim(chbcname(n))//' read error','CHBC FILE ERROR')
        do k = 1 , kz
          do j = 1 , jx
            do i = 1 , iy
              chebdio(j,i,k,n) = xread(j,i,k)
            end do
          end do
        end do
      end do
    end subroutine read_chbc

    subroutine close_chbc
      implicit none
      if ( ichin >= 0 ) then
        istatus = nf90_close(ichin)
        call check_ok(__FILE__,__LINE__, &
              'Error Close CHBC file '//trim(icbcname),'CHBC FILE')
        if ( allocated(chbc_idate) ) deallocate(chbc_idate)
        ichin = -1
      end if
    end subroutine close_chbc

    subroutine check_ok(f,l,m1,mf)
      implicit none
      character(*) , intent(in) :: f, m1 , mf
      integer , intent(in) :: l
      if (istatus /= nf90_noerr) then
        write (6,*) trim(m1)
        write (6,*) nf90_strerror(istatus)
        call fatal(f,l,trim(mf))
      end if
    end subroutine check_ok

end module mod_che_ncio
