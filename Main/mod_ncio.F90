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
module mod_ncio
!
  use mod_runparams
  use mod_cu_interface
  use mod_lm_interface
  use mod_rad_interface , only : iemiss
  use mod_pbl_interface , only : ibltyp
  use mod_mpmessage
  use mod_che_interface
  use mod_memutil
!
  private
!
  public :: ivarname_lookup
  public :: init_mod_ncio , release_mod_ncio
  public :: open_domain , read_domain , read_domain_lake, &
            read_subdomain , read_subdomain_lake,         &
            close_domain
  public :: open_icbc , read_icbc , icbc_search
  public :: prepare_common_out
  public :: writerec_atm , writerec_srf , writerec_sub , &
            writerec_rad , writerec_che , writerec_lak
!
  integer :: idmin , isdmin , ibcin , ncatm , ncsrf , &
             ncsub , ncrad , ncche , nclak
  integer :: istatus
  integer :: ibcrec , ibcnrec
  integer :: iatmrec , isrfrec , isubrec , iradrec , icherec , ilakrec
  integer , dimension(n_atmvar) :: iatmvar
  integer , dimension(n_srfvar) :: isrfvar
  integer , dimension(n_subvar) :: isubvar
  integer , dimension(n_radvar) :: iradvar
  integer , dimension(n_chevar) :: ichevar
  integer , dimension(n_lakvar) :: ilakvar
  character(256) :: dname , sdname , icbcname
  type(rcm_time_and_date) , dimension(:) , allocatable :: icbc_idate
  real(4) , dimension(:) , pointer :: hsigma
  integer , dimension(7) :: icbc_ivar
  logical :: lso4p
  real(8) :: tpd, cfd
  real(8) :: xns2d
  real(4) :: xns2r
  type(rcm_time_and_date) , save :: cordex_refdate

  ! DIM1 is iy ,   DIM2 is jx , DIM3 is time ,       DIM4 is kz
  ! DIM5 is m10 ,  DIM6 is m2 , DIM7 is soil_layer , DIM8 is nv
  ! DIM9 is ntr ,  DIM10 is depth for lake
  integer , dimension(10) :: idims

  integer :: o_is
  integer :: o_ie
  integer :: o_js
  integer :: o_je
  integer :: o_ni
  integer :: o_nj
  integer :: o_isg
  integer :: o_ieg
  integer :: o_jsg
  integer :: o_jeg
  integer :: o_nig
  integer :: o_njg
  integer :: o_nz
  logical :: lwrap , lmaskfill

  real(4) , dimension(:,:) , pointer :: ioxlat
  real(4) , dimension(:,:) , pointer :: ioxlon
  real(4) , dimension(:,:) , pointer :: iotopo
  real(4) , dimension(:,:) , pointer :: iomask
  real(4) , dimension(:,:) , pointer :: iolnds
  real(4) , dimension(:,:) , pointer :: ioxlat_s
  real(4) , dimension(:,:) , pointer :: ioxlon_s
  real(4) , dimension(:,:) , pointer :: iotopo_s
  real(4) , dimension(:,:) , pointer :: iomask_s
  real(4) , dimension(:,:) , pointer :: subio
  real(4) , dimension(:,:,:) , pointer :: dumio
  real(4) , dimension(:,:) , pointer :: sp2d
  real(4) , dimension(:,:) , pointer :: sp2d1
  real(4) , dimension(:,:,:) , pointer :: atmsrfmask
  real(4) , dimension(:,:) , pointer :: atmsrfsum

  integer , dimension(numbat) :: lak_fbats

  character(128) :: cdum

  data lso4p   /.false./
  data lmaskfill /.false./
  data idmin   /-1/
  data isdmin  /-1/
  data ibcin   /-1/
  data ibcrec  / 1/
  data ibcnrec / 0/
  data ncatm   /-1/
  data iatmrec / 1/
  data ncsrf   /-1/
  data isrfrec / 1/
  data ncsub   /-1/
  data isubrec / 1/
  data ncrad   /-1/
  data iradrec / 1/
  data ncche   /-1/
  data icherec / 1/
  data nclak   /-1/
  data ilakrec / 1/

  data lak_fbats / 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, &
                   1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 /

contains

  subroutine cdumlbcs
    implicit none
    select case (iboudy)
      case(0)
       write (cdum,'(a)') 'Fixed'
      case(1)
       write (cdum,'(a)') 'Relaxation, linear technique'
      case(2)
       write (cdum,'(a)') 'Time-dependent'
      case(3)
       write (cdum,'(a)') 'Time and inflow/outflow dependent'
      case(4)
       write (cdum,'(a)') 'Sponge (Perkey & Kreitzberg, MWR 1976)'
      case(5)
       write (cdum,'(a)') 'Relaxation, exponential technique'
      case default 
       write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumlbcs

  subroutine cdumcums
    implicit none
    select case (icup)
      case(1)
       write (cdum,'(a)') 'Kuo'
      case(2)
       write (cdum,'(a)') 'Grell'
      case(3)
       write (cdum,'(a)') 'Betts-Miller (1986)'
      case(4)
       write (cdum,'(a)') 'Emanuel (1991)'
      case(5)
        write (cdum,'(a)') 'Tiedtke (1986)'
      case(98)
        write (cdum,'(a)') 'Grell over ocean, Emanuel (1991) over land'
      case(99)
        write (cdum,'(a)') 'Emanuel (1991) over ocean, Grell over land'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumcums

  subroutine cdumcumcl
    implicit none
    select case (igcc)
      case(1)
        write (cdum,'(a)') 'Arakawa & Schubert (1974)'
      case(2)
        write (cdum,'(a)') 'Fritsch & Chappell (1980)'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumcumcl

  subroutine cdumpbl
    implicit none
    select case (ibltyp)
      case(0)
        write (cdum,'(a)') 'Frictionless'
      case(1)
        write (cdum,'(a)') 'Holtslag PBL (Holtslag, 1990)'
      case(2)
        write (cdum,'(a)') 'UW PBL (Bretherton and McCaa, 2004)'
      case(99)
        write (cdum,'(a)') 'Holtslag PBL, with UW in diag. mode'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumpbl

  subroutine cdummoist
    implicit none
    select case (ipptls)
      case(1)
        write (cdum,'(a)') 'Explicit moisture (SUBEX; Pal et al 2000)'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdummoist

  subroutine cdumocnflx
    implicit none
    select case (iocnflx)
      case(1)
        write (cdum,'(a)') 'Use BATS1e Monin-Obukhov'
      case(2)
        write (cdum,'(a)') 'Zeng et al (1998)'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumocnflx

  subroutine cdumpgfs
    implicit none
    select case (ipgf)
      case(0)
        write (cdum,'(a)') 'Use full fields'
      case(1)
        write (cdum,'(a)') 'Hydrostatic deduction with perturbation temperature'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumpgfs

  subroutine cdumemiss
    implicit none
    select case (iemiss)
      case(0)
        write (cdum,'(a)') 'No'
      case(1)
        write (cdum,'(a)') 'Yes'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumemiss

  subroutine cdumlakes
    implicit none
    select case (lakemod)
      case(0)
        write (cdum,'(a)') 'No'
      case(1)
        write (cdum,'(a)') 'Yes'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumlakes

  subroutine cdumchems
    implicit none
    select case (ichem)
      case(0)
        write (cdum,'(a)') 'Not active'
      case(1)
        write (cdum,'(a)') 'Active'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumchems

  subroutine cdumdcsst
    implicit none
    select case (idcsst)
      case(0)
        write (cdum,'(a)') 'Not active'
      case(1)
        write (cdum,'(a)') 'Active'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumdcsst

  subroutine cdumseaice
    implicit none
    select case (iseaice)
      case(0)
        write (cdum,'(a)') 'Not active'
      case(1)
        write (cdum,'(a)') 'Active'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumseaice

  subroutine cdumdesseas
    implicit none
    select case (idesseas)
      case(0)
        write (cdum,'(a)') 'Not active'
      case(1)
        write (cdum,'(a)') 'Active'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumdesseas

  function ivarname_lookup(ctype,sname)
    implicit none
    integer :: ivarname_lookup
    character(3) , intent(in) :: ctype
    character(len=*) , intent(in) :: sname
    integer :: i

    ivarname_lookup = -1

    if (ctype == 'ATM') then
      do i = 1 , n_atmvar
        if (sname == atm_variables(i)%vname) then
          ivarname_lookup = i
          exit
        end if
      end do
    else if (ctype == 'SRF') then
      do i = 1 , n_srfvar
        if (sname == srf_variables(i)%vname) then
          ivarname_lookup = i
          exit
        end if
      end do
    else if (ctype == 'SUB') then
      do i = 1 , n_subvar
        if (sname == sub_variables(i)%vname) then
          ivarname_lookup = i
          exit
        end if
      end do
    else if (ctype == 'RAD') then
      do i = 1 , n_radvar
        if (sname == rad_variables(i)%vname) then
          ivarname_lookup = i
          exit
        end if
      end do
    else if (ctype == 'CHE') then
      do i = 1 , n_chevar
        if (sname == che_variables(i)%vname) then
          ivarname_lookup = i
          exit
        end if
      end do
    else if (ctype == 'LAK') then
      do i = 1 , n_lakvar
        if (sname == lak_variables(i)%vname) then
          ivarname_lookup = i
          exit
        end if
      end do
    endif
  end function ivarname_lookup

  subroutine init_mod_ncio(lband)
    implicit none
    logical , intent(in) :: lband
    character(3) :: sbstring
    dname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
    write (sbstring,'(i0.3)') nsg
    sdname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN'//sbstring//'.nc'
    icbcname = trim(dirglob)//pthsep//trim(domname)//'_ICBC.'//'YYYYMMDDHH.nc'

    xns2r = 1.0/real(nnsg)
    xns2d = 1.0D0/dble(nnsg)
    cordex_refdate = 1949120100
    call setcal(cordex_refdate,ical)

    if (lband) then
      o_is = 2
      o_ie = iy-1
      o_js = 1
      o_je = jx
      o_ni = iy-2
      o_nj = jx
      o_isg = nsg+1
      o_ieg = iysg-nsg
      o_jsg = 1
      o_jeg = jxsg
      o_nig = iym2sg
      o_njg = jxsg
      o_nz = kz
      lwrap = .true.
    else
      o_is = 2
      o_ie = iy-1
      o_js = 2
      o_je = jx-1
      o_ni = iy-2
      o_nj = jx-2
      o_isg = nsg+1
      o_ieg = iysg-nsg
      o_jsg = nsg+1
      o_jeg = jxsg-nsg
      o_nig = iym2sg
      o_njg = jxm2sg
      o_nz = kz
      lwrap = .false.
    end if
    call getmem1d(hsigma,1,o_nz,'ncio:hsigma')
    call getmem2d(ioxlat,1,o_nj,1,o_ni,'ncio:ioxlat')
    call getmem2d(ioxlon,1,o_nj,1,o_ni,'ncio:ioxlon')
    call getmem2d(iotopo,1,o_nj,1,o_ni,'ncio:iotopo')
    call getmem2d(iomask,1,o_nj,1,o_ni,'ncio:iomask')
    call getmem2d(iolnds,1,o_nj,1,o_ni,'ncio:iolnds')
    call getmem3d(dumio,1,o_nj,1,o_ni,1,o_nz,'ncio:dumio')
    call getmem2d(sp2d,1,jx,1,iy,'ncio:sp2d')
    call getmem3d(atmsrfmask,1,nnsg,1,o_nj,1,o_ni,'ncio:atmsrfmask')
    call getmem2d(atmsrfsum,1,o_nj,1,o_ni,'ncio:atmsrfsum')
    if (nsg > 1) then
      call getmem2d(ioxlat_s,1,o_njg,1,o_nig,'ncio:ioxlat_s')
      call getmem2d(ioxlon_s,1,o_njg,1,o_nig,'ncio:ioxlon_s')
      call getmem2d(iotopo_s,1,o_njg,1,o_nig,'ncio:iotopo_s')
      call getmem2d(iomask_s,1,o_njg,1,o_nig,'ncio:iomask_s')
      call getmem2d(subio,1,o_njg,1,o_nig,'ncio:subio')
      call getmem2d(sp2d1,1,jxsg,1,iysg,'ncio:sp2d1')
    end if
  end subroutine init_mod_ncio

  subroutine open_domain(dx , sigma)
    use netcdf
    implicit none

    real(8) , intent(out) :: dx
    real(8) , dimension(kzp1) :: sigma

    integer :: ivarid , idimid
    integer :: iyy , jxx , kzz , k
    character(6) :: proj
    real(4) :: dsx , iclat , iclon , ptsp
    real(4) , dimension(kzp1) :: rsdum

    write (aline,*) 'open_domain: READING HEADER FILE:', dname
    call say
    istatus = nf90_open(dname, nf90_nowrite, idmin)
    call check_ok(__FILE__,__LINE__,'Error Opening Domain file '//trim(dname), &
                  'DOMAIN FILE')
    if ( nsg > 1 ) then
      write (aline,*) 'READING HEADER SUBDOMAIN FILE:', sdname
      call say
      istatus = nf90_open(sdname, nf90_nowrite, isdmin)
      call check_ok(__FILE__,__LINE__, &
           'Error Opening SubDomain file '//trim(sdname), 'SUBDOM FILE')
    end if
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

  subroutine read_domain(ht,lnd,xlat,xlon,xmap,dmap,f)
    use netcdf
    implicit none

    real(8) , dimension(iy,jx) , intent(out) :: ht
    real(8) , dimension(iy,jx) , intent(out) :: lnd
    real(8) , dimension(iy,jx) , intent(out) :: xlat
    real(8) , dimension(iy,jx) , intent(out) :: xlon
    real(8) , dimension(iy,jx) , intent(out) :: xmap
    real(8) , dimension(iy,jx) , intent(out) :: dmap
    real(8) , dimension(iy,jx) , intent(out) :: f

    integer :: ivarid

    if (idmin < 0) then
      write (6,*) 'Error : Domain file not in open state'
      call fatal(__FILE__,__LINE__, 'DOMAIN FILE')
    end if

    istatus = nf90_inq_varid(idmin, 'topo', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable topo miss', 'DOMAIN FILE')
    istatus = nf90_get_var(idmin, ivarid, sp2d)
    call check_ok(__FILE__,__LINE__,'Variable topo read error', 'DOMAIN FILE')
    ht = dble(transpose(sp2d))
    iotopo = sp2d(o_js:o_je,o_is:o_ie)
    istatus = nf90_inq_varid(idmin, 'landuse', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable landuse miss','DOMAIN FILE')
    istatus = nf90_get_var(idmin, ivarid, sp2d)
    call check_ok(__FILE__,__LINE__,'Variable landuse read error','DOMAIN FILE')
    lnd = dble(transpose(sp2d))
    iolnds = sp2d(o_js:o_je,o_is:o_ie)
    istatus = nf90_inq_varid(idmin, 'xlat', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable xlat miss', 'DOMAIN FILE')
    istatus = nf90_get_var(idmin, ivarid, sp2d)
    call check_ok(__FILE__,__LINE__,'Variable xlat read error', 'DOMAIN FILE')
    xlat = dble(transpose(sp2d))
    ioxlat = sp2d(o_js:o_je,o_is:o_ie)
    istatus = nf90_inq_varid(idmin, 'xlon', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable xlon miss', 'DOMAIN FILE')
    istatus = nf90_get_var(idmin, ivarid, sp2d)
    call check_ok(__FILE__,__LINE__,'Variable xlon read error', 'DOMAIN FILE')
    xlon = dble(transpose(sp2d))
    ioxlon = sp2d(o_js:o_je,o_is:o_ie)
    istatus = nf90_inq_varid(idmin, 'xmap', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable xmap miss', 'DOMAIN FILE')
    istatus = nf90_get_var(idmin, ivarid, sp2d)
    call check_ok(__FILE__,__LINE__,'Variable xmap read error','DOMAIN FILE')
    xmap = dble(transpose(sp2d))
    istatus = nf90_inq_varid(idmin, 'dmap', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable dmap miss','DOMAIN FILE')
    istatus = nf90_get_var(idmin, ivarid, sp2d)
    call check_ok(__FILE__,__LINE__,'Variable dmap read error', 'DOMAIN FILE')
    dmap = dble(transpose(sp2d))
    istatus = nf90_inq_varid(idmin, 'coriol', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable coriol miss', 'DOMAIN FILE')
    istatus = nf90_get_var(idmin, ivarid, sp2d)
    call check_ok(__FILE__,__LINE__,'Variable coriol read error','DOMAIN FILE')
    f = dble(transpose(sp2d))
    istatus = nf90_inq_varid(idmin, 'mask', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable mask miss', 'DOMAIN FILE')
    istatus = nf90_get_var(idmin, ivarid, sp2d)
    call check_ok(__FILE__,__LINE__,'Variable mask read error','DOMAIN FILE')
    iomask = sp2d(o_js:o_je,o_is:o_ie)
  end subroutine read_domain

  subroutine read_domain_lake(hlake)
    use netcdf
    implicit none

    real(8) , dimension(iy,jx) , intent(out) :: hlake

    integer :: ivarid
    real(4) , dimension(jx,iy) :: sp2d

    if (idmin < 0) then
      write (6,*) 'Error : Domain file not in open state'
      call fatal(__FILE__,__LINE__,'DOMAIN FILE')
    end if

    istatus = nf90_inq_varid(idmin, 'dhlake', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable dhlake miss', 'DOMAIN FILE')
    istatus = nf90_get_var(idmin, ivarid, sp2d)
    call check_ok(__FILE__,__LINE__,'Variable dhlake read error','DOMAIN FILE')
    hlake = dble(transpose(sp2d))
  end subroutine read_domain_lake

  subroutine read_subdomain(ht1,lnd1,xlat1,xlon1)
    use netcdf
    implicit none

    real(8) , dimension(nnsg,iy,jx) , intent(out) :: ht1
    real(8) , dimension(nnsg,iy,jx) , intent(out) :: lnd1
    real(8) , dimension(nnsg,iy,jx) , intent(out) :: xlat1
    real(8) , dimension(nnsg,iy,jx) , intent(out) :: xlon1

    integer :: ivarid
    integer :: i , j , n , ii , jj
    
    if (isdmin < 0) then
      write (6,*) 'Error : Subdom file not in open state'
      call fatal(__FILE__,__LINE__,'SUBDOMAIN FILE')
    end if

    istatus = nf90_inq_varid(isdmin, 'topo', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable topo miss', 'SUBDOMAIN FILE')
    istatus = nf90_get_var(isdmin, ivarid, sp2d1)
    call check_ok(__FILE__,__LINE__,'Variable topo read error','SUBDOMAIN FILE')
    do j = 1 , jxsg
      do i = 1 , iysg
        jj = mod(j,nsg)
        if ( jj == 0 ) jj = nsg
        ii = mod(i,nsg)
        if ( ii == 0 ) ii = nsg
        n = (jj-1)*nsg + ii
        jj = (j+nsg-1)/nsg
        ii = (i+nsg-1)/nsg
        ht1(n,ii,jj) = dble(sp2d1(j,i))*egrav
      end do
    end do
    iotopo_s = sp2d1(o_jsg:o_jeg,o_isg:o_ieg)
    istatus = nf90_inq_varid(isdmin, 'landuse', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable landuse miss','SUBDOMAIN FILE')
    istatus = nf90_get_var(isdmin, ivarid, sp2d1)
    call check_ok(__FILE__,__LINE__,'Variable landuse read error', &
                  'SUBDOMAIN FILE')
    do j = 1 , jxsg
      do i = 1 , iysg
        jj = mod(j,nsg)
        if ( jj == 0 ) jj = nsg
        ii = mod(i,nsg)
        if ( ii == 0 ) ii = nsg
        n = (jj-1)*nsg + ii
        jj = (j+nsg-1)/nsg
        ii = (i+nsg-1)/nsg
        lnd1(n,ii,jj) = dble(sp2d1(j,i))
      end do
    end do
    istatus = nf90_inq_varid(isdmin, 'xlat', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable xlat miss','SUBDOMAIN FILE')
    istatus = nf90_get_var(isdmin, ivarid, sp2d1)
    call check_ok(__FILE__,__LINE__,'Variable xlat read error','SUBDOMAIN FILE')
    do j = 1 , jxsg
      do i = 1 , iysg
        jj = mod(j,nsg)
        if ( jj == 0 ) jj = nsg
        ii = mod(i,nsg)
        if ( ii == 0 ) ii = nsg
        n = (jj-1)*nsg + ii
        jj = (j+nsg-1)/nsg
        ii = (i+nsg-1)/nsg
        xlat1(n,ii,jj) = dble(sp2d1(j,i))
      end do
    end do
    ioxlat_s = sp2d1(o_jsg:o_jeg,o_isg:o_ieg)
    istatus = nf90_inq_varid(isdmin, 'xlon', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable xlon miss','SUBDOMAIN FILE')
    istatus = nf90_get_var(isdmin, ivarid, sp2d1)
    call check_ok(__FILE__,__LINE__,'Variable xlon read error','SUBDOMAIN FILE')
    do j = 1 , jxsg
      do i = 1 , iysg
        jj = mod(j,nsg)
        if ( jj == 0 ) jj = nsg
        ii = mod(i,nsg)
        if ( ii == 0 ) ii = nsg
        n = (jj-1)*nsg + ii
        jj = (j+nsg-1)/nsg
        ii = (i+nsg-1)/nsg
        xlon1(n,ii,jj) = dble(sp2d1(j,i))
      end do
    end do
    ioxlon_s = sp2d1(o_jsg:o_jeg,o_isg:o_ieg)
    istatus = nf90_inq_varid(isdmin, 'mask', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable mask miss','SUBDOMAIN FILE')
    istatus = nf90_get_var(isdmin, ivarid, sp2d1)
    call check_ok(__FILE__,__LINE__,'Variable mask read error','SUBDOMAIN FILE')
    iomask_s = sp2d1(o_jsg:o_jeg,o_isg:o_ieg)
  end subroutine read_subdomain

  subroutine read_subdomain_lake(hlake1)
    use netcdf
    implicit none

    real(8) , dimension(nnsg,iy,jx) , intent(out) :: hlake1

    integer :: ivarid
    integer :: i , j , n , ii , jj
    real(4) , dimension(jxsg,iysg) :: sp2d1
    
    if (isdmin < 0) then
      write (6,*) 'Error : Subdom file not in open state'
      call fatal(__FILE__,__LINE__, 'SUBDOMAIN FILE')
    end if

    istatus = nf90_inq_varid(isdmin, 'dhlake', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable dhlake miss','SUBDOMAIN FILE')
    istatus = nf90_get_var(isdmin, ivarid, sp2d1)
    call check_ok(__FILE__,__LINE__,'Variable dhlake read error', &
                  'SUBDOMAIN FILE')
    do j = 1 , jxsg
      do i = 1 , iysg
        jj = mod(j,nsg)
        if ( jj == 0 ) jj = nsg
        ii = mod(i,nsg)
        if ( ii == 0 ) ii = nsg
        n = (jj-1)*nsg + ii
        jj = (j+nsg-1)/nsg
        ii = (i+nsg-1)/nsg
        hlake1(n,ii,jj) = dble(sp2d1(j,i))
      end do
    end do
  end subroutine read_subdomain_lake

  subroutine close_domain
    use netcdf
    implicit none

    if (idmin >= 0) then
      istatus = nf90_close(idmin)
      call check_ok(__FILE__,__LINE__,'Domain file close error','DOMAIN FILE')
      idmin = -1
    end if
    if ( nsg>1 .and. isdmin >=0 ) then
      istatus = nf90_close(isdmin)
      call check_ok(__FILE__,__LINE__,'SubDomain file close error', &
                   'SUBDOMAIN FILE')
      isdmin = -1
    end if

  end subroutine close_domain

  integer function icbc_search(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    type(rcm_time_interval) :: tdif
    character(len=32) :: appdat1, appdat2
    if (idate > icbc_idate(ibcnrec) .or. idate < icbc_idate(1)) then
      icbc_search = -1
    else
      tdif = idate-icbc_idate(1)
      ibcrec = (idnint(tohours(tdif))/ibdyfrq)+1
      if ( ibcrec < 1 .or. ibcrec > ibcnrec ) then
        appdat1 = tochar(idate)
        write (6,*) 'Record is not found in ICBC file for ',appdat1
        appdat1 = tochar(icbc_idate(1))
        appdat2 = tochar(icbc_idate(ibcnrec))
        write (6,*) 'Range is : ', appdat1, '-', appdat2
        call fatal(__FILE__,__LINE__,'ICBC READ')
      end if
      icbc_search = ibcrec
    end if 
  end function icbc_search

  subroutine open_icbc(idate)
    use netcdf
    type(rcm_time_and_date) , intent(in) :: idate
    character(10) :: ctime
    integer :: idimid , itvar , i , chkdiff
    real(8) , dimension(:) , allocatable :: icbc_nctime
    character(64) :: icbc_timeunits , icbc_timecal
    integer :: iyy , jxx , kzz

    call close_icbc
    write (ctime, '(i10)') toint10(idate)
    icbcname = trim(dirglob)//pthsep//trim(domname)//'_ICBC.'//ctime//'.nc'
    istatus = nf90_open(icbcname, nf90_nowrite, ibcin)
    call check_ok(__FILE__,__LINE__, &
        'Error Opening ICBC file '//trim(icbcname),'ICBC FILE OPEN')
    ibcrec = 1
    ibcnrec = 0
    istatus = nf90_inq_dimid(ibcin, 'iy', idimid)
    call check_ok(__FILE__,__LINE__,'Dimension iy miss', 'ICBC FILE')
    istatus = nf90_inquire_dimension(ibcin, idimid, len=iyy)
    call check_ok(__FILE__,__LINE__,'Dimension iy read error','ICBC FILE')
    istatus = nf90_inq_dimid(ibcin, 'jx', idimid)
    call check_ok(__FILE__,__LINE__,'Dimension jx miss', 'ICBC FILE')
    istatus = nf90_inquire_dimension(ibcin, idimid, len=jxx)
    call check_ok(__FILE__,__LINE__,'Dimension jx read error', 'ICBC FILE')
    istatus = nf90_inq_dimid(ibcin, 'kz', idimid)
    call check_ok(__FILE__,__LINE__,'Dimension kz miss', 'ICBC FILE')
    istatus = nf90_inquire_dimension(ibcin, idimid, len=kzz)
    call check_ok(__FILE__,__LINE__,'Dimension kz read error', 'ICBC FILE')
    if ( iyy /= iy .or. jxx /= jx .or. kzz /= kz ) then
      write (6,*) 'Error: dims from regcm.in and ICBC file differ.'
      write (aline,*) 'Input namelist : IY=', iy , '  JX=',  jx , '  KZ=', kz
      call say
      write (aline,*) 'ICBC file      : IY=', iyy, '  JX=',  jxx, '  KZ=', kzz
      call say
      call fatal(__FILE__,__LINE__,'DIMENSION MISMATCH')
    end if
    istatus = nf90_inq_dimid(ibcin, 'time', idimid)
    call check_ok(__FILE__,__LINE__,'Dimension time miss', 'ICBC FILE')
    istatus = nf90_inquire_dimension(ibcin, idimid, len=ibcnrec)
    call check_ok(__FILE__,__LINE__,'Dimension time read error', 'ICBC FILE')
    if ( ibcnrec < 1 ) then
      write (6,*) 'Time var in ICBC has zero dim.'
      call fatal(__FILE__,__LINE__,'ICBC READ')
    end if
    istatus = nf90_inq_varid(ibcin, 'time', itvar)
    call check_ok(__FILE__,__LINE__,'variable time miss', 'ICBC FILE')
    istatus = nf90_get_att(ibcin, itvar, 'units', icbc_timeunits)
    call check_ok(__FILE__,__LINE__,'variable time units miss','ICBC FILE')
    istatus = nf90_get_att(ibcin, itvar, 'calendar', icbc_timecal)
    call check_ok(__FILE__,__LINE__,'variable time calendar miss','ICBC FILE')
    allocate(icbc_nctime(ibcnrec), stat=istatus)
    if ( istatus /= 0 ) then
      write(6,*) 'Memory allocation error in ICBC for time real values'
      call fatal(__FILE__,__LINE__,'ICBC READ')
    end if
    allocate(icbc_idate(ibcnrec), stat=istatus)
    if ( istatus /= 0 ) then
      write(6,*) 'Memory allocation error in ICBC for time array'
      call fatal(__FILE__,__LINE__,'ICBC READ')
    end if
    istatus = nf90_get_var(ibcin, itvar, icbc_nctime)
    call check_ok(__FILE__,__LINE__,'variable time read error', 'ICBC FILE')
    do i = 1 , ibcnrec
      icbc_idate(i) = timeval2date(icbc_nctime(i), icbc_timeunits, icbc_timecal)
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
    istatus = nf90_inq_varid(ibcin, 'ps', icbc_ivar(1))
    call check_ok(__FILE__,__LINE__,'variable ps miss', 'ICBC FILE')
    istatus = nf90_inq_varid(ibcin, 'ts', icbc_ivar(2))
    call check_ok(__FILE__,__LINE__,'variable ts miss', 'ICBC FILE')
    istatus = nf90_inq_varid(ibcin, 'u', icbc_ivar(3))
    call check_ok(__FILE__,__LINE__,'variable u miss', 'ICBC FILE')
    istatus = nf90_inq_varid(ibcin, 'v', icbc_ivar(4))
    call check_ok(__FILE__,__LINE__,'variable v miss', 'ICBC FILE')
    istatus = nf90_inq_varid(ibcin, 't', icbc_ivar(5))
    call check_ok(__FILE__,__LINE__,'variable t miss', 'ICBC FILE')
    istatus = nf90_inq_varid(ibcin, 'qv', icbc_ivar(6))
    call check_ok(__FILE__,__LINE__,'variable qv miss', 'ICBC FILE')
    istatus = nf90_inq_varid(ibcin, 'so4', icbc_ivar(7))
    if ( istatus == nf90_noerr) then
      lso4p = .true.
    end if
  end subroutine open_icbc

  subroutine read_icbc(ps,ts,u,v,t,qv,so4)
    use netcdf
    implicit none
    real(8) , dimension(iy,kz,jx) , intent(out) :: u
    real(8) , dimension(iy,kz,jx) , intent(out) :: v
    real(8) , dimension(iy,kz,jx) , intent(out) :: t
    real(8) , dimension(iy,kz,jx) , intent(out) :: qv
    real(8) , dimension(iy,kz,jx) , intent(out) :: so4
    real(8) , dimension(iy,jx) , intent(out) :: ps
    real(8) , dimension(iy,jx) , intent(out) :: ts

    integer , dimension(4) :: istart , icount
    real(4) , dimension(jx,iy,kz) :: xread
    integer :: i , j , k

    istart(3) = ibcrec
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = iy
    icount(1) = jx
    istatus = nf90_get_var(ibcin, icbc_ivar(1), xread(:,:,1), & 
                           istart(1:3), icount(1:3))
    call check_ok(__FILE__,__LINE__,'variable ps read error', 'ICBC FILE')
    ps = dble(transpose(xread(:,:,1)))
    istatus = nf90_get_var(ibcin, icbc_ivar(2), xread(:,:,1), & 
                           istart(1:3), icount(1:3))
    call check_ok(__FILE__,__LINE__,'variable ts read error', 'ICBC FILE')
    ts = dble(transpose(xread(:,:,1)))
    istart(4) = ibcrec
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = kz
    icount(2) = iy
    icount(1) = jx
    istatus = nf90_get_var(ibcin, icbc_ivar(3), xread, istart, icount)
    call check_ok(__FILE__,__LINE__,'variable u read error', 'ICBC FILE')
    do k = 1 , kz
      do j = 1 , jx
        do i = 1 , iy
          u(i,k,j) = dble(xread(j,i,k))
        end do
      end do
    end do
    istatus = nf90_get_var(ibcin, icbc_ivar(4), xread, istart, icount)
    call check_ok(__FILE__,__LINE__,'variable v read error', 'ICBC FILE')
    do k = 1 , kz
      do j = 1 , jx
        do i = 1 , iy
          v(i,k,j) = dble(xread(j,i,k))
        end do
      end do
    end do
    istatus = nf90_get_var(ibcin, icbc_ivar(5), xread, istart, icount)
    call check_ok(__FILE__,__LINE__,'variable t read error', 'ICBC FILE')
    do k = 1 , kz
      do j = 1 , jx
        do i = 1 , iy
          t(i,k,j) = dble(xread(j,i,k))
        end do
      end do
    end do
    istatus = nf90_get_var(ibcin, icbc_ivar(6), xread, istart, icount)
    call check_ok(__FILE__,__LINE__,'variable qv read error', 'ICBC FILE')
    do k = 1 , kz
      do j = 1 , jx
        do i = 1 , iy
          qv(i,k,j) = dble(xread(j,i,k))
        end do
      end do
    end do
    if (lso4p) then
      istatus = nf90_get_var(ibcin, icbc_ivar(7), xread, istart, icount)
      call check_ok(__FILE__,__LINE__,'variable so4 read error', 'ICBC FILE')
      do k = 1 , kz
        do j = 1 , jx
          do i = 1 , iy
            so4(i,k,j) = dble(xread(j,i,k))
          end do
        end do
      end do
    else
      so4 = d_zero
    end if
  end subroutine read_icbc

  subroutine close_icbc
    use netcdf
    implicit none
    if (ibcin >= 0) then
      istatus = nf90_close(ibcin)
      call check_ok(__FILE__,__LINE__, &
           'Error Close ICBC file '//trim(icbcname),'ICBC FILE')
      if ( allocated(icbc_idate) ) deallocate(icbc_idate)
      ibcin = -1
    end if
  end subroutine close_icbc

  subroutine close_common(ncid, ctype)
    use netcdf
    implicit none
    integer , intent(inout) :: ncid
    character(3) , intent(in) :: ctype
    if (ncid >= 0) then
      istatus = nf90_close(ncid)
      call check_ok(__FILE__,__LINE__, &
                    'Error Close '//ctype//' file', ctype//' FILE')
      ncid = -1
    end if
  end subroutine close_common

  subroutine prepare_common_out(idate,ctype)
    use netcdf
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    character(3) , intent(in) :: ctype
    character(64) :: title
    character(32) :: fbname , ctime
    character(64) :: cmethodmax , cmethodmin
    character(16) :: fterr
    character(256) :: ofname , history
    integer , dimension(8) :: tvals
    real(4) :: hptop , rdum1
    real(4) , dimension(2) :: trlat , rdum2
    real(4) , dimension(iysg) :: yiy
    real(4) , dimension(jxsg) :: xjx
    real(4) , dimension(ndpmax) :: depth
    integer :: ncid
    integer , dimension(3) :: izvar
    integer , dimension(2) :: ivvar
    integer , dimension(4) :: isrvvar
    integer , dimension(5) :: illtpvar
    integer :: itvar , imapvar , i , j , ibnd
    integer :: ichname , ibin , ichtrsol , ichtrdpv , idubinsiz
    integer , dimension(2) :: inmlen
    integer , dimension(2) :: idpv
    integer , dimension(2) :: ibinsiz

    integer , dimension(5) :: tyx
    integer , dimension(5) :: tzyx
    integer , dimension(5) :: t10yx
    integer , dimension(5) :: t2yx
    integer , dimension(5) :: tlyx
    integer , dimension(5) :: tcyx
    integer , dimension(5) :: tczyx
    integer , dimension(5) :: tdyx

    if (ctype == 'ATM') then
      ncid = ncatm
      title = 'ICTP Regional Climatic model V4 ATM output'
      iatmrec = 1
    else if (ctype == 'SRF') then
      ncid = ncsrf
      title = 'ICTP Regional Climatic model V4 SRF output'
      isrfrec = 1
    else if (ctype == 'SUB') then
      ncid = ncsub
      title = 'ICTP Regional Climatic model V4 SUB output'
      isubrec = 1
    else if (ctype == 'RAD') then
      ncid = ncrad
      title = 'ICTP Regional Climatic model V4 RAD output'
      iradrec = 1
    else if (ctype == 'CHE') then
      ncid = ncche
      title = 'ICTP Regional Climatic model V4 CHE output'
      icherec = 1
    else if (ctype == 'LAK') then
      ncid = nclak
      title = 'ICTP Regional Climatic model V4 LAK output'
      ilakrec = 1
    else
      write (aline,*) 'UNKNOWN IO TYPE : ', ctype
      call say
      write (aline,*) 'NOTHING TO DO'
      call say
      return
    end if

    call close_common(ncid, ctype)

    write (fterr, '(a3,a)') ctype, ' FILE'
    write (fbname,'(a,a,i10)') trim(ctype), '.', toint10(idate)
    ofname = trim(dirout)//pthsep//trim(domname)// &
             '_'//trim(fbname)//'.nc'
    ctime = tochar(cordex_refdate)
    write (aline, *) 'Opening new output file ', trim(ofname)
    call say

#ifdef NETCDF4_HDF5
    istatus = nf90_create(ofname, &
              ior(ior(nf90_clobber,nf90_hdf5),nf90_classic_model),ncid)
#else
    istatus = nf90_create(ofname, nf90_clobber, ncid)
#endif
    call check_ok(__FILE__,__LINE__,('Error create file '//trim(ofname)), fterr)

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
    istatus = nf90_put_att(ncid, nf90_global, 'grid_size_in_meters', ds*d_1000)
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
!         ADD RUN PARAMETERS
!
    istatus = nf90_put_att(ncid, nf90_global, 'model_IPCC_scenario', scenario)
    call check_ok(__FILE__,__LINE__,'Error add scenario', fterr)
    call cdumlbcs
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_boundary_conditions' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add lbcs', fterr)
    call cdumcums
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_cumulous_convection_scheme' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add icup', fterr)
    if (icup == 2 .or. icup == 99 .or. icup == 98) then
      call cdumcumcl
      istatus = nf90_put_att(ncid, nf90_global,  &
            'model_convective_closure_assumption' , trim(cdum))
      call check_ok(__FILE__,__LINE__,'Error add igcc', fterr)
    end if
    call cdumpbl
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_boundary_layer_scheme' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add ibltyp', fterr)
    call cdummoist
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_moist_physics_scheme' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add ipptls', fterr)
    call cdumocnflx
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_ocean_flux_scheme' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add iocnflx', fterr)
    call cdumpgfs
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_pressure_gradient_force_scheme' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add ipgf', fterr)
    call cdumemiss
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_use_emission_factor' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add iemiss', fterr)
    call cdumlakes
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_use_lake_model' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add lakemod', fterr)
    call cdumchems
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_chemistry' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add ichem', fterr)
    call cdumdcsst
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_diurnal_cycle_sst' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add dcsst', fterr)
    call cdumseaice
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_seaice_effect' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add seaice', fterr)
    call cdumdesseas
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_seasonal_desert_albedo_effect' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add desseas', fterr)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_simulation_initial_start' , tochar(globidate1))
    call check_ok(__FILE__,__LINE__,'Error add globidate1', fterr)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_simulation_start' , tochar(idate1))
    call check_ok(__FILE__,__LINE__,'Error add idate1', fterr)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_simulation_expected_end' , tochar(idate2))
    call check_ok(__FILE__,__LINE__,'Error add idate2', fterr)
    if (ifrest) then
      cdum = 'Yes'
    else
      cdum = 'No'
    end if
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_simulation_is_a_restart' , trim(cdum))
    call check_ok(__FILE__,__LINE__,'Error add ifrest', fterr)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_timestep_in_seconds' , dt)
    call check_ok(__FILE__,__LINE__,'Error add dt', fterr)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_timestep_in_minutes_solar_rad_calc' , dtrad)
    call check_ok(__FILE__,__LINE__,'Error add dtrad', fterr)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_timestep_in_seconds_bats_calc' , dtsrf)
    call check_ok(__FILE__,__LINE__,'Error add dtsrf', fterr)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_timestep_in_hours_radiation_calc' , dtabem)
    call check_ok(__FILE__,__LINE__,'Error add dtabem', fterr)
    istatus = nf90_put_att(ncid, nf90_global,  &
            'model_timestep_in_hours_boundary_input' , ibdyfrq)
    call check_ok(__FILE__,__LINE__,'Error add dtabem', fterr)
!
!         ADD DIMENSIONS
!
    if (ctype == 'SUB') then
      istatus = nf90_def_dim(ncid, 'iy', o_nig, idims(2))
      call check_ok(__FILE__,__LINE__,'Error create dim iy', fterr)
      istatus = nf90_def_dim(ncid, 'jx', o_njg, idims(1))
      call check_ok(__FILE__,__LINE__,'Error create dim jx', fterr)
    else
      istatus = nf90_def_dim(ncid, 'iy', o_ni, idims(2))
      call check_ok(__FILE__,__LINE__,'Error create dim iy', fterr)
      istatus = nf90_def_dim(ncid, 'jx', o_nj, idims(1))
      call check_ok(__FILE__,__LINE__,'Error create dim jx', fterr)
    end if
    istatus = nf90_def_dim(ncid, 'time', nf90_unlimited, idims(3))
    call check_ok(__FILE__,__LINE__,'Error create dim time', fterr)
    istatus = nf90_def_dim(ncid, 'kz', kz, idims(4))
    call check_ok(__FILE__,__LINE__,'Error create dim kz', fterr)
!
!         OUT TYPE DEPENDENT DIMENSIONS
!
    if (ctype == 'SRF' .or. ctype == 'SUB') then
      istatus = nf90_def_dim(ncid, 'm10', 1, idims(5))
      call check_ok(__FILE__,__LINE__,'Error create dim m10', fterr)
      istatus = nf90_def_dim(ncid, 'm2', 1, idims(6))
      call check_ok(__FILE__,__LINE__,'Error create dim m2', fterr)
      istatus = nf90_def_dim(ncid, 'soil_layer', 2, idims(7))
      call check_ok(__FILE__,__LINE__,'Error create dim soil_layer', fterr)
      istatus = nf90_def_var(ncid, 'm10', nf90_float, idims(5), isrvvar(1))
      call check_ok(__FILE__,__LINE__,'Error add var m10', fterr)
      istatus = nf90_put_att(ncid, isrvvar(1), 'standard_name', 'altitude')
      call check_ok(__FILE__,__LINE__,'Error add m10 standard_name', fterr)
      istatus = nf90_put_att(ncid, isrvvar(1), 'long_name', &
                         'Convenience 10 m elevation level')
      call check_ok(__FILE__,__LINE__,'Error add m10 long_name', fterr)
      istatus = nf90_put_att(ncid, isrvvar(1), 'positive', 'up')
      call check_ok(__FILE__,__LINE__,'Error add m10 positive', fterr)
      istatus = nf90_put_att(ncid, isrvvar(1), 'units', 'm')
      call check_ok(__FILE__,__LINE__,'Error add m10 units', fterr)
      istatus = nf90_def_var(ncid, 'm2', nf90_float, idims(6), isrvvar(2))
      call check_ok(__FILE__,__LINE__,'Error add var m2', fterr)
      istatus = nf90_put_att(ncid, isrvvar(2), 'standard_name', 'altitude')
      call check_ok(__FILE__,__LINE__,'Error add m2 standard_name', fterr)
      istatus = nf90_put_att(ncid, isrvvar(2), 'long_name', &
                         'Convenience 2 m elevation level')
      call check_ok(__FILE__,__LINE__,'Error add m2 long_name', fterr)
      istatus = nf90_put_att(ncid, isrvvar(2), 'positive', 'up')
      call check_ok(__FILE__,__LINE__, 'Error add m2 positive', fterr)
      istatus = nf90_put_att(ncid, isrvvar(2), 'units', 'm')
      call check_ok(__FILE__,__LINE__, 'Error add m2 units', fterr)
      istatus = nf90_def_var(ncid, 'layer', nf90_float, idims(7), isrvvar(3))
      call check_ok(__FILE__,__LINE__,'Error add var layer', fterr)
      istatus = nf90_put_att(ncid, isrvvar(3), 'standard_name', &
                         'model_level_number')
      call check_ok(__FILE__,__LINE__,'Error add layer standard_name', fterr)
      istatus = nf90_put_att(ncid, isrvvar(3), 'long_name', &
                         'Surface and root zone')
      call check_ok(__FILE__,__LINE__,'Error add layer long_name', fterr)
      istatus = nf90_put_att(ncid, isrvvar(3), 'positive', 'down')
      call check_ok(__FILE__,__LINE__,'Error add layer positive', fterr)
      istatus = nf90_put_att(ncid, isrvvar(3), 'units', '1')
      call check_ok(__FILE__,__LINE__,'Error add layer units', fterr)
    end if
    if (ctype == 'SRF') then
      istatus = nf90_def_dim(ncid, 'nv', 2, idims(8))
      call check_ok(__FILE__,__LINE__,'Error create dim nv', fterr)
    end  if
    if (ctype == 'CHE') then
      istatus = nf90_def_dim(ncid, 'tracer', ntr, idims(9))
      call check_ok(__FILE__,__LINE__,'Error create dim tracer', fterr)
      istatus = nf90_def_dim(ncid, 'dust', nbin, ibin)
      call check_ok(__FILE__,__LINE__,'Error create dim dust', fterr)
      istatus = nf90_def_dim(ncid, 'bnd', 2, ibnd)
      call check_ok(__FILE__,__LINE__,'Error create dim dust', fterr)
      istatus = nf90_def_dim(ncid, 'namelen', 6, inmlen(1))
      call check_ok(__FILE__,__LINE__,'Error create dim namelen', fterr)
      inmlen(2) = idims(9)
      idpv(1) = idims(9)
      idpv(2) = ibnd
      ibinsiz(1) = ibin
      ibinsiz(2) = ibnd
    end if
    if (ctype == 'LAK') then
      istatus = nf90_def_dim(ncid, 'depth', ndpmax, idims(10))
      call check_ok(__FILE__,__LINE__,'Error create dim depth', fterr)
    end if
    istatus = nf90_def_var(ncid, 'rcm_map', nf90_int, varid=imapvar)
    call check_ok(__FILE__,__LINE__,'Error add var rcm_map', fterr)
    if (iproj == 'LAMCON') then
      istatus = nf90_put_att(ncid, imapvar, &
                   'grid_mapping_name', 'lambert_conformal_conic')
      call check_ok(__FILE__,__LINE__, &
                    'Error add rcm_map grid_mapping_name',fterr)
    else if (iproj == 'POLSTR') then
      istatus = nf90_put_att(ncid, imapvar, &
                   'grid_mapping_name', 'stereographic')
      call check_ok(__FILE__,__LINE__, &
                    'Error add rcm_map grid_mapping_name',fterr)
    else if (iproj == 'NORMER') then
      istatus = nf90_put_att(ncid, imapvar, &
                   'grid_mapping_name', 'mercator')
      call check_ok(__FILE__,__LINE__, &
                    'Error add rcm_map grid_mapping_name',fterr)
    else if (iproj == 'ROTMER') then
      istatus = nf90_put_att(ncid, imapvar, &
            'grid_mapping_name', 'rotated_latitude_longitude')
      call check_ok(__FILE__,__LINE__, &
                    'Error add rcm_map grid_mapping_name',fterr)
    end if
    istatus = nf90_put_att(ncid, imapvar, 'grid_size_in_meters', ds*d_1000)
    call check_ok(__FILE__,__LINE__,'Error add rcm_map gridsize', fterr)
    istatus = nf90_put_att(ncid, imapvar,   &
                 'latitude_of_projection_origin', clat)
    call check_ok(__FILE__,__LINE__,'Error add rcm_map clat', fterr)
    istatus = nf90_put_att(ncid, imapvar,   &
                 'longitude_of_projection_origin', clon)
    call check_ok(__FILE__,__LINE__,'Error add rcm_map clon', fterr)
    istatus = nf90_put_att(ncid, imapvar,   &
                 'longitude_of_central_meridian', clon)
    call check_ok(__FILE__,__LINE__,'Error add rcm_map gmtllon', fterr)
    if (iproj == 'ROTMER') then
      istatus = nf90_put_att(ncid, imapvar, 'grid_north_pole_latitude', plat)
      call check_ok(__FILE__,__LINE__,'Error add rcm_map plat', fterr)
      istatus = nf90_put_att(ncid, imapvar, 'grid_north_pole_longitude', plon)
      call check_ok(__FILE__,__LINE__,'Error add rcm_map plon', fterr)
    else if (iproj == 'LAMCON') then
      trlat(1) = real(truelatl)
      trlat(2) = real(truelath)
      istatus = nf90_put_att(ncid, imapvar, 'standard_parallel', trlat)
      call check_ok(__FILE__,__LINE__,'Error add rcm_map truelat', fterr)
    else if (iproj == 'NORMER') then
      istatus = nf90_put_att(ncid, imapvar, 'standard_parallel', clat)
      call check_ok(__FILE__,__LINE__,'Error add rcm_map truelat', fterr)
    else if (iproj == 'POLSTR') then
      trlat(1) = 1.0
      istatus = nf90_put_att(ncid, imapvar, &
                 'scale_factor_at_projection_origin', trlat(1:1))
      call check_ok(__FILE__,__LINE__,'Error add rcm_map scfac', fterr)
    end if
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
    if (ctype == 'LAK') then
      istatus = nf90_def_var(ncid, 'depth', nf90_float, idims(10), izvar(3))
      call check_ok(__FILE__,__LINE__,'Error add var depth', fterr)
      istatus = nf90_put_att(ncid, izvar(3), 'standard_name', 'depth')
      call check_ok(__FILE__,__LINE__,'Error add depth standard_name', fterr)
      istatus = nf90_put_att(ncid, izvar(3), 'long_name', &
                           'Depth below surface')
      call check_ok(__FILE__,__LINE__,'Error add depth long_name', fterr)
      istatus = nf90_put_att(ncid, izvar(3), 'units', 'm')
      call check_ok(__FILE__,__LINE__,'Error add depth units', fterr)
      istatus = nf90_put_att(ncid, izvar(3), 'axis', 'Z')
      call check_ok(__FILE__,__LINE__,'Error add depth axis', fterr)
      istatus = nf90_put_att(ncid, izvar(3), 'positive', 'down')
      call check_ok(__FILE__,__LINE__,'Error add depth positive', fterr)
    end if
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
    if (ctype == 'SRF') then
      istatus = nf90_put_att(ncid, itvar, 'bounds', 'tbnds')
      call check_ok(__FILE__,__LINE__,'Error add time bounds', fterr)
    end if
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

    tyx = (/idims(1),idims(2),idims(3),-1,-1/)
    tzyx = (/idims(1),idims(2),idims(4),idims(3),-1/)
    t10yx = (/idims(1),idims(2),idims(5),idims(3),-1/)
    t2yx = (/idims(1),idims(2),idims(6),idims(3),-1/)
    tlyx = (/idims(1),idims(2),idims(7),idims(3),-1/)
    tcyx = (/idims(1),idims(2),idims(9),idims(3),-1/)
    tczyx = (/idims(1),idims(2),idims(4),idims(9),idims(3)/)
    tdyx = (/idims(1),idims(2),idims(10),idims(3),-1/)

    if (ctype == 'ATM') then
      iatmvar = -1
      iatmvar(1) = itvar
      iatmvar(2) = illtpvar(5)
      call addvara(ncid,ctype,tzyx,.false.,3)
      call addvara(ncid,ctype,tzyx,.false.,4)
      call addvara(ncid,ctype,tzyx,.false.,5)
      call addvara(ncid,ctype,tzyx,.false.,6)
      call addvara(ncid,ctype,tzyx,.false.,7)
      call addvara(ncid,ctype,tzyx,.false.,8)
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        call addvara(ncid,ctype,tzyx,.false.,9)
        call addvara(ncid,ctype,tzyx,.false.,10)
        call addvara(ncid,ctype,tzyx,.false.,11)
      end if
      call addvara(ncid,ctype,tyx,.false.,12)
      call addvara(ncid,ctype,tyx,.false.,13)
      call addvara(ncid,ctype,tyx,.true.,14)
      call addvara(ncid,ctype,tyx,.true.,15)
    else if (ctype == 'SRF') then
      isrfvar = -1
      isrfvar(1) = itvar
      istatus = nf90_def_var(ncid, 'tbnds', nf90_double, &
                             (/idims(8),idims(3)/), isrfvar(2))
      call check_ok(__FILE__,__LINE__,'Error add var tbnds', fterr)
      istatus = nf90_put_att(ncid, isrfvar(2), &
                             'calendar', calstr(idate%calendar))
      call check_ok(__FILE__,__LINE__,'Error add tbnds calendar', fterr)
      istatus = nf90_put_att(ncid, isrfvar(2), 'units', 'hours since '//ctime)
      call check_ok(__FILE__,__LINE__,'Error add tbnds units', fterr)
      isrfvar(3) = illtpvar(5)
      call addvara(ncid,ctype,t10yx,.false.,4)
      call addvara(ncid,ctype,t10yx,.false.,5)
      call addvara(ncid,ctype,tyx,.false.,6)
      call addvara(ncid,ctype,tyx,.false.,7)
      call addvara(ncid,ctype,tyx,.true.,8)
      call addvara(ncid,ctype,t2yx,.false.,9)
      call addvara(ncid,ctype,t2yx,.false.,10)
      call addvara(ncid,ctype,tlyx,.true.,11)
      call addvara(ncid,ctype,tyx,.false.,12)
      call addvara(ncid,ctype,tyx,.false.,13)
      call addvara(ncid,ctype,tyx,.true.,14)
      call addvara(ncid,ctype,tyx,.true.,15)
      call addvara(ncid,ctype,tyx,.false.,16)
      call addvara(ncid,ctype,tyx,.false.,17)
      call addvara(ncid,ctype,tyx,.false.,18)
      call addvara(ncid,ctype,tyx,.false.,19)
      call addvara(ncid,ctype,tyx,.false.,20)
      call addvara(ncid,ctype,tyx,.false.,21)
      call addvara(ncid,ctype,tyx,.false.,22)
      write (cmethodmax, '(a,i3,a)') 'time: maximum (interval: ', &
                     idnint(srffrq) , ' hours)'
      write (cmethodmin, '(a,i3,a)') 'time: minimum (interval: ', &
                     idnint(srffrq) , ' hours)'
      call addvara(ncid,ctype,tyx,.false.,23)
      istatus = nf90_put_att(ncid, isrfvar(23), 'cell_methods', cmethodmax)
      call check_ok(__FILE__,__LINE__,'Error add tgmax cell_methods', fterr)
      call addvara(ncid,ctype,tyx,.false.,24)
      istatus = nf90_put_att(ncid, isrfvar(24), 'cell_methods', cmethodmin)
      call check_ok(__FILE__,__LINE__,'Error add tgmin cell_methods', fterr)
      call addvara(ncid,ctype,t2yx,.false.,25)
      istatus = nf90_put_att(ncid, isrfvar(25), 'cell_methods', cmethodmax)
      call check_ok(__FILE__,__LINE__,'Error add t2max cell_methods', fterr)
      call addvara(ncid,ctype,t2yx,.false.,26)
      istatus = nf90_put_att(ncid, isrfvar(26), 'cell_methods', cmethodmin)
      call check_ok(__FILE__,__LINE__,'Error add t2min cell_methods', fterr)
      call addvara(ncid,ctype,t10yx,.false.,27)
      istatus = nf90_put_att(ncid, isrfvar(27), 'cell_methods', cmethodmax)
      call check_ok(__FILE__,__LINE__,'Error add w10max cell_methods', fterr)
      call addvara(ncid,ctype,tyx,.false.,28)
      istatus = nf90_put_att(ncid, isrfvar(28), 'cell_methods', cmethodmin)
      call check_ok(__FILE__,__LINE__,'Error add ps_min cell_methods', fterr)
      call addvara(ncid,ctype,tyx,.false.,29)
      call addvara(ncid,ctype,tyx,.false.,30)
      if ( iseaice == 1 .or. lakemod == 1 ) then
        srf_variables(31)%enabled = .true.
        call addvara(ncid,ctype,tyx,.false.,31)
      end if
    else if (ctype == 'SUB') then
      isubvar = -1
      isubvar(1) = itvar
      isubvar(2) = illtpvar(5)
      call addvara(ncid,ctype,t10yx,.false.,3)
      call addvara(ncid,ctype,t10yx,.false.,4)
      call addvara(ncid,ctype,tyx,.false.,5)
      call addvara(ncid,ctype,tyx,.false.,6)
      call addvara(ncid,ctype,tyx,.true.,7)
      call addvara(ncid,ctype,t2yx,.false.,8)
      call addvara(ncid,ctype,t2yx,.false.,9)
      call addvara(ncid,ctype,tlyx,.true.,10)
      call addvara(ncid,ctype,tyx,.false.,11)
      call addvara(ncid,ctype,tyx,.false.,12)
      call addvara(ncid,ctype,tyx,.true.,13)
      call addvara(ncid,ctype,tyx,.true.,14)
      call addvara(ncid,ctype,tyx,.false.,15)
      call addvara(ncid,ctype,tyx,.false.,16)
    else if (ctype == 'RAD') then
      iradvar = -1
      iradvar(1) = itvar
      iradvar(2) = illtpvar(5)
      call addvara(ncid,ctype,tzyx,.false.,3)
      call addvara(ncid,ctype,tzyx,.false.,4)
      call addvara(ncid,ctype,tzyx,.false.,5)
      call addvara(ncid,ctype,tzyx,.false.,6)
      call addvara(ncid,ctype,tyx,.false.,7)
      call addvara(ncid,ctype,tyx,.false.,8)
      call addvara(ncid,ctype,tyx,.false.,9)
      call addvara(ncid,ctype,tyx,.false.,10)
      call addvara(ncid,ctype,tyx,.false.,11)
      call addvara(ncid,ctype,tyx,.false.,12)
      call addvara(ncid,ctype,tyx,.false.,13)
      call addvara(ncid,ctype,tyx,.false.,14)
      call addvara(ncid,ctype,tyx,.false.,15)
    else if (ctype == 'CHE') then
      istatus = nf90_def_var(ncid, 'chtrname', nf90_char, &
                             inmlen, ichname)
      call check_ok(__FILE__,__LINE__,'Error add var chtrname', fterr)
      istatus = nf90_def_var(ncid, 'chtrsol', nf90_double, &
                             idims(9), ichtrsol)
      call check_ok(__FILE__,__LINE__,'Error add var chtrsol', fterr)
      istatus = nf90_def_var(ncid, 'chtrdpv', nf90_double, &
                             idpv, ichtrdpv)
      call check_ok(__FILE__,__LINE__,'Error add var chtrdpv', fterr)
      istatus = nf90_def_var(ncid, 'dustbinsiz', nf90_double, &
                             ibinsiz, idubinsiz)
      call check_ok(__FILE__,__LINE__,'Error add var dustbinsiz', fterr)
      ichevar = -1
      ichevar(1) = itvar
      ichevar(2) = illtpvar(5)
      call addvara(ncid,ctype,tczyx,.false.,3)
      call addvara(ncid,ctype,tzyx,.false.,4)
      call addvara(ncid,ctype,tzyx,.false.,5)
      call addvara(ncid,ctype,tzyx,.false.,6)
      call addvara(ncid,ctype,tcyx,.false.,7)
      call addvara(ncid,ctype,tcyx,.false.,8)
      call addvara(ncid,ctype,tcyx,.false.,9)
      call addvara(ncid,ctype,tcyx,.false.,10)
      call addvara(ncid,ctype,tcyx,.false.,11)
      call addvara(ncid,ctype,tcyx,.false.,12)
      call addvara(ncid,ctype,tcyx,.false.,13)
      call addvara(ncid,ctype,tyx,.false.,14)
      call addvara(ncid,ctype,tyx,.false.,15)
      call addvara(ncid,ctype,tyx,.false.,16)
      call addvara(ncid,ctype,tyx,.false.,17)
    else if (ctype == 'LAK') then
      ilakvar = -1
      ilakvar(1) = itvar
      ilakvar(2) = illtpvar(5)
      call addvara(ncid,ctype,tyx,.false.,3)
      call addvara(ncid,ctype,tyx,.false.,4)
      call addvara(ncid,ctype,tyx,.true.,5)
      call addvara(ncid,ctype,tyx,.false.,6)
      call addvara(ncid,ctype,tyx,.false.,7)
      call addvara(ncid,ctype,tyx,.false.,8)
      call addvara(ncid,ctype,tyx,.false.,9)
      call addvara(ncid,ctype,tyx,.false.,10)
      call addvara(ncid,ctype,tyx,.false.,11)
      call addvara(ncid,ctype,tyx,.false.,12)
      call addvara(ncid,ctype,tyx,.true.,13)
      call addvara(ncid,ctype,tyx,.true.,14)
      call addvara(ncid,ctype,tyx,.true.,15)
      call addvara(ncid,ctype,tdyx,.true.,16)
    end if

    istatus = nf90_enddef(ncid)
    call check_ok(__FILE__,__LINE__,'Error End Definitions NetCDF output',fterr)

    istatus = nf90_put_var(ncid, izvar(1), hsigma)
    call check_ok(__FILE__,__LINE__,'Error var sigma write', fterr)
    hptop = real(ptop*d_10)
    istatus = nf90_put_var(ncid, izvar(2), hptop)
    call check_ok(__FILE__,__LINE__,'Error var ptop write', fterr)
    if (ctype == 'LAK') then
      do i = 1 , ndpmax
        depth(i) = real(i)
      end do
      istatus = nf90_put_var(ncid, izvar(3), depth)
      call check_ok(__FILE__,__LINE__,'Error var depth write', fterr)
    end if
    if (ctype == 'SUB') then
      yiy(1) = -real((dble((o_nig-1)-1)/2.0D0)*ds)
      xjx(1) = -real((dble((o_njg-1)-1)/2.0D0)*ds)
      do i = 2 , o_nig
        yiy(i) = yiy(i-1)+real(ds)
      end do
      do j = 2 , o_njg
        xjx(j) = xjx(j-1)+real(ds)
      end do
      istatus = nf90_put_var(ncid, ivvar(1), yiy(1:o_nig))
      call check_ok(__FILE__,__LINE__,'Error var iy write', fterr)
      istatus = nf90_put_var(ncid, ivvar(2), xjx(1:o_njg))
      call check_ok(__FILE__,__LINE__,'Error var jx write', fterr)
      istatus = nf90_put_var(ncid, illtpvar(1), ioxlat_s)
      call check_ok(__FILE__,__LINE__,'Error var xlat write', fterr)
      istatus = nf90_put_var(ncid, illtpvar(2), ioxlon_s)
      call check_ok(__FILE__,__LINE__,'Error var xlon write', fterr)
      istatus = nf90_put_var(ncid, illtpvar(3), iotopo_s)
      call check_ok(__FILE__,__LINE__,'Error var topo write', fterr)
      istatus = nf90_put_var(ncid, illtpvar(4), iomask_s)
      call check_ok(__FILE__,__LINE__,'Error var mask write', fterr)
    else
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
    end if
    if (ctype == 'SRF' .or. ctype == 'SUB') then
      rdum1 = 10.0
      istatus = nf90_put_var(ncid, isrvvar(1), rdum1)
      call check_ok(__FILE__,__LINE__,'Error var m10 write', fterr)
      rdum1 = 2.0
      istatus = nf90_put_var(ncid, isrvvar(2), rdum1)
      call check_ok(__FILE__,__LINE__,'Error var m2 write', fterr)
      rdum2(1) = 0.0
      rdum2(2) = 1.0
      istatus = nf90_put_var(ncid, isrvvar(3), rdum2)
      call check_ok(__FILE__,__LINE__,'Error var layer write', fterr)
    end if
    if (ctype == 'CHE') then
      istatus = nf90_put_var(ncid, ichname, chtrname)
      call check_ok(__FILE__,__LINE__,'Error var chtrname write', fterr)
      istatus = nf90_put_var(ncid, ichtrsol, chtrsol)
      call check_ok(__FILE__,__LINE__,'Error var chtrsol write', fterr)
      istatus = nf90_put_var(ncid, ichtrdpv, chtrdpv)
      call check_ok(__FILE__,__LINE__,'Error var chtrdpv write', fterr)
      istatus = nf90_put_var(ncid, idubinsiz, dustbsiz)
      call check_ok(__FILE__,__LINE__,'Error var dustbsiz write', fterr)
    end if

    if ( debug_level > 2 ) then
      istatus = nf90_sync(ncid)
      call check_ok(__FILE__,__LINE__,'Error initial sync', fterr)
    end if

    if (ctype == 'ATM') then
      ncatm = ncid
    else if (ctype == 'SRF') then
      ncsrf = ncid
    else if (ctype == 'SUB') then
      ncsub = ncid
    else if (ctype == 'RAD') then
      ncrad = ncid
    else if (ctype == 'CHE') then
      ncche = ncid
    else if (ctype == 'LAK') then
      nclak = ncid
    end if
  end subroutine prepare_common_out

  subroutine addvara(ncid,ctype,idims,lmiss,nvar)
    use netcdf
    implicit none
    integer , intent(in) :: ncid
    character(3) , intent(in) :: ctype
    integer , dimension(5) , intent(in) :: idims
    logical , intent(in) :: lmiss
    integer , intent(in) :: nvar

    character(len=8)   :: vname
    character(len=128) :: vst , vln
    character(len=16)  :: vuni
    logical :: lreq
    integer :: ivar

    integer :: i , ndims

    ndims = 0
    do i = 1 , 5
      if (idims(i) > 0) ndims = ndims+1
    end do

    select case (ctype)
      case ('ATM')
        vname = atm_variables(nvar)%vname
        vst   = atm_variables(nvar)%vstd_name
        vln   = atm_variables(nvar)%vdesc
        vuni  = atm_variables(nvar)%vunit
        lreq  = atm_variables(nvar)%enabled
      case ('SRF')
        vname = srf_variables(nvar)%vname
        vst   = srf_variables(nvar)%vstd_name
        vln   = srf_variables(nvar)%vdesc
        vuni  = srf_variables(nvar)%vunit
        lreq  = srf_variables(nvar)%enabled
      case ('SUB')
        vname = sub_variables(nvar)%vname
        vst   = sub_variables(nvar)%vstd_name
        vln   = sub_variables(nvar)%vdesc
        vuni  = sub_variables(nvar)%vunit
        lreq  = sub_variables(nvar)%enabled
      case ('LAK')
        vname = lak_variables(nvar)%vname
        vst   = lak_variables(nvar)%vstd_name
        vln   = lak_variables(nvar)%vdesc
        vuni  = lak_variables(nvar)%vunit
        lreq  = lak_variables(nvar)%enabled
      case ('RAD')
        vname = rad_variables(nvar)%vname
        vst   = rad_variables(nvar)%vstd_name
        vln   = rad_variables(nvar)%vdesc
        vuni  = rad_variables(nvar)%vunit
        lreq  = rad_variables(nvar)%enabled
      case ('CHE')
        vname = che_variables(nvar)%vname
        vst   = che_variables(nvar)%vstd_name
        vln   = che_variables(nvar)%vdesc
        vuni  = che_variables(nvar)%vunit
        lreq  = che_variables(nvar)%enabled
      case default
        call fatal(__FILE__,__LINE__,ctype//': Not defined output type')
    end select

    if (lreq) then
      cdum = vname
      istatus = nf90_def_var(ncid, cdum, nf90_float, idims(1:ndims), ivar)
      call check_ok(__FILE__,__LINE__,'Error add var '//vname, ctype//' FILE')
#ifdef NETCDF4_HDF5
      istatus = nf90_def_var_deflate(ncid, ivar, 1, 1, 9)
      call check_ok(__FILE__,__LINE__,'Error setting deflate on var '//vname, &
                    ctype//' FILE')
#endif
      cdum = vst
      istatus = nf90_put_att(ncid, ivar, 'standard_name', cdum)
      call check_ok(__FILE__,__LINE__,'Error add '//vname//' standard_name', &
                    ctype//' FILE')
      cdum = vln
      istatus = nf90_put_att(ncid, ivar, 'long_name', cdum)
      call check_ok(__FILE__,__LINE__,'Error add '//vname//' long_name', &
                    ctype//' FILE')
      cdum = vuni
      istatus = nf90_put_att(ncid, ivar, 'units', cdum)
      call check_ok(__FILE__,__LINE__,'Error add '//vname//' units', &
                    ctype//' FILE')
      istatus = nf90_put_att(ncid, ivar, 'coordinates', 'xlat xlon')
      call check_ok(__FILE__,__LINE__,'Error add '//vname//' coord', &
                    ctype//' FILE')
      istatus = nf90_put_att(ncid, ivar, 'grid_mapping', 'rcm_map')
      call check_ok(__FILE__,__LINE__,'Error add '//vname//' grid_mapping', &
                    ctype//' FILE')
      if (lmiss) then
        istatus = nf90_put_att(ncid, ivar, '_FillValue', smissval)
        call check_ok(__FILE__,__LINE__,'Error add '//vname//' coord', &
                      ctype//' FILE')
      end if
      select case (ctype)
        case ('ATM')
          iatmvar(nvar) = ivar
        case ('SRF')
          isrfvar(nvar) = ivar
        case ('SUB')
          isubvar(nvar) = ivar
        case ('LAK')
          ilakvar(nvar) = ivar
        case ('RAD')
          iradvar(nvar) = ivar
        case ('CHE')
          ichevar(nvar) = ivar
        case default
          call fatal(__FILE__,__LINE__,ctype//': Not defined output type')
      end select
    end if
  end subroutine addvara

  subroutine writerec_srf(nx, ny, numbat, fbat, mask , idate)
    use netcdf
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer , intent(in) :: nx , ny , numbat
    real(4) , dimension(nx,ny,numbat) , intent(in) :: fbat
    integer , dimension(iym1,jxm1) , intent(in) :: mask
    integer :: ivar
    integer :: n
    integer , dimension(4) :: istart , icount
    real(8) , dimension(2) :: nctime
    type(rcm_time_interval) :: tdif
    logical :: lskip
    character(len=36) :: ctime

    if (nx /= o_nj .or. ny /= o_ni) then
      write (6,*) 'Error writing record on SRF file'
      write (6,*) 'Expecting layers ', o_nj, 'x', o_ni
      write (6,*) 'Got layers       ', nx, 'x', ny
      call fatal(__FILE__,__LINE__,'DIMENSION MISMATCH')
    end if
    ctime = tochar(idate)

    istart(2) = isrfrec
    istart(1) = 1
    icount(2) = 1
    icount(1) = 2
    tdif = idate-cordex_refdate
    nctime(2) = tohours(tdif)
    nctime(1) = nctime(2) - srffrq
    istatus = nf90_put_var(ncsrf, isrfvar(1), nctime(2:2), &
                           istart(2:2), icount(2:2))
    call check_ok(__FILE__,__LINE__,'Error writing itime '//ctime, 'SRF FILE')
    istatus = nf90_put_var(ncsrf, isrfvar(2), nctime, &
                           istart(1:2), icount(1:2))
    call check_ok(__FILE__,__LINE__,'Error writing tbnds '//ctime, 'SRF FILE')

    ivar = 3
    lskip = .false.
    do n = 1 , numbat
      if (lskip) then
        lskip = .false.
        cycle
      end if
      if ( srf_variables(ivar)%enabled ) then
        if (ivar == ivarname_lookup('SRF', 'u10m')   .or. &
            ivar == ivarname_lookup('SRF', 'v10m')   .or. &
            ivar == ivarname_lookup('SRF', 'w10max') .or. &
            ivar == ivarname_lookup('SRF', 't2m')    .or. &
            ivar == ivarname_lookup('SRF', 'q2m')    .or. &
            ivar == ivarname_lookup('SRF', 't2max')  .or. &
            ivar == ivarname_lookup('SRF', 't2min')) then
          istart(4) = isrfrec
          istart(3) = 1
          istart(2) = 1
          istart(1) = 1
          icount(4) = 1
          icount(3) = 1
          icount(2) = o_ni
          icount(1) = o_nj
          istatus = nf90_put_var(ncsrf, isrfvar(ivar), &
                          fbat(:,:,n), istart, icount)
          call check_ok(__FILE__,__LINE__, &
                        'Error writing '//srf_variables(ivar)%vname// &
                        ' at '//ctime, 'SRF FILE')
        else if (ivar == ivarname_lookup('SRF', 'smw')) then
          istart(4) = isrfrec
          istart(3) = 1
          istart(2) = 1
          istart(1) = 1
          icount(4) = 1
          icount(3) = 1
          icount(2) = o_ni
          icount(1) = o_nj
          istatus = nf90_put_var(ncsrf, isrfvar(ivar), & 
                          fbat(:,:,n), istart, icount)
          call check_ok(__FILE__,__LINE__, &
                        'Error writing '//srf_variables(ivar)%vname// &
                        ' at '//ctime, 'SRF FILE')
          istart(3) = 2
          istatus = nf90_put_var(ncsrf, isrfvar(ivar), & 
                          fbat(:,:,n+1), istart, icount)
          call check_ok(__FILE__,__LINE__, &
                        'Error writing '//srf_variables(ivar)%vname// &
                        ' at '//ctime, 'SRF FILE')
        else
          istart(3) = isrfrec
          istart(2) = 1
          istart(1) = 1
          icount(3) = 1
          icount(2) = o_ni
          icount(1) = o_nj
          istatus = nf90_put_var(ncsrf, isrfvar(ivar), & 
                   fbat(:,:,n), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__, &
                        'Error writing '//srf_variables(ivar)%vname// &
                        ' at '//ctime, 'SRF FILE')
        end if
      end if
      if (ivar == ivarname_lookup('SRF', 'smw')) then
        lskip = .true.
      end if
      ivar = ivar + 1
    end do

    if ( iseaice == 1 .or. lakemod == 1 ) then
      dumio(:,:,1) = real(transpose(mask(o_is:o_ie,o_js:o_je)))
      istart(3) = isrfrec
      istart(2) = 1
      istart(1) = 1
      icount(3) = 1
      icount(2) = o_ni
      icount(1) = o_nj
      istatus = nf90_put_var(ncsrf, isrfvar(31), & 
               dumio(:,:,1), istart(1:3), icount(1:3))
      call check_ok(__FILE__,__LINE__, &
                    'Error writing '//srf_variables(31)%vname// &
                    ' at '//ctime, 'SRF FILE')
    end if

    if ( debug_level > 2 ) then
      istatus = nf90_sync(ncsrf)
      call check_ok(__FILE__,__LINE__,'Error sync at '//ctime, 'SRF FILE')
    end if
    isrfrec = isrfrec + 1
  end subroutine writerec_srf

  subroutine writerec_sub(nx, ny, ns, nsub, fsub, idate)
    use netcdf
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer , intent(in) :: nx , ny , ns , nsub
    real(4) , dimension(ns*ns,nx/ns,ny/ns,nsub) , intent(in) :: fsub
    integer :: ivar
    integer :: n , nxb , nyb
    integer , dimension(4) :: istart , icount
    real(8) , dimension(1) :: nctime
    type(rcm_time_interval) :: tdif
    character(len=36) :: ctime
    logical :: lskip

    nxb = o_njg / nsg
    nyb = o_nig / nsg

    if (nx /= o_njg .or. ny /= o_nig .or. ns /= nsg) then
      write (6,*) 'Error writing record on SUB file'
      write (6,*) 'Expecting layers ', nsg, 'x', o_njg, 'x', o_nig
      write (6,*) 'Got layers       ', ns, 'x', nx, 'x', ny
      call fatal(__FILE__,__LINE__,'DIMENSION MISMATCH')
    end if

    ctime = tochar(idate)

    istart(1) = isubrec
    icount(1) = 1
    tdif = idate-cordex_refdate
    nctime(1) = tohours(tdif)
    istatus = nf90_put_var(ncsub, isubvar(1), nctime, istart(1:1), icount(1:1))
    call check_ok(__FILE__,__LINE__,'Error writing itime '//ctime, 'SUB FILE')
    ivar = 2
    lskip = .false.
    do n = 1 , nsub
      if (lskip) then
        lskip = .false.
        cycle
      end if
      if ( sub_variables(ivar)%enabled ) then
        call reorder(fsub,subio,nxb,nyb,nsg,nsub,n)
        if (ivar == ivarname_lookup('SUB', 'u10m')   .or. &
            ivar == ivarname_lookup('SUB', 'v10m')   .or. &
            ivar == ivarname_lookup('SUB', 't2m')    .or. &
            ivar == ivarname_lookup('SUB', 'q2m') ) then
          istart(4) = isubrec
          istart(3) = 1
          istart(2) = 1
          istart(1) = 1
          icount(4) = 1
          icount(3) = 1
          icount(2) = o_nig
          icount(1) = o_njg
          istatus = nf90_put_var(ncsub, isubvar(ivar), subio, istart, icount)
          call check_ok(__FILE__,__LINE__, &
                        'Error writing '//sub_variables(ivar)%vname// &
                        ' at '//ctime, 'SUB FILE')
        else if (ivar == ivarname_lookup('SUB', 'smw')) then
          istart(4) = isubrec
          istart(3) = 1
          istart(2) = 1
          istart(1) = 1
          icount(4) = 1
          icount(3) = 1
          icount(2) = o_nig
          icount(1) = o_njg
          istatus = nf90_put_var(ncsub, isubvar(ivar), subio, istart, icount)
          call check_ok(__FILE__,__LINE__, &
                        'Error writing '//sub_variables(ivar)%vname// &
                        ' at '//ctime, 'SUB FILE')
          istart(3) = 2
          call reorder(fsub,subio,nxb,nyb,nsg,nsub,n+1)
          istatus = nf90_put_var(ncsub, isubvar(ivar), subio, istart, icount)
          call check_ok(__FILE__,__LINE__, &
                        'Error writing '//sub_variables(ivar)%vname// &
                        ' at '//ctime, 'SUB FILE')
        else
          istart(3) = isubrec
          istart(2) = 1
          istart(1) = 1
          icount(3) = 1
          icount(2) = o_nig
          icount(1) = o_njg
          istatus = nf90_put_var(ncsub, isubvar(ivar), & 
                                 subio, istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__, &
                        'Error writing '//sub_variables(ivar)%vname// &
                        ' at '//ctime, 'SUB FILE')
        end if
      end if
      if (ivar == ivarname_lookup('SUB', 'smw')) then
        lskip = .true.
      end if
      ivar = ivar + 1
    end do
    if ( debug_level > 2 ) then
      istatus = nf90_sync(ncsub)
      call check_ok(__FILE__,__LINE__,'Error sync at '//ctime, 'SUB FILE')
    end if
    isubrec = isubrec + 1
  end subroutine writerec_sub

  subroutine reorder(fdp,fsp,nx,ny,nz,nn,n)
    implicit none
    integer :: ny , nx , nz , nn , n
    real(4) , dimension(nz*nz,nx,ny,nn) :: fdp
    real(4) , dimension(nx*nz,ny*nz) :: fsp
    intent (in) fdp , ny , nx , nz , nn , n
    intent (out) fsp

    integer :: i , ii , j , jj , k
!
    do j = 1 , nx*nz
      do i = 1 , ny*nz
        jj = mod(j,nz)
        if ( jj == 0 ) jj = nz
        ii = mod(i,nz)
        if ( ii == 0 ) ii = nz
        k = (jj-1)*nz + ii
        jj = (j+nz-1)/nz
        ii = (i+nz-1)/nz
        fsp(j,i) = fdp(k,jj,ii,n)
      end do
    end do
  end subroutine reorder

  subroutine writerec_rad(nx,ny,nz,nrad3d,nrad2d,frad3d,frad2d,ps,idate)
    use netcdf
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer , intent(in) :: nx , ny , nz , nrad3d , nrad2d
    real(4) , dimension(nx,ny,nz,nrad3d) , intent(in) :: frad3d
    real(4) , dimension(nx,ny,nrad2d) , intent(in) :: frad2d
    real(4) , dimension(nx,ny) , intent(in) :: ps
    integer :: ivar
    integer :: n
    integer , dimension(4) :: istart , icount
    real(8) , dimension(1) :: nctime
    type(rcm_time_interval) :: tdif
    character(len=36) :: ctime

    if (nx /= o_nj .or. ny /= o_ni .or. nz /= o_nz) then
      write (6,*) 'Error writing record on RAD file'
      write (6,*) 'Expecting layers ', o_nz, 'x', o_nj, 'x', o_ni
      write (6,*) 'Got layers       ', nz, 'x', nx, 'x', ny
      call fatal(__FILE__,__LINE__,'DIMENSION MISMATCH')
    end if

    ctime = tochar(idate)

    istart(1) = iradrec
    icount(1) = 1
    tdif = idate-cordex_refdate
    nctime(1) = tohours(tdif)
    istatus = nf90_put_var(ncrad, iradvar(1), nctime, istart(1:1), icount(1:1))
    call check_ok(__FILE__,__LINE__,'Error writing itime '//ctime, 'RAD FILE')

    istart(3) = iradrec
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = o_ni
    icount(1) = o_nj
    istatus = nf90_put_var(ncrad, iradvar(2), ps, istart(1:3), icount(1:3))
    call check_ok(__FILE__,__LINE__,'Error writing ps at '//ctime, 'RAD FILE')

    ivar = 3
    do n = 1 , nrad3d
      if ( rad_variables(ivar)%enabled ) then
        istart(4) = iradrec
        istart(3) = 1
        istart(2) = 1
        istart(1) = 1
        icount(4) = 1
        icount(3) = o_nz
        icount(2) = o_ni
        icount(1) = o_nj
        istatus = nf90_put_var(ncrad, iradvar(ivar), &
                               frad3d(:,:,:,n), istart, icount)
        call check_ok(__FILE__,__LINE__, &
                      'Error writing '//rad_variables(ivar)%vname// &
                      ' at '//ctime, 'RAD FILE')
      end if
      ivar = ivar + 1
    end do
    do n = 1 , nrad2d
      if ( rad_variables(ivar)%enabled ) then
        istart(3) = iradrec
        istart(2) = 1
        istart(1) = 1
        icount(3) = 1
        icount(2) = o_ni
        icount(1) = o_nj
        istatus = nf90_put_var(ncrad, iradvar(ivar), & 
                               frad2d(:,:,n), istart(1:3), icount(1:3))
        call check_ok(__FILE__,__LINE__, &
                      'Error writing '//rad_variables(ivar)%vname// &
                      ' at '//ctime, 'RAD FILE')
      end if
      ivar = ivar + 1
    end do

    if ( debug_level > 2 ) then
      istatus = nf90_sync(ncrad)
      call check_ok(__FILE__,__LINE__,'Error sync at '//ctime, 'RAD FILE')
    end if
    iradrec = iradrec + 1
  end subroutine writerec_rad

  subroutine writerec_atm(nx, ny, nnx, nny, nz, ns, u, v, omega,    &
                          t, qv, qc, tke , kth , kzm , ps, rc, rnc, &
                          tgb, swt, rno, mask, idate)
    use netcdf
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer , intent(in) :: nx , nnx, nny , ny , ns , nz
    real(8) , dimension(ny,nz,nx) , intent(in) :: u
    real(8) , dimension(ny,nz,nx) , intent(in) :: v
    real(8) , dimension(ny,nz,nx) , intent(in) :: omega
    real(8) , dimension(ny,nz,nx) , intent(in) :: t
    real(8) , dimension(ny,nz,nx) , intent(in) :: qv
    real(8) , dimension(ny,nz,nx) , intent(in) :: qc
    real(8) , dimension(ny,nz+1,nx) , intent(in) :: tke
    real(8) , dimension(ny,nz+1,nx) , intent(in) :: kth
    real(8) , dimension(ny,nz+1,nx) , intent(in) :: kzm
    real(8) , dimension(ny,nx) , intent(in) :: ps
    real(8) , dimension(ny,nx) , intent(in) :: rc
    real(8) , dimension(ny,nx) , intent(in) :: rnc
    real(8) , dimension(ns,nny,nnx) , intent(in) :: tgb
    real(8) , dimension(ns,nny,nnx) , intent(in) :: swt
    real(8) , dimension(ns,nny,nnx) , intent(in) :: rno
    integer , dimension(ns,nny,nnx) , intent(in) :: mask
    integer :: i , j , n , ip1 , ip2 , jp1 , jp2 , k
    integer , dimension(4) :: istart , icount
    real(8) , dimension(1) :: nctime
    type(rcm_time_interval) :: tdif
    character(len=36) :: ctime

    if (nx < o_nj .or. ny < o_ni .or. nz > o_nz .or. &
        nnx < o_nj .or. nny < o_ni) then
      write (6,*) 'Error writing record on ATM file'
      write (6,*) 'Expecting layers ', o_nz, 'x', o_nj, 'x', o_ni
      write (6,*) 'Got layers 0     ', nz, 'x', nx, 'x', ny
      write (6,*) 'Got layers 1     ', nz, 'x', nnx, 'x', nny
      call fatal(__FILE__,__LINE__,'DIMENSION MISMATCH')
    end if

    ctime = tochar(idate)

    if (.not. lmaskfill) then
      do n = 1 , ns
        atmsrfmask(n,:,:) = real(transpose(mask(n,o_is:o_ie,o_js:o_je)))
      end do
      where ( atmsrfmask > 0.5 )
        atmsrfmask = 1.0
      elsewhere
        atmsrfmask = 0.0
      end where
      atmsrfsum = sum(atmsrfmask,dim=1)
      lmaskfill = .true.
    end if

    istart(1) = iatmrec
    icount(1) = 1
    tdif = idate-cordex_refdate
    nctime(1) = tohours(tdif)
    istatus = nf90_put_var(ncatm, iatmvar(1), nctime, istart(1:1), icount(1:1))
    call check_ok(__FILE__,__LINE__,'Error writing itime '//ctime, 'ATM FILE')

    dumio(:,:,1) = real((transpose(ps(o_is:o_ie,o_js:o_je)) + ptop)*d_10)
    istart(3) = iatmrec
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = o_ni
    icount(1) = o_nj
    istatus = nf90_put_var(ncatm, iatmvar(2), &
                           dumio(:,:,1), istart(1:3), icount(1:3))
    call check_ok(__FILE__,__LINE__,'Error writing ps at '//ctime, 'ATM FILE')

    istart(4) = iatmrec
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = o_nz
    icount(2) = o_ni
    icount(1) = o_nj

    if ( atm_variables(3)%enabled ) then
      do k = 1 , o_nz
        do i = 1 , o_ni
          ip1 = i+1
          ip2 = i+2
          do j = 1 , o_nj
            if (.not. lwrap) then
              jp1 = j+1
              jp2 = j+2
            else
              jp1 = j
              jp2 = j+1
              if (j == o_nj) jp2 = 1
            end if
            dumio(j,i,k) = real(((u(ip1,k,jp1)+u(ip1,k,jp2) + &
                                  u(ip2,k,jp1)+u(ip2,k,jp2))*d_rfour) / &
                                  ps(ip1,jp1))
          end do
        end do
      end do
      istatus = nf90_put_var(ncatm, iatmvar(3), dumio, istart, icount)
      call check_ok(__FILE__,__LINE__, &
                    'Error writing '//atm_variables(3)%vname// &
                    ' at '//ctime, 'ATM FILE')
    end if

    if ( atm_variables(4)%enabled ) then
      do k = 1 , o_nz
        do i = 1 , o_ni
          ip1 = i+1
          ip2 = i+2
          do j = 1 , o_nj
            if (.not. lwrap) then
              jp1 = j+1
              jp2 = j+2
            else
              jp1 = j
              jp2 = j+1
              if (j == o_nj) jp2 = 1
            end if
            dumio(j,i,k) = real(((v(ip1,k,jp1)+v(ip1,k,jp2) + &
                                  v(ip2,k,jp1)+v(ip2,k,jp2))*d_rfour) / &
                                  ps(ip1,jp1))
          end do
        end do
      end do
      istatus = nf90_put_var(ncatm, iatmvar(4), dumio, istart, icount)
      call check_ok(__FILE__,__LINE__,&
                    'Error writing '//atm_variables(4)%vname//' at '//ctime, &
                    'ATM FILE')
    end if

    if ( atm_variables(5)%enabled ) then
      do k = 1 , o_nz
        do i = 1 , o_ni
          ip1 = i+1
          do j = 1 , o_nj
            if (.not. lwrap) then
              jp1 = j+1
            else
              jp1 = j
            end if
            dumio(j,i,k) = real(omega(ip1,k,jp1)*d_10)
          end do
        end do
      end do
      istatus = nf90_put_var(ncatm, iatmvar(5), dumio, istart, icount)
      call check_ok(__FILE__,__LINE__,&
                    'Error writing '//atm_variables(5)%vname//' at '//ctime, &
                    'ATM FILE')
    end if

    if ( atm_variables(6)%enabled ) then
      do k = 1 , o_nz
        do i = 1 , o_ni
          ip1 = i+1
          do j = 1 , o_nj
            if (.not. lwrap) then
              jp1 = j+1
            else
              jp1 = j
            end if
            dumio(j,i,k) = real(t(ip1,k,jp1)/ps(ip1,jp1))
          end do
        end do
      end do
      istatus = nf90_put_var(ncatm, iatmvar(6), dumio, istart, icount)
      call check_ok(__FILE__,__LINE__,&
                    'Error writing '//atm_variables(6)%vname//' at '//ctime, &
                    'ATM FILE')
    end if

    if ( atm_variables(7)%enabled ) then
      dumio = 0.0
      do k = 1 , o_nz
        do i = 1 , o_ni
          ip1 = i+1
          do j = 1 , o_nj
            if (.not. lwrap) then
              jp1 = j+1
            else
              jp1 = j
            end if
            if (qv(ip1,k,jp1) > dlowval) then
              dumio(j,i,k) = real(qv(ip1,k,jp1)/ps(ip1,jp1))
            end if
          end do
        end do
      end do
      istatus = nf90_put_var(ncatm, iatmvar(7), dumio, istart, icount)
      call check_ok(__FILE__,__LINE__,&
                    'Error writing '//atm_variables(7)%vname//' at '//ctime, &
                    'ATM FILE')
    end if

    if ( atm_variables(8)%enabled ) then
      dumio = 0.0
      do k = 1 , o_nz
        do i = 1 , o_ni
          ip1 = i+1
          do j = 1 , o_nj
            if (.not. lwrap) then
              jp1 = j+1
            else
              jp1 = j
            end if
            if (qc(ip1,k,jp1) > dlowval) then
              dumio(j,i,k) = real(qc(ip1,k,jp1)/ps(ip1,jp1))
            end if
          end do
        end do
      end do
      istatus = nf90_put_var(ncatm, iatmvar(8), dumio, istart, icount)
      call check_ok(__FILE__,__LINE__,&
                    'Error writing '//atm_variables(8)%vname//' at '//ctime, &
                    'ATM FILE')
    end if

    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      if ( atm_variables(9)%enabled ) then
        dumio = 0.0
        do k = 1 , o_nz
          do i = 1 , o_ni
            ip1 = i+1
            do j = 1 , o_nj
              if (.not. lwrap) then
                jp1 = j+1
              else
                jp1 = j
              end if
              if (qc(ip1,k,jp1) > dlowval) then
                dumio(j,i,k) = real(tke(ip1,k,jp1)/ps(ip1,jp1))
              end if
            end do
          end do
        end do
        istatus = nf90_put_var(ncatm, iatmvar(9), dumio, istart, icount)
        call check_ok(__FILE__,__LINE__,&
                      'Error writing '//atm_variables(9)%vname//' at '//ctime, &
                      'ATM FILE')
      end if
      if ( atm_variables(10)%enabled ) then
        dumio = 0.0
        do k = 1 , o_nz
          do i = 1 , o_ni
            ip1 = i+1
            do j = 1 , o_nj
              if (.not. lwrap) then
                jp1 = j+1
              else
                jp1 = j
              end if
              if (qc(ip1,k,jp1) > dlowval) then
                dumio(j,i,k) = real(kth(ip1,k,jp1)/ps(ip1,jp1))
              end if
            end do
          end do
        end do
        istatus = nf90_put_var(ncatm, iatmvar(10), dumio, istart, icount)
        call check_ok(__FILE__,__LINE__,&
                'Error writing '//atm_variables(10)%vname//' at '//ctime, &
                'ATM FILE')
      end if
      if ( atm_variables(11)%enabled ) then
        dumio = 0.0
        do k = 1 , o_nz
          do i = 1 , o_ni
            ip1 = i+1
            do j = 1 , o_nj
              if (.not. lwrap) then
                jp1 = j+1
              else
                jp1 = j
              end if
              if (qc(ip1,k,jp1) > dlowval) then
                dumio(j,i,k) = real(kzm(ip1,k,jp1)/ps(ip1,jp1))
              end if
            end do
          end do
        end do
        istatus = nf90_put_var(ncatm, iatmvar(11), dumio, istart, icount)
        call check_ok(__FILE__,__LINE__,&
                'Error writing '//atm_variables(11)%vname//' at '//ctime, &
                'ATM FILE')
      end if
    end if

    istart(3) = iatmrec
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = o_ni
    icount(1) = o_nj

    if ( atm_variables(12)%enabled ) then
      dumio(:,:,1) = 0.0
      where (transpose(rc(o_is:o_ie,o_js:o_je)) > dlowval)
        dumio(:,:,1) = real(transpose(rc(o_is:o_ie,o_js:o_je)))
      end where
      where (transpose(rnc(o_is:o_ie,o_js:o_je)) > dlowval)
        dumio(:,:,1) = dumio(:,:,1) + real(transpose(rnc(o_is:o_ie,o_js:o_je)))
      end where
      dumio(:,:,1) = dumio(:,:,1)*real(tpd)
      istatus = nf90_put_var(ncatm, iatmvar(12), &
                             dumio(:,:,1), istart(1:3), icount(1:3))
      call check_ok(__FILE__,__LINE__,&
                    'Error writing '//atm_variables(12)%vname//' at '//ctime, &
                    'ATM FILE')
    end if

    if ( atm_variables(13)%enabled ) then
      dumio(:,:,1) = real(transpose(sum(tgb(:,o_is:o_ie,o_js:o_je), dim=1)*xns2d))
      istatus = nf90_put_var(ncatm, iatmvar(13), & 
                             dumio(:,:,1), istart(1:3), icount(1:3))
      call check_ok(__FILE__,__LINE__,&
                    'Error writing '//atm_variables(13)%vname//' at '//ctime, &
                    'ATM FILE')
    end if

    if ( atm_variables(14)%enabled ) then
      dumio(:,:,1) = 0.0
      do n = 1 , ns
        where (atmsrfmask(n,:,:) > 0)
          dumio(:,:,1) = dumio(:,:,1) + &
                           real(transpose(swt(n,o_is:o_ie,o_js:o_je)))
        end where
      end do
      where (atmsrfsum > 0)
        dumio(:,:,1) = dumio(:,:,1)/amax1(atmsrfsum/2.0,1.0)
      elsewhere
        dumio(:,:,1) = -1.E34
      end where
      istatus = nf90_put_var(ncatm, iatmvar(14), & 
                             dumio(:,:,1), istart(1:3), icount(1:3))
      call check_ok(__FILE__,__LINE__, &
                    'Error writing '//atm_variables(14)%vname//' at '//ctime, &
                    'ATM FILE')
    end if

    if ( atm_variables(15)%enabled ) then
      dumio(:,:,1) = 0.0
      do n = 1 , ns
        where (atmsrfmask(n,:,:) > 0)
          dumio(:,:,1) = dumio(:,:,1) + &
                           real(transpose(rno(n,o_is:o_ie,o_js:o_je)))
        end where
      end do
      where (atmsrfsum > 0)
        dumio(:,:,1) = dumio(:,:,1)/amax1(atmsrfsum/2.0,1.0)
      elsewhere
        dumio(:,:,1) = -1.E34
      end where
      istatus = nf90_put_var(ncatm, iatmvar(15), & 
                             dumio(:,:,1), istart(1:3), icount(1:3))
      call check_ok(__FILE__,__LINE__, &
                    'Error writing '//atm_variables(15)%vname//' at '//ctime, &
                    'ATM FILE')
    end if

    if ( debug_level > 2 ) then
      istatus = nf90_sync(ncatm)
      call check_ok(__FILE__,__LINE__,'Error sync at '//ctime, 'ATM FILE')
    end if
    iatmrec = iatmrec + 1
  end subroutine writerec_atm

  subroutine writerec_che(nx, ny, nnx, nny, nz, nt, chia, aerext, &
                          aerssa, aerasp, dtrace, wdlsc, wdcvc,   &
                          ddsfc, wxsg, wxaq, cemtrac, aertarf,    &
                          aersrrf, aertalwrf, aersrlwrf, ps,      &
                          idate)
    use netcdf
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer , intent(in) :: nx , ny , nnx , nny , nz , nt
    real(8) , dimension(iy,kz,jx,nt) , intent(in) :: chia
    real(8) , dimension(nny,nz,nnx) , intent(in) :: aerext
    real(8) , dimension(nny,nz,nnx) , intent(in) :: aerssa
    real(8) , dimension(nny,nz,nnx) , intent(in) :: aerasp
    real(8) , dimension(iy,jx,nt) , intent(in) :: dtrace
    real(8) , dimension(iy,jx,nt) , intent(in) :: wdlsc
    real(8) , dimension(iy,jx,nt) , intent(in) :: wdcvc
    real(8) , dimension(iy,jx,nt) , intent(in) :: ddsfc
    real(8) , dimension(iy,jx,nt) , intent(in) :: wxsg
    real(8) , dimension(iy,jx,nt) , intent(in) :: wxaq
    real(8) , dimension(iy,jx,nt) , intent(in) :: cemtrac
    real(8) , dimension(iy,jx) , intent(in) :: ps
    real(8) , dimension(nny,nnx) , intent(in) :: aertarf
    real(8) , dimension(nny,nnx) , intent(in) :: aersrrf
    real(8) , dimension(nny,nnx) , intent(in) :: aertalwrf
    real(8) , dimension(nny,nnx) , intent(in) :: aersrlwrf
    integer :: n , k
    integer , dimension(5) :: istart , icount
    real(8) , dimension(1) :: nctime
    type(rcm_time_interval) :: tdif
    character(len=36) :: ctime

    if (nx < o_nj .or. ny < o_ni .or. nz > o_nz) then
      write (6,*) 'Error writing record on CHE file'
      write (6,*) 'Expecting layers ', o_nz, 'x', o_nj, 'x', o_ni
      write (6,*) 'Got layers       ', nz, 'x', nx, 'x', ny
      call fatal(__FILE__,__LINE__,'DIMENSION MISMATCH')
    end if

    ctime = tochar(idate)

    istart(1) = icherec
    icount(1) = 1
    tdif = idate-cordex_refdate
    nctime(1) = tohours(tdif)
    istatus = nf90_put_var(ncche, ichevar(1), nctime, istart(1:1), icount(1:1))
    call check_ok(__FILE__,__LINE__,'Error writing itime '//ctime, 'CHE FILE')

    dumio(:,:,1) = real(transpose(ps(o_is:o_ie,o_js:o_je)+ptop)*d_10)
    istart(3) = icherec
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = o_ni
    icount(1) = o_nj
    istatus = nf90_put_var(ncche, ichevar(2), &
                           dumio(:,:,1), istart(1:3), icount(1:3))
    call check_ok(__FILE__,__LINE__,'Error writing ps at '//ctime, 'CHE FILE')

    if ( che_variables(3)%enabled ) then
      do n =  1 , nt
        istart(5) = icherec
        istart(4) = n
        istart(3) = 1
        istart(2) = 1
        istart(1) = 1
        icount(5) = 1
        icount(4) = 1
        icount(3) = o_nz
        icount(2) = o_ni
        icount(1) = o_nj
        do k = 1 , nz
          dumio(:,:,k) = real(transpose(chia(o_is:o_ie,k,o_js:o_je,n) / &
                                     ps(o_is:o_ie,o_js:o_je)))
        end do
        istatus = nf90_put_var(ncche, ichevar(3), dumio, istart, icount)
        call check_ok(__FILE__,__LINE__, &
                      'Error writing '//che_variables(3)%vname//' at '//ctime, &
                      'CHE FILE')
      end do
    end if

    istart(4) = icherec
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = o_nz
    icount(2) = o_ni
    icount(1) = o_nj

    if ( che_variables(4)%enabled ) then
      do k = 1 , nz
        dumio(:,:,k) = real(transpose(aerext(o_is:o_ie,k,o_js:o_je)))
      end do
      istatus = nf90_put_var(ncche, ichevar(4), dumio, istart(1:4), icount(1:4))
      call check_ok(__FILE__,__LINE__, &
                    'Error writing '//che_variables(4)%vname//' at '//ctime, &
                    'CHE FILE')
    end if

    if ( che_variables(5)%enabled ) then
      do k = 1 , nz
        dumio(:,:,k) = real(transpose(aerssa(o_is:o_ie,k,o_js:o_je)))
      end do
      istatus = nf90_put_var(ncche, ichevar(5), dumio, istart(1:4), icount(1:4))
      call check_ok(__FILE__,__LINE__, &
                    'Error writing '//che_variables(5)%vname//' at '//ctime, &
                    'CHE FILE')
    end if

    if ( che_variables(6)%enabled ) then
      do k = 1 , nz
        dumio(:,:,k) = real(transpose(aerasp(o_is:o_ie,k,o_js:o_je)))
      end do
      istatus = nf90_put_var(ncche, ichevar(6), dumio, istart(1:4), icount(1:4))
      call check_ok(__FILE__,__LINE__, &
                    'Error writing '//che_variables(6)%vname//' at '//ctime, &
                    'CHE FILE')
    end if

    do n = 1 , nt
      istart(4) = icherec
      istart(3) = n
      istart(2) = 1
      istart(1) = 1
      icount(4) = 1
      icount(3) = 1
      icount(2) = o_ni
      icount(1) = o_nj
      if ( che_variables(7)%enabled ) then
        dumio(:,:,1) = real(transpose(dtrace(o_is:o_ie,o_js:o_je,n)))
        istatus = nf90_put_var(ncche, ichevar(7), &
                               dumio(:,:,1), istart(1:4), icount(1:4))
        call check_ok(__FILE__,__LINE__, &
                      'Error writing '//che_variables(7)%vname//' at '//ctime, &
                      'CHE FILE')
      end if
      if ( che_variables(8)%enabled ) then
        dumio(:,:,1) = real(transpose(wdlsc(o_is:o_ie,o_js:o_je,n))*cfd)
        istatus = nf90_put_var(ncche, ichevar(8), &
                               dumio(:,:,1), istart(1:4), icount(1:4))
        call check_ok(__FILE__,__LINE__, &
                      'Error writing '//che_variables(8)%vname//' at '//ctime, &
                      'CHE FILE')
      end if
      if ( che_variables(9)%enabled ) then
        dumio(:,:,1) = real(transpose(wdcvc(o_is:o_ie,o_js:o_je,n))*cfd)
        istatus = nf90_put_var(ncche, ichevar(9), &
                               dumio(:,:,1), istart(1:4), icount(1:4))
        call check_ok(__FILE__,__LINE__, &
                      'Error writing '//che_variables(9)%vname//' at '//ctime, &
                      'CHE FILE')
      end if
      if ( che_variables(10)%enabled ) then
        dumio(:,:,1) = real(transpose(ddsfc(o_is:o_ie,o_js:o_je,n))*cfd)
        istatus = nf90_put_var(ncche, ichevar(10), &
                               dumio(:,:,1), istart(1:4), icount(1:4))
        call check_ok(__FILE__,__LINE__, &
                      'Error writing '//che_variables(10)%vname// &
                      ' at '//ctime, 'CHE FILE')
      end if
      if ( che_variables(11)%enabled ) then
        dumio(:,:,1) = real(transpose(wxsg(o_is:o_ie,o_js:o_je,n))*cfd)
        istatus = nf90_put_var(ncche, ichevar(11), &
                               dumio(:,:,1), istart(1:4), icount(1:4))
        call check_ok(__FILE__,__LINE__, &
                      'Error writing '//che_variables(11)%vname// &
                      ' at '//ctime, 'CHE FILE')
      end if
      if ( che_variables(12)%enabled ) then
        dumio(:,:,1) = real(transpose(wxaq(o_is:o_ie,o_js:o_je,n))*cfd)
        istatus = nf90_put_var(ncche, ichevar(12), &
                               dumio(:,:,1), istart(1:4), icount(1:4))
        call check_ok(__FILE__,__LINE__, &
                      'Error writing '//che_variables(12)%vname// &
                      ' at '//ctime, 'CHE FILE')
      end if
      if ( che_variables(13)%enabled ) then
        dumio(:,:,1) = real(transpose(cemtrac(o_is:o_ie,o_js:o_je,n))*cfd)
        istatus = nf90_put_var(ncche, ichevar(13), &
                               dumio(:,:,1), istart(1:4), icount(1:4))
        call check_ok(__FILE__,__LINE__, &
                      'Error writing '//che_variables(13)%vname// &
                      ' at '//ctime, 'CHE FILE')
      end if
    end do

    istart(3) = icherec
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = o_ni
    icount(1) = o_nj

    if ( che_variables(14)%enabled ) then
      dumio(:,:,1) = real(transpose(aertarf(o_is:o_ie,o_js:o_je)))
      istatus = nf90_put_var(ncche, ichevar(14), &
                             dumio(:,:,1), istart(1:3), icount(1:3))
      call check_ok(__FILE__,__LINE__, &
                    'Error writing '//che_variables(14)%vname//' at '//ctime, &
                    'CHE FILE')
    end if
    if ( che_variables(15)%enabled ) then
      dumio(:,:,1) = real(transpose(aersrrf(o_is:o_ie,o_js:o_je)))
      istatus = nf90_put_var(ncche, ichevar(15), &
                             dumio(:,:,1), istart(1:3), icount(1:3))
      call check_ok(__FILE__,__LINE__, &
                    'Error writing '//che_variables(15)%vname//' at '//ctime, &
                    'CHE FILE')
    end if
    if ( che_variables(16)%enabled ) then
      dumio(:,:,1) = real(transpose(aertalwrf(o_is:o_ie,o_js:o_je)))
      istatus = nf90_put_var(ncche, ichevar(16), &
                             dumio(:,:,1), istart(1:3), icount(1:3))
      call check_ok(__FILE__,__LINE__, &
                    'Error writing '//che_variables(16)%vname//' at '//ctime, &
                    'CHE FILE')
    end if
    if ( che_variables(17)%enabled ) then
      dumio(:,:,1) = real(transpose(aersrlwrf(o_is:o_ie,o_js:o_je)))
      istatus = nf90_put_var(ncche, ichevar(17), &
                             dumio(:,:,1), istart(1:3), icount(1:3))
      call check_ok(__FILE__,__LINE__, &
                    'Error writing '//che_variables(17)%vname//' at '//ctime, &
                    'CHE FILE')
    end if

    if ( debug_level > 2 ) then
      istatus = nf90_sync(ncche)
      call check_ok(__FILE__,__LINE__,'Error sync at '//ctime, 'CHE FILE')
    end if
    icherec = icherec + 1
  end subroutine writerec_che

  subroutine writerec_lak(nx,ny,numbat,fbat,evl,aveice,hsnow,tlake,idate)
    use netcdf
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer , intent(in) :: nx , ny , numbat
    real(4) , dimension(nx,ny,numbat) , intent(in) :: fbat
    real(8) , dimension(nnsg,iym1,jx) , intent(in) :: evl
    real(8) , dimension(nnsg,iym1,jx) , intent(in) :: aveice
    real(8) , dimension(nnsg,iym1,jx) , intent(in) :: hsnow
    real(8) , dimension(ndpmax,nnsg,iym1,jx) , intent(in) :: tlake
    integer :: ivar
    integer :: n
    integer , dimension(4) :: istart , icount
    real(8) , dimension(1) :: nctime
    type(rcm_time_interval) :: tdif
    character(len=36) :: ctime

    if (nx /= o_nj .or. ny /= o_ni) then
      write (6,*) 'Error writing record on LAK file'
      write (6,*) 'Expecting layers ', o_nj, 'x', o_ni
      write (6,*) 'Got layers       ', nx, 'x', ny
      call fatal(__FILE__,__LINE__,'DIMENSION MISMATCH')
    end if

    ctime = tochar(idate)

    istart(1) = ilakrec
    icount(1) = 1
    tdif = idate-cordex_refdate
    nctime(1) = tohours(tdif)
    istatus = nf90_put_var(nclak, ilakvar(1), nctime, istart(1:1), icount(1:1))
    call check_ok(__FILE__,__LINE__,'Error writing itime '//ctime, 'LAK FILE')

    ivar = 2
    do n = 1 , numbat
      if (lak_fbats(n) == 0) cycle
      istart(3) = ilakrec
      istart(2) = 1
      istart(1) = 1
      icount(3) = 1
      icount(2) = o_ni
      icount(1) = o_nj
      istatus = nf90_put_var(nclak, ilakvar(ivar), & 
                      fbat(:,:,n), istart(1:3), icount(1:3))
      call check_ok(__FILE__,__LINE__, &
                    'Error writing '//lak_variables(ivar)%vname// &
                    ' at '//ctime, 'LAK FILE')
      ivar = ivar + 1
    end do

    ! Add lake model output
    dumio(:,:,1) =  real(transpose(sum(evl(:,o_is:o_ie,o_js:o_je),1))*xns2d)
    istatus = nf90_put_var(nclak, ilakvar(ivar), & 
                    dumio(:,:,1), istart(1:3), icount(1:3))
    call check_ok(__FILE__,__LINE__, &
                  'Error writing '//lak_variables(ivar)%vname// &
                  ' at '//ctime, 'LAK FILE')
    ivar = ivar + 1
    dumio(:,:,1) = real(transpose(sum(aveice(:,o_is:o_ie,o_js:o_je),1))*xns2d)
    istatus = nf90_put_var(nclak, ilakvar(ivar), & 
                    dumio(:,:,1), istart(1:3), icount(1:3))
    call check_ok(__FILE__,__LINE__, &
                  'Error writing '//lak_variables(ivar)%vname// &
                  ' at '//ctime, 'LAK FILE')
    ivar = ivar + 1
    dumio(:,:,1) = real(transpose(sum(hsnow(:,o_is:o_ie,o_js:o_je),1))*xns2d)
    istatus = nf90_put_var(nclak, ilakvar(ivar), & 
             dumio(:,:,1), istart(1:3), icount(1:3))
    call check_ok(__FILE__,__LINE__, &
                  'Error writing '//lak_variables(ivar)%vname// &
                  ' at '//ctime, 'LAK FILE')
    ivar = ivar + 1
    do n = 1 , ndpmax
      istart(4) = ilakrec
      istart(3) = n
      istart(2) = 1
      istart(1) = 1
      icount(4) = 1
      icount(3) = 1
      icount(2) = o_ni
      icount(1) = o_nj
      dumio(:,:,1) = real(transpose( &
                  sum(tlake(n,:,o_is:o_ie,o_js:o_je),1))*xns2d)
      where (iolnds == 14)
        dumio(:,:,1) = dumio(:,:,1) + real(tzero)
      end where
      istatus = nf90_put_var(nclak, ilakvar(ivar), dumio(:,:,1), istart, icount)
      call check_ok(__FILE__,__LINE__, &
                    'Error writing '//lak_variables(ivar)%vname// &
                  ' at '//ctime, 'LAK FILE')
    end do

    if ( debug_level > 2 ) then
      istatus = nf90_sync(nclak)
      call check_ok(__FILE__,__LINE__,'Error sync at '//ctime, 'LAK FILE')
    end if
    ilakrec = ilakrec + 1
  end subroutine writerec_lak

  subroutine check_ok(f,l,m1,mf)
    use netcdf
    implicit none
    character(*) , intent(in) :: f, m1 , mf
    integer , intent(in) :: l
    if (istatus /= nf90_noerr) then 
      write (6,*) trim(m1)
      write (6,*) nf90_strerror(istatus)
      call fatal(f,l,trim(mf))
    end if
  end subroutine check_ok

  subroutine release_mod_ncio
    implicit none
    call close_domain
    call close_icbc
    call close_common(ncatm,'ATM')
    call close_common(ncsrf,'SRF')
    call close_common(ncsub,'SUB')
    call close_common(ncrad,'RAD')
    call close_common(ncche,'CHE')
    call close_common(nclak,'LAK')
  end subroutine release_mod_ncio

end module mod_ncio
