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
        use mod_message
        use mod_date
        use mod_trachem
!
        private
!
        public :: init_mod_ncio , release_mod_ncio
        public :: open_domain , read_domain , read_domain_lake, &
                  read_subdomain , read_subdomain_lake,         &
                  read_texture , close_domain
        public :: open_icbc , read_icbc , icbc_search
        public :: read_aerosol
        public :: prepare_common_out
        public :: writerec_atm , writerec_srf , writerec_sub , &
                  writerec_rad , writerec_che , writerec_lak
!
        integer , parameter :: n_atmvar = 12
        integer , parameter :: n_srfvar = 31
        integer , parameter :: n_subvar = 16
        integer , parameter :: n_radvar = 15
        integer , parameter :: n_chevar = 17
        integer , parameter :: n_lakvar = 16

        integer :: idmin , isdmin , ibcin , ncatm , ncsrf , &
                   ncsub , ncrad , ncche , nclak
        integer :: istatus
        integer :: ibcrec , ibcnrec
        integer :: iatmrec , isrfrec , isubrec , iradrec , &
                   icherec , ilakrec
        integer , dimension(n_atmvar) :: iatmvar
        integer , dimension(n_srfvar) :: isrfvar
        integer , dimension(n_subvar) :: isubvar
        integer , dimension(n_radvar) :: iradvar
        integer , dimension(n_chevar) :: ichevar
        integer , dimension(n_lakvar) :: ilakvar
        character(256) :: dname , sdname , aername , icbcname
        integer , dimension(:) , allocatable :: icbc_idate
        real(4) , dimension(:) , allocatable :: hsigma
        integer , dimension(7) :: icbc_ivar
        logical :: lso4p
        real(8) :: rpt, tpd, cfd
        real(8) :: xns2r

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

        real(4) , dimension(:,:) , allocatable :: ioxlat
        real(4) , dimension(:,:) , allocatable :: ioxlon
        real(4) , dimension(:,:) , allocatable :: iotopo
        real(4) , dimension(:,:) , allocatable :: iomask
        real(4) , dimension(:,:) , allocatable :: iolnds
        real(4) , dimension(:,:) , allocatable :: ioxlat_s
        real(4) , dimension(:,:) , allocatable :: ioxlon_s
        real(4) , dimension(:,:) , allocatable :: iotopo_s
        real(4) , dimension(:,:) , allocatable :: iomask_s
        real(4) , dimension(:,:) , allocatable :: subio
        real(4) , dimension(:,:,:) , allocatable :: dumio
        real(4) , dimension(:,:) , allocatable :: sp2d
        real(4) , dimension(:,:) , allocatable :: sp2d1
        real(4) , dimension(:,:,:) , allocatable :: atmsrfmask
        real(4) , dimension(:,:) , allocatable :: atmsrfsum

        character(len=8) , dimension(n_atmvar) :: atm_names
        character(len=8) , dimension(n_srfvar) :: srf_names
        character(len=8) , dimension(n_subvar) :: sub_names
        character(len=8) , dimension(n_radvar) :: rad_names
        character(len=8) , dimension(n_chevar) :: che_names
        character(len=8) , dimension(n_lakvar) :: lak_names

        integer , dimension(numbat) :: lak_fbats

        character(64) :: cdum

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
        data atm_names / 'time', 'ps', 'u', 'v', 'omega', 't',         &
                         'qv', 'qc', 'tpr', 'tgb', 'swt', 'rno' /
        data srf_names / 'time', 'tbnds', 'ps', 'u10m', 'v10m',        &
                         'uvdrag', 'tg', 'tlef', 't2m', 'q2m', 'smw',  &
                         'tpr', 'evp', 'runoff', 'scv', 'sena', 'flw', &
                         'fsw', 'flwd', 'sina', 'prcv', 'zpbl',        &
                         'tgmax', 'tgmin', 't2max', 't2min', 'w10max', &
                         'ps_min' , 'aldirs' , 'aldifs' , 'seaice' /
        data sub_names / 'time', 'ps', 'u10m', 'v10m', 'uvdrag', 'tg', &
                         'tlef', 't2m' , 'q2m' , 'smw', 'tpr' , 'evp', &
                         'runoff', 'scv', 'sena', 'prcv' /
        data rad_names / 'time', 'ps', 'cld', 'clwp', 'qrs', 'qrl',    &
                         'frsa', 'frla', 'clrst', 'clrss', 'clrlt',    &
                         'clrls', 'solin', 'sabtp', 'firtp' /
        data che_names / 'time', 'ps', 'trac', 'aext8', 'assa8',       &
                         'agfu8', 'colb', 'wdlsc', 'wdcvc', 'sdrdp',   &
                         'xgasc', 'xaquc', 'emiss', 'acstoarf',        &
                         'acstsrrf', 'acstalrf', 'acssrlrf' /
        data lak_names / 'time', 'ps', 'tg' , 'tpr' , 'scv', 'sena' ,  &
                         'flw', 'fsw', 'fld' , 'sina' , 'aldirs' ,     &
                         'aldifs' , 'evl' , 'aveice' , 'hsnow',        &
                         'tlake' /

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
            case(98)
             write (cdum,'(a)') &
               'Grell over ocean, Emanuel (1991) over land'
            case(99)
             write (cdum,'(a)') &
               'Emanuel (1991) over ocean, Grell over land'
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
            case default 
             write (cdum,'(a)') 'Unknown or not specified'
          end select
        end subroutine cdumpbl

        subroutine cdummoist
          implicit none
          select case (ipptls)
            case(1)
             write (cdum,'(a)') &
                'Explicit moisture (SUBEX; Pal et al 2000)'
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
             write (cdum,'(a)') &
               'Hydrostatic deduction with perturbation temperature'
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
              if (sname == atm_names(i)) then
                ivarname_lookup = i
                exit
              end if
            end do
          else if (ctype == 'SRF') then
            do i = 1 , n_srfvar
              if (sname == srf_names(i)) then
                ivarname_lookup = i
                exit
              end if
            end do
          else if (ctype == 'SUB') then
            do i = 1 , n_subvar
              if (sname == sub_names(i)) then
                ivarname_lookup = i
                exit
              end if
            end do
          else if (ctype == 'RAD') then
            do i = 1 , n_radvar
              if (sname == rad_names(i)) then
                ivarname_lookup = i
                exit
              end if
            end do
          else if (ctype == 'CHE') then
            do i = 1 , n_chevar
              if (sname == che_names(i)) then
                ivarname_lookup = i
                exit
              end if
            end do
          else if (ctype == 'LAK') then
            do i = 1 , n_lakvar
              if (sname == lak_names(i)) then
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
          sdname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN'// &
              &    sbstring//'.nc'
          aername = trim(dirglob)//pthsep//trim(domname)//'_AERO.nc'
          icbcname = trim(dirglob)//pthsep//trim(domname)//'_ICBC.'// &
              &      'YYYYMMDDHH.nc'

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
          xns2r = 1.0D0/dble(nnsg)
          allocate(hsigma(o_nz))
          allocate(ioxlat(o_nj,o_ni))
          allocate(ioxlon(o_nj,o_ni))
          allocate(iotopo(o_nj,o_ni))
          allocate(iomask(o_nj,o_ni))
          allocate(iolnds(o_nj,o_ni))
          allocate(dumio(o_nj,o_ni,o_nz))
          allocate(sp2d(jx,iy))
          allocate(atmsrfmask(nnsg,o_nj,o_ni))
          allocate(atmsrfsum(o_nj,o_ni))
          hsigma(:) = 0.0
          ioxlat(:,:) = 0.0
          ioxlon(:,:) = 0.0
          iotopo(:,:) = 0.0
          iomask(:,:) = 0.0
          iolnds(:,:) = 0.0
          dumio(:,:,:) = 0.0
          sp2d(:,:) = 0.0
          atmsrfmask(:,:,:) = 0.0
          atmsrfsum(:,:) = 0.0
          if (nsg > 1) then
            allocate(ioxlat_s(o_njg,o_nig))
            allocate(ioxlon_s(o_njg,o_nig))
            allocate(iotopo_s(o_njg,o_nig))
            allocate(iomask_s(o_njg,o_nig))
            allocate(subio(o_njg,o_nig))
            allocate(sp2d1(jxsg,iysg))
            ioxlat_s(:,:) = 0.0
            ioxlon_s(:,:) = 0.0
            iotopo_s(:,:) = 0.0
            iotopo_s(:,:) = 0.0
            iomask_s(:,:) = 0.0
            subio(:,:) = 0.0
            sp2d1(:,:) = 0.0
          end if
        end subroutine init_mod_ncio

        subroutine open_domain(r8pt , dx , sigma)
          use netcdf
          implicit none

          real(8) , intent(out) :: dx
          real(8) , intent(out) :: r8pt
          real(8) , dimension(kzp1) :: sigma

          integer :: ivarid , idimid
          integer :: iyy , jxx , kzz , k
          character(6) :: proj
          real(4) :: dsx , iclat , iclon , ptsp
          real(4) , dimension(kzp1) :: rsdum

          write (aline,*) 'open_domain: READING HEADER FILE:', dname
          call say
          istatus = nf90_open(dname, nf90_nowrite, idmin)
          call check_ok('Error Opening Domain file '//trim(dname), &
                        'CANNOT OPEN DOMAIN FILE')
          if ( nsg > 1 ) then
            write (aline,*) 'READING HEADER SUBDOMAIN FILE:', sdname
            call say
            istatus = nf90_open(sdname, nf90_nowrite, isdmin)
            call check_ok('Error Opening SubDomain file '//trim(sdname),&
                          'CANNOT OPEN SUBDOM FILE')
          end if
          istatus = nf90_inq_dimid(idmin, 'iy', idimid)
          call check_ok('Dimension iy missing', 'DOMAIN FILE ERROR')
          istatus = nf90_inquire_dimension(idmin, idimid, len=iyy)
          call check_ok('Dimension iy read error', 'DOMAIN FILE ERROR')
          istatus = nf90_inq_dimid(idmin, 'jx', idimid)
          call check_ok('Dimension jx missing', 'DOMAIN FILE ERROR')
          istatus = nf90_inquire_dimension(idmin, idimid, len=jxx)
          call check_ok('Dimension jx read error', 'DOMAIN FILE ERROR')
          istatus = nf90_inq_dimid(idmin, 'kz', idimid)
          call check_ok('Dimension kz missing', 'DOMAIN FILE ERROR')
          istatus = nf90_inquire_dimension(idmin, idimid, len=kzz)
          call check_ok('Dimension kz read error', 'DOMAIN FILE ERROR')
          istatus = nf90_inq_varid(idmin, 'ptop', ivarid)
          call check_ok('Variable ptop missing', 'DOMAIN FILE ERROR')
          istatus = nf90_get_var(idmin, ivarid, ptsp)
          call check_ok('Variable ptop read error', 'DOMAIN FILE ERROR')
          istatus = nf90_get_att(idmin, nf90_global, 'projection', proj)
          call check_ok('Attribute projection missing', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_get_att(idmin, nf90_global, &
                       &         'grid_size_in_meters', dsx)
          call check_ok('Attribute gridsize missing', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_get_att(idmin, nf90_global, &
                     &         'latitude_of_projection_origin', iclat)
          call check_ok('Attribute clat missing', 'DOMAIN FILE ERROR')
          istatus = nf90_get_att(idmin, nf90_global, &
                      &         'longitude_of_projection_origin', iclon)
          call check_ok('Attribute clon missing', 'DOMAIN FILE ERROR')
!
!         Consistency Check
!
          if ( iyy /= iy .or. jxx /= jx .or. kzz /= kzp1 ) then
            write (6,*) 'Error: dims from regcm.in and DOMAIN file ', &
                        'differ.'
            write (aline,*) 'Input namelist : IY=' , iy , &
                   & '  JX=' ,  jx , '  KZ=' , kz
            call say
            write (aline,*) 'DOMAIN file    : IY=' , iyy , &
                   & '  JX=' ,  jxx , '  KZ=' , kzz-1
            call say
            call fatal(__FILE__,__LINE__,'DIMENSION MISMATCH')
          end if
          if (dabs(dble(ptsp*d_r10)-dble(ptop)) > 0.001D+00) then
            write (6,*) 'Error: ptop from regcm.in and DOMAIN file ', &
                        'differ.'
            write (6,*) 'Input namelist = ', ptop
            write (6,*) 'DOMAIN file    = ', ptsp*d_r10
            call fatal(__FILE__,__LINE__, 'DOMAIN ptop ERROR')
          end if
          if (proj /= iproj) then
            write (6,*) 'Error: proj from regcm.in and DOMAIN file ', &
                        'differ.'
            write (6,*) 'Input namelist = ', iproj
            write (6,*) 'DOMAIN file    = ', proj
            call fatal(__FILE__,__LINE__, 'DOMAIN proj ERROR')
          end if
          if (dabs(dble(dsx*d_r1000)-dble(ds)) > 0.001D+00) then
            write (6,*) 'Error: ds from regcm.in and DOMAIN file ', &
                        'differ.'
            write (6,*) 'Input namelist = ', ds
            write (6,*) 'DOMAIN file    = ', dsx*d_r1000
            call fatal(__FILE__,__LINE__, 'DOMAIN ds ERROR')
          end if
          if (dabs(dble(iclat)-dble(clat)) > 0.001D+00) then
            write (6,*) 'Error: clat from regcm.in and DOMAIN file ', &
                        'differ.'
            write (6,*) 'Input namelist = ', clat
            write (6,*) 'DOMAIN file    = ', iclat
            call fatal(__FILE__,__LINE__, 'DOMAIN clat ERROR')
          end if
          if (dabs(dble(iclon)-dble(clon)) > 0.001D+00) then
            write (6,*) 'Error: clon from regcm.in and DOMAIN file ', &
                        'differ.'
            write (6,*) 'Input namelist = ', clon
            write (6,*) 'DOMAIN file    = ', iclon
            call fatal(__FILE__,__LINE__, 'DOMAIN clon ERROR')
          end if
!
!         Assign values in the top data modules
!
          r8pt = dble(ptsp)*d_r10
          rpt = ptop
          tpd = houpd/tapfrq
          cfd = houpd/chemfrq
          dx = dble(dsx)
          istatus = nf90_inq_varid(idmin, 'sigma', ivarid)
          call check_ok('Variable sigma missing', 'DOMAIN FILE ERROR')
          istatus = nf90_get_var(idmin, ivarid, rsdum)
          call check_ok('Variable sigma read error','DOMAIN FILE ERROR')
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
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if

          istatus = nf90_inq_varid(idmin, 'topo', ivarid)
          call check_ok('Variable topo missing', 'DOMAIN FILE ERROR')
          istatus = nf90_get_var(idmin, ivarid, sp2d)
          call check_ok('Variable topo read error', 'DOMAIN FILE ERROR')
          ht = dble(transpose(sp2d))
          iotopo = sp2d(o_js:o_je,o_is:o_ie)
          istatus = nf90_inq_varid(idmin, 'landuse', ivarid)
          call check_ok('Variable landuse missing','DOMAIN FILE ERROR')
          istatus = nf90_get_var(idmin, ivarid, sp2d)
          call check_ok('Variable landuse read error', &
                        'DOMAIN FILE ERROR')
          lnd = dble(transpose(sp2d))
          iolnds = sp2d(o_js:o_je,o_is:o_ie)
          istatus = nf90_inq_varid(idmin, 'xlat', ivarid)
          call check_ok('Variable xlat missing', 'DOMAIN FILE ERROR')
          istatus = nf90_get_var(idmin, ivarid, sp2d)
          call check_ok('Variable xlat read error', 'DOMAIN FILE ERROR')
          xlat = dble(transpose(sp2d))
          ioxlat = sp2d(o_js:o_je,o_is:o_ie)
          istatus = nf90_inq_varid(idmin, 'xlon', ivarid)
          call check_ok('Variable xlon missing', 'DOMAIN FILE ERROR')
          istatus = nf90_get_var(idmin, ivarid, sp2d)
          call check_ok('Variable xlon read error', 'DOMAIN FILE ERROR')
          xlon = dble(transpose(sp2d))
          ioxlon = sp2d(o_js:o_je,o_is:o_ie)
          istatus = nf90_inq_varid(idmin, 'xmap', ivarid)
          call check_ok('Variable xmap missing', 'DOMAIN FILE ERROR')
          istatus = nf90_get_var(idmin, ivarid, sp2d)
          call check_ok('Variable xmap read error','DOMAIN FILE ERROR')
          xmap = dble(transpose(sp2d))
          istatus = nf90_inq_varid(idmin, 'dmap', ivarid)
          call check_ok('Variable dmap missing','DOMAIN FILE ERROR')
          istatus = nf90_get_var(idmin, ivarid, sp2d)
          call check_ok('Variable dmap read error', 'DOMAIN FILE ERROR')
          dmap = dble(transpose(sp2d))
          istatus = nf90_inq_varid(idmin, 'coriol', ivarid)
          call check_ok('Variable coriol missing', 'DOMAIN FILE ERROR')
          istatus = nf90_get_var(idmin, ivarid, sp2d)
          call check_ok('Variable coriol read error', &
                        'DOMAIN FILE ERROR')
          f = dble(transpose(sp2d))
          istatus = nf90_inq_varid(idmin, 'mask', ivarid)
          call check_ok('Variable mask missing', 'DOMAIN FILE ERROR')
          istatus = nf90_get_var(idmin, ivarid, sp2d)
          call check_ok('Variable mask read error', &
                        'DOMAIN FILE ERROR')
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
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if

          istatus = nf90_inq_varid(idmin, 'dhlake', ivarid)
          call check_ok('Variable dhlake missing', 'DOMAIN FILE ERROR')
          istatus = nf90_get_var(idmin, ivarid, sp2d)
          call check_ok('Variable dhlake read error', &
                        'DOMAIN FILE ERROR')
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
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if

          istatus = nf90_inq_varid(isdmin, 'topo', ivarid)
          call check_ok('Variable topo missing', 'SUBDOMAIN FILE ERROR')
          istatus = nf90_get_var(isdmin, ivarid, sp2d1)
          call check_ok('Variable topo read error', &
                        'SUBDOMAIN FILE ERROR')
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
          call check_ok('Variable landuse missing', &
                        'SUBDOMAIN FILE ERROR')
          istatus = nf90_get_var(isdmin, ivarid, sp2d1)
          call check_ok('Variable landuse read error', &
                        'SUBDOMAIN FILE ERROR')
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
          call check_ok('Variable xlat missing','SUBDOMAIN FILE ERROR')
          istatus = nf90_get_var(isdmin, ivarid, sp2d1)
          call check_ok('Variable xlat read error', &
                        'SUBDOMAIN FILE ERROR')
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
          call check_ok('Variable xlon missing','SUBDOMAIN FILE ERROR')
          istatus = nf90_get_var(isdmin, ivarid, sp2d1)
          call check_ok('Variable xlon read error', &
                        'SUBDOMAIN FILE ERROR')
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
          call check_ok('Variable mask missing', &
                        'SUBDOMAIN FILE ERROR')
          istatus = nf90_get_var(isdmin, ivarid, sp2d1)
          call check_ok('Variable mask read error', &
                        'SUBDOMAIN FILE ERROR')
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
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if

          istatus = nf90_inq_varid(isdmin, 'dhlake', ivarid)
          call check_ok('Variable dhlake missing', &
                        'SUBDOMAIN FILE ERROR')
          istatus = nf90_get_var(isdmin, ivarid, sp2d1)
          call check_ok('Variable dhlake read error', &
                        'SUBDOMAIN FILE ERROR')
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
            call check_ok('Domain file close error','DOMAIN FILE ERROR')
            idmin = -1
          end if
          if ( nsg>1 .and. isdmin >=0 ) then
            istatus = nf90_close(isdmin)
            call check_ok('SubDomain file close error', &
                        'SUBDOMAIN FILE ERROR')
            isdmin = -1
          end if

        end subroutine close_domain

        subroutine read_texture(nats,texture)
          use netcdf
          implicit none

          integer , intent(in) :: nats
          real(8) , dimension(iy,jx,nats) , intent(out) :: texture

          integer :: ivarid
          integer :: i , j , n
          integer , dimension(3) :: istart , icount
          real(4), dimension(jx,iy) ::  toto

          if (idmin < 0) then
            istatus = nf90_open(dname, nf90_nowrite, idmin)
            call check_ok('Error Opening Domain file '//trim(dname), &
                       &  'DOMAIN FILE OPEN ERROR')
          end if
          istatus = nf90_inq_varid(idmin, 'texture_fraction', ivarid)
          call check_ok('Variable texture_fraction missing', &
                     &  'DOMAIN FILE ERROR')
          istart(2) = 1
          istart(1) = 1
          icount(3) = 1
          icount(2) = iy
          icount(1) = jx
          do n = 1 , nats
            istart(3) = n
            istatus = nf90_get_var(idmin, ivarid, toto, istart, icount)
            call check_ok('Variable texture_frac read error', &
                     &  'DOMAIN FILE ERROR')
            do j = 1 , jx
              do i = 1 , iy
                texture(i,j,n) = dble(toto(j,i))*0.01D0
                if (texture(i,j,n)<d_zero) texture(i,j,n)=d_zero
              end do
            end do
          end do
          call close_domain

        end subroutine read_texture

        subroutine read_aerosol(chtrname,chemsrc)
          use netcdf
          implicit none

          character(5) , dimension(ntr) , intent(in) :: chtrname
          real(8) , dimension(iy,jx,12,ntr) , intent(out) :: chemsrc

          integer :: ncid , ivarid
          real(4) , dimension(jx,iy) :: toto
          character(5) :: aerctl
          integer , dimension(3) :: istart , icount
          integer :: itr , i , j , m

          istatus = nf90_open(aername, nf90_nowrite, ncid)
          call check_ok('Error Opening Aerosol file '//trim(aername), &
                     &  'AEROSOL FILE OPEN ERROR')

          do itr = 1 , ntr
            aerctl = chtrname(itr)
            write (aline, *) itr , aerctl
            call say
            if ( aerctl(1:4) /= 'DUST') then
              if ( aerctl(1:3) == 'SO2' ) then
                if ( aertyp(4:4) == '1' ) then
                  istatus = nf90_inq_varid(ncid, 'so2', ivarid)
                  call check_ok('Variable so2 missing', &
                            &   'AEROSOL FILE ERROR')
                  istatus = nf90_get_var(ncid, ivarid, toto)
                  call check_ok('Variable so2 read error', &
                            &   'AEROSOL FILE ERROR')
                  do m = 1 , 12
                    do j = 1 , jx
                      do i = 1 , iy
                        chemsrc(i,j,m,itr) = dble(toto(j,i))
                      end do
                    end do
                  end do
                end if
                if ( aertyp(5:5) == '1' ) then
                  istatus = nf90_inq_varid(ncid, 'so2_monthly', ivarid)
                  call check_ok('Variable so2_mon missing', &
                            &   'AEROSOL FILE ERROR')
                  istart(1) = 1
                  istart(2) = 1
                  icount(1) = jx
                  icount(2) = iy
                  icount(3) = 1
                  do m = 1 , 12
                    istart(3) = m
                    istatus = nf90_get_var(ncid,ivarid,toto, &
                                        &  istart,icount)
                    call check_ok('Variable so2_mon read err', &
                            &   'AEROSOL FILE ERROR')
                    do j = 1 , jx
                      do i = 1 , iy
                        chemsrc(i,j,m,itr) = chemsrc(i,j,m,itr) + &
                                              & dble(toto(j,i))
                      end do
                    end do
                  end do
                end if
              else if ( aerctl(1:2) == 'BC' ) then
                if ( aertyp(4:4) == '1' ) then
                  istatus = nf90_inq_varid(ncid, 'bc', ivarid)
                  call check_ok('Variable bc missing', &
                            &   'AEROSOL FILE ERROR')
                  istatus = nf90_get_var(ncid, ivarid, toto)
                  call check_ok('Variable bc read error', &
                            &   'AEROSOL FILE ERROR')
                  do m = 1 , 12
                    do j = 1 , jx
                      do i = 1 , iy
                        chemsrc(i,j,m,itr) = dble(toto(j,i))
                      end do
                    end do
                  end do
                end if
                if ( aertyp(5:5) == '1' ) then
                  istatus = nf90_inq_varid(ncid, 'bc_monthly', ivarid)
                  call check_ok('Variable bc_mon missing', &
                            &   'AEROSOL FILE ERROR')
                  istart(1) = 1
                  istart(2) = 1
                  icount(1) = jx
                  icount(2) = iy
                  icount(3) = 1
                  do m = 1 , 12
                    istart(3) = m
                    istatus = nf90_get_var(ncid,ivarid,toto,  &
                                       &   istart,icount)
                    call check_ok('Variable bc_mon read err', &
                            &   'AEROSOL FILE ERROR')
                    do j = 1 , jx
                      do i = 1 , iy
                        chemsrc(i,j,m,itr) = chemsrc(i,j,m,itr) + &
                                       &        dble(toto(j,i))
                      end do
                    end do
                  end do
                end if
              else if ( aerctl(1:2) == 'OC' ) then
                if ( aertyp(4:4) == '1' ) then
                  istatus = nf90_inq_varid(ncid, 'oc', ivarid)
                  call check_ok('Variable oc missing', &
                            &   'AEROSOL FILE ERROR')
                  istatus = nf90_get_var(ncid, ivarid, toto)
                  call check_ok('Variable oc read error', &
                            &   'AEROSOL FILE ERROR')
                  do m = 1 , 12
                    do j = 1 , jx
                      do i = 1 , iy
                        chemsrc(i,j,m,itr) = dble(toto(j,i))
                      end do
                    end do
                  end do
                end if
                if ( aertyp(5:5) == '1' ) then
                  istatus = nf90_inq_varid(ncid, 'oc_monthly', ivarid)
                  call check_ok('Variable oc_mon missing', &
                            &   'AEROSOL FILE ERROR')
                  istart(1) = 1
                  istart(2) = 1
                  icount(1) = jx
                  icount(2) = iy
                  icount(3) = 1
                  do m = 1 , 12
                    istart(3) = m
                    istatus = nf90_get_var(ncid,ivarid,toto, &
                                       &   istart,icount)
                    call check_ok('Variable oc_mon read err', &
                            &   'AEROSOL FILE ERROR')
                    do j = 1 , jx
                      do i = 1 , iy
                        chemsrc(i,j,m,itr) = chemsrc(i,j,m,itr) + &
                                              & dble(toto(j,i))
                      end do
                    end do
                  end do
                end if
              end if
            end if
          end do

          istatus = nf90_close(ncid)
          call check_ok('Error Closing Aerosol file '//trim(aername), &
                      &   'AEROSOL FILE CLOSE ERROR')

        end subroutine read_aerosol

        subroutine open_icbc(idate)
          use netcdf
          integer , intent(in) :: idate
          character(10) :: cdate
          integer :: idimid , itvar , i , chkdiff
          real(8) , dimension(:) , allocatable :: icbc_xtime
          character(64) :: icbc_timeunits
          integer :: iyy , jxx , kzz

          call close_icbc
          write (cdate, '(i10)') idate
          icbcname = trim(dirglob)//pthsep//trim(domname)//'_ICBC.'// &
              &      cdate//'.nc'
          istatus = nf90_open(icbcname, nf90_nowrite, ibcin)
          call check_ok('Error Opening ICBC file '//trim(icbcname), &
                     &  'ICBC FILE OPEN ERROR')
          ibcrec = 1
          ibcnrec = 0
          istatus = nf90_inq_dimid(ibcin, 'iy', idimid)
          call check_ok('Dimension iy missing', 'ICBC FILE ERROR')
          istatus = nf90_inquire_dimension(ibcin, idimid, len=iyy)
          call check_ok('Dimension iy read error','ICBC FILE ERROR')
          istatus = nf90_inq_dimid(ibcin, 'jx', idimid)
          call check_ok('Dimension jx missing', 'ICBC FILE ERROR')
          istatus = nf90_inquire_dimension(ibcin, idimid, len=jxx)
          call check_ok('Dimension jx read error', 'ICBC FILE ERROR')
          istatus = nf90_inq_dimid(ibcin, 'kz', idimid)
          call check_ok('Dimension kz missing', 'ICBC FILE ERROR')
          istatus = nf90_inquire_dimension(ibcin, idimid, len=kzz)
          call check_ok('Dimension kz read error', 'ICBC FILE ERROR')
          if ( iyy /= iy .or. jxx /= jx .or. kzz /= kz ) then
            write (6,*) 'Error: dims from regcm.in and ICBC file ', &
                        'differ.'
            write (aline,*) 'Input namelist : IY=' , iy , &
                   & '  JX=' ,  jx , '  KZ=' , kz
            call say
            write (aline,*) 'ICBC file      : IY=' , iyy , &
                   & '  JX=' ,  jxx , '  KZ=' , kzz
            call say
            call fatal(__FILE__,__LINE__,'DIMENSION MISMATCH')
          end if
          istatus = nf90_inq_dimid(ibcin, 'time', idimid)
          call check_ok('Dimension time missing', 'ICBC FILE ERROR')
          istatus = nf90_inquire_dimension(ibcin, idimid, len=ibcnrec)
          call check_ok('Dimension time read error', 'ICBC FILE ERROR')
          if ( ibcnrec < 1 ) then
            write (6,*) 'Time variable in ICBC has zero dim.'
            call fatal(__FILE__,__LINE__,'ICBC READ ERROR')
          end if
          istatus = nf90_inq_varid(ibcin, 'time', itvar)
          call check_ok('variable time missing', 'ICBC FILE ERROR')
          istatus = nf90_get_att(ibcin, itvar, 'units', icbc_timeunits)
          call check_ok('variable time units missing','ICBC FILE ERROR')
          allocate(icbc_idate(ibcnrec))
          allocate(icbc_xtime(ibcnrec))
          istatus = nf90_get_var(ibcin, itvar, icbc_xtime)
          call check_ok('variable time read error', 'ICBC FILE ERROR')
          do i = 1 , ibcnrec
            icbc_idate(i) = timeval2idate(icbc_xtime(i), icbc_timeunits)
          end do
          if ( ibcnrec > 1 ) then
            chkdiff = idnint(icbc_xtime(2) - icbc_xtime(1))
            if (chkdiff /= ibdyfrq) then
              write (6,*) 'Time variable in ICBC inconsistency.'
              write (6,*) 'Expecting ibdyfrq = ', ibdyfrq
              write (6,*) 'Found     ibdyfrq = ', chkdiff
              call fatal(__FILE__,__LINE__,'ICBC READ ERROR')
            end if
          end if
          deallocate(icbc_xtime)
          istatus = nf90_inq_varid(ibcin, 'ps', icbc_ivar(1))
          call check_ok('variable ps missing', 'ICBC FILE ERROR')
          istatus = nf90_inq_varid(ibcin, 'ts', icbc_ivar(2))
          call check_ok('variable ts missing', 'ICBC FILE ERROR')
          istatus = nf90_inq_varid(ibcin, 'u', icbc_ivar(3))
          call check_ok('variable u missing', 'ICBC FILE ERROR')
          istatus = nf90_inq_varid(ibcin, 'v', icbc_ivar(4))
          call check_ok('variable v missing', 'ICBC FILE ERROR')
          istatus = nf90_inq_varid(ibcin, 't', icbc_ivar(5))
          call check_ok('variable t missing', 'ICBC FILE ERROR')
          istatus = nf90_inq_varid(ibcin, 'qv', icbc_ivar(6))
          call check_ok('variable qv missing', 'ICBC FILE ERROR')
          istatus = nf90_inq_varid(ibcin, 'so4', icbc_ivar(7))
          if ( istatus == nf90_noerr) then
            lso4p = .true.
          end if
        end subroutine open_icbc

        subroutine read_icbc(idate,ps,ts,u,v,t,qv,so4)
          use netcdf
          implicit none
          integer , intent(in) :: idate
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

          if (idate > icbc_idate(ibcnrec) .or. idate < icbc_idate(1)) then
            write (6,*) 'Cannot find ', idate, ' in ICBC file'
            write (6,*) 'Range is : ', icbc_idate(1) , '-', &
                       & icbc_idate(ibcnrec)
            call fatal(__FILE__,__LINE__,'ICBC READ ERROR')
          end if 

          ibcrec = (idatediff(idate, icbc_idate(1)) / ibdyfrq) + 1
          istart(3) = ibcrec
          istart(2) = 1
          istart(1) = 1
          icount(3) = 1
          icount(2) = iy
          icount(1) = jx
          istatus = nf90_get_var(ibcin, icbc_ivar(1), xread(:,:,1), & 
                                 istart(1:3), icount(1:3))
          call check_ok('variable ps read error', 'ICBC FILE ERROR')
          ps = dble(transpose(xread(:,:,1)))
          istatus = nf90_get_var(ibcin, icbc_ivar(2), xread(:,:,1), & 
                                 istart(1:3), icount(1:3))
          call check_ok('variable ts read error', 'ICBC FILE ERROR')
          ts = dble(transpose(xread(:,:,1)))
          istart(4) = ibcrec
          istart(3) = 1
          istart(2) = 1
          istart(1) = 1
          icount(4) = 1
          icount(3) = kz
          icount(2) = iy
          icount(1) = jx
          istatus = nf90_get_var(ibcin, icbc_ivar(3), xread, &
                    &            istart, icount)
          call check_ok('variable u read error', 'ICBC FILE ERROR')
          do k = 1 , kz
            do j = 1 , jx
              do i = 1 , iy
                u(i,k,j) = dble(xread(j,i,k))
              end do
            end do
          end do
          istatus = nf90_get_var(ibcin, icbc_ivar(4), xread, &
                    &            istart, icount)
          call check_ok('variable v read error', 'ICBC FILE ERROR')
          do k = 1 , kz
            do j = 1 , jx
              do i = 1 , iy
                v(i,k,j) = dble(xread(j,i,k))
              end do
            end do
          end do
          istatus = nf90_get_var(ibcin, icbc_ivar(5), xread, &
                    &            istart, icount)
          call check_ok('variable t read error', 'ICBC FILE ERROR')
          do k = 1 , kz
            do j = 1 , jx
              do i = 1 , iy
                t(i,k,j) = dble(xread(j,i,k))
              end do
            end do
          end do
          istatus = nf90_get_var(ibcin, icbc_ivar(6), xread, &
                    &            istart, icount)
          call check_ok('variable qv read error', 'ICBC FILE ERROR')
          do k = 1 , kz
            do j = 1 , jx
              do i = 1 , iy
                qv(i,k,j) = dble(xread(j,i,k))
              end do
            end do
          end do
          if (lso4p) then
            istatus = nf90_get_var(ibcin, icbc_ivar(7), xread, &
                    &              istart, icount)
            call check_ok('variable so4 read error', 'ICBC FILE ERROR')
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
            call check_ok('Error Closing ICBC file '//trim(icbcname), &
                     &  'ICBC FILE ERROR')
            if (allocated(icbc_idate)) deallocate(icbc_idate)
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
            call check_ok('Error Closing '//ctype//' file', &
                     &  ctype//' FILE ERROR')
            ncid = -1
          end if
        end subroutine close_common

        function icbc_search(idate)
          implicit none
          integer :: icbc_search
          integer , intent(in) :: idate
          if (idate > icbc_idate(ibcnrec) .or. idate < icbc_idate(1)) then
            icbc_search = -1
          else
            icbc_search = (idatediff(idate, icbc_idate(1))/ibdyfrq)+1
          end if 
        end function icbc_search

        subroutine prepare_common_out(idate,ctype)
          use netcdf
          implicit none
          integer , intent(in) :: idate
          character(3) , intent(in) :: ctype
          character(64) :: title
          character(32) :: fbname , csdate
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
          integer :: itvar , imapvar , iyy , im , id , ih , i , j , ibnd
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

          write (fterr, '(a3,a)') ctype, ' FILE ERROR'
          write (fbname,'(a,a,i10)') trim(ctype), '.', idate
          ofname = trim(dirout)//pthsep//trim(domname)// &
                &  '_'//trim(fbname)//'.nc'
          call split_idate(idate0,iyy,im,id,ih)
          write (csdate,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a)') &
               & iyy,'-',im,'-',id,' ',ih,':00:00 UTC'

          write (aline, *) 'Opening new output file ', trim(ofname)
          call say

#ifdef NETCDF4_HDF5
          istatus = nf90_create(ofname, &
                    ior(ior(nf90_clobber,nf90_hdf5),nf90_classic_model),&
                    ncid)
#else
          istatus = nf90_create(ofname, nf90_clobber, ncid)
#endif
          call check_ok(('Error creating file '//trim(ofname)), fterr)

          istatus = nf90_put_att(ncid, nf90_global, 'title', title)
          call check_ok('Error adding global title', fterr)
          istatus = nf90_put_att(ncid, nf90_global, 'institution', &
                   & 'ICTP')
          call check_ok('Error adding global institution', fterr)
          istatus = nf90_put_att(ncid, nf90_global, 'source', &
                   & 'RegCM Model '//'SVN_REV'//' simulation output')
          call check_ok('Error adding global source', fterr)
          istatus = nf90_put_att(ncid, nf90_global, 'Conventions', &
                   & 'CF-1.4')
          call check_ok('Error adding global Conventions', fterr)
          call date_and_time(values=tvals)
          write (history,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a)') &
               tvals(1) , '-' , tvals(2) , '-' , tvals(3) , ' ' ,       &
               tvals(5) , ':' , tvals(6) , ':' , tvals(7) ,             &
               ' : Created by RegCM model'
          istatus = nf90_put_att(ncid, nf90_global, 'history', history)
          call check_ok('Error adding global history', fterr)
          istatus = nf90_put_att(ncid, nf90_global, 'references', &
                   & 'http://eforge.escience-lab.org/gf/project/regcm')
          call check_ok('Error adding global references', fterr)
          istatus = nf90_put_att(ncid, nf90_global, 'experiment', &
                   & domname)
          call check_ok('Error adding global experiment', fterr)
          istatus = nf90_put_att(ncid, nf90_global, 'projection', iproj)
          call check_ok('Error adding global projection', fterr)
          if (iproj == 'LAMCON') then
            istatus = nf90_put_att(ncid, nf90_global, &
                         'grid_mapping_name', 'lambert_conformal_conic')
            call check_ok('Error adding global grid_mapping_name',fterr)
          else if (iproj == 'POLSTR') then
            istatus = nf90_put_att(ncid, nf90_global, &
                         'grid_mapping_name', 'stereographic')
            call check_ok('Error adding global grid_mapping_name',fterr)
          else if (iproj == 'NORMER') then
            istatus = nf90_put_att(ncid, nf90_global, &
                         'grid_mapping_name', 'mercator')
            call check_ok('Error adding global grid_mapping_name',fterr)
          else if (iproj == 'ROTMER') then
            istatus = nf90_put_att(ncid, nf90_global, &
                  'grid_mapping_name', 'rotated_latitude_longitude')
            call check_ok('Error adding global grid_mapping_name',fterr)
          end if
          istatus = nf90_put_att(ncid, nf90_global,   &
                   &   'grid_size_in_meters', ds*d_1000)
          call check_ok('Error adding global gridsize', fterr)
          istatus = nf90_put_att(ncid, nf90_global,   &
                   &   'latitude_of_projection_origin', clat)
          call check_ok('Error adding global clat', fterr)
          istatus = nf90_put_att(ncid, nf90_global,   &
                   &   'longitude_of_projection_origin', clon)
          call check_ok('Error adding global clon', fterr)
          istatus = nf90_put_att(ncid, nf90_global,   &
                   &   'longitude_of_central_meridian', clon)
          call check_ok('Error adding global gmtllon', fterr)
          if (iproj == 'ROTMER') then
            istatus = nf90_put_att(ncid, nf90_global, &
                     &   'grid_north_pole_latitude', plat)
            call check_ok('Error adding global plat', fterr)
            istatus = nf90_put_att(ncid, nf90_global, &
                     &   'grid_north_pole_longitude', plon)
            call check_ok('Error adding global plon', fterr)
          else if (iproj == 'LAMCON') then
            trlat(1) = real(truelatl)
            trlat(2) = real(truelath)
            istatus = nf90_put_att(ncid, nf90_global, &
                     &   'standard_parallel', trlat)
            call check_ok('Error adding global truelat', fterr)
          else if (iproj == 'NORMER') then
            istatus = nf90_put_att(ncid, nf90_global, &
                     &   'standard_parallel', clat)
            call check_ok('Error adding global truelat', fterr)
          else if (iproj == 'POLSTR') then
            trlat(1) = 1.0
            istatus = nf90_put_att(ncid, nf90_global, &
                   &   'scale_factor_at_projection_origin', trlat(1:1))
            call check_ok('Error adding global scfac', fterr)
          end if
!
!         ADD RUN PARAMETERS
!
          istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_IPCC_scenario', scenario)
          call check_ok('Error adding global scenario', fterr)
          call cdumlbcs
          istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_boundary_conditions' , trim(cdum))
          call check_ok('Error adding global lbcs', fterr)
          call cdumcums
          istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_cumulous_convection_scheme' , trim(cdum))
          call check_ok('Error adding global icup', fterr)
          if (icup == 2 .or. icup == 99 .or. icup == 98) then
            call cdumcumcl
            istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_convective_closure_assumption' , trim(cdum))
            call check_ok('Error adding global igcc', fterr)
          end if
          call cdumpbl
          istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_boundary_layer_scheme' , trim(cdum))
          call check_ok('Error adding global ibltyp', fterr)
          call cdummoist
          istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_moist_physics_scheme' , trim(cdum))
          call check_ok('Error adding global ipptls', fterr)
          call cdumocnflx
          istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_ocean_flux_scheme' , trim(cdum))
          call check_ok('Error adding global iocnflx', fterr)
          call cdumpgfs
          istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_pressure_gradient_force_scheme' , trim(cdum))
          call check_ok('Error adding global ipgf', fterr)
          call cdumemiss
          istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_use_emission_factor' , trim(cdum))
          call check_ok('Error adding global iemiss', fterr)
          call cdumlakes
          istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_use_lake_model' , trim(cdum))
          call check_ok('Error adding global lakemod', fterr)
          call cdumchems
          istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_chemistry' , trim(cdum))
          call check_ok('Error adding global ichem', fterr)
          call cdumdcsst
          istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_diurnal_cycle_sst' , trim(cdum))
          call check_ok('Error adding global dcsst', fterr)
          call cdumseaice
          istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_seaice_effect' , trim(cdum))
          call check_ok('Error adding global seaice', fterr)
          call cdumdesseas
          istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_seasonal_desert_albedo_effect' , trim(cdum))
          call check_ok('Error adding global desseas', fterr)
          istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_simulation_initial_start' , globidate1)
          call check_ok('Error adding global globidate1', fterr)
          istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_simulation_start' , idate1)
          call check_ok('Error adding global idate1', fterr)
          istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_simulation_expected_end' , idate2)
          call check_ok('Error adding global idate2', fterr)
          if (ifrest) then
            cdum = 'Yes'
          else
            cdum = 'No'
          end if
          istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_simulation_is_a_restart' , trim(cdum))
          call check_ok('Error adding global ifrest', fterr)
          istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_timestep_in_seconds' , dt)
          call check_ok('Error adding global dt', fterr)
          istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_timestep_in_minutes_solar_rad_calc' , radfrq)
          call check_ok('Error adding global radfrq', fterr)
          istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_timestep_in_seconds_bats_calc' , abatm)
          call check_ok('Error adding global abatm', fterr)
          istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_timestep_in_hours_radiation_calc' , abemh)
          call check_ok('Error adding global abemh', fterr)
          istatus = nf90_put_att(ncid, nf90_global,  &
                  'model_timestep_in_hours_boundary_input' , ibdyfrq)
          call check_ok('Error adding global abemh', fterr)
!
!         ADD DIMENSIONS
!
          if (ctype == 'SUB') then
            istatus = nf90_def_dim(ncid, 'iy', o_nig, idims(2))
            call check_ok('Error creating dimension iy', fterr)
            istatus = nf90_def_dim(ncid, 'jx', o_njg, idims(1))
            call check_ok('Error creating dimension jx', fterr)
          else
            istatus = nf90_def_dim(ncid, 'iy', o_ni, idims(2))
            call check_ok('Error creating dimension iy', fterr)
            istatus = nf90_def_dim(ncid, 'jx', o_nj, idims(1))
            call check_ok('Error creating dimension jx', fterr)
          end if
          istatus = nf90_def_dim(ncid, 'time', nf90_unlimited, &
                              &  idims(3))
          call check_ok('Error creating dimension time', fterr)
          istatus = nf90_def_dim(ncid, 'kz', kz, idims(4))
          call check_ok('Error creating dimension kz', fterr)
!
!         OUT TYPE DEPENDENT DIMENSIONS
!
          if (ctype == 'SRF' .or. ctype == 'SUB') then
            istatus = nf90_def_dim(ncid, 'm10', 1, idims(5))
            call check_ok('Error creating dimension m10', fterr)
            istatus = nf90_def_dim(ncid, 'm2', 1, idims(6))
            call check_ok('Error creating dimension m2', fterr)
            istatus = nf90_def_dim(ncid, 'soil_layer', 2, idims(7))
            call check_ok('Error creating dimension soil_layer', fterr)
            istatus = nf90_def_var(ncid, 'm10', nf90_float, idims(5), &
                                   isrvvar(1))
            call check_ok('Error adding variable m10', fterr)
            istatus = nf90_put_att(ncid, isrvvar(1), 'standard_name', &
                            &  'altitude')
            call check_ok('Error adding m10 standard_name', fterr)
            istatus = nf90_put_att(ncid, isrvvar(1), 'long_name', &
                            &  'Convenience 10 m elevation level')
            call check_ok('Error adding m10 long_name', fterr)
            istatus = nf90_put_att(ncid, isrvvar(1), 'positive', 'up')
            call check_ok('Error adding m10 positive', fterr)
            istatus = nf90_put_att(ncid, isrvvar(1), 'units', 'm')
            call check_ok('Error adding m10 units', fterr)
            istatus = nf90_def_var(ncid, 'm2', nf90_float, idims(6), &
                                   isrvvar(2))
            call check_ok('Error adding variable m2', fterr)
            istatus = nf90_put_att(ncid, isrvvar(2), 'standard_name', &
                            &  'altitude')
            call check_ok('Error adding m2 standard_name', fterr)
            istatus = nf90_put_att(ncid, isrvvar(2), 'long_name', &
                            &  'Convenience 2 m elevation level')
            call check_ok('Error adding m2 long_name', fterr)
            istatus = nf90_put_att(ncid, isrvvar(2), 'positive', 'up')
            call check_ok( 'Error adding m2 positive', fterr)
            istatus = nf90_put_att(ncid, isrvvar(2), 'units', 'm')
            call check_ok( 'Error adding m2 units', fterr)
            istatus = nf90_def_var(ncid, 'layer', nf90_float, idims(7), &
                                   isrvvar(3))
            call check_ok('Error adding variable layer', fterr)
            istatus = nf90_put_att(ncid, isrvvar(3), 'standard_name', &
                            &  'model_level_number')
            call check_ok('Error adding layer standard_name', fterr)
            istatus = nf90_put_att(ncid, isrvvar(3), 'long_name', &
                            &  'Surface and root zone')
            call check_ok('Error adding layer long_name', fterr)
            istatus = nf90_put_att(ncid, isrvvar(3), 'positive', 'down')
            call check_ok('Error adding layer positive', fterr)
            istatus = nf90_put_att(ncid, isrvvar(3), 'units', '1')
            call check_ok('Error adding layer units', fterr)
          end if
          if (ctype == 'SRF') then
            istatus = nf90_def_dim(ncid, 'nv', 2, idims(8))
            call check_ok('Error creating dimension nv', fterr)
          end  if
          if (ctype == 'CHE') then
            istatus = nf90_def_dim(ncid, 'tracer', ntr, idims(9))
            call check_ok('Error creating dimension tracer', fterr)
            istatus = nf90_def_dim(ncid, 'dust', nbin, ibin)
            call check_ok('Error creating dimension dust', fterr)
            istatus = nf90_def_dim(ncid, 'bnd', 2, ibnd)
            call check_ok('Error creating dimension dust', fterr)
            istatus = nf90_def_dim(ncid, 'namelen', 6, inmlen(1))
            call check_ok('Error creating dimension namelen', fterr)
            inmlen(2) = idims(9)
            idpv(1) = idims(9)
            idpv(2) = ibnd
            ibinsiz(1) = ibin
            ibinsiz(2) = ibnd
          end if
          if (ctype == 'LAK') then
            istatus = nf90_def_dim(ncid, 'depth', ndpmax, idims(10))
            call check_ok('Error creating dimension depth', fterr)
          end if
          istatus = nf90_def_var(ncid, 'rcm_map', nf90_int, &
                                 varid=imapvar)
          call check_ok('Error adding variable rcm_map', fterr)
          if (iproj == 'LAMCON') then
            istatus = nf90_put_att(ncid, imapvar, &
                         'grid_mapping_name', 'lambert_conformal_conic')
            call check_ok('Error adding rcm_map grid_mapping_name',fterr)
          else if (iproj == 'POLSTR') then
            istatus = nf90_put_att(ncid, imapvar, &
                         'grid_mapping_name', 'stereographic')
            call check_ok('Error adding rcm_map grid_mapping_name',fterr)
          else if (iproj == 'NORMER') then
            istatus = nf90_put_att(ncid, imapvar, &
                         'grid_mapping_name', 'mercator')
            call check_ok('Error adding rcm_map grid_mapping_name',fterr)
          else if (iproj == 'ROTMER') then
            istatus = nf90_put_att(ncid, imapvar, &
                  'grid_mapping_name', 'rotated_latitude_longitude')
            call check_ok('Error adding rcm_map grid_mapping_name',fterr)
          end if
          istatus = nf90_put_att(ncid, imapvar,   &
                   &   'grid_size_in_meters', ds*d_1000)
          call check_ok('Error adding rcm_map gridsize', fterr)
          istatus = nf90_put_att(ncid, imapvar,   &
                   &   'latitude_of_projection_origin', clat)
          call check_ok('Error adding rcm_map clat', fterr)
          istatus = nf90_put_att(ncid, imapvar,   &
                   &   'longitude_of_projection_origin', clon)
          call check_ok('Error adding rcm_map clon', fterr)
          istatus = nf90_put_att(ncid, imapvar,   &
                   &   'longitude_of_central_meridian', clon)
          call check_ok('Error adding rcm_map gmtllon', fterr)
          if (iproj == 'ROTMER') then
            istatus = nf90_put_att(ncid, imapvar, &
                     &   'grid_north_pole_latitude', plat)
            call check_ok('Error adding rcm_map plat', fterr)
            istatus = nf90_put_att(ncid, imapvar, &
                     &   'grid_north_pole_longitude', plon)
            call check_ok('Error adding rcm_map plon', fterr)
          else if (iproj == 'LAMCON') then
            trlat(1) = real(truelatl)
            trlat(2) = real(truelath)
            istatus = nf90_put_att(ncid, imapvar, &
                     &   'standard_parallel', trlat)
            call check_ok('Error adding rcm_map truelat', fterr)
          else if (iproj == 'NORMER') then
            istatus = nf90_put_att(ncid, imapvar, &
                     &   'standard_parallel', clat)
            call check_ok('Error adding rcm_map truelat', fterr)
          else if (iproj == 'POLSTR') then
            trlat(1) = 1.0
            istatus = nf90_put_att(ncid, imapvar, &
                   &   'scale_factor_at_projection_origin', trlat(1:1))
            call check_ok('Error adding rcm_map scfac', fterr)
          end if
          istatus = nf90_def_var(ncid, 'sigma', nf90_float, &
                              &  idims(4), izvar(1))
          call check_ok('Error adding variable sigma', fterr)
          istatus = nf90_put_att(ncid, izvar(1), 'standard_name', &
                            &  'atmosphere_sigma_coordinate')
          call check_ok('Error adding sigma standard_name', fterr)
          istatus = nf90_put_att(ncid, izvar(1), 'long_name', &
                            &  'Sigma at model layers')
          call check_ok('Error adding sigma long_name', fterr)
          istatus = nf90_put_att(ncid, izvar(1), 'units', '1')
          call check_ok('Error adding sigma units', fterr)
          istatus = nf90_put_att(ncid, izvar(1), 'axis', 'Z')
          call check_ok('Error adding sigma axis', fterr)
          istatus = nf90_put_att(ncid, izvar(1), 'positive', 'down')
          call check_ok('Error adding sigma positive', fterr)
          istatus = nf90_put_att(ncid, izvar(1), 'formula_terms',  &
                     &         'sigma: sigma ps: ps ptop: ptop')
          call check_ok('Error adding sigma formula_terms', fterr)
          if (ctype == 'LAK') then
            istatus = nf90_def_var(ncid, 'depth', nf90_float, &
                                &  idims(10), izvar(3))
            call check_ok('Error adding variable depth', fterr)
            istatus = nf90_put_att(ncid, izvar(3), 'standard_name', &
                              &  'depth')
            call check_ok('Error adding depth standard_name', fterr)
            istatus = nf90_put_att(ncid, izvar(3), 'long_name', &
                              &  'Depth below surface')
            call check_ok('Error adding depth long_name', fterr)
            istatus = nf90_put_att(ncid, izvar(3), 'units', 'm')
            call check_ok('Error adding depth units', fterr)
            istatus = nf90_put_att(ncid, izvar(3), 'axis', 'Z')
            call check_ok('Error adding depth axis', fterr)
            istatus = nf90_put_att(ncid, izvar(3), 'positive', 'down')
            call check_ok('Error adding depth positive', fterr)
          end if
          istatus = nf90_def_var(ncid, 'ptop', nf90_float, &
                           &   varid=izvar(2))
          call check_ok('Error adding variable ptop', fterr)
          istatus = nf90_put_att(ncid, izvar(2), 'standard_name',  &
                            &  'air_pressure')
          call check_ok('Error adding ptop standard_name', fterr)
          istatus = nf90_put_att(ncid, izvar(2), 'long_name', &
                            &  'Pressure at model top')
          call check_ok('Error adding ptop long_name', fterr)
          istatus = nf90_put_att(ncid, izvar(2), 'units', 'hPa')
          call check_ok('Error adding ptop units', fterr)
          istatus = nf90_def_var(ncid, 'iy', nf90_float, idims(2), &
                            &  ivvar(1))
          call check_ok('Error adding variable iy', fterr)
          istatus = nf90_put_att(ncid, ivvar(1), 'standard_name',  &
                            &  'projection_y_coordinate')
          call check_ok('Error adding iy standard_name', fterr)
          istatus = nf90_put_att(ncid, ivvar(1), 'long_name', &
                            &  'y-coordinate in Cartesian system')
          call check_ok('Error adding iy long_name', fterr)
          istatus = nf90_put_att(ncid, ivvar(1), 'units', 'km')
          call check_ok('Error adding iy units', fterr)
          istatus = nf90_def_var(ncid, 'jx', nf90_float, idims(1), &
                            &  ivvar(2))
          call check_ok('Error adding variable jx', fterr)
          istatus = nf90_put_att(ncid, ivvar(2), 'standard_name', &
                            &  'projection_x_coordinate')
          call check_ok('Error adding jx standard_name', fterr)
          istatus = nf90_put_att(ncid, ivvar(2), 'long_name', &
                            &  'x-coordinate in Cartesian system')
          call check_ok('Error adding jx long_name', fterr)
          istatus = nf90_put_att(ncid, ivvar(2), 'units', 'km')
          call check_ok('Error adding jx units', fterr)
          istatus = nf90_def_var(ncid, 'xlat', nf90_float, &
                             &   idims(1:2), illtpvar(1))
          call check_ok('Error adding variable xlat', fterr)
          istatus = nf90_put_att(ncid, illtpvar(1), 'standard_name', &
                            &  'latitude')
          call check_ok('Error adding xlat standard_name', fterr)
          istatus = nf90_put_att(ncid, illtpvar(1), 'long_name', &
                            &  'Latitude at cross points')
          call check_ok('Error adding xlat long_name', fterr)
          istatus = nf90_put_att(ncid, illtpvar(1), 'units', &
                            &  'degrees_north')
          call check_ok('Error adding xlat units', fterr)
          istatus = nf90_def_var(ncid, 'xlon', nf90_float, &
                             &   idims(1:2), illtpvar(2))
          call check_ok('Error adding variable xlon', fterr)
          istatus = nf90_put_att(ncid, illtpvar(2), 'standard_name', &
                            &  'longitude')
          call check_ok('Error adding xlon standard_name', fterr)
          istatus = nf90_put_att(ncid, illtpvar(2), 'long_name', &
                            &  'Longitude at cross points')
          call check_ok('Error adding xlon long_name', fterr)
          istatus = nf90_put_att(ncid, illtpvar(2), 'units',  &
                            &  'degrees_east')
          call check_ok('Error adding xlon units', fterr)
          istatus = nf90_def_var(ncid, 'topo', nf90_float, &
                             &   idims(1:2), illtpvar(3))
          call check_ok('Error adding variable topo', fterr)
          istatus = nf90_put_att(ncid, illtpvar(3), 'standard_name', &
                            &  'surface_altitude')
          call check_ok('Error adding topo standard_name', fterr)
          istatus = nf90_put_att(ncid, illtpvar(3), 'long_name',     &
                            &  'Domain surface elevation')
          call check_ok('Error adding topo long_name', fterr)
          istatus = nf90_put_att(ncid, illtpvar(3), 'units', 'm')
          call check_ok('Error adding topo units', fterr)
          istatus = nf90_put_att(ncid, illtpvar(3), 'coordinates', &
                            &  'xlat xlon')
          call check_ok('Error adding topo coordinates', fterr)
          istatus = nf90_put_att(ncid, illtpvar(3), 'grid_mapping', &
                            &  'rcm_map')
          call check_ok('Error adding topo grid_mapping', fterr)
          istatus = nf90_def_var(ncid, 'mask', nf90_float, &
                             &   idims(1:2), illtpvar(4))
          call check_ok('Error adding variable mask', fterr)
          istatus = nf90_put_att(ncid, illtpvar(4), 'standard_name', &
                            &  'landmask')
          call check_ok('Error adding mask standard_name', fterr)
          istatus = nf90_put_att(ncid, illtpvar(4), 'long_name',     &
                            &  'Domain land/ocean mask')
          call check_ok('Error adding mask long_name', fterr)
          istatus = nf90_put_att(ncid, illtpvar(4), 'units', '1')
          call check_ok('Error adding mask units', fterr)
          istatus = nf90_put_att(ncid, illtpvar(4), 'coordinates', &
                            &  'xlat xlon')
          call check_ok('Error adding mask coordinates', fterr)
          istatus = nf90_put_att(ncid, illtpvar(4), 'grid_mapping', &
                            &  'rcm_map')
          call check_ok('Error adding mask grid_mapping', fterr)
          istatus = nf90_def_var(ncid, 'time', nf90_double, &
                               & idims(3:3), itvar)
          call check_ok('Error adding variable time', fterr)
          istatus = nf90_put_att(ncid, itvar, 'standard_name', 'time')
          call check_ok('Error adding time standard_name', fterr)
          istatus = nf90_put_att(ncid, itvar, 'long_name', 'time')
          call check_ok('Error adding time long_name', fterr)
          istatus = nf90_put_att(ncid, itvar, 'calendar', 'standard')
          call check_ok('Error adding time calendar', fterr)
          istatus = nf90_put_att(ncid, itvar, 'units', &
                             &   'hours since '//csdate)
          call check_ok('Error adding time units', fterr)
          if (ctype == 'SRF') then
            istatus = nf90_put_att(ncid, itvar, 'bounds', 'tbnds')
            call check_ok('Error adding time bounds', fterr)
          end if
          istatus = nf90_def_var(ncid, 'ps', nf90_float, &
                             &   idims(1:3), illtpvar(5))
          call check_ok('Error adding variable ps', fterr)
          istatus = nf90_put_att(ncid, illtpvar(5), 'standard_name', &
                            &  'surface_air_pressure')
          call check_ok('Error adding ps standard_name', fterr)
          istatus = nf90_put_att(ncid, illtpvar(5), 'long_name',     &
                            &  'Surface pressure')
          call check_ok('Error adding ps long_name', fterr)
          istatus = nf90_put_att(ncid, illtpvar(5), 'units', 'hPa')
          call check_ok('Error adding ps units', fterr)
          istatus = nf90_put_att(ncid, illtpvar(5), 'coordinates', &
                            &  'xlat xlon')
          call check_ok('Error adding ps coordinates', fterr)
          istatus = nf90_put_att(ncid, illtpvar(5), 'grid_mapping', &
                            &  'rcm_map')
          call check_ok('Error adding ps grid_mapping', fterr)

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
            call addvara(ncid,ctype,'u','eastward_wind', &
                'U component (westerly) of wind','m s-1', &
                tzyx,.false.,iatmvar(3))
            call addvara(ncid,ctype,'v','northward_wind', &
                'V component (southerly) of wind','m s-1', &
                tzyx,.false.,iatmvar(4))
            call addvara(ncid,ctype,'omega', &
                'lagrangian_tendency_of_air_pressure', &
                'Pressure velocity','hPa s-1',tzyx,.false.,iatmvar(5))
            call addvara(ncid,ctype,'t','air_temperature', &
                'Temperature','K',tzyx,.false.,iatmvar(6))
            call addvara(ncid,ctype,'qv','humidity_mixing_ratio', &
                'Water vapor mixing ratio','kg kg-1', &
                tzyx,.false.,iatmvar(7))
            call addvara(ncid,ctype,'qc', &
                'cloud_liquid_water_mixing_ratio', &
                'Cloud water mixing ratio','kg kg-1', &
                tzyx,.false.,iatmvar(8))
            call addvara(ncid,ctype,'tpr','precipitation_flux', &
                'Total daily precipitation rate','kg m-2 day-1', &
                tyx,.false.,iatmvar(9))
            call addvara(ncid,ctype,'tgb','soil_temperature', &
                'Lower groud temperature','K', &
                tyx,.false.,iatmvar(10))
            call addvara(ncid,ctype,'swt', &
                'moisture_content_of_soil_layer', &
                'Total soil water','kg m-2',tyx,.true.,iatmvar(11))
            call addvara(ncid,ctype,'rno','runoff_flux', &
                'Runoff accumulated infiltration','kg m-2 day-1', &
                tyx,.true.,iatmvar(12))
          else if (ctype == 'SRF') then
            isrfvar = -1
            isrfvar(1) = itvar
            istatus = nf90_def_var(ncid, 'tbnds', nf90_double, &
                                   (/idims(8),idims(3)/), isrfvar(2))
            call check_ok('Error adding variable tbnds', fterr)
            istatus = nf90_put_att(ncid, isrfvar(2), 'calendar', &
                                   'standard')
            call check_ok('Error adding tbnds calendar', fterr)
            istatus = nf90_put_att(ncid, isrfvar(2), 'units', &
                             &   'hours since '//csdate)
            call check_ok('Error adding tbnds units', fterr)
            isrfvar(3) = illtpvar(5)
            call addvara(ncid,ctype,'u10m','eastward_wind', &
                '10 meters U component (westerly) of wind','m s-1', &
                t10yx,.false.,isrfvar(4))
            call addvara(ncid,ctype,'v10m','northward_wind', &
                '10 meters V component (southerly) of wind','m s-1', &
                t10yx,.false.,isrfvar(5))
            call addvara(ncid,ctype,'uvdrag', &
                'surface_drag_coefficient_in_air', &
                'Surface drag stress','1',tyx,.false.,isrfvar(6))
            call addvara(ncid,ctype,'tg','surface_temperature', &
                'Ground temperature','K',tyx,.false.,isrfvar(7))
            call addvara(ncid,ctype,'tlef','canopy_temperature', &
                'Foliage temperature','K',tyx,.true.,isrfvar(8))
            call addvara(ncid,ctype,'t2m','air_temperature', &
                '2 meters temperature','K',t2yx,.false.,isrfvar(9))
            call addvara(ncid,ctype,'q2m','humidity_mixing_ratio', &
                '2 meters vapour mixing ratio','kg kg-1', &
                t2yx,.false.,isrfvar(10))
            call addvara(ncid,ctype,'smw','soil_moisture_content', &
                'Moisture content','kg kg-1',tlyx,.true.,isrfvar(11))
            call addvara(ncid,ctype,'tpr','precipitation_flux', &
                'Total precipitation','kg m-2 day -1',tyx,&
                .false.,isrfvar(12))
            call addvara(ncid,ctype,'evp','water_evaporation_flux', &
                'Total evapotranspiration','kg m-2 day-1', &
                tyx,.false.,isrfvar(13))
            call addvara(ncid,ctype,'runoff','surface_runoff_flux', &
                'Surface runoff','kg m-2 day-1',tyx,.true.,isrfvar(14))
            call addvara(ncid,ctype,'scv','snowfall_flux', &
                'Snow precipitation','kg m-2 day-1',tyx,.true.,isrfvar(15))
            call addvara(ncid,ctype,'sena', &
                'surface_downward_sensible_heat_flux', &
                'Sensible heat flux','W m-2',tyx,.false.,isrfvar(16))
            call addvara(ncid,ctype,'flw', &
                'net_upward_longwave_flux_in_air', &
                'Net infrared energy flux','W m-2', &
                tyx,.false.,isrfvar(17))
            call addvara(ncid,ctype,'fsw', &
                'surface_downwelling_shortwave_flux_in_air', &
                'Solar absorbed energy flux','W m-2', &
                tyx,.false.,isrfvar(18))
            call addvara(ncid,ctype,'fld', &
                'surface_downwelling_longwave_flux_in_air', &
                'Downward LW flux','W m-2',tyx,.false.,isrfvar(19))
            call addvara(ncid,ctype,'sina', &
              'net_downward_radiative_flux_at_top_of_atmosphere_model', &
              'Incident solar energy flux','W m-2', &
                tyx,.false.,isrfvar(20))
            call addvara(ncid,ctype,'prcv', 'convective_rainfall_flux', &
                'Convective precipitation','kg m-2 day-1', &
                tyx,.false.,isrfvar(21))
            call addvara(ncid,ctype,'zpbl', &
                'atmosphere_boundary_layer_thickness', &
                'PBL layer thickness','m',tyx,.false.,isrfvar(22))
            write (cmethodmax, '(a,i3,a)') 'time: maximum (interval: ', &
                   idint(batfrq) , ' hours)'
            write (cmethodmin, '(a,i3,a)') 'time: minimum (interval: ', &
                   idint(batfrq) , ' hours)'
            call addvara(ncid,ctype,'tgmax','surface_temperature', &
                'Maximum surface temperature','K', &
                tyx,.false.,isrfvar(23))
            istatus = nf90_put_att(ncid, isrfvar(23), 'cell_methods', &
                            &  cmethodmax)
            call check_ok('Error adding tgmax cell_methods', fterr)
            call addvara(ncid,ctype,'tgmin','surface_temperature', &
                'Minimum surface temperature','K', &
                tyx,.false.,isrfvar(24))
            istatus = nf90_put_att(ncid, isrfvar(24), 'cell_methods', &
                            &  cmethodmin)
            call check_ok('Error adding tgmin cell_methods', fterr)
            call addvara(ncid,ctype,'t2max','air_temperature', &
                'Maximum 2 meters temperature','K', &
                t2yx,.false.,isrfvar(25))
            istatus = nf90_put_att(ncid, isrfvar(25), 'cell_methods', &
                            &  cmethodmax)
            call check_ok('Error adding t2max cell_methods', fterr)
            call addvara(ncid,ctype,'t2min','air_temperature', &
                'Minimum 2 meters temperature','K', &
                t2yx,.false.,isrfvar(26))
            istatus = nf90_put_att(ncid, isrfvar(26), 'cell_methods', &
                            &  cmethodmin)
            call check_ok('Error adding t2min cell_methods', fterr)
            call addvara(ncid,ctype,'w10max','wind_speed', &
                'Maximum speed of 10m wind','m s-1', &
                t10yx,.false.,isrfvar(27))
            istatus = nf90_put_att(ncid, isrfvar(27), 'cell_methods', &
                            &  cmethodmax)
            call check_ok('Error adding w10max cell_methods', fterr)
            call addvara(ncid,ctype,'ps_min','air_pressure', &
                'Minimum of surface pressure','hPa', &
                tyx,.false.,isrfvar(28))
            istatus = nf90_put_att(ncid, isrfvar(28), 'cell_methods', &
                            &  cmethodmin)
            call check_ok('Error adding ps_min cell_methods', fterr)
            call addvara(ncid,ctype,'aldirs', &
                'surface_albedo_short_wave_direct', &
                'Surface albedo to direct short wave radiation', &
                '1',tyx,.false.,isrfvar(29))
            call addvara(ncid,ctype,'aldifs', &
                'surface_albedo_short_wave_diffuse', &
                'Surface albedo to diffuse short wave radiation', &
                '1',tyx,.false.,isrfvar(30))
            if ( iseaice == 1 .or. lakemod == 1 ) then
              call addvara(ncid,ctype,'seaice', &
                  'seaice_binary_mask', &
                  'Sea ice mask', &
                  '1',tyx,.false.,isrfvar(31))
            end if
          else if (ctype == 'SUB') then
            isubvar = -1
            isubvar(1) = itvar
            isubvar(2) = illtpvar(5)
            call addvara(ncid,ctype,'u10m','eastward_wind', &
                '10 meters U component (westerly) of wind','m s-1', &
                t10yx,.false.,isubvar(3))
            call addvara(ncid,ctype,'v10m','northward_wind', &
                '10 meters V component (southerly) of wind','m s-1', &
                t10yx,.false.,isubvar(4))
            call addvara(ncid,ctype,'uvdrag', &
                'surface_drag_coefficient_in_air', &
                'Surface drag stress','1', tyx, .false.,isubvar(5))
            call addvara(ncid,ctype,'tg','surface_temperature', &
                'Ground temperature','K', tyx, .false.,isubvar(6))
            call addvara(ncid,ctype,'tlef','canopy_temperature', &
                'Foliage temperature','K', tyx, .true.,isubvar(7))
            call addvara(ncid,ctype,'t2m','air_temperature', &
                '2 meters temperature','K', t2yx, .false.,isubvar(8))
            call addvara(ncid,ctype,'q2m','humidity_mixing_ratio', &
                '2 meters vapour mixing ratio','kg kg-1', &
                t2yx,.false.,isubvar(9))
            call addvara(ncid,ctype,'smw','soil_moisture_content', &
                'Moisture content','kg kg-1', tlyx, .true.,isubvar(10))
            call addvara(ncid,ctype,'tpr','precipitation_flux', &
                'Total precipitation','kg m-2 day-1',tyx, &
                .false.,isubvar(11))
            call addvara(ncid,ctype,'evp','water_evaporation_flux', &
                'Total evapotranspiration','kg m-2 day-1', &
                tyx,.false.,isubvar(12))
            call addvara(ncid,ctype,'runoff','surface_runoff_flux', &
                'Surface runoff','kg m-2 day-1',tyx,.true.,isubvar(13))
            call addvara(ncid,ctype,'scv','snowfall_flux', &
                'Snow precipitation','kg m-2 day-1',tyx,.true.,isubvar(14))
            call addvara(ncid,ctype,'sena', &
                'surface_downward_sensible_heat_flux', &
                'Sensible heat flux','W m-2',tyx,.false.,isubvar(15))
            call addvara(ncid,ctype,'prcv', 'convective_rainfall_flux', &
                'Convective precipitation','kg m-2 day-1', &
                tyx,.false.,isubvar(16))
          else if (ctype == 'RAD') then
            iradvar = -1
            iradvar(1) = itvar
            iradvar(2) = illtpvar(5)
            call addvara(ncid,ctype,'cld', &
                'cloud_area_fraction_in_atmosphere_layer', &
                'Cloud fractional cover','1',tzyx,.false.,iradvar(3))
            call addvara(ncid,ctype,'clwp', &
                'atmosphere_cloud_liquid_water_content', &
                'Cloud liquid water content','g m-2', &
                 tzyx,.false.,iradvar(4))
            call addvara(ncid,ctype,'qrs', &
                'tendency_of_air_temperature_due_to_shortwave_heating', &
                'Solar heating rate','K s-1',tzyx,.false.,iradvar(5))
            call addvara(ncid,ctype,'qrl', &
                'tendency_of_air_temperature_due_to_longwave_heating', &
                'Longwave cooling rate','K s-1',tzyx,.false.,iradvar(6))
            call addvara(ncid,ctype,'frsa', &
                'surface_downwelling_shortwave_flux_in_air', &
                'Surface absorbed solar flux','W m-2', &
                tyx,.false.,iradvar(7))
            call addvara(ncid,ctype,'frla', &
                'downwelling_longwave_flux_in_air', &
                'Longwave cooling of surface flux','W m-2', &
                tyx,.false.,iradvar(8))
            call addvara(ncid,ctype,'clrst', &
                'downwelling_shortwave_flux_in_air_assuming_clear_sky', &
                'clearsky total column absorbed solar flux','W m-2', &
                tyx,.false.,iradvar(9))
            call addvara(ncid,ctype,'clrss', &
                'net_downward_shortwave_flux_in_air_assuming_clear_sky', &
                'clearsky surface absorbed solar flux','W m-2', &
                tyx,.false.,iradvar(10))
            call addvara(ncid,ctype,'clrlt', &
                'toa_net_upward_longwave_flux_assuming_clear_sky', &
                'clearsky net upward LW flux at TOA','W m-2', &
                tyx,.false.,iradvar(11))
            call addvara(ncid,ctype,'clrls', &
                'net_upward_longwave_flux_in_air_assuming_clear_sky', &
                'clearsky LW cooling at surface','W m-2', &
                tyx,.false.,iradvar(12))
            call addvara(ncid,ctype,'solin', &
                'toa_instantaneous_shortwave_forcing', &
                'Instantaneous incident solar','W m-2', &
                tyx,.false.,iradvar(13))
            call addvara(ncid,ctype,'sabtp', &
              'atmosphere_net_rate_of_absorption_of_shortwave_energy', &
              'Total column absorbed solar flux','W m-2', &
                tyx,.false.,iradvar(14))
            call addvara(ncid,ctype,'firtp', &
                'atmosphere_net_rate_of_absorption_of_longwave_energy', &
                'net upward LW flux at TOA','W m-2', &
                tyx,.false.,iradvar(15))
          else if (ctype == 'CHE') then
            istatus = nf90_def_var(ncid, 'chtrname', nf90_char, &
                                   inmlen, ichname)
            call check_ok('Error adding variable chtrname', fterr)
            istatus = nf90_def_var(ncid, 'chtrsol', nf90_double, &
                                   idims(9), ichtrsol)
            call check_ok('Error adding variable chtrsol', fterr)
            istatus = nf90_def_var(ncid, 'chtrdpv', nf90_double, &
                                   idpv, ichtrdpv)
            call check_ok('Error adding variable chtrdpv', fterr)
            istatus = nf90_def_var(ncid, 'dustbinsiz', nf90_double, &
                                   ibinsiz, idubinsiz)
            call check_ok('Error adding variable dustbinsiz', fterr)
            ichevar = -1
            ichevar(1) = itvar
            ichevar(2) = illtpvar(5)
            call addvara(ncid,ctype,'trac', &
                'atmosphere_mixing_ratio_of_tracer', &
                'Tracers mixing ratios','kg kg-1', &
                tczyx,.false.,ichevar(3))
            call addvara(ncid,ctype,'aext8', &
                'aerosol_optical_depth', &
                'aer mix. aod.','1',tzyx,.false.,ichevar(4))
            call addvara(ncid,ctype,'assa8', &
                'aerosol_single_scattering_albedo', &
                'aer mix. sin. scat. alb','1',tzyx,.false.,ichevar(5))
            call addvara(ncid,ctype,'agfu8', &
                'aerosol_asymmetry_parameter', &
                'aer mix. sin. scat. alb','1',tzyx,.false.,ichevar(6))
            call addvara(ncid,ctype,'colb', &
                'instantaneous_column_burden', &
                'columnburden inst','mg m-2',tcyx,.false.,ichevar(7))
            call addvara(ncid,ctype,'wdlsc', &
                'tendency_of_wet_deposition_of_tracer'// &
                '_due_to_large_scale_precipitation', &
                'wet dep lgscale','mg m-2 day-1', &
                tcyx,.false.,ichevar(8))
            call addvara(ncid,ctype,'wdcvc', &
                'tendency_of_wet_deposition_of_tracer'// &
                '_due_to_convective_precipitation', &
                'wet dep convect','mg m-2 day-1', &
                tcyx,.false.,ichevar(9))
            call addvara(ncid,ctype,'sdrdp', &
                'tendency_of_dry_deposition_of_tracer', &
                'surf dry depos','mg m-2 day-1', &
                tcyx,.false.,ichevar(10))
            call addvara(ncid,ctype,'xgasc', &
                'tendency_of_gas_conversion_of_tracer', &
                'chem gas conv','mg m-2 day-1', &
                tcyx,.false.,ichevar(11))
            call addvara(ncid,ctype,'xaquc', &
                'tendency_of_aqueous_conversion_of_tracer', &
                'chem aqu conv','mg m-2 day-1', &
                tcyx,.false.,ichevar(12))
            call addvara(ncid,ctype,'emiss', &
                'tendency_of_surface_emission_of_tracer', &
                'surf emission','mg m-2 day-1', &
                tcyx,.false.,ichevar(13))
            call addvara(ncid,ctype,'acstoarf', &
                'toa_instantaneous_shortwave_radiative_forcing', &
                'TOArad SW forcing av.','W m-2', &
                tyx,.false.,ichevar(14))
            call addvara(ncid,ctype,'acstsrrf', &
                'surface_shortwave_radiative_forcing', &
                'SRFrad SW forcing av.','W m-2', &
                tyx,.false.,ichevar(15))
            call addvara(ncid,ctype,'acstalrf', &
                'toa_longwave_radiative_forcing', &
                'TOArad LW forcing av.','W m-2', &
                tyx,.false.,ichevar(16))
            call addvara(ncid,ctype,'acssrlrf', &
                'surface_longwave_radiative_forcing', &
                'SRFrad LW forcing av.','W m-2', &
                tyx,.false.,ichevar(17))
          else if (ctype == 'LAK') then
            ilakvar = -1
            ilakvar(1) = itvar
            ilakvar(2) = illtpvar(5)
            call addvara(ncid,ctype,'tg','surface_temperature', &
                'Ground temperature','K',tyx,.false.,ilakvar(3))
            call addvara(ncid,ctype,'tpr','precipitation_flux', &
                'Total precipitation','kg m-2 day-1',tyx,.false., &
                ilakvar(4))
            call addvara(ncid,ctype,'scv','snowfall_flux', &
                'Snow precipitation','kg m-2 day-1',tyx,.true., &
                ilakvar(5))
            call addvara(ncid,ctype,'sena', &
                'surface_downward_sensible_heat_flux', &
                'Sensible heat flux','W m-2',tyx,.false.,ilakvar(6))
            call addvara(ncid,ctype,'flw', &
                'net_upward_longwave_flux_in_air', &
                'Net infrared energy flux','W m-2', &
                tyx,.false.,ilakvar(7))
            call addvara(ncid,ctype,'fsw', &
                'surface_downwelling_shortwave_flux_in_air', &
                'Solar absorbed energy flux','W m-2', &
                tyx,.false.,ilakvar(8))
            call addvara(ncid,ctype,'fld', &
                'surface_downwelling_longwave_flux_in_air', &
                'Downward LW flux','W m-2',tyx,.false.,ilakvar(9))
            call addvara(ncid,ctype,'sina', &
              'net_downward_radiative_flux_at_top_of_atmosphere_model', &
              'Incident solar energy flux','W m-2', &
                tyx,.false.,ilakvar(10))
            call addvara(ncid,ctype,'aldirs', &
                'surface_albedo_short_wave_direct', &
                'Surface albedo to direct short wave radiation', &
                '1',tyx,.false.,ilakvar(11))
            call addvara(ncid,ctype,'aldifs', &
                'surface_albedo_short_wave_diffuse', &
                'Surface albedo to diffuse short wave radiation', &
                '1',tyx,.false.,ilakvar(12))
            call addvara(ncid,ctype,'evl', &
                'water_evaporation_flux_where_sea_ice', &
                'Water evaporation','mm sec-1',tyx,.true.,ilakvar(13))
            call addvara(ncid,ctype,'aveice', 'floating_ice_thickness', &
                'Floating ice thickness','mm',tyx,.true.,ilakvar(14))
            call addvara(ncid,ctype,'hsnow', &
                'surface_snow_thickness_where_sea_ice', &
                'Floating snow thickness','mm',tyx,.true.,ilakvar(15))
            call addvara(ncid,ctype,'tlake','temperature', &
                'Lake water temperature','K',tdyx,.true.,ilakvar(16))
          end if

          istatus = nf90_enddef(ncid)
          call check_ok('Error End Definitions NetCDF output', fterr)

          istatus = nf90_put_var(ncid, izvar(1), hsigma)
          call check_ok('Error variable sigma write', fterr)
          hptop = real(ptop*d_10)
          istatus = nf90_put_var(ncid, izvar(2), hptop)
          call check_ok('Error variable ptop write', fterr)
          if (ctype == 'LAK') then
            do i = 1 , ndpmax
              depth(i) = real(i)
            end do
            istatus = nf90_put_var(ncid, izvar(3), depth)
            call check_ok('Error variable depth write', fterr)
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
            call check_ok('Error variable iy write', fterr)
            istatus = nf90_put_var(ncid, ivvar(2), xjx(1:o_njg))
            call check_ok('Error variable jx write', fterr)
            istatus = nf90_put_var(ncid, illtpvar(1), ioxlat_s)
            call check_ok('Error variable xlat write', fterr)
            istatus = nf90_put_var(ncid, illtpvar(2), ioxlon_s)
            call check_ok('Error variable xlon write', fterr)
            istatus = nf90_put_var(ncid, illtpvar(3), iotopo_s)
            call check_ok('Error variable topo write', fterr)
            istatus = nf90_put_var(ncid, illtpvar(4), iomask_s)
            call check_ok('Error variable mask write', fterr)
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
            call check_ok('Error variable iy write', fterr)
            istatus = nf90_put_var(ncid, ivvar(2), xjx(1:o_nj))
            call check_ok('Error variable jx write', fterr)
            istatus = nf90_put_var(ncid, illtpvar(1), ioxlat)
            call check_ok('Error variable xlat write', fterr)
            istatus = nf90_put_var(ncid, illtpvar(2), ioxlon)
            call check_ok('Error variable xlon write', fterr)
            istatus = nf90_put_var(ncid, illtpvar(3), iotopo)
            call check_ok('Error variable topo write', fterr)
            istatus = nf90_put_var(ncid, illtpvar(4), iomask)
            call check_ok('Error variable mask write', fterr)
          end if
          if (ctype == 'SRF' .or. ctype == 'SUB') then
            rdum1 = 10.0
            istatus = nf90_put_var(ncid, isrvvar(1), rdum1)
            call check_ok('Error variable m10 write', fterr)
            rdum1 = 2.0
            istatus = nf90_put_var(ncid, isrvvar(2), rdum1)
            call check_ok('Error variable m2 write', fterr)
            rdum2(1) = 0.0
            rdum2(2) = 1.0
            istatus = nf90_put_var(ncid, isrvvar(3), rdum2)
            call check_ok('Error variable layer write', fterr)
          end if
          if (ctype == 'CHE') then
            istatus = nf90_put_var(ncid, ichname, chtrname)
            call check_ok('Error variable chtrname write', fterr)
            istatus = nf90_put_var(ncid, ichtrsol, chtrsol)
            call check_ok('Error variable chtrsol write', fterr)
            istatus = nf90_put_var(ncid, ichtrdpv, chtrdpv)
            call check_ok('Error variable chtrdpv write', fterr)
            istatus = nf90_put_var(ncid, idubinsiz, dustbsiz)
            call check_ok('Error variable dustbsiz write', fterr)
          end if

          if ( debug_level > 2 ) then
            istatus = nf90_sync(ncid)
            call check_ok('Error initial sync', fterr)
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

        subroutine addvara(ncid,ctype,vname,vst,vln,vuni,idims,lmiss, &
                           ivar)
          use netcdf
          implicit none
          integer , intent(in) :: ncid
          character(3) , intent(in) :: ctype
          character(len=*) , intent(in) :: vname
          character(len=*) , intent(in) :: vst , vln , vuni
          integer , dimension(5) , intent(in) :: idims
          logical , intent(in) :: lmiss
          integer , intent(out) :: ivar

          integer :: i , ndims

          ndims = 0
          do i = 1 , 5
            if (idims(i) > 0) ndims = ndims+1
          end do

          cdum = vname
          istatus = nf90_def_var(ncid, cdum, nf90_float, &
                             &   idims(1:ndims), ivar)
          call check_ok('Error adding variable '//vname, &
                        ctype//' FILE ERROR')
#ifdef NETCDF4_HDF5
          istatus = nf90_def_var_deflate(ncid, ivar, 1, 1, 9)
          call check_ok('Error setting deflate on variable '//vname, &
                        ctype//' FILE ERROR')
#endif
          cdum = vst
          istatus = nf90_put_att(ncid, ivar, 'standard_name', cdum)
          call check_ok('Error adding '//vname//' standard_name', &
                        ctype//' FILE ERROR')
          cdum = vln
          istatus = nf90_put_att(ncid, ivar, 'long_name', cdum)
          call check_ok('Error adding '//vname//' long_name', &
                        ctype//' FILE ERROR')
          cdum = vuni
          istatus = nf90_put_att(ncid, ivar, 'units', cdum)
          call check_ok('Error adding '//vname//' units', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, ivar, 'coordinates', &
                            &  'xlat xlon')
          call check_ok('Error adding '//vname//' coordinates', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, ivar, 'grid_mapping', &
                            &  'rcm_map')
          call check_ok('Error adding '//vname//' grid_mapping', &
                        ctype//' FILE ERROR')
          if (lmiss) then
            istatus = nf90_put_att(ncid, ivar, '_FillValue', &
                              &  smissval)
            call check_ok('Error adding '//vname//' coordinates', &
                          ctype//' FILE ERROR')
          end if
        end subroutine addvara

        subroutine writerec_srf(nx, ny, numbat, fbat, mask , idate)
          use netcdf
          implicit none
          integer , intent(in) :: nx , ny , numbat , idate
          real(4) , dimension(nx,ny,numbat) , intent(in) :: fbat
          integer , dimension(nnsg,iym1,jxm1) , intent(in) :: mask
          integer :: ivar
          integer :: n
          integer , dimension(4) :: istart , icount
          real(8) , dimension(2) :: xtime
          character(len=10) :: ctime
          logical :: lskip

          if (nx /= o_nj .or. ny /= o_ni) then
            write (6,*) 'Error writing record on SRF file'
            write (6,*) 'Expecting layers ', o_nj, 'x', o_ni
            write (6,*) 'Got layers       ', nx, 'x', ny
            call fatal(__FILE__,__LINE__,'DIMENSION MISMATCH')
          end if

          write (ctime, '(i10)') idate

          istart(2) = isrfrec
          istart(1) = 1
          icount(2) = 1
          icount(1) = 2
          xtime(2) = dble(idatediff(idate,idate0))
          xtime(1) = xtime(2) - batfrq
          istatus = nf90_put_var(ncsrf, isrfvar(1), xtime(2:2), &
                                 istart(2:2), icount(2:2))
          call check_ok('Error writing itime '//ctime, &
                      'SRF FILE ERROR')
          istatus = nf90_put_var(ncsrf, isrfvar(2), xtime, &
                                 istart(1:2), icount(1:2))
          call check_ok('Error writing tbnds '//ctime, &
                      'SRF FILE ERROR')

          ivar = 3
          lskip = .false.
          do n = 1 , numbat
            if (lskip) then
              lskip = .false.
              cycle
            end if
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
              call check_ok('Error writing '//srf_names(ivar)// &
                            ' at '//ctime, 'SRF FILE ERROR')
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
              call check_ok('Error writing '//srf_names(ivar)// &
                            ' at '//ctime, 'SRF FILE ERROR')
              istart(3) = 2
              istatus = nf90_put_var(ncsrf, isrfvar(ivar), & 
                              fbat(:,:,n+1), istart, icount)
              call check_ok('Error writing '//srf_names(ivar)// &
                            ' at '//ctime, 'SRF FILE ERROR')
              lskip = .true.
            else
              istart(3) = isrfrec
              istart(2) = 1
              istart(1) = 1
              icount(3) = 1
              icount(2) = o_ni
              icount(1) = o_nj
              istatus = nf90_put_var(ncsrf, isrfvar(ivar), & 
                       fbat(:,:,n), istart(1:3), icount(1:3))
              call check_ok('Error writing '//srf_names(ivar)// &
                            ' at '//ctime, 'SRF FILE ERROR')
            end if
            ivar = ivar + 1
          end do

          if ( iseaice == 1 .or. lakemod == 1 ) then
            dumio(:,:,1) = 0.0
            do n = 1 , nnsg
              dumio(:,:,1) = dumio(:,:,1) + &
                            real(transpose(mask(n,o_is:o_ie,o_js:o_je)))
            end do
            dumio(:,:,1) = dumio(:,:,1) / real(nnsg)
            istart(3) = isrfrec
            istart(2) = 1
            istart(1) = 1
            icount(3) = 1
            icount(2) = o_ni
            icount(1) = o_nj
            istatus = nf90_put_var(ncsrf, isrfvar(31), & 
                     dumio(:,:,1), istart(1:3), icount(1:3))
            call check_ok('Error writing '//srf_names(31)// &
                          ' at '//ctime, 'SRF FILE ERROR')
          end if

          if ( debug_level > 2 ) then
            istatus = nf90_sync(ncsrf)
            call check_ok('Error sync at '//ctime, 'SRF FILE ERROR')
          end if
          isrfrec = isrfrec + 1
        end subroutine writerec_srf

        subroutine writerec_sub(nx, ny, ns, nsub, fsub, idate)
          use netcdf
          implicit none
          integer , intent(in) :: nx , ny , ns , nsub , idate
          real(4) , dimension(ns*ns,nx/ns,ny/ns,nsub) , &
                    intent(in) :: fsub
          integer :: ivar
          integer :: n , nxb , nyb
          integer , dimension(4) :: istart , icount
          real(8) , dimension(1) :: xtime
          character(len=10) :: ctime
          logical :: lskip

          nxb = o_njg / nsg
          nyb = o_nig / nsg

          if (nx /= o_njg .or. ny /= o_nig .or. ns /= nsg) then
            write (6,*) 'Error writing record on SUB file'
            write (6,*) 'Expecting layers ', nsg, 'x', o_njg, 'x', o_nig
            write (6,*) 'Got layers       ', ns, 'x', nx, 'x', ny
            call fatal(__FILE__,__LINE__,'DIMENSION MISMATCH')
          end if

          write (ctime, '(i10)') idate

          istart(1) = isubrec
          icount(1) = 1
          xtime(1) = dble(idatediff(idate,idate0))
          istatus = nf90_put_var(ncsub, isubvar(1), xtime, &
                                 istart(1:1), icount(1:1))
          call check_ok('Error writing itime '//ctime, 'SUB FILE ERROR')
          ivar = 2
          lskip = .false.
          do n = 1 , nsub
            if (lskip) then
              lskip = .false.
              cycle
            end if
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
              istatus = nf90_put_var(ncsub, isubvar(ivar), &
                              subio, istart, icount)
              call check_ok('Error writing '//sub_names(ivar)// &
                            ' at '//ctime, 'SUB FILE ERROR')
            else if (ivar == ivarname_lookup('SUB', 'smw')) then
              istart(4) = isubrec
              istart(3) = 1
              istart(2) = 1
              istart(1) = 1
              icount(4) = 1
              icount(3) = 1
              icount(2) = o_nig
              icount(1) = o_njg
              istatus = nf90_put_var(ncsub, isubvar(ivar), & 
                              subio, istart, icount)
              call check_ok('Error writing '//sub_names(ivar)// &
                            ' at '//ctime, 'SUB FILE ERROR')
              istart(3) = 2
              call reorder(fsub,subio,nxb,nyb,nsg,nsub,n+1)
              istatus = nf90_put_var(ncsub, isubvar(ivar), & 
                              subio, istart, icount)
              call check_ok('Error writing '//sub_names(ivar)// &
                            ' at '//ctime, 'SUB FILE ERROR')
              lskip = .true.
            else
              istart(3) = isubrec
              istart(2) = 1
              istart(1) = 1
              icount(3) = 1
              icount(2) = o_nig
              icount(1) = o_njg
              istatus = nf90_put_var(ncsub, isubvar(ivar), & 
                       subio, istart(1:3), icount(1:3))
              call check_ok('Error writing '//sub_names(ivar)// &
                            ' at '//ctime, 'SUB FILE ERROR')
            end if
            ivar = ivar + 1
          end do
          if ( debug_level > 2 ) then
            istatus = nf90_sync(ncsub)
            call check_ok('Error sync at '//ctime, 'SUB FILE ERROR')
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

        subroutine writerec_rad(nx, ny, nz, nrad3d, nrad2d, frad3d, &
                                frad2d, ps, idate)
          use netcdf
          implicit none
          integer , intent(in) :: nx , ny , nz , nrad3d , nrad2d , idate
          real(4) , dimension(nx,ny,nz,nrad3d) , intent(in) :: frad3d
          real(4) , dimension(nx,ny,nrad2d) , intent(in) :: frad2d
          real(4) , dimension(nx,ny) , intent(in) :: ps
          integer :: ivar
          integer :: n
          integer , dimension(4) :: istart , icount
          real(8) , dimension(1) :: xtime
          character(len=10) :: ctime

          if (nx /= o_nj .or. ny /= o_ni .or. nz /= o_nz) then
            write (6,*) 'Error writing record on RAD file'
            write (6,*) 'Expecting layers ', o_nz, 'x', o_nj, 'x', o_ni
            write (6,*) 'Got layers       ', nz, 'x', nx, 'x', ny
            call fatal(__FILE__,__LINE__,'DIMENSION MISMATCH')
          end if

          write (ctime, '(i10)') idate

          istart(1) = iradrec
          icount(1) = 1
          xtime(1) = dble(idatediff(idate,idate0))
          istatus = nf90_put_var(ncrad, iradvar(1), xtime, &
                                 istart(1:1), icount(1:1))
          call check_ok('Error writing itime '//ctime, 'RAD FILE ERROR')

          istart(3) = iradrec
          istart(2) = 1
          istart(1) = 1
          icount(3) = 1
          icount(2) = o_ni
          icount(1) = o_nj
          istatus = nf90_put_var(ncrad, iradvar(2), &
                                 ps, istart(1:3), icount(1:3))
          call check_ok('Error writing ps at '//ctime, 'RAD FILE ERROR')

          ivar = 3
          do n = 1 , nrad3d
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
            call check_ok('Error writing '//rad_names(ivar)// &
                          ' at '//ctime, 'RAD FILE ERROR')
            ivar = ivar + 1
          end do
          do n = 1 , nrad2d
            istart(3) = iradrec
            istart(2) = 1
            istart(1) = 1
            icount(3) = 1
            icount(2) = o_ni
            icount(1) = o_nj
            istatus = nf90_put_var(ncrad, iradvar(ivar), & 
                     frad2d(:,:,n), istart(1:3), icount(1:3))
            call check_ok('Error writing '//rad_names(ivar)// &
                          ' at '//ctime, 'RAD FILE ERROR')
            ivar = ivar + 1
          end do

          if ( debug_level > 2 ) then
            istatus = nf90_sync(ncrad)
            call check_ok('Error sync at '//ctime, 'RAD FILE ERROR')
          end if
          iradrec = iradrec + 1
        end subroutine writerec_rad

        subroutine writerec_atm(nx, ny, nnx, nny, nz, ns, u, v, omega,  &
                                t, qv, qc, ps, rc, rnc, tgb, swt, rno,  &
                                mask, idate)
          use netcdf
          implicit none
          integer , intent(in) :: nx , nnx, nny , ny , ns , nz , idate
          real(8) , dimension(ny,nz,nx) , intent(in) :: u
          real(8) , dimension(ny,nz,nx) , intent(in) :: v
          real(8) , dimension(ny,nz,nx) , intent(in) :: omega
          real(8) , dimension(ny,nz,nx) , intent(in) :: t
          real(8) , dimension(ny,nz,nx) , intent(in) :: qv
          real(8) , dimension(ny,nz,nx) , intent(in) :: qc
          real(8) , dimension(ny,nx) , intent(in) :: ps
          real(8) , dimension(ny,nx) , intent(in) :: rc
          real(8) , dimension(ny,nx) , intent(in) :: rnc
          real(8) , dimension(ns,nny,nnx) , intent(in) :: tgb
          real(8) , dimension(ns,nny,nnx) , intent(in) :: swt
          real(8) , dimension(ns,nny,nnx) , intent(in) :: rno
          integer , dimension(ns,nny,nnx) , intent(in) :: mask
          integer :: i , j , n , ip1 , ip2 , jp1 , jp2 , k
          integer , dimension(4) :: istart , icount
          real(8) , dimension(1) :: xtime
          character(len=10) :: ctime

          if (nx < o_nj .or. ny < o_ni .or. nz > o_nz .or. &
              nnx < o_nj .or. nny < o_ni) then
            write (6,*) 'Error writing record on ATM file'
            write (6,*) 'Expecting layers ', o_nz, 'x', o_nj, 'x', o_ni
            write (6,*) 'Got layers 0     ', nz, 'x', nx, 'x', ny
            write (6,*) 'Got layers 1     ', nz, 'x', nnx, 'x', nny
            call fatal(__FILE__,__LINE__,'DIMENSION MISMATCH')
          end if

          write (ctime, '(i10)') idate

          if (.not. lmaskfill) then
            do n = 1 , ns
              atmsrfmask(n,:,:) = &
                    real(transpose(mask(n,o_is:o_ie,o_js:o_je)))
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
          xtime(1) = dble(idatediff(idate,idate0))
          istatus = nf90_put_var(ncatm, iatmvar(1), xtime, &
                                 istart(1:1), icount(1:1))
          call check_ok('Error writing itime '//ctime, 'ATM FILE ERROR')

          dumio(:,:,1) = real((transpose(ps(o_is:o_ie,o_js:o_je)) + &
                               rpt)*d_10)
          istart(3) = iatmrec
          istart(2) = 1
          istart(1) = 1
          icount(3) = 1
          icount(2) = o_ni
          icount(1) = o_nj
          istatus = nf90_put_var(ncatm, iatmvar(2), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok('Error writing ps at '//ctime, 'ATM FILE ERROR')

          istart(4) = iatmrec
          istart(3) = 1
          istart(2) = 1
          istart(1) = 1
          icount(4) = 1
          icount(3) = o_nz
          icount(2) = o_ni
          icount(1) = o_nj
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
                     u(ip2,k,jp1)+u(ip2,k,jp2))*d_rfour) / ps(ip1,jp1))
              end do
            end do
          end do
          istatus = nf90_put_var(ncatm, iatmvar(3), &
                                 dumio, istart, icount)
          call check_ok('Error writing '//atm_names(3)// &
                        ' at '//ctime, 'ATM FILE ERROR')
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
                     v(ip2,k,jp1)+v(ip2,k,jp2))*d_rfour) / ps(ip1,jp1))
              end do
            end do
          end do
          istatus = nf90_put_var(ncatm, iatmvar(4), &
                                 dumio, istart, icount)
          call check_ok('Error writing '//atm_names(4)//' at '//ctime, &
                        'ATM FILE ERROR')
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
                end if
                dumio(j,i,k) = real(omega(ip1,k,jp1))
              end do
            end do
          end do
          istatus = nf90_put_var(ncatm, iatmvar(5), &
                                 dumio, istart, icount)
          call check_ok('Error writing '//atm_names(5)//' at '//ctime, &
                        'ATM FILE ERROR')
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
                end if
                dumio(j,i,k) = real(t(ip1,k,jp1)/ps(ip1,jp1))
              end do
            end do
          end do
          istatus = nf90_put_var(ncatm, iatmvar(6), &
                                 dumio, istart, icount)
          call check_ok('Error writing '//atm_names(6)//' at '//ctime, &
                        'ATM FILE ERROR')
          dumio = 0.0
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
                end if
                if (qv(ip1,k,jp1) > dlowval) &
                  dumio(j,i,k) = real(qv(ip1,k,jp1)/ps(ip1,jp1))
              end do
            end do
          end do
          istatus = nf90_put_var(ncatm, iatmvar(7), &
                                 dumio, istart, icount)
          call check_ok('Error writing '//atm_names(7)//' at '//ctime, &
                        'ATM FILE ERROR')
          dumio = 0.0
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
                end if
                if (qc(ip1,k,jp1) > dlowval) &
                  dumio(j,i,k) = real(qc(ip1,k,jp1)/ps(ip1,jp1))
              end do
            end do
          end do
          istatus = nf90_put_var(ncatm, iatmvar(8), &
                                 dumio, istart, icount)
          call check_ok('Error writing '//atm_names(8)//' at '//ctime, &
                        'ATM FILE ERROR')

          istart(3) = iatmrec
          istart(2) = 1
          istart(1) = 1
          icount(3) = 1
          icount(2) = o_ni
          icount(1) = o_nj

          dumio(:,:,1) = 0.0
          where (transpose(rc(o_is:o_ie,o_js:o_je)) > dlowval)
            dumio(:,:,1) = real(transpose(rc(o_is:o_ie,o_js:o_je)))
          end where
          where (transpose(rnc(o_is:o_ie,o_js:o_je)) > dlowval)
            dumio(:,:,1) = dumio(:,:,1) + &
                           real(transpose(rnc(o_is:o_ie,o_js:o_je)))
          end where
          dumio(:,:,1) = dumio(:,:,1)*real(tpd)
          istatus = nf90_put_var(ncatm, iatmvar(9), & 
                     dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok('Error writing '//atm_names(9)//' at '//ctime, &
                        'ATM FILE ERROR')
          dumio(:,:,1) = real(transpose(sum(tgb(:,o_is:o_ie,o_js:o_je), &
                                       dim=1)*xns2r))
          istatus = nf90_put_var(ncatm, iatmvar(10), & 
                     dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok('Error writing '//atm_names(10)//' at '//ctime, &
                        'ATM FILE ERROR')
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
          istatus = nf90_put_var(ncatm, iatmvar(11), & 
                     dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok('Error writing '//atm_names(11)//' at '//ctime, &
                        'ATM FILE ERROR')
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
          istatus = nf90_put_var(ncatm, iatmvar(12), & 
                     dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok('Error writing '//atm_names(12)//' at '//ctime, &
                        'ATM FILE ERROR')

          if ( debug_level > 2 ) then
            istatus = nf90_sync(ncatm)
            call check_ok('Error sync at '//ctime, 'ATM FILE ERROR')
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
          integer , intent(in) :: nx , ny , nnx , nny , nz , nt , idate
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
          real(8) , dimension(1) :: xtime
          character(len=10) :: ctime

          if (nx < o_nj .or. ny < o_ni .or. nz > o_nz) then
            write (6,*) 'Error writing record on CHE file'
            write (6,*) 'Expecting layers ', o_nz, 'x', o_nj, 'x', o_ni
            write (6,*) 'Got layers       ', nz, 'x', nx, 'x', ny
            call fatal(__FILE__,__LINE__,'DIMENSION MISMATCH')
          end if

          write (ctime, '(i10)') idate

          istart(1) = icherec
          icount(1) = 1
          xtime(1) = dble(idatediff(idate,idate0))
          istatus = nf90_put_var(ncche, ichevar(1), xtime, &
                                 istart(1:1), icount(1:1))
          call check_ok('Error writing itime '//ctime, 'CHE FILE ERROR')

          dumio(:,:,1) = real(transpose(ps(o_is:o_ie,o_js:o_je)+rpt)*d_10)
          istart(3) = icherec
          istart(2) = 1
          istart(1) = 1
          icount(3) = 1
          icount(2) = o_ni
          icount(1) = o_nj
          istatus = nf90_put_var(ncche, ichevar(2), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok('Error writing ps at '//ctime, 'CHE FILE ERROR')

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
            istatus = nf90_put_var(ncche, ichevar(3), &
                                 dumio, istart, icount)
            call check_ok('Error writing '//che_names(3)//' at '//ctime,&
                          'CHE FILE ERROR')
          end do

          istart(4) = icherec
          istart(3) = 1
          istart(2) = 1
          istart(1) = 1
          icount(4) = 1
          icount(3) = o_nz
          icount(2) = o_ni
          icount(1) = o_nj
          do k = 1 , nz
            dumio(:,:,k) = real(transpose(aerext(o_is:o_ie,k,o_js:o_je)))
          end do
          istatus = nf90_put_var(ncche, ichevar(4), &
                                 dumio, istart(1:4), icount(1:4))
          call check_ok('Error writing '//che_names(4)//' at '//ctime, &
                        'CHE FILE ERROR')
          do k = 1 , nz
            dumio(:,:,k) = real(transpose(aerssa(o_is:o_ie,k,o_js:o_je)))
          end do
          istatus = nf90_put_var(ncche, ichevar(5), &
                                 dumio, istart(1:4), icount(1:4))
          call check_ok('Error writing '//che_names(5)//' at '//ctime, &
                        'CHE FILE ERROR')
          do k = 1 , nz
            dumio(:,:,k) = real(transpose(aerasp(o_is:o_ie,k,o_js:o_je)))
          end do
          istatus = nf90_put_var(ncche, ichevar(6), &
                                 dumio, istart(1:4), icount(1:4))
          call check_ok('Error writing '//che_names(6)//' at '//ctime, &
                        'CHE FILE ERROR')
          do n = 1 , nt
            istart(4) = icherec
            istart(3) = n
            istart(2) = 1
            istart(1) = 1
            icount(4) = 1
            icount(3) = 1
            icount(2) = o_ni
            icount(1) = o_nj
            dumio(:,:,1) = real(transpose(dtrace(o_is:o_ie,o_js:o_je,n)))
            istatus = nf90_put_var(ncche, ichevar(7), &
                                 dumio(:,:,1), istart(1:4), icount(1:4))
            call check_ok('Error writing '//che_names(7)//' at '//ctime,&
                          'CHE FILE ERROR')
            dumio(:,:,1) = real(transpose(wdlsc(o_is:o_ie,o_js:o_je,n))*cfd)
            istatus = nf90_put_var(ncche, ichevar(8), &
                                 dumio(:,:,1), istart(1:4), icount(1:4))
            call check_ok('Error writing '//che_names(8)//' at '//ctime,&
                          'CHE FILE ERROR')
            dumio(:,:,1) = real(transpose(wdcvc(o_is:o_ie,o_js:o_je,n))*cfd)
            istatus = nf90_put_var(ncche, ichevar(9), &
                                 dumio(:,:,1), istart(1:4), icount(1:4))
            call check_ok('Error writing '//che_names(9)//' at '//ctime,&
                          'CHE FILE ERROR')
            dumio(:,:,1) = real(transpose(ddsfc(o_is:o_ie,o_js:o_je,n))*cfd)
            istatus = nf90_put_var(ncche, ichevar(10), &
                                 dumio(:,:,1), istart(1:4), icount(1:4))
            call check_ok('Error writing '//che_names(10)// &
                          ' at '//ctime, 'CHE FILE ERROR')
            dumio(:,:,1) = real(transpose(wxsg(o_is:o_ie,o_js:o_je,n))*cfd)
            istatus = nf90_put_var(ncche, ichevar(11), &
                                 dumio(:,:,1), istart(1:4), icount(1:4))
            call check_ok('Error writing '//che_names(11)// &
                          ' at '//ctime, 'CHE FILE ERROR')
            dumio(:,:,1) = real(transpose(wxaq(o_is:o_ie,o_js:o_je,n))*cfd)
            istatus = nf90_put_var(ncche, ichevar(12), &
                                 dumio(:,:,1), istart(1:4), icount(1:4))
            call check_ok('Error writing '//che_names(12)// &
                          ' at '//ctime, 'CHE FILE ERROR')
            dumio(:,:,1) = real(transpose(cemtrac(o_is:o_ie,o_js:o_je,n))*cfd)
            istatus = nf90_put_var(ncche, ichevar(13), &
                                 dumio(:,:,1), istart(1:4), icount(1:4))
            call check_ok('Error writing '//che_names(13)// &
                          ' at '//ctime, 'CHE FILE ERROR')
          end do
          istart(3) = icherec
          istart(2) = 1
          istart(1) = 1
          icount(3) = 1
          icount(2) = o_ni
          icount(1) = o_nj
          dumio(:,:,1) = real(transpose(aertarf(o_is:o_ie,o_js:o_je)))
          istatus = nf90_put_var(ncche, ichevar(14), &
                               dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok('Error writing '//che_names(14)//' at '//ctime, &
                        'CHE FILE ERROR')
          dumio(:,:,1) = real(transpose(aersrrf(o_is:o_ie,o_js:o_je)))
          istatus = nf90_put_var(ncche, ichevar(15), &
                               dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok('Error writing '//che_names(15)//' at '//ctime, &
                        'CHE FILE ERROR')
          dumio(:,:,1) = real(transpose(aertalwrf(o_is:o_ie,o_js:o_je)))
          istatus = nf90_put_var(ncche, ichevar(16), &
                               dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok('Error writing '//che_names(16)//' at '//ctime, &
                        'CHE FILE ERROR')
          dumio(:,:,1) = real(transpose(aersrlwrf(o_is:o_ie,o_js:o_je)))
          istatus = nf90_put_var(ncche, ichevar(17), &
                               dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok('Error writing '//che_names(17)//' at '//ctime, &
                        'CHE FILE ERROR')

          if ( debug_level > 2 ) then
            istatus = nf90_sync(ncche)
            call check_ok('Error sync at '//ctime, 'CHE FILE ERROR')
          end if
          icherec = icherec + 1
        end subroutine writerec_che

        subroutine writerec_lak(nx, ny, numbat, fbat, evl, aveice, &
                                hsnow, tlake, idate)
          use netcdf
          implicit none
          integer , intent(in) :: nx , ny , numbat , idate
          real(4) , dimension(nx,ny,numbat) , intent(in) :: fbat
          real(8) , dimension(nnsg,iym1,jx) , intent(in) :: evl
          real(8) , dimension(nnsg,iym1,jx) , intent(in) :: aveice
          real(8) , dimension(nnsg,iym1,jx) , intent(in) :: hsnow
          real(8) , dimension(ndpmax,nnsg,iym1,jx) , intent(in) :: tlake
          integer :: ivar
          integer :: n
          integer , dimension(4) :: istart , icount
          real(8) , dimension(1) :: xtime
          character(len=10) :: ctime

          if (nx /= o_nj .or. ny /= o_ni) then
            write (6,*) 'Error writing record on LAK file'
            write (6,*) 'Expecting layers ', o_nj, 'x', o_ni
            write (6,*) 'Got layers       ', nx, 'x', ny
            call fatal(__FILE__,__LINE__,'DIMENSION MISMATCH')
          end if

          write (ctime, '(i10)') idate

          istart(1) = ilakrec
          icount(1) = 1
          xtime(1) = dble(idatediff(idate,idate0))
          istatus = nf90_put_var(nclak, ilakvar(1), xtime, &
                                 istart(1:1), icount(1:1))
          call check_ok('Error writing itime '//ctime, &
                      'LAK FILE ERROR')

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
            call check_ok('Error writing '//lak_names(ivar)// &
                          ' at '//ctime, 'LAK FILE ERROR')
            ivar = ivar + 1
          end do

          ! Add lake model output
          dumio(:,:,1) =  &
               real(transpose(sum(evl(:,o_is:o_ie,o_js:o_je),1))/ &
               dble(nnsg))
          istatus = nf90_put_var(nclak, ilakvar(ivar), & 
                   dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok('Error writing '//lak_names(ivar)// &
                        ' at '//ctime, 'LAK FILE ERROR')
          ivar = ivar + 1
          dumio(:,:,1) =  &
               real(transpose(sum(aveice(:,o_is:o_ie,o_js:o_je),1))/ &
               dble(nnsg))
          istatus = nf90_put_var(nclak, ilakvar(ivar), & 
                   dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok('Error writing '//lak_names(ivar)// &
                        ' at '//ctime, 'LAK FILE ERROR')
          ivar = ivar + 1
          dumio(:,:,1) =  &
               real(transpose(sum(hsnow(:,o_is:o_ie,o_js:o_je),1))/ &
               dble(nnsg))
          istatus = nf90_put_var(nclak, ilakvar(ivar), & 
                   dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok('Error writing '//lak_names(ivar)// &
                        ' at '//ctime, 'LAK FILE ERROR')
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
            dumio(:,:,1) =  &
               real(transpose(sum(tlake(n,:,o_is:o_ie,o_js:o_je),1))/ &
               dble(nnsg))
            where (iolnds == 14)
              dumio(:,:,1) = dumio(:,:,1) + real(tzero)
            end where
            istatus = nf90_put_var(nclak, ilakvar(ivar), & 
                   dumio(:,:,1), istart, icount)
            call check_ok('Error writing '//lak_names(ivar)// &
                        ' at '//ctime, 'LAK FILE ERROR')
          end do

          if ( debug_level > 2 ) then
            istatus = nf90_sync(nclak)
            call check_ok('Error sync at '//ctime, 'LAK FILE ERROR')
          end if
          ilakrec = ilakrec + 1
        end subroutine writerec_lak

        subroutine check_ok(m1,mf)
          use netcdf
          implicit none
          character(*) :: m1 , mf
          if (istatus /= nf90_noerr) then 
            write (6,*) trim(m1)
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__,trim(mf))
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
          if (allocated(ioxlat)) deallocate(ioxlat)
          if (allocated(ioxlon)) deallocate(ioxlon)
          if (allocated(iotopo)) deallocate(iotopo)
          if (allocated(iomask)) deallocate(iomask)
          if (allocated(iolnds)) deallocate(iolnds)
          if (allocated(hsigma)) deallocate(hsigma)
          if (allocated(ioxlat_s)) deallocate(ioxlat_s)
          if (allocated(ioxlon_s)) deallocate(ioxlon_s)
          if (allocated(iotopo_s)) deallocate(iotopo_s)
          if (allocated(iomask_s)) deallocate(iomask_s)
          if (allocated(atmsrfmask)) deallocate(atmsrfmask)
          if (allocated(atmsrfsum)) deallocate(atmsrfsum)
          if (allocated(dumio)) deallocate(dumio)
          if (allocated(sp2d)) deallocate(sp2d)
          if (allocated(subio)) deallocate(subio)
          if (allocated(sp2d1)) deallocate(sp2d1)

        end subroutine release_mod_ncio

      end module mod_ncio
