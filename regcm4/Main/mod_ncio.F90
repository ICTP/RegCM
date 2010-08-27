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

        integer , private :: idmin , isdmin , ibcin , ncatm , ncsrf , &
                             ncsub , ncrad , ncche
        integer , private :: istatus
        integer , private :: ibcrec , ibcnrec
        integer , private :: iatmrec , isrfrec , isubrec , iradrec , &
                             icherec
        integer , dimension(12) :: iatmvar
        integer , dimension(27) :: isrfvar
        integer , dimension(16) :: isubvar
        integer , dimension(15) :: iradvar
        integer , dimension(17) :: ichevar
        character(256) , private :: dname , sdname , aername , icbcname
        integer , dimension(:) , allocatable , private :: icbc_idate
        real(4) , dimension(:) , allocatable , private :: hsigma
        integer , dimension(7) , private :: icbc_ivar
        logical , private :: lso4p

        ! DIM1 is iy ,   DIM2 is jx , DIM3 is time ,       DIM4 is kz
        ! DIM5 is m10 ,  DIM6 is m2 , DIM7 is soil_layer , DIM8 is nv
        ! DIM9 is tracer
        integer , dimension(9) , private :: idims

        data lso4p   /.false./
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

        real(4) , dimension(:,:) , allocatable , private :: ioxlat
        real(4) , dimension(:,:) , allocatable , private :: ioxlon
        real(4) , dimension(:,:) , allocatable , private :: iotopo
        real(4) , dimension(:,:) , allocatable , private :: ioxlat_s
        real(4) , dimension(:,:) , allocatable , private :: ioxlon_s
        real(4) , dimension(:,:) , allocatable , private :: iotopo_s

      contains

        subroutine init_mod_ncio
          use mod_dynparam
          implicit none
          character(3) :: sbstring
          dname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
          write (sbstring,'(i0.3)') nsg
          sdname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN'// &
              &    sbstring//'.nc'
          aername = trim(dirglob)//pthsep//trim(domname)//'_AERO.nc'
          icbcname = trim(dirglob)//pthsep//trim(domname)//'_ICBC.'// &
              &      'YYYYMMDDHH.nc'
          allocate(ioxlat(jx,iy))
          allocate(ioxlon(jx,iy))
          allocate(iotopo(jx,iy))
          allocate(hsigma(kz))
          if (nsg > 1) then
            allocate(ioxlat_s(jxsg,iysg))
            allocate(ioxlon_s(jxsg,iysg))
            allocate(iotopo_s(jxsg,iysg))
          end if
        end subroutine init_mod_ncio

        subroutine open_domain(r8pt , dx , sigma)
          use mod_dynparam
          use mod_message
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

          write (aline,*) 'READING HEADER FILE:', dname
          call say
          istatus = nf90_open(dname, nf90_nowrite, idmin)
          call check_ok(istatus, &
                        'Error Opening Domain file '//trim(dname), &
                        'CANNOT OPEN DOMAIN FILE')
          if ( nsg.gt.1 ) then
            write (aline,*) 'READING HEADER SUBDOMAIN FILE:', sdname
            call say
            istatus = nf90_open(sdname, nf90_nowrite, isdmin)
            call check_ok(istatus, &
                          'Error Opening SubDomain file '//trim(sdname),&
                          'CANNOT OPEN SUBDOM FILE')
          end if
          istatus = nf90_inq_dimid(idmin, 'iy', idimid)
          call check_ok(istatus, 'Dimension iy missing', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_inquire_dimension(idmin, idimid, len=iyy)
          call check_ok(istatus, 'Dimension iy read error', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_inq_dimid(idmin, 'jx', idimid)
          call check_ok(istatus, 'Dimension jx missing', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_inquire_dimension(idmin, idimid, len=jxx)
          call check_ok(istatus, 'Dimension jx read error', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_inq_dimid(idmin, 'kz', idimid)
          call check_ok(istatus, 'Dimension kz missing', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_inquire_dimension(idmin, idimid, len=kzz)
          call check_ok(istatus, 'Dimension kz read error', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_inq_varid(idmin, 'ptop', ivarid)
          call check_ok(istatus, 'Variable ptop missing', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_get_var(idmin, ivarid, ptsp)
          call check_ok(istatus, 'Variable ptop read error', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_get_att(idmin, nf90_global, 'projection', proj)
          call check_ok(istatus, 'Attribute projection missing', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_get_att(idmin, nf90_global, &
                       &         'grid_size_in_meters', dsx)
          call check_ok(istatus, 'Attribute gridsize missing', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_get_att(idmin, nf90_global, &
                     &         'latitude_of_projection_origin', iclat)
          call check_ok(istatus, 'Attribute clat missing', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_get_att(idmin, nf90_global, &
                      &         'longitude_of_projection_origin', iclon)
          call check_ok(istatus, 'Attribute clon missing', &
                        'DOMAIN FILE ERROR')
!
!         Consistency Check
!
          if ( iyy.ne.iy .or. jxx.ne.jx .or. kzz.ne.kzp1 ) then
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
          if (abs(dble(ptsp/10.0)-dble(ptop)) > 0.001D+00) then
            write (6,*) 'Error: ptop from regcm.in and DOMAIN file ', &
                        'differ.'
            write (6,*) 'Input namelist = ', ptop
            write (6,*) 'DOMAIN file    = ', ptsp/10.0
            call fatal(__FILE__,__LINE__, 'DOMAIN ptop ERROR')
          end if
          if (proj /= iproj) then
            write (6,*) 'Error: proj from regcm.in and DOMAIN file ', &
                        'differ.'
            write (6,*) 'Input namelist = ', iproj
            write (6,*) 'DOMAIN file    = ', proj
            call fatal(__FILE__,__LINE__, 'DOMAIN proj ERROR')
          end if
          if (abs(dble(dsx/1000.0)-dble(ds)) > 0.001D+00) then
            write (6,*) 'Error: ds from regcm.in and DOMAIN file ', &
                        'differ.'
            write (6,*) 'Input namelist = ', ds
            write (6,*) 'DOMAIN file    = ', dsx/1000.0
            call fatal(__FILE__,__LINE__, 'DOMAIN ds ERROR')
          end if
          if (abs(dble(iclat)-dble(clat)) > 0.001D+00) then
            write (6,*) 'Error: clat from regcm.in and DOMAIN file ', &
                        'differ.'
            write (6,*) 'Input namelist = ', clat
            write (6,*) 'DOMAIN file    = ', iclat
            call fatal(__FILE__,__LINE__, 'DOMAIN clat ERROR')
          end if
          if (abs(dble(iclon)-dble(clon)) > 0.001D+00) then
            write (6,*) 'Error: clon from regcm.in and DOMAIN file ', &
                        'differ.'
            write (6,*) 'Input namelist = ', clon
            write (6,*) 'DOMAIN file    = ', iclon
            call fatal(__FILE__,__LINE__, 'DOMAIN clon ERROR')
          end if
!
!         Assign values in the top data modules
!
          r8pt = ptsp/10.0
          dx = dsx
          istatus = nf90_inq_varid(idmin, 'sigma', ivarid)
          call check_ok(istatus, 'Variable sigma missing', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_get_var(idmin, ivarid, rsdum)
          call check_ok(istatus, 'Variable sigma read error', &
                        'DOMAIN FILE ERROR')
          sigma = dble(rsdum)
          do k = 1 , kz
            hsigma(k) = (sigma(k)+sigma(k+1))/2.0
          end do

        end subroutine open_domain

        subroutine read_domain(ht,htsd,lnd,xlat,xlon,xmap,dmap,f,snw)
          use mod_dynparam
          use mod_message
          use netcdf
          implicit none

          real(8) , dimension(iy,jx) , intent(out) :: ht
          real(8) , dimension(iy,jx) , intent(out) :: htsd
          real(8) , dimension(iy,jx) , intent(out) :: lnd
          real(8) , dimension(iy,jx) , intent(out) :: xlat
          real(8) , dimension(iy,jx) , intent(out) :: xlon
          real(8) , dimension(iy,jx) , intent(out) :: xmap
          real(8) , dimension(iy,jx) , intent(out) :: dmap
          real(8) , dimension(iy,jx) , intent(out) :: f
          real(8) , dimension(nnsg,iy,jx) , intent(out) :: snw

          integer :: ivarid , n
          real(4) , dimension(jx,iy) :: sp2d

          if (idmin < 0) then
            write (6,*) 'Error : Domain file not in open state'
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if

          istatus = nf90_inq_varid(idmin, 'topo', ivarid)
          call check_ok(istatus, 'Variable topo missing', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_get_var(idmin, ivarid, iotopo)
          call check_ok(istatus, 'Variable topo read error', &
                        'DOMAIN FILE ERROR')
          ht = transpose(iotopo)
          istatus = nf90_inq_varid(idmin, 'htsd', ivarid)
          call check_ok(istatus, 'Variable htsd missing', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_get_var(idmin, ivarid, sp2d)
          call check_ok(istatus, 'Variable htsd read error', &
                        'DOMAIN FILE ERROR')
          htsd = transpose(sp2d)
          istatus = nf90_inq_varid(idmin, 'landuse', ivarid)
          call check_ok(istatus, 'Variable landuse missing', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_get_var(idmin, ivarid, sp2d)
          call check_ok(istatus, 'Variable landuse read error', &
                        'DOMAIN FILE ERROR')
          lnd = transpose(sp2d)
          istatus = nf90_inq_varid(idmin, 'xlat', ivarid)
          call check_ok(istatus, 'Variable xlat missing', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_get_var(idmin, ivarid, ioxlat)
          call check_ok(istatus, 'Variable xlat read error', &
                        'DOMAIN FILE ERROR')
          xlat = transpose(ioxlat)
          istatus = nf90_inq_varid(idmin, 'xlon', ivarid)
          call check_ok(istatus, 'Variable xlon missing', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_get_var(idmin, ivarid, ioxlon)
          call check_ok(istatus, 'Variable xlon read error', &
                        'DOMAIN FILE ERROR')
          xlon = transpose(ioxlon)
          istatus = nf90_inq_varid(idmin, 'xmap', ivarid)
          call check_ok(istatus, 'Variable xmap missing', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_get_var(idmin, ivarid, sp2d)
          call check_ok(istatus, 'Variable xmap read error', &
                        'DOMAIN FILE ERROR')
          xmap = transpose(sp2d)
          istatus = nf90_inq_varid(idmin, 'dmap', ivarid)
          call check_ok(istatus, 'Variable dmap missing', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_get_var(idmin, ivarid, sp2d)
          call check_ok(istatus, 'Variable dmap read error', &
                        'DOMAIN FILE ERROR')
          dmap = transpose(sp2d)
          istatus = nf90_inq_varid(idmin, 'coriol', ivarid)
          call check_ok(istatus, 'Variable coriol missing', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_get_var(idmin, ivarid, sp2d)
          call check_ok(istatus, 'Variable coriol read error', &
                        'DOMAIN FILE ERROR')
          f = transpose(sp2d)
          istatus = nf90_inq_varid(idmin, 'snowam', ivarid)
          call check_ok(istatus, 'Variable snowam missing', &
                        'DOMAIN FILE ERROR')
          istatus = nf90_get_var(idmin, ivarid, sp2d)
          call check_ok(istatus, 'Variable snowam read error', &
                        'DOMAIN FILE ERROR')
          do n = 1 , nnsg
            snw(n,:,:) = transpose(sp2d)
          end do
        end subroutine read_domain

        subroutine read_subdomain(ht1,lnd1,xlat1,xlon1)
          use mod_dynparam
          use mod_message
          use mod_constants , only : gti
          use netcdf
          implicit none

          real(8) , dimension(nnsg,iy,jx) , intent(out) :: ht1
          real(8) , dimension(nnsg,iy,jx) , intent(out) :: lnd1
          real(8) , dimension(nnsg,iy,jx) , intent(out) :: xlat1
          real(8) , dimension(nnsg,iy,jx) , intent(out) :: xlon1

          integer :: ivarid
          integer :: i , j , n , ii , jj
          real(4) , dimension(jxsg,iysg) :: sp2d1
          
          if (isdmin < 0) then
            write (6,*) 'Error : Subdom file not in open state'
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if

          istatus = nf90_inq_varid(isdmin, 'topo', ivarid)
          call check_ok(istatus, 'Variable topo missing', &
                        'SUBDOMAIN FILE ERROR')
          istatus = nf90_get_var(isdmin, ivarid, iotopo_s)
          call check_ok(istatus, 'Variable topo read error', &
                        'SUBDOMAIN FILE ERROR')
          do j = 1 , jxsg
            do i = 1 , iysg
              jj = mod(j,nsg)
              if ( jj.eq.0 ) jj = nsg
              ii = mod(i,nsg)
              if ( ii.eq.0 ) ii = nsg
              n = (jj-1)*nsg + ii
              jj = (j+nsg-1)/nsg
              ii = (i+nsg-1)/nsg
              ht1(n,ii,jj) = iotopo_s(j,i)*gti
            end do
          end do
          istatus = nf90_inq_varid(isdmin, 'landuse', ivarid)
          call check_ok(istatus, 'Variable landuse missing', &
                        'SUBDOMAIN FILE ERROR')
          istatus = nf90_get_var(isdmin, ivarid, sp2d1)
          call check_ok(istatus, 'Variable landuse read error', &
                        'SUBDOMAIN FILE ERROR')
          do j = 1 , jxsg
            do i = 1 , iysg
              jj = mod(j,nsg)
              if ( jj.eq.0 ) jj = nsg
              ii = mod(i,nsg)
              if ( ii.eq.0 ) ii = nsg
              n = (jj-1)*nsg + ii
              jj = (j+nsg-1)/nsg
              ii = (i+nsg-1)/nsg
              lnd1(n,ii,jj) = sp2d1(j,i)
            end do
          end do
          istatus = nf90_inq_varid(isdmin, 'xlat', ivarid)
          call check_ok(istatus, 'Variable xlat missing', &
                        'SUBDOMAIN FILE ERROR')
          istatus = nf90_get_var(isdmin, ivarid, ioxlat_s)
          call check_ok(istatus, 'Variable xlat read error', &
                        'SUBDOMAIN FILE ERROR')
          do j = 1 , jxsg
            do i = 1 , iysg
              jj = mod(j,nsg)
              if ( jj.eq.0 ) jj = nsg
              ii = mod(i,nsg)
              if ( ii.eq.0 ) ii = nsg
              n = (jj-1)*nsg + ii
              jj = (j+nsg-1)/nsg
              ii = (i+nsg-1)/nsg
              xlat1(n,ii,jj) = ioxlat_s(j,i)
            end do
          end do
          istatus = nf90_inq_varid(isdmin, 'xlon', ivarid)
          call check_ok(istatus, 'Variable xlon missing', &
                        'SUBDOMAIN FILE ERROR')
          istatus = nf90_get_var(isdmin, ivarid, ioxlon_s)
          call check_ok(istatus, 'Variable xlon read error', &
                        'SUBDOMAIN FILE ERROR')
          do j = 1 , jxsg
            do i = 1 , iysg
              jj = mod(j,nsg)
              if ( jj.eq.0 ) jj = nsg
              ii = mod(i,nsg)
              if ( ii.eq.0 ) ii = nsg
              n = (jj-1)*nsg + ii
              jj = (j+nsg-1)/nsg
              ii = (i+nsg-1)/nsg
              xlon1(n,ii,jj) = ioxlon_s(j,i)
            end do
          end do

        end subroutine read_subdomain

        subroutine close_domain
          use mod_dynparam
          use netcdf
          implicit none

          if (idmin >= 0) then
            istatus = nf90_close(idmin)
            call check_ok(istatus, 'Domain file close error', &
                        'DOMAIN FILE ERROR')
            idmin = -1
          end if
          if ( nsg>1 .and. isdmin >=0 ) then
            istatus = nf90_close(isdmin)
            call check_ok(istatus, 'SubDomain file close error', &
                        'SUBDOMAIN FILE ERROR')
            isdmin = -1
          end if

        end subroutine close_domain

        subroutine read_texture(nats,texture)
          use mod_dynparam
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
            call check_ok(istatus, &
                       &  'Error Opening Domain file '//trim(dname), &
                       &  'DOMAIN FILE OPEN ERROR')
          end if
          istatus = nf90_inq_varid(idmin, 'texture_fraction', ivarid)
          call check_ok(istatus, 'Variable texture_fraction missing', &
                     &  'DOMAIN FILE ERROR')
          istart(2) = 1
          istart(1) = 1
          icount(3) = 1
          icount(2) = iy
          icount(1) = jx
          do n = 1 , nats
            istart(3) = n
            istatus = nf90_get_var(idmin, ivarid, toto, istart, icount)
            call check_ok(istatus, 'Variable texture_frac read error', &
                     &  'DOMAIN FILE ERROR')
            do j = 1 , jx
              do i = 1 , iy
                texture(i,j,n) = dble(toto(j,i))*0.01
                if (texture(i,j,n)<0.) texture(i,j,n)=0.
              end do
            end do
          end do
          call close_domain

        end subroutine read_texture

        subroutine read_aerosol(chtrname,chemsrc)
          use mod_dynparam
          use mod_message
          use netcdf
          implicit none

          character(5) , dimension(ntr) , intent(in) :: chtrname
          real(8) , dimension(iy,jx,12,ntr) , intent(out) :: chemsrc

          integer :: ncid , ivarid , istatus
          real(4) , dimension(jx,iy) :: toto
          character(5) :: aerctl
          integer , dimension(3) :: istart , icount
          integer :: itr , i , j , m

          istatus = nf90_open(aername, nf90_nowrite, ncid)
          call check_ok(istatus, &
                     &  'Error Opening Aerosol file '//trim(aername), &
                     &  'AEROSOL FILE OPEN ERROR')

          do itr = 1 , ntr
            aerctl = chtrname(itr)
            write (aline, *) itr , aerctl
            call say
            if ( aerctl(1:4).ne.'DUST') then
              if ( aerctl(1:3).eq.'SO2' ) then
                if ( aertyp(4:4).eq.'1' ) then
                  istatus = nf90_inq_varid(ncid, 'so2', ivarid)
                  call check_ok(istatus, 'Variable so2 missing', &
                            &   'AEROSOL FILE ERROR')
                  istatus = nf90_get_var(ncid, ivarid, toto)
                  call check_ok(istatus, 'Variable so2 read error', &
                            &   'AEROSOL FILE ERROR')
                  do m = 1 , 12
                    do j = 1 , jx
                      do i = 1 , iy
                        chemsrc(i,j,m,itr) = toto(j,i)
                      end do
                    end do
                  end do
                end if
                if ( aertyp(5:5).eq.'1' ) then
                  istatus = nf90_inq_varid(ncid, 'so2_monthly', ivarid)
                  call check_ok(istatus, 'Variable so2_mon missing', &
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
                    call check_ok(istatus, 'Variable so2_mon read err', &
                            &   'AEROSOL FILE ERROR')
                    do j = 1 , jx
                      do i = 1 , iy
                        chemsrc(i,j,m,itr) = chemsrc(i,j,m,itr) + &
                                              & toto(j,i)
                      end do
                    end do
                  end do
                end if
              else if ( aerctl(1:2).eq.'BC' ) then
                if ( aertyp(4:4).eq.'1' ) then
                  istatus = nf90_inq_varid(ncid, 'bc', ivarid)
                  call check_ok(istatus, 'Variable bc missing', &
                            &   'AEROSOL FILE ERROR')
                  istatus = nf90_get_var(ncid, ivarid, toto)
                  call check_ok(istatus, 'Variable bc read error', &
                            &   'AEROSOL FILE ERROR')
                  do m = 1 , 12
                    do j = 1 , jx
                      do i = 1 , iy
                        chemsrc(i,j,m,itr) = toto(j,i)
                      end do
                    end do
                  end do
                end if
                if ( aertyp(5:5).eq.'1' ) then
                  istatus = nf90_inq_varid(ncid, 'bc_monthly', ivarid)
                  call check_ok(istatus, 'Variable bc_mon missing', &
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
                    call check_ok(istatus, 'Variable bc_mon read err', &
                            &   'AEROSOL FILE ERROR')
                    do j = 1 , jx
                      do i = 1 , iy
                        chemsrc(i,j,m,itr) = chemsrc(i,j,m,itr) + &
                                       &        toto(j,i)
                      end do
                    end do
                  end do
                end if
              else if ( aerctl(1:2).eq.'OC' ) then
                if ( aertyp(4:4).eq.'1' ) then
                  istatus = nf90_inq_varid(ncid, 'oc', ivarid)
                  call check_ok(istatus, 'Variable oc missing', &
                            &   'AEROSOL FILE ERROR')
                  istatus = nf90_get_var(ncid, ivarid, toto)
                  call check_ok(istatus, 'Variable oc read error', &
                            &   'AEROSOL FILE ERROR')
                  do m = 1 , 12
                    do j = 1 , jx
                      do i = 1 , iy
                        chemsrc(i,j,m,itr) = toto(j,i)
                      end do
                    end do
                  end do
                end if
                if ( aertyp(5:5).eq.'1' ) then
                  istatus = nf90_inq_varid(ncid, 'oc_monthly', ivarid)
                  call check_ok(istatus, 'Variable oc_mon missing', &
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
                    call check_ok(istatus, 'Variable oc_mon read err', &
                            &   'AEROSOL FILE ERROR')
                    do j = 1 , jx
                      do i = 1 , iy
                        chemsrc(i,j,m,itr) = chemsrc(i,j,m,itr) + &
                                              & toto(j,i)
                      end do
                    end do
                  end do
                end if
              end if
            end if
          end do

          istatus = nf90_close(ncid)
          call check_ok(istatus, &
                      & 'Error Closing Aerosol file '//trim(aername), &
                      &   'AEROSOL FILE CLOSE ERROR')

        end subroutine read_aerosol

        subroutine open_icbc(idate)
          use mod_dynparam
          use mod_message
          use netcdf
          use mod_date , only : timeval2idate
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
          call check_ok(istatus, &
                     &  'Error Opening ICBC file '//trim(icbcname), &
                     &  'ICBC FILE OPEN ERROR')
          ibcrec = 1
          ibcnrec = 0
          istatus = nf90_inq_dimid(ibcin, 'iy', idimid)
          call check_ok(istatus, 'Dimension iy missing', &
                     &  'ICBC FILE ERROR')
          istatus = nf90_inquire_dimension(ibcin, idimid, len=iyy)
          call check_ok(istatus, 'Dimension iy read error', &
                     &  'ICBC FILE ERROR')
          istatus = nf90_inq_dimid(ibcin, 'jx', idimid)
          call check_ok(istatus, 'Dimension jx missing', &
                     &  'ICBC FILE ERROR')
          istatus = nf90_inquire_dimension(ibcin, idimid, len=jxx)
          call check_ok(istatus, 'Dimension jx read error', &
                     &  'ICBC FILE ERROR')
          istatus = nf90_inq_dimid(ibcin, 'kz', idimid)
          call check_ok(istatus, 'Dimension kz missing', &
                     &  'ICBC FILE ERROR')
          istatus = nf90_inquire_dimension(ibcin, idimid, len=kzz)
          call check_ok(istatus, 'Dimension kz read error', &
                     &  'ICBC FILE ERROR')
          if ( iyy.ne.iy .or. jxx.ne.jx .or. kzz.ne.kz ) then
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
          call check_ok(istatus, 'Dimension time missing', &
                     &  'ICBC FILE ERROR')
          istatus = nf90_inquire_dimension(ibcin, idimid, len=ibcnrec)
          call check_ok(istatus, 'Dimension time read error', &
                     &  'ICBC FILE ERROR')
          istatus = nf90_inq_varid(ibcin, 'time', itvar)
          call check_ok(istatus, 'variable time missing', &
                     &  'ICBC FILE ERROR')
          istatus = nf90_get_att(ibcin, itvar, 'units', icbc_timeunits)
          call check_ok(istatus, 'variable time units missing', &
                     &  'ICBC FILE ERROR')
          allocate(icbc_idate(ibcnrec))
          allocate(icbc_xtime(ibcnrec))
          istatus = nf90_get_var(ibcin, itvar, icbc_xtime)
          call check_ok(istatus, 'variable time read error', &
                     &  'ICBC FILE ERROR')
          do i = 1 , ibcnrec
            icbc_idate(i) = timeval2idate(icbc_xtime(i), icbc_timeunits)
          end do
          chkdiff = icbc_xtime(2) - icbc_xtime(1)
          deallocate(icbc_xtime)
          if (chkdiff .ne. ibdyfrq) then
            write (6,*) 'Time variable in ICBC inconsistency.'
            write (6,*) 'Expecting ibdyfrq = ', ibdyfrq
            write (6,*) 'Found     ibdyfrq = ', chkdiff
            call fatal(__FILE__,__LINE__,'ICBC READ ERROR')
          end if
          istatus = nf90_inq_varid(ibcin, 'ps', icbc_ivar(1))
          call check_ok(istatus, 'variable ps missing', &
                     &  'ICBC FILE ERROR')
          istatus = nf90_inq_varid(ibcin, 'ts', icbc_ivar(2))
          call check_ok(istatus, 'variable ts missing', &
                     &  'ICBC FILE ERROR')
          istatus = nf90_inq_varid(ibcin, 'u', icbc_ivar(3))
          call check_ok(istatus, 'variable u missing', &
                     &  'ICBC FILE ERROR')
          istatus = nf90_inq_varid(ibcin, 'v', icbc_ivar(4))
          call check_ok(istatus, 'variable v missing', &
                     &  'ICBC FILE ERROR')
          istatus = nf90_inq_varid(ibcin, 't', icbc_ivar(5))
          call check_ok(istatus, 'variable t missing', &
                     &  'ICBC FILE ERROR')
          istatus = nf90_inq_varid(ibcin, 'qv', icbc_ivar(6))
          call check_ok(istatus, 'variable qv missing', &
                     &  'ICBC FILE ERROR')
          istatus = nf90_inq_varid(ibcin, 'so4', icbc_ivar(7))
          if ( istatus == nf90_noerr) then
            lso4p = .true.
          end if
        end subroutine open_icbc

        subroutine read_icbc(idate,ps,ts,u,v,t,qv,so4)
          use mod_dynparam
          use mod_message
          use netcdf
          use mod_date , only : idatediff
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
          call check_ok(istatus, 'variable ps read error', &
                     &  'ICBC FILE ERROR')
          ps = transpose(xread(:,:,1))
          istatus = nf90_get_var(ibcin, icbc_ivar(2), xread(:,:,1), & 
                                 istart(1:3), icount(1:3))
          call check_ok(istatus, 'variable ts read error', &
                     &  'ICBC FILE ERROR')
          ts = transpose(xread(:,:,1))
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
          call check_ok(istatus, 'variable u read error', &
                     &  'ICBC FILE ERROR')
          do k = 1 , kz
            do j = 1 , jx
              do i = 1 , iy
                u(i,k,j) = xread(j,i,k)
              end do
            end do
          end do
          istatus = nf90_get_var(ibcin, icbc_ivar(4), xread, &
                    &            istart, icount)
          call check_ok(istatus, 'variable v read error', &
                     &  'ICBC FILE ERROR')
          do k = 1 , kz
            do j = 1 , jx
              do i = 1 , iy
                v(i,k,j) = xread(j,i,k)
              end do
            end do
          end do
          istatus = nf90_get_var(ibcin, icbc_ivar(5), xread, &
                    &            istart, icount)
          call check_ok(istatus, 'variable t read error', &
                     &  'ICBC FILE ERROR')
          do k = 1 , kz
            do j = 1 , jx
              do i = 1 , iy
                t(i,k,j) = xread(j,i,k)
              end do
            end do
          end do
          istatus = nf90_get_var(ibcin, icbc_ivar(6), xread, &
                    &            istart, icount)
          call check_ok(istatus, 'variable qv read error', &
                     &  'ICBC FILE ERROR')
          do k = 1 , kz
            do j = 1 , jx
              do i = 1 , iy
                qv(i,k,j) = xread(j,i,k)
              end do
            end do
          end do
          if (lso4p) then
            istatus = nf90_get_var(ibcin, icbc_ivar(7), xread, &
                    &              istart, icount)
            call check_ok(istatus, 'variable so4 read error', &
                     &  'ICBC FILE ERROR')
            do k = 1 , kz
              do j = 1 , jx
                do i = 1 , iy
                  so4(i,k,j) = xread(j,i,k)
                end do
              end do
            end do
          else
            so4 = 0.0D+00
          end if
        end subroutine read_icbc

        subroutine close_icbc
          use mod_dynparam
          use netcdf
          implicit none
          if (ibcin >= 0) then
            istatus = nf90_close(ibcin)
            call check_ok(istatus, &
                     &  'Error Closing ICBC file '//trim(icbcname), &
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
            call check_ok(istatus, &
                     &  'Error Closing '//ctype//' file', &
                     &  ctype//' FILE ERROR')
            ncid = -1
          end if
        end subroutine close_common

        function icbc_search(idate)
          use mod_dynparam
          use mod_date , only : idatediff
          implicit none
          integer :: icbc_search
          integer , intent(in) :: idate
          if (idate > icbc_idate(ibcnrec) .or. idate < icbc_idate(1)) then
            icbc_search = -1
          else
            icbc_search = (idatediff(idate, icbc_idate(1))/ibdyfrq)+1
          end if 
        end function icbc_search

        subroutine prepare_common_io(idate,ctype)
          use mod_dynparam
          use mod_date , only : split_idate
          use mod_message
          use netcdf
          implicit none
          integer , intent(in) :: idate
          character(3) , intent(in) :: ctype
          character(64) :: title
          character(32) :: fbname , csdate , cmethodmax , cmethodmin
          character(256) :: ofname , history
          integer , dimension(8) :: tvals
          real(4) , dimension(2) :: trlat
          real(4) :: hptop
          real(4) , dimension(iy) :: yiy
          real(4) , dimension(jx) :: xjx
          integer :: ncid
          integer , dimension(2) :: izvar
          integer , dimension(2) :: ivvar
          integer , dimension(4) :: isrvvar
          integer , dimension(4) :: illtpvar
          integer :: itvar , iyy , im , id , ih , i , j

          if (ctype == 'ATM') then
            ncid = ncatm
            title = 'ICTP Regional Climatic model V4 ATM output'
          else if (ctype == 'SRF') then
            ncid = ncsrf
            title = 'ICTP Regional Climatic model V4 SRF output'
          else if (ctype == 'SUB') then
            ncid = ncsub
            title = 'ICTP Regional Climatic model V4 SUB output'
          else if (ctype == 'RAD') then
            ncid = ncrad
            title = 'ICTP Regional Climatic model V4 RAD output'
          else if (ctype == 'CHE') then
            ncid = ncche
            title = 'ICTP Regional Climatic model V4 CHE output'
          else
            write (aline,*) 'UNKNOWN IO TYPE : ', ctype
            call say
            write (aline,*) 'NOTHING TO DO'
            call say
            return
          end if

          call close_common(ncid, ctype)

          write (fbname,'(a,a,i10)') trim(ctype), '.', idate
          ofname = trim(dirout)//pthsep//trim(domname)// &
                &  '_'//trim(fbname)//'.nc'
          call split_idate(idate,iyy,im,id,ih)
          write (csdate,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a)') &
               & iyy,'-',im,'-',id,' ',ih,':00:00 UTC'

#ifdef NETCDF4_HDF5
          istatus = nf90_create(ofname, ior(nf90_clobber,nf90_hdf5), &
                              ncid)
#else
          istatus = nf90_create(ofname, nf90_clobber, ncid)
#endif
          call check_ok(istatus,('Error creating file '//trim(ofname)), &
                      ctype//' FILE ERROR')

          istatus = nf90_put_att(ncid, nf90_global, 'title', title)
          call check_ok(istatus, 'Error adding global title', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, nf90_global, 'institution', &
                   & 'ICTP')
          call check_ok(istatus, 'Error adding global institution', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, nf90_global, 'source', &
                   & 'RegCM Model '//SVN_REV//' simulation output')
          call check_ok(istatus, 'Error adding global source', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, nf90_global, 'Conventions', &
                   & 'CF-1.4')
          call check_ok(istatus, 'Error adding global Conventions', &
                        ctype//' FILE ERROR')
          call date_and_time(values=tvals)
          write (history,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a)') &
               tvals(1) , '-' , tvals(2) , '-' , tvals(3) , ' ' ,       &
               tvals(5) , ':' , tvals(6) , ':' , tvals(7) ,             &
               ' : Created by RegCM model'
          istatus = nf90_put_att(ncid, nf90_global, 'history', history)
          call check_ok(istatus, 'Error adding global history', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, nf90_global, 'references', &
                   & 'http://eforge.escience-lab.org/gf/project/regcm')
          call check_ok(istatus, 'Error adding global references', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, nf90_global, 'experiment', &
                   & domname)
          call check_ok(istatus, 'Error adding global experiment', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, nf90_global, 'projection', iproj)
          call check_ok(istatus, 'Error adding global projection', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, nf90_global,   &
                   &   'grid_size_in_meters', ds*1000.0)
          call check_ok(istatus, 'Error adding global gridsize', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, nf90_global,   &
                   &   'latitude_of_projection_origin', clat)
          call check_ok(istatus, 'Error adding global clat', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, nf90_global,   &
                   &   'longitude_of_projection_origin', clon)
          call check_ok(istatus, 'Error adding global clon', &
                        ctype//' FILE ERROR')
          if (iproj == 'ROTMER') then
            istatus = nf90_put_att(ncid, nf90_global, &
                     &   'latitude_of_projection_pole', plat)
            call check_ok(istatus, 'Error adding global plat', &
                          ctype//' FILE ERROR')
            istatus = nf90_put_att(ncid, nf90_global, &
                     &   'longitude_of_projection_pole', plon)
            call check_ok(istatus, 'Error adding global plon', &
                          ctype//' FILE ERROR')
          else if (iproj == 'LAMCON') then
            trlat(1) = truelatl
            trlat(2) = truelath
            istatus = nf90_put_att(ncid, nf90_global, &
                     &   'standard_parallel', trlat)
            call check_ok(istatus, 'Error adding global truelat', &
                          ctype//' FILE ERROR')
          end if
!
!         ADD RUN PARAMETERS
!
        ! TBD
!
!         ADD DIMENSIONS
!
          if (ctype == 'SUB') then
            istatus = nf90_def_dim(ncid, 'iy', iysg, idims(2))
            call check_ok(istatus, 'Error creating dimension iy', &
                          ctype//' FILE ERROR')
            istatus = nf90_def_dim(ncid, 'jx', jxsg, idims(1))
            call check_ok(istatus, 'Error creating dimension jx', &
                          ctype//' FILE ERROR')
          else
            istatus = nf90_def_dim(ncid, 'iy', iy, idims(2))
            call check_ok(istatus, 'Error creating dimension iy', &
                          ctype//' FILE ERROR')
            istatus = nf90_def_dim(ncid, 'jx', jx, idims(1))
            call check_ok(istatus, 'Error creating dimension jx', &
                          ctype//' FILE ERROR')
          end if
          istatus = nf90_def_dim(ncid, 'time', nf90_unlimited, &
                              &  idims(3))
          call check_ok(istatus, 'Error creating dimension time', &
                        ctype//' FILE ERROR')
          istatus = nf90_def_dim(ncid, 'kz', kz, idims(4))
          call check_ok(istatus, 'Error creating dimension kz', &
                        ctype//' FILE ERROR')
!
!         OUT TYPE DEPENDENT DIMENSIONS
!
          if (ctype == 'SRF' .or. ctype == 'SUB') then
            istatus = nf90_def_dim(ncid, 'm10', 1, idims(5))
            call check_ok(istatus, 'Error creating dimension m10', &
                          ctype//' FILE ERROR')
            istatus = nf90_def_dim(ncid, 'm2', 1, idims(6))
            call check_ok(istatus, 'Error creating dimension m2', &
                          ctype//' FILE ERROR')
            istatus = nf90_def_dim(ncid, 'soil_layer', 2, idims(7))
            call check_ok(istatus, &
                      &   'Error creating dimension soil_layer', &
                      &   ctype//' FILE ERROR')
            istatus = nf90_def_var(ncid, 'm10', nf90_float, idims(5), &
                                   isrvvar(1))
            call check_ok(istatus, 'Error adding variable m10', &
                          ctype//' FILE ERROR')
            istatus = nf90_put_att(ncid, isrvvar(1), 'standard_name', &
                            &  'altitude')
            call check_ok(istatus, 'Error adding m10 standard_name', &
                          ctype//' FILE ERROR')
            istatus = nf90_put_att(ncid, isrvvar(1), 'long_name', &
                            &  'Convenience 10 m elevation level')
            call check_ok(istatus, 'Error adding m10 long_name', &
                          ctype//' FILE ERROR')
            istatus = nf90_put_att(ncid, isrvvar(1), 'positive', 'up')
            call check_ok(istatus, 'Error adding m10 positive', &
                          ctype//' FILE ERROR')
            istatus = nf90_put_att(ncid, isrvvar(1), 'units', 'm')
            call check_ok(istatus, 'Error adding m10 units', &
                          ctype//' FILE ERROR')
            istatus = nf90_def_var(ncid, 'm2', nf90_float, idims(6), &
                                   isrvvar(2))
            call check_ok(istatus, 'Error adding variable m2', &
                          ctype//' FILE ERROR')
            istatus = nf90_put_att(ncid, isrvvar(2), 'standard_name', &
                            &  'altitude')
            call check_ok(istatus, 'Error adding m2 standard_name', &
                          ctype//' FILE ERROR')
            istatus = nf90_put_att(ncid, isrvvar(2), 'long_name', &
                            &  'Convenience 2 m elevation level')
            call check_ok(istatus, 'Error adding m2 long_name', &
                          ctype//' FILE ERROR')
            istatus = nf90_put_att(ncid, isrvvar(2), 'positive', 'up')
            call check_ok(istatus, 'Error adding m2 positive', &
                          ctype//' FILE ERROR')
            istatus = nf90_put_att(ncid, isrvvar(2), 'units', 'm')
            call check_ok(istatus, 'Error adding m2 units', &
                          ctype//' FILE ERROR')
            istatus = nf90_def_var(ncid, 'layer', nf90_float, idims(7), &
                                   isrvvar(3))
            call check_ok(istatus, 'Error adding variable layer', &
                          ctype//' FILE ERROR')
            istatus = nf90_put_att(ncid, isrvvar(3), 'standard_name', &
                            &  'model_level_number')
            call check_ok(istatus, 'Error adding layer standard_name', &
                          ctype//' FILE ERROR')
            istatus = nf90_put_att(ncid, isrvvar(3), 'long_name', &
                            &  'Surface and root zone')
            call check_ok(istatus, 'Error adding layer long_name', &
                          ctype//' FILE ERROR')
            istatus = nf90_put_att(ncid, isrvvar(3), 'positive', 'down')
            call check_ok(istatus, 'Error adding layer positive', &
                          ctype//' FILE ERROR')
            istatus = nf90_put_att(ncid, isrvvar(3), 'units', '1')
            call check_ok(istatus, 'Error adding layer units', &
                          ctype//' FILE ERROR')
          end if
          if (ctype == 'SRF') then
            istatus = nf90_def_dim(ncid, 'nv', 2, idims(8))
            call check_ok(istatus, 'Error creating dimension nv', &
                          ctype//' FILE ERROR')
          end  if
          if (ctype == 'CHE') then
            istatus = nf90_def_dim(ncid, 'tracer', ntr, idims(9))
            call check_ok(istatus, 'Error creating dimension tracer', &
                          ctype//' FILE ERROR')
          end if
          istatus = nf90_def_var(ncid, 'sigma', nf90_float, &
                              &  idims(4), izvar(1))
          call check_ok(istatus, 'Error adding variable sigma', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, izvar(1), 'standard_name', &
                            &  'atmosphere_sigma_coordinate')
          call check_ok(istatus, 'Error adding sigma standard_name', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, izvar(1), 'long_name', &
                            &  'Sigma at model layers')
          call check_ok(istatus, 'Error adding sigma long_name', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, izvar(1), 'units', '1')
          call check_ok(istatus, 'Error adding sigma units', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, izvar(1), 'axis', 'Z')
          call check_ok(istatus, 'Error adding sigma axis', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, izvar(1), 'positive', 'down')
          call check_ok(istatus, 'Error adding sigma positive', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, izvar(1), 'formula_terms',  &
                     &         'sigma: sigma ps: ps ptop: ptop')
          call check_ok(istatus, 'Error adding sigma formula_terms', &
                        ctype//' FILE ERROR')
          istatus = nf90_def_var(ncid, 'ptop', nf90_float, &
                           &   varid=izvar(2))
          call check_ok(istatus, 'Error adding variable ptop', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, izvar(2), 'standard_name',  &
                            &  'air_pressure')
          call check_ok(istatus, 'Error adding ptop standard_name', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, izvar(2), 'long_name', &
                            &  'Pressure at model top')
          call check_ok(istatus, 'Error adding ptop long_name', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, izvar(2), 'units', 'hPa')
          call check_ok(istatus, 'Error adding ptop units', &
                        ctype//' FILE ERROR')
          istatus = nf90_def_var(ncid, 'iy', nf90_float, idims(2), &
                            &  ivvar(1))
          call check_ok(istatus, 'Error adding variable iy', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, ivvar(1), 'standard_name',  &
                            &  'projection_y_coordinate')
          call check_ok(istatus, 'Error adding iy standard_name', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, ivvar(1), 'long_name', &
                            &  'y-coordinate in Cartesian system')
          call check_ok(istatus, 'Error adding iy long_name', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, ivvar(1), 'units', 'km')
          call check_ok(istatus, 'Error adding iy units', &
                        ctype//' FILE ERROR')
          istatus = nf90_def_var(ncid, 'jx', nf90_float, idims(1), &
                            &  ivvar(2))
          call check_ok(istatus, 'Error adding variable jx', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, ivvar(2), 'standard_name', &
                            &  'projection_x_coordinate')
          call check_ok(istatus, 'Error adding jx standard_name', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, ivvar(2), 'long_name', &
                            &  'x-coordinate in Cartesian system')
          call check_ok(istatus, 'Error adding jx long_name', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, ivvar(2), 'units', 'km')
          call check_ok(istatus, 'Error adding jx units', &
                        ctype//' FILE ERROR')
          istatus = nf90_def_var(ncid, 'xlat', nf90_float, &
                             &   idims(1:2), illtpvar(1))
          call check_ok(istatus, 'Error adding variable xlat', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, illtpvar(1), 'standard_name', &
                            &  'latitude')
          call check_ok(istatus, 'Error adding xlat standard_name', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, illtpvar(1), 'long_name', &
                            &  'Latitude at cross points')
          call check_ok(istatus, 'Error adding xlat long_name', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, illtpvar(1), 'units', &
                            &  'degrees_north')
          call check_ok(istatus, 'Error adding xlat units', &
                        ctype//' FILE ERROR')
          istatus = nf90_def_var(ncid, 'xlon', nf90_float, &
                             &   idims(1:2), illtpvar(2))
          call check_ok(istatus, 'Error adding variable xlon', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, illtpvar(2), 'standard_name', &
                            &  'longitude')
          call check_ok(istatus, 'Error adding xlon standard_name', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, illtpvar(2), 'long_name', &
                            &  'Longitude at cross points')
          call check_ok(istatus, 'Error adding xlon long_name', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, illtpvar(2), 'units',  &
                            &  'degrees_east')
          call check_ok(istatus, 'Error adding xlon units', &
                        ctype//' FILE ERROR')
          istatus = nf90_def_var(ncid, 'topo', nf90_float, &
                             &   idims(1:2), illtpvar(3))
          call check_ok(istatus, 'Error adding variable topo', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, illtpvar(3), 'standard_name', &
                            &  'surface_altitude')
          call check_ok(istatus, 'Error adding topo standard_name', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, illtpvar(3), 'long_name',     &
                            &  'Domain surface elevation')
          call check_ok(istatus, 'Error adding topo long_name', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, illtpvar(3), 'units', 'm')
          call check_ok(istatus, 'Error adding topo units', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, illtpvar(3), 'coordinates', &
                            &  'xlat xlon')
          call check_ok(istatus, 'Error adding topo coordinates', &
                        ctype//' FILE ERROR')
          istatus = nf90_def_var(ncid, 'time', nf90_double, &
                               & idims(3:3), itvar)
          call check_ok(istatus, 'Error adding variable time', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, itvar, 'units', &
                             &   'hours since '//csdate)
          call check_ok(istatus, 'Error adding time units', &
                        ctype//' FILE ERROR')
          if (ctype == 'SRF') then
            istatus = nf90_put_att(ncid, itvar, 'bounds', 'time_bnds')
            call check_ok(istatus, 'Error adding time bounds', &
                          ctype//' FILE ERROR')
            istatus = nf90_def_var(ncid, 'time_bnds', nf90_double, &
                                   (/idims(3),idims(8)/), isrfvar(2))
            call check_ok(istatus, 'Error adding variable time_bnds', &
                          ctype//' FILE ERROR')
            istatus = nf90_put_att(ncid, isrfvar(2), 'calendar', &
                                   'standard')
            call check_ok(istatus, 'Error adding time_bnds calendar', &
                          ctype//' FILE ERROR')
            istatus = nf90_put_att(ncid, isrfvar(2), 'units', &
                             &   'hours since '//csdate)
            call check_ok(istatus, 'Error adding time_bnds units', &
                          ctype//' FILE ERROR')
          end if
          istatus = nf90_def_var(ncid, 'ps', nf90_float, &
                             &   idims(1:3), illtpvar(4))
          call check_ok(istatus, 'Error adding variable ps', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, illtpvar(4), 'standard_name', &
                            &  'surface_air_pressure')
          call check_ok(istatus, 'Error adding ps standard_name', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, illtpvar(4), 'long_name',     &
                            &  'Surface pressure')
          call check_ok(istatus, 'Error adding ps long_name', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, illtpvar(4), 'units', 'hPa')
          call check_ok(istatus, 'Error adding ps units', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, illtpvar(4), 'coordinates', &
                            &  'xlat xlon')
          call check_ok(istatus, 'Error adding ps coordinates', &
                        ctype//' FILE ERROR')

          if (ctype == 'ATM') then
            iatmvar = -1
            iatmvar(1) = itvar
            iatmvar(2) = illtpvar(4)
            call addvara(ncid,ctype,'u','eastward_wind', &
                'U component (westerly) of wind','m s-1', &
                (/idims(1),idims(2),idims(4),idims(3),-1,-1,-1,-1,-1/), &
                .false.,iatmvar(3))
            call addvara(ncid,ctype,'v','northward_wind', &
                'V component (southerly) of wind','m s-1', &
                (/idims(1),idims(2),idims(4),idims(3),-1,-1,-1,-1,-1/), &
                .false.,iatmvar(4))
            call addvara(ncid,ctype,'omega', &
                'lagrangian_tendency_of_air_pressure', &
                'Pressure velocity','hPa s-1', &
                (/idims(1),idims(2),idims(4),idims(3),-1,-1,-1,-1,-1/), &
                .false.,iatmvar(5))
            call addvara(ncid,ctype,'t','air_temperature', &
                'Temperature','K', &
                (/idims(1),idims(2),idims(4),idims(3),-1,-1,-1,-1,-1/), &
                .false.,iatmvar(6))
            call addvara(ncid,ctype,'qv','humidity_mixing_ratio', &
                'Water vapor mixing ratio','kg kg-1', &
                (/idims(1),idims(2),idims(4),idims(3),-1,-1,-1,-1,-1/), &
                .false.,iatmvar(7))
            call addvara(ncid,ctype,'qc', &
                'cloud_liquid_water_mixing_ratio', &
                'Cloud water mixing ratio','kg kg-1', &
                (/idims(1),idims(2),idims(4),idims(3),-1,-1,-1,-1,-1/), &
                .false.,iatmvar(8))
            call addvara(ncid,ctype,'tpr','precipitation_flux', &
                'Total daily precipitation rate','kg m-2 day-1', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,iatmvar(9))
            call addvara(ncid,ctype,'tgb','soil_temperature', &
                'Lower groud temperature in BATS','K', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,iatmvar(10))
            call addvara(ncid,ctype,'swt', &
                'moisture_content_of_soil_layer', &
                'Total soil water','kg m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .true.,iatmvar(11))
            call addvara(ncid,ctype,'rno','runoff_amount', &
                'Runoff accumulated infiltration','kg m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .true.,iatmvar(12))
          else if (ctype == 'SRF') then
            isrfvar = -1
            isrfvar(1) = itvar
            isrfvar(3) = illtpvar(4)
            call addvara(ncid,ctype,'u10m','eastward_wind', &
                '10 meters U component (westerly) of wind','m s-1', &
                (/idims(1),idims(2),idims(5),idims(3),-1,-1,-1,-1,-1/), &
                .false.,isrfvar(4))
            call addvara(ncid,ctype,'v10m','northward_wind', &
                '10 meters V component (southerly) of wind','m s-1', &
                (/idims(1),idims(2),idims(5),idims(3),-1,-1,-1,-1,-1/), &
                .false.,isrfvar(5))
            call addvara(ncid,ctype,'uvdrag', &
                'surface_drag_coefficient_in_air', &
                'Surface drag stress','1', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,isrfvar(6))
            call addvara(ncid,ctype,'tg','surface_temperature', &
                'Ground temperature','K', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,isrfvar(7))
            call addvara(ncid,ctype,'tlef','canopy_temperature', &
                'Foliage temperature','K', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .true.,isrfvar(8))
            call addvara(ncid,ctype,'t2m','air_temperature', &
                '2 meters temperature','K', &
                (/idims(1),idims(2),idims(6),idims(3),-1,-1,-1,-1,-1/), &
                .false.,isrfvar(9))
            call addvara(ncid,ctype,'q2m','humidity_mixing_ratio', &
                '2 meters vapour mixing ratio','kg kg-1', &
                (/idims(1),idims(2),idims(6),idims(3),-1,-1,-1,-1,-1/), &
                .false.,isrfvar(10))
            call addvara(ncid,ctype,'smw','soil_moisture_content', &
                'Moisture content','kg kg-1', &
                (/idims(1),idims(2),idims(7),idims(3),-1,-1,-1,-1,-1/), &
                .true.,isrfvar(11))
            call addvara(ncid,ctype,'tpr','precipitation_amount', &
                'Total precipitation','kg m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,isrfvar(12))
            call addvara(ncid,ctype,'evp','water_evaporation_amount', &
                'Total evapotranspiration','kg m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,isrfvar(13))
            call addvara(ncid,ctype,'runoff','surface_runoff_flux', &
                'Surface runoff','kg m-2 day-1', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .true.,isrfvar(14))
            call addvara(ncid,ctype,'scv','snowfall_amount', &
                'Snow precipitation','kg m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .true.,isrfvar(15))
            call addvara(ncid,ctype,'sena', &
                'surface_downward_sensible_heat_flux', &
                'Sensible heat flux','W m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,isrfvar(16))
            call addvara(ncid,ctype,'fsw', &
                'surface_downwelling_shortwave_flux_in_air', &
                'Solar absorbed energy flux','W m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,isrfvar(17))
            call addvara(ncid,ctype,'fld', &
                'surface_downwelling_longwave_flux_in_air', &
                'Downward LW flux','W m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,isrfvar(18))
            call addvara(ncid,ctype,'sina', &
              'net_downward_radiative_flux_at_top_of_atmosphere_model', &
              'Incident solar energy flux','W m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,isrfvar(19))
            call addvara(ncid,ctype,'prcv', 'convective_rainfall_flux', &
                'Convective precipitation','kg m-2 day-1', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,isrfvar(20))
            call addvara(ncid,ctype,'zpbl', &
                'atmosphere_boundary_layer_thickness', &
                'PBL layer thickness','m', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,isrfvar(21))
            write (cmethodmax, '(a,i3,a)') 'time: maximum (interval: ', &
                   int(batfrq) , ' hours)'
            write (cmethodmin, '(a,i3,a)') 'time: minimum (interval: ', &
                   int(batfrq) , ' hours)'
            call addvara(ncid,ctype,'tgmax','surface_temperature', &
                'Maximum surface temperature','K', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,isrfvar(22))
            istatus = nf90_put_att(ncid, isrfvar(22), 'cell_methods', &
                            &  cmethodmax)
            call check_ok(istatus, 'Error adding tgmax cell_methods', &
                        ctype//' FILE ERROR')
            call addvara(ncid,ctype,'tgmin','surface_temperature', &
                'Minimum surface temperature','K', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,isrfvar(23))
            istatus = nf90_put_att(ncid, isrfvar(23), 'cell_methods', &
                            &  cmethodmin)
            call check_ok(istatus, 'Error adding tgmin cell_methods', &
                        ctype//' FILE ERROR')
            call addvara(ncid,ctype,'t2max','air_temperature', &
                'Maximum 2 meters temperature','K', &
                (/idims(1),idims(2),idims(6),idims(3),-1,-1,-1,-1,-1/), &
                .false.,isrfvar(24))
            istatus = nf90_put_att(ncid, isrfvar(24), 'cell_methods', &
                            &  cmethodmax)
            call check_ok(istatus, 'Error adding t2max cell_methods', &
                        ctype//' FILE ERROR')
            call addvara(ncid,ctype,'t2min','air_temperature', &
                'Minimum 2 meters temperature','K', &
                (/idims(1),idims(2),idims(6),idims(3),-1,-1,-1,-1,-1/), &
                .false.,isrfvar(25))
            istatus = nf90_put_att(ncid, isrfvar(25), 'cell_methods', &
                            &  cmethodmin)
            call check_ok(istatus, 'Error adding t2min cell_methods', &
                        ctype//' FILE ERROR')
            call addvara(ncid,ctype,'w10max','wind_speed', &
                'Maximum speed of 10m wind','m s-1', &
                (/idims(1),idims(2),idims(5),idims(3),-1,-1,-1,-1,-1/), &
                .false.,isrfvar(26))
            istatus = nf90_put_att(ncid, isrfvar(26), 'cell_methods', &
                            &  cmethodmax)
            call check_ok(istatus, 'Error adding w10max cell_methods', &
                        ctype//' FILE ERROR')
            call addvara(ncid,ctype,'ps_min','air_pressure', &
                'Minimum of surface pressure','hPa', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,isrfvar(27))
            istatus = nf90_put_att(ncid, isrfvar(27), 'cell_methods', &
                            &  cmethodmin)
            call check_ok(istatus, 'Error adding ps_min cell_methods', &
                        ctype//' FILE ERROR')
          else if (ctype == 'SUB') then
            isubvar = -1
            isubvar(1) = itvar
            isubvar(2) = illtpvar(4)
            call addvara(ncid,ctype,'u10m','eastward_wind', &
                '10 meters U component (westerly) of wind','m s-1', &
                (/idims(1),idims(2),idims(5),idims(3),-1,-1,-1,-1,-1/), &
                .false.,isubvar(3))
            call addvara(ncid,ctype,'v10m','northward_wind', &
                '10 meters V component (southerly) of wind','m s-1', &
                (/idims(1),idims(2),idims(5),idims(3),-1,-1,-1,-1,-1/), &
                .false.,isubvar(4))
            call addvara(ncid,ctype,'uvdrag', &
                'surface_drag_coefficient_in_air', &
                'Surface drag stress','1', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,isubvar(5))
            call addvara(ncid,ctype,'tg','surface_temperature', &
                'Ground temperature','K', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,isubvar(6))
            call addvara(ncid,ctype,'tlef','canopy_temperature', &
                'Foliage temperature','K', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .true.,isubvar(7))
            call addvara(ncid,ctype,'t2m','air_temperature', &
                '2 meters temperature','K', &
                (/idims(1),idims(2),idims(6),idims(3),-1,-1,-1,-1,-1/), &
                .false.,isubvar(8))
            call addvara(ncid,ctype,'q2m','humidity_mixing_ratio', &
                '2 meters vapour mixing ratio','kg kg-1', &
                (/idims(1),idims(2),idims(6),idims(3),-1,-1,-1,-1,-1/), &
                .false.,isubvar(9))
            call addvara(ncid,ctype,'smw','soil_moisture_content', &
                'Moisture content','kg kg-1', &
                (/idims(1),idims(2),idims(7),idims(3),-1,-1,-1,-1,-1/), &
                .true.,isubvar(10))
            call addvara(ncid,ctype,'tpr','precipitation_amount', &
                'Total precipitation','kg m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,isubvar(11))
            call addvara(ncid,ctype,'evp','water_evaporation_amount', &
                'Total evapotranspiration','kg m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,isubvar(12))
            call addvara(ncid,ctype,'runoff','surface_runoff_flux', &
                'Surface runoff','kg m-2 day-1', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .true.,isubvar(13))
            call addvara(ncid,ctype,'scv','snowfall_amount', &
                'Snow precipitation','kg m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .true.,isubvar(14))
            call addvara(ncid,ctype,'sena', &
                'surface_downward_sensible_heat_flux', &
                'Sensible heat flux','W m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,isubvar(15))
            call addvara(ncid,ctype,'prcv', 'convective_rainfall_flux', &
                'Convective precipitation','kg m-2 day-1', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,isubvar(16))
          else if (ctype == 'RAD') then
            iradvar = -1
            iradvar(1) = itvar
            iradvar(2) = illtpvar(4)
            call addvara(ncid,ctype,'cld', &
                'cloud_area_fraction_in_atmosphere_layer', &
                'Cloud fractional cover','1', &
                (/idims(1),idims(2),idims(4),idims(3),-1,-1,-1,-1,-1/), &
                .false.,iradvar(3))
            call addvara(ncid,ctype,'clwp', &
                'atmosphere_optical_thickness_due_to_cloud', &
                'Cloud liquid water path','1', &
                (/idims(1),idims(2),idims(4),idims(3),-1,-1,-1,-1,-1/), &
                .false.,iradvar(4))
            call addvara(ncid,ctype,'qrs', &
                'tendency_of_air_temperature_due_to_shortwave_heating', &
                'Solar heating rate','K s-1', &
                (/idims(1),idims(2),idims(4),idims(3),-1,-1,-1,-1,-1/), &
                .false.,iradvar(5))
            call addvara(ncid,ctype,'qrl', &
                'tendency_of_air_temperature_due_to_longwave_heating', &
                'Longwave cooling rate','K s-1', &
                (/idims(1),idims(2),idims(4),idims(3),-1,-1,-1,-1,-1/), &
                .false.,iradvar(6))
            call addvara(ncid,ctype,'frsa', &
                'surface_downwelling_shortwave_flux_in_air', &
                'Surface absorbed solar flux','W m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,iradvar(7))
            call addvara(ncid,ctype,'frla', &
                'downwelling_longwave_flux_in_air', &
                'Longwave cooling of surface flux','W m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,iradvar(8))
            call addvara(ncid,ctype,'clrst', &
                'downwelling_shortwave_flux_in_air_assuming_clear_sky', &
                'clearsky total column absorbed solar flux','W m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,iradvar(9))
            call addvara(ncid,ctype,'clrss', &
                'net_downward_shortwave_flux_in_air_assuming_clear_sky', &
                'clearsky surface absorbed solar flux','W m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,iradvar(10))
            call addvara(ncid,ctype,'clrlt', &
                'toa_net_upward_longwave_flux_assuming_clear_sky', &
                'clearsky net upward LW flux at TOA','W m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,iradvar(11))
            call addvara(ncid,ctype,'clrls', &
                'net_upward_longwave_flux_in_air_assuming_clear_sky', &
                'clearsky LW cooling at surface','W m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,iradvar(12))
            call addvara(ncid,ctype,'solin', &
                'toa_instantaneous_shortwave_forcing', &
                'Instantaneous incident solar','W m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,iradvar(13))
            call addvara(ncid,ctype,'sabtp', &
              'atmosphere_net_rate_of_absorption_of_shortwave_energy', &
              'Total column absorbed solar flux','W m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,iradvar(14))
            call addvara(ncid,ctype,'firtp', &
                'atmosphere_net_rate_of_absorption_of_longwave_energy', &
                'net upward LW flux at TOA','W m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,iradvar(15))
          else if (ctype == 'CHE') then
            ichevar = -1
            ichevar(1) = itvar
            ichevar(2) = illtpvar(4)
            call addvara(ncid,ctype,'trac', &
                'atmosphere_mixing_ratio_of_tracer', &
                'Tracers mixing ratios','kg kg-1', &
          (/idims(1),idims(2),idims(4),idims(9),idims(3),-1,-1,-1,-1/), &
                .false.,ichevar(3))
            call addvara(ncid,ctype,'aext8', &
                'aerosol_extincion_coefficient', &
                'aer mix. ext. coef','1', &
                (/idims(1),idims(2),idims(4),idims(3),-1,-1,-1,-1,-1/), &
                .false.,ichevar(4))
            call addvara(ncid,ctype,'assa8', &
                'aerosol_single_scattering_albedo', &
                'aer mix. sin. scat. alb','1', &
                (/idims(1),idims(2),idims(4),idims(3),-1,-1,-1,-1,-1/), &
                .false.,ichevar(5))
            call addvara(ncid,ctype,'agfu8', &
                'aerosol_asymmetry_parameter', &
                'aer mix. sin. scat. alb','1', &
                (/idims(1),idims(2),idims(4),idims(3),-1,-1,-1,-1,-1/), &
                .false.,ichevar(6))
            call addvara(ncid,ctype,'colb', &
                'instantaneous_deposition_of_tracer', &
                'columnburden inst','mg m-2', &
                (/idims(1),idims(2),idims(9),idims(3),-1,-1,-1,-1,-1/), &
                .false.,ichevar(7))
            call addvara(ncid,ctype,'wdlsc', &
                'tendency_of_wet_deposition_of_tracer'// &
                '_due_to_large_scale_precipitation', &
                'wet dep lgscale','mg m-2 day-1', &
                (/idims(1),idims(2),idims(9),idims(3),-1,-1,-1,-1,-1/), &
                .false.,ichevar(8))
            call addvara(ncid,ctype,'wdcvc', &
                'tendency_of_wet_deposition_of_tracer'// &
                '_due_to_convective_precipitation', &
                'wet dep convect','mg m-2 day-1', &
                (/idims(1),idims(2),idims(9),idims(3),-1,-1,-1,-1,-1/), &
                .false.,ichevar(9))
            call addvara(ncid,ctype,'sdrdp', &
                'tendency_of_dry_deposition_of_tracer', &
                'surf dry depos','mg m-2 day-1', &
                (/idims(1),idims(2),idims(9),idims(3),-1,-1,-1,-1,-1/), &
                .false.,ichevar(10))
            call addvara(ncid,ctype,'xgasc', &
                'tendency_of_gas_conversion_of_tracer', &
                'chem gas conv','mg m-2 day-1', &
                (/idims(1),idims(2),idims(9),idims(3),-1,-1,-1,-1,-1/), &
                .false.,ichevar(11))
            call addvara(ncid,ctype,'xaquc', &
                'tendency_of_aqueous_conversion_of_tracer', &
                'chem aqu conv','mg m-2 day-1', &
                (/idims(1),idims(2),idims(9),idims(3),-1,-1,-1,-1,-1/), &
                .false.,ichevar(12))
            call addvara(ncid,ctype,'emiss', &
                'tendency_of_surface_emission_of_tracer', &
                'surf emission','mg m-2 day-1', &
                (/idims(1),idims(2),idims(9),idims(3),-1,-1,-1,-1,-1/), &
                .false.,ichevar(13))
            call addvara(ncid,ctype,'acstoarf', &
                'toa_instantaneous_shortwave_radiative_forcing', &
                'TOArad SW forcing av.','W m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,ichevar(14))
            call addvara(ncid,ctype,'acstsrrf', &
                'surface_instantaneous_shortwave_radiative_forcing', &
                'SRFrad SW forcing av.','W m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,ichevar(15))
            call addvara(ncid,ctype,'acstalrf', &
                'toa_instantaneous_longwave_radiative_forcing', &
                'TOArad LW forcing av.','W m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,ichevar(16))
            call addvara(ncid,ctype,'acssrlrf', &
                'surface_instantaneous_longwave_radiative_forcing', &
                'SRFrad LW forcing av.','W m-2', &
                (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/), &
                .false.,ichevar(17))
          end if

          istatus = nf90_enddef(ncid)
          call check_ok(istatus, 'Error End Definitions NetCDF output', &
                        ctype//' FILE ERROR')

          istatus = nf90_put_var(ncid, izvar(1), hsigma)
          call check_ok(istatus, 'Error variable sigma write', &
                        ctype//' FILE ERROR')
          hptop = ptop*10.0
          istatus = nf90_put_var(ncid, izvar(2), hptop)
          call check_ok(istatus, 'Error variable ptop write', &
                        ctype//' FILE ERROR')
          yiy(1) = -(dble(iy-1)/2.0) * ds
          xjx(1) = -(dble(jx-1)/2.0) * ds
          do i = 2 , iy
            yiy(i) = yiy(i-1)+ds
          end do
          do j = 2 , jx
            xjx(j) = xjx(j-1)+ds
          end do
          istatus = nf90_put_var(ncid, ivvar(1), yiy)
          call check_ok(istatus, 'Error variable iy write', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_var(ncid, ivvar(2), xjx)
          call check_ok(istatus, 'Error variable jx write', &
                        ctype//' FILE ERROR')
          if (ctype == 'SUB') then
            istatus = nf90_put_var(ncid, illtpvar(1), ioxlat_s)
            call check_ok(istatus, 'Error variable xlat write', &
                        ctype//' FILE ERROR')
            istatus = nf90_put_var(ncid, illtpvar(2), ioxlon_s)
            call check_ok(istatus, 'Error variable xlon write', &
                        ctype//' FILE ERROR')
            istatus = nf90_put_var(ncid, illtpvar(3), iotopo_s)
            call check_ok(istatus, 'Error variable topo write', &
                        ctype//' FILE ERROR')
          else
            istatus = nf90_put_var(ncid, illtpvar(1), ioxlat)
            call check_ok(istatus, 'Error variable xlat write', &
                        ctype//' FILE ERROR')
            istatus = nf90_put_var(ncid, illtpvar(2), ioxlon)
            call check_ok(istatus, 'Error variable xlon write', &
                        ctype//' FILE ERROR')
            istatus = nf90_put_var(ncid, illtpvar(3), iotopo)
            call check_ok(istatus, 'Error variable topo write', &
                        ctype//' FILE ERROR')
          end if
          if (ctype == 'SRF' .or. ctype == 'SUB') then
            istatus = nf90_put_var(ncid, isrvvar(1), 10)
            call check_ok(istatus, 'Error variable m10 write', &
                        ctype//' FILE ERROR')
            istatus = nf90_put_var(ncid, isrvvar(2), 2)
            call check_ok(istatus, 'Error variable m2 write', &
                        ctype//' FILE ERROR')
            istatus = nf90_put_var(ncid, isrvvar(3), (/0,1/))
            call check_ok(istatus, 'Error variable layer write', &
                        ctype//' FILE ERROR')
          end if

        end subroutine prepare_common_io

        subroutine addvara(ncid,ctype,vname,vst,vln,vuni,idims,lmiss, &
                           ivar)
          use netcdf
          implicit none
          integer , intent(in) :: ncid
          character(3) , intent(in) :: ctype
          character(len=*) , intent(in) :: vname
          character(len=*) , intent(in) :: vst , vln , vuni
          integer , dimension(9) , intent(in) :: idims
          logical , intent(in) :: lmiss
          integer , intent(out) :: ivar

          real(4) , parameter :: fillv = -1E+34
          character(64) :: cdum
          integer :: i , ndims

          ndims = 1
          do i = 1 , 9
            if (idims(1) > 0) ndims = ndims+1
          end do

          cdum = vname
          istatus = nf90_def_var(ncid, cdum, nf90_float, &
                             &   idims(1:ndims), ivar)
          call check_ok(istatus, 'Error adding variable '//vname, &
                        ctype//' FILE ERROR')
          cdum = vst
          istatus = nf90_put_att(ncid, ivar, 'standard_name', cdum)
          call check_ok(istatus, &
                        'Error adding '//vname//' standard_name', &
                        ctype//' FILE ERROR')
          cdum = vln
          istatus = nf90_put_att(ncid, ivar, 'long_name', cdum)
          call check_ok(istatus, &
                        'Error adding '//vname//' long_name', &
                        ctype//' FILE ERROR')
          cdum = vuni
          istatus = nf90_put_att(ncid, ivar, 'units', cdum)
          call check_ok(istatus, &
                        'Error adding '//vname//' units', &
                        ctype//' FILE ERROR')
          istatus = nf90_put_att(ncid, ivar, 'coordinates', &
                            &  'xlat xlon')
          call check_ok(istatus, &
                        'Error adding '//vname//' coordinates', &
                        ctype//' FILE ERROR')
          if (lmiss) then
            istatus = nf90_put_att(ncid, ivar, '_FillValue', &
                              &  fillv)
            call check_ok(istatus, &
                          'Error adding '//vname//' coordinates', &
                          ctype//' FILE ERROR')
          end if
        end subroutine addvara

        subroutine check_ok(ierr,m1,mf)
          use mod_message
          use netcdf
          implicit none
          integer , intent(in) :: ierr
          character(*) :: m1 , mf
          if (ierr /= nf90_noerr) then 
            write (6,*) m1
            write (6,*) nf90_strerror(ierr)
            call fatal(__FILE__,__LINE__,mf)
          end if
        end subroutine check_ok

        subroutine release_mod_ncio
          call close_domain
          call close_icbc
          call close_common(ncatm,'ATM')
          call close_common(ncsrf,'SRF')
          call close_common(ncsub,'SUB')
          call close_common(ncrad,'RAD')
          call close_common(ncche,'CHE')
          if (allocated(ioxlat)) deallocate(ioxlat)
          if (allocated(ioxlon)) deallocate(ioxlon)
          if (allocated(iotopo)) deallocate(iotopo)
          if (allocated(hsigma)) deallocate(hsigma)
          if (allocated(ioxlat_s)) deallocate(ioxlat_s)
          if (allocated(ioxlon_s)) deallocate(ioxlon_s)
          if (allocated(iotopo_s)) deallocate(iotopo_s)
        end subroutine release_mod_ncio

      end module mod_ncio
