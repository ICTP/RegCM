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

        use mod_dynparam

        integer , private :: iutin , iutin1
        character(256) , private :: dname , sdname , aername

        data iutin  /-1/
        data iutin1 /-1/

      contains

        subroutine static_path_prepare
          implicit none
          character(3) :: sbstring
          dname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
          write (sbstring,'(i0.3)') nsg
          sdname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN'// &
              &    sbstring//'.nc'
          aername = trim(dirglob)//pthsep//trim(domname)//'_AERO.nc'
        end subroutine static_path_prepare

        subroutine open_domain(r8pt , dx , sigma)

          use mod_message
          use netcdf
          implicit none

          real(8) , intent(out) :: dx
          real(8) , intent(out) :: r8pt
          real(8) , dimension(kzp1) :: sigma

          integer :: istatus
          integer :: ivarid , idimid
          integer :: iyy , jxx , kzz
          character(6) :: proj
          real(4) :: dsx , iclat , iclon , ptsp

          write (aline,*) 'READING HEADER FILE:', dname
          call say
          istatus = nf90_open(dname, nf90_nowrite, iutin)
          if ( istatus /= nf90_noerr) then
            write (6,*) 'Error Opening Domain file ', trim(dname)
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'CANNOT OPEN DOMAIN FILE')
          end if
          if ( nsg.gt.1 ) then
            write (aline,*) 'READING HEADER SUBDOMAIN FILE:', sdname
            call say
            istatus = nf90_open(sdname, nf90_nowrite, iutin1)
            if ( istatus /= nf90_noerr) then
              write (6,*) 'Error Opening SubDomain file ', trim(sdname)
              write (6,*) nf90_strerror(istatus)
              call fatal(__FILE__,__LINE__, 'CANNOT OPEN SUBDOM FILE')
            end if
          end if
          istatus = nf90_inq_dimid(iutin, "iy", idimid)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Dimension iy missing'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          istatus = nf90_inquire_dimension(iutin, idimid, len=iyy)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error dimension iy'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          istatus = nf90_inq_dimid(iutin, "jx", idimid)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Dimension jx missing'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          istatus = nf90_inquire_dimension(iutin, idimid, len=jxx)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error dimension jx'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          istatus = nf90_inq_dimid(iutin, "kz", idimid)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Dimension kz missing'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          istatus = nf90_inquire_dimension(iutin, idimid, len=kzz)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error dimension kz'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          istatus = nf90_inq_varid(iutin, "ptop", ivarid)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error ptop variable undefined'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          istatus = nf90_get_var(iutin, ivarid, ptsp)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error reading ptop variable'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          istatus = nf90_get_att(iutin, nf90_global, 'projection', proj)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error reading projection attribute'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          istatus = nf90_get_att(iutin, nf90_global, &
                       &         'grid_size_in_meters', dsx)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error reading grid_size_in_meters attribute'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          istatus = nf90_get_att(iutin, nf90_global, &
                     &         'latitude_of_projection_origin', iclat)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error reading latitude_of_projection_origin'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          istatus = nf90_get_att(iutin, nf90_global, &
                      &         'longitude_of_projection_origin', iclon)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error reading longitude_of_projection_origin'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
!
!         Consistency Check
!
          if ( iyy.ne.iy .or. jxx.ne.jx .or. kzz.ne.kzp1 ) then
            write (aline,*) 'param:  SET IN regcm.in   :  IY=' , iy , &
                   & '  JX=' ,  jx , '  KX=' , kz
            call say
            write (aline,*) 'param:  SET IN DOMAIN file: IYY=' , iyy ,&
                   & ' JXX=' , jxx , ' KZZ=' , kzz
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
          istatus = nf90_inq_varid(iutin, "sigma", ivarid)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error sigma variable undefined'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          istatus = nf90_get_var(iutin, ivarid, sigma)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error reading sigma variable'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if

        end subroutine open_domain

        subroutine read_domain(ht,htsd,lnd,xlat,xlon,xmap,dmap,f,snw)
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

          integer :: istatus , ivarid , n
          real(4) , dimension(jx,iy) :: sp2d

          if (iutin < 0) then
            write (6,*) 'Error : Domain file not in open state'
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if

          istatus = nf90_inq_varid(iutin, "topo", ivarid)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error topo variable undefined'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          istatus = nf90_get_var(iutin, ivarid, sp2d)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error reading topo variable'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          ht = transpose(sp2d)
          istatus = nf90_inq_varid(iutin, "htsd", ivarid)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error htsd variable undefined'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          istatus = nf90_get_var(iutin, ivarid, sp2d)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error reading htsd variable'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          htsd = transpose(sp2d)
          istatus = nf90_inq_varid(iutin, "landuse", ivarid)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error landuse variable undefined'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          istatus = nf90_get_var(iutin, ivarid, sp2d)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error reading landuse variable'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          lnd = transpose(sp2d)
          istatus = nf90_inq_varid(iutin, "xlat", ivarid)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error xlat variable undefined'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          istatus = nf90_get_var(iutin, ivarid, sp2d)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error reading xlat variable'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          xlat = transpose(sp2d)
          istatus = nf90_inq_varid(iutin, "xlon", ivarid)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error xlon variable undefined'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          istatus = nf90_get_var(iutin, ivarid, sp2d)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error reading xlon variable'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          xlon = transpose(sp2d)
          istatus = nf90_inq_varid(iutin, "xmap", ivarid)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error xmap variable undefined'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          istatus = nf90_get_var(iutin, ivarid, sp2d)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error reading xmap variable'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          xmap = transpose(sp2d)
          istatus = nf90_inq_varid(iutin, "dmap", ivarid)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error dmap variable undefined'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          istatus = nf90_get_var(iutin, ivarid, sp2d)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error reading dmap variable'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          dmap = transpose(sp2d)
          istatus = nf90_inq_varid(iutin, "coriol", ivarid)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error coriol variable undefined'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          istatus = nf90_get_var(iutin, ivarid, sp2d)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error reading coriol variable'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          f = transpose(sp2d)
          istatus = nf90_inq_varid(iutin, "snowam", ivarid)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error snowam variable undefined'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          istatus = nf90_get_var(iutin, ivarid, sp2d)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error reading snowam variable'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          do n = 1 , nnsg
            snw(n,:,:) = transpose(sp2d)
          end do
        end subroutine read_domain

        subroutine read_subdomain(ht1,lnd1,xlat1,xlon1)
          use mod_message
          use mod_constants , only : gti
          use netcdf
          implicit none

          real(8) , dimension(nnsg,iy,jx) , intent(out) :: ht1
          real(8) , dimension(nnsg,iy,jx) , intent(out) :: lnd1
          real(8) , dimension(nnsg,iy,jx) , intent(out) :: xlat1
          real(8) , dimension(nnsg,iy,jx) , intent(out) :: xlon1

          integer :: istatus , ivarid
          integer :: i , j , n , ii , jj
          real(4) , dimension(jxsg,iysg) :: sp2d1
          
          if (iutin1 < 0) then
            write (6,*) 'Error : Subdom file not in open state'
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if

          istatus = nf90_inq_varid(iutin1, "topo", ivarid)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error topo variable undefined'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'SUBDOMAIN FILE ERROR')
          end if
          istatus = nf90_get_var(iutin1, ivarid, sp2d1)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error reading topo variable'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'SUBDOMAIN FILE ERROR')
          end if
          do j = 1 , jxsg
            do i = 1 , iysg
              jj = mod(j,nsg)
              if ( jj.eq.0 ) jj = nsg
              ii = mod(i,nsg)
              if ( ii.eq.0 ) ii = nsg
              n = (jj-1)*nsg + ii
              jj = (j+nsg-1)/nsg
              ii = (i+nsg-1)/nsg
              ht1(n,ii,jj) = sp2d1(j,i)*gti
            end do
          end do
          istatus = nf90_inq_varid(iutin1, "landuse", ivarid)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error landuse variable undefined'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'SUBDOMAIN FILE ERROR')
          end if
          istatus = nf90_get_var(iutin1, ivarid, sp2d1)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error reading landuse variable'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'SUBDOMAIN FILE ERROR')
          end if
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
          istatus = nf90_inq_varid(iutin1, "xlat", ivarid)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error xlat variable undefined'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'SUBDOMAIN FILE ERROR')
          end if
          istatus = nf90_get_var(iutin1, ivarid, sp2d1)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error reading xlat variable'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'SUBDOMAIN FILE ERROR')
          end if
          do j = 1 , jxsg
            do i = 1 , iysg
              jj = mod(j,nsg)
              if ( jj.eq.0 ) jj = nsg
              ii = mod(i,nsg)
              if ( ii.eq.0 ) ii = nsg
              n = (jj-1)*nsg + ii
              jj = (j+nsg-1)/nsg
              ii = (i+nsg-1)/nsg
              xlat1(n,ii,jj) = sp2d1(j,i)
            end do
          end do
          istatus = nf90_inq_varid(iutin1, "xlon", ivarid)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error xlon variable undefined'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'SUBDOMAIN FILE ERROR')
          end if
          istatus = nf90_get_var(iutin1, ivarid, sp2d1)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error reading xlon variable'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'SUBDOMAIN FILE ERROR')
          end if
          do j = 1 , jxsg
            do i = 1 , iysg
              jj = mod(j,nsg)
              if ( jj.eq.0 ) jj = nsg
              ii = mod(i,nsg)
              if ( ii.eq.0 ) ii = nsg
              n = (jj-1)*nsg + ii
              jj = (j+nsg-1)/nsg
              ii = (i+nsg-1)/nsg
              xlon1(n,ii,jj) = sp2d1(j,i)
            end do
          end do

        end subroutine read_subdomain

        subroutine close_domain
          use mod_message
          use netcdf
          implicit none

          integer :: istatus

          if (iutin >= 0) then
            istatus = nf90_close(iutin)
            if (istatus /= nf90_noerr) then
              write (6,*) 'Error File close'
              write (6,*) nf90_strerror(istatus)
              call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
            end if
            iutin = -1
          end if
          if ( nsg>1 .and. iutin1 >=0 ) then
            istatus = nf90_close(iutin1)
            if (istatus /= nf90_noerr) then
              write (6,*) 'Error File close'
              write (6,*) nf90_strerror(istatus)
              call fatal(__FILE__,__LINE__, 'SUBDOMAIN FILE ERROR')
            end if
            iutin1 = -1
          end if

        end subroutine close_domain

        subroutine read_texture(nats,texture)
          use netcdf
          use mod_message
          implicit none

          integer , intent(in) :: nats
          real(8) , dimension(iy,jx,nats) , intent(out) :: texture

          integer :: istatus , ivarid
          integer :: i , j , n
          integer , dimension(3) :: istart , icount
          real(4), dimension(jx,iy) ::  toto

          if (iutin < 0) then
            istatus = nf90_open(dname, nf90_nowrite, iutin)
            if ( istatus /= nf90_noerr) then
              write (6,*) 'Error Opening Domain file ', trim(dname)
              write (6,*) nf90_strerror(istatus)
              call fatal(__FILE__,__LINE__, 'CANNOT OPEN DOMAIN FILE')
            end if
          end if
          istatus = nf90_inq_varid(iutin, "texture_fraction", ivarid)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error texture_fraction variable undefined'
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
          end if
          istart(2) = 1
          istart(1) = 1
          icount(3) = 1
          icount(2) = iy
          icount(1) = jx
          do n = 1 , nats
            istart(3) = n
            istatus = nf90_get_var(iutin, ivarid, toto, istart, icount)
            if (istatus /= nf90_noerr) then
              write (6,*) 'Error reading texture_fraction variable'
              write (6,*) nf90_strerror(istatus)
              call fatal(__FILE__,__LINE__, 'DOMAIN FILE ERROR')
            end if
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
          use netcdf
          use mod_message
          implicit none

          character(5) , dimension(ntr) , intent(in) :: chtrname
          real(8) , dimension(iy,jx,12,ntr) , intent(out) :: chemsrc

          integer :: ncid , ivarid , istatus
          real(4) , dimension(jx,iy) :: toto
          character(5) :: aerctl
          integer , dimension(3) :: istart , icount
          integer :: itr , i , j , m

          istatus = nf90_open(aername, nf90_nowrite, ncid)
          if ( istatus /= nf90_noerr) then
            write (6,*) 'Error Opening Aerosol file ', trim(aername)
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__,'AEROSOL FILE OPEN ERROR')
          end if

          do itr = 1 , ntr
            aerctl = chtrname(itr)
            write (aline, *) itr , aerctl
            call say
            if ( aerctl(1:4).ne.'DUST') then
              if ( aerctl(1:3).eq.'SO2' ) then
                if ( aertyp(4:4).eq.'1' ) then
                  istatus = nf90_inq_varid(ncid, "so2", ivarid)
                  if (istatus /= nf90_noerr) then
                    write (6,*) 'Error so2 variable undefined'
                    write (6,*) nf90_strerror(istatus)
                    call fatal(__FILE__,__LINE__,'AEROSOL SO2 ERROR')
                  end if
                  istatus = nf90_get_var(ncid, ivarid, toto)
                  if (istatus /= nf90_noerr) then
                    write (6,*) 'Error reading so2 variable'
                    write (6,*) nf90_strerror(istatus)
                    call fatal(__FILE__,__LINE__,'AEROSOL READ ERROR')
                  end if
                  do m = 1 , 12
                    do j = 1 , jx
                      do i = 1 , iy
                        chemsrc(i,j,m,itr) = toto(j,i)
                      end do
                    end do
                  end do
                end if
                if ( aertyp(5:5).eq.'1' ) then
                  istatus = nf90_inq_varid(ncid, "so2_monthly", ivarid)
                  if (istatus /= nf90_noerr) then
                    write (6,*) 'Error so2_monthly variable undefined'
                    write (6,*) nf90_strerror(istatus)
                    call fatal(__FILE__,__LINE__,'AEROSOL SO2M ERROR')
                  end if
                  istart(1) = 1
                  istart(2) = 1
                  icount(1) = jx
                  icount(2) = iy
                  icount(3) = 1
                  do m = 1 , 12
                    istart(3) = m
                    istatus = nf90_get_var(ncid,ivarid,toto, &
                                        &  istart,icount)
                    if (istatus /= nf90_noerr) then
                      write (6,*) 'Error reading so2_monthly variable'
                      write (6,*) nf90_strerror(istatus)
                      call fatal(__FILE__,__LINE__,  &
                               & 'AEROSOL SO2M READ ERROR')
                    end if
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
                  istatus = nf90_inq_varid(ncid, "bc", ivarid)
                  if (istatus /= nf90_noerr) then
                    write (6,*) 'Error bc variable undefined'
                    write (6,*) nf90_strerror(istatus)
                    call fatal(__FILE__,__LINE__,'AEROSOL BC ERROR')
                  end if
                  istatus = nf90_get_var(ncid, ivarid, toto)
                  if (istatus /= nf90_noerr) then
                    write (6,*) 'Error reading bc variable'
                    write (6,*) nf90_strerror(istatus)
                    call fatal(__FILE__,__LINE__,'AEROSOL READ ERROR')
                  end if
                  do m = 1 , 12
                    do j = 1 , jx
                      do i = 1 , iy
                        chemsrc(i,j,m,itr) = toto(j,i)
                      end do
                    end do
                  end do
                end if
                if ( aertyp(5:5).eq.'1' ) then
                  istatus = nf90_inq_varid(ncid, "bc_monthly", ivarid)
                  if (istatus /= nf90_noerr) then
                    write (6,*) 'Error bc_monthly variable undefined'
                    write (6,*) nf90_strerror(istatus)
                    call fatal(__FILE__,__LINE__,'AEROSOL BCM ERROR')
                  end if
                  istart(1) = 1
                  istart(2) = 1
                  icount(1) = jx
                  icount(2) = iy
                  icount(3) = 1
                  do m = 1 , 12
                    istart(3) = m
                    istatus = nf90_get_var(ncid,ivarid,toto,  &
                                       &   istart,icount)
                    if (istatus /= nf90_noerr) then
                      write (6,*) 'Error reading bc_monthly variable'
                      write (6,*) nf90_strerror(istatus)
                      call fatal(__FILE__,__LINE__, &
                               & 'AEROSOL READ ERROR')
                    end if
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
                  istatus = nf90_inq_varid(ncid, "oc", ivarid)
                  if (istatus /= nf90_noerr) then
                    write (6,*) 'Error oc variable undefined'
                    write (6,*) nf90_strerror(istatus)
                    call fatal(__FILE__,__LINE__,'AEROSOL OC ERROR')
                  end if
                  istatus = nf90_get_var(ncid, ivarid, toto)
                  if (istatus /= nf90_noerr) then
                    write (6,*) 'Error reading oc variable'
                    write (6,*) nf90_strerror(istatus)
                    call fatal(__FILE__,__LINE__,'AEROSOL READ ERROR')
                  end if
                  do m = 1 , 12
                    do j = 1 , jx
                      do i = 1 , iy
                        chemsrc(i,j,m,itr) = toto(j,i)
                      end do
                    end do
                  end do
                end if
                if ( aertyp(5:5).eq.'1' ) then
                  istatus = nf90_inq_varid(ncid, "oc_monthly", ivarid)
                  if (istatus /= nf90_noerr) then
                    write (6,*) 'Error oc_monthly variable undefined'
                    write (6,*) nf90_strerror(istatus)
                    call fatal(__FILE__,__LINE__,'AEROSOL OCM ERROR')
                  end if
                  istart(1) = 1
                  istart(2) = 1
                  icount(1) = jx
                  icount(2) = iy
                  icount(3) = 1
                  do m = 1 , 12
                    istart(3) = m
                    istatus = nf90_get_var(ncid,ivarid,toto, &
                                       &   istart,icount)
                    if (istatus /= nf90_noerr) then
                      write (6,*) 'Error reading oc_monthly variable'
                      write (6,*) nf90_strerror(istatus)
                      call fatal(__FILE__,__LINE__,  &
                               & 'AEROSOL READ ERROR')
                    end if
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
          if ( istatus /= nf90_noerr) then
            write (6,*) 'Error Closing Aerosol file ', trim(aername)
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__,'AEROSOL FILE CLOSE ERROR')
          end if

        end subroutine read_aerosol

      end module mod_ncio
