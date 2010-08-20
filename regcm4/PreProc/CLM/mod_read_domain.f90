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

    module mod_read_domain

    use netcdf
    use mod_dynparam

    implicit none

    real(4) :: clatx , clonx , dsx , grdfacx , offset ,               &
         & perr , platx , plonx , pmax , ptopx , xscale , xlatmax ,   &
         & xlatmin , xlonmax , xlonmin
    integer :: i , ibigendx , idatex , idin , idout , idy , ierr ,    &
             & ifield , ifld , igradsx , ihr , imap , imo , jotyp ,   &
             & irec , iyr , iyy , j , julnc , jxx , kmax , kzz , l
    character(256) :: terfile
    character(6) :: iprojx

    real(4) , allocatable , dimension(:,:) :: xlat , xlon, xlat_dum, xlon_dum
    real(4) , allocatable , dimension(:) :: xlat1d
    real(4) , allocatable , dimension(:) :: xlon1d
    real(4) , allocatable , dimension(:) :: sigx

    integer :: istatus
    integer :: incin
    integer :: idimid
    integer :: ivarid


    contains

    subroutine read_domain()
      implicit none

        !Open the netcdf file
        call handle_nc_err( nf90_open(terfile, nf90_nowrite, incin),   &
          "Opening", trim(terfile))

        !Read the dimensions from the netcdf file
        call handle_nc_err( nf90_inq_dimid(incin, "iy", idimid),   &
          "Finding","iy")
        call handle_nc_err( nf90_inquire_dimension(incin, idimid, len=iyy),    &
          "Reading","iy")

        call handle_nc_err( nf90_inq_dimid(incin, "jx", idimid),   &
          "Finding","jx")
        call handle_nc_err( nf90_inquire_dimension(incin, idimid, len=jxx),    &
          "Reading","jx")

        call handle_nc_err( nf90_inq_dimid(incin, "kz", idimid),   &
          "Finding","kz")
        call handle_nc_err( nf90_inquire_dimension(incin, idimid, len=kzz),    &
          "Reading","kz")

        !Check for consistency with regcm.in
        if ( iyy/=iy .or. jxx/=jx .or. kzz/=kz+1 ) then
          print * , "DOMAIN.INFO is inconsistent with regcm.in"
          print * , "  namelist  : " , iy , jx , kz
          print * , "  DOMAIN    : " , iyy , jxx , kzz
          stop "Dimension mismatch"
        end if

        !Read ds
        call handle_nc_err(  &
           nf90_get_att(incin, NF90_GLOBAL,"grid_size_in_meters", dsx),   &
          "Reading","ds")
        !Convert from m to km
        dsx = dsx/1000

        !Read clatx
        call handle_nc_err(  &
          nf90_get_att(incin, NF90_GLOBAL,"latitude_of_projection_origin" &
                       ,clatx),&
          "Reading","latitude_of_projection_origin")

        !Read clonx
        call handle_nc_err(  &
          nf90_get_att(incin, NF90_GLOBAL,"longitude_of_projection_origin", &
                       clonx),&
          "Reading","longitude_of_projection_origin")

        !Read iproj
        call handle_nc_err(  &
           nf90_get_att(incin, NF90_GLOBAL,"projection", iprojx), &
          "Reading","projection")

        !Only if using the Rotated Mercator projection, read the poles
        if(iprojx.eq."ROTMER")then
          !Read plat
          call handle_nc_err(  &
             nf90_get_att(incin, NF90_GLOBAL,"latitude_of_projection_pole", &
                          platx), &
            "Reading","latitude_of_projection_pole")
          !Read plon
          call handle_nc_err(  &
             nf90_get_att(incin, NF90_GLOBAL,"longitude_of_projection_pole", &
                          plonx), &
            "Reading","longitude_of_projection_pole")
        else
          platx = clatx
          plonx = clonx
        endif

        !Read grdfacx
        call handle_nc_err(  &
           nf90_get_att(incin, NF90_GLOBAL,"grid_factor",grdfacx), &
          "Reading","projection")

        !Read grdfacx
        call handle_nc_err(  &
           nf90_get_att(incin, NF90_GLOBAL,"grid_factor",grdfacx), &
          "Reading","projection")

        !Assume big endian for now
        !TODO: Endianness should be deleted as a variable when the model deals
        !full in netcdf
        ibigendx = 1

        !Read sigx
        call handle_nc_err(  &
          nf90_inq_varid(incin,"sigma",ivarid), &
          "Finding","sigma")
        call handle_nc_err(  &
          nf90_get_var(incin,ivarid,sigx), &
          "Reading","sigma")

        !Read ptopx
        call handle_nc_err(  &
          nf90_inq_varid(incin,"ptop",ivarid), &
          "Finding","ptop")
        call handle_nc_err(  &
          nf90_get_var(incin,ivarid,ptopx), &
          "Reading","ptop")

        !Read ptopx
        call handle_nc_err(  &
          nf90_inq_varid(incin,"ptop",ivarid), &
          "Finding","ptop")
        call handle_nc_err(  &
          nf90_get_var(incin,ivarid,ptopx), &
          "Reading","ptop")

        !Read xlat
        call handle_nc_err(  &
          nf90_inq_varid(incin,"xlat",ivarid), &
          "Finding","xlat")
        call handle_nc_err(  &
          nf90_get_var(incin,ivarid,xlat_dum), &
          "Reading","xlat")
          
        !Read xlon
        call handle_nc_err(  &
          nf90_inq_varid(incin,"xlon",ivarid), &
          "Finding","xlon")
        call handle_nc_err(  &
          nf90_get_var(incin,ivarid,xlon_dum), &
          "Reading","xlon")

        !Set xlat and xlon, swapping the i/j indicies of what was read in
        do i = 1,iy
        do j = 1,jx
          xlat(i,j) = xlat_dum(j,i)
          xlon(i,j) = xlon_dum(j,i)
        enddo
        enddo

        !Close the netcdf flie
        istatus = nf90_close(incin)
        if (istatus /= nf90_noerr) then
          write (6,*) "Error closing Domain file ", trim(terfile)
          write (6,*) nf90_strerror(istatus)
          stop
        end if

      end subroutine read_domain

      !A routine to simply handle NetCDF errors
      subroutine handle_nc_err(incerr,sAction,sVarname)
        implicit none
        integer,intent(in) :: incerr
        character*(*),intent(in) :: sAction,sVarname

         if (incerr /= nf90_noerr) then
          write (6,*) "Error associated with ", trim(sAction)," ",trim(sVarname)
          write (6,*) nf90_strerror(incerr)
          stop
        end if

      end subroutine handle_nc_err

    end module mod_read_domain
