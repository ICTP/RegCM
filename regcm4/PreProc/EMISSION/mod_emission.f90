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

      module mod_emission

      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: ndims = 3
      integer , parameter :: nrecs = 12
      integer , parameter :: nlvls1 = 1
      integer , parameter :: nlats1 = 360
      integer , parameter :: nlons1 = 720
      integer , parameter :: nlvls2 = 1
      integer , parameter :: nlats2 = 180
      integer , parameter :: nlons2 = 360

      character(*) , parameter :: lvl_name = "level"
      character(*) , parameter :: lat_name = "lat"
      character(*) , parameter :: lon_name = "lon"
      character(*) , parameter :: rec_name = "Time"

      integer , parameter :: ema = 1 , emb = 2
      integer , parameter :: emsrc1 = 1 , emsrc2 = 2

      character(*) , parameter :: invntry1 = "retro" ,                  &
                                & invntry2 = "poet" ,                   &
                                & invntry3 = "gfed" ,                   &
                                & invntry4 = "edgar" ,                  &
                                & invntry5 = "mozart"

      character(*) , parameter ::                                       &
         & emiss_name1  = "emission_flux" ,                             &
         & emiss_name2  = "fire_emis" ,                                 &
         & emiss_name3  = "total" ,                                     &
         & emiss_name4  = "bio" ,                                       &
         & emiss_name5  = "oceans" ,                                    &
         & emiss_name6  = "totbb" ,                                     &
         & emiss_name7  = "bb_emission_flux" ,                          &
         & emiss_name8  = "anthro_emission_flux" ,                      &
         & emiss_name9  = "anthro" ,                                    &
         & emiss_name10 = "biogenic" ,                                  &
         & emiss_name11 = "ocean" ,                                     &
         & emiss_name12 = "bb" ,                                        &
         & units        = "units" ,                                     &
         & emiss_units  = "kg(species)/m2/s " ,                         &
         & lat_units    = "degrees_north" ,                             &
         & lon_units    = "degrees_east"

      character(256) :: emibase
      character(*) , parameter ::                                       &
         & datadir1 = "/EMISSION_INVENTORY/retro/" ,                    &
         & datadir2 = "/EMISSION_INVENTORY/poet/netcdf/" ,              &
         & datadir3 = "/EMISSION_INVENTORY/gfed/" ,                     &
         & datadir4 = "/EMISSION_INVENTORY/edgar/netcdf/" ,             &
         & datadir5 = "/EMISSION_INVENTORY/mozart/"

      logical , parameter :: iretro  = .true.
      logical , parameter :: ipoet   = .true.
      logical , parameter :: igfed   = .false.
      logical , parameter :: iedgar  = .true.
      logical , parameter :: imozart = .true.

      integer , parameter :: nspc1a = 7
      integer , parameter :: nspc1b = 10
      integer , parameter :: nspc2a = 5
      integer , parameter :: nspc2b = 15
      integer , parameter :: nspc3  = 1
      integer , parameter :: nspc4a = 4
      integer , parameter :: nspc4b = 4
      integer , parameter :: nspc5a = 1
      integer , parameter :: nspc5b = 2

      character(30) , dimension(nspc1a) :: ele_retroa
      character(30) , dimension(nspc1b) :: ele_retrob
      character(30) , dimension(nspc2a) :: ele_poeta
      character(30) , dimension(nspc2b) :: ele_poetb
      character(30) , dimension(nspc3)  :: ele_gfed
      character(30) , dimension(nspc4a) :: ele_edgara
      character(30) , dimension(nspc4b) :: ele_edgarb
      character(30) , dimension(nspc5a) :: ele_mozrta
      character(30) , dimension(nspc5b) :: ele_mozrtb

      integer :: ny , nx
      integer :: idate1 , idate2
      real(4) , allocatable , dimension(:,:) :: aermm , xlat , xlon

      contains

      subroutine init_emiss(iny,jnx,iidate1,iidate2,dirglob)
        implicit none
        integer , intent (in) :: iny , jnx , iidate1 , iidate2
        character(*) , intent(in) :: dirglob
        ny = iny
        nx = jnx
        idate1 = iidate1
        idate2 = iidate2
        emibase = trim(dirglob)
        allocate(aermm(ny,nx))
        allocate(xlat(ny,nx))
        allocate(xlon(ny,nx))
      end subroutine

      subroutine free_emiss
        implicit none
        deallocate(aermm)
        deallocate(xlat)
        deallocate(xlon)
      end subroutine

      subroutine reademission(nlats,nlons,nlvls,inv,dr,src,nnn,emiss_in,&
                            & loni,lati,yy,month,recc)
      use netcdf
      implicit none
!
! Dummy arguments
!
      character(3) :: dr
      character(1) :: inv
      integer :: month , nlats , nlons , nlvls , nnn , recc , src , yy
      real(4) , dimension(nlons,nlats,nlvls) :: emiss_in
      real(4) , dimension(nlats) :: lati
      real(4) , dimension(nlons) :: loni
      intent (in) dr , inv , month , nlvls , nnn , src
      intent (inout) recc , yy
!
! Local variables
!
      real(4) :: avo , cfact , mw
      integer :: ayr1 , ayr2 , dt , emiss_varid , emsrc , i , idy1 ,    &
               & idy2 , ihr1 , ihr2 , imo1 , imo2 , nyr1 , nyr2 , j ,   &
               & jul1900 , julday , julncep , lat_varid , lon_varid ,   &
               & ncid , ndims_in , ngatts_in , nspc , nvars_in ,        &
               & retval , spc , unlimdimid_in
      integer , dimension(ndims) :: icount , istart
      character(80) :: datadir
      character(20) :: emiss_name
      character(100) :: file_name , cname
      character(6) :: invntry
      character(30) , dimension(nnn) :: elements
!
!     All information about the netcdf file and parameter inside
!     read.parm
!     This is the name of the data file we will read.
!     FILE_NAME is the full path name
!     NAME is the name of the file itself
!     We are reading 3D data,lat X lon X time
!     The istart and icount arrays will tell the netCDF library where to
!     read our data.
!     In addition to the latitude and longitude dimensions, we will also
!     create latitude and longitude variables which will hold the actual
!     latitudes and longitudes. Since they hold data about the
!     coordinate system, the netCDF term for these is: "coordinate
!     variables."
!     We will read emission fields. In netCDF
!     terminology these are called "variables."
!     Program variables to hold the data we will read in. We will only
!     need enough space to hold one timestep of data; one record.
 
      integer :: year
      character(30) :: element

      namelist /species/ ele_retroa , ele_retrob , ele_poeta ,          &
      & ele_poetb , ele_gfed , ele_edgara , ele_edgarb , ele_mozrta ,   &
      & ele_mozrtb
 
!     OPEN(30,file='namelist.input',status='old')
      rewind (30)
      read (30,nml=species)
 
      avo = 6.022E23
      mw = 30.
      cfact = mw*10./avo
 
      nspc = 0
      if ( inv=='a' .and. dr=='DR1' .and. src==1 ) then
        invntry = invntry1
        datadir = datadir1
        nspc = nspc1a
      else if ( inv=='a' .and. dr=='DR1' .and. src==2 ) then
        invntry = invntry1
        datadir = datadir1
        nspc = nspc1b
      else if ( inv=='b' .and. dr=='DR2' .and. src==1 ) then
        invntry = invntry2
        datadir = datadir2
        nspc = nspc2a
      else if ( inv=='b' .and. dr=='DR2' .and. src==2 ) then
        invntry = invntry2
        datadir = datadir2
        nspc = nspc2b
      else if ( inv=='c' .and. dr=='DR3' ) then
        invntry = invntry3
        datadir = datadir3
        nspc = nspc3
      else if ( inv=='d' .and. dr=='DR4' .and. src==1 ) then
        invntry = invntry4
        datadir = datadir4
        nspc = nspc4a
      else if ( inv=='d' .and. dr=='DR4' .and. src==2 ) then
        invntry = invntry4
        datadir = datadir4
        nspc = nspc4b
      else if ( inv=='m' .and. dr=='DR5' .and. src==1 ) then
        invntry = invntry5
        datadir = datadir5
        nspc = nspc5a
      else if ( inv=='m' .and. dr=='DR5' .and. src==2 ) then
        invntry = invntry5
        datadir = datadir5
        nspc = nspc5b
      else
      end if
      if ( src==1 ) emsrc = emsrc1
      if ( src==2 ) emsrc = emsrc2
      dt = 1
!     WRITE(*,*)'|--------------------------------------------|'
!     WRITE(*,*)'|            WELCOME                         |'
!     WRITE(*,*)'|              TO                            |'
!     WRITE(*,*)'|    EMISSION INVENTORIES DATA RETRIEVAL     |'
!     WRITE(*,*)'|              SYSTEM                        |'
!     WRITE(*,*)'|    ICTP, 2007, ESP(PWCg)                   |'
!     WRITE(*,*)'|____________________________________________|_'
 
!     WRITE(*,*)NLONS,NLATS
!     WRITE(*,*)INVNTRY
!     WRITE(*,*)DATADIR
!---------------------------------------------------------------
      call julian2(idate1,dt,julday,julncep,jul1900,nyr1,imo1,idy1,ihr1,&
                 & ayr1)
 
      call julian2(idate2,dt,julday,julncep,jul1900,nyr2,imo2,idy2,ihr2,&
                 & ayr2)
 
!ah   IF (nyr1 .ge. 1960 .and. nyr2 .le. 2000) THEN
!     WRITE(*,*)'THE DATA will be RETRIEVE FROM RETRO INVENTORY'
!ah   ELSE IF (nyr1 .ge. 2000 .and. nyr2 .le. 2005) THEN
!     WRITE(*,*)'THE Bio. DATA will be RETRIEVE FROM GFED INVENTORY'
!ah   ELSE IF (nyr1 .ge. 1960 .and. nyr2 .le. 2005) THEN
!     WRITE(*,*)'THE DATA will be RETRIEVED FROM RETRO + GFED'
!ah   ELSE
!     WRITE(*,*)'YOUR DATA IS OUTSIDE THE DATA RANGE'
!ah   STOP
!ah   END IF
 
!     Open the files.
      if ( emsrc==1 ) then
        open (11,file='anth_'//invntry,status='old')
        open (12,file='anth_spc_'//invntry,status='old')
        write (*,*) 'YOU CHOOSE ANTHRO EMISSION' , 'anth_'//invntry
      else if ( emsrc==2 ) then
        open (11,file='bio_'//invntry,status='old')
        open (12,file='bio_spc_'//invntry,status='old')
        write (*,*) 'YOU CHOOSE FIRE EMISSION'
      else
      end if
!     REWIND(11)
!     REWIND(12)
      do           !loop over species file
        read (11,'(A)',end=100) cname
        read (12,fmt=99001,end=100) year , element
        do spc = 1 , nspc    !loop over species
 
          file_name = trim(emibase)//trim(datadir)//cname
!         WRITE(*,*)YEAR,ELEMENT,FILE_NAME

          if ( emsrc==1 .and. invntry=='retro' ) then
            elements(spc) = ele_retroa(spc)
            if ( yy>2000 ) yy = 2000
          else if ( emsrc==1 .and. invntry=='poet' ) then
            elements(spc) = ele_poeta(spc)
            if ( yy>2000 ) yy = 2000
          else if ( emsrc==1 .and. invntry=='edgar' ) then
            elements(spc) = ele_edgara(spc)
            yy = 2000
          else if ( emsrc==1 .and. invntry=='mozart' ) then
            elements(spc) = ele_mozrta(spc)
            yy = year
          else if ( emsrc==2 .and. invntry=='gfed' ) then
            elements(spc) = ele_gfed(spc)
          else if ( emsrc==2 .and. invntry=='edgar' ) then
            elements(spc) = ele_edgarb(spc)
            yy = 2000
          else if ( emsrc==2 .and. invntry=='retro' ) then
            elements(spc) = ele_retrob(spc)
            if ( yy>2000 ) yy = 2000
          else if ( emsrc==2 .and. invntry=='poet' ) then
            elements(spc) = ele_poetb(spc)
            if ( yy>=2003 ) yy = 2003
          else if ( emsrc==2 .and. invntry=='mozart' ) then
            elements(spc) = ele_mozrtb(spc)
            yy = year
          else
          end if

!ah       write(*,*)'HHHH',ELEMENT,ELEMENTS(spc),YEAR,yy
 
          if ( (element==elements(spc)) .and. (year==yy) ) then
 
            write (*,*) year , element , file_name
 
            retval = nf90_open(file_name,nf90_nowrite,ncid)
            write (*,*) 'logi'
            if ( retval/=nf90_noerr ) call handle_err(retval)
 
            retval = nf90_inquire(ncid,ndims_in,nvars_in,ngatts_in,     &
                   & unlimdimid_in)
            if ( retval/=nf90_noerr ) call handle_err(retval)
!           WRITE(*,*)ncid, ndims_in, nvars_in, ngatts_in,unlimdimid_in
!           WRITE(*,*)retval,nf90_noerr
 
!           Get the varids of the latitude and longitude coordinate
!           variables.
            retval = nf90_inq_varid(ncid,lat_name,lat_varid)
            if ( retval/=nf90_noerr ) call handle_err(retval)
            retval = nf90_inq_varid(ncid,lon_name,lon_varid)
            if ( retval/=nf90_noerr ) call handle_err(retval)
!-----------Spectial case of anthropogenic RETRO---------
 
            retval = nf90_get_var(ncid,lat_varid,lati)
            retval = nf90_get_var(ncid,lon_varid,loni)
 
            write (*,*) 'koko' , loni(1) , lati(1)
!           Get the varids of the emission in  netCDF variables.
            if ( emsrc==1 .and. invntry=='retro' ) then
              emiss_name = emiss_name1
              write (*,*) emiss_name
              retval = nf90_inq_varid(ncid,emiss_name,emiss_varid)
              if ( retval/=nf90_noerr ) call handle_err(retval)
            else if ( emsrc==2 .and. invntry=='retro' ) then
              emiss_name = emiss_name2
              write (*,*) emiss_name
              retval = nf90_inq_varid(ncid,emiss_name,emiss_varid)
              if ( retval/=nf90_noerr ) call handle_err(retval)
            else
            end if
 
            if ( emsrc==1 .and. invntry=='poet' ) then
              emiss_name = emiss_name3
              write (*,*) emiss_name
              retval = nf90_inq_varid(ncid,emiss_name,emiss_varid)
              if ( retval/=nf90_noerr ) call handle_err(retval)
            end if
 
            if ( emsrc==2 .and. invntry=='poet' ) then
              if ( element(1:4)=='bio_' ) then
                emiss_name = emiss_name4
                write (*,*) emiss_name
              else if ( element(1:2)=='o_' ) then
                emiss_name = emiss_name5
                write (*,*) emiss_name
              else
                emiss_name = emiss_name6
              end if
              retval = nf90_inq_varid(ncid,emiss_name,emiss_varid)
              if ( retval/=nf90_noerr ) call handle_err(retval)
            end if
 
            if ( emsrc==2 .and. invntry=='gfed' ) then
              emiss_name = emiss_name7
              write (*,*) emiss_name
              retval = nf90_inq_varid(ncid,emiss_name,emiss_varid)
              if ( retval/=nf90_noerr ) call handle_err(retval)
            end if
 
            if ( emsrc==1 .and. invntry=='edgar' ) then
              emiss_name = emiss_name8
              write (*,*) emiss_name
              retval = nf90_inq_varid(ncid,emiss_name,emiss_varid)
              if ( retval/=nf90_noerr ) call handle_err(retval)
            end if
 
            if ( emsrc==2 .and. invntry=='edgar' ) then
              emiss_name = emiss_name7
              write (*,*) emiss_name
              retval = nf90_inq_varid(ncid,emiss_name,emiss_varid)
              if ( retval/=nf90_noerr ) call handle_err(retval)
            end if
 
            if ( emsrc==1 .and. invntry=='mozart' ) then
              emiss_name = emiss_name9
              write (*,*) emiss_name
              write (*,*) 'XXXXXXXXXXXXXXXXX'
              retval = nf90_inq_varid(ncid,emiss_name,emiss_varid)
              if ( retval/=nf90_noerr ) call handle_err(retval)
            end if
 
            if ( emsrc==2 .and. invntry=='mozart' ) then
              if ( element(1:4)=='bio_' ) then
                emiss_name = emiss_name10
                write (*,*) emiss_name
              else if ( element(1:2)=='o_' ) then
                emiss_name = emiss_name11
                write (*,*) emiss_name
              else
              end if
              retval = nf90_inq_varid(ncid,emiss_name,emiss_varid)
              if ( retval/=nf90_noerr ) call handle_err(retval)
            end if
 
!           Read 1 record of NLVLS*NLATS*NLONS values, starting at the
!           beginning of the record (the (1, 1, 1, rec) element in the
!           netCDF file).
            icount(1) = nlons
            icount(2) = nlats
            icount(3) = 1
            istart(1) = 1
            istart(2) = 1
 
!           Read the surface emission data from the file, one
!           record at a time.
!-----------------read/write part----------------------------------
            istart(3) = month
            retval = nf90_get_var(ncid,emiss_varid,emiss_in,            &
                     &            istart,icount)
            if ( retval/=nf90_noerr ) call handle_err(retval)
            recc = recc + 1
            write (*,*) year , month , recc , element

            call bilinxo(emiss_in,loni,lati,nlons,nlats,aermm,xlon,xlat,&
                       & ny,nx,1)
 
            if ( invntry=='poet' ) then
              cfact = molwgt(element)*10/avo
              write (*,*) 'kkkkk' , element , molwgt(element) , cfact
              write (15,rec=recc) ((cfact*aermm(i,j),j=1,nx),i=1,ny)
 
            else if ( invntry=='mozart' ) then
              cfact = molwgt(element)*10/avo
              write (*,*) 'kkkkk' , element , molwgt(element) , cfact
              write (15,rec=recc) ((cfact*aermm(i,j),j=1,nx),i=1,ny)
 
            else
              write (15,rec=recc) ((aermm(i,j),j=1,nx),i=1,ny)
            end if
!           write(*,*)AERMM
!-------------------------------------------------------------------
!           Close the file. This frees up any internal netCDF resources
!           associated with the file.
 
            retval = nf90_close(ncid)
            if ( retval/=nf90_noerr ) call handle_err(retval)
            write (*,*) 'finish file==' , cname , ncid
            write (*,*) 'last record==' , recc
 
          end if
             !element/year criteria
 
        end do  !next species
 
      end do    !next file name
 
 
!     If we got this far, everything worked as expected. Yipee!
 100  continue
      print * , '*** SUCCESS reading*** '
      write (*,*) 'koooooo'
      close (11)
      close (12)
      close (15)

99001 format (i4,2x,a30)

      end subroutine reademission

      function molwgt(ele)
      implicit none
!
! Dummy arguments
!
      character(len=*) , intent(in) :: ele
      real(4) :: molwgt
!
      if ( ele=='acet' ) molwgt = 58.08
      if ( ele=='butane' ) molwgt = 58.12
      if ( ele=='butene' ) molwgt = 56.
      if ( ele=='c2h4' ) molwgt = 28.
      if ( ele=='c2h5oh' ) molwgt = 28.0
      if ( ele=='c2h6' ) molwgt = 30.07
      if ( ele=='c3h6' ) molwgt = 42.
      if ( ele=='c3h8' ) molwgt = 44.1
      if ( ele=='ch2o' ) molwgt = 30.0
      if ( ele=='ch3cho' ) molwgt = 44.0
      if ( ele=='ch3oh' ) molwgt = 32.04
      if ( ele=='co' ) molwgt = 28.0
      if ( ele=='mek' ) molwgt = 28.0
      if ( ele=='nox' ) molwgt = 30.0
      if ( ele=='toluene' ) molwgt = 28.0
      if ( ele=='bio_acet' ) molwgt = 58.08
      if ( ele=='bio_c2h4' ) molwgt = 28.0
      if ( ele=='bio_c2h6' ) molwgt = 30.07
      if ( ele=='bio_c3h6' ) molwgt = 44.1
      if ( ele=='bio_c3h8' ) molwgt = 44.1
      if ( ele=='bio_ch3oh' ) molwgt = 32.04
      if ( ele=='bio_co' ) molwgt = 28.0
      if ( ele=='bio_isop' ) molwgt = 68.0
      if ( ele=='bio_nox' ) molwgt = 30.0
      if ( ele=='bio_terp' ) molwgt = 136.
      if ( ele=='o_c2h4' ) molwgt = 28.0
      if ( ele=='o_c2h6' ) molwgt = 30.07
      if ( ele=='o_c3h6' ) molwgt = 42.0
      if ( ele=='o_c3h8' ) molwgt = 44.1
      if ( ele=='o_co' ) molwgt = 28.0
      if ( ele=='NH3' ) molwgt = 17.0
      if ( ele=='bio_NH3' ) molwgt = 17.0
      if ( ele=='o_NH3' ) molwgt = 17.0
 
      end function molwgt
!
      subroutine handle_err(errcode)
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: errcode
!
      print * , 'Error: ' , nf90_strerror(errcode)
      stop 2
      end subroutine handle_err
!
      end module mod_emission
