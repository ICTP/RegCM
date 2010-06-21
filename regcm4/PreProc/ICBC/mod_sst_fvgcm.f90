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

      module mod_sst_fvgcm

      contains

      subroutine sst_fvgcm

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Comments on dataset sources and location:                          c
!                                                                    c
! FVGCM    HadAMH_SST in the original netCDF format.                 c
!          for 'RF'          run, from 1959 to 1991, 385 months      c
!          for 'A2' and 'B2' run, from 2069 to 2101, 385 months      c
!                                                                    c
!          ML= 1 is   0.0; ML= 2 is   1.875; => ML=192 is 358.125E   c
!          NL= 1 is  90.0; ML= 2 is  88.75 ; => ML=145 is -90.       c
!                                                                    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use mod_sst_grid
      use mod_interp , only : bilinx
      use mod_printl
      use netcdf

      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: ilon = 192 , jlat = 145
!
! Local variables
!
      integer :: i , idatef , idateo , it , j , k , ludom , lumax ,     &
             &   mrec , nday , nmo , nyear , nsteps
      real(4) , dimension(jlat) :: lati
      real(4) , dimension(ilon) :: loni
      integer , dimension(20) :: lund
      real(4) , dimension(ilon,jlat) :: temp
      real(4) , dimension(ilon,jlat) :: sst
      integer :: idate , idate0
      logical :: there
      character(256) :: sstfile
!
      if ( ssttyp=='FV_RF' ) then
        inquire (file=trim(inpglob)//'/SST/Sst_1959_1991ref.dat',       &
            &    exist=there)
        if ( .not.there ) print * ,                                     &
            & 'Sst_1959_1991ref.dat is not available under ',           &
            & trim(inpglob),'/SST/'
        open (11,file=trim(inpglob)//'/SST/Sst_1959_1991ref.dat',       &
             &form='unformatted',recl=ilon*jlat*ibyte,access='direct')
      else if ( ssttyp=='FV_A2' ) then
        inquire (file=trim(inpglob)//'/SST/Sst_2069_2101_A2.dat',       &
             &   exist=there)
        if ( .not.there ) print * ,                                     &
            & 'Sst_2069_2101_A2.dat is not available under ',           &
            & trim(inpglob),'/SST/'
        open (11,file=trim(inpglob)//'/SST/Sst_2069_2101_A2.dat',       &
             &form='unformatted',recl=ilon*jlat*ibyte,access='direct')
      else if ( ssttyp=='FV_B2' ) then
        inquire (file=trim(inpglob)//'/SST/Sst_2069_2101_B2.dat',       &
              &  exist=there)
        if ( .not.there ) print * ,                                     &
             & 'Sst_2069_2101_B2.dat is not available under ',          &
             & trim(inpglob),'/SST/'
        open (11,file=trim(inpglob)//'/SST/Sst_2069_2101_B2.dat',       &
             &form='unformatted',recl=ilon*jlat*ibyte,access='direct')
      else
        write (*,*) 'PLEASE SET right SSTTYP in regcm.in'
        write (*,*) 'Supported types are FV_RF FV_A2 FV_B2'
        stop
      end if
      write (sstfile,99001) trim(dirglob), pthsep, trim(domname) ,      &
        &  '_SST.RCM'
      open (21,file=sstfile,form='unformatted',status='replace')
 
      idate = globidate1/10000
      if ( idate-(idate/100)*100==1 ) then
        idate = idate - 89
      else
        idate = idate - 1
      end if
      idateo = idate
      idate0 = idateo*10000 + 100
      idate = globidate2/10000
      if ( idate-(idate/100)*100==12 ) then
        idate = idate + 89
      else
        idate = idate + 1
      end if
      idatef = idate
      print * , globidate1 , globidate2 , idateo , idatef
      write (sstfile,99001) trim(dirglob), pthsep, trim(domname) ,      &
          &  '_RCM_SST.dat'
      open (25,file=sstfile,status='unknown',form='unformatted',        &
          & recl=iy*jx*ibyte,access='direct')
      if ( igrads==1 ) then
        write (sstfile,99001) trim(dirglob), pthsep, trim(domname) ,    &
         &   '_RCM_SST.ctl'
        open (31,file=sstfile,status='replace')
        write (31,'(a,a,a)') 'dset ^',trim(domname),'_RCM_SST.dat'
      end if
      nsteps = (idatef/100-idateo/100)*12 + (idatef-(idatef/100)*100)   &
               & - (idateo-(idateo/100)*100) + 1
      call setup_sstfile(idateo,nsteps)
      mrec = 0
 
!     ******    SET UP LONGITUDES AND LATITUDES FOR SST DATA
      do i = 1 , ilon
        loni(i) = float(i-1)*1.875
      end do
      do j = 1 , jlat
        lati(j) = -90. + 1.25*float(j-1)
      end do
 
!     **  REF  SST DATA, 1.875x1.1.25, AVAILABLE FROM 16/1/1959 TO
!     16/1/1991 ** A2&B2 SST DATA, 1.875x1.1.25, AVAILABLE FROM
!     16/1/2069 TO 16/1/2101
      idate = idateo
      do while ( idate<=idatef )
        nyear = idate/100
        nmo = idate - nyear*100
        nday = 16
        write (*,*) idate*10000 + 100 , idate0
!       IF(IDATE*10000+100.EQ.IDATE0) CALL SST_MN(SSTTYP)
        if ( ssttyp=='FV_RF' ) then
          it = (nyear-1959)*12 + nmo
        else
          it = (nyear-2069)*12 + nmo
        end if
        read (11,rec=it) temp
        do j = 1 , jlat
          do i = 1 , ilon
!           if(temp(I,NLAT+1-J).gt.-9000.0.and.
!           &         temp(I,NLAT+1-J).lt.10000.0) then
!           SST2(I,J)=temp(I,NLAT+1-J)
            if ( temp(i,j)>-9000.0 .and. temp(i,j)<10000.0 ) then
              sst(i,j) = temp(i,j)
            else
              sst(i,j) = -9999.
            end if
          end do
        end do
 
!       ******           PRINT OUT DATA AS A CHECK
        if ( nmo==1 ) call printl(sst,ilon,jlat)
 
        call bilinx(sst,sstmm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
        print * , 'XLON,XLAT,SST=' , xlon(1,1) , xlat(1,1) , sstmm(1,1)
 
        do j = 1 , jx
          do i = 1 , iy
            if ( sstmm(i,j)<-5000. .and.                                &
               & (lu(i,j)>13.5 .and. lu(i,j)<15.5) ) then
              do k = 1 , 20
                lund(k) = 0.0
              end do
              lund(nint(lu(i-1,j-1))) = lund(nint(lu(i-1,j-1))) + 2
              lund(nint(lu(i-1,j))) = lund(nint(lu(i-1,j))) + 3
              lund(nint(lu(i-1,j+1))) = lund(nint(lu(i-1,j+1))) + 2
              lund(nint(lu(i,j-1))) = lund(nint(lu(i,j-1))) + 3
              lund(nint(lu(i,j+1))) = lund(nint(lu(i,j+1))) + 3
              lund(nint(lu(i+1,j-1))) = lund(nint(lu(i+1,j-1))) + 2
              lund(nint(lu(i+1,j))) = lund(nint(lu(i+1,j))) + 3
              lund(nint(lu(i+1,j+1))) = lund(nint(lu(i+1,j+1))) + 2
              ludom = 18
              lumax = 0
              do k = 1 , 20
                if ( k<=13 .or. k>=16 ) then
                  if ( lund(k)>lumax ) then
                    ludom = k
                    lumax = lund(k)
                  end if
                end if
              end do
              lu(i,j) = float(ludom)
              print * , ludom , sstmm(i,j)
            end if
            if ( sstmm(i,j)>-100. ) then
              sstmm(i,j) = sstmm(i,j)
            else
              sstmm(i,j) = -9999.
            end if
          end do
        end do
 
!       ******           WRITE OUT SST DATA ON MM4 GRID
        write (21) nday , nmo , nyear , sstmm
        print * , 'WRITING OUT MM4 SST DATA:' , nmo , nyear
        idate = idate + 1
        if ( nmo==12 ) idate = idate + 88
        mrec = mrec + 1
        write (25,rec=mrec) ((sstmm(i,j),j=1,jx),i=1,iy)
      end do
      write (10,rec=4) ((lu(i,j),j=1,jx),i=1,iy)
 
99001 format (a,a,a,a)
      end subroutine sst_fvgcm
!
      end module mod_sst_fvgcm
