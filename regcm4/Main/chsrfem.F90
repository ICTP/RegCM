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
 
      subroutine chsrfem

      use netcdf
      use mod_dynparam
      use mod_mainchem
      use mod_trachem
      use mod_dust
      use mod_iunits
      use mod_aerosol , only : dgmix , dssamix , dextmix
      use mod_message
#ifdef MPP1
      use mod_mppio
#ifndef IBM
      use mpi
#else 
      include 'mpif.h'
#endif 
#ifdef CLM
!      use surfrdMod , only : clm_getsoitex
!      use clm_varsur, only : clm_soitex
#endif
#endif
      implicit none
!
! Local variables
!
      character(5) :: aerctl
      integer :: i , itr , j , k , l , m , n
      integer :: istatus , ncid , ivarid
      integer , dimension(3) :: istart , icount
      logical :: there
      real(4) , dimension(jx,iy) :: toto
#ifdef MPP1
      integer :: ierr
#endif
!

! fisrt activate dust initialization

     write (aline, *) 'Calling inidust'
     call say
     call inidust

! read the monthly aerosol emission files

#ifdef MPP1
      if ( myid.eq.0 ) then
        if (aertyp(4:5).ne.'00') then
          call inaero
          istatus = nf90_open(ffin, nf90_nowrite, ncid)
          if ( istatus /= nf90_noerr) then
            write (6,*) 'Error Opening Aerosol file ', trim(ffin)
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__,'AEROSOL FILE OPEN ERROR')
          end if
        end if
        do itr = 1 , ntr
          aerctl = chtrname(itr)
          write (aline, *) itr , aerctl
          call say
          if ( aerctl(1:4).ne.'DUST' .and. aertyp(4:5).ne.'00' ) then
            if ( aerctl(1:3).eq.'SO2' ) then
              if ( aertyp(4:4).eq.'1' ) then
                istatus = nf90_inq_varid(ncid, "so2", ivarid)
                if (istatus /= nf90_noerr) then
                  write (6,*) 'Error so2 variable undefined'
                  write (6,*) nf90_strerror(istatus)
                  call fatal(__FILE__,__LINE__,'AEROSOL SO2 FIND ERROR')
                end if
                istatus = nf90_get_var(ncid, ivarid, toto)
                if (istatus /= nf90_noerr) then
                  write (6,*) 'Error reading so2 variable'
                  write (6,*) nf90_strerror(istatus)
                  call fatal(__FILE__,__LINE__,'AEROSOL SO2 READ ERROR')
                end if
                do m = 1 , 12
                  do j = 1 , jx
                    do i = 1 , iy
                      chemsrc_io(i,j,m,itr) = toto(j,i)
                    end do
                  end do
                end do
              else
                do m = 1 , 12
                  do j = 1 , jx
                    do i = 1 , iy
                      chemsrc_io(i,j,m,itr) = 0.0D0
                    end do
                  end do
                end do
              end if
              if ( aertyp(5:5).eq.'1' ) then
                istatus = nf90_inq_varid(ncid, "so2_monthly", ivarid)
                if (istatus /= nf90_noerr) then
                  write (6,*) 'Error so2_monthly variable undefined'
                  write (6,*) nf90_strerror(istatus)
                  call fatal(__FILE__,__LINE__,'AEROSOL VAR FIND ERROR')
                end if
                istart(1) = 1
                istart(2) = 1
                icount(1) = jx
                icount(2) = iy
                icount(3) = 1
                do m = 1 , 12
                  istart(3) = m
                  istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
                  if (istatus /= nf90_noerr) then
                    write (6,*) 'Error reading so2_monthly variable'
                    write (6,*) nf90_strerror(istatus)
                    call fatal(__FILE__,__LINE__,                       &
                             & 'AEROSOL VAR READ ERROR')
                  end if
                  do j = 1 , jx
                    do i = 1 , iy
                      chemsrc_io(i,j,m,itr) = chemsrc_io(i,j,m,itr)     &
                      & + toto(j,i)
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
                  call fatal(__FILE__,__LINE__,'AEROSOL BC FIND ERROR')
                end if
                istatus = nf90_get_var(ncid, ivarid, toto)
                if (istatus /= nf90_noerr) then
                  write (6,*) 'Error reading bc variable'
                  write (6,*) nf90_strerror(istatus)
                  call fatal(__FILE__,__LINE__,'AEROSOL BC READ ERROR')
                end if
                do m = 1 , 12
                  do j = 1 , jx
                    do i = 1 , iy
                      chemsrc_io(i,j,m,itr) = toto(j,i)
                    end do
                  end do
                end do
              else
                do m = 1 , 12
                  do j = 1 , jx
                    do i = 1 , iy
                      chemsrc_io(i,j,m,itr) = 0.0D0
                    end do
                  end do
                end do
              end if
              if ( aertyp(5:5).eq.'1' ) then
                istatus = nf90_inq_varid(ncid, "bc_monthly", ivarid)
                if (istatus /= nf90_noerr) then
                  write (6,*) 'Error bc_monthly variable undefined'
                  write (6,*) nf90_strerror(istatus)
                  call fatal(__FILE__,__LINE__,'AEROSOL VAR FIND ERROR')
                end if
                istart(1) = 1
                istart(2) = 1
                icount(1) = jx
                icount(2) = iy
                icount(3) = 1
                do m = 1 , 12
                  istart(3) = m
                  istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
                  if (istatus /= nf90_noerr) then
                    write (6,*) 'Error reading bc_monthly variable'
                    write (6,*) nf90_strerror(istatus)
                    call fatal(__FILE__,__LINE__,                       &
                             & 'AEROSOL VAR READ ERROR')
                  end if
                  do j = 1 , jx
                    do i = 1 , iy
                      chemsrc_io(i,j,m,itr) = chemsrc_io(i,j,m,itr)     &
                      & + toto(j,i)
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
                  call fatal(__FILE__,__LINE__,'AEROSOL OC FIND ERROR')
                end if
                istatus = nf90_get_var(ncid, ivarid, toto)
                if (istatus /= nf90_noerr) then
                  write (6,*) 'Error reading oc variable'
                  write (6,*) nf90_strerror(istatus)
                  call fatal(__FILE__,__LINE__,'AEROSOL OC READ ERROR')
                end if
                do m = 1 , 12
                  do j = 1 , jx
                    do i = 1 , iy
                      chemsrc_io(i,j,m,itr) = toto(j,i)
                    end do
                  end do
                end do
              else
                do m = 1 , 12
                  do j = 1 , jx
                    do i = 1 , iy
                      chemsrc_io(i,j,m,itr) = 0.0D0
                    end do
                  end do
                end do
              end if
              if ( aertyp(5:5).eq.'1' ) then
                istatus = nf90_inq_varid(ncid, "oc_monthly", ivarid)
                if (istatus /= nf90_noerr) then
                  write (6,*) 'Error oc_monthly variable undefined'
                  write (6,*) nf90_strerror(istatus)
                  call fatal(__FILE__,__LINE__,'AEROSOL VAR FIND ERROR')
                end if
                istart(1) = 1
                istart(2) = 1
                icount(1) = jx
                icount(2) = iy
                icount(3) = 1
                do m = 1 , 12
                  istart(3) = m
                  istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
                  if (istatus /= nf90_noerr) then
                    write (6,*) 'Error reading oc_monthly variable'
                    write (6,*) nf90_strerror(istatus)
                    call fatal(__FILE__,__LINE__,                       &
                             & 'AEROSOL VAR READ ERROR')
                  end if
                  do j = 1 , jx
                    do i = 1 , iy
                      chemsrc_io(i,j,m,itr) = chemsrc_io(i,j,m,itr)     &
                      & + toto(j,i)
                    end do
                  end do
                end do
              end if
            end if
          end if
        end do
        do j = 1 , jx
          do itr = 1 , ntr
            do m = 1 , 12
              do i = 1 , iy
                src_0(i,m,itr,j) = chemsrc_io(i,j,m,itr)
              end do
            end do
          end do
        end do
      end if
!     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call mpi_scatter(src_0(1,1,1,1),iy*12*ntr*jxp,mpi_real8,          &
                     & src0(1,1,1,1), iy*12*ntr*jxp,mpi_real8,          &
                     & 0,mpi_comm_world,ierr)
      do j = 1 , jendl
        do itr = 1 , ntr
          do m = 1 , 12
            do i = 1 , iy
              chemsrc(i,j,m,itr) = src0(i,m,itr,j)
            end do
          end do
        end do
      end do
      if (myid == 0) then
        if (aertyp(4:5).ne.'00') then
          istatus = nf90_close(ncid)
          if ( istatus /= nf90_noerr) then
            write (6,*) 'Error Closing Aerosol file ', trim(ffin)
            write (6,*) nf90_strerror(istatus)
            call fatal(__FILE__,__LINE__,'AEROSOL FILE CLOSE ERROR')
          end if
        end if
      end if
#else
      if (aertyp(4:5).ne.'00') then
        call inaero
        istatus = nf90_open(ffin, nf90_nowrite, ncid)
        if ( istatus /= nf90_noerr) then
          write (6,*) 'Error Opening Aerosol file ', trim(ffin)
          write (6,*) nf90_strerror(istatus)
          call fatal(__FILE__,__LINE__,'AEROSOL FILE OPEN ERROR')
        end if
      end if
      do itr = 1 , ntr
        aerctl = chtrname(itr)
        write (aline, *) itr , aerctl
        call say
        if ( aerctl(1:4).ne.'DUST' .and. aertyp(4:5).ne.'00' ) then
          if ( aerctl(1:3).eq.'SO2' ) then
            if ( aertyp(4:4).eq.'1' ) then
              istatus = nf90_inq_varid(ncid, "so2", ivarid)
              if (istatus /= nf90_noerr) then
                write (6,*) 'Error so2 variable undefined'
                write (6,*) nf90_strerror(istatus)
                call fatal(__FILE__,__LINE__,'AEROSOL SO2 FIND ERROR')
              end if
              istatus = nf90_get_var(ncid, ivarid, toto)
              if (istatus /= nf90_noerr) then
                write (6,*) 'Error reading so2 variable'
                write (6,*) nf90_strerror(istatus)
                call fatal(__FILE__,__LINE__,'AEROSOL SO2 READ ERROR')
              end if
              do m = 1 , 12
                do j = 1 , jx
                  do i = 1 , iy
                    chemsrc(i,j,m,itr) = toto(j,i)
                  end do
                end do
              end do
            else
              do m = 1 , 12
                do j = 1 , jx
                  do i = 1 , iy
                    chemsrc(i,j,m,itr) = 0.0D0
                  end do
                end do
              end do
            end if
            if ( aertyp(5:5).eq.'1' ) then
              istatus = nf90_inq_varid(ncid, "so2_monthly", ivarid)
              if (istatus /= nf90_noerr) then
                write (6,*) 'Error so2_monthly variable undefined'
                write (6,*) nf90_strerror(istatus)
                call fatal(__FILE__,__LINE__,'AEROSOL VAR FIND ERROR')
              end if
              istart(1) = 1
              istart(2) = 1
              icount(1) = jx
              icount(2) = iy
              icount(3) = 1
              do m = 1 , 12
                istart(3) = m
                istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
                if (istatus /= nf90_noerr) then
                  write (6,*) 'Error reading so2_monthly variable'
                  write (6,*) nf90_strerror(istatus)
                  call fatal(__FILE__,__LINE__,'AEROSOL VAR READ ERROR')
                end if
                do j = 1 , jx
                  do i = 1 , iy
                    chemsrc(i,j,m,itr) = chemsrc(i,j,m,itr) + toto(j,i)
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
                call fatal(__FILE__,__LINE__,'AEROSOL BC FIND ERROR')
              end if
              istatus = nf90_get_var(ncid, ivarid, toto)
              if (istatus /= nf90_noerr) then
                write (6,*) 'Error reading bc variable'
                write (6,*) nf90_strerror(istatus)
                call fatal(__FILE__,__LINE__,'AEROSOL BC READ ERROR')
              end if
              do m = 1 , 12
                do j = 1 , jx
                  do i = 1 , iy
                    chemsrc(i,j,m,itr) = toto(j,i)
                  end do
                end do
              end do
            else
              do m = 1 , 12
                do j = 1 , jx
                  do i = 1 , iy
                    chemsrc(i,j,m,itr) = 0.0D0
                  end do
                end do
              end do
            end if
            if ( aertyp(5:5).eq.'1' ) then
              istatus = nf90_inq_varid(ncid, "bc_monthly", ivarid)
              if (istatus /= nf90_noerr) then
                write (6,*) 'Error bc_monthly variable undefined'
                write (6,*) nf90_strerror(istatus)
                call fatal(__FILE__,__LINE__,'AEROSOL VAR FIND ERROR')
              end if
              istart(1) = 1
              istart(2) = 1
              icount(1) = jx
              icount(2) = iy
              icount(3) = 1
              do m = 1 , 12
                istart(3) = m
                istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
                if (istatus /= nf90_noerr) then
                  write (6,*) 'Error reading bc_monthly variable'
                  write (6,*) nf90_strerror(istatus)
                  call fatal(__FILE__,__LINE__,'AEROSOL VAR READ ERROR')
                end if
                do j = 1 , jx
                  do i = 1 , iy
                    chemsrc(i,j,m,itr) = chemsrc(i,j,m,itr) + toto(j,i)
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
                call fatal(__FILE__,__LINE__,'AEROSOL OC FIND ERROR')
              end if
              istatus = nf90_get_var(ncid, ivarid, toto)
              if (istatus /= nf90_noerr) then
                write (6,*) 'Error reading oc variable'
                write (6,*) nf90_strerror(istatus)
                call fatal(__FILE__,__LINE__,'AEROSOL OC READ ERROR')
              end if
              do m = 1 , 12
                do j = 1 , jx
                  do i = 1 , iy
                    chemsrc(i,j,m,itr) = toto(j,i)
                  end do
                end do
              end do
            else
              do m = 1 , 12
                do j = 1 , jx
                  do i = 1 , iy
                    chemsrc(i,j,m,itr) = 0.0D0
                  end do
                end do
              end do
            end if
            if ( aertyp(5:5).eq.'1' ) then
              istatus = nf90_inq_varid(ncid, "oc_monthly", ivarid)
              if (istatus /= nf90_noerr) then
                write (6,*) 'Error oc_monthly variable undefined'
                write (6,*) nf90_strerror(istatus)
                call fatal(__FILE__,__LINE__,'AEROSOL VAR FIND ERROR')
              end if
              istart(1) = 1
              istart(2) = 1
              icount(1) = jx
              icount(2) = iy
              icount(3) = 1
              do m = 1 , 12
                istart(3) = m
                istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
                if (istatus /= nf90_noerr) then
                  write (6,*) 'Error reading oc_monthly variable'
                  write (6,*) nf90_strerror(istatus)
                  call fatal(__FILE__,__LINE__,'AEROSOL VAR READ ERROR')
                end if
                do j = 1 , jx
                  do i = 1 , iy
                    chemsrc(i,j,m,itr) = chemsrc(i,j,m,itr) + toto(j,i)
                  end do
                end do
              end do
            end if
          else
          end if
        end if
      end do
      if (aertyp(4:5).ne.'00') then
        istatus = nf90_close(ncid)
        if ( istatus /= nf90_noerr) then
          write (6,*) 'Error Closing Aerosol file ', trim(ffin)
          write (6,*) nf90_strerror(istatus)
          call fatal(__FILE__,__LINE__,'AEROSOL FILE CLOSE ERROR')
        end if
      end if
#endif
 
!     sulfates sources

      do m = 1 , 12
#ifdef MPP1
        do j = 1 , jendl
#else
        do j = 1 , jx
#endif
          do i = 1 , iy
            if ( iso4.gt.0 ) chemsrc(i,j,m,iso4)                        &
               & = 0.02*chemsrc(i,j,m,iso2)
            if ( iso2.gt.0 ) chemsrc(i,j,m,iso2)                        &
               & = 0.98*chemsrc(i,j,m,iso2)
 
!           partition hydrophilic hydrophonic ( cooke et al.1999)
!           BC
            if ( ibchb.gt.0 .and. ibchl.gt.0 ) then
              chemsrc(i,j,m,ibchl) = 0.2*chemsrc(i,j,m,ibchb)
              chemsrc(i,j,m,ibchb) = 0.8*chemsrc(i,j,m,ibchb)
            end if
!           OC
            if ( iochb.gt.0 .and. iochl.gt.0 ) then
              chemsrc(i,j,m,iochl) = 0.5*chemsrc(i,j,m,iochb)
              chemsrc(i,j,m,iochb) = 0.5*chemsrc(i,j,m,iochb)
            end if
          end do
        end do
      end do
 
!     OPtical properties / internal mixing
      if ( mixtype.eq.2 ) then
#ifdef MPP1
        if ( myid.eq.0 ) then
          inquire (file='optdat.bin',exist=there)
          if ( .not.there ) then
            write (*,*) 'For mixtype=2, optdat.bin is required'
            write (*,*) 'ln -s ../Main/Commons/optdat.bin optdat.bin'
            call fatal(__FILE__,__LINE__,'optdat.bin is required')
          end if
          open (iutopt,file='optdat.bin',form='unformatted',            &
              & recl=4*19*11*11*11*11*ibyte,access='direct')
          read (iutopt,rec=1) ((((((dextmix(i,j,k,l,m,n),i=1,4),j=1,19),&
                            & k=1,11),l=1,11),m=1,11),n=1,11)
          read (iutopt,rec=2) ((((((dssamix(i,j,k,l,m,n),i=1,4),j=1,19),&
                            & k=1,11),l=1,11),m=1,11),n=1,11)
          read (iutopt,rec=3) ((((((dgmix(i,j,k,l,m,n),i=1,4),j=1,19),k=&
                            & 1,11),l=1,11),m=1,11),n=1,11)
          close (iutopt)
        end if
        call mpi_bcast(dextmix,4*19*11*11*11*11,mpi_real,0,             &
                     & mpi_comm_world,ierr)
        call mpi_bcast(dssamix,4*19*11*11*11*11,mpi_real,0,             &
                     & mpi_comm_world,ierr)
        call mpi_bcast(dgmix,4*19*11*11*11*11,mpi_real,0,mpi_comm_world,&
                     & ierr)
#else
        inquire (file='optdat.bin',exist=there)
        if ( .not.there ) then
          write (*,*) 'For mixtype=2, optdat.bin is required'
          write (*,*) 'ln -s ../Main/Commons/optdat.bin optdat.bin'
          call fatal(__FILE__,__LINE__,'optdat.bin is required')
        end if
        open (iutopt,file='optdat.bin',form='unformatted',              &
            & recl=4*19*11*11*11*11*ibyte,access='direct')
        read (iutopt,rec=1) ((((((dextmix(i,j,k,l,m,n),i=1,4),j=1,19),k=&
                          & 1,11),l=1,11),m=1,11),n=1,11)
        read (iutopt,rec=2) ((((((dssamix(i,j,k,l,m,n),i=1,4),j=1,19),k=&
                          & 1,11),l=1,11),m=1,11),n=1,11)
        read (iutopt,rec=3) ((((((dgmix(i,j,k,l,m,n),i=1,4),j=1,19),k=1,&
                          & 11),l=1,11),m=1,11),n=1,11)
        close (iutopt)
#endif
 
!       Check !
 
        do k = 1 , 11
          do l = 1 , 11
            do m = 1 , 11
              do n = 1 , 11
 
                if ( k+l+m+n.eq.14 ) then
                  do i = 1 , 4
                    do j = 1 , 19
 
                      if ( (dextmix(i,j,k,l,m,n).lt.0.) .or.            &
                         & (dextmix(i,j,k,l,m,n).gt.20.) ) then
                        write (aline,*) 'problem in dextmix ' ,         &
                                      & dextmix(i,j,k,l,m,n)
                        call say
                        call fatal(__FILE__,__LINE__,'DETMIX ERROR')
                      end if
 
                      if ( (dssamix(i,j,k,l,m,n).lt.0.) .or.            &
                         & (dssamix(i,j,k,l,m,n).gt.1.) ) then
                        write (aline,*) 'problem in dssamix ' ,         &
                                      & dssamix(i,j,k,l,m,n)
                        call say
                        call fatal(__FILE__,__LINE__,'DSSAMIX ERROR')
                      end if
 
                      if ( (dgmix(i,j,k,l,m,n).lt.0.) .or.              &
                         & (dgmix(i,j,k,l,m,n).gt.1.) ) then
                        write (aline,*) 'problem in dgmix ' ,           &
                                      & dgmix(i,j,k,l,m,n)
                        call say
                        call fatal(__FILE__,__LINE__,'DGMIX ERROR')
                      end if
 
                    end do
                  end do
                end if
 
              end do
            end do
          end do
        end do
 
#ifdef MPP1
        if ( myid.eq.0 ) write (*,*) '! OPDATA CHECKED !'
#else
        write (*,*) '! OPDATA CHECKED !'
#endif
      end if
      end subroutine chsrfem
