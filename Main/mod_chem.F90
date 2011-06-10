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
 
module mod_chem

#ifdef CHEMTEST
  use mod_runparams
  use mod_message
  use mod_mppio
  use mod_memutil

  implicit none

  real(8), pointer, dimension(:,:,:)  :: oh,ho2,o3,no3,h2o2
  real(8), pointer, dimension(:,:,:)  :: oh_io,ho2_io,o3_io,no3_io,h2o2_io
  real(8), pointer, dimension(:,:,:)  :: sav7,sav_7

  integer , private :: iutox = 123
  integer , private :: oxrec


  real(8) , parameter :: convertfact
  contains


  subroutine allocate_mod_chem
    implicit none
    call getmem3d(oh,1,iy,1,kz,1,jxp,'chem:oh')
    call getmem3d(ho2,1,iy,1,kz,1,jxp,'chem:ho2')
    call getmem3d(o3,1,iy,1,kz,1,jxp,'chem:o3')
    call getmem3d(no3,1,iy,1,kz,1,jxp,'chem:no3')
    call getmem3d(h2o2,1,iy,1,kz,1,jxp,'chem:h2o2')
    call getmem3d(sav7,1,iy,1,5*kz,1,jxp,'chem:sav7')
    call getmem3d(sav_7,1,iy,1,5*kz,1,jx,'chem:sav_7')
    call getmem3d(oh_io,1,iy,1,kz,1,jx,'chem:oh_io')
    call getmem3d(ho2_io,1,iy,1,kz,1,jx,'chem:ho2_io')
    call getmem3d(o3_io,1,iy,1,kz,1,jx,'chem:o3_io')
    call getmem3d(no3_io,1,iy,1,kz,1,jx,'chem:no3_io')
    call getmem3d(h2o2_io,1,iy,1,kz,1,jx,'chem:h2o2_io')
  end subroutine allocate_mod_chem
!
! this routine reads the oxidant data written by 
! the oxidant program 
!
  subroutine init_chem
#ifndef IBM
    use mpi
#endif
    implicit none 
#ifdef IBM
    include 'mpif.h'
#endif
    logical :: existing
    character(256) :: finm
    integer :: i , j , k , ierr
    integer :: nxxx , nyyy, kzzz 
    real(4) , dimension(iy,jx) :: io2d
!
    existing = .false.

    if ( myid == 0 ) then
      if (ndate0 == globidate1 .or.  &
         (((ndate0/10000)*100+1)*100 == &
         ((globidate1/10000)*100+1)*100 ) ) then
        write (finm,99001) trim(dirglob),pthsep,trim(domname), &
               '_OXBC',globidate1
      else
        write (finm,99001) trim(dirglob),pthsep,trim(domname), &
               '_OXBC',((ndate0/10000)*100+1)*100
      end if
      inquire (file=finm,exist=existing)
      if (.not.existing) then
        write (aline,*) 'The following OXBC File does not exist:' ,  & 
                        trim(finm), 'please check location'
        call say
        call fatal(__FILE__,__LINE__, 'OXBC FILE NOT FOUND')
      else
        open (iutox,file=finm,form='unformatted',status='old', &
              access='direct',recl=iy*jx*ibyte)
      end if
      oxrec = 0
    end if
    do
      if ( myid == 0 ) then
      oxrec = oxrec + 1
      read (iutox,rec=oxrec) ndate0 , nxxx , nyyy , kzzz
      if ( nyyy /= iy .or. nxxx /= jx .or. kzzz /= kz ) then
        write (aline,*) 'SET IN regcm.param: IY=' , iy , ' JX=' , &
                           jx , ' KX=' , kz
        call say
        write (aline,*) 'SET IN OXBC: NY=' , nyyy , ' NX=' ,      &
                           nxx  , ' NZ=' , kzzz
        call say
        call fatal(__FILE__,__LINE__,                             &
                       'IMPROPER DIMENSION SPECIFICATION')
          end if
          print * , 'READING INITAL CONDITIONS' , ndate0
!!            if ( ndate0 < mdatez(nnnchk) ) then
!!            if ( ndate0 < 0 ) then
!!              print * , ndate0 , mdatez(nnnchk) , nnnchk
!!              print * , 'read in datasets at :' , ndate0
!!              oxrec = oxrec + kz*5
!!              print * , 'Searching for proper date: ' , ndate1 
!!                   mdatez(nnnchk+1)
!!              print * , ndate0 , mdatez(nnnchk)
!!              cycle ! Proper date not found
!!            else if ( ndate0 > mdatez(nnnchk) ) then
!!              write (aline,*) ndate0 , mdatez(nnnchk)
!!              call say
!!              call fatal(__FILE__,__LINE__,                             &
!!                         'DATE IN ICBC FILE EXCEEDED DATE IN RegCM')
!!            end if
        end if
        exit ! Found proper date
    end do
!!!!!!!!!!!!!!!!!!!!!! Start read!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ( myid == 0 ) then
        print * , 'OH ROBERTA'
        do k = kz , 1 , -1
          oxrec = oxrec + 1
          read (iutox,rec=oxrec) ((io2d(i,j),j=1,jx),i=1,iy)
          do j = 1 , jx
            do i = 1 , iy
              oh_io(i,k,j) = dble(io2d(i,j))
            end do
          end do
        end do
        print * , 'HO2'
        do k = kz , 1 , -1
          oxrec = oxrec + 1
          read (iutox,rec=oxrec) ((io2d(i,j),j=1,jx),i=1,iy)
          do j = 1 , jx
            do i = 1 , iy
              ho2_io(i,k,j) = dble(io2d(i,j))
            end do
          end do
        end do
        print * , 'O3'
        do k = kz , 1 , -1
          oxrec = oxrec + 1
          read (iutox,rec=oxrec) ((io2d(i,j),j=1,jx),i=1,iy)
          do j = 1 , jx
            do i = 1 , iy
              o3_io(i,k,j) = dble(io2d(i,j))
            end do
          end do
        end do
        print * , 'NO3'
        do k = kz , 1 , -1
          oxrec = oxrec + 1
          read (iutox,rec=oxrec) ((io2d(i,j),j=1,jx),i=1,iy)
          do j = 1 , jx
            do i = 1 , iy
              no3_io(i,k,j) = dble(io2d(i,j))
            end do
          end do
        end do
        print * , 'H2O2'
        do k = kz , 1 , -1
          oxrec = oxrec + 1
          read (iutox,rec=oxrec) ((io2d(i,j),j=1,jx),i=1,iy)
          do j = 1 , jx
            do i = 1 , iy
              h2o2_io(i,k,j) = dble(io2d(i,j))
            end do
          end do
        end do

        do j = 1 , jx
          do k = 1 , kz
            do i = 1 , iy
              sav_7(i,k,j)      = oh_io(i,k,j)
              sav_7(i,kz+k,j)   = ho2_io(i,k,j)
              sav_7(i,kz*2+k,j) = o3_io(i,k,j)
              sav_7(i,kz*3+k,j) = no3_io(i,k,j)
              sav_7(i,kz*4+k,j) = h2o2_io(i,k,j)
            end do
          end do
        end do  
      end if ! end if myid=0
!
!     Start transmission of data to other processors
!
      call mpi_bcast(ndate1,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_scatter(sav_7(1,1,1),iy*(kz*5)*jxp,mpi_real8,        &
                       sav7(1,1,1), iy*(kz*5)*jxp,mpi_real8,        &
                       0,mpi_comm_world,ierr)

    do j = 1 , jendl
       do k = 1 , kz
          do i = 1 , iy
             oh  (i,k,j) = sav7(i,k,j)  
             ho2 (i,k,j) = sav7(i,kz+k,j)
             o3  (i,k,j) = sav7(i,kz*2+k,j)
             no3 (i,k,j) = sav7(i,kz*3+k,j)
             h2o2(i,k,j) = sav7(i,kz*4+k,j)
          end do
       end do
    end do

99001 format (a,a,a,a,i0.10)
  end subroutine init_chem
!
  subroutine bdyin_chem
#ifndef IBM
    use mpi
#endif
    implicit none 
#ifdef IBM
    include 'mpif.h'
#endif
    logical :: existing
    character(256) :: finm
    integer :: i , j , k , ierr , ierr1
    real(4) , dimension(iy,jx) :: io2d
    real(8) :: dtbdys

    existing = .false.
!
    dtbdys = ibdyfrq*secph
    if ( myid == 0 ) then
!     write(*,*)'SSSSSSSSSSSSSSS',ndate1,mdatez(nnnchk+1)
      do
      oxrec=oxrec+1
      read (iutox,rec=oxrec,iostat=ierr1) ndate1
!       write(*,*)'NNNNNNOXXXXX',ndate1,mdatez(nnnchk)
        if ( ierr1 /= 0 ) then
          close (iutox)
          iutox = iutox + 1
          write (finm,99001) trim(dirglob),pthsep,trim(domname),    &
                             '_OXBC',((ndate1/10000)*100+1)*100
          inquire(file=finm,exist=existing)
          if ( .not.existing ) then
            write (aline,*)                                      &
                 'The following OX IBC File does not exist: ' ,  &
                 trim(finm), 'please check location'
            call say
            call fatal(__FILE__,__LINE__,aline) 
          else 
            open (iutox,file=finm,form='unformatted',status='old', &
                  access='direct',recl=iy*jx*ibyte)
          end if
          oxrec = 0
          print * , 'CHANGING OX BDY UNIT NUMBER:  iutbc=' , iutox
          if ( iutox > 999 ) then
            call fatal(__FILE__,__LINE__,'BDY UNIT MAX EXCEEDED')
          end if
          cycle
        end if
!       if ( ndate1 < mdatez(nnnchk) ) then
!         if ( ndate1 < mdatez(nnnchk) ) then
!           print * , 'Searching for proper date: ' , ndate1 ,      &
!                      mdatez(nnnchk+1)
!           print * , 'read in datasets at :' , ndate0
!           oxrec = oxrec + kz*5
!           cycle
!         end if
!       else if ( ndate1 > mdatez(nnnchk+1) ) then
!         print * , 'DATE IN OX BC FILE EXCEEDED DATE IN RegCM'
!         print * , ndate1 , mdatez(nnnchk+1) , nnnchk + 1
!         call fatal(__FILE__,__LINE__,'ICBC date')
!       end if
        exit
      end do
    end if
!!!!!!!!!!!!!!!!!!!!!! Start read!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ( myid == 0 ) then
      print * , 'OH'
      do k = kz , 1 , -1
        oxrec = oxrec + 1
        read (iutox,rec=oxrec) ((io2d(i,j),j=1,jx),i=1,iy)
        do j = 1 , jx
          do i = 1 , iy
            oh_io(i,k,j) = dble(io2d(i,j))
          end do
        end do
      end do
!     write(*,*)'OH====',oh_io(17,1,32)
      write(*,*)'MAX OH',maxval(oh_io(:,kz,:))
      print * , 'HO2'
      do k = kz , 1 , -1
        oxrec = oxrec + 1
        read (iutox,rec=oxrec) ((io2d(i,j),j=1,jx),i=1,iy)
        do j = 1 , jx
          do i = 1 , iy
            ho2_io(i,k,j) = dble(io2d(i,j))
          end do
        end do
      end do
      print * , 'O3'
      do k = kz , 1 , -1
        oxrec = oxrec + 1
        read (iutox,rec=oxrec) ((io2d(i,j),j=1,jx),i=1,iy)
        do j = 1 , jx
          do i = 1 , iy
            o3_io(i,k,j) = dble(io2d(i,j))
          end do
        end do
      end do
!     write(*,*)'O3====',o3_io(17,kz,32)*convertfac
      write(*,*)'MAX O3',maxval(o3_io(:,kz,:))*convertfac
      print * , 'NO3'
      do k = kz , 1 , -1
        oxrec = oxrec + 1
        read (iutox,rec=oxrec) ((io2d(i,j),j=1,jx),i=1,iy)
        do j = 1 , jx
          do i = 1 , iy
            no3_io(i,k,j) = dble(io2d(i,j))
          end do
        end do
      end do
      print * , 'H2O2'
      do k = kz , 1 , -1
        oxrec = oxrec + 1
        read (iutox,rec=oxrec) ((io2d(i,j),j=1,jx),i=1,iy)
        do j = 1 , jx
          do i = 1 , iy
            h2o2_io(i,k,j) = dble(io2d(i,j))
          end do
        end do
      end do
      write(*,*)'RECORD',oxrec
!     write(*,*)'H2O2==',h2o2_io(17,kz,32)*convertfac
      write(*,*)'MAX H2O2',maxval(h2o2_io(:,kz,:))*convertfac

      do j = 1 , jx
        do k = 1 , kz
          do i = 1 , iy
            sav_7(i,k,j)      = oh_io(i,k,j)
            sav_7(i,kz+k,j)   = ho2_io(i,k,j)
            sav_7(i,kz*2+k,j) = o3_io(i,k,j)
            sav_7(i,kz*3+k,j) = no3_io(i,k,j)
            sav_7(i,kz*4+k,j) = h2o2_io(i,k,j)
          end do
        end do
      end do  
    end if ! end if myid=0
!
!   Start transmission of data to other processors
!
    call mpi_scatter(sav_7(1,1,1),iy*(kz*5)*jxp,mpi_real8,        &
                     sav7(1,1,1), iy*(kz*5)*jxp,mpi_real8,        &
                     0,mpi_comm_world,ierr)

    do j = 1 , jendl
      do k = 1 , kz
        do i = 1 , iy
          oh  (i,k,j) = sav7(i,k,j)  
          ho2 (i,k,j) = sav7(i,kz+k,j)
          o3  (i,k,j) = sav7(i,kz*2+k,j)
          no3 (i,k,j) = sav7(i,kz*3+k,j)
          h2o2(i,k,j) = sav7(i,kz*4+k,j)
        end do
      end do
    end do
    write(*,*)'MAX OH 2',maxval(oh(:,kz,:))
    write(*,*)'MAX O3 2',maxval(o3(:,kz,:))*convertfac
    write(*,*)'MAX H2O2 2',maxval(h2o2(:,kz,:))*convertfac

99001 format (a,a,a,a,i0.10)

  end subroutine bdyin_chem

#endif

end module mod_chem
