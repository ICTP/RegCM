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

      program oxidant

      use mod_dynparam
      use mod_date
      use mod_grid
      use mod_ingrid
      use mod_wrtoxd
      use mod_header
      use mod_oxidant

      implicit none
!
!  Local Variables
!
      integer :: idate , idatef , iday , ifile , imon , imonnew ,   &
                & imonold , isize , iyr , nnn , inmber , numfile
      integer :: nnnend , nstart
      integer :: ierr
      character(256) :: namelistfile, prgname
! oxdfile : input global oxidant file
! finame  : output file name e.g TRS__ICBC2003060100      
      character(256) :: oxdfile , finame 
      call header(0)
!
!
!     Read input global namelist
!
      call getarg(0, prgname)
      call getarg(1, namelistfile)
      call initparam(namelistfile, ierr)
      if ( ierr/=0 ) then
        write ( 6, * ) 'Parameter initialization not completed'
        write ( 6, * ) 'Usage : '
        write ( 6, * ) '          ', trim(prgname), ' regcm.in'
        write ( 6, * ) ' '
        write ( 6, * ) 'Check argument and namelist syntax'
        stop
      end if
!      
      call init_grid(iy,jx,kz)
      call init_outoxd(jx,iy,kz) 
      call initdate_icbc
      call finddate_icbc(nstart,globidate1)
      call finddate_icbc(nnnend,globidate2)

      write (*,*) 'NSTART,NNNEND: ' , nstart , nnnend
      write (*,*) 'IDATE1,IDATE2: ' , globidate1 , globidate2

      call commhead

      imonold = 0
      ifile = 501
      do nnn = nstart , nnnend      !time loop
         idate = mdate(nnn)
         iyr = idate/1000000
         imon = idate/10000 - iyr*100
         if ( nnn==nstart .or.                                &
             & (imon/=imonold .and. nnn<nnnend .and. nnn>nstart) ) then
             iday = idate/100 - iyr*10000 - imon*100
         write(*,*)trim(dirglob)
         write (finame,99002) trim(dirglob), pthsep, trim(domname),   &
         '_OXBC', idate
         if ( nnn>nstart ) then
         call getmozart(idate)
         end if
         imonnew = imon + 1
         if ( imon>=12 ) then
         imonnew = 1
         iyr = iyr + 1
         end if
          idatef = iyr*1000000 + imonnew*10000 + 100
          if ( imon==1 .or. imon==3 .or. imon==5 .or. imon==7 .or.      &
             & imon==8 .or. imon==10 .or. imon==12 ) then
            inmber = (32-iday)*4 + 1
          else if ( imon==4 .or. imon==6 .or. imon==9 .or. imon==11 )   &
                  & then
            inmber = (31-iday)*4 + 1
          else
            if ( mod(iyr,4)==0 ) then
              inmber = (30-iday)*4 + 1
            else
              inmber = (29-iday)*4 + 1
            end if
            if ( mod(iyr,100)==0 ) inmber = (29-iday)*4 + 1
            if ( mod(iyr,400)==0 ) inmber = (30-iday)*4 + 1
          end if
          if ( igrads==1 ) call gradsctl(finame,idate,inmber)
          call fexist(finame)
          open (64,file=finame,form='unformatted',status='replace',    &
           & recl=jx*iy*ibyte,access='direct')
          imonold = imon
          ifile = ifile + 1
          noutrec = 0
        end if
          call getmozart(idate)
                  
   
      end do                        !time loop
      call free_outoxd
      call free_grid
99002 format (a,a,a,a,i10)

      end program oxidant
