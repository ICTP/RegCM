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
 
      subroutine radtap
!
      use mod_dynparam
      use mod_param1
      use mod_param2
      use mod_param3 , only : ptop
      use mod_date
      use mod_iunits
      use mod_outrad
      use mod_radbuf
#ifdef MPP1
      use mod_mppio
#else
      use mod_main , only : psa
#endif
      implicit none
!
! Local variables
!
      integer :: i , j , k , n
      real(4) , dimension(jxm2,iym2) :: radpsa
!
!
!     Radiation resolution and I/O parameters
!
!     compute total radiative heating flux for the surface,
!     converting units from cgs to mks
!
      print * , 'Writing rad fields at ktau = ' , ktau , idatex
      do i = 1 , iym2
        do j = 1 , jxm2
#ifdef MPP1
          radpsa(j,i) = (psa_io(i+1,j+1)+ptop)*10.
#else
          radpsa(j,i) = (psa(i+1,j+1)+ptop)*10.
#endif
        end do
      end do


      if ( iotyp.eq.1 ) then
        do n = 1 , nrad3d
          if ( n.ne.1 ) then
            do k = kz , 1 , -1
              nrcrad = nrcrad + 1
#ifdef MPP1
              write (iutrad,rec=nrcrad) ((frad3d_io(j,i,k,n),j=1,jxm2), &
                                      & i=1,iym2)
#else
              write (iutrad,rec=nrcrad) ((frad3d(j,i,k,n),j=1,jxm2),i=1,&
                                      & iym2)
#endif
            end do
          end if
        end do
        do n = 1 , nrad2d
          if ( n.lt.10 .or. n .eq. 22 ) then
                            ! skip everything from alb (10 to 21)
            nrcrad = nrcrad + 1
#ifdef MPP1
            write (iutrad,rec=nrcrad) ((frad2d_io(j,i,n),j=1,jxm2),i=1, &
                                    & iym2)
#else
            write (iutrad,rec=nrcrad) ((frad2d(j,i,n),j=1,jxm2),i=1,iym2&
                                    & )
#endif
          end if
        end do
        nrcrad = nrcrad + 1
        write (iutrad,rec=nrcrad) ((radpsa(j,i),j=1,jxm2),i=1,iym2)
      else if ( iotyp.eq.2 ) then
        write (iutrad) idatex
        do n = 1 , nrad3d
          if ( n.ne.1 ) then
            do k = kz , 1 , -1
#ifdef MPP1
              write (iutrad) ((frad3d_io(j,i,k,n),j=1,jxm2),i=1,iym2)
#else
              write (iutrad) ((frad3d(j,i,k,n),j=1,jxm2),i=1,iym2)
#endif
            end do
          end if
        end do
        do n = 1 , nrad2d
#ifdef MPP1
          if ( n.lt.10 )                                                &
            write (iutrad) ((frad2d_io(j,i,n),j=1,jxm2),i=1,iym2)
#else
          if ( n.lt.10 )                                                &
            write (iutrad) ((frad2d(j,i,n),j=1,jxm2),i=1,iym2)
#endif
                            ! skip everything from alb (10 to 21)
        end do
        write (iutrad) ((radpsa(j,i),j=1,jxm2),i=1,iym2)
      else
      end if
 
      end subroutine radtap
