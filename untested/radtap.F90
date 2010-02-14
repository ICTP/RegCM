!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      subroutine radtap
!
      use mod_regcm_param
      use mod_param1
      use mod_param2
      use mod_date
      use mod_iunits
      use mod_outrad
      use mod_radbuf
#ifdef MPP1
      use mod_mppio
#endif
      implicit none
!
! Local variables
!
      integer :: i , j , k , n
!
!
!     Radiation resolution and I/O parameters
!
!     compute total radiative heating flux for the surface,
!     converting units from cgs to mks
!
      print * , 'Writing rad fields at ktau = ' , ktau , idatex
      if ( iotyp.eq.1 ) then
        do n = 1 , nrad3d
          if ( n.ne.1 ) then
            do k = kx , 1 , -1
              nrcrad = nrcrad + 1
#ifdef MPP1
              write (iutrad,rec=nrcrad) ((frad3d_io(j,i,k,n),j=1,jxm2), &
                                      & i=1,ixm2)
#else
              write (iutrad,rec=nrcrad) ((frad3d(j,i,k,n),j=1,jxm2),i=1,&
                                      & ixm2)
#endif
            end do
          end if
        end do
        do n = 1 , nrad2d
          if ( n.lt.10 ) then
                            ! skip everything from alb (10 to 21)
            nrcrad = nrcrad + 1
#ifdef MPP1
            write (iutrad,rec=nrcrad) ((frad2d_io(j,i,n),j=1,jxm2),i=1, &
                                    & ixm2)
#else
            write (iutrad,rec=nrcrad) ((frad2d(j,i,n),j=1,jxm2),i=1,ixm2&
                                    & )
#endif
          end if
        end do
      else if ( iotyp.eq.2 ) then
        write (iutrad) idatex
        do n = 1 , nrad3d
          if ( n.ne.1 ) then
            do k = kx , 1 , -1
#ifdef MPP1
              write (iutrad) ((frad3d_io(j,i,k,n),j=1,jxm2),i=1,ixm2)
#else
              write (iutrad) ((frad3d(j,i,k,n),j=1,jxm2),i=1,ixm2)
#endif
            end do
          end if
        end do
        do n = 1 , nrad2d
#ifdef MPP1
          if ( n.lt.10 )                                                &
            write (iutrad) ((frad2d_io(j,i,n),j=1,jxm2),i=1,ixm2)
#else
          if ( n.lt.10 )                                                &
            write (iutrad) ((frad2d(j,i,n),j=1,jxm2),i=1,ixm2)
#endif
                            ! skip everything from alb (10 to 21)
        end do
      else
      end if
 
      end subroutine radtap
