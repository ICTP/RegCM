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
 
      subroutine chemtap

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine writes the model chem                           c
!                                                                     c
!     iutl : is the output unit number for large-domain variables.    c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use mod_regcm_param
      use mod_param2
      use mod_iunits
      use mod_date
      use mod_main
      use mod_mainchem
      use mod_trachem
      implicit none
!
! Local variables
!
      integer :: i , itr , j , k
      real(8) :: mpd
 
      if ( iotyp.eq.2 ) write (iutchem) idatex
 
!     3 d fields
 
      do itr = 1 , ntr
        do k = kx , 1 , -1
          do i = 1 , ix - 2
#ifdef MPP1
            do j = 1 , mjx - 2
              fchem(j,i) = chia_io(i+1,k,j+1,itr)/psa_io(i+1,j+1)
            end do
#else
            do j = 1 , jx - 2
              fchem(j,i) = chia(i+1,k,j+1,itr)/psa(i+1,j+1)
            end do
#endif
          end do
          if ( iotyp.eq.1 ) then
            nrcchem = nrcchem + 1
            write (iutchem,rec=nrcchem) fchem
          else if ( iotyp.eq.2 ) then
            write (iutchem) fchem
          else
          end if
        end do
      end do
 
!     optical properties of the tracer mixing (independant of ntr)
 
      do k = kx , 1 , -1
        do i = 1 , ix - 2
#ifdef MPP1
          do j = 1 , mjx - 2
            fchem(j,i) = aerext_io(i+1,k,j+1)
          end do
#else
          do j = 1 , jx - 2
            fchem(j,i) = aerext(i+1,k,j+1)
          end do
#endif
        end do
        if ( iotyp.eq.1 ) then
          nrcchem = nrcchem + 1
          write (iutchem,rec=nrcchem) fchem
        else if ( iotyp.eq.2 ) then
          write (iutchem) fchem
        else
        end if
      end do
 
      do k = kx , 1 , -1
        do i = 1 , ix - 2
#ifdef MPP1
          do j = 1 , mjx - 2
            fchem(j,i) = aerssa_io(i+1,k,j+1)
          end do
#else
          do j = 1 , jx - 2
            fchem(j,i) = aerssa(i+1,k,j+1)
          end do
#endif
        end do
        if ( iotyp.eq.1 ) then
          nrcchem = nrcchem + 1
          write (iutchem,rec=nrcchem) fchem
        else if ( iotyp.eq.2 ) then
          write (iutchem) fchem
        else
        end if
      end do
 
      do k = kx , 1 , -1
        do i = 1 , ix - 2
#ifdef MPP1
          do j = 1 , mjx - 2
            fchem(j,i) = aerasp_io(i+1,k,j+1)
          end do
#else
          do j = 1 , jx - 2
            fchem(j,i) = aerasp(i+1,k,j+1)
          end do
#endif
        end do
        if ( iotyp.eq.1 ) then
          nrcchem = nrcchem + 1
          write (iutchem,rec=nrcchem) fchem
        else if ( iotyp.eq.2 ) then
          write (iutchem) fchem
        else
        end if
      end do
 
!     2d fields
!     instantaneous colum burden
      do itr = 1 , ntr
!       -------
        do i = 1 , ix - 2
#ifdef MPP1
          do j = 1 , mjx - 2
            fchem(j,i) = dtrace_io(i+1,j+1,itr)
          end do
#else
          do j = 1 , jx - 2
            fchem(j,i) = dtrace(i+1,j+1,itr)
          end do
#endif
        end do
        if ( iotyp.eq.1 ) then
          nrcchem = nrcchem + 1
          write (iutchem,rec=nrcchem) fchem
        else if ( iotyp.eq.2 ) then
          write (iutchem) fchem
        else
        end if
!       -------
!       averaged deposition rates
 
!       the variable are the cumul between two output time steps.
!       in mg.m-2 ( cf tarcbud). We want to get the correspnding average
!       rate between two time steps: so multiply by mpd to get a
!       deposition /emissiom rate in mg.m-2.day ( consistant with
!       rainfall output)
        mpd = 24./chemfrq
!       CARE here CUMUL
!       mpd = 1.
 
        do i = 1 , ix - 2
#ifdef MPP1
          do j = 1 , mjx - 2
            fchem(j,i) = wdlsc_io(i+1,j+1,itr)*mpd
          end do
#else
          do j = 1 , jx - 2
            fchem(j,i) = wdlsc(i+1,j+1,itr)*mpd
          end do
#endif
        end do
        if ( iotyp.eq.1 ) then
          nrcchem = nrcchem + 1
          write (iutchem,rec=nrcchem) fchem
        else if ( iotyp.eq.2 ) then
          write (iutchem) fchem
        else
        end if
!       --------
        do i = 1 , ix - 2
#ifdef MPP1
          do j = 1 , mjx - 2
            fchem(j,i) = wdcvc_io(i+1,j+1,itr)*mpd
          end do
#else
          do j = 1 , jx - 2
            fchem(j,i) = wdcvc(i+1,j+1,itr)*mpd
          end do
#endif
        end do
        if ( iotyp.eq.1 ) then
          nrcchem = nrcchem + 1
          write (iutchem,rec=nrcchem) fchem
        else if ( iotyp.eq.2 ) then
          write (iutchem) fchem
        else
        end if
!       -------------
        do i = 1 , ix - 2
#ifdef MPP1
          do j = 1 , mjx - 2
            fchem(j,i) = ddsfc_io(i+1,j+1,itr)*mpd
          end do
#else
          do j = 1 , jx - 2
            fchem(j,i) = ddsfc(i+1,j+1,itr)*mpd
          end do
#endif
        end do
        if ( iotyp.eq.1 ) then
          nrcchem = nrcchem + 1
          write (iutchem,rec=nrcchem) fchem
        else if ( iotyp.eq.2 ) then
          write (iutchem) fchem
        else
        end if
!----------------------
        do i = 1 , ix - 2
#ifdef MPP1
          do j = 1 , mjx - 2
            fchem(j,i) = wxsg_io(i+1,j+1,itr)*mpd
          end do
#else
          do j = 1 , jx - 2
            fchem(j,i) = wxsg(i+1,j+1,itr)*mpd
          end do
#endif
        end do
        if ( iotyp.eq.1 ) then
          nrcchem = nrcchem + 1
          write (iutchem,rec=nrcchem) fchem
        else if ( iotyp.eq.2 ) then
          write (iutchem) fchem
        else
        end if
!----------------------
        do i = 1 , ix - 2
#ifdef MPP1
          do j = 1 , mjx - 2
            fchem(j,i) = wxaq_io(i+1,j+1,itr)*mpd
          end do
#else
          do j = 1 , jx - 2
            fchem(j,i) = wxaq(i+1,j+1,itr)*mpd
          end do
#endif
        end do
        if ( iotyp.eq.1 ) then
          nrcchem = nrcchem + 1
          write (iutchem,rec=nrcchem) fchem
        else if ( iotyp.eq.2 ) then
          write (iutchem) fchem
        else
        end if
!----------------------
        do i = 1 , ix - 2
#ifdef MPP1
          do j = 1 , mjx - 2
            fchem(j,i) = cemtrac_io(i+1,j+1,itr)*mpd
          end do
#else
          do j = 1 , jx - 2
            fchem(j,i) = cemtrac(i+1,j+1,itr)*mpd
          end do
#endif
        end do
        if ( iotyp.eq.1 ) then
          nrcchem = nrcchem + 1
          write (iutchem,rec=nrcchem) fchem
        else if ( iotyp.eq.2 ) then
          write (iutchem) fchem
        else
        end if
 
!       IF YOU WANT THE INSTANTANEOUS EMISSION IN THE SOURCE
 
!       do i=1,ix-2
!       do j=1,jx-2
!       fchem(j,i) = chemsrc(i+1,j+1,lmonth,itr) * 1.e9
!       unit = microg/m2.s
!       end do
!       end do
!       if (iotyp.eq.1) then
!       nrcchem=nrcchem+1
!       write(iutchem,rec=nrcchem) fchem
!       else if (iotyp.eq.2) then
!       write(iutchem) fchem
!       end if
!       ------------------
!       reinitialisation for avraged dposition rates ( a modifier)
#ifdef MPP1
        do i = 1 , ix - 2
          do j = 1 , mjx - 2
            do k = 1 , kx
              remlsc_io(i,k,j,itr) = 0.
              remcvc_io(i,k,j,itr) = 0.
              rxsg_io(i,k,j,itr) = 0.
              rxsaq1_io(i,k,j,itr) = 0.
              rxsaq2_io(i,k,j,itr) = 0.
            end do
            cemtr_io(i,j,itr) = 0.
            remdrd_io(i,j,itr) = 0.
            wdlsc_io(i,j,itr) = 0.
            wdcvc_io(i,j,itr) = 0.
            ddsfc_io(i,j,itr) = 0.
            wxsg_io(i,j,itr) = 0.
            wxaq_io(i,j,itr) = 0.
            cemtrac_io(i,j,itr) = 0.
          end do
        end do
      end do
#else
        do i = 1 , ix - 2
          do j = 1 , jx - 2
            do k = 1 , kx
              remlsc(i,k,j,itr) = 0.
              remcvc(i,k,j,itr) = 0.
              rxsg(i,k,j,itr) = 0.
              rxsaq1(i,k,j,itr) = 0.
              rxsaq2(i,k,j,itr) = 0.
            end do
            cemtr(i,j,itr) = 0.
            remdrd(i,j,itr) = 0.
            wdlsc(i,j,itr) = 0.
            wdcvc(i,j,itr) = 0.
            ddsfc(i,j,itr) = 0.
            wxsg(i,j,itr) = 0.
            wxaq(i,j,itr) = 0.
            cemtrac(i,j,itr) = 0.
          end do
        end do
      end do
#endif
 
!----total aerosol TOA radiative forcing ( independant of the number of tracer)
 
      do i = 1 , ix - 2
#ifdef MPP1
        do j = 1 , mjx - 2
          fchem(j,i) = aertarf_io(i+1,j+1)
        end do
#else
        do j = 1 , jx - 2
          fchem(j,i) = aertarf(i+1,j+1)
        end do
#endif
      end do
      if ( iotyp.eq.1 ) then
        nrcchem = nrcchem + 1
        write (iutchem,rec=nrcchem) fchem
      else if ( iotyp.eq.2 ) then
        write (iutchem) fchem
      else
      end if
 
 
      do i = 1 , ix - 2
#ifdef MPP1
        do j = 1 , mjx - 2
          fchem(j,i) = aersrrf_io(i+1,j+1)
        end do
#else
        do j = 1 , jx - 2
          fchem(j,i) = aersrrf(i+1,j+1)
        end do
#endif
      end do
      if ( iotyp.eq.1 ) then
        nrcchem = nrcchem + 1
        write (iutchem,rec=nrcchem) fchem
      else if ( iotyp.eq.2 ) then
        write (iutchem) fchem
      else
      end if
 
!     reset rad diag to 0 ( averaged between output time stepsin aerout)
#ifdef MPP1
      do j = 1 , mjx - 1
        do i = 1 , ix - 1
          aertarf_io(i,j) = 0.
          aersrrf_io(i,j) = 0.
        end do
      end do
#else
      do j = 1 , jx - 1
        do i = 1 , ix - 1
          aertarf(i,j) = 0.
          aersrrf(i,j) = 0.
        end do
      end do
#endif

      print * , 'Chem variables written ' , idatex
 
      end subroutine chemtap
