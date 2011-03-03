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
 
      module mod_nudge
!
! Relaxation and Sponge Boundary Conditions
!
      use mod_constants
      use mod_dynparam
      use mod_runparams
      use mod_bdycod
      use mod_main
      use mod_service
!
      private
!
      public :: sponge_p , sponge_t , spongeqv , sponge_u , sponge_v
      public :: nudge_p , nudge_t , nudgeqv , nudge_u , nudge_v
!
      contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     these subroutines apply relaxation boundary conditions to the   c
!     tendency term - xpten.                                          c
!                                                                     c
!     ip    : is the number of slices affected by nudging.            c
!                                                                     c
!     xt    : is the time in minutes for variable "psb".              c
!                                                                     c
!     fcoef : are the coefficients for the newtonian term.            c
!                                                                     c
!     gcoef : are the coefficients for the diffusion term.            c
!                                                                     c
!     xpten : is the tendency calculated from the model.              c
!                                                                     c
!     peb, pwb, pss, pnb : are the observed boundary values           c
!                   on east, west, south, and north boundaries.       c
!                                                                     c
!     pebt, pwbt, psbt, pnbt : are the large-scale or observed        c
!             tendencies at east, west, south, and north boundaries.  c
!                                                                     c
!     psb    : is the variable at tau-1.                               c
!                                                                     c
!     ie = iy, je = jx for dot-point variables.                       c
!     ie = iym1, je = jxm1 for cross-point variables.                 c
!                                                                     c
!     j    : is the j'th slice of the tendency to be adjusted.        c
!     ibdy : type of boundary condition relaxation, 1=linear        c
!              5 = exponential                                        c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine nudge_p(ip,fcoef,gcoef,xt,xpten,j,ibdy)
!
      implicit none
!
      real(8) :: fcoef , gcoef , xt
      integer :: ibdy , ip , j
      real(8) , dimension(iy) :: xpten
      intent (in) fcoef , gcoef , ibdy , ip , j , xt
      intent (inout) xpten
!
      real(8) :: dtb , fcx , fls0 , fls1 , fls2 , fls3 , fls4 , gcx
      integer :: i , ii
#ifndef BAND
      integer :: ibeg , iend , jj , jsls
#else
#ifndef MPP1
      integer :: jp1 , jm1
#endif
#endif
!
#ifdef MPP1
#ifndef BAND
      integer :: jwb , jeb
#endif
#endif
      character (len=50) :: subroutine_name='nudge_p'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
#ifdef BAND
!----------------------------------------------------------------------
!
      dtb = xt*minph
!
!-----determine which relaxation method to use:linear/expon.
!
      if ( ibdy == 1 ) then
!
!---------use linear method
!
#ifdef MPP1
!------interior j slices:
         do i = 2 , ip
            ii = iym1 - i + 1
            fcx = fcoef*xfun(i)
            gcx = gcoef*xfun(i)
!.......south boundary:
            fls0 = (pss(i,j)+dtb*psbt(i,j)) - sps2%ps(i,j)
            fls1 = (pss(i,j-1)+dtb*psbt(i,j-1)) - sps2%ps(i,j-1)
            fls2 = (pss(i,j+1)+dtb*psbt(i,j+1)) - sps2%ps(i,j+1)
            fls3 = (pss(i-1,j)+dtb*psbt(i-1,j)) - sps2%ps(i-1,j)
            fls4 = (pss(i+1,j)+dtb*psbt(i+1,j)) - sps2%ps(i+1,j)
            xpten(i) = xpten(i) + fcx*fls0 -                            &
                    & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
            fls0 = (pnb(i,j)+dtb*pnbt(i,j)) - sps2%ps(ii,j)
            fls1 = (pnb(i,j-1)+dtb*pnbt(i,j-1)) - sps2%ps(ii,j-1)
            fls2 = (pnb(i,j+1)+dtb*pnbt(i,j+1)) - sps2%ps(ii,j+1)
            fls3 = (pnb(i-1,j)+dtb*pnbt(i-1,j)) - sps2%ps(ii-1,j)
            fls4 = (pnb(i+1,j)+dtb*pnbt(i+1,j)) - sps2%ps(ii+1,j)
            xpten(ii) = xpten(ii) + fcx*fls0 -                          &
                     & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
         end do
#else
         jm1 = j-1
         jp1 = j+1
         if(jm1 == 0) jm1 = jx
         if(jp1 == jx+1) jp1 = 1
!------interior j slices:
         do i = 2 , ip
            ii = iym1 - i + 1
            fcx = fcoef*xfun(i)
            gcx = gcoef*xfun(i)
!.......south boundary:
            fls0 = (pss(i,j)+dtb*psbt(i,j)) - sps2%ps(i,j)
            fls1 = (pss(i,jm1)+dtb*psbt(i,jm1)) - sps2%ps(i,jm1)
            fls2 = (pss(i,jp1)+dtb*psbt(i,jp1)) - sps2%ps(i,jp1)
            fls3 = (pss(i-1,j)+dtb*psbt(i-1,j)) - sps2%ps(i-1,j)
            fls4 = (pss(i+1,j)+dtb*psbt(i+1,j)) - sps2%ps(i+1,j)
            xpten(i) = xpten(i) + fcx*fls0 -                            &
                    & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
            fls0 = (pnb(i,j)+dtb*pnbt(i,j)) - sps2%ps(ii,j)
            fls1 = (pnb(i,jm1)+dtb*pnbt(i,jm1)) - sps2%ps(ii,jm1)
            fls2 = (pnb(i,jp1)+dtb*pnbt(i,jp1)) - sps2%ps(ii,jp1)
            fls3 = (pnb(i-1,j)+dtb*pnbt(i-1,j)) - sps2%ps(ii-1,j)
            fls4 = (pnb(i+1,j)+dtb*pnbt(i+1,j)) - sps2%ps(ii+1,j)
            xpten(ii) = xpten(ii) + fcx*fls0 -                          &
                     & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
         end do
#endif
!
      else if ( ibdy == 5 ) then
 
!----------use exponential method
 
#ifdef MPP1
!------interior j slices:
         do i = 2 , ip
            ii = iym1 - i + 1
            fcx = fcoef*xfune(i,kz)
            gcx = gcoef*xfune(i,kz)
!........south boundary:
            fls0 = (pss(i,j)+dtb*psbt(i,j)) - sps2%ps(i,j)
            fls1 = (pss(i,j-1)+dtb*psbt(i,j-1)) - sps2%ps(i,j-1)
            fls2 = (pss(i,j+1)+dtb*psbt(i,j+1)) - sps2%ps(i,j+1)
            fls3 = (pss(i-1,j)+dtb*psbt(i-1,j)) - sps2%ps(i-1,j)
            fls4 = (pss(i+1,j)+dtb*psbt(i+1,j)) - sps2%ps(i+1,j)
            xpten(i) = xpten(i) + fcx*fls0 -                            &
                    & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
            fls0 = (pnb(i,j)+dtb*pnbt(i,j)) - sps2%ps(ii,j)
            fls1 = (pnb(i,j-1)+dtb*pnbt(i,j-1)) - sps2%ps(ii,j-1)
            fls2 = (pnb(i,j+1)+dtb*pnbt(i,j+1)) - sps2%ps(ii,j+1)
            fls3 = (pnb(i-1,j)+dtb*pnbt(i-1,j)) - sps2%ps(ii-1,j)
            fls4 = (pnb(i+1,j)+dtb*pnbt(i+1,j)) - sps2%ps(ii+1,j)
            xpten(ii) = xpten(ii) + fcx*fls0 -                          &
                     & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
         end do
#else
         jm1 = j-1
         jp1 = j+1
         if(jm1 == 0) jm1 = jx
         if(jp1 == jx+1) jp1 = 1
!------interior j slices:
         do i = 2 , ip
            ii = iym1 - i + 1
            fcx = fcoef*xfune(i,kz)
            gcx = gcoef*xfune(i,kz)
!........south boundary:
            fls0 = (pss(i,j)+dtb*psbt(i,j)) - sps2%ps(i,j)
            fls1 = (pss(i,jm1)+dtb*psbt(i,jm1)) - sps2%ps(i,jm1)
            fls2 = (pss(i,jp1)+dtb*psbt(i,jp1)) - sps2%ps(i,jp1)
            fls3 = (pss(i-1,j)+dtb*psbt(i-1,j)) - sps2%ps(i-1,j)
            fls4 = (pss(i+1,j)+dtb*psbt(i+1,j)) - sps2%ps(i+1,j)
            xpten(i) = xpten(i) + fcx*fls0 -                            &
                    & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
            fls0 = (pnb(i,j)+dtb*pnbt(i,j)) - sps2%ps(ii,j)
            fls1 = (pnb(i,jm1)+dtb*pnbt(i,jm1)) - sps2%ps(ii,jm1)
            fls2 = (pnb(i,jp1)+dtb*pnbt(i,jp1)) - sps2%ps(ii,jp1)
            fls3 = (pnb(i-1,j)+dtb*pnbt(i-1,j)) - sps2%ps(ii-1,j)
            fls4 = (pnb(i+1,j)+dtb*pnbt(i+1,j)) - sps2%ps(ii+1,j)
            xpten(ii) = xpten(ii) + fcx*fls0 -                          &
                     & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
         end do
#endif
      end if
#else
!----------------------------------------------------------------------
!
      dtb = xt*minph
#ifdef MPP1
      jsls = j + myid*jxp
      jj = jx - jsls
      if ( jj <= ip ) jsls = jj
      jwb = jsls
      if ( jwb > jxp ) jwb = mod(jwb,jxp)
      if ( jwb == 0 ) jwb = jxp
      if ( myid == nproc-1 ) then
        jeb = jsls
      else
        jeb = jsls + 1
      end if
      if ( jeb > jxp ) jeb = mod(jeb,jxp)
      if ( jeb == 0 ) jeb = jxp
#else
      jsls = j
      jj = jx - jsls
      if ( jj <= ip ) jsls = jj
#endif
!
!-----determine which relaxation method to use:linear/expon.
!
      if ( ibdy == 1 ) then
!
!---------use linear method
!
        if ( jsls > ip ) then
!------interior j slices:
          do i = 2 , ip
            ii = iym1 - i + 1
            fcx = fcoef*xfun(i)
            gcx = gcoef*xfun(i)
!.......south boundary:
            fls0 = (pss(i,j)+dtb*psbt(i,j)) - sps2%ps(i,j)
            fls1 = (pss(i,j-1)+dtb*psbt(i,j-1)) - sps2%ps(i,j-1)
            fls2 = (pss(i,j+1)+dtb*psbt(i,j+1)) - sps2%ps(i,j+1)
            fls3 = (pss(i-1,j)+dtb*psbt(i-1,j)) - sps2%ps(i-1,j)
            fls4 = (pss(i+1,j)+dtb*psbt(i+1,j)) - sps2%ps(i+1,j)
            xpten(i) = xpten(i) + fcx*fls0 -                            &
                    & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
            fls0 = (pnb(i,j)+dtb*pnbt(i,j)) - sps2%ps(ii,j)
            fls1 = (pnb(i,j-1)+dtb*pnbt(i,j-1)) - sps2%ps(ii,j-1)
            fls2 = (pnb(i,j+1)+dtb*pnbt(i,j+1)) - sps2%ps(ii,j+1)
            fls3 = (pnb(i-1,j)+dtb*pnbt(i-1,j)) - sps2%ps(ii-1,j)
            fls4 = (pnb(i+1,j)+dtb*pnbt(i+1,j)) - sps2%ps(ii+1,j)
            xpten(ii) = xpten(ii) + fcx*fls0 -                          &
                     & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
!
        else if ( jsls <= ip ) then
!------east or west boundary slices:
          ibeg = 2
          iend = iym1 - 1
          if ( jsls > 2 ) then
            do i = 2 , jsls - 1
              ii = iym1 - i + 1
              fcx = fcoef*xfun(i)
              gcx = gcoef*xfun(i)
!........south boundary:
              fls0 = (pss(i,j)+dtb*psbt(i,j)) - sps2%ps(i,j)
              fls1 = (pss(i,j-1)+dtb*psbt(i,j-1)) - sps2%ps(i,j-1)
              fls2 = (pss(i,j+1)+dtb*psbt(i,j+1)) - sps2%ps(i,j+1)
              fls3 = (pss(i-1,j)+dtb*psbt(i-1,j)) - sps2%ps(i-1,j)
              fls4 = (pss(i+1,j)+dtb*psbt(i+1,j)) - sps2%ps(i+1,j)
              xpten(i) = xpten(i) + fcx*fls0 -                          &
                      & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!.........north boundary:
              fls0 = (pnb(i,j)+dtb*pnbt(i,j)) - sps2%ps(ii,j)
              fls1 = (pnb(i,j-1)+dtb*pnbt(i,j-1)) - sps2%ps(ii,j-1)
              fls2 = (pnb(i,j+1)+dtb*pnbt(i,j+1)) - sps2%ps(ii,j+1)
              fls3 = (pnb(i-1,j)+dtb*pnbt(i-1,j)) - sps2%ps(ii-1,j)
              fls4 = (pnb(i+1,j)+dtb*pnbt(i+1,j)) - sps2%ps(ii+1,j)
              xpten(ii) = xpten(ii) + fcx*fls0 -                        &
                       & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
            ibeg = jsls
            iend = iym1 - jsls + 1
          end if
!
          if ( jj > ip ) then
!-------west-boundary slice:
            fcx = fcoef*xfun(jsls)
            gcx = gcoef*xfun(jsls)
            do i = ibeg , iend
#ifdef MPP1
              fls0 = (pwb(i,jwb)+dtb*pwbt(i,jwb)) - sps2%ps(i,j)
              fls1 = (pwb(i-1,jwb)+dtb*pwbt(i-1,jwb)) - sps2%ps(i-1,j)
              fls2 = (pwb(i+1,jwb)+dtb*pwbt(i+1,jwb)) - sps2%ps(i+1,j)
              fls3 = (pwb(i,jwb-1)+dtb*pwbt(i,jwb-1)) - sps2%ps(i,j-1)
              fls4 = (pwb(i,jwb+1)+dtb*pwbt(i,jwb+1)) - sps2%ps(i,j+1)
#else
              fls0 = (pwb(i,jsls)+dtb*pwbt(i,jsls)) - sps2%ps(i,j)
              fls1 = (pwb(i-1,jsls)+dtb*pwbt(i-1,jsls)) - sps2%ps(i-1,j)
              fls2 = (pwb(i+1,jsls)+dtb*pwbt(i+1,jsls)) - sps2%ps(i+1,j)
              fls3 = (pwb(i,jsls-1)+dtb*pwbt(i,jsls-1)) - sps2%ps(i,j-1)
              fls4 = (pwb(i,jsls+1)+dtb*pwbt(i,jsls+1)) - sps2%ps(i,j+1)
#endif
              xpten(i) = xpten(i) + fcx*fls0 -                          &
                      & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          else if ( jj <= ip ) then
!-------east-boundary slice:
            fcx = fcoef*xfun(jsls)
            gcx = gcoef*xfun(jsls)
            do i = ibeg , iend
#ifdef MPP1
              fls0 = (peb(i,jeb)+dtb*pebt(i,jeb)) - sps2%ps(i,j)
              fls1 = (peb(i-1,jeb)+dtb*pebt(i-1,jeb)) - sps2%ps(i-1,j)
              fls2 = (peb(i+1,jeb)+dtb*pebt(i+1,jeb)) - sps2%ps(i+1,j)
              fls3 = (peb(i,jeb-1)+dtb*pebt(i,jeb-1)) - sps2%ps(i,j-1)
              fls4 = (peb(i,jeb+1)+dtb*pebt(i,jeb+1)) - sps2%ps(i,j+1)
#else
              fls0 = (peb(i,jsls)+dtb*pebt(i,jsls)) - sps2%ps(i,j)
              fls1 = (peb(i-1,jsls)+dtb*pebt(i-1,jsls)) - sps2%ps(i-1,j)
              fls2 = (peb(i+1,jsls)+dtb*pebt(i+1,jsls)) - sps2%ps(i+1,j)
              fls3 = (peb(i,jsls-1)+dtb*pebt(i,jsls-1)) - sps2%ps(i,j-1)
              fls4 = (peb(i,jsls+1)+dtb*pebt(i,jsls+1)) - sps2%ps(i,j+1)
#endif
              xpten(i) = xpten(i) + fcx*fls0 -                          &
                      & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end if
        end if
!
      else if ( ibdy == 5 ) then
 
!----------use exponential method
 
        if ( jsls > ip ) then
!------interior j slices:
          do i = 2 , ip
            ii = iym1 - i + 1
            fcx = fcoef*xfune(i,kz)
            gcx = gcoef*xfune(i,kz)
!........south boundary:
            fls0 = (pss(i,j)+dtb*psbt(i,j)) - sps2%ps(i,j)
            fls1 = (pss(i,j-1)+dtb*psbt(i,j-1)) - sps2%ps(i,j-1)
            fls2 = (pss(i,j+1)+dtb*psbt(i,j+1)) - sps2%ps(i,j+1)
            fls3 = (pss(i-1,j)+dtb*psbt(i-1,j)) - sps2%ps(i-1,j)
            fls4 = (pss(i+1,j)+dtb*psbt(i+1,j)) - sps2%ps(i+1,j)
            xpten(i) = xpten(i) + fcx*fls0 -                            &
                    & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
            fls0 = (pnb(i,j)+dtb*pnbt(i,j)) - sps2%ps(ii,j)
            fls1 = (pnb(i,j-1)+dtb*pnbt(i,j-1)) - sps2%ps(ii,j-1)
            fls2 = (pnb(i,j+1)+dtb*pnbt(i,j+1)) - sps2%ps(ii,j+1)
            fls3 = (pnb(i-1,j)+dtb*pnbt(i-1,j)) - sps2%ps(ii-1,j)
            fls4 = (pnb(i+1,j)+dtb*pnbt(i+1,j)) - sps2%ps(ii+1,j)
            xpten(ii) = xpten(ii) + fcx*fls0 -                          &
                     & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
!
        else if ( jsls <= ip ) then
!------east or west boundary slices:
          ibeg = 2
          iend = iym1 - 1
          if ( jsls > 2 ) then
            do i = 2 , jsls - 1
              ii = iym1 - i + 1
              fcx = fcoef*xfune(i,kz)
              gcx = gcoef*xfune(i,kz)
!.........south boundary:
              fls0 = (pss(i,j)+dtb*psbt(i,j)) - sps2%ps(i,j)
              fls1 = (pss(i,j-1)+dtb*psbt(i,j-1)) - sps2%ps(i,j-1)
              fls2 = (pss(i,j+1)+dtb*psbt(i,j+1)) - sps2%ps(i,j+1)
              fls3 = (pss(i-1,j)+dtb*psbt(i-1,j)) - sps2%ps(i-1,j)
              fls4 = (pss(i+1,j)+dtb*psbt(i+1,j)) - sps2%ps(i+1,j)
              xpten(i) = xpten(i) + fcx*fls0 -                          &
                      & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!.........north boundary:
              fls0 = (pnb(i,j)+dtb*pnbt(i,j)) - sps2%ps(ii,j)
              fls1 = (pnb(i,j-1)+dtb*pnbt(i,j-1)) - sps2%ps(ii,j-1)
              fls2 = (pnb(i,j+1)+dtb*pnbt(i,j+1)) - sps2%ps(ii,j+1)
              fls3 = (pnb(i-1,j)+dtb*pnbt(i-1,j)) - sps2%ps(ii-1,j)
              fls4 = (pnb(i+1,j)+dtb*pnbt(i+1,j)) - sps2%ps(ii+1,j)
              xpten(ii) = xpten(ii) + fcx*fls0 -                        &
                       & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
            ibeg = jsls
            iend = iym1 - jsls + 1
          end if
!
          if ( jj > ip ) then
!-------west-boundary slice:
            fcx = fcoef*xfune(jsls,kz)
            gcx = gcoef*xfune(jsls,kz)
            do i = ibeg , iend
#ifdef MPP1
              fls0 = (pwb(i,jwb)+dtb*pwbt(i,jwb)) - sps2%ps(i,j)
              fls1 = (pwb(i-1,jwb)+dtb*pwbt(i-1,jwb)) - sps2%ps(i-1,j)
              fls2 = (pwb(i+1,jwb)+dtb*pwbt(i+1,jwb)) - sps2%ps(i+1,j)
              fls3 = (pwb(i,jwb-1)+dtb*pwbt(i,jwb-1)) - sps2%ps(i,j-1)
              fls4 = (pwb(i,jwb+1)+dtb*pwbt(i,jwb+1)) - sps2%ps(i,j+1)
#else
              fls0 = (pwb(i,jsls)+dtb*pwbt(i,jsls)) - sps2%ps(i,j)
              fls1 = (pwb(i-1,jsls)+dtb*pwbt(i-1,jsls)) - sps2%ps(i-1,j)
              fls2 = (pwb(i+1,jsls)+dtb*pwbt(i+1,jsls)) - sps2%ps(i+1,j)
              fls3 = (pwb(i,jsls-1)+dtb*pwbt(i,jsls-1)) - sps2%ps(i,j-1)
              fls4 = (pwb(i,jsls+1)+dtb*pwbt(i,jsls+1)) - sps2%ps(i,j+1)
#endif
              xpten(i) = xpten(i) + fcx*fls0 -                          &
                      & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          else if ( jj <= ip ) then
!-------east-boundary slice:
            fcx = fcoef*xfune(jsls,kz)
            gcx = gcoef*xfune(jsls,kz)
            do i = ibeg , iend
#ifdef MPP1
              fls0 = (peb(i,jeb)+dtb*pebt(i,jeb)) - sps2%ps(i,j)
              fls1 = (peb(i-1,jeb)+dtb*pebt(i-1,jeb)) - sps2%ps(i-1,j)
              fls2 = (peb(i+1,jeb)+dtb*pebt(i+1,jeb)) - sps2%ps(i+1,j)
              fls3 = (peb(i,jeb-1)+dtb*pebt(i,jeb-1)) - sps2%ps(i,j-1)
              fls4 = (peb(i,jeb+1)+dtb*pebt(i,jeb+1)) - sps2%ps(i,j+1)
#else
              fls0 = (peb(i,jsls)+dtb*pebt(i,jsls)) - sps2%ps(i,j)
              fls1 = (peb(i-1,jsls)+dtb*pebt(i-1,jsls)) - sps2%ps(i-1,j)
              fls2 = (peb(i+1,jsls)+dtb*pebt(i+1,jsls)) - sps2%ps(i+1,j)
              fls3 = (peb(i,jsls-1)+dtb*pebt(i,jsls-1)) - sps2%ps(i,j-1)
              fls4 = (peb(i,jsls+1)+dtb*pebt(i,jsls+1)) - sps2%ps(i,j+1)
#endif
              xpten(i) = xpten(i) + fcx*fls0 -                          &
                      & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end if
        end if
!
      end if
#endif
      call time_end(subroutine_name,idindx)
      end subroutine nudge_p
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine nudge_t(ip,fcoef,gcoef,xt,ften,j,ibdy)
!
      implicit none
!
      real(8) :: fcoef , gcoef , xt
      integer :: ibdy , ip , j
      real(8) , dimension(iy,kz) :: ften
      intent (in) fcoef , gcoef , ibdy , ip , j , xt
      intent (inout) ften
!
      real(8) :: dtb , fcx , fls0 , fls1 , fls2 , fls3 , fls4 , gcx
      integer :: i , ii , k
#ifndef BAND
      integer :: ibeg , iend , jj , jsls
#else
#ifndef MPP1
      integer :: jp1 , jm1
#endif
#endif
#ifdef MPP1
#ifndef BAND
      integer :: jwb , jeb
#endif
#endif
      character (len=50) :: subroutine_name='nudge_t'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
#ifdef BAND
!----------------------------------------------------------------------
!
      dtb = xt*minph
!
!-----determine which relaxation method to use:linear/expon.
!
      if ( ibdy == 1 ) then
!
!---------use linear method
!
#ifdef MPP1
!------interior j slices:
         do i = 2 , ip
            ii = iym1 - i + 1
            fcx = fcoef*xfun(i)
            gcx = gcoef*xfun(i)
            do k = 1 , kz
!.......south boundary:
              fls0 = (tsb(i,k,j)+dtb*tsbt(i,k,j)) - atm2%t(i,k,j)
              fls1 = (tsb(i,k,j-1)+dtb*tsbt(i,k,j-1)) - atm2%t(i,k,j-1)
              fls2 = (tsb(i,k,j+1)+dtb*tsbt(i,k,j+1)) - atm2%t(i,k,j+1)
              fls3 = (tsb(i-1,k,j)+dtb*tsbt(i-1,k,j)) - atm2%t(i-1,k,j)
              fls4 = (tsb(i+1,k,j)+dtb*tsbt(i+1,k,j)) - atm2%t(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (tnb(i,k,j)+dtb*tnbt(i,k,j)) - atm2%t(ii,k,j)
              fls1 = (tnb(i,k,j-1)+dtb*tnbt(i,k,j-1)) - atm2%t(ii,k,j-1)
              fls2 = (tnb(i,k,j+1)+dtb*tnbt(i,k,j+1)) - atm2%t(ii,k,j+1)
              fls3 = (tnb(i-1,k,j)+dtb*tnbt(i-1,k,j)) - atm2%t(ii-1,k,j)
              fls4 = (tnb(i+1,k,j)+dtb*tnbt(i+1,k,j)) - atm2%t(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
         end do
#else
         jm1 = j-1
         jp1 = j+1
         if(jm1 == 0) jm1 = jx
         if(jp1 == jx+1) jp1 = 1
!------interior j slices:
         do i = 2 , ip
            ii = iym1 - i + 1
            fcx = fcoef*xfun(i)
            gcx = gcoef*xfun(i)
            do k = 1 , kz
!.......south boundary:
              fls0 = (tsb(i,k,j)+dtb*tsbt(i,k,j)) - atm2%t(i,k,j)
              fls1 = (tsb(i,k,jm1)+dtb*tsbt(i,k,jm1)) - atm2%t(i,k,jm1)
              fls2 = (tsb(i,k,jp1)+dtb*tsbt(i,k,jp1)) - atm2%t(i,k,jp1)
              fls3 = (tsb(i-1,k,j)+dtb*tsbt(i-1,k,j)) - atm2%t(i-1,k,j)
              fls4 = (tsb(i+1,k,j)+dtb*tsbt(i+1,k,j)) - atm2%t(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (tnb(i,k,j)+dtb*tnbt(i,k,j)) - atm2%t(ii,k,j)
              fls1 = (tnb(i,k,jm1)+dtb*tnbt(i,k,jm1)) - atm2%t(ii,k,jm1)
              fls2 = (tnb(i,k,jp1)+dtb*tnbt(i,k,jp1)) - atm2%t(ii,k,jp1)
              fls3 = (tnb(i-1,k,j)+dtb*tnbt(i-1,k,j)) - atm2%t(ii-1,k,j)
              fls4 = (tnb(i+1,k,j)+dtb*tnbt(i+1,k,j)) - atm2%t(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
#endif
      else if ( ibdy == 5 ) then
 
!----------use exponential method
 
#ifdef MPP1
!------interior j slices:
         do i = 2 , ip
            ii = iym1 - i + 1
            do k = 1 , kz
              fcx = fcoef*xfune(i,k)
              gcx = gcoef*xfune(i,k)
!........south boundary:
              fls0 = (tsb(i,k,j)+dtb*tsbt(i,k,j)) - atm2%t(i,k,j)
              fls1 = (tsb(i,k,j-1)+dtb*tsbt(i,k,j-1)) - atm2%t(i,k,j-1)
              fls2 = (tsb(i,k,j+1)+dtb*tsbt(i,k,j+1)) - atm2%t(i,k,j+1)
              fls3 = (tsb(i-1,k,j)+dtb*tsbt(i-1,k,j)) - atm2%t(i-1,k,j)
              fls4 = (tsb(i+1,k,j)+dtb*tsbt(i+1,k,j)) - atm2%t(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (tnb(i,k,j)+dtb*tnbt(i,k,j)) - atm2%t(ii,k,j)
              fls1 = (tnb(i,k,j-1)+dtb*tnbt(i,k,j-1)) - atm2%t(ii,k,j-1)
              fls2 = (tnb(i,k,j+1)+dtb*tnbt(i,k,j+1)) - atm2%t(ii,k,j+1)
              fls3 = (tnb(i-1,k,j)+dtb*tnbt(i-1,k,j)) - atm2%t(ii-1,k,j)
              fls4 = (tnb(i+1,k,j)+dtb*tnbt(i+1,k,j)) - atm2%t(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
         end do
#else
         jm1 = j-1
         jp1 = j+1
         if(jm1 == 0) jm1 = jx
         if(jp1 == jx+1) jp1 = 1
!------interior j slices:
         do i = 2 , ip
            ii = iym1 - i + 1
            do k = 1 , kz
              fcx = fcoef*xfune(i,k)
              gcx = gcoef*xfune(i,k)
!........south boundary:
              fls0 = (tsb(i,k,j)+dtb*tsbt(i,k,j)) - atm2%t(i,k,j)
              fls1 = (tsb(i,k,jm1)+dtb*tsbt(i,k,jm1)) - atm2%t(i,k,jm1)
              fls2 = (tsb(i,k,jp1)+dtb*tsbt(i,k,jp1)) - atm2%t(i,k,jp1)
              fls3 = (tsb(i-1,k,j)+dtb*tsbt(i-1,k,j)) - atm2%t(i-1,k,j)
              fls4 = (tsb(i+1,k,j)+dtb*tsbt(i+1,k,j)) - atm2%t(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (tnb(i,k,j)+dtb*tnbt(i,k,j)) - atm2%t(ii,k,j)
              fls1 = (tnb(i,k,jm1)+dtb*tnbt(i,k,jm1)) - atm2%t(ii,k,jm1)
              fls2 = (tnb(i,k,jp1)+dtb*tnbt(i,k,jp1)) - atm2%t(ii,k,jp1)
              fls3 = (tnb(i-1,k,j)+dtb*tnbt(i-1,k,j)) - atm2%t(ii-1,k,j)
              fls4 = (tnb(i+1,k,j)+dtb*tnbt(i+1,k,j)) - atm2%t(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
         end do
#endif
      end if
#else
!----------------------------------------------------------------------
!
      dtb = xt*minph
#ifdef MPP1
      jsls = j + myid*jxp
      jj = jx - jsls
      if ( jj <= ip ) jsls = jj
      jwb = jsls
      if ( jwb > jxp ) jwb = mod(jwb,jxp)
      if ( jwb == 0 ) jwb = jxp
      if ( myid == nproc-1 ) then
        jeb = jsls
      else
        jeb = jsls + 1
      end if
      if ( jeb > jxp ) jeb = mod(jeb,jxp)
      if ( jeb == 0 ) jeb = jxp
#else
      jsls = j
      jj = jx - jsls
      if ( jj <= ip ) jsls = jj
#endif
!
!-----determine which relaxation method to use:linear/expon.
!
      if ( ibdy == 1 ) then
!
!---------use linear method
!
        if ( jsls > ip ) then
!------interior j slices:
          do i = 2 , ip
            ii = iym1 - i + 1
            fcx = fcoef*xfun(i)
            gcx = gcoef*xfun(i)
            do k = 1 , kz
!.......south boundary:
              fls0 = (tsb(i,k,j)+dtb*tsbt(i,k,j)) - atm2%t(i,k,j)
              fls1 = (tsb(i,k,j-1)+dtb*tsbt(i,k,j-1)) - atm2%t(i,k,j-1)
              fls2 = (tsb(i,k,j+1)+dtb*tsbt(i,k,j+1)) - atm2%t(i,k,j+1)
              fls3 = (tsb(i-1,k,j)+dtb*tsbt(i-1,k,j)) - atm2%t(i-1,k,j)
              fls4 = (tsb(i+1,k,j)+dtb*tsbt(i+1,k,j)) - atm2%t(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (tnb(i,k,j)+dtb*tnbt(i,k,j)) - atm2%t(ii,k,j)
              fls1 = (tnb(i,k,j-1)+dtb*tnbt(i,k,j-1)) - atm2%t(ii,k,j-1)
              fls2 = (tnb(i,k,j+1)+dtb*tnbt(i,k,j+1)) - atm2%t(ii,k,j+1)
              fls3 = (tnb(i-1,k,j)+dtb*tnbt(i-1,k,j)) - atm2%t(ii-1,k,j)
              fls4 = (tnb(i+1,k,j)+dtb*tnbt(i+1,k,j)) - atm2%t(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
!
        else if ( jsls <= ip ) then
!------east or west boundary slices:
          ibeg = 2
          iend = iym1 - 1
          if ( jsls > 2 ) then
            do i = 2 , jsls - 1
              ii = iym1 - i + 1
              fcx = fcoef*xfun(i)
              gcx = gcoef*xfun(i)
              do k = 1 , kz
!........south  boundary:
                fls0 = (tsb(i,k,j)+dtb*tsbt(i,k,j)) - atm2%t(i,k,j)
                fls1 = (tsb(i,k,j-1)+dtb*tsbt(i,k,j-1))-atm2%t(i,k,j-1)
                fls2 = (tsb(i,k,j+1)+dtb*tsbt(i,k,j+1))-atm2%t(i,k,j+1)
                fls3 = (tsb(i-1,k,j)+dtb*tsbt(i-1,k,j))-atm2%t(i-1,k,j)
                fls4 = (tsb(i+1,k,j)+dtb*tsbt(i+1,k,j))-atm2%t(i+1,k,j)
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!.........north boundary:
                fls0 = (tnb(i,k,j)+dtb*tnbt(i,k,j)) - atm2%t(ii,k,j)
                fls1 = (tnb(i,k,j-1)+dtb*tnbt(i,k,j-1))-atm2%t(ii,k,j-1)
                fls2 = (tnb(i,k,j+1)+dtb*tnbt(i,k,j+1))-atm2%t(ii,k,j+1)
                fls3 = (tnb(i-1,k,j)+dtb*tnbt(i-1,k,j))-atm2%t(ii-1,k,j)
                fls4 = (tnb(i+1,k,j)+dtb*tnbt(i+1,k,j))-atm2%t(ii+1,k,j)
                ften(ii,k) = ften(ii,k) + fcx*fls0 -                    &
                           & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
            ibeg = jsls
            iend = iym1 - jsls + 1
          end if
!
          if ( jj > ip ) then
!-------west-boundary slice:
            fcx = fcoef*xfun(jsls)
            gcx = gcoef*xfun(jsls)
            do k = 1 , kz
              do i = ibeg , iend
#ifdef MPP1
                fls0 = (twb(i,k,jwb)+dtb*twbt(i,k,jwb)) - atm2%t(i,k,j)
                fls1 = (twb(i-1,k,jwb)+dtb*twbt(i-1,k,jwb))             &
                     & - atm2%t(i-1,k,j)
                fls2 = (twb(i+1,k,jwb)+dtb*twbt(i+1,k,jwb))             &
                     & - atm2%t(i+1,k,j)
                fls3 = (twb(i,k,jwb-1)+dtb*twbt(i,k,jwb-1))             &
                     & - atm2%t(i,k,j-1)
                fls4 = (twb(i,k,jwb+1)+dtb*twbt(i,k,jwb+1))             &
                     & - atm2%t(i,k,j+1)
#else
                fls0 = (twb(i,k,jsls)+dtb*twbt(i,k,jsls))-atm2%t(i,k,j)
                fls1 = (twb(i-1,k,jsls)+dtb*twbt(i-1,k,jsls))           &
                     & - atm2%t(i-1,k,j)
                fls2 = (twb(i+1,k,jsls)+dtb*twbt(i+1,k,jsls))           &
                     & - atm2%t(i+1,k,j)
                fls3 = (twb(i,k,jsls-1)+dtb*twbt(i,k,jsls-1))           &
                     & - atm2%t(i,k,j-1)
                fls4 = (twb(i,k,jsls+1)+dtb*twbt(i,k,jsls+1))           &
                     & - atm2%t(i,k,j+1)
#endif
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
          else if ( jj <= ip ) then
!-------east-boundary slice:
            fcx = fcoef*xfun(jsls)
            gcx = gcoef*xfun(jsls)
            do k = 1 , kz
              do i = ibeg , iend
#ifdef MPP1
                fls0 = (teb(i,k,jeb)+dtb*tebt(i,k,jeb)) - atm2%t(i,k,j)
                fls1 = (teb(i-1,k,jeb)+dtb*tebt(i-1,k,jeb))             &
                     & - atm2%t(i-1,k,j)
                fls2 = (teb(i+1,k,jeb)+dtb*tebt(i+1,k,jeb))             &
                     & - atm2%t(i+1,k,j)
                fls3 = (teb(i,k,jeb-1)+dtb*tebt(i,k,jeb-1))             &
                     & - atm2%t(i,k,j-1)
                fls4 = (teb(i,k,jeb+1)+dtb*tebt(i,k,jeb+1))             &
                     & - atm2%t(i,k,j+1)
#else
                fls0 = (teb(i,k,jsls)+dtb*tebt(i,k,jsls))-atm2%t(i,k,j)
                fls1 = (teb(i-1,k,jsls)+dtb*tebt(i-1,k,jsls))           &
                     & - atm2%t(i-1,k,j)
                fls2 = (teb(i+1,k,jsls)+dtb*tebt(i+1,k,jsls))           &
                     & - atm2%t(i+1,k,j)
                fls3 = (teb(i,k,jsls-1)+dtb*tebt(i,k,jsls-1))           &
                     & - atm2%t(i,k,j-1)
                fls4 = (teb(i,k,jsls+1)+dtb*tebt(i,k,jsls+1))           &
                     & - atm2%t(i,k,j+1)
#endif
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
          end if
        end if
!
      else if ( ibdy == 5 ) then
 
!----------use exponential method
 
        if ( jsls > ip ) then
!------interior j slices:
          do i = 2 , ip
            ii = iym1 - i + 1
            do k = 1 , kz
              fcx = fcoef*xfune(i,k)
              gcx = gcoef*xfune(i,k)
!........south boundary:
              fls0 = (tsb(i,k,j)+dtb*tsbt(i,k,j)) - atm2%t(i,k,j)
              fls1 = (tsb(i,k,j-1)+dtb*tsbt(i,k,j-1)) - atm2%t(i,k,j-1)
              fls2 = (tsb(i,k,j+1)+dtb*tsbt(i,k,j+1)) - atm2%t(i,k,j+1)
              fls3 = (tsb(i-1,k,j)+dtb*tsbt(i-1,k,j)) - atm2%t(i-1,k,j)
              fls4 = (tsb(i+1,k,j)+dtb*tsbt(i+1,k,j)) - atm2%t(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (tnb(i,k,j)+dtb*tnbt(i,k,j)) - atm2%t(ii,k,j)
              fls1 = (tnb(i,k,j-1)+dtb*tnbt(i,k,j-1)) - atm2%t(ii,k,j-1)
              fls2 = (tnb(i,k,j+1)+dtb*tnbt(i,k,j+1)) - atm2%t(ii,k,j+1)
              fls3 = (tnb(i-1,k,j)+dtb*tnbt(i-1,k,j)) - atm2%t(ii-1,k,j)
              fls4 = (tnb(i+1,k,j)+dtb*tnbt(i+1,k,j)) - atm2%t(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
!
        else if ( jsls <= ip ) then
!------east or west boundary slices:
          ibeg = 2
          iend = iym1 - 1
          if ( jsls > 2 ) then
            do i = 2 , jsls - 1
              ii = iym1 - i + 1
              do k = 1 , kz
                fcx = fcoef*xfune(i,k)
                gcx = gcoef*xfune(i,k)
!.........south boundary:
                fls0 = (tsb(i,k,j)+dtb*tsbt(i,k,j)) - atm2%t(i,k,j)
                fls1 = (tsb(i,k,j-1)+dtb*tsbt(i,k,j-1))-atm2%t(i,k,j-1)
                fls2 = (tsb(i,k,j+1)+dtb*tsbt(i,k,j+1))-atm2%t(i,k,j+1)
                fls3 = (tsb(i-1,k,j)+dtb*tsbt(i-1,k,j))-atm2%t(i-1,k,j)
                fls4 = (tsb(i+1,k,j)+dtb*tsbt(i+1,k,j))-atm2%t(i+1,k,j)
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!.........north boundary:
                fls0 = (tnb(i,k,j)+dtb*tnbt(i,k,j)) - atm2%t(ii,k,j)
                fls1 = (tnb(i,k,j-1)+dtb*tnbt(i,k,j-1))-atm2%t(ii,k,j-1)
                fls2 = (tnb(i,k,j+1)+dtb*tnbt(i,k,j+1))-atm2%t(ii,k,j+1)
                fls3 = (tnb(i-1,k,j)+dtb*tnbt(i-1,k,j))-atm2%t(ii-1,k,j)
                fls4 = (tnb(i+1,k,j)+dtb*tnbt(i+1,k,j))-atm2%t(ii+1,k,j)
                ften(ii,k) = ften(ii,k) + fcx*fls0 -                    &
                           & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
            ibeg = jsls
            iend = iym1 - jsls + 1
          end if
!
          if ( jj > ip ) then
!-------west-boundary slice:
            do k = 1 , kz
              fcx = fcoef*xfune(jsls,k)
              gcx = gcoef*xfune(jsls,k)
              do i = ibeg , iend
#ifdef MPP1
                fls0 = (twb(i,k,jwb)+dtb*twbt(i,k,jwb)) - atm2%t(i,k,j)
                fls1 = (twb(i-1,k,jwb)+dtb*twbt(i-1,k,jwb))             &
                     & - atm2%t(i-1,k,j)
                fls2 = (twb(i+1,k,jwb)+dtb*twbt(i+1,k,jwb))             &
                     & - atm2%t(i+1,k,j)
                fls3 = (twb(i,k,jwb-1)+dtb*twbt(i,k,jwb-1))             &
                     & - atm2%t(i,k,j-1)
                fls4 = (twb(i,k,jwb+1)+dtb*twbt(i,k,jwb+1))             &
                     & - atm2%t(i,k,j+1)
#else
                fls0 = (twb(i,k,jsls)+dtb*twbt(i,k,jsls))-atm2%t(i,k,j)
                fls1 = (twb(i-1,k,jsls)+dtb*twbt(i-1,k,jsls))           &
                     & - atm2%t(i-1,k,j)
                fls2 = (twb(i+1,k,jsls)+dtb*twbt(i+1,k,jsls))           &
                     & - atm2%t(i+1,k,j)
                fls3 = (twb(i,k,jsls-1)+dtb*twbt(i,k,jsls-1))           &
                     & - atm2%t(i,k,j-1)
                fls4 = (twb(i,k,jsls+1)+dtb*twbt(i,k,jsls+1))           &
                     & - atm2%t(i,k,j+1)
#endif
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
          else if ( jj <= ip ) then
!-------east-boundary slice:
            do k = 1 , kz
              fcx = fcoef*xfune(jsls,k)
              gcx = gcoef*xfune(jsls,k)
              do i = ibeg , iend
#ifdef MPP1
                fls0 = (teb(i,k,jeb)+dtb*tebt(i,k,jeb)) - atm2%t(i,k,j)
                fls1 = (teb(i-1,k,jeb)+dtb*tebt(i-1,k,jeb))             &
                     & - atm2%t(i-1,k,j)
                fls2 = (teb(i+1,k,jeb)+dtb*tebt(i+1,k,jeb))             &
                     & - atm2%t(i+1,k,j)
                fls3 = (teb(i,k,jeb-1)+dtb*tebt(i,k,jeb-1))             &
                     & - atm2%t(i,k,j-1)
                fls4 = (teb(i,k,jeb+1)+dtb*tebt(i,k,jeb+1))             &
                     & - atm2%t(i,k,j+1)
#else
                fls0 = (teb(i,k,jsls)+dtb*tebt(i,k,jsls))-atm2%t(i,k,j)
                fls1 = (teb(i-1,k,jsls)+dtb*tebt(i-1,k,jsls))           &
                     & - atm2%t(i-1,k,j)
                fls2 = (teb(i+1,k,jsls)+dtb*tebt(i+1,k,jsls))           &
                     & - atm2%t(i+1,k,j)
                fls3 = (teb(i,k,jsls-1)+dtb*tebt(i,k,jsls-1))           &
                     & - atm2%t(i,k,j-1)
                fls4 = (teb(i,k,jsls+1)+dtb*tebt(i,k,jsls+1))           &
                     & - atm2%t(i,k,j+1)
#endif
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
          end if
        end if
!
      end if
#endif
      call time_end(subroutine_name,idindx)
      end subroutine nudge_t
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine nudgeqv(ip,fcoef,gcoef,xt,ften,j,ibdy)
!
      implicit none
!
      real(8) :: fcoef , gcoef , xt
      integer :: ibdy , ip , j
      real(8) , dimension(iy,kz) :: ften
      intent (in) fcoef , gcoef , ibdy , ip , j , xt
      intent (inout) ften
!
      real(8) :: dtb , fcx , fls0 , fls1 , fls2 , fls3 , fls4 , gcx
      integer :: i , ii , k
#ifndef BAND
      integer :: ibeg , iend , jj , jsls
#else
#ifndef MPP1
      integer :: jp1 , jm1
#endif
#endif
#ifdef MPP1
#ifndef BAND
      integer :: jwb , jeb
#endif
#endif
!
      cHARACTER (len=50) :: subroutine_name='nudgeqv'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
#ifdef BAND
!----------------------------------------------------------------------
!
      dtb = xt*minph
!
!-----determine which relaxation method to use:linear/expon.
!
      if ( ibdy == 1 ) then
!
!---------use linear method
!
#ifdef MPP1
!------interior j slices:
         do i = 2 , ip
            ii = iym1 - i + 1
            fcx = fcoef*xfun(i)
            gcx = gcoef*xfun(i)
            do k = 1 , kz
!.......south boundary:
              fls0 = (qsb(i,k,j)+dtb*qsbt(i,k,j)) - atm2%qv(i,k,j)
              fls1 = (qsb(i,k,j-1)+dtb*qsbt(i,k,j-1)) - atm2%qv(i,k,j-1)
              fls2 = (qsb(i,k,j+1)+dtb*qsbt(i,k,j+1)) - atm2%qv(i,k,j+1)
              fls3 = (qsb(i-1,k,j)+dtb*qsbt(i-1,k,j)) - atm2%qv(i-1,k,j)
              fls4 = (qsb(i+1,k,j)+dtb*qsbt(i+1,k,j)) - atm2%qv(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (qnb(i,k,j)+dtb*qnbt(i,k,j)) - atm2%qv(ii,k,j)
              fls1 = (qnb(i,k,j-1)+dtb*qnbt(i,k,j-1))-atm2%qv(ii,k,j-1)
              fls2 = (qnb(i,k,j+1)+dtb*qnbt(i,k,j+1))-atm2%qv(ii,k,j+1)
              fls3 = (qnb(i-1,k,j)+dtb*qnbt(i-1,k,j))-atm2%qv(ii-1,k,j)
              fls4 = (qnb(i+1,k,j)+dtb*qnbt(i+1,k,j))-atm2%qv(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
         end do
#else
         jm1 = j-1
         jp1 = j+1
         if(jm1 == 0) jm1 = jx
         if(jp1 == jx+1) jp1 = 1
!------interior j slices:
         do i = 2 , ip
            ii = iym1 - i + 1
            fcx = fcoef*xfun(i)
            gcx = gcoef*xfun(i)
            do k = 1 , kz
!.......south boundary:
              fls0 = (qsb(i,k,j)+dtb*qsbt(i,k,j)) - atm2%qv(i,k,j)
              fls1 = (qsb(i,k,jm1)+dtb*qsbt(i,k,jm1)) - atm2%qv(i,k,jm1)
              fls2 = (qsb(i,k,jp1)+dtb*qsbt(i,k,jp1)) - atm2%qv(i,k,jp1)
              fls3 = (qsb(i-1,k,j)+dtb*qsbt(i-1,k,j)) - atm2%qv(i-1,k,j)
              fls4 = (qsb(i+1,k,j)+dtb*qsbt(i+1,k,j)) - atm2%qv(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (qnb(i,k,j)+dtb*qnbt(i,k,j)) - atm2%qv(ii,k,j)
              fls1 = (qnb(i,k,jm1)+dtb*qnbt(i,k,jm1))-atm2%qv(ii,k,jm1)
              fls2 = (qnb(i,k,jp1)+dtb*qnbt(i,k,jp1))-atm2%qv(ii,k,jp1)
              fls3 = (qnb(i-1,k,j)+dtb*qnbt(i-1,k,j))-atm2%qv(ii-1,k,j)
              fls4 = (qnb(i+1,k,j)+dtb*qnbt(i+1,k,j))-atm2%qv(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
         end do
#endif
      else if ( ibdy == 5 ) then
 
!----------use exponential method
 
#ifdef MPP1
!------interior j slices:
         do i = 2 , ip
            ii = iym1 - i + 1
            do k = 1 , kz
              fcx = fcoef*xfune(i,k)
              gcx = gcoef*xfune(i,k)
!........south boundary:
              fls0 = (qsb(i,k,j)+dtb*qsbt(i,k,j)) - atm2%qv(i,k,j)
              fls1 = (qsb(i,k,j-1)+dtb*qsbt(i,k,j-1)) - atm2%qv(i,k,j-1)
              fls2 = (qsb(i,k,j+1)+dtb*qsbt(i,k,j+1)) - atm2%qv(i,k,j+1)
              fls3 = (qsb(i-1,k,j)+dtb*qsbt(i-1,k,j)) - atm2%qv(i-1,k,j)
              fls4 = (qsb(i+1,k,j)+dtb*qsbt(i+1,k,j)) - atm2%qv(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (qnb(i,k,j)+dtb*qnbt(i,k,j)) - atm2%qv(ii,k,j)
              fls1 = (qnb(i,k,j-1)+dtb*qnbt(i,k,j-1))-atm2%qv(ii,k,j-1)
              fls2 = (qnb(i,k,j+1)+dtb*qnbt(i,k,j+1))-atm2%qv(ii,k,j+1)
              fls3 = (qnb(i-1,k,j)+dtb*qnbt(i-1,k,j))-atm2%qv(ii-1,k,j)
              fls4 = (qnb(i+1,k,j)+dtb*qnbt(i+1,k,j))-atm2%qv(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
         end do
#else
         jm1 = j-1
         jp1 = j+1
         if(jm1 == 0) jm1 = jx
         if(jp1 == jx+1) jp1 = 1
!------interior j slices:
         do i = 2 , ip
            ii = iym1 - i + 1
            do k = 1 , kz
              fcx = fcoef*xfune(i,k)
              gcx = gcoef*xfune(i,k)
!........south boundary:
              fls0 = (qsb(i,k,j)+dtb*qsbt(i,k,j)) - atm2%qv(i,k,j)
              fls1 = (qsb(i,k,jm1)+dtb*qsbt(i,k,jm1)) - atm2%qv(i,k,jm1)
              fls2 = (qsb(i,k,jp1)+dtb*qsbt(i,k,jp1)) - atm2%qv(i,k,jp1)
              fls3 = (qsb(i-1,k,j)+dtb*qsbt(i-1,k,j)) - atm2%qv(i-1,k,j)
              fls4 = (qsb(i+1,k,j)+dtb*qsbt(i+1,k,j)) - atm2%qv(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (qnb(i,k,j)+dtb*qnbt(i,k,j)) - atm2%qv(ii,k,j)
              fls1 = (qnb(i,k,jm1)+dtb*qnbt(i,k,jm1))-atm2%qv(ii,k,jm1)
              fls2 = (qnb(i,k,jp1)+dtb*qnbt(i,k,jp1))-atm2%qv(ii,k,jp1)
              fls3 = (qnb(i-1,k,j)+dtb*qnbt(i-1,k,j))-atm2%qv(ii-1,k,j)
              fls4 = (qnb(i+1,k,j)+dtb*qnbt(i+1,k,j))-atm2%qv(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
         end do
#endif
      end if
#else
!----------------------------------------------------------------------
!
      dtb = xt*minph
#ifdef MPP1
      jsls = j + myid*jxp
      jj = jx - jsls
      if ( jj <= ip ) jsls = jj
      jwb = jsls
      if ( jwb > jxp ) jwb = mod(jwb,jxp)
      if ( jwb == 0 ) jwb = jxp
      if ( myid == nproc-1 ) then
        jeb = jsls
      else
        jeb = jsls + 1
      end if
      if ( jeb > jxp ) jeb = mod(jeb,jxp)
      if ( jeb == 0 ) jeb = jxp
#else
      jsls = j
      jj = jx - jsls
      if ( jj <= ip ) jsls = jj
#endif
!
!-----determine which relaxation method to use:linear/expon.
!
      if ( ibdy == 1 ) then
!
!---------use linear method
!
        if ( jsls > ip ) then
!------interior j slices:
          do i = 2 , ip
            ii = iym1 - i + 1
            fcx = fcoef*xfun(i)
            gcx = gcoef*xfun(i)
            do k = 1 , kz
!.......south boundary:
              fls0 = (qsb(i,k,j)+dtb*qsbt(i,k,j)) - atm2%qv(i,k,j)
              fls1 = (qsb(i,k,j-1)+dtb*qsbt(i,k,j-1)) - atm2%qv(i,k,j-1)
              fls2 = (qsb(i,k,j+1)+dtb*qsbt(i,k,j+1)) - atm2%qv(i,k,j+1)
              fls3 = (qsb(i-1,k,j)+dtb*qsbt(i-1,k,j)) - atm2%qv(i-1,k,j)
              fls4 = (qsb(i+1,k,j)+dtb*qsbt(i+1,k,j)) - atm2%qv(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (qnb(i,k,j)+dtb*qnbt(i,k,j)) - atm2%qv(ii,k,j)
              fls1 = (qnb(i,k,j-1)+dtb*qnbt(i,k,j-1))-atm2%qv(ii,k,j-1)
              fls2 = (qnb(i,k,j+1)+dtb*qnbt(i,k,j+1))-atm2%qv(ii,k,j+1)
              fls3 = (qnb(i-1,k,j)+dtb*qnbt(i-1,k,j))-atm2%qv(ii-1,k,j)
              fls4 = (qnb(i+1,k,j)+dtb*qnbt(i+1,k,j))-atm2%qv(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
!
        else if ( jsls <= ip ) then
!------east or west boundary slices:
          ibeg = 2
          iend = iym1 - 1
          if ( jsls > 2 ) then
            do i = 2 , jsls - 1
              ii = iym1 - i + 1
              fcx = fcoef*xfun(i)
              gcx = gcoef*xfun(i)
              do k = 1 , kz
!........south  boundary:
                fls0 = (qsb(i,k,j)+dtb*qsbt(i,k,j)) - atm2%qv(i,k,j)
                fls1 = (qsb(i,k,j-1)+dtb*qsbt(i,k,j-1))-atm2%qv(i,k,j-1)
                fls2 = (qsb(i,k,j+1)+dtb*qsbt(i,k,j+1))-atm2%qv(i,k,j+1)
                fls3 = (qsb(i-1,k,j)+dtb*qsbt(i-1,k,j))-atm2%qv(i-1,k,j)
                fls4 = (qsb(i+1,k,j)+dtb*qsbt(i+1,k,j))-atm2%qv(i+1,k,j)
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!.........north boundary:
                fls0 = (qnb(i,k,j)+dtb*qnbt(i,k,j)) - atm2%qv(ii,k,j)
                fls1 = (qnb(i,k,j-1)+dtb*qnbt(i,k,j-1))- &
                       atm2%qv(ii,k,j-1)
                fls2 = (qnb(i,k,j+1)+dtb*qnbt(i,k,j+1))- &
                       atm2%qv(ii,k,j+1)
                fls3 = (qnb(i-1,k,j)+dtb*qnbt(i-1,k,j))- &
                       atm2%qv(ii-1,k,j)
                fls4 = (qnb(i+1,k,j)+dtb*qnbt(i+1,k,j))- &
                       atm2%qv(ii+1,k,j)
                ften(ii,k) = ften(ii,k) + fcx*fls0 -                    &
                           & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
            ibeg = jsls
            iend = iym1 - jsls + 1
          end if
!
          if ( jj > ip ) then
!-------west-boundary slice:
            fcx = fcoef*xfun(jsls)
            gcx = gcoef*xfun(jsls)
            do k = 1 , kz
              do i = ibeg , iend
#ifdef MPP1
                fls0 = (qwb(i,k,jwb)+dtb*qwbt(i,k,jwb)) - atm2%qv(i,k,j)
                fls1 = (qwb(i-1,k,jwb)+dtb*qwbt(i-1,k,jwb))             &
                     & - atm2%qv(i-1,k,j)
                fls2 = (qwb(i+1,k,jwb)+dtb*qwbt(i+1,k,jwb))             &
                     & - atm2%qv(i+1,k,j)
                fls3 = (qwb(i,k,jwb-1)+dtb*qwbt(i,k,jwb-1))             &
                     & - atm2%qv(i,k,j-1)
                fls4 = (qwb(i,k,jwb+1)+dtb*qwbt(i,k,jwb+1))             &
                     & - atm2%qv(i,k,j+1)
#else
                fls0 = (qwb(i,k,jsls)+dtb*qwbt(i,k,jsls))-atm2%qv(i,k,j)
                fls1 = (qwb(i-1,k,jsls)+dtb*qwbt(i-1,k,jsls))           &
                     & - atm2%qv(i-1,k,j)
                fls2 = (qwb(i+1,k,jsls)+dtb*qwbt(i+1,k,jsls))           &
                     & - atm2%qv(i+1,k,j)
                fls3 = (qwb(i,k,jsls-1)+dtb*qwbt(i,k,jsls-1))           &
                     & - atm2%qv(i,k,j-1)
                fls4 = (qwb(i,k,jsls+1)+dtb*qwbt(i,k,jsls+1))           &
                     & - atm2%qv(i,k,j+1)
#endif
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
          else if ( jj <= ip ) then
!-------east-boundary slice:
            fcx = fcoef*xfun(jsls)
            gcx = gcoef*xfun(jsls)
            do k = 1 , kz
              do i = ibeg , iend
#ifdef MPP1
                fls0 = (qeb(i,k,jeb)+dtb*qebt(i,k,jeb)) - atm2%qv(i,k,j)
                fls1 = (qeb(i-1,k,jeb)+dtb*qebt(i-1,k,jeb))             &
                     & - atm2%qv(i-1,k,j)
                fls2 = (qeb(i+1,k,jeb)+dtb*qebt(i+1,k,jeb))             &
                     & - atm2%qv(i+1,k,j)
                fls3 = (qeb(i,k,jeb-1)+dtb*qebt(i,k,jeb-1))             &
                     & - atm2%qv(i,k,j-1)
                fls4 = (qeb(i,k,jeb+1)+dtb*qebt(i,k,jeb+1))             &
                     & - atm2%qv(i,k,j+1)
#else
                fls0 = (qeb(i,k,jsls)+dtb*qebt(i,k,jsls))-atm2%qv(i,k,j)
                fls1 = (qeb(i-1,k,jsls)+dtb*qebt(i-1,k,jsls))           &
                     & - atm2%qv(i-1,k,j)
                fls2 = (qeb(i+1,k,jsls)+dtb*qebt(i+1,k,jsls))           &
                     & - atm2%qv(i+1,k,j)
                fls3 = (qeb(i,k,jsls-1)+dtb*qebt(i,k,jsls-1))           &
                     & - atm2%qv(i,k,j-1)
                fls4 = (qeb(i,k,jsls+1)+dtb*qebt(i,k,jsls+1))           &
                     & - atm2%qv(i,k,j+1)
#endif
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
          end if
        end if
!
      else if ( ibdy == 5 ) then
 
!----------use exponential method
 
        if ( jsls > ip ) then
!------interior j slices:
          do i = 2 , ip
            ii = iym1 - i + 1
            do k = 1 , kz
              fcx = fcoef*xfune(i,k)
              gcx = gcoef*xfune(i,k)
!........south boundary:
              fls0 = (qsb(i,k,j)+dtb*qsbt(i,k,j)) - atm2%qv(i,k,j)
              fls1 = (qsb(i,k,j-1)+dtb*qsbt(i,k,j-1)) - atm2%qv(i,k,j-1)
              fls2 = (qsb(i,k,j+1)+dtb*qsbt(i,k,j+1)) - atm2%qv(i,k,j+1)
              fls3 = (qsb(i-1,k,j)+dtb*qsbt(i-1,k,j)) - atm2%qv(i-1,k,j)
              fls4 = (qsb(i+1,k,j)+dtb*qsbt(i+1,k,j)) - atm2%qv(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (qnb(i,k,j)+dtb*qnbt(i,k,j)) - atm2%qv(ii,k,j)
              fls1 = (qnb(i,k,j-1)+dtb*qnbt(i,k,j-1))-atm2%qv(ii,k,j-1)
              fls2 = (qnb(i,k,j+1)+dtb*qnbt(i,k,j+1))-atm2%qv(ii,k,j+1)
              fls3 = (qnb(i-1,k,j)+dtb*qnbt(i-1,k,j))-atm2%qv(ii-1,k,j)
              fls4 = (qnb(i+1,k,j)+dtb*qnbt(i+1,k,j))-atm2%qv(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
!
        else if ( jsls <= ip ) then
!------east or west boundary slices:
          ibeg = 2
          iend = iym1 - 1
          if ( jsls > 2 ) then
            do i = 2 , jsls - 1
              ii = iym1 - i + 1
              do k = 1 , kz
                fcx = fcoef*xfune(i,k)
                gcx = gcoef*xfune(i,k)
!.........south boundary:
                fls0 = (qsb(i,k,j)+dtb*qsbt(i,k,j)) - atm2%qv(i,k,j)
                fls1 = (qsb(i,k,j-1)+dtb*qsbt(i,k,j-1))-atm2%qv(i,k,j-1)
                fls2 = (qsb(i,k,j+1)+dtb*qsbt(i,k,j+1))-atm2%qv(i,k,j+1)
                fls3 = (qsb(i-1,k,j)+dtb*qsbt(i-1,k,j))-atm2%qv(i-1,k,j)
                fls4 = (qsb(i+1,k,j)+dtb*qsbt(i+1,k,j))-atm2%qv(i+1,k,j)
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!.........north boundary:
                fls0 = (qnb(i,k,j)+dtb*qnbt(i,k,j)) - atm2%qv(ii,k,j)
                fls1 = (qnb(i,k,j-1)+dtb*qnbt(i,k,j-1)) - &
                       atm2%qv(ii,k,j-1)
                fls2 = (qnb(i,k,j+1)+dtb*qnbt(i,k,j+1)) - &
                       atm2%qv(ii,k,j+1)
                fls3 = (qnb(i-1,k,j)+dtb*qnbt(i-1,k,j)) - &
                       atm2%qv(ii-1,k,j)
                fls4 = (qnb(i+1,k,j)+dtb*qnbt(i+1,k,j)) - &
                       atm2%qv(ii+1,k,j)
                ften(ii,k) = ften(ii,k) + fcx*fls0 -                    &
                           & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
            ibeg = jsls
            iend = iym1 - jsls + 1
          end if
!
          if ( jj > ip ) then
!-------west-boundary slice:
            do k = 1 , kz
              fcx = fcoef*xfune(jsls,k)
              gcx = gcoef*xfune(jsls,k)
              do i = ibeg , iend
#ifdef MPP1
                fls0 = (qwb(i,k,jwb)+dtb*qwbt(i,k,jwb)) - atm2%qv(i,k,j)
                fls1 = (qwb(i-1,k,jwb)+dtb*qwbt(i-1,k,jwb))             &
                     & - atm2%qv(i-1,k,j)
                fls2 = (qwb(i+1,k,jwb)+dtb*qwbt(i+1,k,jwb))             &
                     & - atm2%qv(i+1,k,j)
                fls3 = (qwb(i,k,jwb-1)+dtb*qwbt(i,k,jwb-1))             &
                     & - atm2%qv(i,k,j-1)
                fls4 = (qwb(i,k,jwb+1)+dtb*qwbt(i,k,jwb+1))             &
                     & - atm2%qv(i,k,j+1)
#else
                fls0 = (qwb(i,k,jsls)+dtb*qwbt(i,k,jsls))-atm2%qv(i,k,j)
                fls1 = (qwb(i-1,k,jsls)+dtb*qwbt(i-1,k,jsls))           &
                     & - atm2%qv(i-1,k,j)
                fls2 = (qwb(i+1,k,jsls)+dtb*qwbt(i+1,k,jsls))           &
                     & - atm2%qv(i+1,k,j)
                fls3 = (qwb(i,k,jsls-1)+dtb*qwbt(i,k,jsls-1))           &
                     & - atm2%qv(i,k,j-1)
                fls4 = (qwb(i,k,jsls+1)+dtb*qwbt(i,k,jsls+1))           &
                     & - atm2%qv(i,k,j+1)
#endif
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
          else if ( jj <= ip ) then
!-------east-boundary slice:
            do k = 1 , kz
              fcx = fcoef*xfune(jsls,k)
              gcx = gcoef*xfune(jsls,k)
              do i = ibeg , iend
#ifdef MPP1
                fls0 = (qeb(i,k,jeb)+dtb*qebt(i,k,jeb)) - atm2%qv(i,k,j)
                fls1 = (qeb(i-1,k,jeb)+dtb*qebt(i-1,k,jeb))             &
                     & - atm2%qv(i-1,k,j)
                fls2 = (qeb(i+1,k,jeb)+dtb*qebt(i+1,k,jeb))             &
                     & - atm2%qv(i+1,k,j)
                fls3 = (qeb(i,k,jeb-1)+dtb*qebt(i,k,jeb-1))             &
                     & - atm2%qv(i,k,j-1)
                fls4 = (qeb(i,k,jeb+1)+dtb*qebt(i,k,jeb+1))             &
                     & - atm2%qv(i,k,j+1)
#else
                fls0 = (qeb(i,k,jsls)+dtb*qebt(i,k,jsls))-atm2%qv(i,k,j)
                fls1 = (qeb(i-1,k,jsls)+dtb*qebt(i-1,k,jsls))           &
                     & - atm2%qv(i-1,k,j)
                fls2 = (qeb(i+1,k,jsls)+dtb*qebt(i+1,k,jsls))           &
                     & - atm2%qv(i+1,k,j)
                fls3 = (qeb(i,k,jsls-1)+dtb*qebt(i,k,jsls-1))           &
                     & - atm2%qv(i,k,j-1)
                fls4 = (qeb(i,k,jsls+1)+dtb*qebt(i,k,jsls+1))           &
                     & - atm2%qv(i,k,j+1)
#endif
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
          end if
        end if
!
      end if
#endif
      call time_end(subroutine_name,idindx)
      end subroutine nudgeqv
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine nudge_u(ip,fcoef,gcoef,xt,ften,j,ibdy)
!
      implicit none
!
      real(8) :: fcoef , gcoef , xt
      integer :: ibdy , ip , j
      real(8) , dimension(iy,kz) :: ften
      intent (in) fcoef , gcoef , ibdy , ip , j , xt
      intent (inout) ften
!
      real(8) :: dtb , fcx , fls0 , fls1 , fls2 , fls3 , fls4 , gcx
      integer :: i , ii , k
#ifndef BAND
      integer :: ibeg , iend , jj , jsls
#else
#ifndef MPP1
      integer :: jp1 , jm1
#endif
#endif
#ifdef MPP1
#ifndef BAND
      integer :: jew
#endif
#endif
      character (len=50) :: subroutine_name='nudge_u'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
#ifdef BAND
!----------------------------------------------------------------------
!
      dtb = xt*minph
!
!-----determine which relaxation method to use:linear/expon.
!
      if ( ibdy == 1 ) then
!
!---------use linear method
!
#ifdef MPP1
!------interior j slices:
         do i = 2 , ip
            ii = iy - i + 1
            fcx = fcoef*xfun(i)
            gcx = gcoef*xfun(i)
            do k = 1 , kz
!.......south boundary:
              fls0 = (usb(i,k,j)+dtb*usbt(i,k,j)) - atm2%u(i,k,j)
              fls1 = (usb(i,k,j-1)+dtb*usbt(i,k,j-1)) - atm2%u(i,k,j-1)
              fls2 = (usb(i,k,j+1)+dtb*usbt(i,k,j+1)) - atm2%u(i,k,j+1)
              fls3 = (usb(i-1,k,j)+dtb*usbt(i-1,k,j)) - atm2%u(i-1,k,j)
              fls4 = (usb(i+1,k,j)+dtb*usbt(i+1,k,j)) - atm2%u(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (unb(i,k,j)+dtb*unbt(i,k,j)) - atm2%u(ii,k,j)
              fls1 = (unb(i,k,j-1)+dtb*unbt(i,k,j-1)) - atm2%u(ii,k,j-1)
              fls2 = (unb(i,k,j+1)+dtb*unbt(i,k,j+1)) - atm2%u(ii,k,j+1)
              fls3 = (unb(i-1,k,j)+dtb*unbt(i-1,k,j)) - atm2%u(ii-1,k,j)
              fls4 = (unb(i+1,k,j)+dtb*unbt(i+1,k,j)) - atm2%u(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
         end do
#else
         jm1 = j-1
         jp1 = j+1
         if(jm1 == 0) jm1 = jx
         if(jp1 == jx+1) jp1 = 1
!------interior j slices:
         do i = 2 , ip
            ii = iy - i + 1
            fcx = fcoef*xfun(i)
            gcx = gcoef*xfun(i)
            do k = 1 , kz
!.......south boundary:
              fls0 = (usb(i,k,j)+dtb*usbt(i,k,j)) - atm2%u(i,k,j)
              fls1 = (usb(i,k,jm1)+dtb*usbt(i,k,jm1)) - atm2%u(i,k,jm1)
              fls2 = (usb(i,k,jp1)+dtb*usbt(i,k,jp1)) - atm2%u(i,k,jp1)
              fls3 = (usb(i-1,k,j)+dtb*usbt(i-1,k,j)) - atm2%u(i-1,k,j)
              fls4 = (usb(i+1,k,j)+dtb*usbt(i+1,k,j)) - atm2%u(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (unb(i,k,j)+dtb*unbt(i,k,j)) - atm2%u(ii,k,j)
              fls1 = (unb(i,k,jm1)+dtb*unbt(i,k,jm1)) - atm2%u(ii,k,jm1)
              fls2 = (unb(i,k,jp1)+dtb*unbt(i,k,jp1)) - atm2%u(ii,k,jp1)
              fls3 = (unb(i-1,k,j)+dtb*unbt(i-1,k,j)) - atm2%u(ii-1,k,j)
              fls4 = (unb(i+1,k,j)+dtb*unbt(i+1,k,j)) - atm2%u(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
         end do
#endif
      else if ( ibdy == 5 ) then
 
!----------use exponential method
#ifdef MPP1
!------interior j slices:
         do i = 2 , ip
            ii = iy - i + 1
            do k = 1 , kz
              fcx = fcoef*xfune(i,k)
              gcx = gcoef*xfune(i,k)
!........south boundary:
              fls0 = (usb(i,k,j)+dtb*usbt(i,k,j)) - atm2%u(i,k,j)
              fls1 = (usb(i,k,j-1)+dtb*usbt(i,k,j-1)) - atm2%u(i,k,j-1)
              fls2 = (usb(i,k,j+1)+dtb*usbt(i,k,j+1)) - atm2%u(i,k,j+1)
              fls3 = (usb(i-1,k,j)+dtb*usbt(i-1,k,j)) - atm2%u(i-1,k,j)
              fls4 = (usb(i+1,k,j)+dtb*usbt(i+1,k,j)) - atm2%u(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (unb(i,k,j)+dtb*unbt(i,k,j)) - atm2%u(ii,k,j)
              fls1 = (unb(i,k,j-1)+dtb*unbt(i,k,j-1)) - atm2%u(ii,k,j-1)
              fls2 = (unb(i,k,j+1)+dtb*unbt(i,k,j+1)) - atm2%u(ii,k,j+1)
              fls3 = (unb(i-1,k,j)+dtb*unbt(i-1,k,j)) - atm2%u(ii-1,k,j)
              fls4 = (unb(i+1,k,j)+dtb*unbt(i+1,k,j)) - atm2%u(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
         end do
#else
         jm1 = j-1
         jp1 = j+1
         if(jm1 == 0) jm1 = jx
         if(jp1 == jx+1) jp1 = 1
!------interior j slices:
         do i = 2 , ip
            ii = iy - i + 1
            do k = 1 , kz
              fcx = fcoef*xfune(i,k)
              gcx = gcoef*xfune(i,k)
!........south boundary:
              fls0 = (usb(i,k,j)+dtb*usbt(i,k,j)) - atm2%u(i,k,j)
              fls1 = (usb(i,k,jm1)+dtb*usbt(i,k,jm1)) - atm2%u(i,k,jm1)
              fls2 = (usb(i,k,jp1)+dtb*usbt(i,k,jp1)) - atm2%u(i,k,jp1)
              fls3 = (usb(i-1,k,j)+dtb*usbt(i-1,k,j)) - atm2%u(i-1,k,j)
              fls4 = (usb(i+1,k,j)+dtb*usbt(i+1,k,j)) - atm2%u(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (unb(i,k,j)+dtb*unbt(i,k,j)) - atm2%u(ii,k,j)
              fls1 = (unb(i,k,jm1)+dtb*unbt(i,k,jm1)) - atm2%u(ii,k,jm1)
              fls2 = (unb(i,k,jp1)+dtb*unbt(i,k,jp1)) - atm2%u(ii,k,jp1)
              fls3 = (unb(i-1,k,j)+dtb*unbt(i-1,k,j)) - atm2%u(ii-1,k,j)
              fls4 = (unb(i+1,k,j)+dtb*unbt(i+1,k,j)) - atm2%u(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
         end do
#endif
      end if
#else
!----------------------------------------------------------------------
!
      dtb = xt*minph
#ifdef MPP1
      jsls = j + myid*jxp
      jj = jxp1 - jsls
      if ( jj <= ip ) jsls = jj
      jew = jsls
      if ( jew > jxp ) jew = mod(jsls,jxp)
      if ( jew == 0 ) jew = jxp
#else
      jsls = j
      jj = jxp1 - jsls
      if ( jj <= ip ) jsls = jj
#endif
!
!-----determine which relaxation method to use:linear/expon.
!
      if ( ibdy == 1 ) then
!
!---------use linear method
!
        if ( jsls > ip ) then
!------interior j slices:
          do i = 2 , ip
            ii = iy - i + 1
            fcx = fcoef*xfun(i)
            gcx = gcoef*xfun(i)
            do k = 1 , kz
!.......south boundary:
              fls0 = (usb(i,k,j)+dtb*usbt(i,k,j)) - atm2%u(i,k,j)
              fls1 = (usb(i,k,j-1)+dtb*usbt(i,k,j-1)) - atm2%u(i,k,j-1)
              fls2 = (usb(i,k,j+1)+dtb*usbt(i,k,j+1)) - atm2%u(i,k,j+1)
              fls3 = (usb(i-1,k,j)+dtb*usbt(i-1,k,j)) - atm2%u(i-1,k,j)
              fls4 = (usb(i+1,k,j)+dtb*usbt(i+1,k,j)) - atm2%u(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (unb(i,k,j)+dtb*unbt(i,k,j)) - atm2%u(ii,k,j)
              fls1 = (unb(i,k,j-1)+dtb*unbt(i,k,j-1)) - atm2%u(ii,k,j-1)
              fls2 = (unb(i,k,j+1)+dtb*unbt(i,k,j+1)) - atm2%u(ii,k,j+1)
              fls3 = (unb(i-1,k,j)+dtb*unbt(i-1,k,j)) - atm2%u(ii-1,k,j)
              fls4 = (unb(i+1,k,j)+dtb*unbt(i+1,k,j)) - atm2%u(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
!
        else if ( jsls <= ip ) then
!------east or west boundary slices:
          ibeg = 2
          iend = iym1
          if ( jsls > 2 ) then
            do i = 2 , jsls - 1
              ii = iy - i + 1
              fcx = fcoef*xfun(i)
              gcx = gcoef*xfun(i)
              do k = 1 , kz
!........south  boundary:
                fls0 = (usb(i,k,j)+dtb*usbt(i,k,j)) - atm2%u(i,k,j)
                fls1 = (usb(i,k,j-1)+dtb*usbt(i,k,j-1))- atm2%u(i,k,j-1)
                fls2 = (usb(i,k,j+1)+dtb*usbt(i,k,j+1))- atm2%u(i,k,j+1)
                fls3 = (usb(i-1,k,j)+dtb*usbt(i-1,k,j))- atm2%u(i-1,k,j)
                fls4 = (usb(i+1,k,j)+dtb*usbt(i+1,k,j))- atm2%u(i+1,k,j)
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!.........north boundary:
                fls0 = (unb(i,k,j)+dtb*unbt(i,k,j)) - atm2%u(ii,k,j)
                fls1 = (unb(i,k,j-1)+dtb*unbt(i,k,j-1))-atm2%u(ii,k,j-1)
                fls2 = (unb(i,k,j+1)+dtb*unbt(i,k,j+1))-atm2%u(ii,k,j+1)
                fls3 = (unb(i-1,k,j)+dtb*unbt(i-1,k,j))-atm2%u(ii-1,k,j)
                fls4 = (unb(i+1,k,j)+dtb*unbt(i+1,k,j))-atm2%u(ii+1,k,j)
                ften(ii,k) = ften(ii,k) + fcx*fls0 -                    &
                           & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
            ibeg = jsls
            iend = iy - jsls + 1
          end if
!
          if ( jj > ip ) then
!-------west-boundary slice:
            fcx = fcoef*xfun(jsls)
            gcx = gcoef*xfun(jsls)
            do k = 1 , kz
              do i = ibeg , iend
#ifdef MPP1
                fls0 = (uwb(i,k,jew)+dtb*uwbt(i,k,jew)) - atm2%u(i,k,j)
                fls1 = (uwb(i-1,k,jew)+dtb*uwbt(i-1,k,jew))             &
                     & - atm2%u(i-1,k,j)
                fls2 = (uwb(i+1,k,jew)+dtb*uwbt(i+1,k,jew))             &
                     & - atm2%u(i+1,k,j)
                fls3 = (uwb(i,k,jew-1)+dtb*uwbt(i,k,jew-1))             &
                     & - atm2%u(i,k,j-1)
                fls4 = (uwb(i,k,jew+1)+dtb*uwbt(i,k,jew+1))             &
                     & - atm2%u(i,k,j+1)
#else
                fls0 = (uwb(i,k,jsls)+dtb*uwbt(i,k,jsls))-atm2%u(i,k,j)
                fls1 = (uwb(i-1,k,jsls)+dtb*uwbt(i-1,k,jsls))           &
                     & - atm2%u(i-1,k,j)
                fls2 = (uwb(i+1,k,jsls)+dtb*uwbt(i+1,k,jsls))           &
                     & - atm2%u(i+1,k,j)
                fls3 = (uwb(i,k,jsls-1)+dtb*uwbt(i,k,jsls-1))           &
                     & - atm2%u(i,k,j-1)
                fls4 = (uwb(i,k,jsls+1)+dtb*uwbt(i,k,jsls+1))           &
                     & - atm2%u(i,k,j+1)
#endif
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
          else if ( jj <= ip ) then
!-------east-boundary slice:
            fcx = fcoef*xfun(jsls)
            gcx = gcoef*xfun(jsls)
            do k = 1 , kz
              do i = ibeg , iend
#ifdef MPP1
                fls0 = (ueb(i,k,jew)+dtb*uebt(i,k,jew)) - atm2%u(i,k,j)
                fls1 = (ueb(i-1,k,jew)+dtb*uebt(i-1,k,jew))             &
                     & - atm2%u(i-1,k,j)
                fls2 = (ueb(i+1,k,jew)+dtb*uebt(i+1,k,jew))             &
                     & - atm2%u(i+1,k,j)
                fls3 = (ueb(i,k,jew-1)+dtb*uebt(i,k,jew-1))             &
                     & - atm2%u(i,k,j-1)
                fls4 = (ueb(i,k,jew+1)+dtb*uebt(i,k,jew+1))             &
                     & - atm2%u(i,k,j+1)
#else
                fls0 = (ueb(i,k,jsls)+dtb*uebt(i,k,jsls))-atm2%u(i,k,j)
                fls1 = (ueb(i-1,k,jsls)+dtb*uebt(i-1,k,jsls))           &
                     & - atm2%u(i-1,k,j)
                fls2 = (ueb(i+1,k,jsls)+dtb*uebt(i+1,k,jsls))           &
                     & - atm2%u(i+1,k,j)
                fls3 = (ueb(i,k,jsls-1)+dtb*uebt(i,k,jsls-1))           &
                     & - atm2%u(i,k,j-1)
                fls4 = (ueb(i,k,jsls+1)+dtb*uebt(i,k,jsls+1))           &
                     & - atm2%u(i,k,j+1)
#endif
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
          end if
        end if
!
      else if ( ibdy == 5 ) then
 
!----------use exponential method
 
        if ( jsls > ip ) then
!------interior j slices:
          do i = 2 , ip
            ii = iy - i + 1
            do k = 1 , kz
              fcx = fcoef*xfune(i,k)
              gcx = gcoef*xfune(i,k)
!........south boundary:
              fls0 = (usb(i,k,j)+dtb*usbt(i,k,j)) - atm2%u(i,k,j)
              fls1 = (usb(i,k,j-1)+dtb*usbt(i,k,j-1)) - atm2%u(i,k,j-1)
              fls2 = (usb(i,k,j+1)+dtb*usbt(i,k,j+1)) - atm2%u(i,k,j+1)
              fls3 = (usb(i-1,k,j)+dtb*usbt(i-1,k,j)) - atm2%u(i-1,k,j)
              fls4 = (usb(i+1,k,j)+dtb*usbt(i+1,k,j)) - atm2%u(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (unb(i,k,j)+dtb*unbt(i,k,j)) - atm2%u(ii,k,j)
              fls1 = (unb(i,k,j-1)+dtb*unbt(i,k,j-1)) - atm2%u(ii,k,j-1)
              fls2 = (unb(i,k,j+1)+dtb*unbt(i,k,j+1)) - atm2%u(ii,k,j+1)
              fls3 = (unb(i-1,k,j)+dtb*unbt(i-1,k,j)) - atm2%u(ii-1,k,j)
              fls4 = (unb(i+1,k,j)+dtb*unbt(i+1,k,j)) - atm2%u(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
!
        else if ( jsls <= ip ) then
!------east or west boundary slices:
          ibeg = 2
          iend = iym1
          if ( jsls > 2 ) then
            do i = 2 , jsls - 1
              ii = iy - i + 1
              do k = 1 , kz
                fcx = fcoef*xfune(i,k)
                gcx = gcoef*xfune(i,k)
!.........south boundary:
                fls0 = (usb(i,k,j)+dtb*usbt(i,k,j)) - atm2%u(i,k,j)
                fls1 = (usb(i,k,j-1)+dtb*usbt(i,k,j-1))-atm2%u(i,k,j-1)
                fls2 = (usb(i,k,j+1)+dtb*usbt(i,k,j+1))-atm2%u(i,k,j+1)
                fls3 = (usb(i-1,k,j)+dtb*usbt(i-1,k,j))-atm2%u(i-1,k,j)
                fls4 = (usb(i+1,k,j)+dtb*usbt(i+1,k,j))-atm2%u(i+1,k,j)
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!.........north boundary:
                fls0 = (unb(i,k,j)+dtb*unbt(i,k,j)) - atm2%u(ii,k,j)
                fls1 = (unb(i,k,j-1)+dtb*unbt(i,k,j-1))-atm2%u(ii,k,j-1)
                fls2 = (unb(i,k,j+1)+dtb*unbt(i,k,j+1))-atm2%u(ii,k,j+1)
                fls3 = (unb(i-1,k,j)+dtb*unbt(i-1,k,j))-atm2%u(ii-1,k,j)
                fls4 = (unb(i+1,k,j)+dtb*unbt(i+1,k,j))-atm2%u(ii+1,k,j)
                ften(ii,k) = ften(ii,k) + fcx*fls0 -                    &
                           & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
            ibeg = jsls
            iend = iy - jsls + 1
          end if
!
          if ( jj > ip ) then
!-------west-boundary slice:
            do k = 1 , kz
              fcx = fcoef*xfune(jsls,k)
              gcx = gcoef*xfune(jsls,k)
              do i = ibeg , iend
#ifdef MPP1
                fls0 = (uwb(i,k,jew)+dtb*uwbt(i,k,jew)) - atm2%u(i,k,j)
                fls1 = (uwb(i-1,k,jew)+dtb*uwbt(i-1,k,jew))             &
                     & - atm2%u(i-1,k,j)
                fls2 = (uwb(i+1,k,jew)+dtb*uwbt(i+1,k,jew))             &
                     & - atm2%u(i+1,k,j)
                fls3 = (uwb(i,k,jew-1)+dtb*uwbt(i,k,jew-1))             &
                     & - atm2%u(i,k,j-1)
                fls4 = (uwb(i,k,jew+1)+dtb*uwbt(i,k,jew+1))             &
                     & - atm2%u(i,k,j+1)
#else
                fls0 = (uwb(i,k,jsls)+dtb*uwbt(i,k,jsls))-atm2%u(i,k,j)
                fls1 = (uwb(i-1,k,jsls)+dtb*uwbt(i-1,k,jsls))           &
                     & - atm2%u(i-1,k,j)
                fls2 = (uwb(i+1,k,jsls)+dtb*uwbt(i+1,k,jsls))           &
                     & - atm2%u(i+1,k,j)
                fls3 = (uwb(i,k,jsls-1)+dtb*uwbt(i,k,jsls-1))           &
                     & - atm2%u(i,k,j-1)
                fls4 = (uwb(i,k,jsls+1)+dtb*uwbt(i,k,jsls+1))           &
                     & - atm2%u(i,k,j+1)
#endif
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
          else if ( jj <= ip ) then
!-------east-boundary slice:
            do k = 1 , kz
              fcx = fcoef*xfune(jsls,k)
              gcx = gcoef*xfune(jsls,k)
              do i = ibeg , iend
#ifdef MPP1
                fls0 = (ueb(i,k,jew)+dtb*uebt(i,k,jew)) - atm2%u(i,k,j)
                fls1 = (ueb(i-1,k,jew)+dtb*uebt(i-1,k,jew))             &
                     & - atm2%u(i-1,k,j)
                fls2 = (ueb(i+1,k,jew)+dtb*uebt(i+1,k,jew))             &
                     & - atm2%u(i+1,k,j)
                fls3 = (ueb(i,k,jew-1)+dtb*uebt(i,k,jew-1))             &
                     & - atm2%u(i,k,j-1)
                fls4 = (ueb(i,k,jew+1)+dtb*uebt(i,k,jew+1))             &
                     & - atm2%u(i,k,j+1)
#else
                fls0 = (ueb(i,k,jsls)+dtb*uebt(i,k,jsls))-atm2%u(i,k,j)
                fls1 = (ueb(i-1,k,jsls)+dtb*uebt(i-1,k,jsls))           &
                     & - atm2%u(i-1,k,j)
                fls2 = (ueb(i+1,k,jsls)+dtb*uebt(i+1,k,jsls))           &
                     & - atm2%u(i+1,k,j)
                fls3 = (ueb(i,k,jsls-1)+dtb*uebt(i,k,jsls-1))           &
                     & - atm2%u(i,k,j-1)
                fls4 = (ueb(i,k,jsls+1)+dtb*uebt(i,k,jsls+1))           &
                     & - atm2%u(i,k,j+1)
#endif
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
          end if
        end if
!
      end if
#endif
      call time_end(subroutine_name,idindx)
      end subroutine nudge_u
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine nudge_v(ip,fcoef,gcoef,xt,ften,j,ibdy)
!
      implicit none
!
      real(8) :: fcoef , gcoef , xt
      integer :: ibdy , ip , j
      real(8) , dimension(iy,kz) :: ften
      intent (in) fcoef , gcoef , ibdy , ip , j , xt
      intent (inout) ften
!
      real(8) :: dtb , fcx , fls0 , fls1 , fls2 , fls3 , fls4 , gcx
      integer :: i , ii , k
#ifndef BAND
      integer :: ibeg , iend , jj , jsls
#else
#ifndef MPP1
      integer :: jp1 , jm1
#endif
#endif
#ifdef MPP1
#ifndef BAND
      integer :: jew
#endif
#endif
      character (len=50) :: subroutine_name='nudge_v'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
#ifdef BAND
!----------------------------------------------------------------------
!
      dtb = xt*minph
!
!-----determine which relaxation method to use:linear/expon.
!
      if ( ibdy == 1 ) then
!
!---------use linear method
!
#ifdef MPP1
!------interior j slices:
         do i = 2 , ip
            ii = iy - i + 1
            fcx = fcoef*xfun(i)
            gcx = gcoef*xfun(i)
            do k = 1 , kz
!.......south boundary:
              fls0 = (vsb(i,k,j)+dtb*vsbt(i,k,j)) - atm2%v(i,k,j)
              fls1 = (vsb(i,k,j-1)+dtb*vsbt(i,k,j-1)) - atm2%v(i,k,j-1)
              fls2 = (vsb(i,k,j+1)+dtb*vsbt(i,k,j+1)) - atm2%v(i,k,j+1)
              fls3 = (vsb(i-1,k,j)+dtb*vsbt(i-1,k,j)) - atm2%v(i-1,k,j)
              fls4 = (vsb(i+1,k,j)+dtb*vsbt(i+1,k,j)) - atm2%v(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (vnb(i,k,j)+dtb*vnbt(i,k,j)) - atm2%v(ii,k,j)
              fls1 = (vnb(i,k,j-1)+dtb*vnbt(i,k,j-1)) - atm2%v(ii,k,j-1)
              fls2 = (vnb(i,k,j+1)+dtb*vnbt(i,k,j+1)) - atm2%v(ii,k,j+1)
              fls3 = (vnb(i-1,k,j)+dtb*vnbt(i-1,k,j)) - atm2%v(ii-1,k,j)
              fls4 = (vnb(i+1,k,j)+dtb*vnbt(i+1,k,j)) - atm2%v(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
         end do
#else
         jm1 = j-1
         jp1 = j+1
         if(jm1 == 0) jm1 = jx
         if(jp1 == jx+1) jp1 = 1
!------interior j slices:
         do i = 2 , ip
            ii = iy - i + 1
            fcx = fcoef*xfun(i)
            gcx = gcoef*xfun(i)
            do k = 1 , kz
!.......south boundary:
              fls0 = (vsb(i,k,j)+dtb*vsbt(i,k,j)) - atm2%v(i,k,j)
              fls1 = (vsb(i,k,jm1)+dtb*vsbt(i,k,jm1)) - atm2%v(i,k,jm1)
              fls2 = (vsb(i,k,jp1)+dtb*vsbt(i,k,jp1)) - atm2%v(i,k,jp1)
              fls3 = (vsb(i-1,k,j)+dtb*vsbt(i-1,k,j)) - atm2%v(i-1,k,j)
              fls4 = (vsb(i+1,k,j)+dtb*vsbt(i+1,k,j)) - atm2%v(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (vnb(i,k,j)+dtb*vnbt(i,k,j)) - atm2%v(ii,k,j)
              fls1 = (vnb(i,k,jm1)+dtb*vnbt(i,k,jm1)) - atm2%v(ii,k,jm1)
              fls2 = (vnb(i,k,jp1)+dtb*vnbt(i,k,jp1)) - atm2%v(ii,k,jp1)
              fls3 = (vnb(i-1,k,j)+dtb*vnbt(i-1,k,j)) - atm2%v(ii-1,k,j)
              fls4 = (vnb(i+1,k,j)+dtb*vnbt(i+1,k,j)) - atm2%v(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
         end do
#endif
      else if ( ibdy == 5 ) then
 
!----------use exponential method
#ifdef MPP1
!------interior j slices:
         do i = 2 , ip
            ii = iy - i + 1
            do k = 1 , kz
              fcx = fcoef*xfune(i,k)
              gcx = gcoef*xfune(i,k)
!........south boundary:
              fls0 = (vsb(i,k,j)+dtb*vsbt(i,k,j)) - atm2%v(i,k,j)
              fls1 = (vsb(i,k,j-1)+dtb*vsbt(i,k,j-1)) - atm2%v(i,k,j-1)
              fls2 = (vsb(i,k,j+1)+dtb*vsbt(i,k,j+1)) - atm2%v(i,k,j+1)
              fls3 = (vsb(i-1,k,j)+dtb*vsbt(i-1,k,j)) - atm2%v(i-1,k,j)
              fls4 = (vsb(i+1,k,j)+dtb*vsbt(i+1,k,j)) - atm2%v(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (vnb(i,k,j)+dtb*vnbt(i,k,j)) - atm2%v(ii,k,j)
              fls1 = (vnb(i,k,j-1)+dtb*vnbt(i,k,j-1)) - atm2%v(ii,k,j-1)
              fls2 = (vnb(i,k,j+1)+dtb*vnbt(i,k,j+1)) - atm2%v(ii,k,j+1)
              fls3 = (vnb(i-1,k,j)+dtb*vnbt(i-1,k,j)) - atm2%v(ii-1,k,j)
              fls4 = (vnb(i+1,k,j)+dtb*vnbt(i+1,k,j)) - atm2%v(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
         end do
#else
         jm1 = j-1
         jp1 = j+1
         if(jm1 == 0) jm1 = jx
         if(jp1 == jx+1) jp1 = 1
!------interior j slices:
         do i = 2 , ip
            ii = iy - i + 1
            do k = 1 , kz
              fcx = fcoef*xfune(i,k)
              gcx = gcoef*xfune(i,k)
!........south boundary:
              fls0 = (vsb(i,k,j)+dtb*vsbt(i,k,j)) - atm2%v(i,k,j)
              fls1 = (vsb(i,k,jm1)+dtb*vsbt(i,k,jm1)) - atm2%v(i,k,jm1)
              fls2 = (vsb(i,k,jp1)+dtb*vsbt(i,k,jp1)) - atm2%v(i,k,jp1)
              fls3 = (vsb(i-1,k,j)+dtb*vsbt(i-1,k,j)) - atm2%v(i-1,k,j)
              fls4 = (vsb(i+1,k,j)+dtb*vsbt(i+1,k,j)) - atm2%v(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (vnb(i,k,j)+dtb*vnbt(i,k,j)) - atm2%v(ii,k,j)
              fls1 = (vnb(i,k,jm1)+dtb*vnbt(i,k,jm1)) - atm2%v(ii,k,jm1)
              fls2 = (vnb(i,k,jp1)+dtb*vnbt(i,k,jp1)) - atm2%v(ii,k,jp1)
              fls3 = (vnb(i-1,k,j)+dtb*vnbt(i-1,k,j)) - atm2%v(ii-1,k,j)
              fls4 = (vnb(i+1,k,j)+dtb*vnbt(i+1,k,j)) - atm2%v(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
         end do
#endif
      end if
#else
!----------------------------------------------------------------------
!
      dtb = xt*minph
#ifdef MPP1
      jsls = j + myid*jxp
      jj = jxp1 - jsls
      if ( jj <= ip ) jsls = jj
      jew = jsls
      if ( jew > jxp ) jew = mod(jsls,jxp)
      if ( jew == 0 ) jew = jxp
#else
      jsls = j
      jj = jxp1 - jsls
      if ( jj <= ip ) jsls = jj
#endif
!
!-----determine which relaxation method to use:linear/expon.
!
      if ( ibdy == 1 ) then
!
!---------use linear method
!
        if ( jsls > ip ) then
!------interior j slices:
          do i = 2 , ip
            ii = iy - i + 1
            fcx = fcoef*xfun(i)
            gcx = gcoef*xfun(i)
            do k = 1 , kz
!.......south boundary:
              fls0 = (vsb(i,k,j)+dtb*vsbt(i,k,j)) - atm2%v(i,k,j)
              fls1 = (vsb(i,k,j-1)+dtb*vsbt(i,k,j-1)) - atm2%v(i,k,j-1)
              fls2 = (vsb(i,k,j+1)+dtb*vsbt(i,k,j+1)) - atm2%v(i,k,j+1)
              fls3 = (vsb(i-1,k,j)+dtb*vsbt(i-1,k,j)) - atm2%v(i-1,k,j)
              fls4 = (vsb(i+1,k,j)+dtb*vsbt(i+1,k,j)) - atm2%v(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (vnb(i,k,j)+dtb*vnbt(i,k,j)) - atm2%v(ii,k,j)
              fls1 = (vnb(i,k,j-1)+dtb*vnbt(i,k,j-1)) - atm2%v(ii,k,j-1)
              fls2 = (vnb(i,k,j+1)+dtb*vnbt(i,k,j+1)) - atm2%v(ii,k,j+1)
              fls3 = (vnb(i-1,k,j)+dtb*vnbt(i-1,k,j)) - atm2%v(ii-1,k,j)
              fls4 = (vnb(i+1,k,j)+dtb*vnbt(i+1,k,j)) - atm2%v(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
!
        else if ( jsls <= ip ) then
!------east or west boundary slices:
          ibeg = 2
          iend = iym1
          if ( jsls > 2 ) then
            do i = 2 , jsls - 1
              ii = iy - i + 1
              fcx = fcoef*xfun(i)
              gcx = gcoef*xfun(i)
              do k = 1 , kz
!........south  boundary:
                fls0 = (vsb(i,k,j)+dtb*vsbt(i,k,j)) - atm2%v(i,k,j)
                fls1 = (vsb(i,k,j-1)+dtb*vsbt(i,k,j-1))-atm2%v(i,k,j-1)
                fls2 = (vsb(i,k,j+1)+dtb*vsbt(i,k,j+1))-atm2%v(i,k,j+1)
                fls3 = (vsb(i-1,k,j)+dtb*vsbt(i-1,k,j))-atm2%v(i-1,k,j)
                fls4 = (vsb(i+1,k,j)+dtb*vsbt(i+1,k,j))-atm2%v(i+1,k,j)
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!.........north boundary:
                fls0 = (vnb(i,k,j)+dtb*vnbt(i,k,j)) - atm2%v(ii,k,j)
                fls1 = (vnb(i,k,j-1)+dtb*vnbt(i,k,j-1))-atm2%v(ii,k,j-1)
                fls2 = (vnb(i,k,j+1)+dtb*vnbt(i,k,j+1))-atm2%v(ii,k,j+1)
                fls3 = (vnb(i-1,k,j)+dtb*vnbt(i-1,k,j))-atm2%v(ii-1,k,j)
                fls4 = (vnb(i+1,k,j)+dtb*vnbt(i+1,k,j))-atm2%v(ii+1,k,j)
                ften(ii,k) = ften(ii,k) + fcx*fls0 -                    &
                           & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
            ibeg = jsls
            iend = iy - jsls + 1
          end if
!
          if ( jj > ip ) then
!-------west-boundary slice:
            fcx = fcoef*xfun(jsls)
            gcx = gcoef*xfun(jsls)
            do k = 1 , kz
              do i = ibeg , iend
#ifdef MPP1
                fls0 = (vwb(i,k,jew)+dtb*vwbt(i,k,jew)) - atm2%v(i,k,j)
                fls1 = (vwb(i-1,k,jew)+dtb*vwbt(i-1,k,jew))             &
                     & - atm2%v(i-1,k,j)
                fls2 = (vwb(i+1,k,jew)+dtb*vwbt(i+1,k,jew))             &
                     & - atm2%v(i+1,k,j)
                fls3 = (vwb(i,k,jew-1)+dtb*vwbt(i,k,jew-1))             &
                     & - atm2%v(i,k,j-1)
                fls4 = (vwb(i,k,jew+1)+dtb*vwbt(i,k,jew+1))             &
                     & - atm2%v(i,k,j+1)
#else
                fls0 = (vwb(i,k,jsls)+dtb*vwbt(i,k,jsls))-atm2%v(i,k,j)
                fls1 = (vwb(i-1,k,jsls)+dtb*vwbt(i-1,k,jsls))           &
                     & - atm2%v(i-1,k,j)
                fls2 = (vwb(i+1,k,jsls)+dtb*vwbt(i+1,k,jsls))           &
                     & - atm2%v(i+1,k,j)
                fls3 = (vwb(i,k,jsls-1)+dtb*vwbt(i,k,jsls-1))           &
                     & - atm2%v(i,k,j-1)
                fls4 = (vwb(i,k,jsls+1)+dtb*vwbt(i,k,jsls+1))           &
                     & - atm2%v(i,k,j+1)
#endif
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
          else if ( jj <= ip ) then
!-------east-boundary slice:
            fcx = fcoef*xfun(jsls)
            gcx = gcoef*xfun(jsls)
            do k = 1 , kz
              do i = ibeg , iend
#ifdef MPP1
                fls0 = (veb(i,k,jew)+dtb*vebt(i,k,jew)) - atm2%v(i,k,j)
                fls1 = (veb(i-1,k,jew)+dtb*vebt(i-1,k,jew))             &
                     & - atm2%v(i-1,k,j)
                fls2 = (veb(i+1,k,jew)+dtb*vebt(i+1,k,jew))             &
                     & - atm2%v(i+1,k,j)
                fls3 = (veb(i,k,jew-1)+dtb*vebt(i,k,jew-1))             &
                     & - atm2%v(i,k,j-1)
                fls4 = (veb(i,k,jew+1)+dtb*vebt(i,k,jew+1))             &
                     & - atm2%v(i,k,j+1)
#else
                fls0 = (veb(i,k,jsls)+dtb*vebt(i,k,jsls))-atm2%v(i,k,j)
                fls1 = (veb(i-1,k,jsls)+dtb*vebt(i-1,k,jsls))           &
                     & - atm2%v(i-1,k,j)
                fls2 = (veb(i+1,k,jsls)+dtb*vebt(i+1,k,jsls))           &
                     & - atm2%v(i+1,k,j)
                fls3 = (veb(i,k,jsls-1)+dtb*vebt(i,k,jsls-1))           &
                     & - atm2%v(i,k,j-1)
                fls4 = (veb(i,k,jsls+1)+dtb*vebt(i,k,jsls+1))           &
                     & - atm2%v(i,k,j+1)
#endif
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
          end if
        end if
!
      else if ( ibdy == 5 ) then
 
!----------use exponential method
 
        if ( jsls > ip ) then
!------interior j slices:
          do i = 2 , ip
            ii = iy - i + 1
            do k = 1 , kz
              fcx = fcoef*xfune(i,k)
              gcx = gcoef*xfune(i,k)
!........south boundary:
              fls0 = (vsb(i,k,j)+dtb*vsbt(i,k,j)) - atm2%v(i,k,j)
              fls1 = (vsb(i,k,j-1)+dtb*vsbt(i,k,j-1)) - atm2%v(i,k,j-1)
              fls2 = (vsb(i,k,j+1)+dtb*vsbt(i,k,j+1)) - atm2%v(i,k,j+1)
              fls3 = (vsb(i-1,k,j)+dtb*vsbt(i-1,k,j)) - atm2%v(i-1,k,j)
              fls4 = (vsb(i+1,k,j)+dtb*vsbt(i+1,k,j)) - atm2%v(i+1,k,j)
              ften(i,k) = ften(i,k) + fcx*fls0 -                        &
                        & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
              fls0 = (vnb(i,k,j)+dtb*vnbt(i,k,j)) - atm2%v(ii,k,j)
              fls1 = (vnb(i,k,j-1)+dtb*vnbt(i,k,j-1)) - atm2%v(ii,k,j-1)
              fls2 = (vnb(i,k,j+1)+dtb*vnbt(i,k,j+1)) - atm2%v(ii,k,j+1)
              fls3 = (vnb(i-1,k,j)+dtb*vnbt(i-1,k,j)) - atm2%v(ii-1,k,j)
              fls4 = (vnb(i+1,k,j)+dtb*vnbt(i+1,k,j)) - atm2%v(ii+1,k,j)
              ften(ii,k) = ften(ii,k) + fcx*fls0 -                      &
                         & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
!
        else if ( jsls <= ip ) then
!------east or west boundary slices:
          ibeg = 2
          iend = iym1
          if ( jsls > 2 ) then
            do i = 2 , jsls - 1
              ii = iy - i + 1
              do k = 1 , kz
                fcx = fcoef*xfune(i,k)
                gcx = gcoef*xfune(i,k)
!.........south boundary:
                fls0 = (vsb(i,k,j)+dtb*vsbt(i,k,j)) - atm2%v(i,k,j)
                fls1 = (vsb(i,k,j-1)+dtb*vsbt(i,k,j-1))-atm2%v(i,k,j-1)
                fls2 = (vsb(i,k,j+1)+dtb*vsbt(i,k,j+1))-atm2%v(i,k,j+1)
                fls3 = (vsb(i-1,k,j)+dtb*vsbt(i-1,k,j))-atm2%v(i-1,k,j)
                fls4 = (vsb(i+1,k,j)+dtb*vsbt(i+1,k,j))-atm2%v(i+1,k,j)
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
!.........north boundary:
                fls0 = (vnb(i,k,j)+dtb*vnbt(i,k,j)) - atm2%v(ii,k,j)
                fls1 = (vnb(i,k,j-1)+dtb*vnbt(i,k,j-1))-atm2%v(ii,k,j-1)
                fls2 = (vnb(i,k,j+1)+dtb*vnbt(i,k,j+1))-atm2%v(ii,k,j+1)
                fls3 = (vnb(i-1,k,j)+dtb*vnbt(i-1,k,j))-atm2%v(ii-1,k,j)
                fls4 = (vnb(i+1,k,j)+dtb*vnbt(i+1,k,j))-atm2%v(ii+1,k,j)
                ften(ii,k) = ften(ii,k) + fcx*fls0 -                    &
                           & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
            ibeg = jsls
            iend = iy - jsls + 1
          end if
!
          if ( jj > ip ) then
!-------west-boundary slice:
            do k = 1 , kz
              fcx = fcoef*xfune(jsls,k)
              gcx = gcoef*xfune(jsls,k)
              do i = ibeg , iend
#ifdef MPP1
                fls0 = (vwb(i,k,jew)+dtb*vwbt(i,k,jew)) - atm2%v(i,k,j)
                fls1 = (vwb(i-1,k,jew)+dtb*vwbt(i-1,k,jew))             &
                     & - atm2%v(i-1,k,j)
                fls2 = (vwb(i+1,k,jew)+dtb*vwbt(i+1,k,jew))             &
                     & - atm2%v(i+1,k,j)
                fls3 = (vwb(i,k,jew-1)+dtb*vwbt(i,k,jew-1))             &
                     & - atm2%v(i,k,j-1)
                fls4 = (vwb(i,k,jew+1)+dtb*vwbt(i,k,jew+1))             &
                     & - atm2%v(i,k,j+1)
#else
                fls0 = (vwb(i,k,jsls)+dtb*vwbt(i,k,jsls))-atm2%v(i,k,j)
                fls1 = (vwb(i-1,k,jsls)+dtb*vwbt(i-1,k,jsls))           &
                     & - atm2%v(i-1,k,j)
                fls2 = (vwb(i+1,k,jsls)+dtb*vwbt(i+1,k,jsls))           &
                     & - atm2%v(i+1,k,j)
                fls3 = (vwb(i,k,jsls-1)+dtb*vwbt(i,k,jsls-1))           &
                     & - atm2%v(i,k,j-1)
                fls4 = (vwb(i,k,jsls+1)+dtb*vwbt(i,k,jsls+1))           &
                     & - atm2%v(i,k,j+1)
#endif
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
          else if ( jj <= ip ) then
!-------east-boundary slice:
            do k = 1 , kz
              fcx = fcoef*xfune(jsls,k)
              gcx = gcoef*xfune(jsls,k)
              do i = ibeg , iend
#ifdef MPP1
                fls0 = (veb(i,k,jew)+dtb*vebt(i,k,jew)) - atm2%v(i,k,j)
                fls1 = (veb(i-1,k,jew)+dtb*vebt(i-1,k,jew))             &
                     & - atm2%v(i-1,k,j)
                fls2 = (veb(i+1,k,jew)+dtb*vebt(i+1,k,jew))             &
                     & - atm2%v(i+1,k,j)
                fls3 = (veb(i,k,jew-1)+dtb*vebt(i,k,jew-1))             &
                     & - atm2%v(i,k,j-1)
                fls4 = (veb(i,k,jew+1)+dtb*vebt(i,k,jew+1))             &
                     & - atm2%v(i,k,j+1)
#else
                fls0 = (veb(i,k,jsls)+dtb*vebt(i,k,jsls))-atm2%v(i,k,j)
                fls1 = (veb(i-1,k,jsls)+dtb*vebt(i-1,k,jsls))           &
                     & - atm2%v(i-1,k,j)
                fls2 = (veb(i+1,k,jsls)+dtb*vebt(i+1,k,jsls))           &
                     & - atm2%v(i+1,k,j)
                fls3 = (veb(i,k,jsls-1)+dtb*vebt(i,k,jsls-1))           &
                     & - atm2%v(i,k,j-1)
                fls4 = (veb(i,k,jsls+1)+dtb*vebt(i,k,jsls+1))           &
                     & - atm2%v(i,k,j+1)
#endif
                ften(i,k) = ften(i,k) + fcx*fls0 -                      &
                          & gcx*c203*(fls1+fls2+fls3+fls4-d_four*fls0)
              end do
            end do
          end if
        end if
!
      end if
#endif
      call time_end(subroutine_name,idindx)
      end subroutine nudge_v
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine applies sponge boundary condition to the        c
!     tendency term - ften.                                           c
!                                                                     c
!     ip   : is the number of slices affected by sponge boundary.     c
!                                                                     c
!     wg   : are the weightings.                                      c
!                                                                     c
!     ften : is the tendency calculated from the model.               c
!                                                                     c
!     pebt, pwbt, pnbt, psbt : are the large-scale or observed        c
!            tendencies at east, west, north, and south boundaries.   c
!                                                                     c
!     ie = iy, je = jx for dot-point variables.                       c
!     ie = iym1, je = jxm1 for cross-point variables.                 c
!                                                                     c
!     j    : is the j'th slice of the tendency to be adjusted.        c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine sponge_p(ip,wg,ften,j)
!
      implicit none
!
      integer :: ip , j
      real(8) , dimension(iy) :: ften
      real(8) , dimension(ip) :: wg
      intent (in) ip , j , wg
      intent (inout) ften
!
      integer :: i , ii
#ifndef BAND
      integer :: ibeg , iend , jj , jsls
#endif
#ifdef MPP1
#ifndef BAND
      integer :: jwb , jeb
#endif
#endif
!
      character (len=50) :: subroutine_name='sponge_p'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
#ifdef BAND
!----------------------------------------------------------------------
!
!-----interior j slices:
      do i = 2 , ip
         ii = iy - i
!.......south boundary:
         ften(i) = wg(i)*ften(i) + (d_one-wg(i))*psbt(i,j)
!.......north boundary:
         ften(ii) = wg(i)*ften(ii) + (d_one-wg(i))*pnbt(i,j)
      end do

#else
!----------------------------------------------------------------------
!
#ifdef MPP1
      jsls = j + myid*jxp
      jj = jx - jsls
      if ( jj <= ip ) jsls = jj
      jwb = jsls
      if ( jwb > jxp ) jwb = mod(jwb,jxp)
      if ( jwb == 0 ) jwb = jxp
      if ( myid == nproc-1 ) then
        jeb = jsls
      else
        jeb = jsls + 1
      end if
      if ( jeb > jxp ) jeb = mod(jeb,jxp)
      if ( jeb == 0 ) jeb = jxp
#else
      jsls = j
      jj = jx - jsls
      if ( jj <= ip ) jsls = jj
#endif
!
      if ( jsls > ip ) then
!-----interior j slices:
        do i = 2 , ip
          ii = iy - i
!.......south boundary:
          ften(i) = wg(i)*ften(i) + (d_one-wg(i))*psbt(i,j)
!.......north boundary:
          ften(ii) = wg(i)*ften(ii) + (d_one-wg(i))*pnbt(i,j)
        end do
!
      else if ( jsls <= ip ) then
        ibeg = 2
        iend = iym1 - 1
        if ( jsls > 2 ) then
          do i = 2 , jsls - 1
            ii = iy - i
!........south boundary:
            ften(i) = wg(i)*ften(i) + (d_one-wg(i))*psbt(i,j)
!........north boundary:
            ften(ii) = wg(i)*ften(ii) + (d_one-wg(i))*pnbt(i,j)
          end do
          ibeg = jsls
          iend = iy - jsls
        end if
!
        if ( jj > ip ) then
!------west-boundary slice:
          do i = ibeg , iend
#ifdef MPP1
            if ( jsls <= ip ) ften(i) = wg(jsls)*ften(i) + (d_one-wg(jsls))&
                                      & *pwbt(i,jwb)
#else
            ften(i) = wg(jsls)*ften(i) + (d_one-wg(jsls))*pwbt(i,jsls)
#endif
          end do
        else if ( jj <= ip ) then
!------east-boundary slice:
          do i = ibeg , iend
#ifdef MPP1
            if ( jsls <= ip ) ften(i) = wg(jsls)*ften(i) + (d_one-wg(jsls))&
                                      & *pebt(i,jeb)
#else
            ften(i) = wg(jsls)*ften(i) + (d_one-wg(jsls))*pebt(i,jsls)
#endif
          end do
        end if
!
      end if

#endif
      call time_end(subroutine_name,idindx)
      end subroutine sponge_p
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine sponge_t(ip,wg,ften,j)
!
      implicit none
!
      integer :: ip , j
      real(8) , dimension(iy,kz) :: ften
      real(8) , dimension(ip) :: wg
      intent (in) ip , j , wg
      intent (inout) ften
!
      integer :: i , ii , k
#ifndef BAND
      integer :: ibeg , iend , jj , jsls
#endif
#ifdef MPP1
#ifndef BAND
      integer :: jwb , jeb
#endif
#endif
      character (len=50) :: subroutine_name='sponge_t'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
!
#ifdef BAND
!----------------------------------------------------------------------
!
!-----interior j slices:
      do i = 2 , ip
         ii = iy - i
         do k = 1 , kz
!.......south boundary:
            ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*tsbt(i,k,j)
!.......north boundary:
            ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*tnbt(i,k,j)
         end do
      end do

#else
!----------------------------------------------------------------------
!
#ifdef MPP1
      jsls = j + myid*jxp
      jj = jx - jsls
      if ( jj <= ip ) jsls = jj
      jwb = jsls
      if ( jwb > jxp ) jwb = mod(jwb,jxp)
      if ( jwb == 0 ) jwb = jxp
      if ( myid == nproc-1 ) then
        jeb = jsls
      else
        jeb = jsls + 1
      end if
      if ( jeb > jxp ) jeb = mod(jeb,jxp)
      if ( jeb == 0 ) jeb = jxp
#else
      jsls = j
      jj = jx - jsls
      if ( jj <= ip ) jsls = jj
#endif
!
      if ( jsls > ip ) then
!-----interior j slices:
        do i = 2 , ip
          ii = iy - i
          do k = 1 , kz
!.......south boundary:
            ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*tsbt(i,k,j)
!.......north boundary:
            ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*tnbt(i,k,j)
          end do
        end do
!
      else if ( jsls <= ip ) then
        ibeg = 2
        iend = iym1 - 1
        if ( jsls > 2 ) then
          do i = 2 , jsls - 1
            ii = iy - i
            do k = 1 , kz
!........south boundary:
              ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*tsbt(i,k,j)
!........north boundary:
              ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*tnbt(i,k,j)
            end do
          end do
          ibeg = jsls
          iend = iy - jsls
        end if
!
        if ( jj > ip ) then
!------west-boundary slice:
          do k = 1 , kz
            do i = ibeg , iend
#ifdef MPP1
              if ( jsls <= ip ) ften(i,k) = wg(jsls)*ften(i,k)          &
                 & + (d_one-wg(jsls))*twbt(i,k,jwb)
#else
              ften(i,k) = wg(jsls)*ften(i,k) + (d_one-wg(jsls))            &
                        & *twbt(i,k,jsls)
#endif
            end do
          end do
        else if ( jj <= ip ) then
!------east-boundary slice:
          do k = 1 , kz
            do i = ibeg , iend
#ifdef MPP1
              if ( jsls <= ip ) ften(i,k) = wg(jsls)*ften(i,k)          &
                 & + (d_one-wg(jsls))*tebt(i,k,jeb)
#else
              ften(i,k) = wg(jsls)*ften(i,k) + (d_one-wg(jsls))            &
                        & *tebt(i,k,jsls)
#endif
            end do
          end do
        end if
!
      end if
#endif

      call time_end(subroutine_name,idindx)
      end subroutine sponge_t
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine spongeqv(ip,wg,ften,j)
!
      implicit none
!
      integer :: ip , j
      real(8) , dimension(iy,kz) :: ften
      real(8) , dimension(ip) :: wg
      intent (in) ip , j , wg
      intent (inout) ften
!
      integer :: i , ii , k
#ifndef BAND
      integer :: ibeg , iend , jj , jsls
#endif
#ifdef MPP1
#ifndef BAND
      integer :: jwb , jeb
#endif
#endif
      character (len=50) :: subroutine_name='spongeqv'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
!
#ifdef BAND
!----------------------------------------------------------------------
!
!-----interior j slices:
      do i = 2 , ip
         ii = iy - i
         do k = 1 , kz
!.......south boundary:
            ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*qsbt(i,k,j)
!.......north boundary:
            ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*qnbt(i,k,j)
         end do
      end do
#else
!----------------------------------------------------------------------
!
#ifdef MPP1
      jsls = j + myid*jxp
      jj = jx - jsls
      if ( jj <= ip ) jsls = jj
      jwb = jsls
      if ( jwb > jxp ) jwb = mod(jwb,jxp)
      if ( jwb == 0 ) jwb = jxp
      if ( myid == nproc-1 ) then
        jeb = jsls
      else
        jeb = jsls + 1
      end if
      if ( jeb > jxp ) jeb = mod(jeb,jxp)
      if ( jeb == 0 ) jeb = jxp
#else
      jsls = j
      jj = jx - jsls
      if ( jj <= ip ) jsls = jj
#endif
!
      if ( jsls > ip ) then
!-----interior j slices:
        do i = 2 , ip
          ii = iy - i
          do k = 1 , kz
!.......south boundary:
            ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*qsbt(i,k,j)
!.......north boundary:
            ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*qnbt(i,k,j)
          end do
        end do
!
      else if ( jsls <= ip ) then
        ibeg = 2
        iend = iym1 - 1
        if ( jsls > 2 ) then
          do i = 2 , jsls - 1
            ii = iy - i
            do k = 1 , kz
!........south boundary:
              ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*qsbt(i,k,j)
!........north boundary:
              ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*qnbt(i,k,j)
            end do
          end do
          ibeg = jsls
          iend = iy - jsls
        end if
!
        if ( jj > ip ) then
!------west-boundary slice:
          do k = 1 , kz
            do i = ibeg , iend
#ifdef MPP1
              if ( jsls <= ip ) ften(i,k) = wg(jsls)*ften(i,k)          &
                 & + (d_one-wg(jsls))*qwbt(i,k,jwb)
#else
              ften(i,k) = wg(jsls)*ften(i,k) + (d_one-wg(jsls))            &
                        & *qwbt(i,k,jsls)
#endif
            end do
          end do
        else if ( jj <= ip ) then
!------east-boundary slice:
          do k = 1 , kz
            do i = ibeg , iend
#ifdef MPP1
              if ( jsls <= ip ) ften(i,k) = wg(jsls)*ften(i,k)          &
                 & + (d_one-wg(jsls))*qebt(i,k,jeb)
#else
              ften(i,k) = wg(jsls)*ften(i,k) + (d_one-wg(jsls))            &
                        & *qebt(i,k,jsls)
#endif
            end do
          end do
        end if
!
      end if
#endif
      call time_end(subroutine_name,idindx)
      end subroutine spongeqv
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine sponge_u(ip,wg,ften,j)
!
      implicit none
!
      integer :: ip , j
      real(8) , dimension(iy,kz) :: ften
      real(8) , dimension(ip) :: wg
      intent (in) ip , j , wg
      intent (inout) ften
!
      integer :: i , ii , k
#ifndef BAND
      integer :: ibeg , iend , jj , jsls
#endif
#ifdef MPP1
#ifndef BAND
      integer :: jew
#endif
#endif
      character (len=50) :: subroutine_name='sponge_u'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
!
#ifdef BAND
!----------------------------------------------------------------------
!
!-----interior j slices:
      do i = 2 , ip
         ii = iy - i + 1
         do k = 1 , kz
!.......south boundary:
            ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*usbt(i,k,j)
!.......north boundary:
            ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*unbt(i,k,j)
         end do
      end do
#else
!----------------------------------------------------------------------
!
#ifdef MPP1
      jsls = j + myid*jxp
      jj = jxp1 - jsls
      if ( jj <= ip ) jsls = jj
      jew = jsls
      if ( jew > jxp ) jew = mod(jsls,jxp)
      if ( jew == 0 ) jew = jxp
#else
      jsls = j
      jj = jxp1 - jsls
      if ( jj <= ip ) jsls = jj
#endif
!
      if ( jsls > ip ) then
!-----interior j slices:
        do i = 2 , ip
          ii = iy - i + 1
          do k = 1 , kz
!.......south boundary:
            ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*usbt(i,k,j)
!.......north boundary:
            ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*unbt(i,k,j)
          end do
        end do
!
      else if ( jsls <= ip ) then
        ibeg = 2
        iend = iym1
        if ( jsls > 2 ) then
          do i = 2 , jsls - 1
            ii = iy - i + 1
            do k = 1 , kz
!........south boundary:
              ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*usbt(i,k,j)
!........north boundary:
              ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*unbt(i,k,j)
            end do
          end do
          ibeg = jsls
          iend = iy - jsls + 1
        end if
!
        if ( jj > ip ) then
!------west-boundary slice:
          do k = 1 , kz
            do i = ibeg , iend
#ifdef MPP1
              if ( jsls <= ip ) ften(i,k) = wg(jsls)*ften(i,k)          &
                 & + (d_one-wg(jsls))*uwbt(i,k,jew)
#else
              ften(i,k) = wg(jsls)*ften(i,k) + (d_one-wg(jsls))            &
                        & *uwbt(i,k,jsls)
#endif
            end do
          end do
        else if ( jj <= ip ) then
!------east-boundary slice:
          do k = 1 , kz
            do i = ibeg , iend
#ifdef MPP1
              if ( jsls <= ip ) ften(i,k) = wg(jsls)*ften(i,k)          &
                 & + (d_one-wg(jsls))*uebt(i,k,jew)
#else
              ften(i,k) = wg(jsls)*ften(i,k) + (d_one-wg(jsls))            &
                        & *uebt(i,k,jsls)
#endif
            end do
          end do
        end if
!
      end if
#endif
      call time_end(subroutine_name,idindx)
      end subroutine sponge_u
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine sponge_v(ip,wg,ften,j)
!
      implicit none
!
      integer :: ip , j
      real(8) , dimension(iy,kz) :: ften
      real(8) , dimension(ip) :: wg
      intent (in) ip , j , wg
      intent (inout) ften
!
      integer :: i , ii , k
#ifndef BAND
      integer :: ibeg , iend , jj , jsls
#endif
#ifdef MPP1
#ifndef BAND
      integer :: jew
#endif
#endif
      character (len=50) :: subroutine_name='sponge_v'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
!
#ifdef BAND
!----------------------------------------------------------------------
!
!-----interior j slices:
      do i = 2 , ip
         ii = iy - i + 1
         do k = 1 , kz
!.......south boundary:
            ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*vsbt(i,k,j)
!.......north boundary:
            ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*vnbt(i,k,j)
         end do
      end do
#else
!----------------------------------------------------------------------
!
#ifdef MPP1
      jsls = j + myid*jxp
      jj = jxp1 - jsls
      if ( jj <= ip ) jsls = jj
      jew = jsls
      if ( jew > jxp ) jew = mod(jsls,jxp)
      if ( jew == 0 ) jew = jxp
#else
      jsls = j
      jj = jxp1 - jsls
      if ( jj <= ip ) jsls = jj
#endif
!
      if ( jsls > ip ) then
!-----interior j slices:
        do i = 2 , ip
          ii = iy - i + 1
          do k = 1 , kz
!.......south boundary:
            ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*vsbt(i,k,j)
!.......north boundary:
            ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*vnbt(i,k,j)
          end do
        end do
!
      else if ( jsls <= ip ) then
        ibeg = 2
        iend = iym1
        if ( jsls > 2 ) then
          do i = 2 , jsls - 1
            ii = iy - i + 1
            do k = 1 , kz
!........south boundary:
              ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*vsbt(i,k,j)
!........north boundary:
              ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*vnbt(i,k,j)
            end do
          end do
          ibeg = jsls
          iend = iy - jsls + 1
        end if
!
        if ( jj > ip ) then
!------west-boundary slice:
          do k = 1 , kz
            do i = ibeg , iend
#ifdef MPP1
              if ( jsls <= ip ) ften(i,k) = wg(jsls)*ften(i,k)          &
                 & + (d_one-wg(jsls))*vwbt(i,k,jew)
#else
              ften(i,k) = wg(jsls)*ften(i,k) + (d_one-wg(jsls))            &
                        & *vwbt(i,k,jsls)
#endif
            end do
          end do
        else if ( jj <= ip ) then
!------east-boundary slice:
          do k = 1 , kz
            do i = ibeg , iend
#ifdef MPP1
              if ( jsls <= ip ) ften(i,k) = wg(jsls)*ften(i,k)          &
                 & + (d_one-wg(jsls))*vebt(i,k,jew)
#else
              ften(i,k) = wg(jsls)*ften(i,k) + (d_one-wg(jsls))            &
                        & *vebt(i,k,jsls)
#endif
            end do
          end do
        end if
!
      end if
#endif

      call time_end(subroutine_name,idindx)
      end subroutine sponge_v
!
      function xfun(mm)
        implicit none
        real(8) :: xfun
        integer , intent(in) :: mm
        xfun = dble(nspgd-mm)/(dble(nspgd)-d_two)
      end function xfun
      function xfune(mm,kk)
        implicit none
        real(8) :: xfune
        integer , intent(in) :: mm , kk
        xfune = dexp(-dble(mm-2)/anudg(kk))
        end function xfune
!
      end module mod_nudge
