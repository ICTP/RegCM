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
  private
!
  public :: sponge_p , sponge_t , spongeqv , sponge_u , sponge_v
  public :: nudge_p
  public :: nudge
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
!     xt    : is the time in seconds for variable "psb".              c
!                                                                     c
!     fcoef : are the coefficients for the newtonian term.            c
!                                                                     c
!     gcoef : are the coefficients for the diffusion term.            c
!                                                                     c
!     xpten : is the tendency calculated from the model.              c
!                                                                     c
!     peb, pwb, psb, pnb : are the observed boundary values           c
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
  use mod_bdycod
  use mod_atm_interface
  use mod_runparams
  use mod_service
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
  integer :: ibeg , iend , jj , jsls , jwb , jeb
#endif
  character (len=64) :: subroutine_name='nudge_p'
  integer :: idindx=0
!
  call time_begin(subroutine_name,idindx)
!
!----------------------------------------------------------------------
!
  dtb = xt
#ifdef BAND
!
!-----determine which relaxation method to use:linear/expon.
!
  if ( ibdy == 1 ) then
!
!---------use linear method
!
!------interior j slices:
     do i = 2 , ip
        ii = iym1 - i + 1
        fcx = fcoef*xfun(i)
        gcx = gcoef*xfun(i)
!.......south boundary:
        fls0 = (psb(i,j)+dtb*psbt(i,j)) - sps2%ps(i,j)
        fls1 = (psb(i,j-1)+dtb*psbt(i,j-1)) - sps2%ps(i,j-1)
        fls2 = (psb(i,j+1)+dtb*psbt(i,j+1)) - sps2%ps(i,j+1)
        fls3 = (psb(i-1,j)+dtb*psbt(i-1,j)) - sps2%ps(i-1,j)
        fls4 = (psb(i+1,j)+dtb*psbt(i+1,j)) - sps2%ps(i+1,j)
        xpten(i) = xpten(i) + fcx*fls0 -                            &
                & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
        fls0 = (pnb(i,j)+dtb*pnbt(i,j)) - sps2%ps(ii,j)
        fls1 = (pnb(i,j-1)+dtb*pnbt(i,j-1)) - sps2%ps(ii,j-1)
        fls2 = (pnb(i,j+1)+dtb*pnbt(i,j+1)) - sps2%ps(ii,j+1)
        fls3 = (pnb(i-1,j)+dtb*pnbt(i-1,j)) - sps2%ps(ii-1,j)
        fls4 = (pnb(i+1,j)+dtb*pnbt(i+1,j)) - sps2%ps(ii+1,j)
        xpten(ii) = xpten(ii) + fcx*fls0 -                          &
                 & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
     end do
!
  else if ( ibdy == 5 ) then
 
!----------use exponential method
 
!------interior j slices:
     do i = 2 , ip
        ii = iym1 - i + 1
        fcx = fcoef*xfune(i,kz)
        gcx = gcoef*xfune(i,kz)
!........south boundary:
        fls0 = (psb(i,j)+dtb*psbt(i,j)) - sps2%ps(i,j)
        fls1 = (psb(i,j-1)+dtb*psbt(i,j-1)) - sps2%ps(i,j-1)
        fls2 = (psb(i,j+1)+dtb*psbt(i,j+1)) - sps2%ps(i,j+1)
        fls3 = (psb(i-1,j)+dtb*psbt(i-1,j)) - sps2%ps(i-1,j)
        fls4 = (psb(i+1,j)+dtb*psbt(i+1,j)) - sps2%ps(i+1,j)
        xpten(i) = xpten(i) + fcx*fls0 -                            &
                & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
        fls0 = (pnb(i,j)+dtb*pnbt(i,j)) - sps2%ps(ii,j)
        fls1 = (pnb(i,j-1)+dtb*pnbt(i,j-1)) - sps2%ps(ii,j-1)
        fls2 = (pnb(i,j+1)+dtb*pnbt(i,j+1)) - sps2%ps(ii,j+1)
        fls3 = (pnb(i-1,j)+dtb*pnbt(i-1,j)) - sps2%ps(ii-1,j)
        fls4 = (pnb(i+1,j)+dtb*pnbt(i+1,j)) - sps2%ps(ii+1,j)
        xpten(ii) = xpten(ii) + fcx*fls0 -                          &
                 & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
     end do
  end if
#else
!----------------------------------------------------------------------
!
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
        fls0 = (psb(i,j)+dtb*psbt(i,j)) - sps2%ps(i,j)
        fls1 = (psb(i,j-1)+dtb*psbt(i,j-1)) - sps2%ps(i,j-1)
        fls2 = (psb(i,j+1)+dtb*psbt(i,j+1)) - sps2%ps(i,j+1)
        fls3 = (psb(i-1,j)+dtb*psbt(i-1,j)) - sps2%ps(i-1,j)
        fls4 = (psb(i+1,j)+dtb*psbt(i+1,j)) - sps2%ps(i+1,j)
        xpten(i) = xpten(i) + fcx*fls0 -                            &
                & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
        fls0 = (pnb(i,j)+dtb*pnbt(i,j)) - sps2%ps(ii,j)
        fls1 = (pnb(i,j-1)+dtb*pnbt(i,j-1)) - sps2%ps(ii,j-1)
        fls2 = (pnb(i,j+1)+dtb*pnbt(i,j+1)) - sps2%ps(ii,j+1)
        fls3 = (pnb(i-1,j)+dtb*pnbt(i-1,j)) - sps2%ps(ii-1,j)
        fls4 = (pnb(i+1,j)+dtb*pnbt(i+1,j)) - sps2%ps(ii+1,j)
        xpten(ii) = xpten(ii) + fcx*fls0 -                          &
                 & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
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
          fls0 = (psb(i,j)+dtb*psbt(i,j)) - sps2%ps(i,j)
          fls1 = (psb(i,j-1)+dtb*psbt(i,j-1)) - sps2%ps(i,j-1)
          fls2 = (psb(i,j+1)+dtb*psbt(i,j+1)) - sps2%ps(i,j+1)
          fls3 = (psb(i-1,j)+dtb*psbt(i-1,j)) - sps2%ps(i-1,j)
          fls4 = (psb(i+1,j)+dtb*psbt(i+1,j)) - sps2%ps(i+1,j)
          xpten(i) = xpten(i) + fcx*fls0 -                          &
                  & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
!.........north boundary:
          fls0 = (pnb(i,j)+dtb*pnbt(i,j)) - sps2%ps(ii,j)
          fls1 = (pnb(i,j-1)+dtb*pnbt(i,j-1)) - sps2%ps(ii,j-1)
          fls2 = (pnb(i,j+1)+dtb*pnbt(i,j+1)) - sps2%ps(ii,j+1)
          fls3 = (pnb(i-1,j)+dtb*pnbt(i-1,j)) - sps2%ps(ii-1,j)
          fls4 = (pnb(i+1,j)+dtb*pnbt(i+1,j)) - sps2%ps(ii+1,j)
          xpten(ii) = xpten(ii) + fcx*fls0 -                        &
                   & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
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
          fls0 = (pwb(i,jwb)+dtb*pwbt(i,jwb)) - sps2%ps(i,j)
          fls1 = (pwb(i-1,jwb)+dtb*pwbt(i-1,jwb)) - sps2%ps(i-1,j)
          fls2 = (pwb(i+1,jwb)+dtb*pwbt(i+1,jwb)) - sps2%ps(i+1,j)
          fls3 = (pwb(i,jwb-1)+dtb*pwbt(i,jwb-1)) - sps2%ps(i,j-1)
          fls4 = (pwb(i,jwb+1)+dtb*pwbt(i,jwb+1)) - sps2%ps(i,j+1)
          xpten(i) = xpten(i) + fcx*fls0 -                          &
                  & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      else if ( jj <= ip ) then
!-------east-boundary slice:
        fcx = fcoef*xfun(jsls)
        gcx = gcoef*xfun(jsls)
        do i = ibeg , iend
          fls0 = (peb(i,jeb)+dtb*pebt(i,jeb)) - sps2%ps(i,j)
          fls1 = (peb(i-1,jeb)+dtb*pebt(i-1,jeb)) - sps2%ps(i-1,j)
          fls2 = (peb(i+1,jeb)+dtb*pebt(i+1,jeb)) - sps2%ps(i+1,j)
          fls3 = (peb(i,jeb-1)+dtb*pebt(i,jeb-1)) - sps2%ps(i,j-1)
          fls4 = (peb(i,jeb+1)+dtb*pebt(i,jeb+1)) - sps2%ps(i,j+1)
          xpten(i) = xpten(i) + fcx*fls0 -                          &
                  & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
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
        fls0 = (psb(i,j)+dtb*psbt(i,j)) - sps2%ps(i,j)
        fls1 = (psb(i,j-1)+dtb*psbt(i,j-1)) - sps2%ps(i,j-1)
        fls2 = (psb(i,j+1)+dtb*psbt(i,j+1)) - sps2%ps(i,j+1)
        fls3 = (psb(i-1,j)+dtb*psbt(i-1,j)) - sps2%ps(i-1,j)
        fls4 = (psb(i+1,j)+dtb*psbt(i+1,j)) - sps2%ps(i+1,j)
        xpten(i) = xpten(i) + fcx*fls0 -                            &
                & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
        fls0 = (pnb(i,j)+dtb*pnbt(i,j)) - sps2%ps(ii,j)
        fls1 = (pnb(i,j-1)+dtb*pnbt(i,j-1)) - sps2%ps(ii,j-1)
        fls2 = (pnb(i,j+1)+dtb*pnbt(i,j+1)) - sps2%ps(ii,j+1)
        fls3 = (pnb(i-1,j)+dtb*pnbt(i-1,j)) - sps2%ps(ii-1,j)
        fls4 = (pnb(i+1,j)+dtb*pnbt(i+1,j)) - sps2%ps(ii+1,j)
        xpten(ii) = xpten(ii) + fcx*fls0 -                          &
                 & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
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
          fls0 = (psb(i,j)+dtb*psbt(i,j)) - sps2%ps(i,j)
          fls1 = (psb(i,j-1)+dtb*psbt(i,j-1)) - sps2%ps(i,j-1)
          fls2 = (psb(i,j+1)+dtb*psbt(i,j+1)) - sps2%ps(i,j+1)
          fls3 = (psb(i-1,j)+dtb*psbt(i-1,j)) - sps2%ps(i-1,j)
          fls4 = (psb(i+1,j)+dtb*psbt(i+1,j)) - sps2%ps(i+1,j)
          xpten(i) = xpten(i) + fcx*fls0 -                          &
                  & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
!.........north boundary:
          fls0 = (pnb(i,j)+dtb*pnbt(i,j)) - sps2%ps(ii,j)
          fls1 = (pnb(i,j-1)+dtb*pnbt(i,j-1)) - sps2%ps(ii,j-1)
          fls2 = (pnb(i,j+1)+dtb*pnbt(i,j+1)) - sps2%ps(ii,j+1)
          fls3 = (pnb(i-1,j)+dtb*pnbt(i-1,j)) - sps2%ps(ii-1,j)
          fls4 = (pnb(i+1,j)+dtb*pnbt(i+1,j)) - sps2%ps(ii+1,j)
          xpten(ii) = xpten(ii) + fcx*fls0 -                        &
                   & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
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
          fls0 = (pwb(i,jwb)+dtb*pwbt(i,jwb)) - sps2%ps(i,j)
          fls1 = (pwb(i-1,jwb)+dtb*pwbt(i-1,jwb)) - sps2%ps(i-1,j)
          fls2 = (pwb(i+1,jwb)+dtb*pwbt(i+1,jwb)) - sps2%ps(i+1,j)
          fls3 = (pwb(i,jwb-1)+dtb*pwbt(i,jwb-1)) - sps2%ps(i,j-1)
          fls4 = (pwb(i,jwb+1)+dtb*pwbt(i,jwb+1)) - sps2%ps(i,j+1)
          xpten(i) = xpten(i) + fcx*fls0 -                          &
                  & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      else if ( jj <= ip ) then
!-------east-boundary slice:
        fcx = fcoef*xfune(jsls,kz)
        gcx = gcoef*xfune(jsls,kz)
        do i = ibeg , iend
          fls0 = (peb(i,jeb)+dtb*pebt(i,jeb)) - sps2%ps(i,j)
          fls1 = (peb(i-1,jeb)+dtb*pebt(i-1,jeb)) - sps2%ps(i-1,j)
          fls2 = (peb(i+1,jeb)+dtb*pebt(i+1,jeb)) - sps2%ps(i+1,j)
          fls3 = (peb(i,jeb-1)+dtb*pebt(i,jeb-1)) - sps2%ps(i,j-1)
          fls4 = (peb(i,jeb+1)+dtb*pebt(i,jeb+1)) - sps2%ps(i,j+1)
          xpten(i) = xpten(i) + fcx*fls0 -                          &
                  & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
    end if
!
  end if
#endif
  call time_end(subroutine_name,idindx)
  end subroutine nudge_p
!
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
  use mod_bdycod
  use mod_atm_interface
  use mod_runparams
  use mod_service
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
#ifndef BAND
  integer :: jwb , jeb
#endif
!
  character (len=64) :: subroutine_name='sponge_p'
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
        if ( jsls <= ip ) then
          ften(i) = wg(jsls)*ften(i) + (d_one-wg(jsls))*pwbt(i,jwb)
        end if
      end do
    else if ( jj <= ip ) then
!------east-boundary slice:
      do i = ibeg , iend
        if ( jsls <= ip ) then
          ften(i) = wg(jsls)*ften(i) + (d_one-wg(jsls))*pebt(i,jeb)
        end if
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
  use mod_bdycod
  use mod_atm_interface
  use mod_runparams
  use mod_service
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
#ifndef BAND
  integer :: jwb , jeb
#endif
  character (len=64) :: subroutine_name='sponge_t'
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
        ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*xtb%sbt(i,k,j)
!.......north boundary:
        ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*xtb%nbt(i,k,j)
     end do
  end do

#else
!----------------------------------------------------------------------
!
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
!
  if ( jsls > ip ) then
!-----interior j slices:
    do i = 2 , ip
      ii = iy - i
      do k = 1 , kz
!.......south boundary:
        ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*xtb%sbt(i,k,j)
!.......north boundary:
        ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*xtb%nbt(i,k,j)
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
          ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*xtb%sbt(i,k,j)
!........north boundary:
          ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*xtb%nbt(i,k,j)
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
          if ( jsls <= ip ) then
            ften(i,k) = wg(jsls)*ften(i,k)+(d_one-wg(jsls))*xtb%wbt(i,k,jwb)
          end if
        end do
      end do
    else if ( jj <= ip ) then
!------east-boundary slice:
      do k = 1 , kz
        do i = ibeg , iend
          if ( jsls <= ip ) then
            ften(i,k) = wg(jsls)*ften(i,k)+(d_one-wg(jsls))*xtb%ebt(i,k,jeb)
          end if
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
  use mod_bdycod
  use mod_atm_interface
  use mod_runparams
  use mod_service
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
  integer :: jwb , jeb
#endif
  character (len=64) :: subroutine_name='spongeqv'
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
        ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*xqb%sbt(i,k,j)
!.......north boundary:
        ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*xqb%nbt(i,k,j)
     end do
  end do
#else
!----------------------------------------------------------------------
!
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
!
  if ( jsls > ip ) then
!-----interior j slices:
    do i = 2 , ip
      ii = iy - i
      do k = 1 , kz
!.......south boundary:
        ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*xqb%sbt(i,k,j)
!.......north boundary:
        ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*xqb%nbt(i,k,j)
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
          ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*xqb%sbt(i,k,j)
!........north boundary:
          ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*xqb%nbt(i,k,j)
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
          if ( jsls <= ip ) then
            ften(i,k) = wg(jsls)*ften(i,k)+(d_one-wg(jsls))*xqb%wbt(i,k,jwb)
          end if
        end do
      end do
    else if ( jj <= ip ) then
!------east-boundary slice:
      do k = 1 , kz
        do i = ibeg , iend
          if ( jsls <= ip ) then
            ften(i,k) = wg(jsls)*ften(i,k)+(d_one-wg(jsls))*xqb%ebt(i,k,jeb)
          end if
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
  use mod_bdycod
  use mod_atm_interface
  use mod_runparams
  use mod_service
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
  integer :: jew
#endif
  character (len=64) :: subroutine_name='sponge_u'
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
        ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*xub%sbt(i,k,j)
!.......north boundary:
        ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*xub%nbt(i,k,j)
     end do
  end do
#else
!----------------------------------------------------------------------
!
  jsls = j + myid*jxp
  jj = jxp1 - jsls
  if ( jj <= ip ) jsls = jj
  jew = jsls
  if ( jew > jxp ) jew = mod(jsls,jxp)
  if ( jew == 0 ) jew = jxp
!
  if ( jsls > ip ) then
!-----interior j slices:
    do i = 2 , ip
      ii = iy - i + 1
      do k = 1 , kz
!.......south boundary:
        ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*xub%sbt(i,k,j)
!.......north boundary:
        ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*xub%nbt(i,k,j)
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
          ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*xub%sbt(i,k,j)
!........north boundary:
          ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*xub%nbt(i,k,j)
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
          if ( jsls <= ip ) then
            ften(i,k) = wg(jsls)*ften(i,k)+(d_one-wg(jsls))*xub%wbt(i,k,jew)
          end if
        end do
      end do
    else if ( jj <= ip ) then
!------east-boundary slice:
      do k = 1 , kz
        do i = ibeg , iend
          if ( jsls <= ip ) then
            ften(i,k) = wg(jsls)*ften(i,k)+(d_one-wg(jsls))*xub%ebt(i,k,jew)
          end if
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
  use mod_bdycod
  use mod_atm_interface
  use mod_runparams
  use mod_service
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
  integer :: jew
#endif
  character (len=64) :: subroutine_name='sponge_v'
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
        ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*xvb%sbt(i,k,j)
!.......north boundary:
        ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*xvb%nbt(i,k,j)
     end do
  end do
#else
!----------------------------------------------------------------------
!
  jsls = j + myid*jxp
  jj = jxp1 - jsls
  if ( jj <= ip ) jsls = jj
  jew = jsls
  if ( jew > jxp ) jew = mod(jsls,jxp)
  if ( jew == 0 ) jew = jxp
!
  if ( jsls > ip ) then
!-----interior j slices:
    do i = 2 , ip
      ii = iy - i + 1
      do k = 1 , kz
!.......south boundary:
        ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*xvb%sbt(i,k,j)
!.......north boundary:
        ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*xvb%nbt(i,k,j)
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
          ften(i,k) = wg(i)*ften(i,k) + (d_one-wg(i))*xvb%sbt(i,k,j)
!........north boundary:
          ften(ii,k) = wg(i)*ften(ii,k) + (d_one-wg(i))*xvb%nbt(i,k,j)
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
          if ( jsls <= ip ) then
            ften(i,k) = wg(jsls)*ften(i,k)+(d_one-wg(jsls))*xvb%wbt(i,k,jew)
          end if
        end do
      end do
    else if ( jj <= ip ) then
!------east-boundary slice:
      do k = 1 , kz
        do i = ibeg , iend
          if ( jsls <= ip ) then
            ften(i,k) = wg(jsls)*ften(i,k)+(d_one-wg(jsls))*xvb%ebt(i,k,jew)
          end if
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
    use mod_runparams
    implicit none
    real(8) :: xfun
    integer , intent(in) :: mm
    xfun = dble(nspgd-mm)/dble(nspgd-2)
  end function xfun
!
  function xfune(mm,kk)
    use mod_runparams
    implicit none
    real(8) :: xfune
    integer , intent(in) :: mm , kk
    xfune = dexp(-dble(mm-2)/anudg(kk))
  end function xfune
!
  subroutine nudge(ldot,ip,fcoef,gcoef,xt,f,ften,j,ibdy,bnd)
!
  use mod_runparams
  use mod_service
  use mod_atm_interface , only : vbound
  implicit none
!
  logical , intent(in) :: ldot ! Dot flag
  integer , intent(in) :: ibdy , ip , j
  real(8) , intent(in) :: fcoef , gcoef , xt
  real(8) , intent(in) , dimension(iy,kz,-1:jxp+2) :: f
  type(vbound) , intent(in) :: bnd
  real(8) , intent(inout) , dimension(iy,kz,jxp) :: ften
!
  real(8) :: fcx , fls0 , fls1 , fls2 , fls3 , fls4 , gcx
  integer :: i , ido , ii , k
#ifndef BAND
  integer :: ibeg , iend , jj , jsls , jwb , jeb , jew
#endif
  character (len=64) :: subroutine_name='nudge'
  integer :: idindx=0
!
  call time_begin(subroutine_name,idindx)
!
!----------------------------------------------------------------------
!
  ido = 0
  if ( ldot ) then
    ido = 1
  end if

#ifdef BAND
!
!-----determine which relaxation method to use:linear/expon.
!
  if ( ibdy == 1 ) then
!
!---------use linear method
!
!------interior j slices:
     do i = 2 , ip
        ii = iym1 + ido - i + 1
        fcx = fcoef*xfun(i)
        gcx = gcoef*xfun(i)
        do k = 1 , kz
!.......south boundary:
          fls0 = (bnd%sb(i,k,j)+xt*bnd%sbt(i,k,j)) - f(i,k,j)
          fls1 = (bnd%sb(i,k,j-1)+xt*bnd%sbt(i,k,j-1)) - f(i,k,j-1)
          fls2 = (bnd%sb(i,k,j+1)+xt*bnd%sbt(i,k,j+1)) - f(i,k,j+1)
          fls3 = (bnd%sb(i-1,k,j)+xt*bnd%sbt(i-1,k,j)) - f(i-1,k,j)
          fls4 = (bnd%sb(i+1,k,j)+xt*bnd%sbt(i+1,k,j)) - f(i+1,k,j)
          ften(i,k,j) = ften(i,k,j) + fcx*fls0 - &
                    & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
          fls0 = (bnd%nb(i,k,j)+xt*bnd%nbt(i,k,j)) - f(ii,k,j)
          fls1 = (bnd%nb(i,k,j-1)+xt*bnd%nbt(i,k,j-1)) - f(ii,k,j-1)
          fls2 = (bnd%nb(i,k,j+1)+xt*bnd%nbt(i,k,j+1)) - f(ii,k,j+1)
          fls3 = (bnd%nb(i-1,k,j)+xt*bnd%nbt(i-1,k,j)) - f(ii-1,k,j)
          fls4 = (bnd%nb(i+1,k,j)+xt*bnd%nbt(i+1,k,j)) - f(ii+1,k,j)
          ften(ii,k,j) = ften(ii,k,j) + fcx*fls0 - &
                     & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
     end do

  else if ( ibdy == 5 ) then
 
!----------use exponential method
 
!------interior j slices:
     do i = 2 , ip
        ii = iym1 + ido - i + 1
        do k = 1 , kz
          fcx = fcoef*xfune(i,k)
          gcx = gcoef*xfune(i,k)
!........south boundary:
          fls0 = (bnd%sb(i,k,j)+xt*bnd%sbt(i,k,j)) - f(i,k,j)
          fls1 = (bnd%sb(i,k,j-1)+xt*bnd%sbt(i,k,j-1)) - f(i,k,j-1)
          fls2 = (bnd%sb(i,k,j+1)+xt*bnd%sbt(i,k,j+1)) - f(i,k,j+1)
          fls3 = (bnd%sb(i-1,k,j)+xt*bnd%sbt(i-1,k,j)) - f(i-1,k,j)
          fls4 = (bnd%sb(i+1,k,j)+xt*bnd%sbt(i+1,k,j)) - f(i+1,k,j)
          ften(i,k,j) = ften(i,k,j) + fcx*fls0 - &
                    & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
          fls0 = (bnd%nb(i,k,j)+xt*bnd%nbt(i,k,j)) - f(ii,k,j)
          fls1 = (bnd%nb(i,k,j-1)+xt*bnd%nbt(i,k,j-1)) - f(ii,k,j-1)
          fls2 = (bnd%nb(i,k,j+1)+xt*bnd%nbt(i,k,j+1)) - f(ii,k,j+1)
          fls3 = (bnd%nb(i-1,k,j)+xt*bnd%nbt(i-1,k,j)) - f(ii-1,k,j)
          fls4 = (bnd%nb(i+1,k,j)+xt*bnd%nbt(i+1,k,j)) - f(ii+1,k,j)
          ften(ii,k,j) = ften(ii,k,j) + fcx*fls0 - &
                     & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
     end do
  end if
#else
!----------------------------------------------------------------------
!
  jsls = j + myid*jxp
  if ( ldot ) then
    jj = jxp1 - jsls
    if ( jj <= ip ) jsls = jj
    jew = jsls
    if ( jew > jxp ) jew = mod(jsls,jxp)
    if ( jew == 0 ) jew = jxp
    jwb = jew
    jeb = jew
  else
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
  end if
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
        ii = iym1 + ido - i + 1
        fcx = fcoef*xfun(i)
        gcx = gcoef*xfun(i)
        do k = 1 , kz
!.......south boundary:
          fls0 = (bnd%sb(i,k,j)+xt*bnd%sbt(i,k,j)) - f(i,k,j)
          fls1 = (bnd%sb(i,k,j-1)+xt*bnd%sbt(i,k,j-1)) - f(i,k,j-1)
          fls2 = (bnd%sb(i,k,j+1)+xt*bnd%sbt(i,k,j+1)) - f(i,k,j+1)
          fls3 = (bnd%sb(i-1,k,j)+xt*bnd%sbt(i-1,k,j)) - f(i-1,k,j)
          fls4 = (bnd%sb(i+1,k,j)+xt*bnd%sbt(i+1,k,j)) - f(i+1,k,j)
          ften(i,k,j) = ften(i,k,j) + fcx*fls0 - &
                    & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
          fls0 = (bnd%nb(i,k,j)+xt*bnd%nbt(i,k,j)) - f(ii,k,j)
          fls1 = (bnd%nb(i,k,j-1)+xt*bnd%nbt(i,k,j-1)) - f(ii,k,j-1)
          fls2 = (bnd%nb(i,k,j+1)+xt*bnd%nbt(i,k,j+1)) - f(ii,k,j+1)
          fls3 = (bnd%nb(i-1,k,j)+xt*bnd%nbt(i-1,k,j)) - f(ii-1,k,j)
          fls4 = (bnd%nb(i+1,k,j)+xt*bnd%nbt(i+1,k,j)) - f(ii+1,k,j)
          ften(ii,k,j) = ften(ii,k,j) + fcx*fls0 - &
                     & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end do
!
    else if ( jsls <= ip ) then
!------east or west boundary slices:
      ibeg = 2
      iend = iym1 + ido - 1
      if ( jsls > 2 ) then
        do i = 2 , jsls - 1
          ii = iym1 + ido - i + 1
          fcx = fcoef*xfun(i)
          gcx = gcoef*xfun(i)
          do k = 1 , kz
!........south  boundary:
            fls0 = (bnd%sb(i,k,j)+xt*bnd%sbt(i,k,j)) - f(i,k,j)
            fls1 = (bnd%sb(i,k,j-1)+xt*bnd%sbt(i,k,j-1)) - f(i,k,j-1)
            fls2 = (bnd%sb(i,k,j+1)+xt*bnd%sbt(i,k,j+1)) - f(i,k,j+1)
            fls3 = (bnd%sb(i-1,k,j)+xt*bnd%sbt(i-1,k,j)) - f(i-1,k,j)
            fls4 = (bnd%sb(i+1,k,j)+xt*bnd%sbt(i+1,k,j)) - f(i+1,k,j)
            ften(i,k,j) = ften(i,k,j) + fcx*fls0 - &
                      & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
!.........north boundary:
            fls0 = (bnd%nb(i,k,j)+xt*bnd%nbt(i,k,j)) - f(ii,k,j)
            fls1 = (bnd%nb(i,k,j-1)+xt*bnd%nbt(i,k,j-1)) - f(ii,k,j-1)
            fls2 = (bnd%nb(i,k,j+1)+xt*bnd%nbt(i,k,j+1)) - f(ii,k,j+1)
            fls3 = (bnd%nb(i-1,k,j)+xt*bnd%nbt(i-1,k,j)) - f(ii-1,k,j)
            fls4 = (bnd%nb(i+1,k,j)+xt*bnd%nbt(i+1,k,j)) - f(ii+1,k,j)
            ften(ii,k,j) = ften(ii,k,j) + fcx*fls0 - &
                       & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
        end do
        ibeg = jsls
        iend = iym1 + ido - jsls + 1
      end if
!
      if ( jj > ip ) then
!-------west-boundary slice:
        fcx = fcoef*xfun(jsls)
        gcx = gcoef*xfun(jsls)
        do k = 1 , kz
          do i = ibeg , iend
            fls0 = (bnd%wb(i,k,jwb)+xt*bnd%wbt(i,k,jwb)) - f(i,k,j)
            fls1 = (bnd%wb(i-1,k,jwb)+xt*bnd%wbt(i-1,k,jwb)) - f(i-1,k,j)
            fls2 = (bnd%wb(i+1,k,jwb)+xt*bnd%wbt(i+1,k,jwb)) - f(i+1,k,j)
            fls3 = (bnd%wb(i,k,jwb-1)+xt*bnd%wbt(i,k,jwb-1)) - f(i,k,j-1)
            fls4 = (bnd%wb(i,k,jwb+1)+xt*bnd%wbt(i,k,jwb+1)) - f(i,k,j+1)
            ften(i,k,j) = ften(i,k,j) + fcx*fls0 - &
                      & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
        end do
      else if ( jj <= ip ) then
!-------east-boundary slice:
        fcx = fcoef*xfun(jsls)
        gcx = gcoef*xfun(jsls)
        do k = 1 , kz
          do i = ibeg , iend
            fls0 = (bnd%eb(i,k,jeb)+xt*bnd%ebt(i,k,jeb)) - f(i,k,j)
            fls1 = (bnd%eb(i-1,k,jeb)+xt*bnd%ebt(i-1,k,jeb)) - f(i-1,k,j)
            fls2 = (bnd%eb(i+1,k,jeb)+xt*bnd%ebt(i+1,k,jeb)) - f(i+1,k,j)
            fls3 = (bnd%eb(i,k,jeb-1)+xt*bnd%ebt(i,k,jeb-1)) - f(i,k,j-1)
            fls4 = (bnd%eb(i,k,jeb+1)+xt*bnd%ebt(i,k,jeb+1)) - f(i,k,j+1)
            ften(i,k,j) = ften(i,k,j) + fcx*fls0 -  &
                      & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
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
        ii = iym1 + ido - i + 1
        do k = 1 , kz
          fcx = fcoef*xfune(i,k)
          gcx = gcoef*xfune(i,k)
!........south boundary:
          fls0 = (bnd%sb(i,k,j)+xt*bnd%sbt(i,k,j)) - f(i,k,j)
          fls1 = (bnd%sb(i,k,j-1)+xt*bnd%sbt(i,k,j-1)) - f(i,k,j-1)
          fls2 = (bnd%sb(i,k,j+1)+xt*bnd%sbt(i,k,j+1)) - f(i,k,j+1)
          fls3 = (bnd%sb(i-1,k,j)+xt*bnd%sbt(i-1,k,j)) - f(i-1,k,j)
          fls4 = (bnd%sb(i+1,k,j)+xt*bnd%sbt(i+1,k,j)) - f(i+1,k,j)
          ften(i,k,j) = ften(i,k,j) + fcx*fls0 - &
                    & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
          fls0 = (bnd%nb(i,k,j)+xt*bnd%nbt(i,k,j)) - f(ii,k,j)
          fls1 = (bnd%nb(i,k,j-1)+xt*bnd%nbt(i,k,j-1)) - f(ii,k,j-1)
          fls2 = (bnd%nb(i,k,j+1)+xt*bnd%nbt(i,k,j+1)) - f(ii,k,j+1)
          fls3 = (bnd%nb(i-1,k,j)+xt*bnd%nbt(i-1,k,j)) - f(ii-1,k,j)
          fls4 = (bnd%nb(i+1,k,j)+xt*bnd%nbt(i+1,k,j)) - f(ii+1,k,j)
          ften(ii,k,j) = ften(ii,k,j) + fcx*fls0 - &
                     & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end do
!
    else if ( jsls <= ip ) then
!------east or west boundary slices:
      ibeg = 2
      iend = iym1 + ido - 1
      if ( jsls > 2 ) then
        do i = 2 , jsls - 1
          ii = iym1 + ido - i + 1
          do k = 1 , kz
            fcx = fcoef*xfune(i,k)
            gcx = gcoef*xfune(i,k)
!.........south boundary:
            fls0 = (bnd%sb(i,k,j)+xt*bnd%sbt(i,k,j)) - f(i,k,j)
            fls1 = (bnd%sb(i,k,j-1)+xt*bnd%sbt(i,k,j-1)) - f(i,k,j-1)
            fls2 = (bnd%sb(i,k,j+1)+xt*bnd%sbt(i,k,j+1)) - f(i,k,j+1)
            fls3 = (bnd%sb(i-1,k,j)+xt*bnd%sbt(i-1,k,j)) - f(i-1,k,j)
            fls4 = (bnd%sb(i+1,k,j)+xt*bnd%sbt(i+1,k,j)) - f(i+1,k,j)
            ften(i,k,j) = ften(i,k,j) + fcx*fls0 - &
                      & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
!.........north boundary:
            fls0 = (bnd%nb(i,k,j)+xt*bnd%nbt(i,k,j)) - f(ii,k,j)
            fls1 = (bnd%nb(i,k,j-1)+xt*bnd%nbt(i,k,j-1)) - f(ii,k,j-1)
            fls2 = (bnd%nb(i,k,j+1)+xt*bnd%nbt(i,k,j+1)) - f(ii,k,j+1)
            fls3 = (bnd%nb(i-1,k,j)+xt*bnd%nbt(i-1,k,j)) - f(ii-1,k,j)
            fls4 = (bnd%nb(i+1,k,j)+xt*bnd%nbt(i+1,k,j)) - f(ii+1,k,j)
            ften(ii,k,j) = ften(ii,k,j) + fcx*fls0 - &
                       & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
        end do
        ibeg = jsls
        iend = iym1 + ido - jsls + 1
      end if
!
      if ( jj > ip ) then
!-------west-boundary slice:
        do k = 1 , kz
          fcx = fcoef*xfune(jsls,k)
          gcx = gcoef*xfune(jsls,k)
          do i = ibeg , iend
            fls0 = (bnd%wb(i,k,jwb)+xt*bnd%wbt(i,k,jwb)) - f(i,k,j)
            fls1 = (bnd%wb(i-1,k,jwb)+xt*bnd%wbt(i-1,k,jwb)) - f(i-1,k,j)
            fls2 = (bnd%wb(i+1,k,jwb)+xt*bnd%wbt(i+1,k,jwb)) - f(i+1,k,j)
            fls3 = (bnd%wb(i,k,jwb-1)+xt*bnd%wbt(i,k,jwb-1)) - f(i,k,j-1)
            fls4 = (bnd%wb(i,k,jwb+1)+xt*bnd%wbt(i,k,jwb+1)) - f(i,k,j+1)
            ften(i,k,j) = ften(i,k,j) + fcx*fls0 - &
                      & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
        end do
      else if ( jj <= ip ) then
!-------east-boundary slice:
        do k = 1 , kz
          fcx = fcoef*xfune(jsls,k)
          gcx = gcoef*xfune(jsls,k)
          do i = ibeg , iend
            fls0 = (bnd%eb(i,k,jeb)+xt*bnd%ebt(i,k,jeb)) - f(i,k,j)
            fls1 = (bnd%eb(i-1,k,jeb)+xt*bnd%ebt(i-1,k,jeb)) - f(i-1,k,j)
            fls2 = (bnd%eb(i+1,k,jeb)+xt*bnd%ebt(i+1,k,jeb)) - f(i+1,k,j)
            fls3 = (bnd%eb(i,k,jeb-1)+xt*bnd%ebt(i,k,jeb-1)) - f(i,k,j-1)
            fls4 = (bnd%eb(i,k,jeb+1)+xt*bnd%ebt(i,k,jeb+1)) - f(i,k,j+1)
            ften(i,k,j) = ften(i,k,j) + fcx*fls0 - &
                      & gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
        end do
      end if
    end if
!
  end if
#endif
  call time_end(subroutine_name,idindx)

  end subroutine nudge
!
end module mod_nudge
