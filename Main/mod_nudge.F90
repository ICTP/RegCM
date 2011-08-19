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

  interface nudge
    module procedure nudge3d , nudge2d
  end interface nudge
!
  public :: nudge
!
  contains
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
     ften(i) = wg(i)*ften(i) + (d_one-wg(i))*xpsb%sbt(i,j)
!.......north boundary:
     ften(ii) = wg(i)*ften(ii) + (d_one-wg(i))*xpsb%nbt(i,j)
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
      ften(i) = wg(i)*ften(i) + (d_one-wg(i))*xpsb%sbt(i,j)
!.......north boundary:
      ften(ii) = wg(i)*ften(ii) + (d_one-wg(i))*xpsb%nbt(i,j)
    end do
!
  else if ( jsls <= ip ) then
    ibeg = 2
    iend = iym1 - 1
    if ( jsls > 2 ) then
      do i = 2 , jsls - 1
        ii = iy - i
!........south boundary:
        ften(i) = wg(i)*ften(i) + (d_one-wg(i))*xpsb%sbt(i,j)
!........north boundary:
        ften(ii) = wg(i)*ften(ii) + (d_one-wg(i))*xpsb%nbt(i,j)
      end do
      ibeg = jsls
      iend = iy - jsls
    end if
!
    if ( jj > ip ) then
!------west-boundary slice:
      do i = ibeg , iend
        if ( jsls <= ip ) then
          ften(i) = wg(jsls)*ften(i) + (d_one-wg(jsls))*xpsb%wbt(i,jwb)
        end if
      end do
    else if ( jj <= ip ) then
!------east-boundary slice:
      do i = ibeg , iend
        if ( jsls <= ip ) then
          ften(i) = wg(jsls)*ften(i) + (d_one-wg(jsls))*xpsb%ebt(i,jeb)
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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     these subroutines apply relaxation boundary conditions to the   c
!     tendency term - ften - of variable f                            c
!                                                                     c
!     ldot  : logical dot (u,v) / cross (t,q,p) flag                  c
!                                                                     c
!     ip    : is the number of slices affected by nudging.            c
!                                                                     c
!     xt    : is the time in seconds for variable f                   c
!                                                                     c
!     fcoef : are the coefficients for the newtonian term.            c
!                                                                     c
!     gcoef : are the coefficients for the diffusion term.            c
!                                                                     c
!     ften  : is the tendency calculated from the model.              c
!                                                                     c
!     j     : is the j'th slice of the tendency to be adjusted.       c
!                                                                     c
!     nk    : is the number of vertical level to be adjusted.         c
!                                                                     c
!     ibdy  : type of boundary condition relaxation, 1=linear         c
!              5 = exponential                                        c
!                                                                     c
!     bnd   : Boundary condition data structure                       c
!             2D or 3D (managed by interface declaration)             c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine nudge3d(ldot,ip,fcoef,gcoef,xt,f,ften,j,nk,ibdy,bnd)
!
  use mod_runparams
  use mod_service
  use mod_atm_interface , only : v3dbound
  implicit none
!
  logical , intent(in) :: ldot ! Dot flag
  integer , intent(in) :: ibdy , nk , ip , j
  real(8) , intent(in) :: fcoef , gcoef , xt
  real(8) , intent(in) , dimension(iy,nk,-1:jxp+2) :: f
  type(v3dbound) , intent(in) :: bnd
  real(8) , intent(inout) , dimension(iy,nk,jxp) :: ften
!
  real(8) :: fcx , fls0 , fls1 , fls2 , fls3 , fls4 , gcx
  integer :: i , ido , ii , k
#ifndef BAND
  integer :: ibeg , iend , jj , jsls , jwb , jeb , jew
#endif
  character (len=64) :: subroutine_name='nudge3d'
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
        do k = 1 , nk
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
        do k = 1 , nk
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
        do k = 1 , nk
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
          do k = 1 , nk
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
        do k = 1 , nk
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
        do k = 1 , nk
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
        do k = 1 , nk
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
          do k = 1 , nk
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
        do k = 1 , nk
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
        do k = 1 , nk
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

  end subroutine nudge3d
!
! ###################################################################
!
  subroutine nudge2d(ldot,ip,fcoef,gcoef,xt,f,ften,j,nk,ibdy,bnd)
!
  use mod_runparams
  use mod_service
  use mod_atm_interface , only : v2dbound
  implicit none
!
  logical , intent(in) :: ldot ! Dot flag
  integer , intent(in) :: ibdy , nk , ip , j
  real(8) , intent(in) :: fcoef , gcoef , xt
  real(8) , intent(in) , dimension(iy,-1:jxp+2) :: f
  type(v2dbound) , intent(in) :: bnd
  real(8) , intent(inout) , dimension(iy,jxp) :: ften
!
  real(8) :: fcx , fls0 , fls1 , fls2 , fls3 , fls4 , gcx
  integer :: i , ido , ii , k
#ifndef BAND
  integer :: ibeg , iend , jj , jsls , jwb , jeb , jew
#endif
  character (len=64) :: subroutine_name='nudge2d'
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
!.......south boundary:
        fls0 = (bnd%sb(i,j)+xt*bnd%sbt(i,j)) - f(i,j)
        fls1 = (bnd%sb(i,j-1)+xt*bnd%sbt(i,j-1)) - f(i,j-1)
        fls2 = (bnd%sb(i,j+1)+xt*bnd%sbt(i,j+1)) - f(i,j+1)
        fls3 = (bnd%sb(i-1,j)+xt*bnd%sbt(i-1,j)) - f(i-1,j)
        fls4 = (bnd%sb(i+1,j)+xt*bnd%sbt(i+1,j)) - f(i+1,j)
        ften(i,j) = ften(i,j) + fcx*fls0 - &
                      gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
        fls0 = (bnd%nb(i,j)+xt*bnd%nbt(i,j)) - f(ii,j)
        fls1 = (bnd%nb(i,j-1)+xt*bnd%nbt(i,j-1)) - f(ii,j-1)
        fls2 = (bnd%nb(i,j+1)+xt*bnd%nbt(i,j+1)) - f(ii,j+1)
        fls3 = (bnd%nb(i-1,j)+xt*bnd%nbt(i-1,j)) - f(ii-1,j)
        fls4 = (bnd%nb(i+1,j)+xt*bnd%nbt(i+1,j)) - f(ii+1,j)
        ften(ii,j) = ften(ii,j) + fcx*fls0 - &
                       gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
     end do

  else if ( ibdy == 5 ) then
 
!----------use exponential method
 
!------interior j slices:
     do i = 2 , ip
        ii = iym1 + ido - i + 1
        fcx = fcoef*xfune(i,kz)
        gcx = gcoef*xfune(i,kz)
!........south boundary:
        fls0 = (bnd%sb(i,j)+xt*bnd%sbt(i,j)) - f(i,j)
        fls1 = (bnd%sb(i,j-1)+xt*bnd%sbt(i,j-1)) - f(i,j-1)
        fls2 = (bnd%sb(i,j+1)+xt*bnd%sbt(i,j+1)) - f(i,j+1)
        fls3 = (bnd%sb(i-1,j)+xt*bnd%sbt(i-1,j)) - f(i-1,j)
        fls4 = (bnd%sb(i+1,j)+xt*bnd%sbt(i+1,j)) - f(i+1,j)
        ften(i,j) = ften(i,j) + fcx*fls0 - &
                      gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
        fls0 = (bnd%nb(i,j)+xt*bnd%nbt(i,j)) - f(ii,j)
        fls1 = (bnd%nb(i,j-1)+xt*bnd%nbt(i,j-1)) - f(ii,j-1)
        fls2 = (bnd%nb(i,j+1)+xt*bnd%nbt(i,j+1)) - f(ii,j+1)
        fls3 = (bnd%nb(i-1,j)+xt*bnd%nbt(i-1,j)) - f(ii-1,j)
        fls4 = (bnd%nb(i+1,j)+xt*bnd%nbt(i+1,j)) - f(ii+1,j)
        ften(ii,j) = ften(ii,j) + fcx*fls0 - &
                       gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
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
!.......south boundary:
        fls0 = (bnd%sb(i,j)+xt*bnd%sbt(i,j)) - f(i,j)
        fls1 = (bnd%sb(i,j-1)+xt*bnd%sbt(i,j-1)) - f(i,j-1)
        fls2 = (bnd%sb(i,j+1)+xt*bnd%sbt(i,j+1)) - f(i,j+1)
        fls3 = (bnd%sb(i-1,j)+xt*bnd%sbt(i-1,j)) - f(i-1,j)
        fls4 = (bnd%sb(i+1,j)+xt*bnd%sbt(i+1,j)) - f(i+1,j)
        ften(i,j) = ften(i,j) + fcx*fls0 - &
                      gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
        fls0 = (bnd%nb(i,j)+xt*bnd%nbt(i,j)) - f(ii,j)
        fls1 = (bnd%nb(i,j-1)+xt*bnd%nbt(i,j-1)) - f(ii,j-1)
        fls2 = (bnd%nb(i,j+1)+xt*bnd%nbt(i,j+1)) - f(ii,j+1)
        fls3 = (bnd%nb(i-1,j)+xt*bnd%nbt(i-1,j)) - f(ii-1,j)
        fls4 = (bnd%nb(i+1,j)+xt*bnd%nbt(i+1,j)) - f(ii+1,j)
        ften(ii,j) = ften(ii,j) + fcx*fls0 - &
                       gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
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
!........south  boundary:
          fls0 = (bnd%sb(i,j)+xt*bnd%sbt(i,j)) - f(i,j)
          fls1 = (bnd%sb(i,j-1)+xt*bnd%sbt(i,j-1)) - f(i,j-1)
          fls2 = (bnd%sb(i,j+1)+xt*bnd%sbt(i,j+1)) - f(i,j+1)
          fls3 = (bnd%sb(i-1,j)+xt*bnd%sbt(i-1,j)) - f(i-1,j)
          fls4 = (bnd%sb(i+1,j)+xt*bnd%sbt(i+1,j)) - f(i+1,j)
          ften(i,j) = ften(i,j) + fcx*fls0 - &
                        gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
!.........north boundary:
          fls0 = (bnd%nb(i,j)+xt*bnd%nbt(i,j)) - f(ii,j)
          fls1 = (bnd%nb(i,j-1)+xt*bnd%nbt(i,j-1)) - f(ii,j-1)
          fls2 = (bnd%nb(i,j+1)+xt*bnd%nbt(i,j+1)) - f(ii,j+1)
          fls3 = (bnd%nb(i-1,j)+xt*bnd%nbt(i-1,j)) - f(ii-1,j)
          fls4 = (bnd%nb(i+1,j)+xt*bnd%nbt(i+1,j)) - f(ii+1,j)
          ften(ii,j) = ften(ii,j) + fcx*fls0 - &
                         gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
        ibeg = jsls
        iend = iym1 + ido - jsls + 1
      end if
!
      if ( jj > ip ) then
!-------west-boundary slice:
        fcx = fcoef*xfun(jsls)
        gcx = gcoef*xfun(jsls)
        do i = ibeg , iend
          fls0 = (bnd%wb(i,jwb)+xt*bnd%wbt(i,jwb)) - f(i,j)
          fls1 = (bnd%wb(i-1,jwb)+xt*bnd%wbt(i-1,jwb)) - f(i-1,j)
          fls2 = (bnd%wb(i+1,jwb)+xt*bnd%wbt(i+1,jwb)) - f(i+1,j)
          fls3 = (bnd%wb(i,jwb-1)+xt*bnd%wbt(i,jwb-1)) - f(i,j-1)
          fls4 = (bnd%wb(i,jwb+1)+xt*bnd%wbt(i,jwb+1)) - f(i,j+1)
          ften(i,j) = ften(i,j) + fcx*fls0 - &
                        gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      else if ( jj <= ip ) then
!-------east-boundary slice:
        fcx = fcoef*xfun(jsls)
        gcx = gcoef*xfun(jsls)
        do i = ibeg , iend
          fls0 = (bnd%eb(i,jeb)+xt*bnd%ebt(i,jeb)) - f(i,j)
          fls1 = (bnd%eb(i-1,jeb)+xt*bnd%ebt(i-1,jeb)) - f(i-1,j)
          fls2 = (bnd%eb(i+1,jeb)+xt*bnd%ebt(i+1,jeb)) - f(i+1,j)
          fls3 = (bnd%eb(i,jeb-1)+xt*bnd%ebt(i,jeb-1)) - f(i,j-1)
          fls4 = (bnd%eb(i,jeb+1)+xt*bnd%ebt(i,jeb+1)) - f(i,j+1)
          ften(i,j) = ften(i,j) + fcx*fls0 -  &
                        gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
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
        fcx = fcoef*xfune(i,kz)
        gcx = gcoef*xfune(i,kz)
!........south boundary:
        fls0 = (bnd%sb(i,j)+xt*bnd%sbt(i,j)) - f(i,j)
        fls1 = (bnd%sb(i,j-1)+xt*bnd%sbt(i,j-1)) - f(i,j-1)
        fls2 = (bnd%sb(i,j+1)+xt*bnd%sbt(i,j+1)) - f(i,j+1)
        fls3 = (bnd%sb(i-1,j)+xt*bnd%sbt(i-1,j)) - f(i-1,j)
        fls4 = (bnd%sb(i+1,j)+xt*bnd%sbt(i+1,j)) - f(i+1,j)
        ften(i,j) = ften(i,j) + fcx*fls0 - &
                      gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
!........north boundary:
        fls0 = (bnd%nb(i,j)+xt*bnd%nbt(i,j)) - f(ii,j)
        fls1 = (bnd%nb(i,j-1)+xt*bnd%nbt(i,j-1)) - f(ii,j-1)
        fls2 = (bnd%nb(i,j+1)+xt*bnd%nbt(i,j+1)) - f(ii,j+1)
        fls3 = (bnd%nb(i-1,j)+xt*bnd%nbt(i-1,j)) - f(ii-1,j)
        fls4 = (bnd%nb(i+1,j)+xt*bnd%nbt(i+1,j)) - f(ii+1,j)
        ften(ii,j) = ften(ii,j) + fcx*fls0 - &
                       gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
      end do
!
    else if ( jsls <= ip ) then
!------east or west boundary slices:
      ibeg = 2
      iend = iym1 + ido - 1
      if ( jsls > 2 ) then
        do i = 2 , jsls - 1
          ii = iym1 + ido - i + 1
          fcx = fcoef*xfune(i,kz)
          gcx = gcoef*xfune(i,kz)
!.........south boundary:
          fls0 = (bnd%sb(i,j)+xt*bnd%sbt(i,j)) - f(i,j)
          fls1 = (bnd%sb(i,j-1)+xt*bnd%sbt(i,j-1)) - f(i,j-1)
          fls2 = (bnd%sb(i,j+1)+xt*bnd%sbt(i,j+1)) - f(i,j+1)
          fls3 = (bnd%sb(i-1,j)+xt*bnd%sbt(i-1,j)) - f(i-1,j)
          fls4 = (bnd%sb(i+1,j)+xt*bnd%sbt(i+1,j)) - f(i+1,j)
          ften(i,j) = ften(i,j) + fcx*fls0 - &
                        gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
!.........north boundary:
          fls0 = (bnd%nb(i,j)+xt*bnd%nbt(i,j)) - f(ii,j)
          fls1 = (bnd%nb(i,j-1)+xt*bnd%nbt(i,j-1)) - f(ii,j-1)
          fls2 = (bnd%nb(i,j+1)+xt*bnd%nbt(i,j+1)) - f(ii,j+1)
          fls3 = (bnd%nb(i-1,j)+xt*bnd%nbt(i-1,j)) - f(ii-1,j)
          fls4 = (bnd%nb(i+1,j)+xt*bnd%nbt(i+1,j)) - f(ii+1,j)
          ften(ii,j) = ften(ii,j) + fcx*fls0 - &
                         gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
        ibeg = jsls
        iend = iym1 + ido - jsls + 1
      end if
!
      if ( jj > ip ) then
!-------west-boundary slice:
        fcx = fcoef*xfune(jsls,kz)
        gcx = gcoef*xfune(jsls,kz)
        do i = ibeg , iend
          fls0 = (bnd%wb(i,jwb)+xt*bnd%wbt(i,jwb)) - f(i,j)
          fls1 = (bnd%wb(i-1,jwb)+xt*bnd%wbt(i-1,jwb)) - f(i-1,j)
          fls2 = (bnd%wb(i+1,jwb)+xt*bnd%wbt(i+1,jwb)) - f(i+1,j)
          fls3 = (bnd%wb(i,jwb-1)+xt*bnd%wbt(i,jwb-1)) - f(i,j-1)
          fls4 = (bnd%wb(i,jwb+1)+xt*bnd%wbt(i,jwb+1)) - f(i,j+1)
          ften(i,j) = ften(i,j) + fcx*fls0 - &
                        gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      else if ( jj <= ip ) then
!-------east-boundary slice:
        fcx = fcoef*xfune(jsls,kz)
        gcx = gcoef*xfune(jsls,kz)
        do i = ibeg , iend
          fls0 = (bnd%eb(i,jeb)+xt*bnd%ebt(i,jeb)) - f(i,j)
          fls1 = (bnd%eb(i-1,jeb)+xt*bnd%ebt(i-1,jeb)) - f(i-1,j)
          fls2 = (bnd%eb(i+1,jeb)+xt*bnd%ebt(i+1,jeb)) - f(i+1,j)
          fls3 = (bnd%eb(i,jeb-1)+xt*bnd%ebt(i,jeb-1)) - f(i,j-1)
          fls4 = (bnd%eb(i,jeb+1)+xt*bnd%ebt(i,jeb+1)) - f(i,j+1)
          ften(i,j) = ften(i,j) + fcx*fls0 - &
                        gcx*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
    end if
!
  end if
#endif
  call time_end(subroutine_name,idindx)

  end subroutine nudge2d
!
end module mod_nudge
