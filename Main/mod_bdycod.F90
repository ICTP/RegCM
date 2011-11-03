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

module mod_bdycod
!
! Storage and subroutines for input of boundary values and tendencies
! of p*u, p*v, p*t, p*qv and p*, and outmost 2 slices of u and v for
! large domain.
! Relaxation and Sponge Boundary Conditions routines
!
  use mod_dynparam
  use mod_mppparam
  use mod_memutil
  use mod_atm_interface
  use mod_che_interface
  use mod_pbl_interface , only : set_tke_bc
  use mod_lm_interface
  use mod_mpmessage 
  use mod_ncio
  use mod_mppio
  use mod_service
!
  private
!
  public :: allocate_mod_bdycon , bdyin , bdyval
!
#ifndef BAND
  public :: uj1 , uj2 , ujlx , ujl
  public :: vj1 , vj2 , vjlx , vjl
#endif
  public :: ui1 , ui2 , uilx , uil
  public :: vi1 , vi2 , vilx , vil
  public :: ts0 , so0
  public :: ts1 ! FOR DCSST
!
  real(8) , pointer , dimension(:,:) :: ui1 , ui2 , uil , uilx , &
                                        vi1 , vi2 , vil , vilx
#ifndef BAND
  real(8) , pointer , dimension(:,:) :: uj1 , uj2 , ujl , ujlx , &
                                        vj1 , vj2 , vjl , vjlx
#endif
  real(8) , pointer , dimension(:,:,:) :: so0 , so1
  real(8) , pointer , dimension(:,:) :: ts0 , ts1
!
  interface nudge
    module procedure nudge3d , nudge2d
  end interface nudge
!
  interface sponge
    module procedure sponge3d , sponge2d
  end interface sponge
!
  public :: sponge , nudge
!
  contains
!
  subroutine allocate_mod_bdycon
    implicit none
    character (len=64) :: subroutine_name='allocate_mod_bdycon'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)

    call getmem2d(ts0,1,iy,1,jxp,'bdycon:ts0')
    call getmem2d(ts1,1,iy,1,jxp,'bdycon:ts1')
!
    call getmem3d(so0,1,iy,1,kz,1,jxp,'bdycon:so0')
    call getmem3d(so1,1,iy,1,kz,1,jxp,'bdycon:so1')
    call getmem2d(ui1,1,kz,0,jxp+1,'bdycon:ui1')
    call getmem2d(ui2,1,kz,0,jxp+1,'bdycon:ui2')
    call getmem2d(uil,1,kz,0,jxp+1,'bdycon:uil')
    call getmem2d(uilx,1,kz,0,jxp+1,'bdycon:uilx')
    call getmem2d(vi1,1,kz,0,jxp+1,'bdycon:vi1')
    call getmem2d(vi2,1,kz,0,jxp+1,'bdycon:vi2')
    call getmem2d(vil,1,kz,0,jxp+1,'bdycon:vil')
    call getmem2d(vilx,1,kz,0,jxp+1,'bdycon:vilx')
!
#ifndef BAND
    call getmem2d(uj1,1,iy,1,kz,'bdycon:uj1')
    call getmem2d(uj2,1,iy,1,kz,'bdycon:uj2')
    call getmem2d(ujl,1,iy,1,kz,'bdycon:ujl')
    call getmem2d(ujlx,1,iy,1,kz,'bdycon:ujlx')
    call getmem2d(vj1,1,iy,1,kz,'bdycon:vj1')
    call getmem2d(vj2,1,iy,1,kz,'bdycon:vj2')
    call getmem2d(vjl,1,iy,1,kz,'bdycon:vjl')
    call getmem2d(vjlx,1,iy,1,kz,'bdycon:vjlx')
#endif
    call time_end(subroutine_name,idindx)
  end subroutine allocate_mod_bdycon
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! this subroutine reads in the boundary conditions.
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine bdyin
!
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
!
    integer :: i , j , k , nn , nnb , mmrec
    integer :: ierr , ndeb , ndwb , nxeb , nxwb
#ifndef BAND
    integer :: nkk
#endif
    character(len=32) :: appdat
    real(8) , dimension(iy,jxp) :: psdot , tdum
    integer :: n
    character (len=64) :: subroutine_name='bdyin'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
    if ( myid == 0 ) then
      if ( ehso4 ) then
        do k = 1 , kz
          do j = 1 , jendl
            do i = 1 , iy
              sulfate(i,k,j) = so0(i,k,j)
            end do
          end do
        end do
      end if
      bdydate2 = bdydate2 + intbdy
      write (6,'(a,i10)') 'SEARCH BC data for ', toint10(bdydate2)
      mmrec = icbc_search(bdydate2)
      if (mmrec < 0) then
        call open_icbc(monfirst(bdydate2))
        mmrec = icbc_search(bdydate2)
        if (mmrec < 0) then
          appdat = tochar(bdydate2)
          call fatal(__FILE__,__LINE__,'ICBC for '//appdat//' not found')
        end if
      end if
      call read_icbc(ps1_io,ts1_io,ub1_io,vb1_io,tb1_io,qb1_io,so1_io)
      ps1_io = ps1_io*d_r10
      do j = 1 , jx
        do k = 1 , kz
          do i = 1 , iy
            sav_0(i,k,j)      = ub1_io(i,k,j)
            sav_0(i,kz+k,j)   = vb1_io(i,k,j)
            sav_0(i,kz*2+k,j) = qb1_io(i,k,j)
            sav_0(i,kz*3+k,j) = tb1_io(i,k,j)
          end do
        end do
        do i = 1 , iy
          sav_0(i,kz*4+1,j) = ps1_io(i,j)
          sav_0(i,kz*4+2,j) = ts1_io(i,j)
        end do
      end do
      if ( ehso4 ) then
        do j = 1 , jx
          do k = 1 , kz
            do i = 1 , iy
              sav_0s(i,k,j) = so1_io(i,k,j)
            end do
          end do
        end do
      end if
    end if
!
    call date_bcast(bdydate2,0,mycomm,ierr)
!
    call mpi_scatter(sav_0,iy*(kz*4+2)*jxp,mpi_real8,        &
                     sav0, iy*(kz*4+2)*jxp,mpi_real8,        &
                     0,mycomm,ierr)
    do j = 1 , jendl
      do k = 1 , kz
        do i = 1 , iy
          xub%b1(i,k,j) = sav0(i,k,j)
          xvb%b1(i,k,j) = sav0(i,kz+k,j)
          xqb%b1(i,k,j) = sav0(i,kz*2+k,j)
          xtb%b1(i,k,j) = sav0(i,kz*3+k,j)
        end do
      end do
      do i = 1 , iy
        xpsb%b1(i,j) = sav0(i,kz*4+1,j)
        ts1(i,j) = sav0(i,kz*4+2,j)
      end do
    end do
    if ( ehso4 ) then
      call mpi_scatter(sav_0s,iy*kz*jxp,mpi_real8,  &
                       sav0s, iy*kz*jxp,mpi_real8,  &
                       0,mycomm,ierr)
      do j = 1 , jendl
        do k = 1 , kz
          do i = 1 , iy
            so1(i,k,j) = sav0s(i,k,j)
          end do
        end do
      end do
    end if
!
!   Convert surface pressure to pstar
!
    do j = 1 , jendl
      do i = 1 , iy
        xpsb%b1(i,j) = xpsb%b1(i,j) - ptop
      end do
    end do
!
!  this routine determines p(.) from p(x) by a 4-point
!  interpolation. on the x-grid, a p(x) point outside the grid
!  domain is assumed to satisfy p(0,j)=p(1,j); p(iy,j)=p(iym1,j);
!  and similarly for the i's.
!
    call mpi_sendrecv(xpsb%b1(:,jxp),iy,mpi_real8,ieast,1,   &
                      xpsb%b1(:,0),  iy,mpi_real8,iwest,1,   &
                      mycomm,mpi_status_ignore,ierr)
    do j = jbegin , jendx
      do i = 2 , iym1
        psdot(i,j) = d_rfour*(xpsb%b1(i,j)+xpsb%b1(i-1,j) + &
                              xpsb%b1(i,j-1)+xpsb%b1(i-1,j-1))
      end do
    end do
#ifdef BAND
    do j = jbegin , jendx
      psdot(1,j)  = d_half*(xpsb%b1(1,j)+xpsb%b1(1,j-1))
      psdot(iy,j) = d_half*(xpsb%b1(iym1,j)+xpsb%b1(iym1,j-1))
    end do
#else
!
    do i = 2 , iym1
      if ( myid == 0 ) then
        psdot(i,1) = d_half*(xpsb%b1(i,1)+xpsb%b1(i-1,1))
      end if
      if ( myid == nproc-1 ) then
        psdot(i,jendl) = d_half*(xpsb%b1(i,jendx)+xpsb%b1(i-1,jendx))
      end if
    end do
!
    do j = jbegin , jendx
      psdot(1,j)  = d_half*(xpsb%b1(1,j)+xpsb%b1(1,j-1))
      psdot(iy,j) = d_half*(xpsb%b1(iym1,j)+xpsb%b1(iym1,j-1))
    end do
!
    if ( myid == 0 ) then
      psdot(1,1)  = xpsb%b1(1,1)
      psdot(iy,1) = xpsb%b1(iym1,1)
    end if
    if ( myid == nproc-1 ) then
      psdot(1,jendl)  = xpsb%b1(1,jendx)
      psdot(iy,jendl) = xpsb%b1(iym1,jendx)
    end if
#endif
!
!   Couple pressure u,v,t,q
!
    do k = 1 , kz
      do j = 1 , jendl
        do i = 1 , iy
          xub%b1(i,k,j) = xub%b1(i,k,j)*psdot(i,j)
          xvb%b1(i,k,j) = xvb%b1(i,k,j)*psdot(i,j)
          xtb%b1(i,k,j) = xtb%b1(i,k,j)*xpsb%b1(i,j)
          xqb%b1(i,k,j) = xqb%b1(i,k,j)*xpsb%b1(i,j)
        end do
      end do
    end do
!
!   compute boundary conditions for p*:
!
#ifdef BAND
    nxwb = 0
    nxeb = 0
#else
    if ( nspgx <= jxp ) then
      nxwb = nspgx
    else
      nkk = nspgx/jxp
      if ( nspgx == nkk*jxp ) then
        nxwb = jxp
      else
        nxwb = nspgx - nkk*jxp
      end if
    end if
    if ( nxwb+myid*jxp > nspgx ) then
      nxwb = 0
    else if ( nxwb+myid*jxp < nspgx ) then
      nxwb = jxp
    end if
!
    if ( nspgx <= jxp-1 ) then
      nxeb = nspgx
    else
      nkk = (nspgx-jxp+1)/jxp
      if ( (nspgx-jxp+1) == nkk*jxp ) then
        nxeb = jxp
      else
        nxeb = (nspgx-jxp+1) - nkk*jxp
      end if
    end if
    if ( jxm1-(myid*jxp+jxp-nxeb) > nspgx ) then
      nxeb = 0
    else if ( jxm1-(myid*jxp+jxp-nxeb) < nspgx ) then
      nxeb = min0(jendx,jxp)
    end if
    do nn = 1 , nxwb
      do i = 1 , iym1
        xpsb%wb(i,nn) = xpsb%b0(i,nn)
        xpsb%wbt(i,nn) = (xpsb%b1(i,nn)-xpsb%b0(i,nn))/dtbdys
      end do
    end do
    do nn = 1 , nxeb
      nnb = min0(jendx,jxp) - nn + 1
      do i = 1 , iym1
        xpsb%eb(i,nn) = xpsb%b0(i,nnb)
        xpsb%ebt(i,nn) = (xpsb%b1(i,nnb)-xpsb%b0(i,nnb))/dtbdys
      end do
    end do
#endif
    do nn = 1 , nspgx
      nnb = iym1 - nn + 1
      do j = 1 , jendx
        xpsb%nb(nn,j) = xpsb%b0(nnb,j)
        xpsb%sb(nn,j) = xpsb%b0(nn,j)
        xpsb%nbt(nn,j) = (xpsb%b1(nnb,j)-xpsb%b0(nnb,j))/dtbdys
        xpsb%sbt(nn,j) = (xpsb%b1(nn,j)-xpsb%b0(nn,j))/dtbdys
      end do
    end do
!
!   compute boundary conditions for p*u and p*v:
!
#ifdef BAND
    ndwb = 0
    ndeb = 0
#else
    if ( nspgd <= jxp ) then
      ndwb = nspgd
    else
      nkk = nspgd/jxp
      if ( nspgd == nkk*jxp ) then
        ndwb = jxp
      else
        ndwb = nspgd - nkk*jxp
      end if
    end if
    if ( ndwb+myid*jxp > nspgd ) then
      ndwb = 0
    else if ( ndwb+myid*jxp < nspgd ) then
      ndwb = jxp
    end if
!
    if ( nspgd <= jendl ) then
      ndeb = nspgd
    else
      nkk = nspgd/jxp
      if ( nspgd == nkk*jxp ) then
        ndeb = jxp
      else
        ndeb = nspgd - nkk*jxp
      end if
    end if
    if ( jx-(myid*jxp+jxp-ndeb) > nspgd ) then
      ndeb = 0
    else if ( jx-(myid*jxp+jxp-ndeb) < nspgd ) then
      ndeb = jxp
    end if
    do nn = 1 , ndwb
      do k = 1 , kz
        do i = 1 , iy
          xub%wb(i,k,nn)  = xub%b0(i,k,nn)
          xvb%wb(i,k,nn)  = xvb%b0(i,k,nn)
          xub%wbt(i,k,nn) = (xub%b1(i,k,nn)-xub%b0(i,k,nn))/dtbdys
          xvb%wbt(i,k,nn) = (xvb%b1(i,k,nn)-xvb%b0(i,k,nn))/dtbdys
        end do
      end do
    end do
    do nn = 1 , ndeb
      nnb = min0(jendl,jxp) - nn + 1
      do k = 1 , kz
        do i = 1 , iy
          xub%eb(i,k,nn)  = xub%b0(i,k,nnb)
          xvb%eb(i,k,nn)  = xvb%b0(i,k,nnb)
          xub%ebt(i,k,nn) = (xub%b1(i,k,nnb)-xub%b0(i,k,nnb))/dtbdys
          xvb%ebt(i,k,nn) = (xvb%b1(i,k,nnb)-xvb%b0(i,k,nnb))/dtbdys
        end do
      end do
    end do
#endif
    do nn = 1 , nspgd
      nnb = iy - nn + 1
      do k = 1 , kz
        do j = 1 , jendl
          xub%nb(nn,k,j)  = xub%b0(nnb,k,j)
          xub%sb(nn,k,j)  = xub%b0(nn,k,j)
          xvb%nb(nn,k,j)  = xvb%b0(nnb,k,j)
          xvb%sb(nn,k,j)  = xvb%b0(nn,k,j)
          xub%nbt(nn,k,j) = (xub%b1(nnb,k,j)-xub%b0(nnb,k,j))/dtbdys
          xub%sbt(nn,k,j) = (xub%b1(nn,k,j)-xub%b0(nn,k,j))/dtbdys
          xvb%nbt(nn,k,j) = (xvb%b1(nnb,k,j)-xvb%b0(nnb,k,j))/dtbdys
          xvb%sbt(nn,k,j) = (xvb%b1(nn,k,j)-xvb%b0(nn,k,j))/dtbdys
        end do
      end do
    end do
!
!   compute boundary conditions for p*t and p*qv:
!
#ifndef BAND
    do nn = 1 , nxwb
      do k = 1 , kz
        do i = 1 , iym1
          xtb%wb(i,k,nn)  = xtb%b0(i,k,nn)
          xqb%wb(i,k,nn)  = xqb%b0(i,k,nn)
          xtb%wbt(i,k,nn) = (xtb%b1(i,k,nn)-xtb%b0(i,k,nn))/dtbdys
          xqb%wbt(i,k,nn) = (xqb%b1(i,k,nn)-xqb%b0(i,k,nn))/dtbdys
        end do
      end do
    end do
    do nn = 1 , nxeb
      nnb = min0(jendx,jxp) - nn + 1
      do k = 1 , kz
        do i = 1 , iym1
          xtb%eb(i,k,nn)  = xtb%b0(i,k,nnb)
          xqb%eb(i,k,nn)  = xqb%b0(i,k,nnb)
          xtb%ebt(i,k,nn) = (xtb%b1(i,k,nnb)-xtb%b0(i,k,nnb))/dtbdys
          xqb%ebt(i,k,nn) = (xqb%b1(i,k,nnb)-xqb%b0(i,k,nnb))/dtbdys
        end do
      end do
    end do
#endif
    do nn = 1 , nspgx
      nnb = iym1 - nn + 1
      do k = 1 , kz
        do j = 1 , jendx
          xtb%nb(nn,k,j)  = xtb%b0(nnb,k,j)
          xtb%sb(nn,k,j)  = xtb%b0(nn,k,j)
          xqb%nb(nn,k,j)  = xqb%b0(nnb,k,j)
          xqb%sb(nn,k,j)  = xqb%b0(nn,k,j)
          xtb%nbt(nn,k,j) = (xtb%b1(nnb,k,j)-xtb%b0(nnb,k,j))/dtbdys
          xtb%sbt(nn,k,j) = (xtb%b1(nn,k,j)-xtb%b0(nn,k,j))/dtbdys
          xqb%nbt(nn,k,j) = (xqb%b1(nnb,k,j)-xqb%b0(nnb,k,j))/dtbdys
          xqb%sbt(nn,k,j) = (xqb%b1(nn,k,j)-xqb%b0(nn,k,j))/dtbdys
        end do
      end do
    end do

    if ( myid == 0 ) then
      write (6,'(a,i10,a,i10)') 'READY  BC from     ' , &
            toint10(bdydate1) , ' to ' , toint10(bdydate2)
    end if

    bdydate1 = bdydate2
    call date_bcast(bdydate1,0,mycomm,ierr)

    do j = 1 , jendx
      do i = 1 , iym1
        tdum(i,j) = ts1(i,j)
      end do
    end do
    do k = 1 , kz
      do j = 1 , jendl
        do i = 1 , iy
          xub%b0(i,k,j) = xub%b1(i,k,j)
          xvb%b0(i,k,j) = xvb%b1(i,k,j)
          xqb%b0(i,k,j) = xqb%b1(i,k,j)
          xtb%b0(i,k,j) = xtb%b1(i,k,j)
        end do
      end do
    end do
    do j = 1 , jendl
      do i = 1 , iy
        xpsb%b0(i,j) = xpsb%b1(i,j)
        ts0(i,j) = ts1(i,j)
      end do
    end do
    if ( ehso4 ) then
      do k = 1 , kz
        do j = 1 , jendl
          do i = 1 , iy
            so0(i,k,j) = so1(i,k,j)
          end do
        end do
      end do
    end if
   
    if ( idatex < bdydate1 ) then
      do j = 1 , jendx
        do i = 1 , iym1
          if ( iswater(mddom%lndcat(i,j)) ) then
            if (idcsst == 1) then
              sts1%tg(i,j) = tdum(i,j) + dtskin(i,j)
              sts2%tg(i,j) = tdum(i,j) + dtskin(i,j)
            else
              sts1%tg(i,j) = tdum(i,j)
              sts2%tg(i,j) = tdum(i,j)
            end if
            if ( iseaice == 1 ) then
              if ( lakemod == 1 .and. islake(mddom%lndcat(i,j)) ) cycle
!
              if ( tdum(i,j) <= icetemp ) then
                sts1%tg(i,j) = icetemp
                sts2%tg(i,j) = icetemp
                tdum(i,j) = icetemp
                ldmsk(i,j) = 2
                do n = 1, nnsg
                  ocld2d(n,i,j) = 2
                  sice2d(n,i,j) = d_1000
                  scv2d(n,i,j) = d_zero
                end do
              else
                sts1%tg(i,j) = tdum(i,j)
                sts2%tg(i,j) = tdum(i,j)
                ldmsk(i,j) = 0
                do n = 1, nnsg
                  ocld2d(n,i,j) = 0
                  sice2d(n,i,j) = d_zero
                  scv2d(n,i,j)  = d_zero
                end do
              end if
            end if
          end if
        end do
      end do
    end if
    call time_end(subroutine_name,idindx)
  end subroutine bdyin
!
! This subroutine sets the boundary values of u and v according
! to the boundary conditions specified.
!
!     ib = 0 : fixed
!        = 1 : relaxation, linear technique
!        = 2 : time dependent
!        = 3 : time dependent and inflow/outflow dependent
!        = 4 : sponge
!        = 5 : relaxation, exponential technique
!
!     dtb    : elapsed time from the initial boundary values.
!
  subroutine bdyuv(ib,dtb)
!
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
!
    real(8) :: dtb
    integer :: ib
    intent (in) dtb , ib
!
    integer :: i , j , k
    integer :: ierr
    character (len=64) :: subroutine_name='bdyuv'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
!   compute the p* at dot points:
!
    call mpi_sendrecv(sps1%ps(1,jxp),iy,mpi_real8,ieast,1,   &
                      sps1%ps(1,0),  iy,mpi_real8,iwest,1,   &
                      mycomm,mpi_status_ignore,ierr)
!
!   interior points:
!
    do j = jbegin , jendx
      do i = 2 , iym1
        sps1%pdot(i,j) = d_rfour*(sps1%ps(i,j)  +sps1%ps(i-1,j) + &
                                  sps1%ps(i,j-1)+sps1%ps(i-1,j-1))
      end do
    end do
!
!   east and west boundaries:
!
#ifndef BAND
    do i = 2 , iym1
      if ( myid == 0 ) then
        sps1%pdot(i,1) = d_half*(sps1%ps(i,1)+sps1%ps(i-1,1))
      end if
      if ( myid == nproc-1 ) then
        sps1%pdot(i,jendl) = d_half*(sps1%ps(i,jendx)+sps1%ps(i-1,jendx))
      end if
    end do
#endif
!
!   north and south boundaries:
!
    do j = jbegin , jendx
      sps1%pdot(1,j)  = d_half*(sps1%ps(1,j)   +sps1%ps(1,j-1))
      sps1%pdot(iy,j) = d_half*(sps1%ps(iym1,j)+sps1%ps(iym1,j-1))
    end do
!
!   corner points:
!
#ifndef BAND
    if ( myid == 0 ) then
      sps1%pdot(1,1)  = sps1%ps(1,1)
      sps1%pdot(iy,1) = sps1%ps(iym1,1)
    end if
    if ( myid == nproc-1 ) then
      sps1%pdot(1,jendl)  = sps1%ps(1,jendx)
      sps1%pdot(iy,jendl) = sps1%ps(iym1,jendx)
    end if
#endif
!
!   interior silces:
!
    do k = 1 , kz
#ifndef BAND
!
      do i = 2 , iym1
        if ( myid == 0 ) then
          uj2(i,k) = atm1%u(i,k,2)/sps1%pdot(i,2)
          vj2(i,k) = atm1%v(i,k,2)/sps1%pdot(i,2)
        end if
        if ( myid == nproc-1 ) then
          ujlx(i,k) = atm1%u(i,k,jendx)/sps1%pdot(i,jendx)
          vjlx(i,k) = atm1%v(i,k,jendx)/sps1%pdot(i,jendx)
        end if
      end do
#endif
!
#ifdef BAND
      do j = 1 , jendl
#else
      do j = jbegin , jendx
#endif
        ui2(k,j)  = atm1%u(2,k,j)/sps1%pdot(2,j)
        vi2(k,j)  = atm1%v(2,k,j)/sps1%pdot(2,j)
        uilx(k,j) = atm1%u(iym1,k,j)/sps1%pdot(iym1,j)
        vilx(k,j) = atm1%v(iym1,k,j)/sps1%pdot(iym1,j)
      end do
    end do
!
!   boundary silces:
!
    if ( ib == 0 ) then
!
!     fixed boundary conditions:
!
      do k = 1 , kz
#ifndef BAND
!
!       west (j = 1) and east (j = jx) boundaries:
!
        do i = 1 , iy
          if ( myid == 0 ) then
            uj1(i,k) = xub%wb(i,k,1)/sps1%pdot(i,1)
            vj1(i,k) = xvb%wb(i,k,1)/sps1%pdot(i,1)
          end if
          if ( myid == nproc-1 ) then
            ujl(i,k) = xub%eb(i,k,1)/sps1%pdot(i,jendl)
            vjl(i,k) = xvb%eb(i,k,1)/sps1%pdot(i,jendl)
          end if
        end do
#endif
!
!       south (i = 1) and north (i = iy) boundaries:
!
        do j = 1 , jendl
          ui1(k,j) = xub%sb(1,k,j)/sps1%pdot(1,j)
          vi1(k,j) = xvb%sb(1,k,j)/sps1%pdot(1,j)
          uil(k,j) = xub%nb(1,k,j)/sps1%pdot(iy,j)
          vil(k,j) = xvb%nb(1,k,j)/sps1%pdot(iy,j)
        end do
      end do
    else
!
!     time-dependent boundary conditions:
!
      do k = 1 , kz
#ifndef BAND
!
!       west (j = 1) and east (j = jx) boundaries:
!
        do i = 1 , iy
          if ( myid == 0 ) then
            uj1(i,k) = (xub%wb(i,k,1)+dtb*xub%wbt(i,k,1))/sps1%pdot(i,1)
            vj1(i,k) = (xvb%wb(i,k,1)+dtb*xvb%wbt(i,k,1))/sps1%pdot(i,1)
          end if
          if ( myid == nproc-1 ) then
            ujl(i,k) = (xub%eb(i,k,1)+dtb*xub%ebt(i,k,1))/sps1%pdot(i,jendl)
            vjl(i,k) = (xvb%eb(i,k,1)+dtb*xvb%ebt(i,k,1))/sps1%pdot(i,jendl)
          end if
        end do
#endif
!
!       south (i = 1) and north (i = iy) boundaries:
!
        do j = 1 , jendl
          ui1(k,j) = (xub%sb(1,k,j)+dtb*xub%sbt(1,k,j))/sps1%pdot(1,j)
          vi1(k,j) = (xvb%sb(1,k,j)+dtb*xvb%sbt(1,k,j))/sps1%pdot(1,j)
          uil(k,j) = (xub%nb(1,k,j)+dtb*xub%nbt(1,k,j))/sps1%pdot(iy,j)
          vil(k,j) = (xvb%nb(1,k,j)+dtb*xvb%nbt(1,k,j))/sps1%pdot(iy,j)
        end do
      end do
!
    end if
!
!   fill up the interior silces:
!
#ifndef BAND
    do k = 1 , kz
      if ( myid == 0 ) then
        uj2(1,k)  = ui1(k,2)
        uj2(iy,k) = uil(k,2)
        ui2(k,1)  = uj1(2,k)
        uilx(k,1) = uj1(iym1,k)
        vj2(1,k)  = vi1(k,2)
        vj2(iy,k) = vil(k,2)
        vi2(k,1)  = vj1(2,k)
        vilx(k,1) = vj1(iym1,k)
      end if
      if ( myid == nproc-1 ) then
        ujlx(1,k)     = ui1(k,jendx)
        ujlx(iy,k)    = uil(k,jendx)
        ui2(k,jendl)  = ujl(2,k)
        uilx(k,jendl) = ujl(iym1,k)
        vjlx(1,k)     = vi1(k,jendx)
        vjlx(iy,k)    = vil(k,jendx)
        vi2(k,jendl)  = vjl(2,k)
        vilx(k,jendl) = vjl(iym1,k)
      end if
    end do
#endif
#ifndef BAND
    if ( myid /= nproc-1 ) then
#endif
      do k = 1 , kz
        var1snd(k,1) = ui1(k,jxp)
        var1snd(k,2) = vi1(k,jxp)
        var1snd(k,3) = ui2(k,jxp)
        var1snd(k,4) = vi2(k,jxp)
        var1snd(k,5) = uilx(k,jxp)
        var1snd(k,6) = vilx(k,jxp)
        var1snd(k,7) = uil(k,jxp)
        var1snd(k,8) = vil(k,jxp)
      end do
#ifndef BAND
    end if
#endif
    call mpi_sendrecv(var1snd(1,1),kz*8,mpi_real8,ieast,1,    &
                      var1rcv(1,1),kz*8,mpi_real8,iwest,1,    &
                      mycomm,mpi_status_ignore,ierr)
#ifndef BAND
    if ( myid /= 0 ) then
#endif
      do k = 1 , kz
        ui1(k,0)  = var1rcv(k,1)
        vi1(k,0)  = var1rcv(k,2)
        ui2(k,0)  = var1rcv(k,3)
        vi2(k,0)  = var1rcv(k,4)
        uilx(k,0) = var1rcv(k,5)
        vilx(k,0) = var1rcv(k,6)
        uil(k,0)  = var1rcv(k,7)
        vil(k,0)  = var1rcv(k,8)
      end do
#ifndef BAND
    end if
#endif
!
#ifndef BAND
    if ( myid /= 0 ) then
#endif
      do k = 1 , kz
        var1snd(k,1) = ui1(k,1)
        var1snd(k,2) = vi1(k,1)
        var1snd(k,3) = ui2(k,1)
        var1snd(k,4) = vi2(k,1)
        var1snd(k,5) = uilx(k,1)
        var1snd(k,6) = vilx(k,1)
        var1snd(k,7) = uil(k,1)
        var1snd(k,8) = vil(k,1)
      end do
#ifndef BAND
    end if
#endif
    call mpi_sendrecv(var1snd(1,1),kz*8,mpi_real8,iwest,2,     &
                      var1rcv(1,1),kz*8,mpi_real8,ieast,2,     &
                      mycomm,mpi_status_ignore,ierr)
#ifndef BAND
    if ( myid /= nproc-1 ) then
#endif
      do k = 1 , kz
        ui1(k,jxp+1)  = var1rcv(k,1)
        vi1(k,jxp+1)  = var1rcv(k,2)
        ui2(k,jxp+1)  = var1rcv(k,3)
        vi2(k,jxp+1)  = var1rcv(k,4)
        uilx(k,jxp+1) = var1rcv(k,5)
        vilx(k,jxp+1) = var1rcv(k,6)
        uil(k,jxp+1)  = var1rcv(k,7)
        vil(k,jxp+1)  = var1rcv(k,8)
      end do
#ifndef BAND
    end if
#endif

    call time_end(subroutine_name,idindx)
  end subroutine bdyuv
!
! This subroutine sets the boundary values for p*, p*u, p*v,
! p*t, p*qv, p*qc, and p*qr.
!
!     ---the boundary values of p*u and p*v are extrapolated from
!        the interior points.
!
!     ---the boundary values of p* and p*t are specified.
!
!     ---the boundary values of p*qv, p*qc, and p*qr depend on
!        inflow/outflow conditions, if iboudy = 3 or 4.
!
!     xt     : is the time in seconds the variables xxa represent.
!
!     iexec  : = 1 ; represents this subroutine is called for the
!                    first time in this forecast run.
!              > 1 ; represents subsequent calls to this subroutine.
!
  subroutine bdyval(xt,iexec)
!
    implicit none
!
    integer :: iexec
    real(8) :: xt
    intent (in) iexec , xt
!
    real(8) :: chix , chix1 , chix2 , dtb , qcx , qcx2 , qvx , qvx1 , qvx2 , vavg
    integer :: itr , j , k
#ifndef BAND
    integer :: i
    real(8) :: uavg
#endif
    character (len=64) :: subroutine_name='bdyval'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
!   Fill up the boundary value for xxb variables from xxa variables:
!   if this subroutine is called for the first time, this part
!   shall be skipped.
!
    if ( iexec /= 1 ) then
!
!     for p*:
!
#ifndef BAND
      do i = 1 , iym1
        if ( myid == 0 )       sps2%ps(i,1)     = sps1%ps(i,1)
        if ( myid == nproc-1 ) sps2%ps(i,jendx) = sps1%ps(i,jendx)
      end do
#endif
      do j = jbegin , jendm
        sps2%ps(1,j)    = sps1%ps(1,j)
        sps2%ps(iym1,j) = sps1%ps(iym1,j)
      end do
!
!     for p*u and p*v:
!
      do k = 1 , kz
#ifndef BAND
        do i = 1 , iy
          if ( myid == 0 ) then
            atm2%u(i,k,1) = atm1%u(i,k,1)/mddom%msfd(i,1)
            atm2%v(i,k,1) = atm1%v(i,k,1)/mddom%msfd(i,1)
          end if
          if ( myid == nproc-1 ) then
            atm2%u(i,k,jendl) = atm1%u(i,k,jendl)/mddom%msfd(i,jendl)
            atm2%v(i,k,jendl) = atm1%v(i,k,jendl)/mddom%msfd(i,jendl)
          end if
        end do
#endif
        do j = jbegin , jendx
          atm2%u(1,k,j)  = atm1%u(1,k,j)/mddom%msfd(1,j)
          atm2%u(iy,k,j) = atm1%u(iy,k,j)/mddom%msfd(iy,j)
          atm2%v(1,k,j)  = atm1%v(1,k,j)/mddom%msfd(1,j)
          atm2%v(iy,k,j) = atm1%v(iy,k,j)/mddom%msfd(iy,j)
        end do
      end do
!
!     for p*t:
!
      do k = 1 , kz
#ifndef BAND
        do i = 1 , iym1
          if ( myid == 0 )       atm2%t(i,k,1)     = atm1%t(i,k,1)
          if ( myid == nproc-1 ) atm2%t(i,k,jendx) = atm1%t(i,k,jendx)
        end do
#endif
        do j = jbegin , jendm
          atm2%t(1,k,j)    = atm1%t(1,k,j)
          atm2%t(iym1,k,j) = atm1%t(iym1,k,j)
        end do
      end do
!
!     for p*qv:
!
      do k = 1 , kz
#ifndef BAND
        do i = 1 , iym1
          if ( myid == 0 )       atm2%qv(i,k,1)     = atm1%qv(i,k,1)
          if ( myid == nproc-1 ) atm2%qv(i,k,jendx) = atm1%qv(i,k,jendx)
        end do
#endif
        do j = jbegin , jendm
          atm2%qv(1,k,j)    = atm1%qv(1,k,j)
          atm2%qv(iym1,k,j) = atm1%qv(iym1,k,j)
        end do
      end do
!
      if ( ichem == 1 ) then
!
!       for p*chi (tracers)
!
        do itr = 1 , ntr
          do k = 1 , kz
#ifndef BAND
            do i = 1 , iym1
              if ( myid == 0 )       chib(i,k,1,itr)     = chia(i,k,1,itr)
              if ( myid == nproc-1 ) chib(i,k,jendx,itr) = chia(i,k,jendx,itr)
            end do
#endif
            do j = jbegin , jendm
              chib(1,k,j,itr)    = chia(1,k,j,itr)
              chib(iym1,k,j,itr) = chia(iym1,k,j,itr)
            end do
          end do
        end do
      end if
!
!     for p*qc:
!
      do k = 1 , kz
#ifndef BAND
        do i = 1 , iym1
          if ( myid == 0 )       atm2%qc(i,k,1)     = atm1%qc(i,k,1)
          if ( myid == nproc-1 ) atm2%qc(i,k,jendx) = atm1%qc(i,k,jendx)
        end do
#endif
        do j = jbegin , jendm
          atm2%qc(1,k,j)    = atm1%qc(1,k,j)
          atm2%qc(iym1,k,j) = atm1%qc(iym1,k,j)
        end do
      end do
!
    end if  !end if (iexec /= 1) test
!
!   Compute the boundary values for xxa variables:
!
!   compute the time interval for boundary tendency:
!
    dtb = xt
    if ( nbdytime == 0 .and. ktau > 0 ) then
      dtb = dtbdys
    end if
!
!   set boundary values for p*:
!   set boundary conditions for p*u and p*v:
!
    if ( .not.(iexec == 1 .and. ifrest) ) then
!
      if ( iboudy == 0 ) then
!
!       fixed boundary conditions:
!
#ifndef BAND
        do i = 1 , iym1
          if ( myid == 0 )       sps1%ps(i,1)     = xpsb%wb(i,1)
          if ( myid == nproc-1 ) sps1%ps(i,jendx) = xpsb%eb(i,1)
        end do
#endif
        do j = jbegin , jendm
          sps1%ps(1,j)    = xpsb%sb(1,j)
          sps1%ps(iym1,j) = xpsb%nb(1,j)
        end do
!
        do k = 1 , kz
#ifndef BAND
          do i = 1 , iy
            if ( myid == 0 ) then
              atm1%u(i,k,1) = xub%wb(i,k,1)
              atm1%v(i,k,1) = xvb%wb(i,k,1)
            end if
            if ( myid == nproc-1 ) then
              atm1%u(i,k,jendl) = xub%eb(i,k,1)
              atm1%v(i,k,jendl) = xvb%eb(i,k,1)
            end if
          end do
#endif
          do j = jbegin , jendx
            atm1%u(1,k,j)  = xub%sb(1,k,j)
            atm1%u(iy,k,j) = xub%nb(1,k,j)
            atm1%v(1,k,j)  = xvb%sb(1,k,j)
            atm1%v(iy,k,j) = xvb%nb(1,k,j)
          end do
        end do
      end if
!
!     time-dependent boundary conditions:
!
#ifndef BAND
      do i = 1 , iym1
        if ( myid == 0 ) then
          sps1%ps(i,1) = xpsb%wb(i,1) + dtb*xpsb%wbt(i,1)
        end if
        if ( myid == nproc-1 ) then
          sps1%ps(i,jendx) = xpsb%eb(i,1) + dtb*xpsb%ebt(i,1)
        end if
      end do
#endif
      do j = jbegin , jendm
        sps1%ps(1,j)    = xpsb%sb(1,j) + dtb*xpsb%sbt(1,j)
        sps1%ps(iym1,j) = xpsb%nb(1,j) + dtb*xpsb%nbt(1,j)
      end do
!
      do k = 1 , kz
#ifndef BAND
        do i = 1 , iy
          if ( myid == 0 ) then
            atm1%u(i,k,1) = xub%wb(i,k,1) + dtb*xub%wbt(i,k,1)
            atm1%v(i,k,1) = xvb%wb(i,k,1) + dtb*xvb%wbt(i,k,1)
          end if
          if ( myid == nproc-1 ) then
            atm1%u(i,k,jendl) = xub%eb(i,k,1) + dtb*xub%ebt(i,k,1)
            atm1%v(i,k,jendl) = xvb%eb(i,k,1) + dtb*xvb%ebt(i,k,1)
          end if
        end do
#endif
        do j = jbegin , jendx
          atm1%u(1,k,j)  = xub%sb(1,k,j) + dtb*xub%sbt(1,k,j)
          atm1%u(iy,k,j) = xub%nb(1,k,j) + dtb*xub%nbt(1,k,j)
          atm1%v(1,k,j)  = xvb%sb(1,k,j) + dtb*xvb%sbt(1,k,j)
          atm1%v(iy,k,j) = xvb%nb(1,k,j) + dtb*xvb%nbt(1,k,j)
        end do
      end do
    end if
!
!   get boundary values of u and v:
!
    call bdyuv(iboudy,dtb)
!
    if ( iexec == 1 .and. ifrest ) return
!
!     set boundary values for p*t:
!     set boundary values for p*qv:
!
    if ( iboudy == 0 ) then
!
!     fixed boundary conditions:
!
      do k = 1 , kz
#ifndef BAND
        do i = 1 , iym1
          if ( myid == 0 )       atm1%t(i,k,1)     = xtb%wb(i,k,1)
          if ( myid == nproc-1 ) atm1%t(i,k,jendx) = xtb%eb(i,k,1)
        end do
#endif
        do j = jbegin , jendm
          atm1%t(1,k,j)    = xtb%sb(1,k,j)
          atm1%t(iym1,k,j) = xtb%nb(1,k,j)
        end do
      end do
      do k = 1 , kz
#ifndef BAND
        do i = 1 , iym1
          if ( myid == 0 )       atm1%qv(i,k,1)     = xqb%wb(i,k,1)
          if ( myid == nproc-1 ) atm1%qv(i,k,jendx) = xqb%eb(i,k,1)
        end do
#endif
        do j = jbegin , jendm
          atm1%qv(1,k,j)    = xqb%sb(1,k,j)
          atm1%qv(iym1,k,j) = xqb%nb(1,k,j)
        end do
      end do
    end if
!
!   time-dependent boundary conditions:
!
    do k = 1 , kz
#ifndef BAND
      do i = 1 , iym1
        if ( myid == 0 ) then
          atm1%t(i,k,1) = xtb%wb(i,k,1) + dtb*xtb%wbt(i,k,1)
        end if
        if ( myid == nproc-1 ) then
          atm1%t(i,k,jendx) = xtb%eb(i,k,1) + dtb*xtb%ebt(i,k,1)
        end if
      end do
#endif
      do j = jbegin , jendm
        atm1%t(1,k,j)    = xtb%sb(1,k,j) + dtb*xtb%sbt(1,k,j)
        atm1%t(iym1,k,j) = xtb%nb(1,k,j) + dtb*xtb%nbt(1,k,j)
      end do
    end do
    do k = 1 , kz
#ifndef BAND
      do i = 1 , iym1
        if ( myid == 0 ) then
          atm1%qv(i,k,1) = xqb%wb(i,k,1) + dtb*xqb%wbt(i,k,1)
        end if
        if ( myid == nproc-1 ) then
          atm1%qv(i,k,jendx) = xqb%eb(i,k,1) + dtb*xqb%ebt(i,k,1)
        end if
      end do
#endif
      do j = jbegin , jendm
        atm1%qv(1,k,j)    = xqb%sb(1,k,j) + dtb*xqb%sbt(1,k,j)
        atm1%qv(iym1,k,j) = xqb%nb(1,k,j) + dtb*xqb%nbt(1,k,j)
      end do
    end do
!
    if ( iboudy == 3 .or. iboudy == 4 ) then
!
!     determine boundary values depends on inflow/outflow:
!
      do k = 1 , kz
#ifndef BAND
!
!       west boundary:
!
        if ( myid == 0 ) then
          do i = 1 , iym1
            qvx1 = atm1%qv(i,k,1)/sps1%ps(i,1)
            qvx2 = atm1%qv(i,k,2)/sps1%ps(i,2)
            uavg = uj1(i,k) + uj1(i+1,k) + uj2(i,k) + uj2(i+1,k)
            if ( uavg >= d_zero ) then
              qvx = qvx1
            else
              qvx = qvx2
            end if
            atm1%qv(i,k,1) = qvx*sps1%ps(i,1)
          end do
        end if
!
!       east boundary:
!
        if ( myid == nproc-1 ) then
          do i = 1 , iym1
            qvx1 = atm1%qv(i,k,jendx)/sps1%ps(i,jendx)
            qvx2 = atm1%qv(i,k,jendm)/sps1%ps(i,jendm)
            uavg = ujlx(i,k) + ujlx(i+1,k) + ujl(i,k) + ujl(i+1,k)
            if ( uavg < d_zero ) then
              qvx = qvx1
            else
              qvx = qvx2
            end if
            atm1%qv(i,k,jendx) = qvx*sps1%ps(i,jendx)
          end do
        end if
#endif
!
!       south boundary:
!
        do j = jbegin , jendm
          qvx1 = atm1%qv(1,k,j)/sps1%ps(1,j)
          qvx2 = atm1%qv(2,k,j)/sps1%ps(2,j)
          vavg = vi1(k,j) + vi1(k,j+1) + vi2(k,j) + vi2(k,j+1)
          if ( vavg >= d_zero ) then
            qvx = qvx1
          else
            qvx = qvx2
          end if
          atm1%qv(1,k,j) = qvx*sps1%ps(1,j)
        end do
!
!       north boundary:
!
        do j = jbegin , jendm
          qvx1 = atm1%qv(iym1,k,j)/sps1%ps(iym1,j)
          qvx2 = atm1%qv(iym2,k,j)/sps1%ps(iym2,j)
          vavg = vilx(k,j) + vilx(k,j+1) + vil(k,j) + vil(k,j+1)
          if ( vavg < d_zero ) then
            qvx = qvx1
          else
            qvx = qvx2
          end if
          atm1%qv(iym1,k,j) = qvx*sps1%ps(iym1,j)
        end do
!
      end do
    end if !end if (iboudy == 3.or.4) test
!
!   set boundary values for p*qc and p*qr:
!   *** note ***
!   for large domain, we assume the boundary tendencies are not available.
!
!   if the boundary values and tendencies are not available,
!   determine boundary values depends on inflow/outflow:
!   inflow  : set it equal to zero.
!   outflow : get from interior point.
!
    do k = 1 , kz
#ifndef BAND
!
!     west boundary:
!
      if ( myid == 0 ) then
        do i = 1 , iym1
          qcx2 = atm1%qc(i,k,2)/sps1%ps(i,2)
          uavg = uj1(i,k) + uj1(i+1,k) + uj2(i,k) + uj2(i+1,k)
          if ( uavg >= d_zero ) then
            qcx = d_zero
          else
            qcx = qcx2
          end if
          atm1%qc(i,k,1) = qcx*sps1%ps(i,1)
        end do
      end if
!
!     east boundary:
!
      if ( myid == nproc-1 ) then
        do i = 1 , iym1
          qcx2 = atm1%qc(i,k,jendm)/sps1%ps(i,jendm)
          uavg = ujlx(i,k) + ujlx(i+1,k) + ujl(i,k) + ujl(i+1,k)
          if ( uavg < d_zero ) then
            qcx = d_zero
          else
            qcx = qcx2
          end if
          atm1%qc(i,k,jendx) = qcx*sps1%ps(i,jendx)
        end do
      end if
#endif
!
!     south boundary:
!
      do j = jbegin , jendm
        qcx2 = atm1%qc(2,k,j)/sps1%ps(2,j)
        vavg = vi1(k,j) + vi1(k,j+1) + vi2(k,j) + vi2(k,j+1)
        if ( vavg >= d_zero ) then
          qcx = d_zero
        else
          qcx = qcx2
        end if
        atm1%qc(1,k,j) = qcx*sps1%ps(1,j)
      end do
!
!     north boundary:
!
      do j = jbegin , jendm
        qcx2 = atm1%qc(iym2,k,j)/sps1%ps(iym2,j)
        vavg = vilx(k,j) + vilx(k,j+1) + vil(k,j) + vil(k,j+1)
        if ( vavg < d_zero ) then
          qcx = d_zero
        else
          qcx = qcx2
        end if
        atm1%qc(iym1,k,j) = qcx*sps1%ps(iym1,j)
      end do
!
    end do
   
    if ( ichem == 1 ) then
!   
!     add tracer bc's
!
      do itr = 1 , ntr
        do k = 1 , kz
#ifndef BAND
!
!         west  boundary:
!
          if ( myid == 0 ) then
            do i = 1 , iym1
              chix1 = chia(i,k,1,itr)/sps1%ps(i,1)
              chix2 = chia(i,k,2,itr)/sps1%ps(i,2)
              uavg = uj1(i,k) + uj1(i+1,k) + uj2(i,k) + uj2(i+1,k)
              if ( uavg >= d_zero ) then
                chix = chix1
              else
                chix = chix2
              end if
              chia(i,k,1,itr) = chix*sps1%ps(i,1)
            end do
          end if
!
!         east  boundary:
!
          if ( myid == nproc-1 ) then
            do i = 1 , iym1
              chix1 = chia(i,k,jendx,itr)/sps1%ps(i,jendx)
              chix2 = chia(i,k,jendm,itr)/sps1%ps(i,jendm)
              uavg = ujlx(i,k) + ujlx(i+1,k) + ujl(i,k) + ujl(i+1,k)
              if ( uavg < d_zero ) then
                chix = chix1
              else
                chix = chix2
              end if
              chia(i,k,jendx,itr) = chix*sps1%ps(i,jendx)
            end do
          end if
#endif
!
!         south boundary:
!
          do j = jbegin , jendm
            chix1 = chia(1,k,j,itr)/sps1%ps(1,j)
            chix2 = chia(2,k,j,itr)/sps1%ps(2,j)
            vavg = vi1(k,j) + vi1(k,j+1) + vi2(k,j) + vi2(k,j+1)
            if ( vavg >= d_zero ) then
              chix = chix1
            else
              chix = chix2
            end if
            chia(1,k,j,itr) = chix*sps1%ps(1,j)
          end do
!
!         north boundary:
!
          do j = jbegin , jendm
            chix1 = chia(iym1,k,j,itr)/sps1%ps(iym1,j)
            chix2 = chia(iym2,k,j,itr)/sps1%ps(iym2,j)
            vavg = vilx(k,j) + vilx(k,j+1) + vil(k,j) + vil(k,j+1)
            if ( vavg < d_zero ) then
              chix = chix1
            else
              chix = chix2
            end if
            chia(iym1,k,j,itr) = chix*sps1%ps(iym1,j)
          end do
        end do
      end do
    end if
!
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      call set_tke_bc(atm1,atm2)
    end if
!
    call time_end(subroutine_name,idindx)
!
  end subroutine bdyval
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
! this subroutine applies sponge boundary condition to the        c
! tendency term - ften.                                           c
!                                                                 c
! ldot  : logical dot (u,v) / cross (t,q,p) flag                  c
!                                                                 c
! ip    : is the number of slices affected by nudging.            c
!                                                                 c
! wg   : are the weightings.                                      c
!                                                                 c
! ften  : is the tendency calculated from the model.              c
!                                                                 c
! j     : is the j'th slice of the tendency to be adjusted.       c
!                                                                 c
! nk    : is the number of vertical level to be adjusted.         c
!                                                                 c
! bnd   : Boundary condition data structure                       c
!         2D or 3D (managed by interface declaration)             c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine sponge3d(ldot,ip,wg,ften,j,nk,bnd)
!
    implicit none
!
    logical , intent(in) :: ldot
    integer , intent(in) :: ip , j , nk
    real(8) , intent(in) , dimension(ip) :: wg
    type(v3dbound) , intent(in) :: bnd
    real(8) , intent(inout) , dimension(iy,nk,jxp) :: ften
!
    integer :: i , ii , k , ido
#ifndef BAND
    integer :: ibeg , iend , jj , jsls
#endif
#ifndef BAND
    integer :: jwb , jeb , jew
#endif
    character (len=64) :: subroutine_name='sponge3d'
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
!-----interior j slices:
    do i = 2 , ip
       ii = iy - i + ido
       do k = 1 , nk
!.......south boundary:
          ften(i,k,j) = wg(i)*ften(i,k,j) + (d_one-wg(i))*bnd%sbt(i,k,j)
!.......north boundary:
          ften(ii,k,j) = wg(i)*ften(ii,k,j) + (d_one-wg(i))*bnd%nbt(i,k,j)
       end do
    end do

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
    if ( jsls > ip ) then
!-----interior j slices:
      do i = 2 , ip
        ii = iy - i + ido
        do k = 1 , nk
!.......south boundary:
          ften(i,k,j) = wg(i)*ften(i,k,j) + (d_one-wg(i))*bnd%sbt(i,k,j)
!.......north boundary:
          ften(ii,k,j) = wg(i)*ften(ii,k,j) + (d_one-wg(i))*bnd%nbt(i,k,j)
        end do
      end do
!
    else if ( jsls <= ip ) then
      ibeg = 2
      iend = iym1 - 1 + ido
      if ( jsls > 2 ) then
        do i = 2 , jsls - 1
          ii = iy - i + ido
          do k = 1 , nk
!........south boundary:
            ften(i,k,j) = wg(i)*ften(i,k,j) + (d_one-wg(i))*bnd%sbt(i,k,j)
!........north boundary:
            ften(ii,k,j) = wg(i)*ften(ii,k,j) + (d_one-wg(i))*bnd%nbt(i,k,j)
          end do
        end do
        ibeg = jsls
        iend = iy - jsls + ido
      end if
!
      if ( jj > ip ) then
!------west-boundary slice:
        do k = 1 , nk
          do i = ibeg , iend
            if ( jsls <= ip ) then
              ften(i,k,j) = wg(jsls)*ften(i,k,j) + &
                            (d_one-wg(jsls))*bnd%wbt(i,k,jwb)
            end if
          end do
        end do
      else if ( jj <= ip ) then
!------east-boundary slice:
        do k = 1 , nk
          do i = ibeg , iend
            if ( jsls <= ip ) then
              ften(i,k,j) = wg(jsls)*ften(i,k,j) + &
                            (d_one-wg(jsls))*bnd%ebt(i,k,jeb)
            end if
          end do
        end do
      end if
!
    end if
#endif

    call time_end(subroutine_name,idindx)
  end subroutine sponge3d
!
  subroutine sponge2d(ldot,ip,wg,ften,j,nk,bnd)
!
    implicit none
!
    logical , intent(in) :: ldot
    integer , intent(in) :: ip , j , nk
    real(8) , intent(in) , dimension(ip) :: wg
    type(v2dbound) , intent(in) :: bnd
    real(8) , intent(inout) , dimension(iy,jxp) :: ften
!
    integer :: i , ii , k , ido
#ifndef BAND
    integer :: ibeg , iend , jj , jsls
#endif
#ifndef BAND
    integer :: jwb , jeb , jew
#endif
    character (len=64) :: subroutine_name='sponge2d'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
!----------------------------------------------------------------------
!
    k = nk ! not to have warnings
    ido = 0
    if ( ldot ) then
      ido = 1
    end if
#ifdef BAND
!
!-----interior j slices:
    do i = 2 , ip
       ii = iy - i + ido
!.......south boundary:
        ften(i,j) = wg(i)*ften(i,j) + (d_one-wg(i))*bnd%sbt(i,j)
!.......north boundary:
        ften(ii,j) = wg(i)*ften(ii,j) + (d_one-wg(i))*bnd%nbt(i,j)
    end do

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
    if ( jsls > ip ) then
!-----interior j slices:
      do i = 2 , ip
        ii = iy - i + ido
!.......south boundary:
        ften(i,j) = wg(i)*ften(i,j) + (d_one-wg(i))*bnd%sbt(i,j)
!.......north boundary:
        ften(ii,j) = wg(i)*ften(ii,j) + (d_one-wg(i))*bnd%nbt(i,j)
      end do
!
    else if ( jsls <= ip ) then
      ibeg = 2
      iend = iym1 - 1 + ido
      if ( jsls > 2 ) then
        do i = 2 , jsls - 1
          ii = iy - i + ido
!........south boundary:
          ften(i,j) = wg(i)*ften(i,j) + (d_one-wg(i))*bnd%sbt(i,j)
!........north boundary:
          ften(ii,j) = wg(i)*ften(ii,j) + (d_one-wg(i))*bnd%nbt(i,j)
        end do
        ibeg = jsls
        iend = iy - jsls + ido
      end if
!
      if ( jj > ip ) then
!------west-boundary slice:
        do i = ibeg , iend
          if ( jsls <= ip ) then
            ften(i,j) = wg(jsls)*ften(i,j)+(d_one-wg(jsls))*bnd%wbt(i,jwb)
          end if
        end do
      else if ( jj <= ip ) then
!------east-boundary slice:
        do i = ibeg , iend
          if ( jsls <= ip ) then
            ften(i,j) = wg(jsls)*ften(i,j)+(d_one-wg(jsls))*bnd%ebt(i,jeb)
          end if
        end do
      end if
!
    end if
#endif

    call time_end(subroutine_name,idindx)
  end subroutine sponge2d
!
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
  function xfun(mm,ldot)
    implicit none
    real(8) :: xfun
    integer , intent(in) :: mm
    logical , intent(in) :: ldot
    if ( ldot ) then
      xfun = dble(nspgd-mm)/dble(nspgd-2)
    else
      xfun = dble(nspgx-mm)/dble(nspgx-2)
    end if
  end function xfun
!
  function xfune(mm,kk)
    implicit none
    real(8) :: xfune
    integer , intent(in) :: mm , kk
    xfune = dexp(-dble(mm-2)/anudg(kk))
  end function xfune
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
! These subroutines apply relaxation boundary conditions to the   c
! tendency term - ften - of variable f                            c
!                                                                 c
! ldot  : logical dot (u,v) / cross (t,q,p) flag                  c
!                                                                 c
! ip    : is the number of slices affected by nudging.            c
!                                                                 c
! xt    : is the time in seconds for variable f                   c
!                                                                 c
! fcoef : are the coefficients for the newtonian term.            c
!                                                                 c
! gcoef : are the coefficients for the diffusion term.            c
!                                                                 c
! ften  : is the tendency calculated from the model.              c
!                                                                 c
! j     : is the j'th slice of the tendency to be adjusted.       c
!                                                                 c
! nk    : is the number of vertical level to be adjusted.         c
!                                                                 c
! ibdy  : type of boundary condition relaxation, 1=linear         c
!         5 = exponential                                         c
!                                                                 c
! bnd   : Boundary condition data structure                       c
!         2D or 3D (managed by interface declaration)             c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine nudge3d(ldot,ip,fcoef,gcoef,xt,f,ften,j,nk,ibdy,bnd)
!
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
          fcx = fcoef*xfun(i,ldot)
          gcx = gcoef*xfun(i,ldot)
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
          fcx = fcoef*xfun(i,ldot)
          gcx = gcoef*xfun(i,ldot)
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
            fcx = fcoef*xfun(i,ldot)
            gcx = gcoef*xfun(i,ldot)
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
          fcx = fcoef*xfun(jsls,ldot)
          gcx = gcoef*xfun(jsls,ldot)
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
          fcx = fcoef*xfun(jsls,ldot)
          gcx = gcoef*xfun(jsls,ldot)
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
    k = nk ! not to have warnings
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
          fcx = fcoef*xfun(i,ldot)
          gcx = gcoef*xfun(i,ldot)
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
          fcx = fcoef*xfun(i,ldot)
          gcx = gcoef*xfun(i,ldot)
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
            fcx = fcoef*xfun(i,ldot)
            gcx = gcoef*xfun(i,ldot)
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
          fcx = fcoef*xfun(jsls,ldot)
          gcx = gcoef*xfun(jsls,ldot)
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
          fcx = fcoef*xfun(jsls,ldot)
          gcx = gcoef*xfun(jsls,ldot)
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
end module mod_bdycod
