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
!
  use mod_runparams
  use mod_memutil
  use mod_main
  use mod_che_main
  use mod_bats
  use mod_message 
  use mod_ncio
  use mod_cvaria
  use mod_mppio
  use mod_tcm_interface
!
  private
!
  public :: allocate_mod_bdycon , bdyin , bdyval
#ifndef BAND
  public :: uj1 , uj2 , ujlx , ujl
  public :: vj1 , vj2 , vjlx , vjl
#endif
  public :: ui1 , ui2 , uilx , uil
  public :: vi1 , vi2 , vilx , vil
  public :: ub0 , vb0 , qb0 , tb0 , ps0 , ts0 , so0
#ifndef BAND
  public :: peb , pebt , pwb , pwbt
  public :: ueb , uebt , veb , vebt , uwb , uwbt , vwb , vwbt
#endif
  public :: pnb , pnbt , psb , psbt
  public :: unb , unbt , vnb , vnbt , usb , usbt , vsb , vsbt
#ifndef BAND
  public :: teb , tebt , qeb , qebt , twb , twbt , qwb , qwbt
#endif
  public :: tnb , tnbt , qnb , qnbt , tsb , tsbt , qsb , qsbt
  public :: ts1 ! FOR DCSST
!
#ifndef BAND
  real(8) , pointer , dimension(:,:) :: uj1 , uj2 , ujl , ujlx ,  &
         vj1 , vj2 , vjl , vjlx
#endif
  real(8) , pointer , dimension(:,:) :: ps0 , ps1
  real(8) , pointer , dimension(:,:,:) :: qb0 , qb1 , so0 , so1 , &
         tb0 , tb1 , ub0 , ub1 , vb0 , vb1
  real(8) , pointer , dimension(:,:) :: ts0 , ts1
!
#ifndef BAND
  real(8) , pointer , dimension(:,:) :: peb , pebt , pwb , pwbt
#endif
  real(8) , pointer , dimension(:,:) :: pnb , pnbt , psbt , psb
#ifndef BAND
  real(8) , pointer , dimension(:,:,:) :: qeb , qebt , qwb , qwbt , &
        teb , tebt , twb , twbt , ueb , uebt , uwb , uwbt , veb ,   &
        vebt , vwb , vwbt
#endif
  real(8) , pointer , dimension(:,:,:) :: qnb , qnbt , qsb ,    &
        qsbt , tnb , tnbt , tsb , tsbt
  real(8) , pointer , dimension(:,:) :: ui1 , ui2 , uil , uilx ,&
        vi1 , vi2 , vil , vilx
  real(8) , pointer , dimension(:,:,:) :: unb , unbt , usb ,    &
        usbt , vnb , vnbt , vsb , vsbt
!
  contains
!
  subroutine allocate_mod_bdycon
    implicit none
    character (len=50) :: subroutine_name='allocate_mod_bdycon'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)

    call getmem2d(ps0,1,iy,0,jxp+1,'bdycon:ps0')
    call getmem2d(ps1,1,iy,0,jxp+1,'bdycon:ps1')
    call getmem2d(ts0,1,iy,1,jxp,'bdycon:ts0')
    call getmem2d(ts1,1,iy,1,jxp,'bdycon:ts1')
!
    call getmem3d(qb0,1,iy,1,kz,1,jxp,'bdycon:qb0')
    call getmem3d(qb1,1,iy,1,kz,1,jxp,'bdycon:qb1')
    call getmem3d(so0,1,iy,1,kz,1,jxp,'bdycon:so0')
    call getmem3d(so1,1,iy,1,kz,1,jxp,'bdycon:so1')
    call getmem3d(tb0,1,iy,1,kz,1,jxp,'bdycon:tb0')
    call getmem3d(tb1,1,iy,1,kz,1,jxp,'bdycon:tb1')
    call getmem3d(ub0,1,iy,1,kz,1,jxp,'bdycon:ub0')
    call getmem3d(ub1,1,iy,1,kz,1,jxp,'bdycon:ub1')
    call getmem3d(vb0,1,iy,1,kz,1,jxp,'bdycon:vb0')
    call getmem3d(vb1,1,iy,1,kz,1,jxp,'bdycon:vb1')
#ifndef BAND
    call getmem2d(peb,1,iy,0,jxp+1,'bdycon:peb')
    call getmem2d(pebt,1,iy,0,jxp+1,'bdycon:pebt')
    call getmem2d(pwb,1,iy,0,jxp+1,'bdycon:pwb')
    call getmem2d(pwbt,1,iy,0,jxp+1,'bdycon:pwbt')
#endif
    call getmem2d(pnb,1,nspgx,0,jxp+1,'bdycon:pnb')
    call getmem2d(pnbt,1,nspgx,0,jxp+1,'bdycon:pnbt')
    call getmem2d(psb,1,nspgx,0,jxp+1,'bdycon:psb')
    call getmem2d(psbt,1,nspgx,0,jxp+1,'bdycon:psbt')
#ifndef BAND
    call getmem3d(qeb,1,iy,1,kz,0,jxp+1,'bdycon:qeb')
    call getmem3d(qebt,1,iy,1,kz,0,jxp+1,'bdycon:qebt')
    call getmem3d(qwb,1,iy,1,kz,0,jxp+1,'bdycon:qwb')
    call getmem3d(qwbt,1,iy,1,kz,0,jxp+1,'bdycon:qwbt')
    call getmem3d(teb,1,iy,1,kz,0,jxp+1,'bdycon:teb')
    call getmem3d(tebt,1,iy,1,kz,0,jxp+1,'bdycon:tebt')
    call getmem3d(twb,1,iy,1,kz,0,jxp+1,'bdycon:twb')
    call getmem3d(twbt,1,iy,1,kz,0,jxp+1,'bdycon:twbt')
    call getmem3d(ueb,1,iy,1,kz,0,jxp+1,'bdycon:ueb')
    call getmem3d(uebt,1,iy,1,kz,0,jxp+1,'bdycon:uebt')
    call getmem3d(uwb,1,iy,1,kz,0,jxp+1,'bdycon:uwb')
    call getmem3d(uwbt,1,iy,1,kz,0,jxp+1,'bdycon:uwbt')
    call getmem3d(veb,1,iy,1,kz,0,jxp+1,'bdycon:veb')
    call getmem3d(vebt,1,iy,1,kz,0,jxp+1,'bdycon:vebt')
    call getmem3d(vwb,1,iy,1,kz,0,jxp+1,'bdycon:vwb')
    call getmem3d(vwbt,1,iy,1,kz,0,jxp+1,'bdycon:vwbt')
#endif
    call getmem3d(qnb,1,nspgx,1,kz,0,jxp+1,'bdycon:qnb')
    call getmem3d(qnbt,1,nspgx,1,kz,0,jxp+1,'bdycon:qnbt')
    call getmem3d(qsb,1,nspgx,1,kz,0,jxp+1,'bdycon:qsb')
    call getmem3d(qsbt,1,nspgx,1,kz,0,jxp+1,'bdycon:qsbt')
    call getmem3d(tnb,1,nspgx,1,kz,0,jxp+1,'bdycon:tnb')
    call getmem3d(tnbt,1,nspgx,1,kz,0,jxp+1,'bdycon:tnbt')
    call getmem3d(tsb,1,nspgx,1,kz,0,jxp+1,'bdycon:tsb')
    call getmem3d(tsbt,1,nspgx,1,kz,0,jxp+1,'bdycon:tsbt')
!
    call getmem2d(ui1,1,kz,0,jxp+1,'bdycon:ui1')
    call getmem2d(ui2,1,kz,0,jxp+1,'bdycon:ui2')
    call getmem2d(uil,1,kz,0,jxp+1,'bdycon:uil')
    call getmem2d(uilx,1,kz,0,jxp+1,'bdycon:uilx')
    call getmem2d(vi1,1,kz,0,jxp+1,'bdycon:vi1')
    call getmem2d(vi2,1,kz,0,jxp+1,'bdycon:vi2')
    call getmem2d(vil,1,kz,0,jxp+1,'bdycon:vil')
    call getmem2d(vilx,1,kz,0,jxp+1,'bdycon:vilx')
!
    call getmem3d(unb,1,nspgd,1,kz,0,jxp+1,'bdycon:unb')
    call getmem3d(unbt,1,nspgd,1,kz,0,jxp+1,'bdycon:unbt')
    call getmem3d(usb,1,nspgd,1,kz,0,jxp+1,'bdycon:usb')
    call getmem3d(usbt,1,nspgd,1,kz,0,jxp+1,'bdycon:usbt')
    call getmem3d(vnb,1,nspgd,1,kz,0,jxp+1,'bdycon:vnb')
    call getmem3d(vnbt,1,nspgd,1,kz,0,jxp+1,'bdycon:vnbt')
    call getmem3d(vsb,1,nspgd,1,kz,0,jxp+1,'bdycon:vsb')
    call getmem3d(vsbt,1,nspgd,1,kz,0,jxp+1,'bdycon:vsbt')
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
    real(8) , dimension(iy,jxp) :: psdot , tdum
    integer :: n
    character (len=50) :: subroutine_name='bdyin'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
    if ( dabs(xtime) > 0.0001D0 ) return
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
      write (6,'(a,i10)') 'SEARCH BC data for ', bdydate2%toidate()
      mmrec = icbc_search(bdydate2)
      if (mmrec < 0) then
        call open_icbc(monfirst(bdydate2))
        mmrec = icbc_search(bdydate2)
        if (mmrec < 0) then
          call fatal(__FILE__,__LINE__,'ICBC for '//bdydate2%tostring()//' not found')
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
    call bdydate2%broadcast(0,mpi_comm_world,ierr)
!
    call mpi_scatter(sav_0,iy*(kz*4+2)*jxp,mpi_real8,        &
                     sav0, iy*(kz*4+2)*jxp,mpi_real8,        &
                     0,mpi_comm_world,ierr)
    do j = 1 , jendl
      do k = 1 , kz
        do i = 1 , iy
          ub1(i,k,j) = sav0(i,k,j)
          vb1(i,k,j) = sav0(i,kz+k,j)
          qb1(i,k,j) = sav0(i,kz*2+k,j)
          tb1(i,k,j) = sav0(i,kz*3+k,j)
        end do
      end do
      do i = 1 , iy
        ps1(i,j) = sav0(i,kz*4+1,j)
        ts1(i,j) = sav0(i,kz*4+2,j)
      end do
    end do
    if ( ehso4 ) then
      call mpi_scatter(sav_0s,iy*kz*jxp,mpi_real8,  &
                       sav0s, iy*kz*jxp,mpi_real8,  &
                       0,mpi_comm_world,ierr)
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
        ps1(i,j) = ps1(i,j) - r8pt
      end do
    end do
!
!  this routine determines p(.) from p(x) by a 4-point
!  interpolation. on the x-grid, a p(x) point outside the grid
!  domain is assumed to satisfy p(0,j)=p(1,j); p(iy,j)=p(iym1,j);
!  and similarly for the i's.
!
    call mpi_sendrecv(ps1(1,jxp),iy,mpi_real8,ieast,1,   &
                      ps1(1,0),  iy,mpi_real8,iwest,1,   &
                      mpi_comm_world,mpi_status_ignore,ierr)
    do j = jbegin , jendx
      do i = 2 , iym1
        psdot(i,j) = d_rfour*(ps1(i,j)+ps1(i-1,j)+ps1(i,j-1)+ps1(i-1,j-1))
      end do
    end do
#ifdef BAND
    do j = jbegin , jendx
      psdot(1,j)  = d_half*(ps1(1,j)+ps1(1,j-1))
      psdot(iy,j) = d_half*(ps1(iym1,j)+ps1(iym1,j-1))
    end do
#else
!
    do i = 2 , iym1
      if ( myid == 0 )       psdot(i,1)     = d_half*(ps1(i,1)+ps1(i-1,1))
      if ( myid == nproc-1 ) psdot(i,jendl) = d_half*(ps1(i,jendx)+ps1(i-1,jendx))
    end do
!
    do j = jbegin , jendx
      psdot(1,j)  = d_half*(ps1(1,j)+ps1(1,j-1))
      psdot(iy,j) = d_half*(ps1(iym1,j)+ps1(iym1,j-1))
    end do
!
    if ( myid == 0 ) then
      psdot(1,1)  = ps1(1,1)
      psdot(iy,1) = ps1(iym1,1)
    end if
    if ( myid == nproc-1 ) then
      psdot(1,jendl)  = ps1(1,jendx)
      psdot(iy,jendl) = ps1(iym1,jendx)
    end if
#endif
!
!   Couple pressure u,v,t,q
!
    do k = 1 , kz
      do j = 1 , jendl
        do i = 1 , iy
          ub1(i,k,j) = ub1(i,k,j)*psdot(i,j)
          vb1(i,k,j) = vb1(i,k,j)*psdot(i,j)
          tb1(i,k,j) = tb1(i,k,j)*ps1(i,j)
          qb1(i,k,j) = qb1(i,k,j)*ps1(i,j)
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
        pwb(i,nn) = ps0(i,nn)
        pwbt(i,nn) = (ps1(i,nn)-ps0(i,nn))/dtbdys
      end do
    end do
    do nn = 1 , nxeb
      nnb = min0(jendx,jxp) - nn + 1
      do i = 1 , iym1
        peb(i,nn) = ps0(i,nnb)
        pebt(i,nn) = (ps1(i,nnb)-ps0(i,nnb))/dtbdys
      end do
    end do
#endif
    do nn = 1 , nspgx
      nnb = iym1 - nn + 1
      do j = 1 , jendx
        pnb(nn,j) = ps0(nnb,j)
        psb(nn,j) = ps0(nn,j)
        pnbt(nn,j) = (ps1(nnb,j)-ps0(nnb,j))/dtbdys
        psbt(nn,j) = (ps1(nn,j)-ps0(nn,j))/dtbdys
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
          uwb(i,k,nn)  = ub0(i,k,nn)
          vwb(i,k,nn)  = vb0(i,k,nn)
          uwbt(i,k,nn) = (ub1(i,k,nn)-ub0(i,k,nn))/dtbdys
          vwbt(i,k,nn) = (vb1(i,k,nn)-vb0(i,k,nn))/dtbdys
        end do
      end do
    end do
    do nn = 1 , ndeb
      nnb = min0(jendl,jxp) - nn + 1
      do k = 1 , kz
        do i = 1 , iy
          ueb(i,k,nn)  = ub0(i,k,nnb)
          veb(i,k,nn)  = vb0(i,k,nnb)
          uebt(i,k,nn) = (ub1(i,k,nnb)-ub0(i,k,nnb))/dtbdys
          vebt(i,k,nn) = (vb1(i,k,nnb)-vb0(i,k,nnb))/dtbdys
        end do
      end do
    end do
#endif
    do nn = 1 , nspgd
      nnb = iy - nn + 1
      do k = 1 , kz
        do j = 1 , jendl
          unb(nn,k,j)  = ub0(nnb,k,j)
          usb(nn,k,j)  = ub0(nn,k,j)
          vnb(nn,k,j)  = vb0(nnb,k,j)
          vsb(nn,k,j)  = vb0(nn,k,j)
          unbt(nn,k,j) = (ub1(nnb,k,j)-ub0(nnb,k,j))/dtbdys
          usbt(nn,k,j) = (ub1(nn,k,j)-ub0(nn,k,j))/dtbdys
          vnbt(nn,k,j) = (vb1(nnb,k,j)-vb0(nnb,k,j))/dtbdys
          vsbt(nn,k,j) = (vb1(nn,k,j)-vb0(nn,k,j))/dtbdys
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
          twb(i,k,nn)  = tb0(i,k,nn)
          qwb(i,k,nn)  = qb0(i,k,nn)
          twbt(i,k,nn) = (tb1(i,k,nn)-tb0(i,k,nn))/dtbdys
          qwbt(i,k,nn) = (qb1(i,k,nn)-qb0(i,k,nn))/dtbdys
        end do
      end do
    end do
    do nn = 1 , nxeb
      nnb = min0(jendx,jxp) - nn + 1
      do k = 1 , kz
        do i = 1 , iym1
          teb(i,k,nn)  = tb0(i,k,nnb)
          qeb(i,k,nn)  = qb0(i,k,nnb)
          tebt(i,k,nn) = (tb1(i,k,nnb)-tb0(i,k,nnb))/dtbdys
          qebt(i,k,nn) = (qb1(i,k,nnb)-qb0(i,k,nnb))/dtbdys
        end do
      end do
    end do
#endif
    do nn = 1 , nspgx
      nnb = iym1 - nn + 1
      do k = 1 , kz
        do j = 1 , jendx
          tnb(nn,k,j)  = tb0(nnb,k,j)
          tsb(nn,k,j)  = tb0(nn,k,j)
          qnb(nn,k,j)  = qb0(nnb,k,j)
          qsb(nn,k,j)  = qb0(nn,k,j)
          tnbt(nn,k,j) = (tb1(nnb,k,j)-tb0(nnb,k,j))/dtbdys
          tsbt(nn,k,j) = (tb1(nn,k,j)-tb0(nn,k,j))/dtbdys
          qnbt(nn,k,j) = (qb1(nnb,k,j)-qb0(nnb,k,j))/dtbdys
          qsbt(nn,k,j) = (qb1(nn,k,j)-qb0(nn,k,j))/dtbdys
        end do
      end do
    end do

    if ( myid == 0 ) then
      write (6,'(a,i10,a,i10)') 'READY  BC from     ' , &
            bdydate1%toidate() , ' to ' , bdydate2%toidate()
    end if

    bdydate1 = bdydate2
    call bdydate1%broadcast(0,mpi_comm_world,ierr)

    do j = 1 , jendx
      do i = 1 , iym1
        tdum(i,j) = ts1(i,j)
      end do
    end do
    do k = 1 , kz
      do j = 1 , jendl
        do i = 1 , iy
          ub0(i,k,j) = ub1(i,k,j)
          vb0(i,k,j) = vb1(i,k,j)
          qb0(i,k,j) = qb1(i,k,j)
          tb0(i,k,j) = tb1(i,k,j)
        end do
      end do
    end do
    do j = 1 , jendl
      do i = 1 , iy
        ps0(i,j) = ps1(i,j)
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
    character (len=50) :: subroutine_name='bdyuv'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
!   compute the p* at dot points:
!
    call mpi_sendrecv(sps1%ps(1,jxp),iy,mpi_real8,ieast,1,   &
                      sps1%ps(1,0),  iy,mpi_real8,iwest,1,   &
                      mpi_comm_world,mpi_status_ignore,ierr)
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
            uj1(i,k) = uwb(i,k,1)/sps1%pdot(i,1)
            vj1(i,k) = vwb(i,k,1)/sps1%pdot(i,1)
          end if
          if ( myid == nproc-1 ) then
            ujl(i,k) = ueb(i,k,1)/sps1%pdot(i,jendl)
            vjl(i,k) = veb(i,k,1)/sps1%pdot(i,jendl)
          end if
        end do
#endif
!
!       south (i = 1) and north (i = iy) boundaries:
!
        do j = 1 , jendl
          ui1(k,j) = usb(1,k,j)/sps1%pdot(1,j)
          vi1(k,j) = vsb(1,k,j)/sps1%pdot(1,j)
          uil(k,j) = unb(1,k,j)/sps1%pdot(iy,j)
          vil(k,j) = vnb(1,k,j)/sps1%pdot(iy,j)
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
            uj1(i,k) = (uwb(i,k,1)+dtb*uwbt(i,k,1))/sps1%pdot(i,1)
            vj1(i,k) = (vwb(i,k,1)+dtb*vwbt(i,k,1))/sps1%pdot(i,1)
          end if
          if ( myid == nproc-1 ) then
            ujl(i,k) = (ueb(i,k,1)+dtb*uebt(i,k,1))/sps1%pdot(i,jendl)
            vjl(i,k) = (veb(i,k,1)+dtb*vebt(i,k,1))/sps1%pdot(i,jendl)
          end if
        end do
#endif
!
!       south (i = 1) and north (i = iy) boundaries:
!
        do j = 1 , jendl
          ui1(k,j) = (usb(1,k,j)+dtb*usbt(1,k,j))/sps1%pdot(1,j)
          vi1(k,j) = (vsb(1,k,j)+dtb*vsbt(1,k,j))/sps1%pdot(1,j)
          uil(k,j) = (unb(1,k,j)+dtb*unbt(1,k,j))/sps1%pdot(iy,j)
          vil(k,j) = (vnb(1,k,j)+dtb*vnbt(1,k,j))/sps1%pdot(iy,j)
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
                      mpi_comm_world,mpi_status_ignore,ierr)
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
                      mpi_comm_world,mpi_status_ignore,ierr)
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
    character (len=50) :: subroutine_name='bdyval'
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
    if ( dabs(xt) < 0.00001D0 .and. ktau > 0 ) then
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
          if ( myid == 0 )       sps1%ps(i,1)     = pwb(i,1)
          if ( myid == nproc-1 ) sps1%ps(i,jendx) = peb(i,1)
        end do
#endif
        do j = jbegin , jendm
          sps1%ps(1,j)    = psb(1,j)
          sps1%ps(iym1,j) = pnb(1,j)
        end do
!
        do k = 1 , kz
#ifndef BAND
          do i = 1 , iy
            if ( myid == 0 ) then
              atm1%u(i,k,1) = uwb(i,k,1)
              atm1%v(i,k,1) = vwb(i,k,1)
            end if
            if ( myid == nproc-1 ) then
              atm1%u(i,k,jendl) = ueb(i,k,1)
              atm1%v(i,k,jendl) = veb(i,k,1)
            end if
          end do
#endif
          do j = jbegin , jendx
            atm1%u(1,k,j)  = usb(1,k,j)
            atm1%u(iy,k,j) = unb(1,k,j)
            atm1%v(1,k,j)  = vsb(1,k,j)
            atm1%v(iy,k,j) = vnb(1,k,j)
          end do
        end do
      end if
!
!     time-dependent boundary conditions:
!
#ifndef BAND
      do i = 1 , iym1
        if ( myid == 0 )       sps1%ps(i,1)     = pwb(i,1) + dtb*pwbt(i,1)
        if ( myid == nproc-1 ) sps1%ps(i,jendx) = peb(i,1) + dtb*pebt(i,1)
      end do
#endif
      do j = jbegin , jendm
        sps1%ps(1,j)    = psb(1,j) + dtb*psbt(1,j)
        sps1%ps(iym1,j) = pnb(1,j) + dtb*pnbt(1,j)
      end do
!
      do k = 1 , kz
#ifndef BAND
        do i = 1 , iy
          if ( myid == 0 ) then
            atm1%u(i,k,1) = uwb(i,k,1) + dtb*uwbt(i,k,1)
            atm1%v(i,k,1) = vwb(i,k,1) + dtb*vwbt(i,k,1)
          end if
          if ( myid == nproc-1 ) then
            atm1%u(i,k,jendl) = ueb(i,k,1) + dtb*uebt(i,k,1)
            atm1%v(i,k,jendl) = veb(i,k,1) + dtb*vebt(i,k,1)
          end if
        end do
#endif
        do j = jbegin , jendx
          atm1%u(1,k,j)  = usb(1,k,j) + dtb*usbt(1,k,j)
          atm1%u(iy,k,j) = unb(1,k,j) + dtb*unbt(1,k,j)
          atm1%v(1,k,j)  = vsb(1,k,j) + dtb*vsbt(1,k,j)
          atm1%v(iy,k,j) = vnb(1,k,j) + dtb*vnbt(1,k,j)
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
          if ( myid == 0 )       atm1%t(i,k,1)     = twb(i,k,1)
          if ( myid == nproc-1 ) atm1%t(i,k,jendx) = teb(i,k,1)
        end do
#endif
        do j = jbegin , jendm
          atm1%t(1,k,j)    = tsb(1,k,j)
          atm1%t(iym1,k,j) = tnb(1,k,j)
        end do
      end do
      do k = 1 , kz
#ifndef BAND
        do i = 1 , iym1
          if ( myid == 0 )       atm1%qv(i,k,1)     = qwb(i,k,1)
          if ( myid == nproc-1 ) atm1%qv(i,k,jendx) = qeb(i,k,1)
        end do
#endif
        do j = jbegin , jendm
          atm1%qv(1,k,j)    = qsb(1,k,j)
          atm1%qv(iym1,k,j) = qnb(1,k,j)
        end do
      end do
    end if
!
!   time-dependent boundary conditions:
!
    do k = 1 , kz
#ifndef BAND
      do i = 1 , iym1
        if ( myid == 0 )       atm1%t(i,k,1)     = twb(i,k,1) + dtb*twbt(i,k,1)
        if ( myid == nproc-1 ) atm1%t(i,k,jendx) = teb(i,k,1) + dtb*tebt(i,k,1)
      end do
#endif
      do j = jbegin , jendm
        atm1%t(1,k,j)    = tsb(1,k,j) + dtb*tsbt(1,k,j)
        atm1%t(iym1,k,j) = tnb(1,k,j) + dtb*tnbt(1,k,j)
      end do
    end do
    do k = 1 , kz
#ifndef BAND
      do i = 1 , iym1
        if ( myid == 0 )       atm1%qv(i,k,1)     = qwb(i,k,1) + dtb*qwbt(i,k,1)
        if ( myid == nproc-1 ) atm1%qv(i,k,jendx) = qeb(i,k,1) + dtb*qebt(i,k,1)
      end do
#endif
      do j = jbegin , jendm
        atm1%qv(1,k,j)    = qsb(1,k,j) + dtb*qsbt(1,k,j)
        atm1%qv(iym1,k,j) = qnb(1,k,j) + dtb*qnbt(1,k,j)
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
      call set_tke_bc
    end if
!
    call time_end(subroutine_name,idindx)
!
  end subroutine bdyval
!
end module mod_bdycod
