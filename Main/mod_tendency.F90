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
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
module mod_tendency

  use mod_runparams
  use mod_mppparam
  use mod_mpmessage
  use mod_service
  use mod_memutil
  use mod_atm_interface
  use mod_che_interface
  use mod_cu_interface
  use mod_lm_interface
  use mod_rad_interface
  use mod_pbl_interface
  use mod_bdycod
  use mod_che_bdyco
  use mod_precip
  use mod_slice
  use mod_sun
  use mod_diagnosis
  use mod_advection
  use mod_diffusion
  use mod_mppio
#ifdef CLM
  use mod_clm
  use mod_mtrxclm
  use clm_varsur
#endif

  private

  public :: allocate_mod_tend , tend

  real(dp) , pointer , dimension(:,:,:) :: divl
  real(dp) , pointer , dimension(:,:,:) :: ttld , xkc , xkcf , td , phi
  real(dp) , pointer , dimension(:,:,:) :: ps4
  real(dp) , pointer , dimension(:,:,:) :: ps_4 
  real(dp) , pointer , dimension(:,:) :: psc , psd , pten
  real(dp) , pointer , dimension(:,:,:) :: wrkkuo1
  real(dp) , pointer , dimension(:,:,:) :: wrkkuo2

  integer(8) , parameter :: irep = 50

  contains

  subroutine allocate_mod_tend
    implicit none

    call getmem3d(divl,1,jxp,icross1,icross2,1,kz,'tendency:divl')
    call getmem3d(ttld,1,jxp,icross1,icross2,1,kz,'tend:ttld')
    call getmem3d(phi,0,jxp,idot1,idot2,1,kz,'tendency:phi')
    call getmem3d(xkc,1,jxp,idot1,idot2,1,kz,'tendency:xkc')
    call getmem3d(xkcf,1,jxp,idot1,idot2,1,kzp1,'tendency:xkcf')
    call getmem3d(ps4,1,jxp,icross1,icross2,1,4,'tendency:ps4')
    call getmem3d(ps_4,1,jx,icross1,icross2,1,4,'tendency:ps_4')
    call getmem3d(td,1,jxp,icross1,icross2,1,kz,'tendency:td')
    call getmem2d(psc,1,jxp,icross1,icross2,'tendency:psc')
    call getmem2d(psd,0,jxp+1,icross1,icross2,'tendency:psd')
    call getmem2d(pten,1,jxp,icross1,icross2,'tendency:pten')
    if ( icup == 1 ) then
      call getmem3d(wrkkuo1,0,jxp+1,icross1,icross2,1,kz,'tendency:wrkkuo1')
      call getmem3d(wrkkuo2,0,jxp+1,icross1,icross2,1,kz,'tendency:wrkkuo2')
    end if
  end subroutine allocate_mod_tend

  subroutine tend

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the tendencies of the prognostic       c
!     variables p*, u, v, and t.                                      c
!                                                                     c
!     p*u, p*v, p*t ,p*qv, and p*qc stored in main common block.      c
!                                                                     c
!     all the two-dimension arrays stored in main common block.       c
!                                                                     c
!     east/west boundary conditions stored in common block /bdycod/ . c
!                                                                     c
!     north/south boundary conditions stored in common block          c
!              /bdycod/.                                              c
!                                                                     c
!     all the integers stored in common block /param1/.               c
!                                                                     c
!     all the constants stored in common block /param1/.              c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
#ifndef IBM
    use mpi
#else
    include 'mpif.h'
#endif
    implicit none
!
    real(dp) :: cell , chias , chibs , dudx , dudy , dvdx , dvdy , &
               psasum , pt2bar , pt2tot , ptnbar , maxv ,          &
               ptntot , qcas , qcbs , qvas , qvbs , rovcpm ,       &
               rtbar , sigpsa , tv , tv1 , tv2 , tv3 , tv4 , tva , &
               tvavg , tvb , tvc , xmsf , xtm1 , theta , eccf,sod
    real(dp) , pointer , dimension(:,:,:) :: spchiten , spchi , spchia , &
                                            spchib3d
    integer :: i , iptn , itr , j , k , lev , n , ii , jj , kk
    integer :: ierr , icons_mpi
    logical :: loutrad , labsem
    character (len=32) :: appdat
    character (len=64) :: subroutine_name='tend'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
    !
    ! Calculate eccentricity factor for radiation calculations
    !
#ifdef CLM
    eccf  = r2ceccf
#else
    calday = yeardayfrac(idatex)
    theta = twopi*calday/dayspy
    eccf = 1.000110D0 + 0.034221D0*dcos(theta) +  &
           0.001280D0 * dsin(theta) + &
           0.000719D0 * dcos(d_two*theta) + &
           0.000077D0 * dsin(d_two*theta)
#endif
!
!   multiply ua and va by inverse of mapscale factor at dot point:
!
    do k = 1 , kz
      do i = ide1 , ide2
        do j = jde1 , jde2
          atm1%u(j,i,k) = atm1%u(j,i,k)*mddom%msfd(j,i)
          atm1%v(j,i,k) = atm1%v(j,i,k)*mddom%msfd(j,i)
        end do
      end do
    end do

    call deco1_exchange_left(sfs%psa,1,ice1,ice2)
    call deco1_exchange_right(sfs%psa,1,ice1,ice2)
    call psc2psd(sfs%psb,psdot)
    !
    ! Internal U,V points
    !
    do k = 1 , kz
      do i = idii1 , idii2
        do j = jdii1 , jdii2
          xmsf = mddom%msfd(j,i)
          atmx%u(j,i,k) = atm1%u(j,i,k)/(psdot(j,i)*xmsf)
          atmx%v(j,i,k) = atm1%v(j,i,k)/(psdot(j,i)*xmsf)
        end do
      end do
    end do
    !
    ! Boundary points
    !
    if ( ma%hasleft ) then
      do k = 1 , kz
        do i = idi1 , idi2
          atmx%u(jdi1,i,k) = wui(i,k)
          atmx%v(jdi1,i,k) = wvi(i,k)
        end do
      end do
      do k = 1 , kz
        do i = idi1 , idi2
          atmx%u(jde1,i,k) = wue(i,k)
          atmx%v(jde1,i,k) = wve(i,k)
        end do
      end do
      ! inflow/outflow dependence
      if ( iboudy == 3 .or. iboudy == 4 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            if ( atmx%u(jde1,i,k) <= d_zero ) then
              atmx%u(jde1,i,k) = atmx%u(jdi1,i,k)
              atmx%v(jde1,i,k) = atmx%v(jdi1,i,k)
            end if
          end do
        end do
      end if
    end if
    if ( ma%hasright ) then
      do k = 1 , kz
        do i = idi1 , idi2
          atmx%u(jdi2,i,k) = eui(i,k)
          atmx%v(jdi2,i,k) = evi(i,k)
        end do
      end do
      do k = 1 , kz
        do i = idi1 , idi2
          atmx%u(jde2,i,k) = eue(i,k)
          atmx%v(jde2,i,k) = eve(i,k)
        end do
      end do
      ! inflow/outflow dependence
      if ( iboudy == 3 .or. iboudy == 4 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            if ( atmx%u(jde2,i,k) >= d_zero ) then
              atmx%u(jde2,i,k) = atmx%u(jdi2,i,k)
              atmx%v(jde2,i,k) = atmx%v(jdi2,i,k)
            end if
          end do
        end do
      end if
    end if
    if ( ma%hasbottom ) then
      do k = 1 , kz
        do j = jdi1 , jdi2
          atmx%u(j,idi1,k) = sui(j,k)
          atmx%v(j,idi1,k) = svi(j,k)
        end do
      end do
      do k = 1 , kz
        do j = jde1 , jde2
          atmx%u(j,ide1,k) = sue(j,k)
          atmx%v(j,ide1,k) = sve(j,k)
        end do
      end do
      if ( iboudy == 3 .or. iboudy == 4 ) then
        ! inflow/outflow dependence
        do k = 1 , kz
          do j = jde1 , jde2
            if ( atmx%v(j,ide1,k) >= d_zero ) then
              atmx%v(j,ide1,k) = atmx%v(j,idi1,k)
              atmx%u(j,ide1,k) = atmx%u(j,idi1,k)
            end if
          end do
        end do
      end if
    end if
    if ( ma%hastop ) then
      do k = 1 , kz
        do j = jdi1 , jdi2
          atmx%u(j,idi2,k) = nui(j,k)
          atmx%v(j,idi2,k) = nvi(j,k)
        end do
      end do
      do k = 1 , kz
        do j = jde1 , jde2
          atmx%u(j,ide2,k) = nue(j,k)
          atmx%v(j,ide2,k) = nve(j,k)
        end do
      end do
      if ( iboudy == 3 .or. iboudy == 4 ) then
        ! inflow/outflow dependence
        do k = 1 , kz
          do j = jde1 , jde2
            if ( atmx%v(j,ide2,k) <= d_zero ) then
              atmx%v(j,ide2,k) = atmx%v(j,idi2,k)
              atmx%u(j,ide2,k) = atmx%u(j,idi2,k)
            end if
          end do
        end do
      end if
    end if
   !
   ! T , QV , QC decouple
   !
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          atmx%t(j,i,k)  = atm1%t(j,i,k)/sfs%psa(j,i)
          atmx%qv(j,i,k) = atm1%qv(j,i,k)/sfs%psa(j,i)
          atmx%qc(j,i,k) = atm1%qc(j,i,k)/sfs%psa(j,i)
        end do
      end do
    end do
    !
    ! call tracer decoupling routine for multiple (ntr) species
    !
    if ( ichem == 1 ) then
      do n = 1 , ntr
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              chi(j,i,k,n) = chia(j,i,k,n)/sfs%psa(j,i)
            end do
          end do
        end do
      end do
    end if
!
!=======================================================================
!
    call deco1_exchange_left(sfs%psb,1,ice1,ice2)
    call deco1_exchange_right(sfs%psb,1,ice1,ice2)
!
    call deco1_exchange_left(atm1%u,1,ide1,ide2,1,kz)
    call deco1_exchange_right(atm1%u,1,ide1,ide2,1,kz)
    call deco1_exchange_left(atm1%v,1,ide1,ide2,1,kz)
    call deco1_exchange_right(atm1%v,1,ide1,ide2,1,kz)
!
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      call deco1_exchange_left(atm1%tke,1,ice1,ice2,1,kz)
      call deco1_exchange_right(atm1%tke,1,ice1,ice2,1,kz)
    end if
!
    call deco1_exchange_left(atm2%u,1,ide1,ide2,1,kz)
    call deco1_exchange_right(atm2%u,1,ide1,ide2,1,kz)
    call deco1_exchange_left(atm2%v,1,ide1,ide2,1,kz)
    call deco1_exchange_right(atm2%v,1,ide1,ide2,1,kz)
    call deco1_exchange_left(atm2%t,1,ice1,ice2,1,kz)
    call deco1_exchange_right(atm2%t,1,ice1,ice2,1,kz)
    call deco1_exchange_left(atm2%qv,1,ice1,ice2,1,kz)
    call deco1_exchange_right(atm2%qv,1,ice1,ice2,1,kz)
!
    call deco1_exchange_left(atmx%u,1,ide1,ide2,1,kz)
    call deco1_exchange_right(atmx%u,1,ide1,ide2,1,kz)
    call deco1_exchange_left(atmx%v,1,ide1,ide2,1,kz)
    call deco1_exchange_right(atmx%v,1,ide1,ide2,1,kz)
    call deco1_exchange_left(atmx%t,1,ice1,ice2,1,kz)
    call deco1_exchange_right(atmx%t,1,ice1,ice2,1,kz)
    call deco1_exchange_left(atmx%qv,1,ice1,ice2,1,kz)
    call deco1_exchange_right(atmx%qv,1,ice1,ice2,1,kz)
    call deco1_exchange_left(atmx%qc,1,ice1,ice2,1,kz)
    call deco1_exchange_right(atmx%qc,1,ice1,ice2,1,kz)
!
    if ( ichem == 1 ) then
      call deco1_exchange_left(chi,1,ice1,ice2,1,kz,1,ntr)
      call deco1_exchange_right(chi,1,ice1,ice2,1,kz,1,ntr)
      call deco1_exchange_left(chib,1,ice1,ice2,1,kz,1,ntr)
      call deco1_exchange_right(chib,1,ice1,ice2,1,kz,1,ntr)
    end if
!
!   Calculate pdot
!
    call deco1_exchange_left(sfs%psb,1,ice1,ice2)
    call deco1_exchange_right(sfs%psb,1,ice1,ice2)
    call psc2psd(sfs%psb,psdot)
!
!=======================================================================
!
    call mkslice

#ifdef CLM
    if ( init_grid ) then
      call initclm(ifrest,idate1,idate2,dx,dtrad,dtsrf)
      init_grid = .false.
    end if
#endif
!
!=======================================================================
!
    call deco1_exchange_left(atms%ubd3d,2,ide1,ide2,1,kz)
    call deco1_exchange_right(atms%ubd3d,2,ide1,ide2,1,kz)
    call deco1_exchange_left(atms%vbd3d,2,ide1,ide2,1,kz)
    call deco1_exchange_right(atms%vbd3d,2,ide1,ide2,1,kz)
    call deco1_exchange_left(atms%tb3d,2,ice1,ice2,1,kz)
    call deco1_exchange_right(atms%tb3d,2,ice1,ice2,1,kz)
    call deco1_exchange_left(atms%ubx3d,2,ice1,ice2,1,kz)
    call deco1_exchange_right(atms%ubx3d,2,ice1,ice2,1,kz)
    call deco1_exchange_left(atms%vbx3d,2,ice1,ice2,1,kz)
    call deco1_exchange_right(atms%vbx3d,2,ice1,ice2,1,kz)
    call deco1_exchange_left(atms%qvb3d,2,ice1,ice2,1,kz)
    call deco1_exchange_right(atms%qvb3d,2,ice1,ice2,1,kz)
    call deco1_exchange_left(atms%qcb3d,2,ice1,ice2,1,kz)
    call deco1_exchange_right(atms%qcb3d,2,ice1,ice2,1,kz)

    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      call deco1_exchange_left(atm2%tke,2,ice1,ice2,1,kz)
      call deco1_exchange_right(atm2%tke,2,ice1,ice2,1,kz)
    end if

    if ( ichem == 1 ) then
      call deco1_exchange_left(atms%chib3d,2,ice1,ice2,1,kz,1,ntr)
      call deco1_exchange_right(atms%chib3d,2,ice1,ice2,1,kz,1,ntr)
    end if
!
!**********************************************************************
!
!     compute the pressure tendency:
!
    pten(:,:)   = d_zero
    qdot(:,:,:) = d_zero
    omega(:,:,:) = d_zero

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          divl(j,i,k) = (atm1%u(j+1,i+1,k)+atm1%u(j+1,i,k)- &
                         atm1%u(j,i+1,k)  -atm1%u(j,i,k)) + &
                        (atm1%v(j+1,i+1,k)+atm1%v(j,i+1,k)- &
                         atm1%v(j+1,i,k)  -atm1%v(j,i,k))
          pten(j,i) = pten(j,i) - divl(j,i,k)*dsigma(k) / &
                      (dx2*mddom%msfx(j,i)*mddom%msfx(j,i))
        end do
      end do
    end do
!
!   compute vertical sigma-velocity (qdot):
!
    do k = 2 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          qdot(j,i,k) = qdot(j,i,k-1) - (pten(j,i)+divl(j,i,k-1) / &
                 (dx2*mddom%msfx(j,i)*mddom%msfx(j,i)))*dsigma(k-1)/sfs%psa(j,i)
         end do
      end do
    end do
    call deco1_exchange_left(qdot,1,ide1,ide2,1,kz)
    call deco1_exchange_right(qdot,1,ide1,ide2,1,kz)
!
!   compute omega
!
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          omega(j,i,k) = d_half*sfs%psa(j,i)* &
                     (qdot(j,i,k+1)+qdot(j,i,k))+a(k)*(pten(j,i)+   &
                     ((atmx%u(j,i,k)    +atmx%u(j,i+1,k)+           &
                       atmx%u(j+1,i+1,k)+atmx%u(j+1,i,k))*          &
                     (sfs%psa(j+1,i)-sfs%psa(j-1,i))+               &
                     (atmx%v(j,i,k)    +atmx%v(j,i+1,k)+            &
                      atmx%v(j+1,i+1,k)+atmx%v(j+1,i,k))*           &
                     (sfs%psa(j,i+1)-sfs%psa(j,i-1)))/              &
                     (dx8*mddom%msfx(j,i)))
        end do
      end do
    end do

    if ( ichem == 1 ) then
#ifndef BAND
      call deco1_exchange_left(chieb,1,ice1,ice2,1,kz,1,ntr)
      call deco1_exchange_right(chieb,1,ice1,ice2,1,kz,1,ntr)
      call deco1_exchange_left(chiebt,1,ice1,ice2,1,kz,1,ntr)
      call deco1_exchange_right(chiebt,1,ice1,ice2,1,kz,1,ntr)
      call deco1_exchange_left(chiwb,1,ice1,ice2,1,kz,1,ntr)
      call deco1_exchange_right(chiwb,1,ice1,ice2,1,kz,1,ntr)
      call deco1_exchange_left(chiwbt,1,ice1,ice2,1,kz,1,ntr)
      call deco1_exchange_right(chiwbt,1,ice1,ice2,1,kz,1,ntr)
#endif
      call deco1_exchange_left(chinb,1,1,nspgx,1,kz,1,ntr)
      call deco1_exchange_right(chinb,1,1,nspgx,1,kz,1,ntr)
      call deco1_exchange_left(chinbt,1,1,nspgx,1,kz,1,ntr)
      call deco1_exchange_right(chinbt,1,1,nspgx,1,kz,1,ntr)
      call deco1_exchange_left(chisb,1,1,nspgx,1,kz,1,ntr)
      call deco1_exchange_right(chisb,1,1,nspgx,1,kz,1,ntr)
      call deco1_exchange_left(chisbt,1,1,nspgx,1,kz,1,ntr)
      call deco1_exchange_right(chisbt,1,1,nspgx,1,kz,1,ntr)
    end if
!
    if ( iboudy == 4 ) then
      call sponge(ba_cr,xpsb,pten)
    else if ( iboudy == 1 .or. iboudy == 5 ) then
      xtm1 = xbctime - dtsec
      if ( nbdytime == 0 .and. ktau /= 0 ) xtm1 = -dtsec
      call nudge(ba_cr,xtm1,sfs%psb,iboudy,xpsb,pten)
    end if
    !
    ! psc : forecast pressure
    ! psd : weighted p* (psd)
    !
    do i = ice1 , ice2
      do j = jce1 , jce2
        psd(j,i) = sfs%psa(j,i)
      end do
    end do
    call deco1_exchange_left(psd,1,ice1,ice2)

    do i = ici1 , ici2
      do j = jci1 , jci2
        psc(j,i) = sfs%psb(j,i) + pten(j,i)*dt
      end do
    end do
    if ( ma%hasleft ) then
      do i = ici1 , ici2
        psc(jce1,i) = sfs%psb(jce1,i) + xpsb%bt(jce1,i)*dt
      end do
    end if
    if ( ma%hasright ) then
      do i = ici1 , ici2
        psc(jce2,i) = sfs%psb(jce2,i) + xpsb%bt(jce2,i)*dt
      end do
    end if
    if ( ma%hasbottom ) then
      do j = jce1 , jce2
        psc(j,ice1) = sfs%psb(j,ice1) + xpsb%bt(j,ice1)*dt
      end do
    end if
    if ( ma%hastop ) then
      do j = jce1 , jce2
        psc(j,ice2) = sfs%psb(j,ice2) + xpsb%bt(j,ice2)*dt
      end do
    end if
!
!   compute bleck (1977) noise parameters:
!
    do j = jci1 , jci2
      do i = ici1 , ici2
        ps4(j,i,1) = pten(j,i)
        ps4(j,i,2) = psc(j,i)
        ps4(j,i,3) = sfs%psb(j,i)
        ps4(j,i,4) = sfs%psa(j,i)
      end do
    end do

    call deco1_gather(ps4,ps_4,jci1,jci2,ici1,ici2,1,4)

    if ( ktau /= 0 ) then
      iptn = 0
      ptntot = d_zero
      pt2tot = d_zero
      if ( myid == 0 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            iptn = iptn + 1
            ptntot = ptntot + dabs(ps_4(j,i,1))
            pt2tot = pt2tot +                       &
                     dabs((ps_4(j,i,2)+ps_4(j,i,3)- &
                           d_two*ps_4(j,i,4))/(dt*dt*d_rfour))
          end do
        end do
      end if
      call mpi_bcast(iptn,1,mpi_integer,0,mycomm,ierr)
      call mpi_bcast(ptntot,1,mpi_real8,0,mycomm,ierr)
      call mpi_bcast(pt2tot,1,mpi_real8,0,mycomm,ierr)
    end if
!
!   calculate solar zenith angle
!
    if ( ktau == 0 .or. ichem == 1 .or. &
      mod(ktau+1,ntsrf) == 0 .or. mod(ktau+1,ntrad) == 0 ) then
      call zenitm(coszrs,jci1,jci2,ici1,ici2)
    end if

    xkc(:,:,:) = d_zero 
    !
    ! No diffusion of TKE on lower boundary (kzp1)
    !
    xkcf(:,:,:) = d_zero 
!
!     compute the horizontal diffusion coefficient and stored in xkc:
!     the values are calculated at cross points, but they also used
!     for dot-point variables.
!
    do k = 2 , kz
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          dudx = atm2%u(j+1,i,k) + atm2%u(j+1,i+1,k) - &
                 atm2%u(j,i,k)   - atm2%u(j,i+1,k)
          dvdx = atm2%v(j+1,i,k) + atm2%v(j+1,i+1,k) - &
                 atm2%v(j,i,k)   - atm2%v(j,i+1,k)
          dudy = atm2%u(j,i+1,k) + atm2%u(j+1,i+1,k) - &
                 atm2%u(j,i,k)   - atm2%u(j+1,i,k)
          dvdy = atm2%v(j,i+1,k) + atm2%v(j+1,i+1,k) - &
                 atm2%v(j,i,k)   - atm2%v(j+1,i,k)
          cell = (xkhz*hgfact(j,i)+c200 * &
                  dsqrt((dudx-dvdy)*(dudx-dvdy)+(dvdx+dudy)*(dvdx+dudy)))
          xkc(j,i,k) = dmin1(cell,xkhmax)
          if ( k > 1 ) then
            ! TAO: Interpolate the diffusion coefficients to full levels
            ! for use in diffusion of TKE
            cell = twt(k,1)*xkc(j,i,k) + twt(k,2)*xkc(j,i,k-1)
            xkcf(j,i,k) = dmin1(nuk*cell,xkhmax)
          else
            ! TAO: Multiply the horizontal diffusion coefficient by
            ! nuk for TKE.  Without this multiplication, it appears
            ! that TKE does not diffuse fast enough, and stabilities
            ! appear in the TKE field.  While this is different from
            ! Bretherton's treatment, it is consistent with the
            ! scaling of the vertical TKE diffusivity.
            xkcf(j,i,1) = nuk*xkc(j,i,1)
          end if
        end do
      end do
    end do
!
!   compute the temperature tendency:
!
    aten%t(:,:,:) = d_zero
    aten%qv(:,:,:) = d_zero
    aten%qc(:,:,:) = d_zero
!
!   compute the horizontal advection term:
!
    call hadv(cross,aten%t,atmx%t,kz,1)
!
!   compute the vertical advection term:
!
    if ( ibltyp /= 2 .and. ibltyp /= 99 ) then
      call vadv(cross,aten%t,atm1%t,kz,1)
    else
      if ( iuwvadv == 1 ) then
        call vadv(cross,aten%t,atm1%t,kz,6)
      else
        call vadv(cross,aten%t,atm1%t,kz,1)
      end if
    end if
!
!   compute the adiabatic term:
!
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          rovcpm = rgas/(cpd*(d_one+0.8D0*(atmx%qv(j,i,k))))
          tv = atmx%t(j,i,k)*(d_one+ep1*(atmx%qv(j,i,k)))
          aten%t(j,i,k) = aten%t(j,i,k) + (omega(j,i,k)*rovcpm*tv) / &
                          (ptop/sfs%psa(j,i)+a(k))
        end do
      end do
    end do
!
!   compute the diffusion term for t and store in difft:
!
    adf%difft(:,:,:) = d_zero
    adf%diffq(:,:,:) = d_zero
!
    call diffu_x(adf%difft,atms%tb3d,sfs%psb,xkc,kz)
!
!   compute the moisture tendencies:
!
!   icup = 1 : kuo-anthes cumulus parameterizaion scheme
!   icup = 2 : grell cumulus paramterization scheme
!   icup = 3 : betts-miller (1986)
!   icup = 4 : emanuel (1991)
!   icup = 5 : tiedtke (1986)
!   icup = 99: grell over land, emanuel over ocean
!   icup = 98: emanuel over land, grell over ocean
!
    call hadv(cross,aten%qv,atmx%qv,kz,1)
    if ( icup /= 1 ) then
      if ( ibltyp /= 2 .and. ibltyp /= 99 ) then
        call vadv(cross,aten%qv,atm1%qv,kz,2)
      else
        if ( iuwvadv == 1 ) then
          call vadv(cross,aten%qv,atm1%qv,kz,6)
        else
          call vadv(cross,aten%qv,atm1%qv,kz,2)
        end if
      end if
    end if
    !
    ! Zero out radiative clouds
    !
    cldfra(:,:,:) = d_zero
    cldlwc(:,:,:) = d_zero

    if ( icup == 1 ) then
      call cupara(jci1,jci2,ici1,ici2,ktau)
    end if
    if ( icup == 2 .or. icup == 99 .or. icup == 98 ) then
      call cuparan(jci1,jci2,ici1,ici2,ktau)
    end if
    if ( icup == 3 ) then
      call bmpara(jci1,jci2,ici1,ici2,ktau)
    end if
    if ( icup == 4 .or. icup == 99 .or. icup == 98 ) then
      call cupemandrv(jci1,jci2,ici1,ici2,ktau)
    end if
    if ( icup == 5 ) then
      call tiedtkedrv(jci1,jci2,ici1,ici2,ktau)
    end if

    if ( ipptls == 1 ) then
      call hadv(cross,aten%qc,atmx%qc,kz,2)
      if ( ibltyp /= 2 .and. ibltyp /= 99 ) then
        call vadv(cross,aten%qc,atm1%qc,kz,5)
      else
        if ( iuwvadv == 1 ) then
          call vadv(cross,aten%qc,atm1%qc,kz,6)
        else
          call vadv(cross,aten%qc,atm1%qc,kz,5)
        end if
      end if
      call pcp(jci1,jci2,ici1,ici2)
      call cldfrac(jci1,jci2,ici1,ici2)
!
!     need also to set diffq to 0 here before calling diffut
!
      adf%diffq(:,:,:) = d_zero
 
!     compute the diffusion terms:
!     the diffusion term for qv is stored in diffq. before
!     completing aten%qv computation, do not use diffq for other
!     purpose.
!
      call diffu_x(adf%diffq,atms%qvb3d,sfs%psb,xkc,kz)
      call diffu_x(aten%qc,atms%qcb3d,sfs%psb,xkc,kz)
    end if
!
    if ( ichem == 1 ) then
      !
      ! TRANSPORT OF TRACERS : initialize tracer tendencies
      !
      chiten(:,:,:,:) = d_zero
      !
      ! horizontal and vertical advection
      !
      do itr = 1 , ntr
        ! Here assignpnt does not work with gfortran with a sliced array.
        ! Doing explicit work on bounds.
        spchiten                      => chiten(:,:,:,itr)
        spchi(lbound(chi,1):,1:,1:)   => chi(:,:,:,itr)
        spchia(lbound(chia,1):,1:,1:) => chia(:,:,:,itr)
        spchib3d(lbound(chib3d,1):,1:,1:) => chib3d(:,:,:,itr)

        call hadv(cross,spchiten,spchi,kz,2)

        if ( icup /= 1 ) then
          if ( ibltyp /= 2 .and. ibltyp /= 99 ) then
            call vadv(cross,spchiten,spchia,kz,5)
          else
            if ( iuwvadv == 1 ) then
              call vadv(cross,spchiten,spchia,kz,6)
            else
              call vadv(cross,spchiten,spchia,kz,5)
            end if
          end if
        end if

        ! horizontal diffusion: initialize scratch vars to 0.
        ! need to compute tracer tendencies due to diffusion

        call diffu_x(spchiten,spchib3d,sfs%psb,xkc,kz)

      end do ! end tracer loop
      !
      ! Compute chemistry tendencies (other yhan transport)
      !
      sod = dble(idatex%second_of_day)
      call tractend2(jci1,jci2,ici1,ici2,ktau,xyear,xmonth,xday,calday,sod)
      !
    end if ! ichem
!
!----------------------------------------------------------------------
!  compute the pbl fluxes:
!  the diffusion and pbl tendencies of t and qv are stored in
!  difft and diffq.
!
    aten%u(:,:,:) = d_zero
    aten%v(:,:,:) = d_zero
!
!   calculate albedo
!
    if ( ktau == 0 .or. mod(ktau+1,ntrad) == 0 ) then
#ifdef CLM
      call albedoclm(xmonth,jci1,jci2,ici1,ici2)
#else
      call albedobats(xmonth,jci1,jci2,ici1,ici2)
#endif
    end if
!
!   call radiative transfer package
!
    if ( ktau == 0 .or. mod(ktau+1,ntrad) == 0 ) then
      loutrad = (ktau == 0 .or. mod(ktau+1,krad) == 0)
      if ( irrtm == 1 ) then
        call rrtmg_driver(jci1,jci2,ici1,ici2,xyear,eccf,loutrad)
      else
        labsem = (ktau == 0 .or. mod(ktau+1,ntabem) == 0)
        if ( labsem .and. myid == 0 ) then
          print *, 'Doing emission/absorbtion calculation...'
        end if
        call colmod3(jci1,jci2,ici1,ici2,xyear,eccf,loutrad,labsem)
      end if
    end if
 
#ifndef CLM
!
!   call mtrxbats for surface physics calculations
!
    if ( ktau == 0 .or. mod(ktau+1,ntsrf) == 0 ) then
      dtbat = dt*d_half*dble(ntsrf)
      if ( ktau == 0 ) dtbat = dt
      call mtrxbats(jci1,jci2,ici1,ici2,ktau)
    end if
#endif
 
#ifdef CLM
!
!   call mtrxclm for surface physics calculations
!
    if ( ktau == 0 .or. mod(ktau+1,ntrad) == 0 ) then
      r2cdoalb = .true.
    else
      r2cdoalb = .false.
    end if
    if ( ktau == 0 .or. mod(ktau+1,ntsrf) == 0 ) then
      ! Timestep used is the same as for bats
      if ( ktau == 0 ) then
        r2cnstep = 0
      else
        r2cnstep = (ktau+1)/ntsrf
      end if
      dtbat = dt*d_half*ntsrf
      ! CLM j loop is in mtrxclm
      call mtrxclm(ktau)
    end if
#endif

    if ( icup == 1 ) then
      wrkkuo1(1:jxp,:,:) = rsheat(:,:,:)
      wrkkuo2(1:jxp,:,:) = rswat(:,:,:)
      call deco1_exchange_left(wrkkuo1,1,ice1,ice2,1,kz)
      call deco1_exchange_right(wrkkuo1,1,ice1,ice2,1,kz)
      call deco1_exchange_left(wrkkuo2,1,ice1,ice2,1,kz)
      call deco1_exchange_right(wrkkuo2,1,ice1,ice2,1,kz)
      call htdiff(wrkkuo1,wrkkuo2,dxsq,akht1,jci1,jci2,ici1,ici2)
    end if
#ifndef BAND
    ! diagnostic on total evaporation
    if (debug_level > 2) call conqeva
#endif
!
!   Call medium resolution PBL
!
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      ! Call the Grenier and Bretherton (2001) / Bretherton (2004) TCM
      call uwtcm
      call uvcross2dot(uwten%u,uwten%v,aten%u,aten%v)
      call get_data_from_tcm(uwstateb,uwten,aten,atm1,atm2,.true.)
    end if
    if ( ibltyp == 1 .or. ibltyp == 99 ) then
      ! Call the Holtslag PBL
      call deco1_exchange_left(sfs%psb,1,ice1,ice2)
      call deco1_exchange_right(sfs%psb,1,ice1,ice2)
      call psc2psd(sfs%psb,psdot)
      call deco1_exchange_left(sfs%uvdrag,1,ice1,ice2)
      call holtbl(jci1,jci2,ici1,ici2)
    end if

    if ( ibltyp == 99 ) then
      call check_conserve_qt(holtten%qv,holtten%qc,uwten,uwstateb,kz)
      adf%diffq = adf%diffq + holtten%qv
      aten%qc = aten%qc + holtten%qc
    end if
!
!   add ccm radiative transfer package-calculated heating rates to
!   temperature tendency
!
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1, jci2
          ! heating rate in deg/sec
          aten%t(j,i,k) = aten%t(j,i,k) + sfs%psb(j,i)*heatrt(j,i,k)
        end do
      end do
    end do
!
!   add horizontal diffusion and pbl tendencies for t and qv to aten%t
!   and aten%qv for calculating condensational term in subroutine
!   "condtq".
!
    if ( ipptls == 1 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            aten%t(j,i,k) = aten%t(j,i,k) + adf%difft(j,i,k)
          end do
        end do
      end do
!
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            aten%qv(j,i,k) = aten%qv(j,i,k) + adf%diffq(j,i,k)
          end do
        end do
      end do
!
!     compute the condensation and precipitation terms for explicit
!     moisture scheme:
!
      call condtq(jci1,jci2,ici1,ici2,psc)
    end if
!
!   subtract horizontal diffusion and pbl tendencies from aten%t and
!   aten%qv for appling the sponge boundary conditions on t and qv:
!
    if ( iboudy == 4 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            aten%t(j,i,k) = aten%t(j,i,k) - adf%difft(j,i,k)
          end do
        end do
      end do
      call sponge(kz,ba_cr,xtb,aten%t)
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            aten%t(j,i,k) = aten%t(j,i,k) + adf%difft(j,i,k)
          end do
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            aten%qv(j,i,k) = aten%qv(j,i,k) - adf%diffq(j,i,k)
          end do
        end do
      end do
      call sponge(kz,ba_cr,xqb,aten%qv)
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            aten%qv(j,i,k) = aten%qv(j,i,k) + adf%diffq(j,i,k)
          end do
        end do
      end do
    end if
!
!   apply the nudging boundary conditions:
!
    if ( iboudy == 1 .or. iboudy == 5 ) then
      xtm1 = xbctime - dtsec
      if ( nbdytime == 0 .and. ktau /= 0 ) xtm1 = -dtsec
      call nudge(kz,ba_cr,xtm1,atm2%t,iboudy,xtb,aten%t)
      call nudge(kz,ba_cr,xtm1,atm2%qv,iboudy,xqb,aten%qv)
    end if

    if ( ichem == 1 ) then
      ! keep nudge_chi for now 
      if ( iboudy == 1 .or. iboudy == 5 ) then
        xtm1 = xbctime - dtsec
        if ( nbdytime == 0 .and. ktau /= 0 ) xtm1 = -dtsec
        do j = jci1 , jci2
          call nudge_chi(nspgx-1,fnudge,gnudge,xtm1,chiten(j,:,:,:),j,iboudy)
        end do
      end if
    end if
!
!   forecast t, qv, and qc at tau+1:
!
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          atmc%t(j,i,k) = atm2%t(j,i,k) + dt*aten%t(j,i,k)
          atmc%qv(j,i,k) = atm2%qv(j,i,k) + dt*aten%qv(j,i,k)
        end do
      end do
    end do
    do k = 1 , kz
      do i = icii1 , icii2
        do j = jcii1 , jcii2
          atmc%qc(j,i,k) = atm2%qc(j,i,k) + dt*aten%qc(j,i,k)
        end do
      end do
    end do
!
!   forecast tracer chi at at tau+1:
!
    if ( ichem == 1 ) then
      do itr = 1 , ntr
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              chic(j,i,k,itr) = chib(j,i,k,itr) + dt*chiten(j,i,k,itr)
            end do
          end do
        end do
      end do
    end if
!
!   compute weighted p*t (td) for use in ssi:
!
    if ( ipgf == 1 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            tvc = atmc%t(j,i,k)*(d_one+ep1*(atmc%qv(j,i,k))/psc(j,i))
            tva = atm1%t(j,i,k)*(d_one+ep1*(atmx%qv(j,i,k)))
            tvb = atm2%t(j,i,k)*(d_one+ep1* &
                                 (atm2%qv(j,i,k))/sfs%psb(j,i))
            td(j,i,k) = alpha*(tvc+tvb) + beta*tva
            ttld(j,i,k) = td(j,i,k) - psd(j,i) * &
                      t00pg*((a(k)*psd(j,i)+ptop)/p00pg)**pgfaa1
          end do
        end do
      end do
      if ( ma%hasleft ) then
        do k = 1 , kz
          do i = ici1 , ici2
            td(jce1,i,k) = atm1%t(jce1,i,k)*(d_one+ep1*(atmx%qv(jce1,i,k)))
            ttld(jce1,i,k) = td(jce1,i,k) - sfs%psa(jce1,i) * &
                        t00pg*((a(k)*sfs%psa(jce1,i)+ptop)/p00pg)**pgfaa1
          end do
        end do
      end if
      if ( ma%hasright ) then
        do k = 1 , kz
          do i = ici1 , ici2
            td(jce2,i,k) = atm1%t(jce2,i,k)*(d_one+ep1*(atmx%qv(jce2,i,k)))
            ttld(jce2,i,k) = td(jce2,i,k) - sfs%psa(jce2,i) * &
                     t00pg*((a(k)*sfs%psa(jce2,i)+ptop)/p00pg)**pgfaa1
          end do
        end do
      end if
      if ( ma%hasbottom ) then
        do k = 1 , kz
          do j = jce1 , jce2
            td(j,ice1,k) = atm1%t(j,ice1,k)*(d_one+ep1*(atmx%qv(j,ice1,k)))
            ttld(j,ice1,k) = td(j,ice1,k) - sfs%psa(j,ice1) * &
                     t00pg*((a(k)*sfs%psa(j,ice1)+ptop)/p00pg)**pgfaa1
          end do
        end do
      end if
      if ( ma%hastop ) then
        do k = 1 , kz
          do j = jce1 , jce2
            td(j,ice2,k) = atm1%t(j,ice2,k)*(d_one+ep1*(atmx%qv(j,ice2,k)))
            ttld(j,ice2,k) = td(j,ice2,k) - sfs%psa(j,ice2) * &
                     t00pg*((a(k)*sfs%psa(j,ice2)+ptop)/p00pg)**pgfaa1
          end do
        end do
      end if
    else if ( ipgf == 0 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            tvc = atmc%t(j,i,k)*(d_one+ep1*(atmc%qv(j,i,k))/psc(j,i))
            tva = atm1%t(j,i,k)*(d_one+ep1*(atmx%qv(j,i,k)))
            tvb = atm2%t(j,i,k)*(d_one+ep1* &
                 (atm2%qv(j,i,k))/sfs%psb(j,i))
            td(j,i,k) = alpha*(tvc+tvb) + beta*tva
          end do
        end do
      end do
      if ( ma%hasleft ) then
        do k = 1 , kz
          do i = ici1 , ici2
            td(jce1,i,k) = atm1%t(jce1,i,k)*(d_one+ep1*(atmx%qv(jce1,i,k)))
          end do
        end do
      end if
      if ( ma%hasright ) then
        do k = 1 , kz
          do i = ici1 , ici2
            td(jce2,i,k) = atm1%t(jce2,i,k)*(d_one+ep1*(atmx%qv(jce2,i,k)))
          end do
        end do
      end if
      if ( ma%hasbottom ) then
        do k = 1 , kz
          do j = jce1 , jce2
            td(j,ice1,k) = atm1%t(j,ice1,k)*(d_one+ep1*(atmx%qv(j,ice1,k)))
          end do
        end do
      end if
      if ( ma%hastop ) then
        do k = 1 , kz
          do j = jce1 , jce2
            td(j,ice2,k) = atm1%t(j,ice2,k)*(d_one+ep1*(atmx%qv(j,ice2,k)))
          end do
        end do
      end if
    end if
!
!----------------------------------------------------------------------
!   compute the u and v tendencies:
!
!   compute the diffusion terms:
!   put diffusion and pbl tendencies of u and v in difuu and difuv.
!
    do k = 1 , kz
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          adf%difuu(j,i,k) = aten%u(j,i,k)
          adf%difuv(j,i,k) = aten%v(j,i,k)
        end do
      end do
    end do
!
    call diffu_d(adf%difuu,atms%ubd3d,psdot,mddom%msfd,xkc,1)
    call diffu_d(adf%difuv,atms%vbd3d,psdot,mddom%msfd,xkc,1)
!
!   compute the horizontal advection terms for u and v:
!
    aten%u(:,:,:) = d_zero
    aten%v(:,:,:) = d_zero
!
    call hadv(dot,aten%u,atmx%u,kz,3)
    call hadv(dot,aten%v,atmx%v,kz,3)
!
!   compute coriolis terms:
!
    do k = 1 , kz
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          aten%u(j,i,k) = aten%u(j,i,k) + &
                       mddom%coriol(j,i)*atm1%v(j,i,k)/mddom%msfd(j,i)
          aten%v(j,i,k) = aten%v(j,i,k) - &
                       mddom%coriol(j,i)*atm1%u(j,i,k)/mddom%msfd(j,i)
        end do
      end do
    end do
!
!   compute pressure gradient terms:
!
    if ( ipgf == 1 ) then
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            psasum = psd(j,i) + psd(j,i-1) + psd(j-1,i) + psd(j-1,i-1)
            sigpsa = psasum
            tv1 = atmx%t(j-1,i-1,k)*(d_one+ep1*(atmx%qv(j-1,i-1,k)))
            tv2 = atmx%t(j-1,i,k)*(d_one+ep1*(atmx%qv(j-1,i,k)))
            tv3 = atmx%t(j,i-1,k)*(d_one+ep1*(atmx%qv(j,i-1,k)))
            tv4 = atmx%t(j,i,k)*(d_one+ep1*(atmx%qv(j,i,k)))
            rtbar = tv1 + tv2 + tv3 + tv4 - d_four*t00pg*             &
                    ((a(k)*psasum*d_rfour+ptop)/p00pg)**pgfaa1
            rtbar = rgas*rtbar*sigpsa/16.0D0
            aten%u(j,i,k) = aten%u(j,i,k) - rtbar * &
                  (dlog(d_half*(psd(j,i)+psd(j,i-1))*a(k)+ptop) -     &
                   dlog(d_half*(psd(j-1,i)+psd(j-1,i-1))*a(k)+ptop))/ &
                   (dx*mddom%msfd(j,i))
            aten%v(j,i,k) = aten%v(j,i,k) - rtbar * &
                  (dlog(d_half*(psd(j,i)+psd(j-1,i))*a(k)+ptop) -     &
                   dlog(d_half*(psd(j-1,i-1)+psd(j,i-1))*a(k)+ptop))/ &
                   (dx*mddom%msfd(j,i))
          end do
        end do
      end do
    else if ( ipgf == 0 ) then
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            psasum = psd(j,i) + psd(j,i-1) + psd(j-1,i) + psd(j-1,i-1)
            sigpsa = psasum
            tv1 = atmx%t(j-1,i-1,k)*(d_one+ep1*(atmx%qv(j-1,i-1,k)))
            tv2 = atmx%t(j-1,i,k)*(d_one+ep1*(atmx%qv(j-1,i,k)))
            tv3 = atmx%t(j,i-1,k)*(d_one+ep1*(atmx%qv(j,i-1,k)))
            tv4 = atmx%t(j,i,k)*(d_one+ep1*(atmx%qv(j,i,k)))
            rtbar = rgas*(tv1+tv2+tv3+tv4)*sigpsa/16.0D0
            aten%u(j,i,k) = aten%u(j,i,k) - rtbar * &
                   (dlog(d_half*(psd(j,i)+psd(j,i-1))*a(k)+ptop) -    &
                    dlog(d_half*(psd(j-1,i)+psd(j-1,i-1))*a(k)+ptop))/&
                    (dx*mddom%msfd(j,i))
            aten%v(j,i,k) = aten%v(j,i,k) - rtbar *                   &
                   (dlog(d_half*(psd(j,i)+psd(j-1,i))*a(k)+ptop) -    &
                    dlog(d_half*(psd(j-1,i-1)+psd(j,i-1))*a(k)+ptop))/&
                    (dx*mddom%msfd(j,i))
          end do
        end do
      end do
    end if
!
!   compute geopotential height at half-k levels, cross points:
!
    if ( ipgf == 1 ) then
      do i = ice1 , ice2
        do j = jce1 , jce2
          tv = (ttld(j,i,kz)/psd(j,i))/(d_one+atmx%qc(j,i,kz)/ &
                                       (d_one+atmx%qv(j,i,kz)))
          phi(j,i,kz) = mddom%ht(j,i) + &
                   rgas*t00pg/pgfaa1*((psd(j,i)+ptop)/p00pg)**pgfaa1
          phi(j,i,kz) = phi(j,i,kz) - rgas * &
                  tv*dlog((a(kz)+ptop/psd(j,i))/(d_one+ptop/psd(j,i)))
        end do
      end do
      do k = 1 , kzm1
        lev = kz - k
        do i = ice1 , ice2
          do j = jce1 , jce2
            tvavg = ((ttld(j,i,lev)*dsigma(lev)+ttld(j,i,lev+1)* &
                    dsigma(lev+1))/(psd(j,i)*(dsigma(lev)+       &
                    dsigma(lev+1))))/(d_one+atmx%qc(j,i,lev)/    &
                    (d_one+atmx%qv(j,i,lev)))
            phi(j,i,lev) = phi(j,i,lev+1) - rgas *    &
                   tvavg*dlog((a(lev)+ptop/psd(j,i))/ &
                             (a(lev+1)+ptop/psd(j,i)))
          end do
        end do
      end do
    else if ( ipgf == 0 ) then
      do i = ice1 , ice2
        do j = jce1 , jce2
          tv = (td(j,i,kz)/psd(j,i))/(d_one+atmx%qc(j,i,kz)/  &
               (d_one+atmx%qv(j,i,kz)))
          phi(j,i,kz) = mddom%ht(j,i) - rgas * &
               tv*dlog((a(kz)+ptop/psd(j,i))/(d_one+ptop/psd(j,i)))
        end do
      end do
      do k = 1 , kzm1
        lev = kz - k
        do i = ice1 , ice2
          do j = jce1 , jce2
            tvavg = ((td(j,i,lev)*dsigma(lev)+td(j,i,lev+1)*   &
                    dsigma(lev+1))/(psd(j,i)*(dsigma(lev)+     &
                    dsigma(lev+1))))/(d_one+atmx%qc(j,i,lev)/  &
                    (d_one+atmx%qv(j,i,lev)))
            phi(j,i,lev) = phi(j,i,lev+1) - rgas *    &
                   tvavg*dlog((a(lev)+ptop/psd(j,i))  &
                           /(a(lev+1)+ptop/psd(j,i)))
          end do
        end do
      end do
    end if
!
    call deco1_exchange_left(phi,1,ide1,ide2,1,kz)
!
!   compute the geopotential gradient terms:
!
    do k = 1 , kz
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          aten%u(j,i,k) = aten%u(j,i,k) -                              &
               (psd(j-1,i-1)+psd(j-1,i)+psd(j,i-1)+psd(j,i)) *         &
               (phi(j,i,k)+phi(j,i-1,k)-phi(j-1,i,k)-phi(j-1,i-1,k)) / &
               (dx8*mddom%msfd(j,i))
          aten%v(j,i,k) = aten%v(j,i,k) -                              &
               (psd(j-1,i-1)+psd(j-1,i)+psd(j,i-1)+psd(j,i)) *         &
               (phi(j,i,k)+phi(j-1,i,k)-phi(j,i-1,k)-phi(j-1,i-1,k)) / &
               (dx8*mddom%msfd(j,i))
        end do
      end do
    end do
!
!   compute the vertical advection terms:
!
    call vadv(dot,aten%u,atm1%u,kz,4)
    call vadv(dot,aten%v,atm1%v,kz,4)
!
!   apply the sponge boundary condition on u and v:
!
    if ( iboudy == 4 ) then
      call sponge(kz,ba_dt,xub,aten%u)
      call sponge(kz,ba_dt,xvb,aten%v)
    end if
!
!   apply the nudging boundary conditions:
!
    if ( iboudy == 1 .or. iboudy == 5 ) then
      xtm1 = xbctime - dtsec
      if ( nbdytime == 0 .and. ktau /= 0 ) xtm1 = -dtsec
      call nudge(kz,ba_dt,xtm1,atm2%u,iboudy,xub,aten%u)
      call nudge(kz,ba_dt,xtm1,atm2%v,iboudy,xvb,aten%v)
    end if
!
!   add the diffusion and pbl tendencies to aten%u and aten%v:
!
    do k = 1 , kz
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          aten%u(j,i,k) = aten%u(j,i,k) + adf%difuu(j,i,k)
          aten%v(j,i,k) = aten%v(j,i,k) + adf%difuv(j,i,k)
        end do
      end do
    end do
!
!   forecast p*u and p*v at tau+1:
!
    do k = 1 , kz
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          atmc%u(j,i,k) = atm2%u(j,i,k) + dt*aten%u(j,i,k)
          atmc%v(j,i,k) = atm2%v(j,i,k) + dt*aten%v(j,i,k)
        end do
      end do
    end do
!
!   Couple TKE to ps for use in vertical advection
!
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            uwstatea%tkeps(j,i,k) = atm1%tke(j,i,k)*sfs%psa(j,i)
            uwstatea%advtke(j,i,k) = d_zero
          end do
        end do
      end do
      do i = ice1 , ice2
        do j = jce1 , jce2
          uwstatea%tkeps(j,i,kz+1) = atm1%tke(j,i,kz+1)*sfs%psa(j,i)
        end do
      end do
    end if

    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      do j = jci1 , jci2
        ! Calculate the horizontal advective tendency for TKE
        call hadvtke(uwstatea,atm1,twt,dx4,j)
        ! Calculate the vertical advective tendency for TKE
        call vadvtke(uwstatea,qdot,j,2)
        ! Calculate the horizontal, diffusive tendency for TKE
      end do
      call diffu_x(uwstatea%advtke,atms%tkeb3d,sfs%psb,xkcf,kzp1)
    end if
!
!   store the xxa variables in xxb and xxc in xxa:
!   perform time smoothing operations.
!
    do k = 1 , kz
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          atm2%u(j,i,k) = omuhf*atm1%u(j,i,k)/mddom%msfd(j,i) + &
                          gnuhf*(atm2%u(j,i,k)+atmc%u(j,i,k))
          atm2%v(j,i,k) = omuhf*atm1%v(j,i,k)/mddom%msfd(j,i) + &
                          gnuhf*(atm2%v(j,i,k)+atmc%v(j,i,k))
          atm1%u(j,i,k) = atmc%u(j,i,k)
          atm1%v(j,i,k) = atmc%v(j,i,k)
          ! TAO: Once the full loop above is completed, update the TKE
          ! tendency if the UW PBL is running.  NOTE!!! Do not try to
          ! combine these loops with the above loop Advection MUST be
          ! done in a loop separate from the updates.  (I lost 3 days
          ! of working to disocover that this is a problem because I
          ! thought it would be clever to combine loops--TAO)
          if ( ibltyp == 2 .or. ibltyp == 99 ) then
            ! Add the advective tendency to the TKE tendency calculated
            ! by the UW TKE routine
             aten%tke(j,i,k) = aten%tke(j,i,k) + &
                               uwstatea%advtke(j,i,k)/sfs%psa(j,i)
             ! Do a filtered time integration
             atmc%tke(j,i,k) = max(tkemin,atm2%tke(j,i,k) + &
                               dttke*aten%tke(j,i,k))
             atm2%tke(j,i,k) = max(tkemin,omuhf*atm1%tke(j,i,k) + &
                               gnuhf*(atm2%tke(j,i,k) + atmc%tke(j,i,k)))
             atm1%tke(j,i,k) = atmc%tke(j,i,k)
          end if ! TKE tendency update
        end do
      end do
    end do
!
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          atm2%t(j,i,k) = omuhf*atm1%t(j,i,k) + &
                          gnuhf*(atm2%t(j,i,k)+atmc%t(j,i,k))
          atm1%t(j,i,k) = atmc%t(j,i,k)
          qvas = atmc%qv(j,i,k)
          if ( qvas < dlowval ) qvas = minqx
          qvbs = omuhf*atm1%qv(j,i,k) + gnuhf*(atm2%qv(j,i,k)+atmc%qv(j,i,k))
          if ( qvbs < dlowval ) qvbs = minqx
          atm2%qv(j,i,k) = qvbs
          atm1%qv(j,i,k) = qvas
          qcas = atmc%qc(j,i,k)
          if ( qcas < dlowval ) qcas = d_zero
          qcbs = omu*atm1%qc(j,i,k) + gnu*(atm2%qc(j,i,k)+atmc%qc(j,i,k))
          if ( qcbs < dlowval ) qcbs = d_zero
          atm2%qc(j,i,k) = qcbs
          atm1%qc(j,i,k) = qcas
        end do
      end do
    end do
    if ( ichem == 1 ) then
      do itr = 1 , ntr
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              chias = chic(j,i,k,itr)
              if ( chias < dlowval ) chias = dlowval
              chibs = omu*chia(j,i,k,itr)                             &
                      + gnu*(chib(j,i,k,itr)+chic(j,i,k,itr))
              if ( chibs < dlowval ) chibs = dlowval
              chib(j,i,k,itr) = chibs
              chia(j,i,k,itr) = chias
            end do
          end do
        end do
      end do
    end if

    do i = ici1 , ici2
      do j = jci1 , jci2
        sfs%psb(j,i) = omuhf*sfs%psa(j,i) + gnuhf*(sfs%psb(j,i)+psc(j,i))
        sfs%psa(j,i) = psc(j,i)
      end do
    end do
!
!   increment elapsed forecast time:
!
    ktau = ktau + 1
    xbctime = xbctime + dtsec
    nbdytime = nbdytime + ntsec
    idatex = idatex + intmdl

    if ( mod(ktau,khour) == 0 ) then
      call split_idate(idatex,xyear,xmonth,xday,xhour)
    end if
    if ( mod(ktau,kbdy) == 0 ) then
      nbdytime = 0
      xbctime = d_zero
    end if

    if ( ktau == 2 ) then
      dt = dt2
      dtcum  = dt2
      dtche  = dt2
      dtpbl  = dt2
      rdtpbl = d_one/dt2
      dttke  = dt2
    end if
    !
    ! compute the amounts advected through the lateral boundaries:
    ! *** note *** we must calculate the amounts advected through
    ! the lateral boundaries before updating the values
    ! at boundary slices.
    !
#ifndef BAND
    if (debug_level > 2) then
      call conadv
      if ( ichem == 1 ) call tracdiag(xkc)
    end if
#endif
    !
    ! do cumulus transport of tracers
    !
    if ( ichem == 1 .and. ichcumtra == 1 ) call cumtran
#ifndef BAND
    ! 
    ! trace the mass conservation of dry air and water substance:
    !
    if (debug_level > 2) call conmas
#endif
    !
    ! budgets for tracers
    !
    if ( ichem == 1 ) then
      call tracbud
#ifndef BAND
      if (debug_level > 2) call contrac
#endif
    end if
    !
    ! Print out noise parameter
    !
    if ( ktau > 1 ) then
      ptnbar = ptntot/dble(iptn)
      pt2bar = pt2tot/dble(iptn)
      icons_mpi = 0
      call mpi_allreduce(total_precip_points,icons_mpi,1,mpi_integer, &
                         mpi_sum,mycomm,ierr)
      ! Added a check for nan... The following inequality is wanted.
      if ((ptnbar /= ptnbar) .or. &
         ((ptnbar > d_zero) .eqv. (ptnbar <= d_zero))) then
        maxv = maxval(aten%t)
        if ( maxv > d_one ) then
          print *, 'MAXVAL ATEN T :', maxv
          maxv = maxv - 0.001D0
          do kk = 1 , kz
            do ii = ici1 , ici2
              do jj = jci1 , jci2
                if ( aten%t(jj,ii,kk) > maxv ) then
                  print *, 'II :', ii , ', JJ :', myid*jxp+jj , ', KK :', kk
                end if
              end do
            end do
          end do
        end if
        if ( myid == 0 ) then
          write (*,*) 'WHUUUUBBBASAAAGASDDWD!!!!!!!!!!!!!!!!'
          write (*,*) 'No more atmosphere here....'
          write (*,*) 'CFL violation detected, so model STOP'
          write (*,*) '#####################################'
          write (*,*) '#            DECREASE DT !!!!       #'
          write (*,*) '#####################################'
          call fatal(__FILE__,__LINE__,'CFL VIOLATION')
        end if
      end if
      if ( myid == 0 ) then
        if ( mod(ktau,irep) == 0 ) then
          appdat = tochar(idatex)
          write(6,99001) appdat , ktau , ptnbar , pt2bar , icons_mpi
        end if
      end if

99001 format (a23,', ktau = ',i10, ' :  1st, 2nd time deriv of ps = ',2E12.5, &
             ',  no. of points w/convection = ',i7)
    end if
!
    call time_end(subroutine_name,idindx)
!
  end subroutine tend
!
end module mod_tendency
