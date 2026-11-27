!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_massck
  !
  ! Daily mass of air and water of the hydrostatic and the non hydrostatic
  ! sigma cores (debug_level > 0). The MOLOCH core has its own budget.
  !
  use mod_intkinds, only : ik4
  use mod_realkinds, only : rk16, rk8
  use mod_stdio, only : stdout
  use mod_dynparam, only : ici1, ici2, jci1, jci2, kz, myid, nspgx
  use mod_constants, only : d_zero, egrav, regrav
  use mod_runparams, only : alarm_day, dsigma, dtsec, dxsq, iqlst, rcmtimer, &
    syncro_dbg
  use mod_mppparam, only : italk, ma, sumall
  use mod_atm_interface, only : atm1, crrate, mddom, ncrrate, sfs

  implicit none

  private

  public :: massck

#ifndef QUAD_PRECISION
  integer, parameter :: wrkp = rk8
#else
  integer, parameter :: wrkp = rk16
#endif

  real(wrkp), parameter :: q_zero = 0.0_wrkp
  real(wrkp), public :: dryini = q_zero
  real(wrkp), public :: watini = q_zero
  real(rk8), public :: dryerror = d_zero
  real(rk8), public :: waterror = d_zero

  ! Daily sums for the report (kg), and masses at the start of the day
  real(wrkp) :: scrai = 0.0_wrkp
  real(wrkp) :: sevap = 0.0_wrkp
  real(wrkp) :: dry0 = -1.0_wrkp
  real(wrkp) :: wat0 = -1.0_wrkp
  real(wrkp) :: all0 = -1.0_wrkp
  real(wrkp) :: area = -1.0_wrkp
  real(wrkp) :: allarea = -1.0_wrkp

  contains
  !
  ! Daily report of the mass of air (hPa) and water (mm) of the points
  ! inside the lateral relaxation zone (nspgx points from the boundary,
  ! changed only by the dynamics and the physics), with the evaporation
  ! and the precipitation, and of the air of all the internal points.
  ! The cell area is (dx/m)**2, the pressures are in cb.
  !
  subroutine massck
    implicit none
    real(wrkp) :: tdrym, tqmass, tcrai, tqeva, tarea, tall, tallarea
    real(wrkp) :: drymass, qmass, craim, evapm, allmass, csa
    real(rk8) :: c
    integer(ik4) :: i, j, k, n, i1, i2, j1, j2

    j1 = jci1
    j2 = jci2
    i1 = ici1
    i2 = ici2
    if ( ma%has_bdyleft ) j1 = jci1 + nspgx
    if ( ma%has_bdyright ) j2 = jci2 - nspgx
    if ( ma%has_bdybottom ) i1 = ici1 + nspgx
    if ( ma%has_bdytop ) i2 = ici2 - nspgx
    tdrym = q_zero
    tqmass = q_zero
    tcrai = q_zero
    tqeva = q_zero
    tarea = q_zero
    tall = q_zero
    tallarea = q_zero
    do i = ici1, ici2
      do j = jci1, jci2
        csa = dxsq/(mddom%msfx(j,i)*mddom%msfx(j,i))
        tallarea = tallarea + csa
        tall = tall + sfs%psa(j,i)*csa*1000.0_wrkp*regrav
        if ( i < i1 .or. i > i2 .or. j < j1 .or. j > j2 ) cycle
        tarea = tarea + csa
        tdrym = tdrym + sfs%psa(j,i)*csa*1000.0_wrkp*regrav
        tcrai = tcrai + (crrate(j,i)+ncrrate(j,i))*csa*dtsec
        tqeva = tqeva + sfs%qfx(j,i)*csa*dtsec
        do n = 1, iqlst
          do k = 1, kz
            tqmass = tqmass + atm1%qx(j,i,k,n) * &
              dsigma(k)*csa*1000.0_wrkp*regrav
          end do
        end do
      end do
    end do
    call sumall(tdrym,drymass)
    call sumall(tqmass,qmass)
    call sumall(tall,allmass)
    if ( area < d_zero ) then
      call sumall(tarea,area)
      call sumall(tallarea,allarea)
    end if
    call sumall(tcrai,craim)
    call sumall(tqeva,evapm)

    if ( myid == italk ) then
      if ( dry0 < d_zero .or. rcmtimer%start( ) ) then
        dry0 = drymass
        wat0 = qmass
        all0 = allmass
        scrai = q_zero
        sevap = q_zero
      else
        scrai = scrai + craim
        sevap = sevap + evapm
      end if
      if ( .not. rcmtimer%start( ) .and. &
           ( alarm_day%act( ) .or. syncro_dbg%act( ) .or. &
             rcmtimer%next_is_endtime ) ) then
        c = real(egrav,rk8)/100.0_rk8
        write(stdout,'(a,a)') ' MASS BUDGET (hPa/day) to ', &
          trim(rcmtimer%str( ))
        write(stdout,'(a,f10.4)') '   mass (hPa)          = ', &
          real(allmass/allarea,rk8)*c
        write(stdout,'(a,f10.4)') '   change              = ', &
          real((allmass-all0)/allarea,rk8)*c
        write(stdout,'(a,f10.4)') '   mass inside the relaxation zone = ', &
          real(drymass/area,rk8)*c
        write(stdout,'(a,f10.4)') '   change              = ', &
          real((drymass-dry0)/area,rk8)*c
        write(stdout,'(a,a)') &
          ' WATER BUDGET inside the relaxation zone (mm/day) to ', &
          trim(rcmtimer%str( ))
        write(stdout,'(a,f10.4)') '   water (mm)          = ', &
          real(qmass/area,rk8)
        write(stdout,'(a,f10.4)') '   change              = ', &
          real((qmass-wat0)/area,rk8)
        write(stdout,'(a,f10.4)') '   evaporation         = ', &
          real(sevap/area,rk8)
        write(stdout,'(a,f10.4)') '   precipitation       = ', &
          real(scrai/area,rk8)
        dry0 = drymass
        wat0 = qmass
        all0 = allmass
        scrai = q_zero
        sevap = q_zero
      end if
      ! Saved in the restart file
      dryini = drymass
      watini = qmass
      dryerror = d_zero
      waterror = d_zero
    end if
  end subroutine massck

end module mod_massck

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
