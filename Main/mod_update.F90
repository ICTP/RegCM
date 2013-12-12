!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     This file is part of ICTP RegCM.
!
!     ICTP RegCM is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     ICTP RegCM is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::      
!
      module mod_update
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
       use mod_intkinds, only : ik4
       use mod_realkinds, only : rk8
       use mod_regcm_types
       use mod_memutil
!
      implicit none
      private

      type(imp_data), public :: importFields
      type(exp_data), public :: exportFields
!
      integer(ik4), pointer :: ldmskb(:,:)
      integer(ik4), pointer :: wetdry(:,:)
!
      real(rk8), parameter :: zeroval = 0.0d20
      real(rk8), parameter :: missing_r8 = 1.0d20
      real(rk8), parameter :: tol = missing_r8/2.0d0
!
!-----------------------------------------------------------------------
!     Public subroutines 
!-----------------------------------------------------------------------
!
      public :: RCM_Get
      public :: RCM_Put
      public :: RCM_Allocate
!
!-----------------------------------------------------------------------
!     Module constants 
!-----------------------------------------------------------------------
!
      real(rk8), parameter :: beta = 1.25 ! gustiness coeff
      real(rk8), parameter :: von  = 0.4  ! von Karman constant
      real(rk8), parameter :: fdg  = 1.00 ! ratio of thermal to wind von Karman
      real(rk8), parameter :: tdk  = 273.16
      real(rk8), parameter :: grav = 9.82 ! accel of earth grav
!
      contains
!
      subroutine RCM_Allocate()
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use mod_bats_common, only : cplmsk
      use mod_atm_interface , only : mddom
      use mod_dynparam, only : ice1, ice2, jce1, jce2
      use mod_dynparam, only : ici1, ici2, jci1, jci2
!
      implicit none
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: i, j
      real(rk8), parameter :: initval = 1.0d20
      real(rk8), parameter :: zeroval = 0.0d20
!
!-----------------------------------------------------------------------
!     Allocate arrays 
!-----------------------------------------------------------------------
!
      call getmem2d(exportFields%psfc,jce1,jce2,ice1,ice2,'cpl:psfc')
      call getmem2d(exportFields%tsfc,jce1,jce2,ice1,ice2,'cpl:tsfc')
      call getmem2d(exportFields%qsfc,jce1,jce2,ice1,ice2,'cpl:qsfc')
      call getmem2d(exportFields%swrd,jce1,jce2,ice1,ice2,'cpl:swrd')
      call getmem2d(exportFields%lwrd,jce1,jce2,ice1,ice2,'cpl:lwrd')
      call getmem2d(exportFields%dlwr,jce1,jce2,ice1,ice2,'cpl:dlwr')
      call getmem2d(exportFields%lhfx,jce1,jce2,ice1,ice2,'cpl:lhfx')
      call getmem2d(exportFields%shfx,jce1,jce2,ice1,ice2,'cpl:shfx')
      call getmem2d(exportFields%prec,jce1,jce2,ice1,ice2,'cpl:prec')
      call getmem2d(exportFields%wndu,jce1,jce2,ice1,ice2,'cpl:wndu')
      call getmem2d(exportFields%wndv,jce1,jce2,ice1,ice2,'cpl:wndv')
      call getmem2d(exportFields%rnof,jci1,jci2,ici1,ici2,'cpl:rnof')
      call getmem2d(exportFields%snof,jci1,jci2,ici1,ici2,'cpl:snof')
      call getmem2d(exportFields%taux,jce1,jce2,ice1,ice2,'cpl:taux')
      call getmem2d(exportFields%tauy,jce1,jce2,ice1,ice2,'cpl:tauy')
      call getmem2d(exportFields%wspd,jce1,jce2,ice1,ice2,'cpl:wspd')
      call getmem2d(exportFields%nflx,jce1,jce2,ice1,ice2,'cpl:nflx')
      call getmem2d(exportFields%sflx,jce1,jce2,ice1,ice2,'cpl:sflx')
      call getmem2d(exportFields%snow,jce1,jce2,ice1,ice2,'cpl:snow')
      call getmem2d(exportFields%dswr,jce1,jce2,ice1,ice2,'cpl:dswr')
!
      call getmem2d(importFields%sst,jce1,jce2,ice1,ice2,'cpl:sst')
      call getmem2d(importFields%sit,jce1,jce2,ice1,ice2,'cpl:sit')
      call getmem2d(importFields%msk,jce1,jce2,ice1,ice2,'cpl:msk')
!
      call getmem2d(ldmskb,jci1,jci2,ici1,ici2,'cpl:ldmsk')
      call getmem2d(wetdry,jci1,jci2,ici1,ici2,'cpl:wetdry')
!
!-----------------------------------------------------------------------
!     Initialize arrays 
!-----------------------------------------------------------------------
!
      do i = ici1, ici2
        do j = jci1, jci2
          exportFields%psfc(j,i) = initval
          exportFields%tsfc(j,i) = initval
          exportFields%qsfc(j,i) = initval
          exportFields%swrd(j,i) = initval
          exportFields%lwrd(j,i) = initval
          exportFields%dlwr(j,i) = initval
          exportFields%lhfx(j,i) = initval
          exportFields%shfx(j,i) = initval
          exportFields%prec(j,i) = initval
          exportFields%wndu(j,i) = initval
          exportFields%wndv(j,i) = initval
          exportFields%rnof(j,i) = zeroval
          exportFields%snof(j,i) = zeroval
          exportFields%taux(j,i) = zeroval
          exportFields%tauy(j,i) = zeroval
          exportFields%wspd(j,i) = zeroval
          exportFields%nflx(j,i) = zeroval
          exportFields%sflx(j,i) = zeroval
          exportFields%snow(j,i) = zeroval
          exportFields%dswr(j,i) = zeroval
!
          importFields%sst(j,i) = initval
          importFields%sit(j,i) = initval
          importFields%msk(j,i) = initval 
!
          ldmskb(j,i) = mddom%ldmsk(j,i)
          wetdry(j,i) = 0
        end do
      end do
!
      end subroutine RCM_Allocate
!
      subroutine RCM_Get(localPet)
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use mod_constants
      use mod_atm_interface, only : sfs, mddom , mdsub
      use mod_bats_common, only : cplmsk, ntcpl, tgbrd
      use mod_bats_common, only : sfice, sent, sncv, tgrd
      use mod_dynparam, only : ice1, ice2, jce1, jce2
      use mod_dynparam, only : ici1, ici2, jci1, jci2, nnsg
      use mod_dynparam, only : global_cross_istart, global_cross_jstart
      use mod_runparams, only : iswater, ktau
!
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      integer, intent(in) :: localPet
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: i, j, ii, jj, n
      logical :: flag
      real(rk8) :: toth
      real(rk8), parameter :: iceminh = d_10
      real(rk8), parameter :: href = d_two * iceminh
      real(rk8), parameter :: steepf = 1.0D0
      integer(ik4) :: ix, jy, imin, imax, jmin, jmax, srad, hveg(22)
!
!-----------------------------------------------------------------------
!     Retrieve information from OCN component
!-----------------------------------------------------------------------
!
      do i = ici1, ici2
      ii = global_cross_istart+i-1
      do j = jci1, jci2
      jj = global_cross_jstart+j-1
!        
      if (iswater(mddom%lndcat(j,i))) then
!
!-----------------------------------------------------------------------
!     Update: Sea Surface Temperature (SST)
!-----------------------------------------------------------------------
!
      if (importFields%sst(j,i) .lt. tol) then
        ! create fixed coupling mask
        if (mod(ktau+1, ntcpl*2) == 0) then
          cplmsk(j,i) = 1
        end if
!        
        sfs%tga(j,i) = importFields%sst(j,i)
        sfs%tgb(j,i) = importFields%sst(j,i)
        sfs%tgbb(j,i) = importFields%sst(j,i)
        tgrd(:,j,i)  = importFields%sst(j,i)
        tgbrd(:,j,i) = importFields%sst(j,i)
      end if
!
!-----------------------------------------------------------------------
!     Update: Mask and land-use type (based on dynamic wet-dry)
!-----------------------------------------------------------------------
!
!      if (importFields%msk(j,i) .lt. tol .and. ldmskb(j,i) == 0) then
!        if (importFields%msk(j,i) .lt. 1.0) then
!          flag = .false.
!          if (mddom%ldmsk(j,i) == 0 .or. mddom%ldmsk(j,i) == 2) flag = .true.
!          ! set land-sea mask
!          mddom%ldmsk(j,i) = 1
!          do n = 1, nnsg
!            mdsub%ldmsk(n,j,i) = mddom%ldmsk(j,i)
!          end do
!          ! count land-use type in a specified search radius (srad)
!          srad = 10
!          jmin = j-srad
!          if (j-srad < jci1) jmin = jci1
!          jmax = j+srad
!          if (j+srad > jci2) jmax = jci2
!          imin = i-srad
!          if (i-srad < ici1) imin = ici1
!          imax = i+srad
!          if (i+srad > ici2) imax = ici2
!          hveg = 0
!          do ix = imin, imax
!            do jy = jmin, jmax
!              do n = 1, nnsg
!                hveg(mdsub%iveg(n,jy,ix)) = hveg(mdsub%iveg(n,jy,ix))+1
!              end do
!            end do
!          end do        
!          hveg(14) = 0
!          hveg(15) = 0
!          ! set array to store change
!          wetdry(j,i) = 1
!          ! write debug info
!          if (flag) then
!            write(*,20) jj, ii, 'water', 'land ',                       &
!                        mddom%ldmsk(j,i)
!          end if
!        else
!          if (mddom%ldmsk(j,i) == 1 .and. wetdry(j,i) == 1) then
!            flag = .false.
!            if (ldmskb(j,i) /= mddom%ldmsk(j,i)) flag = .true.
!            ! set land-sea mask to its original value
!            mddom%ldmsk(j,i) = ldmskb(j,i)             
!            do n = 1, nnsg
!              mdsub%ldmsk(n,j,i) = mddom%ldmsk(j,i)
!            end do
!            ! set array to store change
!            wetdry(j,i) = 0
!            ! write debug info
!            if (flag) then
!              write(*,20) jj, ii, 'land ', 'water',                     &
!                        mddom%ldmsk(j,i)
!            end if
!          end if
!        end if
!      end if
!
!-----------------------------------------------------------------------
!     Update: Sea-ice, mask and land-use type (based on sea-ice module) 
!-----------------------------------------------------------------------
! 
      if (importFields%sit(j,i) < tol .and. mddom%ldmsk(j,i) /= 1) then
        if (importFields%sit(j,i) > iceminh) then
          flag = .false.
          if (mddom%ldmsk(j,i) == 0) flag = .true.
          ! set land-sea mask
          mddom%ldmsk(j,i) = 2
          do n = 1, nnsg
            mdsub%ldmsk(n,j,i) = 2
            ! set sea ice thikness (in mm)
            sfice(n,j,i) = importFields%sit(j,i) 
          end do
          ! write debug info
          if (flag) then
            write(*,30) jj, ii, 'water', 'ice  ',                       &
                        mddom%ldmsk(j,i), sfice(1,j,i)
          end if
        else
          if (ldmskb(j,i) == 0 .and. mddom%ldmsk(j,i) == 2) then
            ! reduce to one tenth surface ice: it should melt away
            do n = 1, nnsg
              ! check that sea ice is melted or not
              if (sfice(n,j,i) <= iceminh) then
                if (ldmskb(j,i) /= mddom%ldmsk(j,i)) flag = .true.
                ! set land-sea mask to its original value
                mddom%ldmsk(j,i) = ldmskb(j,i)
                mdsub%ldmsk(n,j,i) = ldmskb(j,i)
                ! set land-use type to its original value
                ! set sea ice thikness (in mm)
                sfice(n,j,i) = d_zero 
              else
                flag = .false.
              end if
            end do
            ! write debug info
            if (flag) then
              write(*,40) jj, ii, 'ice  ', 'water',                     &
                          mddom%ldmsk(j,i), sfice(1,j,i)
            end if 
          end if           
        end if
      end if 
      end if      
!
      end do
      end do
!
!-----------------------------------------------------------------------
!     Format definition 
!-----------------------------------------------------------------------
!
 20   format(' ATM land-sea mask is changed at (',I3,',',I3,') : ',     &
             A5,' --> ',A5,' [',I2,']')
 30   format(' ATM sea-ice is formed at (',I3,',',I3,') : ',            &
             A5,' --> ',A5,' ['I2,' - ',F12.4,']')
 40   format(' ATM sea-ice is melted at (',I3,',',I3,') : ',       &
             A5,' --> ',A5,' [',I2,' - ',F12.4,']')
!
      end subroutine RCM_Get
!
      subroutine RCM_Put(localPet)
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use mod_constants
      use mod_runparams, only : kday, ktau, rnsrf_for_day
      use mod_dynparam, only : ici1, ici2, jci1, jci2, nnsg, ptop
      use mod_bats_common, only : sfps, t2m, q2m, evpr, sent, prcp, &
                                  u10m, v10m, srnof, trnof, rdnnsg,  &
                                  dailyrnf, taux, tauy, sncv, runoffcount
      use mod_atm_interface, only : flw, flwd, fsw , mddom , solvs, &
                                    solvsd
!
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      integer, intent(in) :: localPet
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: i, j 
!
!-----------------------------------------------------------------------
!     Send information to OCN component
!-----------------------------------------------------------------------
!
      do i = ici1, ici2
        do j = jci1, jci2
          exportFields%psfc(j,i) = (sfps(j,i)+ptop)*d_10
          exportFields%tsfc(j,i) = t2m(1,j,i)
          exportFields%qsfc(j,i) = q2m(1,j,i)
          exportFields%swrd(j,i) = fsw(j,i)
          exportFields%lwrd(j,i) = flw(j,i) 
          exportFields%dlwr(j,i) = flwd(j,i)
          exportFields%lhfx(j,i) = evpr(1,j,i)*wlhv
          exportFields%shfx(j,i) = sent(1,j,i)
          exportFields%prec(j,i) = prcp(1,j,i)
          exportFields%wndu(j,i) = u10m(1,j,i)
          exportFields%wndv(j,i) = v10m(1,j,i)
          exportFields%taux(j,i) = taux(1,j,i)
          exportFields%tauy(j,i) = tauy(1,j,i)
          exportFields%wspd(j,i) = dsqrt(u10m(1,j,i)**2+v10m(1,j,i)**2)
          exportFields%nflx(j,i) = fsw(j,i)-evpr(1,j,i)*wlhv-sent(1,j,i)-flw(j,i)
          exportFields%sflx(j,i) = evpr(1,j,i)-prcp(1,j,i)
          exportFields%snow(j,i) = sncv(1,j,i)
          exportFields%dswr(j,i) = solvsd(j,i)+solvs(j,i) 
        end do
      end do
!
!-----------------------------------------------------------------------
!     Send information to RTM component
!-----------------------------------------------------------------------
!
      if (mod(ktau+1,kday) == 0) then
        where (mddom%ldmsk > 0)
          exportFields%rnof = dailyrnf(:,:,1)/runoffcount
          exportFields%snof = dailyrnf(:,:,2)/runoffcount
        else where
          exportFields%rnof = zeroval
          exportFields%snof = zeroval
        end where
        dailyrnf(:,:,:) = zeroval
        runoffcount = zeroval
      end if
!
      end subroutine RCM_Put
!
      end module mod_update
