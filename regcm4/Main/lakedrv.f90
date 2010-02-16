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
 
      subroutine lakedrv(jslc)
      use mod_regcm_param
      use mod_param1
      use mod_iunits
      use mod_bats
      use mod_date
      implicit none
!
! Dummy arguments
!
      integer :: jslc
      intent (in) jslc
!
! Local variables
!
      real(8) :: aveice , dayl , evl , flw , fsw , hlat , hsen , prec , &
               & ql , snow , tgl , tl , vl , zl
      integer :: ilake , n
!
      do ilake = 1 , ixm1
        do n = 1 , nnsg
          if ( veg2d1(n,ilake,jslc).eq.14. ) then
            dayl = (nnnnnn-nstrt0)/4. + (xtime+dtmin)/1440.
            tl = ts1d(n,ilake)
            vl = dsqrt(us1d(ilake)**2+vs1d(ilake)**2)
            zl = z1d(n,ilake)
            ql = qs1d(n,ilake)
            fsw = fsw1d(ilake)
            flw = -1.*flw1d(ilake)
            prec = prca2d(ilake,jslc)      !  units of prec = mm
            hsen = -1.*sent1d(n,ilake)
            hlat = -1.*evpr1d(n,ilake)
            call lake(iutlak,dayl,dtlake,tl,vl,zl,ql,fsw,flw,hsen,hlat, &
                    & tgl,evl,prec,aveice,snow)
 
            tg1d(n,ilake) = tgl
            tgb1d(n,ilake) = tgl
            if ( aveice.le.10. ) then
              ldoc1d(n,ilake) = 0.
              sice1d(n,ilake) = 0.
              scv1d(n,ilake) = 0.
              sag1d(n,ilake) = 0.
            else
              ldoc1d(n,ilake) = 2.
              sice1d(n,ilake) = aveice  !  units of ice = mm
              scv1d(n,ilake) = snow   !  units of snow = mm h2o
            end if
 
          end if
        end do
      end do
 
      end subroutine lakedrv
