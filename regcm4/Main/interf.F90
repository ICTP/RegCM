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
 
      subroutine interf(innex,j)

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  this subroutine interfaces mm42d and bats variables
!
!  innex = 1 ,   mm42d --> bats
!  innex = 2 ,   bats --> mm42d
!
      use mod_regcm_param
      use mod_param1 , only : dtbat , dtmin
      use mod_param2
      use mod_param3 , only : ptop
      use mod_main
      use mod_pbldim
      use mod_slice
      use mod_bats
      use mod_constants , only : tau1 , zlnd , zoce , zsno , rgti ,     &
                 & rgas , tmelt
      use mod_date , only : jyear , jyear0 , jyearr , ntime , ktau ,    &
                 & ktaur
      implicit none
!
! Dummy arguments
!
      integer :: innex , j
      intent (in) innex , j
!
! Local variables
!
      real(8) :: amxtem , facb , facs , fact , factuv , facv , fracb ,  &
               & fracs , fracv , hl , mmpd , rh0 , satvp , sfac ,       &
               & solvt , wpm2
      integer :: i , n , nnn
      real(4) :: real_4
!
!     ******    fbat contains nummap+1 fields to be written out to
!     iutbat ******    fields are written out in following order:
!     1.  anemon west  wind (41)
!     2.  anemom south wind (42)
!     3.  drag- surface stress, in si (48)
!     4.  ground temp (4)
!     5.  temp of foliage (7)
!     6.  anemom temp (5)
!     7.  anemom spec. humidity (15)
!     8.  upper layer soil water (26)
!     9.  root zone soil water (25)
!     10.  accum precip (19)
!     11.  accum evap (39)
!     12.  accum surf runoff (12)
!     13.  snow depth in mm h2o (30)
!     14.  accum sensible heat (40)
!     15.  accum net ir (38)
!     16.  accum net solar abs (37)
!     17.  accum downward ir
!     18.  accum solar incident at surface
!     19.  convective precipitation
!     20.  surface pressure
!     21.  pbl height
 
      if ( innex.eq.1 ) then
        do i = np1 , nbmax
          do n = 1 , nnsg
            p1d0(n,i) = (psb(i,j)+ptop)*1000.
            z1d(n,i) = za(i,kx,j)
            ts1d0(n,i) = thx3d(i,kx,j)
            qs1d0(n,i) = qvb3d(i,kx,j)/(1.+qvb3d(i,kx,j))
            qs1d(n,i) = qs1d0(n,i)
 
            hl = 597.3 - .566*(ts1d0(n,i)-tmelt)
            satvp = 6.11*dexp(9.045*hl*(1./tmelt-1./ts1d0(n,i)))
            rh0 = dmax1(qs1d0(n,i)/(.622*satvp/(p1d0(n,i)*0.01-satvp)), &
                & 0.D0)
 
            ts1d(n,i) = ts1d0(n,i) - 6.5E-3*rgti*(ht1(n,i,j)-ht(i,j))
            p1d(n,i) = p1d0(n,i)*(ts1d(n,i)/ts1d0(n,i))
 
            hl = 597.3 - .566*(ts1d(n,i)-tmelt)
            satvp = 6.11*dexp(9.045*hl*(1./tmelt-1./ts1d(n,i)))
            qs1d(n,i) = dmax1(rh0*.622*satvp/(p1d(n,i)*0.01-satvp),0.D0)
 
            tg1d(n,i) = tg2d(n,i,j)
            rhs1d(n,i) = p1d(n,i)/(rgas*ts1d(n,i))
            prcp1d(n,i) = pptnc(i,j) + pptc(i,j)
!
!           quantities stored on 2d surface array for bats use only
!
            tgb1d(n,i) = tgb2d(n,i,j)
            taf1d(n,i) = taf2d(n,i,j)
            tlef1d(n,i) = tlef2d(n,i,j)
            tsw1d(n,i) = swt2d(n,i,j)
            rsw1d(n,i) = srw2d(n,i,j)
            ssw1d(n,i) = ssw2d(n,i,j)
            ldew1d(n,i) = dew2d(n,i,j)
            sag1d(n,i) = sag2d(n,i,j)
            scv1d(n,i) = scv2d(n,i,j)
            sice1d(n,i) = sice2d(n,i,j)
            gwet1d(n,i) = gwet2d(n,i,j)
            sent1d(n,i) = hfx(i,j)
            evpr1d(n,i) = qfx(i,j)
            ldoc1d(n,i) = ocld2d(n,i,j)
            ircp1d(n,i) = ircp2d(n,i,j)
            lveg(n,i) = nint(veg2d1(n,i,j))
            amxtem = dmax1(298.-tgb1d(n,i),0.D0)
            sfac = 1. - dmax1(0.D0,1.-0.0016*amxtem**2)
            if ( lveg(n,i).eq.0 ) then
              veg1d(n,i) = 0.
            else
              veg1d(n,i) = vegc(lveg(n,i)) - seasf(lveg(n,i))*sfac
            end if
            emiss_1d(n,i) = emiss2d(n,i,j)
          end do
 
          rh0 = 0.0D0
          do n = 1 , nnsg
            rh0 = rh0 + (qs1d(n,i)-qs1d0(n,i))
          end do
          rh0 = rh0/nnsg
          do n = 1 , nnsg
            qs1d(n,i) = dmax1(qs1d(n,i)-rh0,0.0D0)
          end do
 
          us1d(i) = ubx3d(i,kx,j)
          vs1d(i) = vbx3d(i,kx,j)
          fsw1d(i) = fsw2d(i,j)
          flw1d(i) = flw2d(i,j)
          solis(i) = sol2d(i,j)
          sabveg(i) = sabv2d(i,j)
          solvt = solvd2d(i,j) + solvs2d(i,j)
          if ( solvt.gt.0.0 ) then
            fracd(i) = solvd2d(i,j)/solvt
          else
            fracd(i) = 0.2
          end if
          czen(i) = dmax1(coszrs(i),0.D0)
        end do
        npts = nbmax
 
      else if ( innex.eq.2 ) then
 
        do i = np1 , nbmax
          uvdrag(i,j) = 0.0
          hfx(i,j) = 0.0
          qfx(i,j) = 0.0
          tgb(i,j) = 0.0
          tga(i,j) = 0.0
          tgbb(i,j) = 0.0
!chem2
          ssw2da(i,j) = 0.0
          sdeltk2d(i,j) = 0.0
          sdelqk2d(i,j) = 0.0
          sfracv2d(i,j) = 0.0
          sfracb2d(i,j) = 0.0
          sfracs2d(i,j) = 0.0
          svegfrac2d(i,j) = 0.0
!chem2_
          do n = 1 , nnsg
            uvdrag(i,j) = uvdrag(i,j) + drag1d(n,i)
            hfx(i,j) = hfx(i,j) + sent1d(n,i)
            qfx(i,j) = qfx(i,j) + evpr1d(n,i)
            tgb(i,j) = tgb(i,j) + tg1d(n,i)
            tga(i,j) = tga(i,j) + tg1d(n,i)
!chem2
            ssw2da(i,j) = ssw2da(i,j) + ssw1d(n,i)
            sdeltk2d(i,j) = sdeltk2d(i,j) + delt1d(n,i)
            sdelqk2d(i,j) = sdelqk2d(i,j) + delq1d(n,i)
            sfracv2d(i,j) = sfracv2d(i,j) + sigf(n,i)
            sfracb2d(i,j) = sfracb2d(i,j) + (1.-veg1d(n,i))             &
                          & *(1.-scvk(n,i))
            sfracs2d(i,j) = sfracs2d(i,j) + veg1d(n,i)*wt(n,i)          &
                          & + (1.-veg1d(n,i))*scvk(n,i)
            svegfrac2d(i,j) = svegfrac2d(i,j) + veg1d(n,i)
!chem2_
 
            if ( iocnflx.eq.1 .or.                                      &
               & (iocnflx.eq.2 .and. ocld2d(n,i,j).ge.0.5) ) then
              tgbb(i,j) = tgbb(i,j)                                     &
                        & + ((1.-veg1d(n,i))*tg1d(n,i)**4+veg1d(n,i)    &
                        & *tlef1d(n,i)**4)**0.25
            else
              tgbb(i,j) = tgbb(i,j) + tg1d(n,i)
            end if
            if ( ocld2d(n,i,j).lt.0.5 ) then
              ssw1d(n,i) = -1.E34
              rsw1d(n,i) = -1.E34
              tsw1d(n,i) = -1.E34
              rno1d(n,i) = -1.E34
              rnos1d(n,i) = -1.E34
              scv1d(n,i) = -1.E34
            end if
          end do
          uvdrag(i,j) = uvdrag(i,j)/float(nnsg)
          hfx(i,j) = hfx(i,j)/float(nnsg)
          qfx(i,j) = qfx(i,j)/float(nnsg)
          tgb(i,j) = tgb(i,j)/float(nnsg)
          tga(i,j) = tga(i,j)/float(nnsg)
          tgbb(i,j) = tgbb(i,j)/float(nnsg)
!chem2
          ssw2da(i,j) = ssw2da(i,j)/float(nnsg)
          sdeltk2d(i,j) = sdeltk2d(i,j)/float(nnsg)
          sdelqk2d(i,j) = sdelqk2d(i,j)/float(nnsg)
          sfracv2d(i,j) = sfracv2d(i,j)/float(nnsg)
          sfracb2d(i,j) = sfracb2d(i,j)/float(nnsg)
          sfracs2d(i,j) = sfracs2d(i,j)/float(nnsg)
          svegfrac2d(i,j) = svegfrac2d(i,j)/float(nnsg)
!chem2_
          do n = 1 , nnsg
            snowc(n,i,j) = scv1d(n,i)
            tg2d(n,i,j) = tg1d(n,i)
            tgb2d(n,i,j) = tgb1d(n,i)
            taf2d(n,i,j) = taf1d(n,i)
            tlef2d(n,i,j) = tlef1d(n,i)
            swt2d(n,i,j) = tsw1d(n,i)
            srw2d(n,i,j) = rsw1d(n,i)
            ssw2d(n,i,j) = ssw1d(n,i)
            dew2d(n,i,j) = ldew1d(n,i)
            sag2d(n,i,j) = sag1d(n,i)
            scv2d(n,i,j) = scv1d(n,i)
            sice2d(n,i,j) = sice1d(n,i)
            gwet2d(n,i,j) = gwet1d(n,i)
            ocld2d(n,i,j) = ldoc1d(n,i)
            ircp2d(n,i,j) = ircp1d(n,i)
            evpa2d(n,i,j) = evpa2d(n,i,j) + dtbat*evpr1d(n,i)
            sena2d(n,i,j) = sena2d(n,i,j) + dtbat*sent1d(n,i)
            if ( rnos2d(n,i,j).gt.-1.E10 .and. rnos1d(n,i).gt.-1.E10 )  &
               & then
              rnos2d(n,i,j) = rnos2d(n,i,j) + rnos1d(n,i)/tau1*dtbat
            else
              rnos2d(n,i,j) = -1.E34
            end if
            if ( rno2d(n,i,j).gt.-1.E10 .and. rnos1d(n,i)               &
               & .gt.-1.E10 .and. rno1d(n,i).gt.-1.E10 ) then
              rno2d(n,i,j) = rno2d(n,i,j) + (rno1d(n,i)-rnos1d(n,i))    &
                           & /tau1*dtbat
            else
              rno2d(n,i,j) = -1.E34
            end if
          end do
 
!
!         quantities stored on 2d surface array for bats use only
!
          prca2d(i,j) = prca2d(i,j) + dtbat*pptc(i,j)
          prnca2d(i,j) = prnca2d(i,j) + dtbat*pptnc(i,j)
          flwa2d(i,j) = flwa2d(i,j) + dtbat*flw1d(i)
          flwda2d(i,j) = flwda2d(i,j) + dtbat*flwd2d(i,j)
          fswa2d(i,j) = fswa2d(i,j) + dtbat*fsw1d(i)
          svga2d(i,j) = svga2d(i,j) + dtbat*sabveg(i)
          sina2d(i,j) = sina2d(i,j) + dtbat*sinc2d(i,j)
          pptnc(i,j) = 0.
          pptc(i,j) = 0.
        end do
        do i = np1 , nbmax
#ifdef MPP1
          u10m_o(j,i-1) = 0.0
          v10m_o(j,i-1) = 0.0
          tg_o(j,i-1) = 0.0
          t2m_o(j,i-1) = 0.0
          do n = 1 , nnsg
            if ( ocld2d(n,i,j).ge.0.5 ) then
              fracv = sigf(n,i)
              fracb = (1.-veg1d(n,i))*(1.-scvk(n,i))
              fracs = veg1d(n,i)*wt(n,i) + (1.-veg1d(n,i))*scvk(n,i)
              facv = dlog(z1(n,i)/2.)/dlog(z1(n,i)/rough(lveg(n,i)))
              facb = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zlnd)
              facs = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zsno)
              fact = fracv*facv + fracb*facb + fracs*facs
              facv = dlog(z1(n,i)/10.)/dlog(z1(n,i)/rough(lveg(n,i)))
              facb = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zlnd)
              facs = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zsno)
              factuv = fracv*facv + fracb*facb + fracs*facs
              u10m1d(n,i) = us1d(i)*(1.-factuv)
              v10m1d(n,i) = vs1d(i)*(1.-factuv)
              t2m_1d(n,i) = ts1d(n,i) - delt1d(n,i)*fact
            else if ( iocnflx.eq.1 ) then
              fact = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zoce)
              factuv = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zoce)
              u10m1d(n,i) = us1d(i)*(1.-factuv)
              v10m1d(n,i) = vs1d(i)*(1.-factuv)
              t2m_1d(n,i) = ts1d(n,i) - delt1d(n,i)*fact
            else
            end if
            tg_s(n,j,i-1) = tg1d(n,i)
            u10m_s(n,j,i-1) = u10m1d(n,i)
            v10m_s(n,j,i-1) = v10m1d(n,i)
            t2m_s(n,j,i-1) = t2m_1d(n,i)
 
            u10m_o(j,i-1) = u10m_o(j,i-1) + u10m1d(n,i)
            v10m_o(j,i-1) = v10m_o(j,i-1) + v10m1d(n,i)
            t2m_o(j,i-1) = t2m_o(j,i-1) + t2m_1d(n,i)
            tg_o(j,i-1) = tg_o(j,i-1) + tg1d(n,i)
          end do
          u10m_o(j,i-1) = u10m_o(j,i-1)/float(nnsg)
          v10m_o(j,i-1) = v10m_o(j,i-1)/float(nnsg)
          t2m_o(j,i-1) = t2m_o(j,i-1)/float(nnsg)
          tg_o(j,i-1) = tg_o(j,i-1)/float(nnsg)
 
          tgmx_o(j,i-1) = amax1(tgmx_o(j,i-1),tg_o(j,i-1))
          tgmn_o(j,i-1) = amin1(tgmn_o(j,i-1),tg_o(j,i-1))
          t2mx_o(j,i-1) = amax1(t2mx_o(j,i-1),t2m_o(j,i-1))
          t2mn_o(j,i-1) = amin1(t2mn_o(j,i-1),t2m_o(j,i-1))
          w10x_o(j,i-1) = amax1(w10x_o(j,i-1),sqrt(u10m_o(j,i-1)**2+    &
                        & v10m_o(j,i-1)**2))
          real_4 = (psb(i,j)+ptop)*10.
          psmn_o(j,i-1) = amin1(psmn_o(j,i-1),real_4)
#else
          u10m_o(j-1,i-1) = 0.0
          v10m_o(j-1,i-1) = 0.0
          tg_o(j-1,i-1) = 0.0
          t2m_o(j-1,i-1) = 0.0
          do n = 1 , nnsg
            if ( ocld2d(n,i,j).ge.0.5 ) then
              fracv = sigf(n,i)
              fracb = (1.-veg1d(n,i))*(1.-scvk(n,i))
              fracs = veg1d(n,i)*wt(n,i) + (1.-veg1d(n,i))*scvk(n,i)
              facv = dlog(z1(n,i)/2.)/dlog(z1(n,i)/rough(lveg(n,i)))
              facb = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zlnd)
              facs = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zsno)
              fact = fracv*facv + fracb*facb + fracs*facs
              facv = dlog(z1(n,i)/10.)/dlog(z1(n,i)/rough(lveg(n,i)))
              facb = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zlnd)
              facs = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zsno)
              factuv = fracv*facv + fracb*facb + fracs*facs
              u10m1d(n,i) = us1d(i)*(1.-factuv)
              v10m1d(n,i) = vs1d(i)*(1.-factuv)
              t2m_1d(n,i) = ts1d(n,i) - delt1d(n,i)*fact
            else if ( iocnflx.eq.1 ) then
              fact = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zoce)
              factuv = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zoce)
              u10m1d(n,i) = us1d(i)*(1.-factuv)
              v10m1d(n,i) = vs1d(i)*(1.-factuv)
              t2m_1d(n,i) = ts1d(n,i) - delt1d(n,i)*fact
            else
            end if
            tg_s(n,j-1,i-1) = tg1d(n,i)
            u10m_s(n,j-1,i-1) = u10m1d(n,i)
            v10m_s(n,j-1,i-1) = v10m1d(n,i)
            t2m_s(n,j-1,i-1) = t2m_1d(n,i)

            u10m_o(j-1,i-1) = u10m_o(j-1,i-1) + u10m1d(n,i)
            v10m_o(j-1,i-1) = v10m_o(j-1,i-1) + v10m1d(n,i)
            t2m_o(j-1,i-1) = t2m_o(j-1,i-1) + t2m_1d(n,i)
            tg_o(j-1,i-1) = tg_o(j-1,i-1) + tg1d(n,i)
          end do
          u10m_o(j-1,i-1) = u10m_o(j-1,i-1)/float(nnsg)
          v10m_o(j-1,i-1) = v10m_o(j-1,i-1)/float(nnsg)
          t2m_o(j-1,i-1) = t2m_o(j-1,i-1)/float(nnsg)
          tg_o(j-1,i-1) = tg_o(j-1,i-1)/float(nnsg)
          tgmx_o(j-1,i-1) = amax1(tgmx_o(j-1,i-1),tg_o(j-1,i-1))
          tgmn_o(j-1,i-1) = amin1(tgmn_o(j-1,i-1),tg_o(j-1,i-1))
          t2mx_o(j-1,i-1) = amax1(t2mx_o(j-1,i-1),t2m_o(j-1,i-1))
          t2mn_o(j-1,i-1) = amin1(t2mn_o(j-1,i-1),t2m_o(j-1,i-1))
          w10x_o(j-1,i-1) = amax1(w10x_o(j-1,i-1),sqrt(u10m_o(j-1,i-1)**&
                          & 2+v10m_o(j-1,i-1)**2))
          real_4 = (psb(i,j)+ptop)*10.
          psmn_o(j-1,i-1) = amin1(psmn_o(j-1,i-1),real_4)
#endif
        end do

        if ( mod(ntime+nint(dtmin*60.),kbats).eq.0 .or.                 &
           & (jyear.eq.jyearr .and. ktau.eq.ktaur) ) then
          if ( jyear.eq.jyear0 .and. ktau.le.1 ) then
            mmpd = 86400./dtbat
            wpm2 = 1./dtbat
          else if ( jyear.eq.jyear0 .and. dble(ktau*dtmin)              &
                  & .le.batfrq*60.+0.01 ) then
            mmpd = 24./(batfrq-dtmin/60.)
            wpm2 = 1./((batfrq-dtmin/60.)*3600.)
          else
            mmpd = 24./batfrq
            wpm2 = 1./(batfrq*3600.)
          end if
          do i = np1 , nbmax
#ifdef MPP1
            drag_o(j,i-1) = 0.0
            q2m_o(j,i-1) = 0.0
            evpa_o(j,i-1) = 0.0
            sena_o(j,i-1) = 0.0
            do n = 1 , nnsg
              if ( ocld2d(n,i,j).ge.0.5 ) then
                fracv = sigf(n,i)
                fracb = (1.-veg1d(n,i))*(1.-scvk(n,i))
                fracs = veg1d(n,i)*wt(n,i) + (1.-veg1d(n,i))*scvk(n,i)
                facv = dlog(z1(n,i)/2.)/dlog(z1(n,i)/rough(lveg(n,i)))
                facb = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zlnd)
                facs = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zsno)
                fact = fracv*facv + fracb*facb + fracs*facs
                q2m_1d(n,i) = qs1d(n,i) - delq1d(n,i)*fact
              else if ( iocnflx.eq.1 ) then
                fact = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zoce)
                q2m_1d(n,i) = qs1d(n,i) - delq1d(n,i)*fact
              else
              end if
              q2m_s(n,j,i-1) = q2m_1d(n,i)
              drag_s(n,j,i-1) = drag1d(n,i)
              evpa_s(n,j,i-1) = evpa2d(n,i,j)*mmpd
              sena_s(n,j,i-1) = sena2d(n,i,j)*wpm2
              tpr_s(n,j,i-1) = (prnca2d(i,j)+prca2d(i,j))*mmpd
              prcv_s(n,j,i-1) = prca2d(i,j)*mmpd
              ps_s(n,j,i-1) = p1d(n,i)*0.01
 
              q2m_o(j,i-1) = q2m_o(j,i-1) + q2m_1d(n,i)
              drag_o(j,i-1) = drag_o(j,i-1) + drag1d(n,i)
              evpa_o(j,i-1) = evpa_o(j,i-1) + evpa2d(n,i,j)
              sena_o(j,i-1) = sena_o(j,i-1) + sena2d(n,i,j)
            end do
            tpr_o(j,i-1) = (prnca2d(i,j)+prca2d(i,j))*mmpd
            q2m_o(j,i-1) = q2m_o(j,i-1)/float(nnsg)
            drag_o(j,i-1) = drag_o(j,i-1)/float(nnsg)
            evpa_o(j,i-1) = evpa_o(j,i-1)/float(nnsg)*mmpd
            sena_o(j,i-1) = sena_o(j,i-1)/float(nnsg)*wpm2
            flwa_o(j,i-1) = flwa2d(i,j)*wpm2
            fswa_o(j,i-1) = fswa2d(i,j)*wpm2
            flwd_o(j,i-1) = flwda2d(i,j)*wpm2
            sina_o(j,i-1) = sina2d(i,j)*wpm2
            prcv_o(j,i-1) = prca2d(i,j)*mmpd
            ps_o(j,i-1) = (psb(i,j)+ptop)*10.
            zpbl_o(j,i-1) = zpbl(i,j)
 
            tlef_o(j,i-1) = 0.0
            ssw_o(j,i-1) = 0.0
            rsw_o(j,i-1) = 0.0
            rnos_o(j,i-1) = 0.0
            scv_o(j,i-1) = 0.0
            nnn = 0
            do n = 1 , nnsg
              if ( ocld2d(n,i,j).ge.0.5 ) then
                tlef_o(j,i-1) = tlef_o(j,i-1) + tlef1d(n,i)
                ssw_o(j,i-1) = ssw_o(j,i-1) + ssw1d(n,i)
                rsw_o(j,i-1) = rsw_o(j,i-1) + rsw1d(n,i)
                rnos_o(j,i-1) = rnos_o(j,i-1) + rnos2d(n,i,j)
                scv_o(j,i-1) = scv_o(j,i-1) + scv1d(n,i)
                tlef_s(n,j,i-1) = tlef1d(n,i)
                ssw_s(n,j,i-1) = ssw1d(n,i)
                rsw_s(n,j,i-1) = rsw1d(n,i)
                rnos_s(n,j,i-1) = rnos2d(n,i,j)*mmpd
                scv_s(n,j,i-1) = scv1d(n,i)
                nnn = nnn + 1
              else
                tlef_s(n,j,i-1) = -1.E34
                ssw_s(n,j,i-1) = -1.E34
                rsw_s(n,j,i-1) = -1.E34
                rnos_s(n,j,i-1) = -1.E34
                scv_s(n,j,i-1) = -1.E34
              end if
            end do
            if ( nnn.ge.max0(nnsg/2,1) ) then
              tlef_o(j,i-1) = tlef_o(j,i-1)/float(nnn)
              ssw_o(j,i-1) = ssw_o(j,i-1)/float(nnn)
              rsw_o(j,i-1) = rsw_o(j,i-1)/float(nnn)
              rnos_o(j,i-1) = rnos_o(j,i-1)/float(nnn)*mmpd
              scv_o(j,i-1) = scv_o(j,i-1)/float(nnn)
            else
              tlef_o(j,i-1) = -1.E34
              ssw_o(j,i-1) = -1.E34
              rsw_o(j,i-1) = -1.E34
              rnos_o(j,i-1) = -1.E34
              scv_o(j,i-1) = -1.E34
            end if
#else
            drag_o(j-1,i-1) = 0.0
            q2m_o(j-1,i-1) = 0.0
            evpa_o(j-1,i-1) = 0.0
            sena_o(j-1,i-1) = 0.0
            do n = 1 , nnsg
              if ( ocld2d(n,i,j).ge.0.5 ) then
                fracv = sigf(n,i)
                fracb = (1.-veg1d(n,i))*(1.-scvk(n,i))
                fracs = veg1d(n,i)*wt(n,i) + (1.-veg1d(n,i))*scvk(n,i)
                facv = dlog(z1(n,i)/2.)/dlog(z1(n,i)/rough(lveg(n,i)))
                facb = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zlnd)
                facs = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zsno)
                fact = fracv*facv + fracb*facb + fracs*facs
                q2m_1d(n,i) = qs1d(n,i) - delq1d(n,i)*fact
              else if ( iocnflx.eq.1 ) then
                fact = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zoce)
                q2m_1d(n,i) = qs1d(n,i) - delq1d(n,i)*fact
              else
              end if
              q2m_s(n,j-1,i-1) = q2m_1d(n,i)
              drag_s(n,j-1,i-1) = drag1d(n,i)
              evpa_s(n,j-1,i-1) = evpa2d(n,i,j)*mmpd
              sena_s(n,j-1,i-1) = sena2d(n,i,j)*wpm2
              tpr_s(n,j-1,i-1) = (prnca2d(i,j)+prca2d(i,j))*mmpd
              prcv_s(n,j-1,i-1) = prca2d(i,j)*mmpd
              ps_s(n,j-1,i-1) = p1d(n,i)*0.01

              q2m_o(j-1,i-1) = q2m_o(j-1,i-1) + q2m_1d(n,i)
              drag_o(j-1,i-1) = drag_o(j-1,i-1) + drag1d(n,i)
              evpa_o(j-1,i-1) = evpa_o(j-1,i-1) + evpa2d(n,i,j)
              sena_o(j-1,i-1) = sena_o(j-1,i-1) + sena2d(n,i,j)
            end do
            tpr_o(j-1,i-1) = (prnca2d(i,j)+prca2d(i,j))*mmpd
            q2m_o(j-1,i-1) = q2m_o(j-1,i-1)/float(nnsg)
            drag_o(j-1,i-1) = drag_o(j-1,i-1)/float(nnsg)
            evpa_o(j-1,i-1) = evpa_o(j-1,i-1)/float(nnsg)*mmpd
            sena_o(j-1,i-1) = sena_o(j-1,i-1)/float(nnsg)*wpm2
            flwa_o(j-1,i-1) = flwa2d(i,j)*wpm2
            fswa_o(j-1,i-1) = fswa2d(i,j)*wpm2
            flwd_o(j-1,i-1) = flwda2d(i,j)*wpm2
            sina_o(j-1,i-1) = sina2d(i,j)*wpm2
            prcv_o(j-1,i-1) = prca2d(i,j)*mmpd
            ps_o(j-1,i-1) = (psb(i,j)+ptop)*10.
            zpbl_o(j-1,i-1) = zpbl(i,j)

            tlef_o(j-1,i-1) = 0.0
            ssw_o(j-1,i-1) = 0.0
            rsw_o(j-1,i-1) = 0.0
            rnos_o(j-1,i-1) = 0.0
            scv_o(j-1,i-1) = 0.0
            nnn = 0
            do n = 1 , nnsg
              if ( ocld2d(n,i,j).ge.0.5 ) then
                tlef_o(j-1,i-1) = tlef_o(j-1,i-1) + tlef1d(n,i)
                ssw_o(j-1,i-1) = ssw_o(j-1,i-1) + ssw1d(n,i)
                rsw_o(j-1,i-1) = rsw_o(j-1,i-1) + rsw1d(n,i)
                rnos_o(j-1,i-1) = rnos_o(j-1,i-1) + rnos2d(n,i,j)
                scv_o(j-1,i-1) = scv_o(j-1,i-1) + scv1d(n,i)
                tlef_s(n,j-1,i-1) = tlef1d(n,i)
                ssw_s(n,j-1,i-1) = ssw1d(n,i)
                rsw_s(n,j-1,i-1) = rsw1d(n,i)
                rnos_s(n,j-1,i-1) = rnos2d(n,i,j)*mmpd
                scv_s(n,j-1,i-1) = scv1d(n,i)
                nnn = nnn + 1
              else
                tlef_s(n,j-1,i-1) = -1.E34
                ssw_s(n,j-1,i-1) = -1.E34
                rsw_s(n,j-1,i-1) = -1.E34
                rnos_s(n,j-1,i-1) = -1.E34
                scv_s(n,j-1,i-1) = -1.E34
              end if
            end do
            if ( nnn.ge.max0(nnsg/2,1) ) then
              tlef_o(j-1,i-1) = tlef_o(j-1,i-1)/float(nnn)
              ssw_o(j-1,i-1) = ssw_o(j-1,i-1)/float(nnn)
              rsw_o(j-1,i-1) = rsw_o(j-1,i-1)/float(nnn)
              rnos_o(j-1,i-1) = rnos_o(j-1,i-1)/float(nnn)*mmpd
              scv_o(j-1,i-1) = scv_o(j-1,i-1)/float(nnn)
            else
              tlef_o(j-1,i-1) = -1.E34
              ssw_o(j-1,i-1) = -1.E34
              rsw_o(j-1,i-1) = -1.E34
              rnos_o(j-1,i-1) = -1.E34
              scv_o(j-1,i-1) = -1.E34
            end if
#endif
 
!           ******    reset accumulation arrays to zero
            do n = 1 , nnsg
              evpa2d(n,i,j) = 0.
              rnos2d(n,i,j) = 0.
              sena2d(n,i,j) = 0.
            end do
            prnca2d(i,j) = 0.
            prca2d(i,j) = 0.
            flwa2d(i,j) = 0.
            flwda2d(i,j) = 0.
            fswa2d(i,j) = 0.
            svga2d(i,j) = 0.
            sina2d(i,j) = 0.
          end do
        end if
!
      else                 !end innex test
      end if
      end subroutine interf
