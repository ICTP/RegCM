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
 
  module mod_che_tend
!
! Tendency and budget for tracer transport and chemicals
!
  use mod_constants
  use mod_dynparam
  use mod_che_common
  use mod_che_indices
  use mod_che_param
  use mod_che_sox
  use mod_che_drydep
  use mod_che_wetdep
  use mod_che_emission
! use mod_chem_sox
! use mod_sea_salt
! use mod_chem_emis
  use mod_che_dust
  use mod_che_seasalt
  use mod_che_carbonaer
  use mod_che_mppio
  private

  public :: tractend2 , tracbud

  contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     This subroutine computes the tendencies for tracer transport and
!     chemistry
!
!     j:             index of j slice in current computation
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine tractend2(jstart,jend,istart,iend,ktau,lmonth,xkc)
      implicit none
!
      integer , intent(in) :: jstart , jend , istart , iend , lmonth
      integer(8) , intent(in) :: ktau
      real(8) , pointer , dimension(:,:,:) , intent(in):: xkc
!
      real(8) :: agct , ak00t , ak0tm , akval , clmin , facb , facs , &
                 fact , facv , pres10 , qsat10 , remcum , satvp ,     &
                 shu10 , u10 , v10 , chias , chibs

      real(8) , dimension(ntr) :: agingtend
      real(8) , dimension(iy,kz) :: wk, rho , settend , &
                                    ttb, wl, fracloud, fracum , prec

      integer :: i , j , ibin , itr , k , kk , kb , kdwd
      integer , dimension(iy) :: ivegcov

      real(8) , dimension(iy,kz,ntr) :: pdepv
      real(8) , dimension(iy,ntr) :: ddepa


      real(8) , dimension(iy) :: psurf , rh10 , soilw , srad ,  &
          temp10 , tsurf , vegfrac , wid10 , zeff , ustar

      real(8) , dimension(iy,nbin) :: rsfrow
      real(8), dimension(iy,sbin) :: seasalt_flx
      real(8), dimension(iy,ntr) :: drydepvg
      real(8) , dimension(ntr) :: wetrem , wetrem_cvc
!
      integer :: igaschem !!!PROVISOIRE
!
!**************************************************************************
!     A : PRELIMINARY CALCULATIONS
!*************************************************************************
! 
      do j = jstart , jend
        rho = d_zero
        wl = d_zero
        ttb = d_zero
        prec = d_zero
        fracloud = d_zero
        fracum = d_zero
        psurf = d_zero
        igaschem = 1

!       the unit: rho - kg/m3, wl - g/m3
        do k = 1 , kz
          do i = istart , iend
!           rho(i,k) = (sps2%ps(i,j)*a(k)+r8pt)* &
!      what the hell   1000./287./atm2%t(i,k,j)*sps2%ps(i,j)
            rho(i,k) = crhob3d(i,k,j)
            wl(i,k) = cqcb3d(i,k,j)*rho(i,k)
            ttb(i,k) = ctb3d(i,k,j)
!           precipiation rate is a rquired variable for deposition routines.
!           It is directly taken as rembc (saved in precip routine) in mm/hr !!
            prec(i,k) = crembc(i,k) / 3600.D0 !passed in mm/s  
          end do
        end do
!       cloud fractionnal cover for wet deposition
!       large scale : fracloud, calculated from fcc coming from pcp.f
!       cumulus scale : fracum, calculated from the total cloud fraction
!       (as defined for the radiation scheme in cldfrac.f routine)
        do i = istart , iend
          do k = 1 , kz
            fracloud(i,k) = cfcc(i,k,j)
            fracum(i,k) = d_zero
          end do
          if ( kcumtop(j,i) > 0 ) then
            do kk = kcumtop(j,i) , kz
              fracum(i,kk) = ccldfra(j,i,kk) - fracloud(i,kk)
            end do
          end if
        end do
!
!       variables used for natural fluxes and deposition velocities 
! 
        ivegcov=0   
        do i = istart , iend
          ivegcov(i) = cveg2d(i,j)
          psurf(i) = cpsb(i,j) * 1.0D3 + ptop
 
!         method based on bats diagnostic in routine interf.
 
          if ( (ivegcov(i) /= 0) ) then
            facv = dlog(cza(i,kz,j)/d_10) / &
                   dlog(cza(i,kz,j)/crough(ivegcov(i)))
            facb = dlog(cza(i,kz,j)/d_10)/dlog(cza(i,kz,j)/zlnd)
            facs = dlog(cza(i,kz,j)/d_10)/dlog(cza(i,kz,j)/zsno)
 
!           fact = csfracv2d(i,j)*facv 
            fact = cvegfrac(i,j) * facv + (d_one-cvegfrac(i,j)) * facb
!           FAB REVOIR CETTE partie et definir interface pour sfracs,sfracv
!                   + sfracb2d(i,j)*facb + sfracs2d(i,j)*facs
! 
!           grid level effective roughness lenght (linear averaging for now)
!           zeff(i) = rough(ivegcov(i))*sfracv2d(i,j) + &
!                     zlnd * sfracb2d(i,j) + zsno * sfracs2d(i,j)
            zeff(i) = crough(ivegcov(i))*cvegfrac(i,j) + &
                      zlnd*(d_one-cvegfrac(i,j))
          else
!           water surface
            fact = dlog(cza(i,kz,j)/d_10)/dlog(cza(i,kz,j)/zoce)
            zeff(i) = zoce
          end if
!         10 m wind
          u10 = (cubx3d(i,kz,j))*(1-fact)
          v10 = (cvbx3d(i,kz,j))*(1-fact)
          wid10(i) = sqrt(u10**2+v10**2)
!         10 m air temperature
          temp10(i) = ttb(i,kz) - csdeltk2d(i,j)*fact
!         specific  humidity at 10m
          shu10 = cqvb3d(i,kz,j)/ &
                  (d_one+cqvb3d(i,kz,j))-csdelqk2d(i,j)*fact
!         back to mixing ratio
          shu10 = shu10/(1-shu10)
!         saturation mixing ratio at 10m
          if ( temp10(i) > tzero ) then
            satvp = svp1*1.0D3*dexp(svp2*(temp10(i)-tzero)/(temp10(i)-svp3))
          else
            satvp = svp4*1.0D3*dexp(svp5-svp6/temp10(i))
          end if
          pres10 = psurf(i) - 98.0D0
          qsat10 = ep2*satvp/(pres10-satvp)
!         relative humidity at 10m
          rh10(i) = d_zero
          if ( qsat10 > d_zero ) rh10(i) = shu10/qsat10
!
!         friction velocity ( from uvdrag so updtaed at  bats or clm frequency )
!
!!$   FAB AFIXER         ustar(i) = sqrt ( sfsta%uvdrag(i,j)             *           &
!!$           &           sqrt ( (atm2%u(i,kz,j)/sps2%ps(i,j) )**2   +     &
!!$           &                  (atm2%v(i,kz,j)/sps2%ps(i,j) )**2 ) /     &
!!$           &                   rho(i,kz) )
 
!           soil wetness
 
          soilw(i) = cssw2da(i,j)/cdepuv(idnint(clndcat(i,j)))/(2650.0D0 * &
                (d_one-cxmopor(ciexsol(idnint(clndcat(i,j))))))
!         fraction of vegetation
          vegfrac(i) = cvegfrac(i,j)
!         surface temperature
!         over land recalculated from the BATS  deltk air/ surface
!         temperature account for a composite temperature between
!         bare ground and vegetation
          if ( ivegcov(i) /= 0 ) then
            tsurf(i) = ttb(i,kz) - csdeltk2d(i,j)
          else
!           ocean temperature in this case
            tsurf(i) = ctg(i,j)
          end if
 
!        aborbed solar radiation (for stb criteria used to calculate
!        aerodynamic resistance)
 
         srad(i) = csol2d(i,j)
 
        end do
!
!       END of preliminary calculations)
!
!*****************************************************************
! B :CALCULATION OF TRACER TENDENCY (except full gas phase chemistry solver)
!*****************************************************************
!
!       SOX CHEMSITRY ( from offline oxidant) 
!
        if ( igaschem == 0 ) then
!FAB :    regler le probleme de gas vs aerosol only
          if (iso2 > 0 .and. iso4 >0.) then
            call chemsox(j,wl,fracloud,fracum,rho,ttb)
          end if
        end if
!
!       aging of carboneaceous aerosols
!
        if ( (ibchb > 0 .and. ibchl > 0 ) .or. &
             (iochb > 0 .and. iochl > 0) ) then
          call aging_carb(j)
        end if

        ! NATURAL EMISSIONS FLUX and tendencies  (dust -sea salt)       

        if ( idust(1) > 0 ) then
          call sfflux(iy,2,iym2,j,ivegcov,vegfrac,ustar, &
                      zeff,soilw,wid10,rho(:,kz),dustbsiz,rsfrow)     
        end if
!       if ( isslt(1) > 0 ) call sea_salt(j,wid10,ivegcov,seasalt_flx)
!
!       update emission tendencies from inventories

        call emis_tend(ktau,j,lmonth)
!
!       aerosol settling and drydep 
!       include calculation of dry dep/settling velocities and 
!       updating tendencies
!
        pdepv = d_zero
        ddepa = d_zero
        if ( idust(1) > 0 ) then
          call drydep_aero(j,nbin,idust,rhodust,ivegcov,ttb,rho,hlev,psurf, &
                           temp10,tsurf,srad,rh10,wid10,zeff,dustbsiz,      &
                           pdepv,ddepa)
        end if

!       if (isslt(1) >0 ) then
!         call drydep_aero(j,sbin,isslt,rhosslt,ivegcov,ttb,rho,hlev,psurf, &
!                          temp10,tsurf,srad,rh10,wid10,zeff,ssltbsiz,      &
!                          pdepv,ddepa)
!       end if 

        if ( icarb(1) > 0 ) then
          ibin = count( icarb > 0 ) 
          call drydep_aero(j,ibin,icarb(1:ibin),rhooc,ivegcov,ttb,rho,hlev, &
                           psurf,temp10,tsurf,srad,rh10,wid10,zeff,         &
                           carbsiz(1:ibin,:),pdepv,ddepa)
        end if 
!!$
!       GAS phase dry deposition velocity + tendencies
!       option compatible with BATS and CLM
!!$
        if ( igaschem == 1 ) then
          call drydep_gas(j,ivegcov,rh10,srad,tsurf,prec(:,kz),temp10,  &
                          wid10,zeff,drydepvg)
        end if
!!$
!       WET deposition (rainout and washout) for aerosol
!!$
        if ( idust(1) > 0 ) then
          call wetdepa(j,nbin,idust,dustbsiz,rhodust,ttb,wl,fracloud, &
                       fracum,psurf,hlev,rho,prec,pdepv)  
        end if

!       if ( isslt(1) > 0 )  then   
!         call wetdepa(j,sbin,isslt,ssltbsiz,rhosslt,ttb,wl,fracloud, &
!                      fracum,psurf,hlev,rho, prec, pdepv )  
!       end if
!       if ( icarb(1) > 0 )  then   
!         ibin = count( icarb > 0 ) 
!         call wetdepa(j,ibin,icarb(1:ibin),carbsiz(1:ibin,:),rhobchl, &
!                      ttb,wl,fracloud,fracum,psurf,hlev,rho,prec,pdepv)  
!       end if
!
!!$
!       Wet Deposition for gasphase species 
!!$
        if ( igaschem == 1 ) then
          call sethet(j,cza(:,:,j),cht(:,j),ttb,checum,cremrat, &
                      chevap,dtche,rho,chib(:,:,j,:),iym3,cpsb(2:iym2,j))
        end if
      end do
    end subroutine tractend2
!
    subroutine conv_trans
      implicit none
!
!!$      if ( ichcumtra.eq.2 ) then
!!$        do k = 2 , kz
!!$          do i = 2 , iym2
!!$            wk(i,k) = (d_one/sps1%ps(i,j))                                 &
!!$                    & *(twt(k,1)*chib(i,k,j,itr)+twt(k,2)*chib(i,k-1,j, &
!!$                    & itr))
!!$ 
!!$            cutend_up(i,k) = 0.
!!$            cutend_dwd(i,k) = 0.
!!$          end do
!!$        end do
!!$ 
!!$        do i = 2 , iym2
!!$ 
!!$          if ( icumtop(i,j) /= 0 ) then
!!$ 
!!$            kt = max0(icumtop(i,j),3)
!!$            kb = icumbot(i,j)
!!$            kdwd = icumdwd(i,j)
!!$ 
!!$!           cutend(i,kt) =  mflx(i) * g * 1.0d-3*
!!$!           &               (wk(i,kb)-wk(i,kt))/(sigma(kb)-sigma(kt))
!!$ 
!!$!           transport linked to updraft
!!$!           betwwen kt et kdwd , the tendancy is averaged (mixing)
!!$ 
!!$            if ( kdwd < kt ) then
!!$              write (aline, *) 'Problem in tractend2 !'
!!$              call say
!!$            end if
!!$            do k = kt , kdwd
!!$              cutend_up(i,k) = mflx(i,1)*egrav*1.0D-3*wk(i,kb)             &
!!$                             & /(sigma(kdwd)-sigma(kt))
!!$            end do
!!$ 
!!$            cutend_up(i,kb) = -mflx(i,1)*egrav*1.0D-3*wk(i,kb)/(dsigma(kb))
!!$!           transport linked to downdraft
!!$ 
!!$            cutend_dwd(i,kdwd) = -mflx(i,2)*egrav*1.0D-3*wk(i,kdwd)        &
!!$                               & /(dsigma(kdwd))
!!$ 
!!$            cutend_dwd(i,kz) = +mflx(i,2)*egrav*1.0D-3*wk(i,kdwd)          &
!!$                             & /(dsigma(kz))
!!$ 
!!$            do k = kt , kz
!!$              chiten(i,k,j,itr) = chiten(i,k,j,itr) + cutend_up(i,k)    &
!!$                                & + cutend_dwd(i,k)
!!$            end do
!!$          end if
!!$        end do
!!$      end if
!!$ 
      end subroutine conv_trans

      subroutine tracbud
      implicit none
!
! Local variables
!
      integer :: i , itr , j , k
!
      do itr = 1 , ntr
        do j = 1 , jendx
          do i = 1 , iym1
            dtrace(i,j,itr) = 0.0
            wdlsc(i,j,itr) = 0.0
            wdcvc(i,j,itr) = 0.0
            wxsg(i,j,itr) = 0.0
            wxaq(i,j,itr) = 0.0
            ddsfc(i,j,itr) = 0.0
          end do
        end do
      end do
 
!-----tracers (unit = kg):
      do itr = 1 , ntr
        do j = 1 , jendx
          do i = 1 , iym1
            do k = 1 , kz
              dtrace(i,j,itr) = dtrace(i,j,itr) + chia(i,k,j,itr)       &
                              & *cdsigma(k)
              wdlsc(i,j,itr) = wdlsc(i,j,itr) + remlsc(i,k,j,itr)       &
                             & *cdsigma(k)
              wdcvc(i,j,itr) = wdcvc(i,j,itr) + remcvc(i,k,j,itr)       &
                             & *cdsigma(k)
              wxsg(i,j,itr) = wxsg(i,j,itr) + rxsg(i,k,j,itr)*cdsigma(k)
!             sum ls and conv contribution
              wxaq(i,j,itr) = wxaq(i,j,itr)                             &
                            & + (rxsaq1(i,k,j,itr)+rxsaq2(i,k,j,itr))   &
                            & *cdsigma(k)
            end do
            ddsfc(i,j,itr) = ddsfc(i,j,itr) + remdrd(i,j,itr)*cdsigma(kz)
!           Source cumulated diag(care the unit are alredy .m-2)
            cemtrac(i,j,itr) = cemtr(i,j,itr)
          end do
        end do
      end do

      do itr = 1 , ntr
        do j = 1 , jendx
          do i = 1 , iym1
            dtrace(i,j,itr) = 1.D6*dtrace(i,j,itr)*d_1000*regrav
                                                        ! unit: mg/m2
            wdlsc(i,j,itr) = 1.D6*wdlsc(i,j,itr)*d_1000*regrav
            wdcvc(i,j,itr) = 1.D6*wdcvc(i,j,itr)*d_1000*regrav
            ddsfc(i,j,itr) = 1.D6*ddsfc(i,j,itr)*d_1000*regrav
            wxsg(i,j,itr) = 1.D6*wxsg(i,j,itr)*d_1000*regrav
            wxaq(i,j,itr) = 1.D6*wxaq(i,j,itr)*d_1000*regrav
!           emtrac isbuilt from chsurfem so just need the 1e6*dt/2
!           factor to to pass im mg/m2
            cemtrac(i,j,itr) = 1.D6*cemtrac(i,j,itr)
          end do
        end do
      end do

      end subroutine tracbud
!
      end module mod_che_tend
