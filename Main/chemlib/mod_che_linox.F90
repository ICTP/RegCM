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

module mod_che_linox
  !
  ! Tracer convective transport
  !
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_che_common
  use mod_dynparam

  implicit none

  private

  public :: linox_em , allocate_mod_che_linox

  real(rkx) , pointer, dimension(:,:) :: flashrate

contains

  subroutine allocate_mod_che_linox
    implicit none
    if ( ichem == 1 ) then
       call getmem2d(flashrate,jce1,jce2,ice1,ice2,'che_linox:flashrate')
    end if
  end subroutine allocate_mod_che_linox

  subroutine linox_em(j,ivegcov)

    implicit none 

    integer(ik4) , intent(in) :: j
    integer(ik4) , dimension(ici1:ici2), intent(in) :: ivegcov
    integer(ik4) :: i,k,kcumzero,kcumm10
    real(rkx), dimension(kz) :: znox_prod_ic, znox_prod_cg, amass 
    real(rkx):: zsum_lmass_ic ,zsum_lmass_cg, &
         pic_rate, pcg_rate, cloud_top_height,zthickice,zbeta, fac 



    do i = ici1,ici2 

       ! cloud top heigh from convection scheme
       if (kcumtop(j,i) <= 0) cycle
       cloud_top_height  =  cza(j,i,kcumtop(j,i)) + 0.5 * cdzq(j,i,kcumtop(j,i))
       
       ! Lightning frequency for continentaland maritime  clouds as a function of cloud
       ! top height : Price and Rind 

       if (ivegcov(i) > 0 ) then 
          ! Flashes/min over land boxes
          flashrate (j,i)   = 3.44e-5 * ( ( cloud_top_height * 1.e-3 )**4.9  )
       else if (ivegcov(i) == 0) then 
          flashrate (j,i)  = 6.4e-4  * ( ( cloud_top_height  * 1e-3 )**1.73 )
       end if
       ! flash per second parameterization Price and Rind
       flashrate(j,i) = flashrate(j,i)/60.
       ! apply a scale factor ( Price and Rind, 1994) 
       fac  = dx**2 * (mathpi/180.)**2 / (earthrad**2 *cos(cxlat(j,i) * mathpi/180.)) 
       flashrate(j,i) = flashrate(j,i)*0.97241*exp(0.048203*fac) 

       ! IC / CG flash ratio in function of icy cloud thickness  
       kcumzero =0
       icelev:  do  k = kcumbot(j,i), kcumtop(j,i),-1
          if (ctb3d(j,i,k) < tzero ) then 
             kcumzero = k 
             exit icelev
          end if
       end do icelev

       ! if convective cloud is not icy, cycle to the next i column 
       if (kcumzero ==0 .or. kcumzero == kcumtop(j,i) ) cycle 

       kcumm10 =0
       minus10:  do  k = kcumzero, kcumtop(j,i),-1
          if (ctb3d(j,i,k) < tzero - 10. ) then 
             kcumm10 = k 
             exit minus10
          end if
       end do minus10


       zthickice =  (cloud_top_height -  cza(j,i,kcumzero) ) *1.e-3

       zbeta = (((0.021*zthickice -  0.648) &
            *zthickice   +  7.493) &
            *zthickice  - 36.540) &
            *zthickice + 63.090
       zbeta = min( 48.7, max( 0.19, zbeta ) ) ! 5.5km < zthickice < 14km
       !
      
       znox_prod_ic(:) = 0.0
       znox_prod_cg(:) = 0.0
       amass(:) = 0.0
       zsum_lmass_ic = 0.0
       zsum_lmass_cg = 0.0

       ! IC_NOx production rate
       !

               
       pic_rate = (zbeta / (1.0+zbeta))*flashrate(j,i)
       do k = kcumzero, kcumtop(j,i), -1
          amass(k) = crhob3d(j,i,k)*  cdzq(j,i,k) * dx**2
          zsum_lmass_ic   = zsum_lmass_ic + amass(k)
          znox_prod_ic(k) = pic_rate * 6.7e25 * amass(k)
       end do
       znox_prod_ic(:) = znox_prod_ic(:) / zsum_lmass_ic
       !
       ! cg_nox production rate
       ! ic can happen without cg
       if (kcumm10 > 0) then  
          pcg_rate  = (1.0 / (1.0+zbeta))*flashrate(j,i)
          do k = kz, kcumm10, -1             
             amass(k) = crhob3d(j,i,k)*  cdzq(j,i,k) * dx**2
             zsum_lmass_cg   = zsum_lmass_cg + amass(k)
             znox_prod_cg(k) = pcg_rate * 6.7e26 * amass(k)
          end do
          znox_prod_cg(:)  = znox_prod_cg(:)  / zsum_lmass_cg
       end if

       ! Update units (molecules NO /s => kg.kg-1.s-1)
       ! update emission tendencies and diag
       ! 
       !               ------------------------------------------------
       !!                  
       where(amass(:) > 0.0)
          znox_prod_ic(:) = znox_prod_ic(:) * 30.e-3 / (navgdr * amass(:))  
          znox_prod_cg(:) =  znox_prod_cg(:)* 30.e-3 / (navgdr * amass(:))  
       end where


       chiten(j,i,:,ino) =   chiten(j,i,:,ino) + (znox_prod_ic(:) + znox_prod_cg(:)) * cpsb(j,i)

       if ( ichdiag == 1 ) then
          cemisdiag(j,i,:,ino) = cemisdiag(j,i,:,ino ) + &
               (znox_prod_ic(:) +  znox_prod_cg(:)) * cfdout
       end if

       !
    end do ! en loop on i 


  end subroutine linox_em

end module mod_che_linox
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
