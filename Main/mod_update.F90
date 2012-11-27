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
      use ESMF
!
      use mod_couplerr
      use mod_bats_common
!
      implicit none
      private
!
!-----------------------------------------------------------------------
!     Public subroutines 
!-----------------------------------------------------------------------
!
      public :: RCM_PutExportData
      public :: RCM_GetImportData
!
      contains
!
      subroutine RCM_PutExportData (localPet)
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use mod_mppparam, only : ma
      use mod_constants, only : wlhv
      use mod_dynparam, only : ice1, ice2, jce1, jce2
      use mod_dynparam, only : ici1, ici2, jci1, jci2, kz
      use mod_dynparam, only : global_istart, global_jstart
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
      integer :: i, j, ii, jj, k, n, rc
      character (len=40) :: name
      character (len=100) :: outfile
!
!-----------------------------------------------------------------------
!     Initialize the import and export fields 
!-----------------------------------------------------------------------
!
      do n = 1, nNest(Iatmos)
!
!-----------------------------------------------------------------------
!     Create export state fields 
!-----------------------------------------------------------------------
!
      do k = 1, size(models(Iatmos)%dataExport(:,n), dim=1)
!
!-----------------------------------------------------------------------
!     Debug: write size of pointers    
!-----------------------------------------------------------------------
!
      name = trim(adjustl(models(Iatmos)%dataExport(k,n)%name))
      if (cpl_dbglevel > 1) then
        write(*,50) localPet, 0, adjustl("PTR/ATM/EXP/"//name),         &
                    lbound(models(Iatmos)%dataExport(k,n)%ptr, dim=1),  &
                    ubound(models(Iatmos)%dataExport(k,n)%ptr, dim=1),  &
                    lbound(models(Iatmos)%dataExport(k,n)%ptr, dim=2),  &
                    ubound(models(Iatmos)%dataExport(k,n)%ptr, dim=2)
        write(*,50) localPet, 0, adjustl("DAT/ATM/EXP/"//name),         &
                    lbound(fbat(:,:,ps_o), dim=1),                      &
                    ubound(fbat(:,:,ps_o), dim=1),                      &
                    lbound(fbat(:,:,ps_o), dim=2),                      &
                    ubound(fbat(:,:,ps_o), dim=2)
        write(*,50) localPet, 0, adjustl("IND/ATM/EXP/"//name),         &
                    jci1, jci2, ici1, ici2
      end if
!
!-----------------------------------------------------------------------
!     Fill pointers with data 
!-----------------------------------------------------------------------
!
      select case (trim(adjustl(name)))
!
!-----------------------------------------------------------------------
!     Surface atmospheric pressure (mb)
!-----------------------------------------------------------------------
!
      case ('Pair')
      do i = ici1, ici2
        do j = jci1, jci2
        ii = global_istart+i-1
        jj = global_jstart+j-1
        models(Iatmos)%dataExport(k,n)%ptr(ii,jj) = fbat(j,i,ps_o)
        end do
      end do
!
!-----------------------------------------------------------------------
!     Surface (2m) air temperature (K)
!-----------------------------------------------------------------------
!
      case ('Tair')
      do i = ici1, ici2
        do j = jci1, jci2
        ii = global_istart+i-1
        jj = global_jstart+j-1
        models(Iatmos)%dataExport(k,n)%ptr(ii,jj) = fbat(j,i,t2m_o)
        end do
      end do
!
!-----------------------------------------------------------------------
!     Surface (2m) specific humidity (kg/kg)
!-----------------------------------------------------------------------
!           
      case ('Qair')
      do i = ici1, ici2
        do j = jci1, jci2
        ii = global_istart+i-1
        jj = global_jstart+j-1
        models(Iatmos)%dataExport(k,n)%ptr(ii,jj) = fbat(j,i,q2m_o)
        end do
      end do
!          
!-----------------------------------------------------------------------
!     Shortwave radiation (W/m2) 
!-----------------------------------------------------------------------
! 
      case ('Swrad')
      do i = ici1, ici2
        do j = jci1, jci2
        ii = global_istart+i-1
        jj = global_jstart+j-1
        models(Iatmos)%dataExport(k,n)%ptr(ii,jj) = fbat(j,i,fswa_o)
        end do
      end do
!
!-----------------------------------------------------------------------
!     Downward longwave radiation (W/m2) 
!-----------------------------------------------------------------------
! 
      case ('Lwrad_down')
      do i = ici1, ici2
        do j = jci1, jci2
        ii = global_istart+i-1
        jj = global_jstart+j-1
        models(Iatmos)%dataExport(k,n)%ptr(ii,jj) = fbat(j,i,flwd_o)
        end do
      end do
!          
!-----------------------------------------------------------------------
!     Net longwave radiation (W/m2) 
!-----------------------------------------------------------------------
! 
      case ('Lwrad')
      do i = ici1, ici2
        do j = jci1, jci2
        ii = global_istart+i-1
        jj = global_jstart+j-1
        models(Iatmos)%dataExport(k,n)%ptr(ii,jj) = fbat(j,i,flwa_o)
        end do
      end do
!          
!-----------------------------------------------------------------------
!     Sensible heat flux (W/m2) 
!-----------------------------------------------------------------------
! 
      case ('Shflx')
      do i = ici1, ici2
        do j = jci1, jci2
        ii = global_istart+i-1
        jj = global_jstart+j-1
        models(Iatmos)%dataExport(k,n)%ptr(ii,jj) = fbat(j,i,sena_o)
        end do
      end do
!          
!-----------------------------------------------------------------------
!     Latent heat flux (W/m2) 
!-----------------------------------------------------------------------
! 
      case ('Lhflx')
      do i = ici1, ici2
        do j = jci1, jci2
        ii = global_istart+i-1
        jj = global_jstart+j-1
        models(Iatmos)%dataExport(k,n)%ptr(ii,jj) = fbat(j,i,evpa_o)*   &
                                                    wlhv*day2s
        end do
      end do
!
!-----------------------------------------------------------------------
!     Precipitation (m/s)
!-----------------------------------------------------------------------
! 
      case ('Rain')
      do i = ici1, ici2
        do j = jci1, jci2
        ii = global_istart+i-1
        jj = global_jstart+j-1
        models(Iatmos)%dataExport(k,n)%ptr(ii,jj) = fbat(j,i,tpr_o)*    &
                                                    day2s
        end do
      end do
!          
!-----------------------------------------------------------------------
!     Surface (10m) U-wind speed (m/s)
!-----------------------------------------------------------------------
! 
      case ('Uwind')
      do i = ici1, ici2
        do j = jci1, jci2
        ii = global_istart+i-1
        jj = global_jstart+j-1
        models(Iatmos)%dataExport(k,n)%ptr(ii,jj) = fbat(j,i,u10m_o)
        end do
      end do
!          
!-----------------------------------------------------------------------
!     Surface (10 m) V-wind speed (m/s)
!-----------------------------------------------------------------------
! 
      case ('Vwind')
      do i = ici1, ici2
        do j = jci1, jci2
        ii = global_istart+i-1
        jj = global_jstart+j-1
        models(Iatmos)%dataExport(k,n)%ptr(ii,jj) = fbat(j,i,v10m_o)
        end do
      end do
!          
!-----------------------------------------------------------------------
!     Net freshwater flux (m/s) 
!-----------------------------------------------------------------------
! 
      case ('EminP')
      do i = ici1, ici2
        do j = jci1, jci2
        ii = global_istart+i-1
        jj = global_jstart+j-1
        models(Iatmos)%dataExport(k,n)%ptr(ii,jj) = (fbat(j,i,evpa_o)-  &
                                                    fbat(j,i,tpr_o))*   &
                                                    day2s
        end do
      end do
!          
!-----------------------------------------------------------------------
!     Net heat flux (W/m2) 
!-----------------------------------------------------------------------
!
      case ('NHeat')
      do i = ici1, ici2
        do j = jci1, jci2
        ii = global_istart+i-1
        jj = global_jstart+j-1
        models(Iatmos)%dataExport(k,n)%ptr(ii,jj) = fbat(j,i,fswa_o)-   &
                                                    fbat(j,i,sena_o)-   &
                                                    fbat(j,i,evpa_o)*   &
                                                    wlhv*day2s-         &
                                                    fbat(j,i,flwa_o)
        end do
      end do
!          
!-----------------------------------------------------------------------
!     Surface (10m) U-wind stress (Pa) 
!-----------------------------------------------------------------------
!
      case ('Ustr')
      do i = ici1, ici2
        do j = jci1, jci2
        ii = global_istart+i-1
        jj = global_jstart+j-1
        models(Iatmos)%dataExport(k,n)%ptr(ii,jj) = taux(1,j,i)
        end do
      end do
!          
!-----------------------------------------------------------------------
!     Surface (10m) V-wind stress (Pa) 
!-----------------------------------------------------------------------
!
      case ('Vstr')
      do i = ici1, ici2
        do j = jci1, jci2
        ii = global_istart+i-1
        jj = global_jstart+j-1
        models(Iatmos)%dataExport(k,n)%ptr(ii,jj) = tauy(1,j,i)
        end do
      end do
      end select
!
!-----------------------------------------------------------------------
!     Fill domain boundaries with data   
!-----------------------------------------------------------------------
!
      if (ma%has_bdytop) then
        ii = global_istart+ici2-1
        models(Iatmos)%dataExport(k,n)%ptr(ii+1,:) =                    &
                       models(Iatmos)%dataExport(k,n)%ptr(ii,:)
        models(Iatmos)%dataExport(k,n)%ptr(ii+2,:) =                    &
                       models(Iatmos)%dataExport(k,n)%ptr(ii,:)
      end if
!
      if (ma%has_bdyright) then
        jj = global_jstart+jci2-1
        models(Iatmos)%dataExport(k,n)%ptr(:,jj+1) =                    &
                       models(Iatmos)%dataExport(k,n)%ptr(:,jj)
        models(Iatmos)%dataExport(k,n)%ptr(:,jj+2) =                    &
                       models(Iatmos)%dataExport(k,n)%ptr(:,jj)
      end if
!
      if (ma%has_bdybottom) then
        ii = global_istart+ici1-1
        models(Iatmos)%dataExport(k,n)%ptr(ii-1,:) =                    &
                       models(Iatmos)%dataExport(k,n)%ptr(ii,:)
      end if
!
      if (ma%has_bdyleft) then
        jj = global_jstart+jci1-1
        models(Iatmos)%dataExport(k,n)%ptr(:,jj-1) =                    &
                       models(Iatmos)%dataExport(k,n)%ptr(:,jj)
      end if
!
!-----------------------------------------------------------------------
!     Debug: write field to ASCII file    
!-----------------------------------------------------------------------
!
      if (cpl_dbglevel > 3) then
      write(outfile,                                                    &
            fmt='(A10,"_",A,"_",I4,"-",I2.2,"-",                        &
            I2.2,"_",I2.2,"_",I2.2,".txt")')                            &
            'atm_export',                                               &
            trim(adjustl(name)),                                        &
            models(Iatmos)%time%year,                                   &
            models(Iatmos)%time%month,                                  &
            models(Iatmos)%time%day,                                    &
            models(Iatmos)%time%hour,                                   &
            localPet
!
      open (unit=99, file = trim(outfile)) 
      call print_matrix_r8(models(Iatmos)%dataExport(k,n)%ptr,          &
                           1, 1, localPet, 99, "PTR/ATM/EXP")
      close(99)
      end if
!
!-----------------------------------------------------------------------
!     Debug: write field to file
!-----------------------------------------------------------------------
!
      if (cpl_dbglevel > 2) then
      write(outfile,                                                    &
            fmt='(A10,"_",A,"_",I4,"-",I2.2,"-",I2.2,"_",I2.2,".nc")')  &
            'atm_export',                                               &
            trim(adjustl(name)),                                        &
            models(Iatmos)%time%year,                                   &
            models(Iatmos)%time%month,                                  &
            models(Iatmos)%time%day,                                    &
            models(Iatmos)%time%hour
!
      call ESMF_FieldWrite(models(Iatmos)%dataExport(k,n)%field,        &
                           trim(adjustl(outfile)),                      &
                           rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
!
      end do
      end do
!
!-----------------------------------------------------------------------
!     Format definition 
!-----------------------------------------------------------------------
!     
 50   format(" PET(",I3,") - DE(",I2,") - ", A20, " : ", 4I8)
!
!-----------------------------------------------------------------------
!     Set return flag to success
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine RCM_PutExportData
!
      subroutine RCM_GetImportData (localPet)
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use mod_atm_interface, only : sfs
      use mod_runparams, only : idcsst 
      use mod_dynparam, only : global_istart, global_jstart
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
      character (len=40) :: name
      character (len=100) :: outfile
      integer :: i, j, ii, jj, k, m, n, p, rc, localDECount
      real*8 :: hice, toth, tol, scale_factor, add_offset
#ifdef ROMSICE
      ! minimum ice depth in mm: less that this is removed
      real(dp) , parameter :: iceminh = d_10
      ! reference hgt in mm for latent heat removal from ice
      real(dp) , parameter :: href = d_two * iceminh
      ! steepness factor of latent heat removal
      real(dp) , parameter :: steepf = 1.0D0  ! Tuning needed
#endif
!
      real(ESMF_KIND_R8), pointer :: ptr(:,:)
!
!-----------------------------------------------------------------------
!     Loop over number of nested/composed meshs 
!-----------------------------------------------------------------------
!
      do p = 1, nNest(Iatmos)
!
!-----------------------------------------------------------------------
!     Get import fields 
!-----------------------------------------------------------------------
!
      do k = 1, size(models(Iatmos)%dataImport(:,p), dim=1)
      name = models(Iatmos)%dataImport(k,p)%name
!
!-----------------------------------------------------------------------
!     Get number of local DEs
!-----------------------------------------------------------------------
! 
      call ESMF_GridGet (models(Iatmos)%grid(p),                        &
                         localDECount=localDECount,                     &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get pointer
!-----------------------------------------------------------------------
! 
      do m = 0, localDECount-1
      call ESMF_FieldGet (models(Iatmos)%dataImport(k,p)%field,         &
                          localDE=m,                                    &
                          farrayPtr=ptr,                                &
                          rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Debug: write size of pointers    
!-----------------------------------------------------------------------
!
      if (cpl_dbglevel > 1) then
        write(*,60) localPet, m, adjustl("PTR/ATM/IMP/"//name),         &
                    lbound(ptr, dim=1), ubound(ptr, dim=1),             &
                    lbound(ptr, dim=2), ubound(ptr, dim=2)
      end if
!
!-----------------------------------------------------------------------
!     Put data to RCM variable
!-----------------------------------------------------------------------
!
      tol = MISSING_R8/2.0d0
      scale_factor = models(Iatmos)%dataImport(k,p)%scale_factor
      add_offset = models(Iatmos)%dataImport(k,p)%add_offset
!
      select case (trim(adjustl(name)))
      case('SST')      
      if (cpl_dbglevel > 1) then
        write(*,60) localPet, m, adjustl("DAT/ATM/IMP/"//name),         &
                    lbound(sfs%tga, dim=1), ubound(sfs%tga, dim=1),     &
                    lbound(sfs%tga, dim=2), ubound(sfs%tga, dim=2)
      end if
!
      if ( idcsst == 1 ) then
      do i = ici1, ici2
        do j = jci1, jci2
          ii = global_istart+i-1
          jj = global_jstart+j-1
          do n = 1 , nnsg
            if (iveg1(n,j,i) == 12 .or.                                 &
                iveg1(n,j,i) == 14 .or.                                 &
                iveg1(n,j,i) == 15) then
               if (ptr(ii,jj) .lt. tol) then
                  sfs%tga(j,i) = (ptr(ii,jj)*scale_factor)+             &
                                 add_offset+dtskin(j,i)
                  sfs%tgb(j,i) = (ptr(ii,jj)*scale_factor)+             &
                                 add_offset+dtskin(j,i)
               end if
             end if
           end do
         end do
      end do
      else
      do i = ici1, ici2
        do j = jci1, jci2
          ii = global_istart+i-1
          jj = global_jstart+j-1
          do n = 1 , nnsg
            if (iveg1(n,j,i) == 12 .or.                                 &
                iveg1(n,j,i) == 14 .or.                                 &
                iveg1(n,j,i) == 15) then
              if (ptr(ii,jj) .lt. tol) then 
                 sfs%tga(j,i) = (ptr(ii,jj)*scale_factor)+add_offset
                 sfs%tgb(j,i) = (ptr(ii,jj)*scale_factor)+add_offset 
              end if
            end if
          end do
        end do
      end do
      end if
#ifdef ROMSICE
      case('Hice')      
      if (cpl_dbglevel > 1) then
        write(*,60) localPet, m, adjustl("DAT/ATM/IMP/"//name),         &
                    lbound(sfs%hfx, dim=1), ubound(sfs%hfx, dim=1),     &
                    lbound(sfs%hfx, dim=2), ubound(sfs%hfx, dim=2)
      end if

      do i = ici1, ici2
        do j = jci1, jci2
          ii = global_istart+i-1
          jj = global_jstart+j-1
          do n = 1 , nnsg
            if (iveg1(n,j,i) == 12 .or.                                 &
                iveg1(n,j,i) == 14 .or.                                 &
                iveg1(n,j,i) == 15) then
              if (ptr(ii,jj) .lt. tol) then
                hice = (ptr(ii,jj)*scale_factor)+add_offset  
                if (hice .gt. iceminh) then
                  ! change land-use type as ice covered
                  ldmsk1(n,j,i) = 2
                  lveg(n,j,i) = 12
                  ! reduce sensible heat flux for ice presence                    
                  toth = hice+fbat(j,i,scv_o) 
                  if ( toth > href ) then
                    sfs%hfx(j,i) = sfs%hfx(j,i)*(href/toth)**steepf
                  end if
                else
                  ! change land-use type to its original value
                  ldmsk1(n,j,i) = 0
                  lveg(n,j,i) = iveg1(n,j,i)
                end if
              end if
            end if
          end do
        end do
      end do
#endif
      end select
      end do
!
!-----------------------------------------------------------------------
!     Debug: write field to ASCII file    
!-----------------------------------------------------------------------
!
      if (cpl_dbglevel > 3) then
      write(outfile,                                                    &
            fmt='(A10,"_",A,"_",I4,"-",I2.2,"-",                        &
            I2.2,"_",I2.2,"_",I2.2,".txt")')                            &
            'atm_import',                                               &
            trim(adjustl(name)),                                        &
            models(Iatmos)%time%year,                                   &
            models(Iatmos)%time%month,                                  &
            models(Iatmos)%time%day,                                    &
            models(Iatmos)%time%hour,                                   &
            localPet
!
      open (unit=99, file = trim(outfile))
      call print_matrix_r8(ptr, 1, 1, localPet, 99, "PTR/ATM/IMP")
      close(99)
      end if
!
!-----------------------------------------------------------------------
!     Debug: write field to NetCDF file
!-----------------------------------------------------------------------
!
      if (cpl_dbglevel > 2) then
      write(outfile,                                                    &
            fmt='(A10,"_",A,"_",I4,"-",I2.2,"-",I2.2,"_",I2.2,".nc")')  &
            'atm_import',                                               &
            trim(adjustl(name)),                                        &
            models(Iatmos)%time%year,                                   &
            models(Iatmos)%time%month,                                  &
            models(Iatmos)%time%day,                                    &
            models(Iatmos)%time%hour
!
      call ESMF_FieldWrite(models(Iatmos)%dataImport(k,p)%field,        &
                           trim(adjustl(outfile)),                      &
                           rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
!
!-----------------------------------------------------------------------
!     Nullify pointer to make sure that it does not point on a random 
!     part in the memory 
!-----------------------------------------------------------------------
!
      if (associated(ptr)) then
        nullify(ptr)
      end if
      end do
      end do
!
!-----------------------------------------------------------------------
!     Format definition 
!-----------------------------------------------------------------------
!
 60   format(" PET(",I3,") - DE(",I2,") - ", A20, " : ", 4I8)
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine RCM_GetImportData
!
      end module mod_update
