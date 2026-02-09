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

module mod_che_mppio
!
  use mod_intkinds
  use mod_realkinds
  use mod_mppparam
  use mod_runparams
  use mod_dynparam
  use mod_memutil
  use mod_mpmessage
  use mod_che_param
  use mod_che_common
  use mod_che_species

  implicit none (type, external)

  private

  real(rkx), pointer, contiguous, public, dimension(:,:,:,:) :: rainout_io, washout_io
  real(rkx), pointer, contiguous, public, dimension(:,:,:) :: remdrd_io, convpr_io
  real(rkx), pointer, contiguous, public, dimension(:,:,:) :: duflux_io, voflux_io
  real(rkx), pointer, contiguous, public, dimension(:,:) :: ssw2da_io, sdelt_io,   &
                                         sdelq_io, sfracv2d_io, &
                                         sfracb2d_io, sfracs2d_io, &
                                         svegfrac2d_io
  real(rkx), pointer, contiguous, public, dimension(:,:,:,:) :: trac_io
  real(rkx), pointer, contiguous, public, dimension(:,:,:,:) :: chia_io, chib_io
  real(rkx), pointer, contiguous, public, dimension(:,:,:,:) :: chemall_io
  real(rkx), pointer, contiguous, public, dimension(:,:,:,:) :: taucldsp_io

  public :: allocate_mod_che_mppio

  contains
    !
    ! This routines allocate all the arrays contained in the module
    !
    subroutine allocate_mod_che_mppio
      implicit none (type, external)

      if ( ichem == 1 ) then

        if ( do_parallel_save ) then
          if ( idynamic == 3 ) then
            call getmem(trac_io,jce1,jce2,ice1,ice2, &
                          1,kz,1,ntr,'che_mppio:trac_io')
          else
            call getmem(chia_io,jce1,jce2,ice1,ice2, &
                          1,kz,1,ntr,'che_mppio:chia_io')
            call getmem(chib_io,jce1,jce2,ice1,ice2, &
                          1,kz,1,ntr,'che_mppio:chib_io')
          end if
          call getmem(convpr_io,jci1,jci2,ici1,ici2, &
                        1,kz,'che_mppio:convpr_io')
          call getmem(rainout_io,jci1,jci2,ici1,ici2, &
                        1,kz,1,ntr,'che_mppio:rainout_io')
          call getmem(washout_io,jci1,jci2,ici1,ici2, &
                        1,kz,1,ntr,'che_mppio:washout_io')
          call getmem(remdrd_io,jci1,jci2,ici1,ici2, &
                        1,ntr,'che_mppio:remdrd_io')
          call getmem(ssw2da_io,jci1,jci2,ici1,ici2, &
                        'che_mppio:ssw2da_io')
#ifdef CLM45
          call getmem(duflux_io,jci1,jci2,ici1,ici2, &
                        1,4,'che_mppio:duflux_io')
          call getmem(voflux_io,jci1,jci2,ici1,ici2, &
                        1,ntr,'che_mppio:voflux_io')
#else
          call getmem(sdelq_io,jci1,jci2,ici1,ici2, &
                        'che_mppio:sdelq_io')
          call getmem(sdelt_io,jci1,jci2,ici1,ici2, &
                        'che_mppio:sdelt_io')
          call getmem(svegfrac2d_io,jci1,jci2,ici1,ici2, &
                        'che_mppio:svegfrac2d_io')
#endif
          call getmem(sfracb2d_io,jci1,jci2,ici1,ici2, &
                        'che_mppio:sfracb2d_io')
          call getmem(sfracs2d_io,jci1,jci2,ici1,ici2, &
                        'che_mppio:sfracs2d_io')
          call getmem(sfracv2d_io,jci1,jci2,ici1,ici2, &
                        'che_mppio:sfracv2d_io')
          if ( igaschem == 1 .and. ichsolver > 0 ) then
            call getmem(chemall_io,jci1,jci2,ici1,ici2, &
                          1,kz,1,totsp,'che_mppio:chemall_io')
            call getmem(taucldsp_io,jci1,jci2,ici1,ici2, &
                          0,kz,1,nspi,'che_mppio:taucldsp_io')
          end if

          return

        end if

        if ( myid == iocpu ) then
          call getmem(convpr_io,jcross1,jcross2,icross1,icross2, &
                        1,kz,'che_mppio:convpr_io')
          call getmem(rainout_io,jcross1,jcross2,icross1,icross2, &
                        1,kz,1,ntr,'che_mppio:rainout_io')
          call getmem(washout_io,jcross1,jcross2,icross1,icross2, &
                        1,kz,1,ntr,'che_mppio:washout_io')
          call getmem(remdrd_io,jcross1,jcross2,icross1,icross2, &
                        1,ntr,'che_mppio:remdrd_io')

          call getmem(ssw2da_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:ssw2da_io')
#ifdef CLM45
          call getmem(duflux_io,jcross1,jcross2,icross1,icross2, &
                        1,4,'che_mppio:duflux_io')
          call getmem(voflux_io,jcross1,jcross2,icross1,icross2, &
                        1,ntr,'che_mppio:voflux_io')
#else
          call getmem(sdelq_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:sdelq_io')
          call getmem(sdelt_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:sdelt_io')
          call getmem(svegfrac2d_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:svegfrac2d_io')
#endif
          call getmem(sfracb2d_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:sfracb2d_io')
          call getmem(sfracs2d_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:sfracs2d_io')
          call getmem(sfracv2d_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:sfracv2d_io')

          if ( idynamic == 3 ) then
            call getmem(trac_io,jcross1,jcross2,icross1,icross2, &
                          1,kz,1,ntr,'che_mppio:trac_io')
          else
            call getmem(chia_io,jcross1,jcross2,icross1,icross2, &
                          1,kz,1,ntr,'che_mppio:chia_io')
            call getmem(chib_io,jcross1,jcross2,icross1,icross2, &
                          1,kz,1,ntr,'che_mppio:chib_io')
          end if
          if ( igaschem == 1 .and. ichsolver > 0 ) then
            call getmem(chemall_io,jcross1,jcross2,icross1,icross2, &
                          1,kz,1,totsp,'che_mppio:chemall_io')
            call getmem(taucldsp_io,jcross1,jcross2,icross1,icross2, &
                          0,kz,1,nspi,'che_mppio:taucldsp_io')
          end if
        end if
      end if
    end subroutine allocate_mod_che_mppio
!
end module mod_che_mppio
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
