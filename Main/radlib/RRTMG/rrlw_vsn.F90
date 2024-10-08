      module rrlw_vsn

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_lw version information

! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
!hnamrtm :character:
!hnamini :character:
!hnamcld :character:
!hnamclc :character:
!hnamrtr :character:
!hnamrtx :character:
!hnamrtc :character:
!hnamset :character:
!hnamtau :character:
!hnamatm :character:
!hnamutl :character:
!hnamext :character:
!hnamkg  :character:
!
! hvrrtm :character:
! hvrini :character:
! hvrcld :character:
! hvrclc :character:
! hvrrtr :character:
! hvrrtx :character:
! hvrrtc :character:
! hvrset :character:
! hvrtau :character:
! hvratm :character:
! hvrutl :character:
! hvrext :character:
! hvrkg  :character:
!------------------------------------------------------------------

      character(len=18) hvrrtm,hvrini,hvrcld,hvrclc,hvrrtr,hvrrtx, &
                   hvrrtc,hvrset,hvrtau,hvratm,hvrutl,hvrext
      character(len=20) hnamrtm,hnamini,hnamcld,hnamclc,hnamrtr,hnamrtx, &
                   hnamrtc,hnamset,hnamtau,hnamatm,hnamutl,hnamext

      character(len=18) hvrkg
      character(len=20) hnamkg

      end module rrlw_vsn


! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
