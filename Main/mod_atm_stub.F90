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

module mod_atm_stub

  use mod_dynparam
  use mod_stdio
  use mod_runparams
  use mod_mppparam
  use mod_mpmessage
  use mod_service
  use mod_memutil
  use mod_regcm_types

  implicit none

  private

  type(domain) , public :: mddom
  type(domain_subgrid) , public :: mdsub
  type(surfstate) , public :: sfs

  integer(ik4) :: ix1 , ix2 , jx1 , jx2
  integer(ik4) :: id1 , id2 , jd1 , jd2

  public :: allocate_mod_atm_interface
  public :: setup_model_indexes

  contains

    subroutine setup_model_indexes
      implicit none
      ma%jbl1 = 1
      ma%jbl2 = 2
      ma%jbl3 = 3
      ma%jbl4 = 4
      ma%jbl6 = 6
      ma%jbr1 = 1
      ma%jbr2 = 2
      ma%jbr3 = 3
      ma%jbr4 = 4
      ma%jbr6 = 6
      ma%ibt1 = 1
      ma%ibt2 = 2
      ma%ibt3 = 3
      ma%ibt4 = 4
      ma%ibt6 = 6
      ma%ibb1 = 1
      ma%ibb2 = 2
      ma%ibb3 = 3
      ma%ibb4 = 4
      ma%ibb6 = 6
      if ( ma%has_bdyleft ) then
        ma%jbl1 = 0
        ma%jbl2 = 0
        ma%jbl3 = 0
        ma%jbl4 = 0
        ma%jbl6 = 0
      end if
      if ( ma%has_bdyright ) then
        ma%jbr1 = 0
        ma%jbr2 = 0
        ma%jbr3 = 0
        ma%jbr4 = 0
        ma%jbr6 = 0
      end if
      if ( ma%has_bdytop ) then
        ma%ibt1 = 0
        ma%ibt2 = 0
        ma%ibt3 = 0
        ma%ibt4 = 0
        ma%ibt6 = 0
      end if
      if ( ma%has_bdybottom ) then
        ma%ibb1 = 0
        ma%ibb2 = 0
        ma%ibb3 = 0
        ma%ibb4 = 0
        ma%ibb6 = 0
      end if
      jde1  = global_dot_jstart
      jdi1  = global_dot_jstart
      jdii1 = global_dot_jstart
      jde2  = global_dot_jend
      jdi2  = global_dot_jend
      jdii2 = global_dot_jend
      ide1  = global_dot_istart
      idi1  = global_dot_istart
      idii1 = global_dot_istart
      ide2  = global_dot_iend
      idi2  = global_dot_iend
      idii2 = global_dot_iend
      if ( ma%has_bdyleft ) then
        jdi1 = jde1 + 1
        jdii1 = jde1 + 2
      end if
      if ( ma%has_bdyright ) then
        jdi2 = jde2 - 1
        jdii2 = jde2 - 2
      end if
      if ( ma%has_bdybottom ) then
        idi1 = ide1 + 1
        idii1 = ide1 + 2
      end if
      if ( ma%has_bdytop ) then
        idi2 = ide2 - 1
        idii2 = ide2 - 2
      end if
      jce1  = global_cross_jstart
      jci1  = global_cross_jstart
      jcii1 = global_cross_jstart
      jce2  = global_cross_jend
      jci2  = global_cross_jend
      jcii2 = global_cross_jend
      ice1  = global_cross_istart
      ici1  = global_cross_istart
      icii1 = global_cross_istart
      ice2  = global_cross_iend
      ici2  = global_cross_iend
      icii2 = global_cross_iend
      if ( ma%has_bdyleft ) then
        jci1 = jce1 + 1
        jcii1 = jce1 + 2
      end if
      if ( ma%has_bdyright ) then
        jci2 = jce2 - 1
        jcii2 = jce2 - 2
      end if
      if ( ma%has_bdybottom ) then
        ici1 = ice1 + 1
        icii1 = ice1 + 2
      end if
      if ( ma%has_bdytop ) then
        ici2 = ice2 - 1
        icii2 = ice2 - 2
      end if
      idi1ga = idi1 - ma%ibb1
      idi2ga = idi2 + ma%ibt1
      jdi1ga = jdi1 - ma%jbl1
      jdi2ga = jdi2 + ma%jbr1
      ide1ga = ide1 - ma%ibb1
      ide2ga = ide2 + ma%ibt1
      jde1ga = jde1 - ma%jbl1
      jde2ga = jde2 + ma%jbr1
      ici1ga = ici1 - ma%ibb1
      ici2ga = ici2 + ma%ibt1
      jci1ga = jci1 - ma%jbl1
      jci2ga = jci2 + ma%jbr1
      ice1ga = ice1 - ma%ibb1
      ice2ga = ice2 + ma%ibt1
      jce1ga = jce1 - ma%jbl1
      jce2ga = jce2 + ma%jbr1

      idi1gb = idi1 - ma%ibb2
      idi2gb = idi2 + ma%ibt2
      jdi1gb = jdi1 - ma%jbl2
      jdi2gb = jdi2 + ma%jbr2
      ide1gb = ide1 - ma%ibb2
      ide2gb = ide2 + ma%ibt2
      jde1gb = jde1 - ma%jbl2
      jde2gb = jde2 + ma%jbr2
      ici1gb = ici1 - ma%ibb2
      ici2gb = ici2 + ma%ibt2
      jci1gb = jci1 - ma%jbl2
      jci2gb = jci2 + ma%jbr2
      ice1gb = ice1 - ma%ibb2
      ice2gb = ice2 + ma%ibt2
      jce1gb = jce1 - ma%jbl2
      jce2gb = jce2 + ma%jbr2

      idi1gc = idi1 - ma%ibb3
      idi2gc = idi2 + ma%ibt3
      jdi1gc = jdi1 - ma%jbl3
      jdi2gc = jdi2 + ma%jbr3
      ide1gc = ide1 - ma%ibb3
      ide2gc = ide2 + ma%ibt3
      jde1gc = jde1 - ma%jbl3
      jde2gc = jde2 + ma%jbr3
      ici1gc = ici1 - ma%ibb3
      ici2gc = ici2 + ma%ibt3
      jci1gc = jci1 - ma%jbl3
      jci2gc = jci2 + ma%jbr3
      ice1gc = ice1 - ma%ibb3
      ice2gc = ice2 + ma%ibt3
      jce1gc = jce1 - ma%jbl3
      jce2gc = jce2 + ma%jbr3

      idi1sl = idi1gb - ma%ibb2
      idi2sl = idi2gb + ma%ibt2
      jdi1sl = jdi1gb - ma%jbl2
      jdi2sl = jdi2gb + ma%jbr2
      ide1sl = ide1gb - ma%ibb2
      ide2sl = ide2gb + ma%ibt2
      jde1sl = jde1gb - ma%jbl2
      jde2sl = jde2gb + ma%jbr2
      ici1sl = ici1gb - ma%ibb2
      ici2sl = ici2gb + ma%ibt2
      jci1sl = jci1gb - ma%jbl2
      jci2sl = jci2gb + ma%jbr2
      ice1sl = ice1gb - ma%ibb2
      ice2sl = ice2gb + ma%ibt2
      jce1sl = jce1gb - ma%jbl2
      jce2sl = jce2gb + ma%jbr2

      jde1sg = (jde1-1)*nsg+1
      jde2sg = jde2*nsg
      ide1sg = (ide1-1)*nsg+1
      ide2sg = ide2*nsg

#ifdef DEBUG
      write(ndebug,*) 'TOPLEFT     = ', ma%topleft
      write(ndebug,*) 'TOP         = ', ma%top
      write(ndebug,*) 'TOPRIGHT    = ', ma%topright
      write(ndebug,*) 'RIGHT       = ', ma%right
      write(ndebug,*) 'BOTTOMRIGHT = ', ma%bottomright
      write(ndebug,*) 'BOTTOM      = ', ma%bottom
      write(ndebug,*) 'BOTTOMLEFT  = ', ma%bottomleft
      write(ndebug,*) 'LEFT        = ', ma%left
      write(ndebug,*) 'DOTPEXTJI1 : ', jde1 , jde2 , ide1 , ide2
      write(ndebug,*) 'DOTPINTJI1 : ', jdi1 , jdi2 , idi1 , idi2
      write(ndebug,*) 'DOTPINTJI2 : ', jdii1 , jdii2 , idii1 , idii2
      write(ndebug,*) 'CRXPEXTJI1 : ', jce1 , jce2 , ice1 , ice2
      write(ndebug,*) 'CRXPINTJI1 : ', jci1 , jci2 , ici1 , ici2
      write(ndebug,*) 'CRXPINTJI2 : ', jcii1 , jcii2 , icii1 , icii2
      write(ndebug,*) 'TOPBDY   : ', ma%has_bdytop
      write(ndebug,*) 'BTMBDY   : ', ma%has_bdybottom
      write(ndebug,*) 'RGTBDY   : ', ma%has_bdyright
      write(ndebug,*) 'LFTBDY   : ', ma%has_bdyleft
      flush(ndebug)
#endif
    end subroutine setup_model_indexes

    subroutine allocate_domain(dom)
      implicit none
      type(domain) , intent(out) :: dom
      call getmem2d(dom%ht,jde1gb,jde2gb,ide1gb,ide2gb,'storage:ht')
      call getmem2d(dom%lndcat,jde1,jde2,ide1,ide2,'storage:lndcat')
      call getmem2d(dom%lndtex,jde1,jde2,ide1,ide2,'storage:lndtex')
      call getmem2d(dom%xlat,jde1ga,jde2ga,ide1ga,ide2ga,'storage:xlat')
      call getmem2d(dom%xlon,jde1ga,jde2ga,ide1ga,ide2ga,'storage:xlon')
      call getmem2d(dom%mask,jde1,jde2,ide1,ide2,'storage:mask')
      call getmem2d(dom%dlat,jde1,jde2,ide1,ide2,'storage:dlat')
      call getmem2d(dom%dlon,jde1,jde2,ide1,ide2,'storage:dlon')
      if ( idynamic == 3 ) then
        call getmem2d(dom%clv,jce1,jce2,ide1ga,ide2ga,'storage:clv')
        call getmem2d(dom%fmyu,jce1,jce2,ice1,ice2,'storage:fmyu')
        call getmem2d(dom%hx,jdi1ga,jdi2ga,ice1,ice2,'storage:hx')
        call getmem2d(dom%hy,jce1,jce2,idi1ga,idi2ga,'storage:hy')
      end if
      call getmem2d(dom%msfx,jd1,jd2,id1,id2,'storage:msfx')
      call getmem2d(dom%msfd,jd1,jd2,id1,id2,'storage:msfd')
      call getmem2d(dom%coriol,jde1,jde2,ide1,ide2,'storage:f')
      call getmem2d(dom%snowam,jde1,jde2,ide1,ide2,'storage:snowam')
      call getmem2d(dom%smoist,jde1,jde2,ide1,ide2,'storage:smoist')
      call getmem3d(dom%rmoist,jde1,jde2,ide1,ide2, &
                    1,num_soil_layers,'storage:rmoist')
      call getmem2d(dom%ldmsk,jci1,jci2,ici1,ici2,'storage:ldmsk')
      call getmem2d(dom%iveg,jci1,jci2,ici1,ici2,'storage:iveg')
      call getmem2d(dom%itex,jci1,jci2,ici1,ici2,'storage:itex')
      call getmem2d(dom%xmsf,jdi1,jdi2,idi1,idi2,'storage:xmsf')
      call getmem2d(dom%dmsf,jdi1,jdi2,idi1,idi2,'storage:dmsf')
      if ( lakemod == 1 ) then
        call getmem2d(dom%dhlake,jde1,jde2,ide1,ide2,'storage:dhlake')
      end if
      if ( idynamic == 2 ) then
        call getmem2d(dom%ef,jdi1ga,jdi2ga,idi1ga,idi2ga,'storage:ef')
        call getmem2d(dom%ddx,jdi1ga,jdi2ga,idi1ga,idi2ga,'storage:ddx')
        call getmem2d(dom%ddy,jdi1ga,jdi2ga,idi1ga,idi2ga,'storage:ddy')
        call getmem2d(dom%ex,jci1,jci2,ici1,ici2,'storage:ex')
        call getmem2d(dom%crx,jci1,jci2,ici1,ici2,'storage:crx')
        call getmem2d(dom%cry,jci1,jci2,ici1,ici2,'storage:cry')
        call getmem2d(dom%dmdy,jdi1,jdi2,idi1,idi2,'storage:dmdy')
        call getmem2d(dom%dmdx,jdi1,jdi2,idi1,idi2,'storage:dmdx')
      end if
    end subroutine allocate_domain

    subroutine allocate_domain_subgrid(sub)
      implicit none
      type(domain_subgrid) , intent(out) :: sub
      call getmem3d(sub%ht,1,nnsg,jde1,jde2,ide1,ide2,'storage:ht')
      call getmem3d(sub%lndcat,1,nnsg,jde1,jde2,ide1,ide2,'storage:lndcat')
      call getmem3d(sub%lndtex,1,nnsg,jde1,jde2,ide1,ide2,'storage:lndtex')
      call getmem3d(sub%xlat,1,nnsg,jde1,jde2,ide1,ide2,'storage:xlat')
      call getmem3d(sub%xlon,1,nnsg,jde1,jde2,ide1,ide2,'storage:xlon')
      call getmem3d(sub%mask,1,nnsg,jde1,jde2,ide1,ide2,'storage:xlon')
      call getmem3d(sub%ldmsk,1,nnsg,jci1,jci2,ici1,ici2,'storage:ldmsk')
      call getmem3d(sub%iveg,1,nnsg,jci1,jci2,ici1,ici2,'storage:iveg')
      call getmem3d(sub%itex,1,nnsg,jci1,jci2,ici1,ici2,'storage:itex')
      if ( lakemod == 1 ) then
        call getmem3d(sub%dhlake,1,nnsg,jde1,jde2,ide1,ide2,'storage:dhlake')
      end if
    end subroutine allocate_domain_subgrid

    subroutine allocate_surfstate(sfs)
      implicit none
      type(surfstate) , intent(out) :: sfs
      call getmem2d(sfs%psa,jce1ga,jce2ga,ice1ga,ice2ga,'surf:psa')
      call getmem2d(sfs%psdota,jde1ga,jde2ga,ide1ga,ide2ga,'surf:psdota')
      call getmem2d(sfs%psb,jx1,jx2,ix1,ix2,'surf:psb')
      call getmem2d(sfs%psdotb,jd1,jd2,id1,id2,'surf:psdotb')
      call getmem2d(sfs%psc,jce1,jce2,ice1,ice2,'surf:psc')
      call getmem2d(sfs%tg,jci1,jci2,ici1,ici2,'surf:tg')
      call getmem2d(sfs%hfx,jci1,jci2,ici1,ici2,'surf:hfx')
      call getmem2d(sfs%qfx,jci1,jci2,ici1,ici2,'surf:qfx')
      call getmem2d(sfs%rainc,jci1,jci2,ici1,ici2,'surf:rainc')
      call getmem2d(sfs%rainnc,jci1,jci2,ici1,ici2,'surf:rainnc')
      if ( ipptls > 1 ) then
        call getmem2d(sfs%snownc,jci1,jci2,ici1,ici2,'surf:snownc')
      end if
      call getmem2d(sfs%tgbb,jci1,jci2,ici1,ici2,'surf:tgbb')
      call getmem2d(sfs%uvdrag,jci1,jci2,ici1,ici2,'surf:uvdrag')
      call getmem2d(sfs%zo,jci1,jci2,ici1,ici2,'surf:zo')
      call getmem2d(sfs%ram1,jci1,jci2,ici1,ici2,'surf:ram1')
      call getmem2d(sfs%rah1,jci1,jci2,ici1,ici2,'surf:rah1')
      call getmem2d(sfs%br,jci1,jci2,ici1,ici2,'surf:br')
      call getmem2d(sfs%q2m,jci1,jci2,ici1,ici2,'surf:q2m')
      call getmem2d(sfs%ustar,jci1,jci2,ici1,ici2,'surf:ustar')
      call getmem2d(sfs%w10m,jci1,jci2,ici1,ici2,'surf:w10m')
      call getmem2d(sfs%u10m,jci1,jci2,ici1,ici2,'surf:u10m')
      call getmem2d(sfs%v10m,jci1,jci2,ici1,ici2,'surf:v10m')
      if ( ibltyp == 4 ) then
        call getmem2d(sfs%uz0,jci1,jci2,ici1,ici2,'surf:uz0')
        call getmem2d(sfs%vz0,jci1,jci2,ici1,ici2,'surf:vz0')
        call getmem2d(sfs%thz0,jci1,jci2,ici1,ici2,'surf:thz0')
        call getmem2d(sfs%qz0,jci1,jci2,ici1,ici2,'surf:qz0')
      end if
    end subroutine allocate_surfstate

    subroutine allocate_mod_atm_interface
      implicit none

      id1 = ide1ga
      id2 = ide2ga
      jd1 = jde1ga
      jd2 = jde2ga
      ix1 = ice1ga
      ix2 = ice2ga
      jx1 = jce1ga
      jx2 = jce2ga

      call allocate_domain(mddom)
      call allocate_domain_subgrid(mdsub)
      call allocate_surfstate(sfs)

    end subroutine allocate_mod_atm_interface

end module mod_atm_stub

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
