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

module mod_pbl_interface

  use m_realkinds
  use mod_service
  use mod_constants
  use mod_atm_interface , only : atmstate , diffx , slice , surfstate , &
                   surfpstate , surftstate , domain , uvcross2dot
  use mod_pbl_common
  use mod_pbl_holtbl
  use mod_pbl_uwtcm

  public

  contains

  subroutine init_pbl(atm2,atms,aten,holtten,uwten,adf,heatrt,chib,chiten,    &
                      remdrd,sps2,sts2,sfsta,mddom,ldmsk,a,sigma,dsigma,ptop, &
                      chtrdpv,chtrname,ichem,ichdrydepo,dt)
    implicit none
    integer , intent(in) :: ichem , ichdrydepo
    type (atmstate) , intent(in) :: atm2 , aten , holtten , uwten
    type (slice) , intent(in) :: atms
    type (diffx) , intent(in) :: adf
    type (domain) , intent(in) :: mddom
    type (surfstate) , intent(in) :: sfsta
    type (surftstate) , intent(in) :: sts2
    type (surfpstate) , intent(in) :: sps2
    real(dp) , pointer , dimension(:,:,:) :: heatrt
    real(dp) , pointer , dimension(:,:,:,:) :: chib
    real(dp) , pointer , dimension(:,:,:,:) :: chiten
    real(dp) , pointer , dimension(:,:,:) :: remdrd
    integer , pointer , dimension(:,:) :: ldmsk
    real(dp) , pointer , dimension(:) :: a
    real(dp) , pointer , dimension(:) :: sigma
    real(dp) , pointer , dimension(:) :: dsigma
    real(dp) :: dt , ptop
    real(dp) , pointer , dimension(:,:) :: chtrdpv
    character(len=5) , pointer , dimension(:) :: chtrname

    ptp = ptop
    if ( ichem == 1 )      lchem = .true.
    if ( ichdrydepo == 1 ) lchdrydepo = .true.

    itcmstart = 2
    itcmend = iy-1 
    if ( myid == 0 ) then
      jtcmstart = 2
    else
      jtcmstart = 1
    end if

    if ( myid == nproc-1 ) then
      jtcmend = jxp-1 
    else
      jtcmend = jxp 
    end if

    dttke = dt
    dtpbl = dt
    tkemin = 1.0D-8

    call assignpnt(aten%u,uten)
    call assignpnt(aten%v,vten)
    call assignpnt(aten%t,tten)
    call assignpnt(aten%tke,tketen)
    call assignpnt(aten%qv,qvten)
    call assignpnt(aten%qc,qcten)
    call assignpnt(uwten%u,uuwten)
    call assignpnt(uwten%v,vuwten)
    call assignpnt(uwten%t,tuwten)
    call assignpnt(uwten%tke,tkeuwten)
    call assignpnt(uwten%qv,qvuwten)
    call assignpnt(uwten%qc,qcuwten)
    call assignpnt(atms%ubx3d,uatm)
    call assignpnt(atms%vbx3d,vatm)
    call assignpnt(atms%ubd3d,udatm)
    call assignpnt(atms%vbd3d,vdatm)
    call assignpnt(atms%tb3d,tatm)
    call assignpnt(atms%qvb3d,qvatm)
    call assignpnt(atms%qcb3d,qcatm)
    call assignpnt(atm2%tke,tkeatm)
    call assignpnt(atms%thx3d,thxatm)
    call assignpnt(adf%difft,difft)
    call assignpnt(adf%diffq,diffq)
    call assignpnt(heatrt,radheatrt)
    call assignpnt(holtten%qv,diagqv)
    call assignpnt(holtten%qc,diagqc)
    call assignpnt(chib,chmx)
    call assignpnt(chiten,chten)
    call assignpnt(remdrd,drmr)
    call assignpnt(sps2%ps,sfcps)
    call assignpnt(sps2%pdot,sfcpd)
    call assignpnt(sts2%tg,tg)
    call assignpnt(sfsta%qfx,qfx)
    call assignpnt(sfsta%hfx,hfx)
    call assignpnt(sfsta%zpbl,zpbl)
    call assignpnt(sfsta%uvdrag,uvdrag)
    call assignpnt(mddom%coriol,coriolis)
    call assignpnt(mddom%msfx,mapfcx)
    call assignpnt(ldmsk,landmsk)
    call assignpnt(a,hlev)
    call assignpnt(sigma,flev)
    call assignpnt(dsigma,dlev)
    call assignpnt(chtrdpv,depvel)
    chname => chtrname

  end subroutine init_pbl

  subroutine get_data_from_tcm(tcmstate,tcmtend,aten,atm1,atm2,bRegridWinds)
    implicit none 
    type(atmstate) , intent(inout) :: tcmtend , aten , atm1 , atm2
    type(tcm_state) , intent(inout) :: tcmstate
    logical , intent(in) :: bRegridWinds

    ! Don't update the model variables if we are the diagnostic mode
    ! (Holtslag running, and UW updating tke)
    if ( ibltyp /= 99 ) then
      !
      ! Put the t and qv tendencies in to difft and diffq for
      ! application of the sponge boundary conditions (see mod_tendency)
      !
      diffq(itcmstart:itcmend,:,jtcmstart:jtcmend) =  &
                  diffq(itcmstart:itcmend,:,jtcmstart:jtcmend) +  &
                  tcmtend%qv(itcmstart:itcmend,:,jtcmstart:jtcmend)
      difft(itcmstart:itcmend,:,jtcmstart:jtcmend) =  &
                  difft(itcmstart:itcmend,:,jtcmstart:jtcmend) +  &
                  tcmtend%t(itcmstart:itcmend,:,jtcmstart:jtcmend)

      ! Put the cloud water tendency in aten
      aten%qc(itcmstart:itcmend,:,jtcmstart:jtcmend) =  &
                 aten%qc(itcmstart:itcmend,:,jtcmstart:jtcmend) +   &
                 tcmtend%qc(itcmstart:itcmend,:,jtcmstart:jtcmend)

      if ( bRegridWinds ) then
        !
        ! If the TCM calculations were all done on the cross grid, then
        ! the u and v tendencies need to be regridded to the dot grid
        !
        call uvcross2dot(tcmtend,aten)
      else
        ! Otherwise, simply set the appropriate variables in aten
        aten%u(itcmstart:itcmend,:,jtcmstart:jtcmend) = &
                 aten%u(itcmstart:itcmend,:,jtcmstart:jtcmend) +    &
                 tcmtend%u(itcmstart:itcmend,:,jtcmstart:jtcmend)
        aten%v(itcmstart:itcmend,:,jtcmstart:jtcmend) = &
                 aten%v(itcmstart:itcmend,:,jtcmstart:jtcmend) +    &
                 tcmtend%v(itcmstart:itcmend,:,jtcmstart:jtcmend)
      end if

      zpbl(itcmstart:itcmend,jtcmstart:jtcmend) = &
                tcmstate%zpbl(itcmstart:itcmend,jtcmstart:jtcmend)
    end if

!   !
!   ! Interpolate kzm and kth from the interfaces to the midpoints
!   ! if the diffusivities are to be used in the holtslag model
!   !
!   if ( ibltyp == 99 ) then
!     do k = 1 , kz
!       tcmstate%kzm(:,k,:) = sqrt(tcmstate%kzm(:,k,:)*tcmstate%kzm(:,k+1,:))
!       tcmstate%kth(:,k,:) = sqrt(tcmstate%kth(:,k,:)*tcmstate%kth(:,k+1,:))
!     end do
!   end if
!   !
!   ! Shift kth and kzm for output if the UW model is running
!   !
!   if ( ibltyp == 2 ) then
!     do k = 1 , kz
!       tcmstate%kzm(:,k,:) = tcmstate%kzm(:,k+1,:)
!       tcmstate%kth(:,k,:) = tcmstate%kth(:,k+1,:)
!     end do
!   end if

    aten%tke(itcmstart:itcmend,:,jtcmstart:jtcmend) = &
                tcmtend%tke(itcmstart:itcmend,:,jtcmstart:jtcmend)
    !
    ! Set the surface tke (diagnosed)
    !
    atm1%tke(itcmstart:itcmend,kzp1,jtcmstart:jtcmend) =  &
               tcmstate%srftke(itcmstart:itcmend,jtcmstart:jtcmend)
    atm2%tke(itcmstart:itcmend,kzp1,jtcmstart:jtcmend) =  &
               tcmstate%srftke(itcmstart:itcmend,jtcmstart:jtcmend)

!   tcmtend%qc = 0.0d0
!   tcmtend%qv = 0.0d0
!   tcmtend%t = 0.0d0
!   tcmtend%u = 0.0d0
!   tcmtend%v = 0.0d0
!   tcmtend%tke = 0.0d0


  end subroutine get_data_from_tcm

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     This subroutine computes the horizontal flux-divergence terms   c
!     for tke.  Second-order difference is used.                      c
!                                                                     c
!     dxx    : is the horizontal distance.                            c
!     j      : is the j'th slice of f anf ften.                       c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine hadvtke(tcmstate,atm,twt,dxx,j)
    implicit none
    integer , intent(in) :: j
    real(dp) , intent(in) :: dxx
    type(atmstate) , intent(in) :: atm
    type(tcm_state) , intent(inout) :: tcmstate
    real(dp) , pointer , dimension(:,:) :: twt
    character (len=64) :: subroutine_name='hadvtke'
    integer :: idindx = 0
    integer :: i , k , jm1 , jp1
!
    call time_begin(subroutine_name,idindx)

    jm1 = j - 1
    jp1 = j + 1

    do k = 2 , kz
      do i = 2 , iym1

        !
        ! Interpoalte the winds to the full sigma levels
        ! while the advection term is calculated
        !
        tcmstate%advtke(i,k,j)= tcmstate%advtke(i,k,j)  &
           -(((atm%u(i+1,k-1,jp1)+atm%u(i,k-1,jp1))*twt(k,2)  &
             +(atm%u(i+1,k,jp1)  +atm%u(i,k,jp1))*twt(k,1)) &
             *( atm%tke(i,k,j)+atm%tke(i,k,jp1))  &
            -((atm%u(i+1,k-1,j)+atm%u(i,k-1,j))*twt(k,2)  &
             +(atm%u(i+1,k,j)  +atm%u(i,k,j))*twt(k,1)) &
             *( atm%tke(i,k,j)+atm%tke(i,k,jm1))  &
            +((atm%v(i+1,k-1,j)+atm%v(i+1,k-1,jp1))*twt(k,2)  &
             +(atm%v(i+1,k,j)  +atm%v(i+1,k,jp1))*twt(k,1)) &
             *( atm%tke(i,k,j)+atm%tke(i+1,k,j))  &
            -((atm%v(i,k-1,j)+atm%v(i,k-1,jp1))*twt(k,2)  &
             +(atm%v(i,k,j)  +atm%v(i,k,jp1))*twt(k,1)) &
             *( atm%tke(i,k,j)+atm%tke(i-1,k,j))) &
             /(dxx*mapfcx(i,j)*mapfcx(i,j))

!TAO Debug:
!       tcmstate%advtke(i,k,j)= tcmstate%advtke(i,k,j)  + 1d-8*atm%tke(i,k,jp1)

      end do
    end do

    call time_end(subroutine_name,idindx)

  end subroutine hadvtke

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the vertical flux-divergence terms.    c
!                                                                     c
!     j      : jth slice of variable fa.                              c
!                                                                     c
!     ind = 1 : Bretherton's vertical advection method                c
!           2 : Alternate vertical advection method (unknown origin)  c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine vadvtke(tcmstate,qdot,j,ind)
    implicit none
    integer , intent(in) :: j , ind
    type(tcm_state) :: tcmstate
    real(dp) , dimension(:,:,:) , intent(in) :: qdot
    real(dp) , dimension(iy,kz) :: dotqdot , ftmp
    character (len=64) :: subroutine_name='vadvtke'
    integer :: idindx = 0
!
    integer :: i , k

    call time_begin(subroutine_name,idindx)
        
    !
    ! Use Bretherton's method for tke advection
    !
    if ( ind == 1 ) then
      do k = 1 , kz
        do i = 2 , iym1
          dotqdot(i,k) = (qdot(i,k,j) + qdot(i,k+1,j))
          ftmp(i,k) = 0.5D0*(tcmstate%tkeps(i,k,j) + tcmstate%tkeps(i,k+1,j))
        end do
      end do

      do k = 2 , kz
        do i = 2 , iym1
          tcmstate%advtke(i,k,j) = tcmstate%advtke(i,k,j) - &
                                    (dotqdot(i,k)*ftmp(i,k) -  &
                                    dotqdot(i,k-1)*ftmp(i,k-1))/ &
                                    (dlev(k)+dlev(k-1))
        end do
      end do
    !
    ! Use an alternative method (this came from where?)
    !
    else
      do i = 2 , iym1
        tcmstate%advtke(i,1,j)=tcmstate%advtke(i,1,j)-  &
                              qdot(i,2,j)*tcmstate%tkeps(i,2,j)/dlev(1)
      end do
      do k = 2 , kzm1
        do i = 2 , iym1
          tcmstate%advtke(i,k,j)=tcmstate%advtke(i,k,j) &
                              -(qdot(i,k+1,j)*tcmstate%tkeps(i,k+1,j)   &
                                -qdot(i,k,j)*tcmstate%tkeps(i,k,j))/dlev(k)
        end do
      end do
      do i = 2 , iym1
        tcmstate%advtke(i,kz,j)=tcmstate%advtke(i,kz,j)+  &
                              qdot(i,kz,j)*tcmstate%tkeps(i,kz,j)/dlev(kz)
      end do
    end if

    call time_end(subroutine_name,idindx)

  end subroutine vadvtke

  subroutine set_tke_bc(atm1,atm2)
    implicit none
    type(atmstate) , intent(inout) :: atm1 , atm2
    !
    ! Set tke boundary conditions
    !
    if ( myid == 0 ) then
      atm1%tke(:,:,0:1) = tkemin ! East boundary
      atm2%tke(:,:,-1:1) = tkemin ! East boundary
    end if
    if ( myid == nproc-1 ) then
      atm1%tke(:,:,jxp:jxp+1) = tkemin ! West boundary
      atm2%tke(:,:,jxp:jxp+2) = tkemin ! West boundary
    end if

    atm1%tke(1,:,:) = tkemin  !North boundary
    atm2%tke(1,:,:) = tkemin  !North boundary
    atm1%tke(iy,:,:) = tkemin  !South boundary
    atm2%tke(iy,:,:) = tkemin  !South boundary
    !
    ! End set tke boundary conditions
    !
  end subroutine set_tke_bc

  subroutine check_conserve_qt(rcmqvten,rcmqcten,tcmtend,tcmstate,kmax)
    implicit none
    type(atmstate) , intent(in) :: tcmtend
    type(tcm_state) :: tcmstate
    real(dp) , dimension(:,:,:) , intent(in) :: rcmqvten , rcmqcten
    integer , intent(in) :: kmax
    real(dp) , dimension(kmax) :: rho1d , rhobydpdz1d
    real(dp) :: qwtcm , qwrcm , qwanom , dtops , xps , ps2 , dza
    integer :: i , j , k

    do i = 2 , iym1
      do j = jbegin , jendx

        do k = 1 , kzm1
          xps = (hlev(k)*sfcps(i,j)+ptp)
          ps2 = (hlev(k+1)*sfcps(i,j)+ptp)
          dza = za(i,k,j) - za(i,k+1,j)
          rhobydpdz1d(k) = d_1000*(ps2-xps)/(egrav*dza)

        end do
        rhobydpdz1d(kz) = rhox2d(i,j)
        dtops = dtpbl/sfcps(i,j)

        rho1d = d_1000*(hlev*sfcps(i,j) + ptp) / &
                    ( rgas * tatm(i,:,j) *  &
                    (d_one + ep1* qvatm(i,:,j) - qcatm(i,:,j))  )


        qwtcm = sum((tcmtend%qv(i,:,j) + tcmtend%qc(i,:,j))  &
                      *rho1d*dzq(i,:,j))*dtops
        qwrcm = sum((rcmqvten(i,:,j) + rcmqcten(i,:,j))   &
                      *rho1d*dzq(i,:,j))*dtops
!                     *rhobydpdz1d*dzq(i,:,j))*dtops

!       qwanom = qwtcm - qwrcm
!       qwanom = qwrcm-qfx(i,j)*dt
        qwanom = qwtcm - dtpbl*qfx(i,j)
        tcmstate%kzm(i,1,j) = qwanom

      end do
    end do
  end subroutine check_conserve_qt

end module mod_pbl_interface
