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
 
module mod_cu_em
!
! Kerry Emanuel Convective scheme
!
  use mod_dynparam
  use mod_memutil
  use mod_cu_common
!
  private
!
  public :: allocate_mod_cu_em , cupemandrv
  public :: minsig , elcrit , tlcrit , entp , sigd , sigs ,    &
            omtrain , omtsnow , coeffr , coeffs , cu , betae , &
            dtmax , alphae , damp , minorig
!
  real(dp) :: alphae , betae , coeffr , coeffs , cu , damp , dtmax , &
             elcrit , entp , minsig , omtrain , omtsnow , sigd ,    &
             sigs , tlcrit
!
  integer :: minorig
  real(dp) , parameter :: cl = 2500.0D0
  real(dp) , parameter :: cpvmcl = cl - cpv
  real(dp) , parameter :: mincbmf = 1.0D-30
!
  real(dp) , public , pointer , dimension(:,:) :: cbmf2d
!
  contains
!
  subroutine allocate_mod_cu_em
    implicit none
    call getmem2d(cbmf2d,jce1,jce2,ice1,ice2,'mod_cu_em:cbmf2d')
  end subroutine allocate_mod_cu_em
!
!
! **********************************************
! **** Driver for Emanuel Convection Scheme ****
! **********************************************
!
  subroutine cupemandrv(jstart,jend,istart,iend,ktau)
! 
    implicit none
!
    integer , intent(in) :: jstart , jend , istart , iend
    integer(8) , intent(in) :: ktau
!
    integer , parameter :: ntra = 0
!
    real(dp) :: akclth , cbmf , pret , qprime , tprime , wd , prainx
    real(dp) , dimension(kz) :: fq , ft , fu , fv , pcup , qcup ,      &
                               qscup , tcup , ucup , vcup
    real(dp) , dimension(kz,1) :: ftra , tra
    integer :: i , j , k , iflag , kbase , kclth , kk , ktop
    real(dp) , dimension(kzp1) :: phcup
!
    total_precip_points = 0
    do i = istart , iend
      do j = jstart , jend
        if ( icup /= 4 ) then
          if ( cucontrol(j,i) /= 4 ) cycle
        end if
        do k = 1 , kz
          kk = kzp1 - k
          tcup(k) = tas(j,i,kk)                         ! [k]
          qcup(k) = qvas(j,i,kk)/(d_one+qvas(j,i,kk))   ! [kg/kg]
          qscup(k) = qsas(j,i,kk)/(d_one+qsas(j,i,kk))  ! [kg/kg]
          ucup(k) = uas(j,i,kk)                         ! [m/s]
          vcup(k) = vas(j,i,kk)                         ! [m/s]
          pcup(k) = pas(j,i,kk)*d_10                    ! [hPa]
        end do
        do k = 1 , kzp1
          kk = kzp1 - k + 1
          phcup(k) = (flev(kk)*sfcps(j,i)+ptop)*d_10 ! [hPa]
        end do
        cbmf = cbmf2d(j,i)                              ! [(kg/m**2)/s]
   
        call cupeman(tcup,qcup,qscup,ucup,vcup,tra,pcup,phcup,kz,kzp1,  &
                     kzm1,ntra,iflag,ft,fq,fu,fv,ftra,pret,wd,    &
                     tprime,qprime,cbmf,kbase,ktop)
   
        cbmf2d(j,i) = cbmf
   
        ! iflag=0: No moist convection; atmosphere stable or surface
        !          temperature < 250K or surface humidity is negative.
        ! iflag=1: Moist convection occurs.
        ! iflag=2: No moist convection: lifted condensation level above 200 mb.
        ! iflag=3: No moist convection: cloud base higher than level kzm2.
        ! iflag=4: Moist convection occurs, but CFL condition on the
        !          subsidence warming is violated. (Does not terminate scheme.)
        if ( iflag == 1 .or. iflag == 4 ) then ! If moist convection
!         if ( iflag == 4 ) then               ! If CFL violation
!           print*,'EMAN CFL VIOLATION: ',i,j,cbmf
!         end if
   
          ! Tendencies
          do k = 1 , kz
            kk = kzp1 - k
            tten(j,i,kk) = ft(k)*sfcps(j,i) + tten(j,i,kk)
            qvten(j,i,kk) = fq(k)/(d_one-fq(k))* &
                              sfcps(j,i)+qvten(j,i,kk)
            ! There is a bit of an inconsistency here...  The wind
            ! tendencies from convection are on cross points, but the
            ! model wants them on dot points.
            uten(j,i,kk) = fu(k)*sfcps(j,i) + uten(j,i,kk)
            vten(j,i,kk) = fv(k)*sfcps(j,i) + vten(j,i,kk)
          end do
   
          ! Cloud fraction and cloud water
          kclth = ktop - kbase + 1
          akclth = d_one/dble(kclth)
          do k = kbase , ktop
            kk = kzp1 - k
            rcldlwc(j,i,kk) = cllwcv
            rcldfra(j,i,kk) = d_one - (d_one-clfrcv)**akclth
          end do
   
          ! Precipitation
          prainx = pret*dtmdl
          if ( prainx > dlowval ) then
            rainc(j,i)  = rainc(j,i)  + prainx  ! mm
            if ( ktau == 0 .and. debug_level > 2 ) then
              lmpcpc(j,i) = lmpcpc(j,i) + pret
            else
              lmpcpc(j,i) = lmpcpc(j,i) + pret*aprdiv
            end if
            total_precip_points = total_precip_points + 1
          end if
        end if
      end do
    end do
! 
  end subroutine cupemandrv
!
!**************************************************************************
!****                       subroutine cupeman (formerly convect)     *****
!****                          version 4.3C                           *****
!****                          20 May, 2002                           *****
!****                          Kerry Emanuel                          *****
!**************************************************************************
!
!-----------------------------------------------------------------------------
!    *** on input:      ***
!
!     t:   array of absolute temperature (k) of dimension nd, with first
!           index corresponding to lowest model level. note that this array
!           will be altered by the subroutine if dry convective adjustment
!           occurs and if ipbl is not equal to 0.
!
!     q:   array of specific humidity (gm/gm) of dimension nd, with first
!            index corresponding to lowest model level. must be defined
!            at same grid levels as t. note that this array will be altered
!            if dry convective adjustment occurs and if ipbl is not equal to 0.
!
!     qs:  array of saturation specific humidity of dimension nd, with first
!            index corresponding to lowest model level. must be defined
!            at same grid levels as t. note that this array will be altered
!            if dry convective adjustment occurs and if ipbl is not equal to 0.
!
!     u:   array of zonal wind velocity (m/s) of dimension nd, witth first
!            index corresponding with the lowest model level. defined at
!            same levels as t. note that this array will be altered if
!            dry convective adjustment occurs and if ipbl is not equal to 0.
!
!     v:   same as u but for meridional velocity.
!
!     tra: array of passive tracer mixing ratio, of dimensions (nd,ntra),
!            where ntra is the number of different tracers. if no
!            convective tracer transport is needed, define a dummy
!            input array of dimension (nd,1). tracers are defined at
!            same vertical levels as t. note that this array will be altered
!            if dry convective adjustment occurs and if ipbl is not equal to 0.
!
!     p:   array of pressure (mb) of dimension nd, with first
!            index corresponding to lowest model level. must be defined
!            at same grid levels as t.
!
!     ph:  array of pressure (mb) of dimension nd+1, with first index
!            corresponding to lowest level. these pressures are defined at
!            levels intermediate between those of p, t, q and qs. the first
!            value of ph should be greater than (i.e. at a lower level than)
!            the first value of the array p.
!
!     nd:  the dimension of the arrays t,q,qs,p,ph,ft and fq
!
!     nl:  the maximum number of levels to which convection can
!            penetrate, plus 1.
!            nl must be less than or equal to nd-1.
!
!     ntra:the number of different tracers. if no tracer transport
!            is needed, set this equal to 1. (on most compilers, setting
!            ntra to 0 will bypass tracer calculation, saving some cpu.)
!
!----------------------------------------------------------------------------
!    ***   on output:         ***
!
!     iflag: an output integer whose value denotes the following:
!
!                value                        interpretation
!                -----                        --------------
!                  0               no moist convection; atmosphere is not
!                                  unstable, or surface temperature is less
!                                  than 250 k or surface specific humidity
!                                  is non-positive.
!
!                  1               moist convection occurs.
!
!                  2               no moist convection: lifted condensation
!                                  level is above the 200 mb level.
!
!                  3               no moist convection: cloud base is higher
!                                  then the level nl-1.
!
!                  4               moist convection occurs, but a cfl condition
!                                  on the subsidence warming is violated. this
!                                  does not cause mod_the scheme to terminate.
!
!     ft:   array of temperature tendency (k/s) of dimension nd, defined at same
!             grid levels as t, q, qs and p.
!
!     fq:   array of specific humidity tendencies ((gm/gm)/s) of dimension nd,
!             defined at same grid levels as t, q, qs and p.
!
!     fu:   array of forcing of zonal velocity (m/s^2) of dimension nd,
!             defined at same grid levels as t.
!
!     fv:   same as fu, but for forcing of meridional velocity.
!
!     ftra: array of forcing of tracer content, in tracer mixing ratio per
!             second, defined at same levels as t. dimensioned (nd,ntra).
!
!     precip: scalar convective precipitation rate (mm/s).
!
!     wd:    a convective downdraft velocity scale. for use in surface
!             flux parameterizations. see convect.ps file for details.
!
!     tprime: a convective downdraft temperature perturbation scale (k).
!              for use in surface flux parameterizations. see convect.ps
!              file for details.
!
!     qprime: a convective downdraft specific humidity
!              perturbation scale (gm/gm).
!              for use in surface flux parameterizations. see convect.ps
!              file for details.
!
!     cbmf:   the cloud base mass flux ((kg/m**2)/s). this scalar value must
!              be stored by the calling program and returned to convect at
!              its next call. that is, the value of cbmf must be "remembered"
!              by the calling program between calls to convect.
!
!------------------------------------------------------------------------------
!
!    ***  the parameter na should in general be greater than   ***
!    ***                or equal to  nd + 1                    ***
!
!------------------------------------------------------------------------------
! modifications for regcm:
!   1. units for precipitation were change from mm/day to mm/s
!   2. the thermodynamic constants were made consistent with those
!      of regcm.
!   3. dependence of latent heat of vaporization on tempertature
!      removed because it is neglected in regcm. this is done by
!      setting cpv equal to cl and setting wlhv to the regcm value.
!   4. added cloud base (icb) and cloud top (inb) to the output to
!      compute the cloud fraction and cloud liquid water content.
!   5. each variable is now explicitly declared.  that is, the
!      "implicit none" option was added.
!   6. the value minorig is increased because mod_the thickness of the
!      lowest layer(s) is(are) too small. if the thickness is too
!      small for the given timestep, the mass of the layer is likely
!      to be evacuated.
!   7. a maximum value to the cloud base mass flux has been added.
!
  subroutine cupeman(t,q,qs,u,v,tra,p,ph,nd,na,nl,ntra,iflag,ft,fq,fu,fv, &
                     ftra,precip,wd,tprime,qprime,cbmf,icb,inb)
!
    implicit none
!
    real(dp) :: cbmf , precip , qprime , tprime , wd
    integer :: icb , iflag , inb , na , nd , nl , ntra
    real(dp) , dimension(nd) :: fq , ft , fu , fv , p , ph , q , qs ,  &
                               t , u , v
    real(dp) , dimension(nd,1) :: ftra , tra
    intent (in) na , ntra , ph , p , nd , nl
    intent (out) tprime , wd
    intent (inout) cbmf , fq , ft , ftra , fu , fv , icb , iflag ,    &
                   inb , precip , q , qprime , qs , t , tra , u , v
!
    real(dp) :: a2 , ad , afac , ahm , ahmax , ahmin , alt , altem ,   &
               alv , alvnew , am , amp1 , anum , asij , asum , awat , &
               b6 , bf2 , bsum , by , byp , c6 , cape , capem ,       &
               cbmfold , chi , coeff , cpinv , cwat , damps , dbo ,   &
               dbosum , defrac , dei , delm , delp , delt0 , delti ,  &
               denom , dhdp , dphinv , dpinv , dtma , dtmnx , dtpbl , &
               elacrit , ents , epmax , fac , fqold , frac , ftold ,  &
               ftraold , fuold , fvold , plcl , qnew , qp1 , qsm ,    &
               qstm , qti , rat , rdcp , revap , rh , rm , scrit , sigt
    real(dp) , dimension(na) :: clw , cpn , ep , evap , gz , h , hm ,  &
                               hp , lv , lvcp , m , mp , qp , sigp ,  &
                               th , told , tp , tratm , tv , tvp ,    &
                               up , vp , water , wt
    real(dp) , dimension(na,na) :: elij , ment , qent , sij , uent , vent
    integer :: i , ihmin , inb1 , ipbl , j , jc , jn , jtt , k , nk
    integer , dimension(na) :: nent
    real(dp) :: sjmax , sjmin , smid , smin , stemp , tc , tca ,       &
               thbar , tnew , traav , tvaplcl , tvpplcl , tvx , tvy , &
               uav , um , vav , vm , wdtrain , x
    real(dp) , dimension(na,na,ntra) :: traent
    real(dp) , dimension(na,ntra) :: trap
!
!   specify switches                        
!
!   ipbl: set to zero to bypass dry adiabatic adjustment      
!         any other value results in dry adiabatic adjustment   
!         (zero value recommended for use in models with boundary layer schemes)
!
!   minorig: lowest level from which convection may originate 
!            (should be first model level at which t is defined   
!            for models using bulk pbl schemes; otherwise, it should
!            be the first model level at which t is defined above
!            the surface layer) 
!
    ipbl = 0
!
!   assign values of thermodynamic constants, gravity, and liquid 
!   water density.  these should be consistent with those used in
!   calling program
!   note: these are also specified in subroutine tlift
!
    delti = d_one/dtcum
!
!   initialize output arrays and parameters
!
    do i = 1 , nd
      ft(i) = d_zero
      fq(i) = d_zero
      fu(i) = d_zero
      fv(i) = d_zero
      do j = 1 , ntra
        ftra(i,j) = d_zero
      end do
    end do
    do i = 1 , nl + 1
      rdcp = (rgas*(d_one-q(i))+q(i)*rwat)/(cpd*(d_one-q(i))+q(i)*cpv)
      th(i) = t(i)*(d_1000/p(i))**rdcp
    end do
    precip = d_zero
    wd = d_zero
    tprime = d_zero
    qprime = d_zero
    iflag = 0
!
    if ( ipbl /= 0 ) then
!
!     perform dry adiabatic adjustment
!
      jc = 0
      do i = nl - 1 , 1 , -1
        jn = 0
        asum = th(i)*(d_one+q(i)*rgowi-q(i))
        do j = i + 1 , nl
          asum = asum + th(j)*(d_one+q(j)*rgowi-q(j))
          thbar = asum/dble(j+1-i)
          if ( (th(j)*(d_one+q(j)*rgowi-q(j))) < thbar ) jn = j
        end do
        if ( i == 1 ) jn = max0(jn,2)
        if ( jn /= 0 ) then
          do
            ahm = d_zero
            rm = d_zero
            um = d_zero
            vm = d_zero
            do k = 1 , ntra
              tratm(k) = d_zero
            end do
            do j = i , jn
              ahm = ahm + (cpd*(d_one-q(j))+q(j)*cpv)*t(j)*(ph(j)-ph(j+1))
              rm = rm + q(j)*(ph(j)-ph(j+1))
              um = um + u(j)*(ph(j)-ph(j+1))
              vm = vm + v(j)*(ph(j)-ph(j+1))
              do k = 1 , ntra
                tratm(k) = tratm(k) + tra(j,k)*(ph(j)-ph(j+1))
              end do
            end do
            dphinv = d_one/(ph(i)-ph(jn+1))
            rm = rm*dphinv
            um = um*dphinv
            vm = vm*dphinv
            do k = 1 , ntra
              tratm(k) = tratm(k)*dphinv
            end do
            a2 = d_zero
            do j = i , jn
              q(j) = rm
              u(j) = um
              v(j) = vm
              do k = 1 , ntra
                tra(j,k) = tratm(k)
              end do
              rdcp = (rgas*(d_one-q(j))+q(j)*rwat)/(cpd*(d_one-q(j))+q(j)*cpv)
              x = (0.001D0*p(j))**rdcp
              told(j) = t(j)
              t(j) = x
              a2 = a2 + (cpd*(d_one-q(j))+q(j)*cpv)*x*(ph(j)-ph(j+1))
            end do
            do j = i , jn
              th(j) = ahm/a2
              t(j) = t(j)*th(j)
              tc = told(j) - tzero
              alv = wlhv - cpvmcl*tc
              qs(j) = qs(j) + qs(j)*(d_one+qs(j)*(rgowi-d_one)) * &
                      alv*(t(j)-told(j))/(rwat*told(j)*told(j))
            end do
            if ( ((th(jn+1)*(d_one+q(jn+1)*rgowi-q(jn+1))) <  &
                  (th(jn)*(d_one+q(jn)*rgowi-q(jn)))) ) then
              jn = jn + 1
              cycle
            end if
            if ( i == 1 ) jc = jn
            exit
          end do
        end if
      end do
!
!     remove any supersaturation that results from adjustment
!
      if ( jc > 1 ) then
        do j = 1 , jc
          if ( qs(j) < q(j) ) then
            alv = wlhv - cpvmcl*(t(j)-tzero)
            tnew = t(j) + alv*(q(j)-qs(j))/(cpd*(d_one-q(j))+cl*q(j)+qs(j) * &
                    (cpv-cl+alv*alv/(rwat*t(j)*t(j))))
            alvnew = wlhv - cpvmcl*(tnew-tzero)
            qnew = (alv*q(j)-(tnew-t(j))*(cpd*(d_one-q(j))+cl*q(j)))/alvnew
!rcm        precip = precip+24.*3600.*1.0e5*(ph(j)-ph(j+1))*  ! mm/d
            precip = precip + 1.0D5*(ph(j)-ph(j+1)) * &
                       (q(j)-qnew)*regrav/(dtcum*d_1000)  ! mm/s
            t(j) = tnew
            q(j) = qnew
            qs(j) = qnew
          end if
        end do
      end if
!
    end if
!
!   calculate arrays of geopotential, heat capacity and static energy
!
    gz(1) = d_zero
    cpn(1) = cpd*(d_one-q(1)) + q(1)*cpv
    h(1) = t(1)*cpn(1)
    lv(1) = wlhv - cpvmcl*(t(1)-tzero)
    hm(1) = lv(1)*q(1)
    tv(1) = t(1)*(d_one+q(1)*rgowi-q(1))
    ahmin = 1.0D12
    ihmin = nl
    do i = 2 , nl + 1
      tvx = t(i)*(d_one+q(i)*rgowi-q(i))
      tvy = t(i-1)*(d_one+q(i-1)*rgowi-q(i-1))
      gz(i) = gz(i-1) + (rgas*d_half)*(tvx+tvy)*(p(i-1)-p(i))/ph(i)
      cpn(i) = cpd*(d_one-q(i)) + cpv*q(i)
      h(i) = t(i)*cpn(i) + gz(i)
      lv(i) = wlhv - cpvmcl*(t(i)-tzero)
      hm(i) = (cpd*(d_one-q(i))+cl*q(i))*(t(i)-t(1))+lv(i)*q(i)+gz(i)
      tv(i) = t(i)*(d_one+q(i)*rgowi-q(i))
!
!     find level of minimum moist static energy
!
      if ( i >= minorig .and. hm(i) < ahmin .and. hm(i) < hm(i-1) ) then
        ahmin = hm(i)
        ihmin = i
      end if
    end do
    ihmin = min0(ihmin,nl-1)
!
!   find that model level below the level of minimum moist
!   static energy that has the maximum value of moist static
!   energy
!
    ahmax = d_zero
    nk = nl
    do i = minorig , ihmin
      if ( hm(i) > ahmax ) then
        nk = i
        ahmax = hm(i)
      end if
    end do
!
!   check whether parcel level temperature and specific humidity are reasonable
!   skip convection if hm increases monotonically upward
!
    if ( t(nk) < 250.0D0 .or. q(nk) <= d_zero .or. &
         ihmin == (nl-1) ) then
      iflag = 0
      cbmf = d_zero
      return
    end if
!
!   calculate lifted condensation level of air at parcel origin level
!   (within 0.2% of formula of bolton, mon. wea. rev.,1980)
!
    rh = q(nk)/qs(nk)
    chi = t(nk)/(1669.0D0-122.0D0*rh-t(nk))
    plcl = p(nk)*(rh**chi)
    if ( plcl < 200.0D0 .or. plcl >= 2000.0D0 ) then
      iflag = 2
      cbmf = d_zero
      return
    end if
!
!   calculate first level above lcl (=icb)
!
    icb = nl - 1
    do i = nk + 1 , nl
      if ( p(i) < plcl ) icb = min0(icb,i)
    end do
    if ( icb >= (nl-1) ) then
      iflag = 3
      cbmf = d_zero
      return
    end if
!
!   find temperature up through icb and test for instability     
!
!   subroutine tlift calculates part of the lifted parcel virtual
!   temperature, the actual temperature and the adiabatic
!   liquid water content
!
    call tlift(p,t,q,qs,gz,icb,nk,tvp,tp,clw,nd,nl,1)
    do i = nk , icb
      tvp(i) = tvp(i) - tp(i)*q(nk)
    end do
!
!   if there was no convection at last time step and parcel
!   is stable at icb then skip rest of calculation     
!
    if ( dabs(cbmf) < mincbmf .and. tvp(icb) <= (tv(icb)-dtmax) ) then
      iflag = 0
      return
    end if
!
!   if this point is reached, moist convective adjustment is necessary
!
    if ( iflag /= 4 ) iflag = 1
!
!   find the rest of the lifted parcel temperatures
!
    call tlift(p,t,q,qs,gz,icb,nk,tvp,tp,clw,nd,nl,2)
!
!   set the precipitation efficiencies and the fraction of
!   precipitation falling outside of cloud
!   these may be functions of tp(i), p(i) and clw(i)
!
    do i = 1 , nk
      ep(i) = d_zero
      sigp(i) = sigs
    end do
    do i = nk + 1 , nl
      tca = tp(i) - tzero
      if ( tca >= d_zero ) then
        elacrit = elcrit
      else
        elacrit = elcrit*(d_one-tca/tlcrit)
      end if
      elacrit = dmax1(elacrit,d_zero)
      epmax = 0.999D0
      ep(i) = epmax*(d_one-elacrit/dmax1(clw(i),1.0D-8))
      ep(i) = dmax1(ep(i),d_zero)
      ep(i) = dmin1(ep(i),epmax)
      sigp(i) = sigs
    end do
!
!   calculate virtual temperature and lifted parcel
!   virtual temperature
!
    do i = icb + 1 , nl
      tvp(i) = tvp(i) - tp(i)*q(nk)
    end do
    tvp(nl+1) = tvp(nl) - (gz(nl+1)-gz(nl))*rcpd
!
!   now initialize various arrays used in the computations
!
    do i = 1 , nl + 1
      hp(i) = h(i)
      nent(i) = 0
      water(i) = d_zero
      evap(i) = d_zero
      wt(i) = omtsnow
      mp(i) = d_zero
      m(i) = d_zero
      lvcp(i) = lv(i)/cpn(i)
      do j = 1 , nl + 1
        qent(i,j) = q(j)
        elij(i,j) = d_zero
        ment(i,j) = d_zero
        sij(i,j) = d_zero
        uent(i,j) = u(j)
        vent(i,j) = v(j)
        do k = 1 , ntra
          traent(i,j,k) = tra(j,k)
        end do
      end do
    end do
    qp(1) = q(1)
    up(1) = u(1)
    vp(1) = v(1)
    do i = 1 , ntra
      trap(1,i) = tra(1,i)
    end do
    do i = 2 , nl + 1
      qp(i) = q(i-1)
      up(i) = u(i-1)
      vp(i) = v(i-1)
      do j = 1 , ntra
        trap(i,j) = tra(i-1,j)
      end do
    end do
!
!   find the first model level (inb1) above the parcel's
!   highest level of neutral buoyancy
!   and the highest level of positive cape (inb)
!
    cape = d_zero
    capem = d_zero
    inb = icb + 1
    inb1 = inb
    byp = d_zero
    do i = icb + 1 , nl - 1
      by = (tvp(i)-tv(i))*(ph(i)-ph(i+1))/p(i)
      cape = cape + by
      if ( by >= d_zero ) inb1 = i + 1
      if ( cape > d_zero ) then
        inb = i + 1
        byp = (tvp(i+1)-tv(i+1))*(ph(i+1)-ph(i+2))/p(i+1)
        capem = cape
      end if
    end do
    inb = max0(inb,inb1)
    cape = capem + byp
    defrac = capem - cape
    defrac = dmax1(defrac,0.001D0)
    frac = -cape/defrac
    frac = dmin1(frac,d_one)
    frac = dmax1(frac,d_zero)
!
!   calculate liquid water static energy of lifted parcel
!
    do i = icb , inb
      hp(i) = h(nk) + (lv(i)+(cpd-cpv)*t(i))*ep(i)*clw(i)
    end do
!
!   calculate cloud base mass flux and rates of mixing, m(i),
!   at each model level                    
!
    dbosum = d_zero
!
!   interpolate difference between lifted parcel and
!   environmental temperatures to lifted condensation level
!
    tvpplcl = tvp(icb-1) - rgas*tvp(icb-1)*(p(icb-1)-plcl)/(cpn(icb-1)*p(icb-1))
    tvaplcl = tv(icb) + (tvp(icb)-tvp(icb+1))*(plcl-p(icb))/(p(icb)-p(icb+1))
    dtpbl = d_zero
    do i = nk , icb - 1
      dtpbl = dtpbl + (tvp(i)-tv(i))*(ph(i)-ph(i+1))
    end do
    dtpbl = dtpbl/(ph(nk)-ph(icb))
    dtmnx = tvpplcl - tvaplcl + dtmax + dtpbl
    dtma = dtmnx
!
!   adjust cloud base mass flux
!
    cbmfold = cbmf
    delt0 = 300.0D0
    damps = damp*dtcum/delt0
    cbmf = (d_one-damps)*cbmf + 0.1D0*alphae*dtma
    cbmf = dmax1(cbmf,d_zero)
!
!   if cloud base mass flux is zero, skip rest of calculation
!
    if ( dabs(cbmf) < mincbmf .and. dabs(cbmfold) < mincbmf ) return
!
!   calculate rates of mixing, m(i)
!
    m(icb) = d_zero
    do i = icb + 1 , inb
      k = min0(i,inb1)
      dbo = dabs(tv(k)-tvp(k)) + entp*0.02D0*(ph(k)-ph(k+1))
      dbosum = dbosum + dbo
      m(i) = cbmf*dbo
    end do
    do i = icb + 1 , inb
      m(i) = m(i)/dbosum
    end do
!
!   calculate entrained air mass flux (ment), total water mixing ratio (qent),
!   total condensed water (elij), and mixing fraction (sij)         
!
    do i = icb + 1 , inb
      qti = q(nk) - ep(i)*clw(i)
      do j = icb , inb
        bf2 = d_one + lv(j)*lv(j)*qs(j)/(rwat*t(j)*t(j)*cpd)
        anum = h(j) - hp(i) + (cpv-cpd)*t(j)*(qti-q(j))
        denom = h(i) - hp(i) + (cpd-cpv)*(q(i)-qti)*t(j)
        dei = denom
        if ( dabs(dei) < 0.01D0 ) dei = 0.01D0
        sij(i,j) = anum/dei
        sij(i,i) = d_one
        altem = sij(i,j)*q(i) + (d_one-sij(i,j))*qti - qs(j)
        altem = altem/bf2
        cwat = clw(j)*(d_one-ep(j))
        stemp = sij(i,j)
        if ( (stemp < d_zero .or. stemp > d_one .or. &
              altem > cwat) .and. j > i ) then
          anum = anum - lv(j)*(qti-qs(j)-cwat*bf2)
          denom = denom + lv(j)*(q(i)-qti)
          if ( dabs(denom) < 0.01D0 ) denom = 0.01D0
          sij(i,j) = anum/denom
          altem = sij(i,j)*q(i) + (d_one-sij(i,j))*qti - qs(j)
          altem = altem - (bf2-d_one)*cwat
        end if
        if ( sij(i,j) > d_zero .and. sij(i,j) < 0.9D0 ) then
          qent(i,j) = sij(i,j)*q(i) + (d_one-sij(i,j))*qti
          uent(i,j) = sij(i,j)*u(i) + (d_one-sij(i,j))*u(nk)
          vent(i,j) = sij(i,j)*v(i) + (d_one-sij(i,j))*v(nk)
          do k = 1 , ntra
            traent(i,j,k) = sij(i,j)*tra(i,k) + (d_one-sij(i,j))*tra(nk,k)
          end do
          elij(i,j) = altem
          elij(i,j) = dmax1(d_zero,elij(i,j))
          ment(i,j) = m(i)/(d_one-sij(i,j))
          nent(i) = nent(i) + 1
        end if
        sij(i,j) = dmax1(d_zero,sij(i,j))
        sij(i,j) = dmin1(d_one,sij(i,j))
      end do
!
!     if no air can entrain at level i assume that updraft detrains
!     at that level and calculate detrained air flux and properties
!
      if ( nent(i) == 0 ) then
        ment(i,i) = m(i)
        qent(i,i) = q(nk) - ep(i)*clw(i)
        uent(i,i) = u(nk)
        vent(i,i) = v(nk)
        do j = 1 , ntra
          traent(i,i,j) = tra(nk,j)
        end do
        elij(i,i) = clw(i)
        sij(i,i) = d_one
      end if
    end do
    sij(inb,inb) = d_one
!
!   normalize entrained air mass fluxes to represent equal
!   probabilities of mixing
!
    do i = icb + 1 , inb
      if ( nent(i) /= 0 ) then
        qp1 = q(nk) - ep(i)*clw(i)
        anum = h(i) - hp(i) - lv(i)*(qp1-qs(i))
        denom = h(i) - hp(i) + lv(i)*(q(i)-qp1)
        if ( dabs(denom) < 0.01D0 ) denom = 0.01D0
        scrit = anum/denom
        alt = qp1 - qs(i) + scrit*(q(i)-qp1)
        if ( alt < d_zero ) scrit = d_one
        scrit = dmax1(scrit,d_zero)
        asij = d_zero
        smin = d_one
        do j = icb , inb
          if ( sij(i,j) > d_zero .and. sij(i,j) < 0.9D0 ) then
            if ( j > i ) then
              smid = dmin1(sij(i,j),scrit)
              sjmax = smid
              sjmin = smid
              if ( smid < smin .and. sij(i,j+1) < smid ) then
                smin = smid
                sjmax = dmin1(sij(i,j+1),sij(i,j),scrit)
                sjmin = dmax1(sij(i,j-1),sij(i,j))
                sjmin = dmin1(sjmin,scrit)
              end if
            else
              sjmax = dmax1(sij(i,j+1),scrit)
              smid = dmax1(sij(i,j),scrit)
              sjmin = d_zero
              if ( j > 1 ) sjmin = sij(i,j-1)
              sjmin = dmax1(sjmin,scrit)
            end if
            delp = dabs(sjmax-smid)
            delm = dabs(sjmin-smid)
            asij = asij + (delp+delm)*(ph(j)-ph(j+1))
            ment(i,j) = ment(i,j)*(delp+delm)*(ph(j)-ph(j+1))
          end if
        end do
        asij = dmax1(1.0D-21,asij)
        asij = d_one/asij
        do j = icb , inb
          ment(i,j) = ment(i,j)*asij
        end do
        bsum = d_zero
        do j = icb , inb
          bsum = bsum + ment(i,j)
        end do
        if ( bsum < 1.0D-18 ) then
          nent(i) = 0
          ment(i,i) = m(i)
          qent(i,i) = q(nk) - ep(i)*clw(i)
          uent(i,i) = u(nk)
          vent(i,i) = v(nk)
          do j = 1 , ntra
            traent(i,i,j) = tra(nk,j)
          end do
          elij(i,i) = clw(i)
          sij(i,i) = d_one
        end if
      end if
    end do
!
!   check whether ep(inb)=0, if so, skip precipitating
!   downdraft calculation
!
    if ( ep(inb) >= 0.0001D0 ) then
!
!     integrate liquid water equation to find condensed water
!     and condensed water flux
!
      jtt = 2
!
!     begin downdraft loop                   
!
      do i = inb , 1 , -1
!
!       calculate detrained precipitation           
!
        wdtrain = egrav*ep(i)*m(i)*clw(i)
        if ( i > 1 ) then
          do j = 1 , i - 1
            awat = elij(j,i) - (d_one-ep(i))*clw(i)
            awat = dmax1(d_zero,awat)
            wdtrain = wdtrain + egrav*awat*ment(j,i)
          end do
        end if
!
!       find rain water and evaporation using provisional
!       estimates of qp(i)and qp(i-1)
!
!       value of terminal velocity and coefficient of evaporation for snow
!
        coeff = coeffs
        wt(i) = omtsnow
!
!       value of terminal velocity and coefficient of evaporation for rain
!
        if ( t(i) > tzero) then
          coeff = coeffr
          wt(i) = omtrain
        end if
        qsm = (q(i)+qp(i+1))*d_half
        afac = coeff*ph(i)*(qs(i)-qsm)/(1.0D4+2.0D3*ph(i)*qs(i))
        afac = dmax1(afac,d_zero)
        sigt = sigp(i)
        sigt = dmax1(d_zero,sigt)
        sigt = dmin1(d_one,sigt)
        b6 = d_100*(ph(i)-ph(i+1))*sigt*afac/wt(i)
        c6 = (water(i+1)*wt(i+1)+wdtrain/sigd)/wt(i)
        revap = (-b6+dsqrt(b6*b6+d_four*c6))*d_half
        evap(i) = sigt*afac*revap
        water(i) = revap*revap
!
!       calculate precipitating downdraft mass flux under
!       hydrostatic approximation
!
        if ( i /= 1 ) then
          dhdp = (h(i)-h(i-1))/(p(i-1)-p(i))
          dhdp = dmax1(dhdp,d_10)
          mp(i) = d_100*regrav*lv(i)*sigd*evap(i)/dhdp
          mp(i) = dmax1(mp(i),d_zero)
!
!         add small amount of inertia to downdraft             
!
          fac = 20.0D0/(ph(i-1)-ph(i))
          mp(i) = (fac*mp(i+1)+mp(i))/(d_one+fac)
!
!         force mp to decrease linearly to zero             
!         between about 950 mb and the surface          
!
          if ( p(i) > (0.949D0*p(1)) ) then
            jtt = max0(jtt,i)
            mp(i) = mp(jtt)*(p(1)-p(i))/(p(1)-p(jtt))
          end if
        end if
!
!       find mixing ratio of precipitating downdraft
!
        if ( i /= inb ) then
          if ( i == 1 ) then
            qstm = qs(1)
          else
            qstm = qs(i-1)
          end if
          if ( mp(i) > mp(i+1) ) then
            rat = mp(i+1)/mp(i)
            qp(i) = qp(i+1)*rat + q(i)*(d_one-rat) + &
                    d_100*regrav*sigd*(ph(i)-ph(i+1))*(evap(i)/mp(i))
            up(i) = up(i+1)*rat + u(i)*(d_one-rat)
            vp(i) = vp(i+1)*rat + v(i)*(d_one-rat)
            do j = 1 , ntra
              trap(i,j) = trap(i+1,j)*rat + trap(i,j)*(d_one-rat)
            end do
          else if ( mp(i+1) > d_zero ) then
            qp(i) = (gz(i+1)-gz(i)+qp(i+1)*(lv(i+1)+t(i+1)*(cl-cpd)) + &
                     cpd*(t(i+1)-t(i)))/(lv(i)+t(i)*(cl-cpd))
            up(i) = up(i+1)
            vp(i) = vp(i+1)
            do j = 1 , ntra
              trap(i,j) = trap(i+1,j)
            end do
          end if
          qp(i) = dmin1(qp(i),qstm)
          qp(i) = dmax1(qp(i),d_zero)
        end if
      end do
!
!     calculate surface precipitation in mm/s
!
!     precip = precip+wt(1)*sigd*water(1)*3600.*24000./(d_1000*g)   ! mm/d
      precip = precip + wt(1)*sigd*water(1)/egrav ! mm/s
    end if
!
!
!   calculate downdraft velocity scale and surface temperature and
!   water vapor fluctuations         
!
    wd = betae*dabs(mp(icb))*0.01D0*rgas*t(icb)/(sigd*p(icb))
    qprime = (qp(1)-q(1))*d_half
    tprime = wlhv*qprime*rcpd
!
!   calculate tendencies of lowest level potential temperature and mixing ratio
!
    dpinv = 0.01D0/(ph(1)-ph(2))
    am = d_zero
    if ( nk == 1 ) then
      do k = 2 , inb
        am = am + m(k)
      end do
    end if
    if ( (d_two*egrav*dpinv*am) >= delti ) iflag = 4
    ft(1) = ft(1) + egrav*dpinv*am*(t(2)-t(1)+(gz(2)-gz(1))/cpn(1))
    ft(1) = ft(1) - lvcp(1)*sigd*evap(1)
    ft(1) = ft(1) + sigd*wt(2)*(cl-cpd)*water(2)*(t(2)-t(1))*dpinv/cpn(1)
    fq(1) = fq(1) + egrav*mp(2)*(qp(2)-q(1))*dpinv + sigd*evap(1)
    fq(1) = fq(1) + egrav*am*(q(2)-q(1))*dpinv
    fu(1) = fu(1) + egrav*dpinv*(mp(2)*(up(2)-u(1))+am*(u(2)-u(1)))
    fv(1) = fv(1) + egrav*dpinv*(mp(2)*(vp(2)-v(1))+am*(v(2)-v(1)))
    do j = 1 , ntra
      ftra(1,j) = ftra(1,j) + egrav*dpinv * &
               (mp(2)*(trap(2,j)-tra(1,j)) + am*(tra(2,j)-tra(1,j)))
    end do
    do j = 2 , inb
      fq(1) = fq(1) + egrav*dpinv*ment(j,1)*(qent(j,1)-q(1))
      fu(1) = fu(1) + egrav*dpinv*ment(j,1)*(uent(j,1)-u(1))
      fv(1) = fv(1) + egrav*dpinv*ment(j,1)*(vent(j,1)-v(1))
      do k = 1 , ntra
        ftra(1,k) = ftra(1,k) + egrav*dpinv*ment(j,1)*(traent(j,1,k)-tra(1,k))
      end do
    end do
!
!   calculate tendencies of potential temperature and mixing ratio
!   at levels above the lowest level    
!   first find the net saturated updraft and downdraft mass fluxes
!   through each level          
!
    do i = 2 , inb
      dpinv = 0.01D0/(ph(i)-ph(i+1))
      cpinv = d_one/cpn(i)
      amp1 = d_zero
      ad = d_zero
      if ( i >= nk ) then
        do k = i + 1 , inb + 1
          amp1 = amp1 + m(k)
        end do
      end if
      do k = 1 , i
        do j = i + 1 , inb + 1
          amp1 = amp1 + ment(k,j)
        end do
      end do
      if ( (d_two*egrav*dpinv*amp1) >= delti ) iflag = 4
      do k = 1 , i - 1
        do j = i , inb
          ad = ad + ment(j,k)
        end do
      end do
      ft(i) = ft(i) + egrav*dpinv*(amp1*(t(i+1)-t(i)+ &
              (gz(i+1)-gz(i))*cpinv)-ad*(t(i)-t(i-1)+ &
              (gz(i)-gz(i-1))*cpinv)) - sigd*lvcp(i)*evap(i)
      ft(i) = ft(i) + egrav*dpinv*ment(i,i) * &
              (hp(i)-h(i)+t(i)*(cpv-cpd)*(q(i)-qent(i,i)))*cpinv
      ft(i) = ft(i) + sigd*wt(i+1)*(cl-cpd)*water(i+1)*(t(i+1)-t(i))*dpinv*cpinv
      fq(i) = fq(i) + egrav*dpinv*(amp1*(q(i+1)-q(i))-ad*(q(i)-q(i-1)))
      fu(i) = fu(i) + egrav*dpinv*(amp1*(u(i+1)-u(i))-ad*(u(i)-u(i-1)))
      fv(i) = fv(i) + egrav*dpinv*(amp1*(v(i+1)-v(i))-ad*(v(i)-v(i-1)))
      do k = 1 , ntra
        ftra(i,k) = ftra(i,k) + egrav*dpinv * &
             (amp1*(tra(i+1,k)-tra(i,k)) - ad*(tra(i,k)-tra(i-1,k)))
      end do
      do k = 1 , i - 1
        awat = elij(k,i) - (d_one-ep(i))*clw(i)
        awat = dmax1(awat,d_zero)
        fq(i) = fq(i) + egrav*dpinv*ment(k,i)*(qent(k,i)-awat-q(i))
        fu(i) = fu(i) + egrav*dpinv*ment(k,i)*(uent(k,i)-u(i))
        fv(i) = fv(i) + egrav*dpinv*ment(k,i)*(vent(k,i)-v(i))
        do j = 1 , ntra
          ftra(i,j) = ftra(i,j)+egrav*dpinv*ment(k,i)*(traent(k,i,j)-tra(i,j))
        end do
      end do
      do k = i , inb
        fq(i) = fq(i) + egrav*dpinv*ment(k,i)*(qent(k,i)-q(i))
        fu(i) = fu(i) + egrav*dpinv*ment(k,i)*(uent(k,i)-u(i))
        fv(i) = fv(i) + egrav*dpinv*ment(k,i)*(vent(k,i)-v(i))
        do j = 1 , ntra
          ftra(i,j) = ftra(i,j) + egrav*dpinv*ment(k,i)*(traent(k,i,j)-tra(i,j))
        end do
      end do
      fq(i) = fq(i) + sigd*evap(i) + egrav * &
              (mp(i+1)*(qp(i+1)-q(i)) - mp(i)*(qp(i)-q(i-1)))*dpinv
      fu(i) = fu(i) + egrav * &
              (mp(i+1)*(up(i+1)-u(i)) - mp(i)*(up(i)-u(i-1)))*dpinv
      fv(i) = fv(i) + egrav * &
              (mp(i+1)*(vp(i+1)-v(i)) - mp(i)*(vp(i)-v(i-1)))*dpinv
      do j = 1 , ntra
        ftra(i,j) = ftra(i,j) + egrav*dpinv * &
               (mp(i+1)*(trap(i+1,j)-tra(i,j)) - mp(i)*(trap(i,j)-trap(i-1,j)))
      end do
    end do
!
!   adjust tendencies at top of convection layer to reflect
!   actual position of the level zero cape
!
    fqold = fq(inb)
    fq(inb) = fq(inb)*(d_one-frac)
    fq(inb-1) = fq(inb-1) + frac*fqold * &
                ((ph(inb)-ph(inb+1))/(ph(inb-1)-ph(inb)))*lv(inb)/lv(inb-1)
    ftold = ft(inb)
    ft(inb) = ft(inb)*(d_one-frac)
    ft(inb-1) = ft(inb-1) + frac*ftold * &
                ((ph(inb)-ph(inb+1))/(ph(inb-1)-ph(inb)))*cpn(inb)/cpn(inb-1)
    fuold = fu(inb)
    fu(inb) = fu(inb)*(d_one-frac)
    fu(inb-1) = fu(inb-1) + frac*fuold * &
                ((ph(inb)-ph(inb+1))/(ph(inb-1)-ph(inb)))
    fvold = fv(inb)
    fv(inb) = fv(inb)*(d_one-frac)
    fv(inb-1) = fv(inb-1) + frac*fvold * &
                ((ph(inb)-ph(inb+1))/(ph(inb-1)-ph(inb)))
    do k = 1 , ntra
      ftraold = ftra(inb,k)
      ftra(inb,k) = ftra(inb,k)*(d_one-frac)
      ftra(inb-1,k) = ftra(inb-1,k) + frac*ftraold * &
                      (ph(inb)-ph(inb+1))/(ph(inb-1)-ph(inb))
    end do
!
!   very slightly adjust tendencies to force exact
!   enthalpy, momentum and tracer conservation
!
    ents = d_zero
    uav = d_zero
    vav = d_zero
    do i = 1 , inb
      ents = ents + (cpn(i)*ft(i)+lv(i)*fq(i))*(ph(i)-ph(i+1))
      uav = uav + fu(i)*(ph(i)-ph(i+1))
      vav = vav + fv(i)*(ph(i)-ph(i+1))
    end do
    ents = ents/(ph(1)-ph(inb+1))
    uav = uav/(ph(1)-ph(inb+1))
    vav = vav/(ph(1)-ph(inb+1))
    do i = 1 , inb
      ft(i) = ft(i) - ents/cpn(i)
      fu(i) = (d_one-cu)*(fu(i)-uav)
      fv(i) = (d_one-cu)*(fv(i)-vav)
    end do
    do k = 1 , ntra
      traav = d_zero
      do i = 1 , inb
        traav = traav + ftra(i,k)*(ph(i)-ph(i+1))
      end do
      traav = traav/(ph(1)-ph(inb+1))
      do i = 1 , inb
        ftra(i,k) = ftra(i,k) - traav
      end do
    end do
!
  end subroutine cupeman
!
! Calculate lifting level temperature
!
  subroutine tlift(p,t,q,qs,gz,icb,nk,tvp,tpk,clw,nd,nl,kk)

    implicit none
!
    integer :: icb , kk , nd , nk , nl
    real(dp) , dimension(nd) :: clw , gz , p , q , qs , t , tpk , tvp
    intent (in) gz , icb , kk , nd , nk , nl , p , q , qs , t
    intent (out) tvp
    intent (inout) clw , tpk
!
    real(dp) :: ah0 , ahg , alv , cpinv , cpp , denom , es , &
               qg , rg , s , tc , tg
    integer :: i , j , nsb , nst
!
!   calculate certain parcel quantities, including static energy
!
    ah0 = (cpd*(d_one-q(nk))+cl*q(nk))*t(nk) + q(nk) * &
          (wlhv-cpvmcl*(t(nk)-tzero)) + gz(nk)
    cpp = cpd*(d_one-q(nk)) + q(nk)*cpv
    cpinv = d_one/cpp
!
    if ( kk == 1 ) then
!
!     calculate lifted parcel quantities below cloud base
!
      do i = 1 , icb - 1
        clw(i) = d_zero
      end do
      do i = nk , icb - 1
        tpk(i) = t(nk) - (gz(i)-gz(nk))*cpinv
        tvp(i) = tpk(i)*(d_one+q(nk)*rgowi)
      end do
    end if
!
!   find lifted parcel quantities above cloud base
!
    nst = icb
    nsb = icb
    if ( kk == 2 ) then
      nst = nl
      nsb = icb + 1
    end if
    do i = nsb , nst
      tg = t(i)
      qg = qs(i)
      alv = wlhv - cpvmcl*(t(i)-tzero)
      do j = 1 , 2
        s = cpd + alv*alv*qg/(rwat*t(i)*t(i))
        s = d_one/s
        ahg = cpd*tg + (cl-cpd)*q(nk)*t(i) + alv*qg + gz(i)
        tg = tg + s*(ah0-ahg)
        tg = dmax1(tg,35.0D0)
        tc = tg - tzero
        denom = 243.5D0 + tc
        if ( tc >= d_zero ) then
          es = 6.112D0*dexp(17.67D0*tc/denom)
        else
          es = dexp(23.33086D0-6111.72784D0/tg+0.15215D0*dlog(tg))
        end if
        qg = rgow*es/(p(i)-es*(d_one-rgow))
      end do
      alv = wlhv - cpvmcl*(t(i)-tzero)
      tpk(i) = (ah0-(cl-cpd)*q(nk)*t(i)-gz(i)-alv*qg)*rcpd
      clw(i) = q(nk) - qg
      clw(i) = dmax1(d_zero,clw(i))
      rg = qg/(d_one-q(nk))
      tvp(i) = tpk(i)*(d_one+rg*rgowi)
    end do
!
  end subroutine tlift
!
end module mod_cu_em
