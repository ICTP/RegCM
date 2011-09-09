       subroutine chemistry(jj,chemin,chemox,taa,psaa,zena, &
                            ktau,idatein,TOD)

       use mod_chem_interface, only : chemall,nvar, taucld
       use mod_indices
       use mod_dynparam
       use mod_constants
       use mod_runparams, only : dtchem
       use mpi
       implicit none

       include 'boxvars.EXT'        ! Box model variables
       include 'chemmech.EXT'       ! Mechanism variables for chemistry
       include 'chemvars.EXT'       ! Chem variables and indices -
      INTEGER     ,INTENT(IN) :: jj
      INTEGER     ,INTENT(IN) :: ktau
      INTEGER     ,INTENT(IN) :: idatein
      REAL(KIND=8),INTENT(IN) :: TOD    !abt added for time of day

      REAL(KIND=8),DIMENSION(2:iym2,1:kz,nvar),INTENT(IN )  :: chemin
      REAL(KIND=8),DIMENSION(2:iym2,1:kz,nvar),INTENT(OUT ) :: chemox
      REAL(KIND=8),DIMENSION(2:iym2,1:kz),INTENT(IN)      :: taa,psaa
      REAL(kind=8),DIMENSION(2:iym2),INTENT(IN)             :: zena
! LOCAL VARIABLES
      INTEGER       :: ixi,klev,ierr
      REAL(KIND=8),DIMENSION(2:iym2,1:kz)      :: cfactor
      REAL(kind=8), PARAMETER :: kb=1.380658E-19
      REAL(KIND=8),DIMENSION(2:iym2,1:kz,1:56)  :: jphoto
!ah  variable for photolysis
      REAL(kind=8),DIMENSION(1:kz)            :: taucab,taucbl
      REAL(kind=8),DIMENSION(1:iy,1:kz,1:jxp) :: taucld2
      REAL(kind=8),DIMENSION(1:iy,0:kz+1)     :: psaa2,taa2
      REAL(kind=8),DIMENSION(1:kz)            :: hcbl,hcab
      REAL(kind=8)                            :: levav
      INTEGER                                 :: kbl,kab,ll
!ah_


      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      time      = 900.
      IDATE     = idatein
      ktaubx    = ktau
      xhour  = TOD        !abt added for time of day
      c_numitr  = 20

      ! initialize jphoto to zero
       jphoto(:,:,:) = 0.0

      ! Reorder from top-down to bottom-up
       DO klev= 1,kz
       DO ixi = 2,iym2
           ll = (kz+1)-klev
           taucld2(ixi,ll,jj) = taucld(ixi,klev,jj)
           psaa2(ixi,ll)      = psaa(ixi,klev)
           taa2(ixi,ll)       = taa(ixi,klev)
       END DO
       END DO


      DO klev=kz,1,-1
      DO ixi=2,iym2

          ll  = (kz+1)-klev
      altmid(1) = psaa2(ixi,klev)
      temp(1)   = taa2(ixi,klev)
      zenith    = zena(ixi) 
      cfactor(ixi,klev)   = psaa2(ixi,klev)*10./(kb*taa2(ixi,klev))
      dens(1)   = cfactor(ixi,klev)
         taucab(ll)  = 0.0
         hcab(ll)    = 0.0
         taucbl(ll)  = 0.0
         hcbl(ll)    = 0.0
         IF(ll .lt. 18)THEN
           DO kab=ll+1,18
             taucab(ll)   = taucab(ll)+taucld2(ixi,kab,jj)
             levav        = 0.5*(psaa2(ixi,kab)+psaa2(ixi,kab+1))
             hcab(ll)     = hcab(ll)  + taucld2(ixi,kab,jj)*levav
           END DO
         END IF
         taucab(ll) = taucab(ll) +taucld2(ixi,ll,jj)*0.5
         levav        = 0.5*(psaa2(ixi,ll)+psaa2(ixi,ll+1))
         hcab(ll)   = hcab(ll)   + taucld2(ixi,ll,jj)*0.5*levav
         IF(taucab(ll) .eq. 0.0 )THEN
           hcab(ll) = 0.0
         ELSE
           hcab(ll)   = hcab(ll)/taucab(ll)
         END IF
         DEPTHA   = taucab(ll)
         ALTABOVE = hcab(ll)
         IF(ll .gt. 1)THEN
           DO kbl=1,ll-1
             taucbl(ll)   = taucbl(ll)+taucld2(ixi,kbl,jj)
             levav        = 0.5*(psaa2(ixi,kbl)+psaa2(ixi,kbl-1))
             hcbl(ll)     = hcbl(ll)  +taucld2(ixi,kbl,jj)*levav
           END DO
         END IF
         taucbl(ll)   = taucbl(ll) + taucld2(ixi,ll,jj)*0.5
         levav        = 0.5*(psaa2(ixi,ll)+psaa2(ixi,ll-1))
         hcbl(ll)   = hcbl(ll)   + taucld2(ixi,ll,jj)*0.5*levav
         IF(taucbl(ll) .eq. 0.0 )THEN
           hcbl(ll) = 0.0
         ELSE
           hcbl(ll)   = hcbl(ll)/taucbl(ll)
         END IF
         DEPTHB     = taucbl(ll)
         ALTBELOW   = hcbl(ll)
         kmax      = 1





      do ic=1,nvar
      xr(1,ic) = 0.0
      end do

      do ic=1,nvar
      xr(1,ic) = chemall(ixi,klev,jj,ic) 
      end do

      xh2o           = chemin(ixi,klev,ind_H2O)
      xr(1,ind_H2O)  = xh2o
      xr(1,ind_H2)   = 0.1E13
      xr(1,ind_O3)   = chemin(ixi,klev,ind_O3) !0.094E+13
      xr(1,ind_NO2)  = chemin(ixi,klev,ind_NO2) !0.300E+09
      xr(1,ind_NO)   = chemin(ixi,klev,ind_NO) !0.300E+09
      xr(1,ind_CO)   = chemin(ixi,klev,ind_CO) !0.100E+13
      xr(1,ind_H2O2) = chemin(ixi,klev,ind_H2O2) !0.200E+11
      xr(1,ind_HNO3) = chemin(ixi,klev,ind_HNO3) !0.200E+10
      xr(1,ind_N2O5) = chemin(ixi,klev,ind_N2O5) !0.100E+08
      xr(1,ind_SO2)  = chemin(ixi,klev,ind_SO2)  !0.200E+11
      xr(1,ind_SULF) = chemin(ixi,klev,ind_SULF) !0.200E+11
      xr(1,ind_DMS)  = chemin(ixi,klev,ind_DMS) !0.200E+10
      xr(1,ind_HCHO) = chemin(ixi,klev,ind_HCHO) !0.200E+10
      xr(1,ind_ALD2) = chemin(ixi,klev,ind_ALD2) !0.200E+10
      xr(1,ind_ISOP) = chemin(ixi,klev,ind_ISOP) !0.500E+10
      xr(1,ind_C2H6) = chemin(ixi,klev,ind_C2H6) !0.500E+10
      xr(1,ind_PAR)  = chemin(ixi,klev,ind_PAR) !0.200E+08
      xr(1,ind_ACET) = chemin(ixi,klev,ind_ACET) !0.200E+10
      xr(1,ind_MOH)  = chemin(ixi,klev,ind_MOH) !0.200E+10
      xr(1,ind_PRPE) = chemin(ixi,klev,ind_PRPE) !0.200E+09
      xr(1,ind_BUTE) = chemin(ixi,klev,ind_BUTE) !0.200E+07
      xr(1,ind_TOLU) = chemin(ixi,klev,ind_TOLU) !0.200E+07
      xr(1,ind_XYLE) = chemin(ixi,klev,ind_XYLE) !0.000E+10
      xr(1,ind_CH4)  = chemin(ixi,klev,ind_CH4)
      xr(1,ind_PAN)  = chemin(ixi,klev,ind_PAN) !0.750E+10
      xr(1,ind_ETHE) = chemin(ixi,klev,ind_ETHE) !10.200E+09

      call chemmain

      do ic=1,nvar
      chemall(ixi,klev,jj,ic) =xr(1,ic)
      end do
       ! Store photolysis rates
          DO ic=1,56
          jphoto(ixi,klev,ic) = c_jval(1,ic)
          END DO

       chemox(ixi,klev,ind_O3)    =   xr(1,ind_O3)
       chemox(ixi,klev,ind_NO2)   =   xr(1,ind_NO2)
       chemox(ixi,klev,ind_NO)    =   xr(1,ind_NO)
       chemox(ixi,klev,ind_CO)    =   xr(1,ind_CO)
       chemox(ixi,klev,ind_H2O2)  =   xr(1,ind_H2O2)
       chemox(ixi,klev,ind_HNO3)  =   xr(1,ind_HNO3)
       chemox(ixi,klev,ind_N2O5)  =   xr(1,ind_N2O5)
       chemox(ixi,klev,ind_SO2)   =   xr(1,ind_SO2)
       chemox(ixi,klev,ind_SULF)  =   xr(1,ind_SULF)
       chemox(ixi,klev,ind_DMS)   =   xr(1,ind_DMS)
       chemox(ixi,klev,ind_HCHO)  =   xr(1,ind_HCHO)
       chemox(ixi,klev,ind_ALD2)  =   xr(1,ind_ALD2)
       chemox(ixi,klev,ind_ISOP)  =   xr(1,ind_ISOP)
       chemox(ixi,klev,ind_C2H6)  =   xr(1,ind_C2H6)
       chemox(ixi,klev,ind_PAR)   =   xr(1,ind_PAR)
       chemox(ixi,klev,ind_TOLU)  =   xr(1,ind_TOLU)
       chemox(ixi,klev,ind_XYLE)  =   xr(1,ind_XYLE)
       chemox(ixi,klev,ind_ETHE)  =   xr(1,ind_ETHE)
       chemox(ixi,klev,ind_PAN)   =   xr(1,ind_PAN)
       chemox(ixi,klev,ind_CH4)   =   xr(1,ind_CH4)
       chemox(ixi,klev,ind_PRPE)  =   xr(1,ind_PRPE)
       chemox(ixi,klev,ind_BUTE)  =   xr(1,ind_BUTE)
       chemox(ixi,klev,ind_MOH)   =   xr(1,ind_MOH)
       chemox(ixi,klev,ind_ACET)  =   xr(1,ind_ACET)


      END DO  !ixi loop
      END DO  !klev loop



!        IF(jj .eq. 1 ) THEN
!        write(*,555)ktau*100,minval(chemall(2:iym2,kz,jj,ind_O3))&
!                         ,maxval(chemall(2:iym2,kz,jj,ind_O3))   &
!                         ,maxval(chemall(2:iym2,kz,jj,ind_OH)) &
!!                        ,maxval(chemall(2:iym2,kz,jj,ind_HO2))  &
!                       ,maxval(chemall(2:iym2,kz,jj,ind_ro2))  &
!                         ,minval(jphoto(2:iym2,kz,4))                    &
!                          ,myid
!         END IF

555     format(I5,5x,5(1pe10.3,2x),I4,1x,I2)

      end subroutine chemistry

