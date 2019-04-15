!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or
!    modify
!    it under the terms of the GNU General Public License as
!    published by
!    the Free Software Foundation, either version 3 of the
!    License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty
!    of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public
!    License
!    along with ICTP RegCM.  If not, see
!    <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_che_common_isorropia

  ! ------------
  ! Terms of Use
  ! ------------
  !
  ! These codes not to be included in any commercial package, or used for any
  ! commercial applications (or for profit) without prior authorization from the
  ! code authors (CF, AN, SP, CP).
  ! The code is to be used for educational or non-profit purposes only.
  ! Any other usage must first have authorization from the authors
  ! ISORROPIA/ISORROPIA-II cannot be modified in any way without the author's
  ! consent.
  ! No portion of the ISORROPIA/ISORROPIA-II source code can be used in other
  ! codes without the author's consent.
  ! The codes are provided as-is, and the authors have no liability from its
  ! usage.
  ! Usage of the model, for any purpose (educational, research, or other) must
  ! acknowledge the usage of these codes, i.e.
  ! Links to the ISORROPIA webpage must be provided
  ! (http://nenes.eas.gatech.edu/ISORROPIA) where users can download the latest
  ! version of the code, as well as manuals and other materials.
  ! The main ISORROPIA reference (Nenes et al., Aquatic Geochemistry, 1998) must
  ! be cited in all publications and documentation.
  ! The main ISORROPIA-II reference (Fountoukis and Nenes, Atmos.Chem.Phys.,
  ! 2007) must be cited in all publications and documentation.
  ! If ISORROPIA/ISORROPIA-II is to be included within another model, some
  ! kind of agreement is required (such as an e-mail confirmation) that
  ! subsequent users will abide to the terms as outlined here.
  ! Documentation and publications using ISORROPIA/ISORROPIA-II should cite the
  ! ISORROPIA webpage.
  !

  use mod_realkinds
  use mod_intkinds

  implicit none

  !=======================================================================
  ! *** ISORROPIA CODE II
  ! *** INCLUDE FILE 'ISRPIA.INC'
  ! *** THIS FILE CONTAINS THE DECLARATIONS OF THE GLOBAL CONSTANTS
  !     AND VARIABLES.
  !
  ! *** COPYRIGHT 1996-2008 , UNIVERSITY OF MIAMI , CARNEGIE MELLON UNIVERSITY ,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  ! *** UPDATED BY CHRISTOS FOUNTOUKIS
  !
  !=======================================================================

  ! Integer Parameters

  integer(ik4) , parameter :: ncomp  =   8
  integer(ik4) , parameter :: nions  =  10
  integer(ik4) , parameter :: ngasaq =   3
  integer(ik4) , parameter :: nslds  =  19
  integer(ik4) , parameter :: npair  =  23
  integer(ik4) , parameter :: nzsr   = 100
  integer(ik4) , parameter :: nerrmx =  25

  !=================================================================
  ! Input variables
  !=================================================================

  integer(ik4) :: iprob , metstbl
  integer(ik4) :: nadj
  real(rkx) , dimension(ncomp) :: w , waer
  real(rkx) :: temp , rh
!      common /inpt/ w,waer,temp,rh,iprob,metstbl,nadj

  !=================================================================
  ! Water activities of pure salt solutions
  !=================================================================

  real(rkx) , dimension(nzsr) :: awas , awss , awac , awsc
  real(rkx) , dimension(nzsr) :: awan , awsn , awsb , awab
  real(rkx) , dimension(nzsr) :: awsa , awlc , awcs , awcn
  real(rkx) , dimension(nzsr) :: awcc , awps , awpb , awpn
  real(rkx) , dimension(nzsr) :: awpc , awms , awmn , awmc

!      common /zsr/ awas,awss,awac,awsc,awan,awsn,awsb,awab,
!     &             awsa,awlc,awcs,awcn,awcc,awps,awpb,awpn,
!     &             awpc,awms,awmn,awmc

  !=================================================================
  ! Deliquescence relative humidities
  !=================================================================

  integer(ik4) wftyp
  real(rkx) :: drh2so4 ,  drnh42s4 , drnahso4 , drnacl ,   drnano3
  real(rkx) :: drna2so4 , drnh4hs4 , drlc ,     drnh4no3 , drnh4cl
  real(rkx) :: drcaso4 ,  drcano32 , drcacl2 ,  drk2so4 ,  drkhso4
  real(rkx) :: drkno3 , drkcl , drmgso4 , drmgno32 , drmgcl2

!      common /drh / drh2so4,drnh42s4,drnahso4,drnacl,drnano3,
!     &              drna2so4,drnh4hs4,drlc,drnh4no3,drnh4cl,
!     &              drcaso4,drcano32,drcacl2,drk2so4,drkhso4,
!     &              drkno3,drkcl,drmgso4,drmgno32,drmgcl2

  real(rkx) ::  drmlcab , drmlcas , drmasan , drmg1 , drmg2
  real(rkx) ::  drmg3 , drmh1 , drmh2 , drmi1 , drmi2
  real(rkx) ::  drmi3 , drmq1 , drmr1 , drmr2 , drmr3
  real(rkx) ::  drmr4 , drmr5 , drmr6 , drmr7 , drmr8
  real(rkx) ::  drmr9 , drmr10 , drmr11 , drmr12 , drmr13

!      common /mdrh/ drmlcab,drmlcas,drmasan,drmg1,drmg2,
!     &              drmg3,drmh1,drmh2,drmi1,drmi2,
!     &              drmi3,drmq1,drmr1,drmr2,drmr3,
!     &              drmr4,drmr5,drmr6,drmr7,drmr8,
!     &              drmr9,drmr10,drmr11,drmr12,drmr13,
!     &              wftyp

  real(rkx) :: drmo1 , drmo2 , drmo3 , drml1 , drml2
  real(rkx) :: drml3 , drmm1 , drmm2 , drmp1 , drmp2
  real(rkx) :: drmp3 , drmp4 , drmp5 , drmv1

!      common /mdrh2/ drmo1,drmo2,drmo3,drml1,drml2,
!     &               drml3,drmm1,drmm2,drmp1,drmp2,
!     &               drmp3,drmp4,drmp5,drmv1

  !=================================================================
  ! Variables for liquid aerosol phase
  !=================================================================

  real(rkx) , dimension(nions) :: molal
  real(rkx) , dimension(npair) :: molalr , m0
  real(rkx) , dimension(nions) :: z
  real(rkx) , dimension(npair) :: zz
  real(rkx) :: epsact
  real(rkx) , dimension(npair) :: gamou , gamin
  real(rkx) , dimension(ngasaq) :: gasaq
  real(rkx) , dimension(npair) :: gama
  real(rkx) :: coh , chno3 , chcl , water
  integer(ik4) :: iacalc
  real(rkx) :: ionic
  logical :: calaou , calain , frst , dryf

!      common /ions/ molal,molalr,gama,zz,
!     &              z,gamou,gamin,m0,
!     &              gasaq,epsact,coh,chno3,
!     &              chcl,water,ionic,iacalc,
!     &              frst,calain,calaou,dryf

  !=================================================================
  ! Variables for solid aerosol phase
  !=================================================================

  real(rkx) :: ch2so4 ,  cnh42s4 , cnh4hs4 , cnacl ,   cna2so4
  real(rkx) :: cnano3 ,  cnh4no3 , cnh4cl ,  cnahso4 , clc , ccaso4
  real(rkx) :: ccano32 , ccacl2 ,  ck2so4 ,  ckhso4 ,  ckno3 , ckcl
  real(rkx) :: cmgso4 ,  cmgno32 , cmgcl2

!      common /salt/ ch2so4,cnh42s4,cnh4hs4,cnacl,cna2so4,
!     &              cnano3,cnh4no3,cnh4cl,cnahso4,clc,ccaso4,
!     &              ccano32,ccacl2,ck2so4,ckhso4,ckno3,ckcl,
!     &              cmgso4,cmgno32,cmgcl2

  !=================================================================
  ! Variables for gas phase
  !=================================================================

  real(rkx) :: gnh3 , ghno3 , ghcl

!      common /gas / gnh3,ghno3,ghcl

  !=================================================================
  ! Equilibrium constants
  !=================================================================

  real(rkx) :: xk1 , xk2 , xk3 , xk4 , xk5 , xk6 , xk7 , xk8 , xk9 , xk10
  real(rkx) :: xk11 , xk12 , xk13 , xk14 , xkw , xk21 , xk22
  real(rkx) :: xk31 , xk32 , xk41 , xk42
  real(rkx) :: xk15 , xk16 , xk17 , xk18 , xk19 , xk20 , xk23
  real(rkx) :: xk24 , xk25

!      common /equk/ xk1,xk2,xk3,xk4,xk5,xk6,xk7,xk8,xk9,xk10,
!     &              xk11,xk12,xk13,xk14,xkw,xk21,xk22,xk31,xk32,xk41,
!     &              xk42,xk15,xk16,xk17,xk18,xk19,xk20,xk23,
!     &              xk24,xk25,xk26,xk27

  !=================================================================
  ! Molecular Weights
  !=================================================================

  real(rkx) :: r
  real(rkx) , dimension(nions) :: imw
  real(rkx) , dimension(ncomp) :: wmw
  real(rkx) , dimension(npair) :: smw

!      common /othr/ r,imw,wmw,smw

  !=================================================================
  ! Solution/info variables
  !=================================================================

  character(len=15) :: scase
  real(rkx) :: sulratw , sulrat ,  sodrat , so4rat ,  crnarat , crrat

!      common /case/ sulratw,sulrat,sodrat,so4rat,crnarat,crrat,scase

  real(rkx) :: eps
  integer(ik4) :: maxit , nsweep , ndiv , iclact

!      common /soln/ eps,maxit,nsweep,ndiv,iclact

  !=================================================================
  ! Error system
  !=================================================================

  character(len=40) , dimension(nerrmx) :: errmsg
  integer(ik4) , dimension(nerrmx) :: errstk
  integer(ik4) :: nofer
  logical :: stkofl

!      common /eror/ stkofl,nofer,errstk,errmsg

  !=================================================================
  ! Generic Variables
  !=================================================================

  character(len=15) :: version
  real(rkx) :: great , tiny1 , tiny2 , zero , one

!      common /cgen/ great,tiny,tiny2,zero,one,version

end module mod_che_common_isorropia

module mod_che_common_solut
  use mod_realkinds
  use mod_intkinds
  real(rkx) :: chi1 , chi2 , chi3 , chi4 , chi5 , chi6 , chi7 , chi8
  real(rkx) :: chi9 , chi10 , chi11 , chi12 , chi13 , chi14 , chi15
  real(rkx) :: chi16 , chi17 , psi1 , psi2 , psi3 , psi4 , psi5 , psi6
  real(rkx) :: psi7 , psi8 , psi9 , psi10 , psi11 , psi12 , psi13
  real(rkx) :: psi14 , psi15 , psi16 , psi17 , a1 , a2 , a3 , a4 , a5 , a6
  real(rkx) :: a7 , a8 , a9 , a10 , a11 , a12 , a13 , a14 , a15 , a16 , a17
end module mod_che_common_solut

module mod_che_common_caseg
  use mod_realkinds
  use mod_intkinds
  real(rkx) :: chi1 , chi2 , chi3 , chi4 , chi5 , chi6 , lamda
  real(rkx) :: psi1 , psi2 , psi3 , psi4 , psi5 , psi6 , psi7
  real(rkx) :: a1 , a2 , a3 , a4 , a5 , a6 , a7
end module mod_che_common_caseg

module mod_che_common_casej
  use mod_realkinds
  use mod_intkinds
  real(rkx) :: chi1 , chi2 , chi3 , lamda , kapa , psi1 , psi2 , psi3
  real(rkx) :: a1 , a2 , a3
end module mod_che_common_casej

module mod_che_common_casek
  use mod_realkinds
  use mod_intkinds
  real(rkx) :: chi1 , chi2 , chi3 , chi4 , lamda , kapa , psi1 , psi2 , psi3
  real(rkx) :: a1 , a2 , a3 , a4
end module mod_che_common_casek

module mod_che_common_caseo
  use mod_realkinds
  use mod_intkinds
  real(rkx) :: chi1 , chi2 , chi3 , chi4 , chi5 , chi6 , chi7 , chi8
  real(rkx) :: chi9 , lamda , psi1 , psi2 , psi3 , psi4 , psi5
  real(rkx) :: psi6 , psi7 , psi8 , psi9 , a1 , a2 , a3 , a4
  real(rkx) :: a5 , a6 , a7 , a8 , a9
end module mod_che_common_caseo
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
