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
  real(rk8) , dimension(ncomp) :: w , waer
  real(rk8) :: temp , rh
!      common /inpt/ w,waer,temp,rh,iprob,metstbl,nadj

  !=================================================================
  ! Water activities of pure salt solutions
  !=================================================================

  real(rk8) , dimension(nzsr) :: awas , awss , awac , awsc
  real(rk8) , dimension(nzsr) :: awan , awsn , awsb , awab
  real(rk8) , dimension(nzsr) :: awsa , awlc , awcs , awcn
  real(rk8) , dimension(nzsr) :: awcc , awps , awpb , awpn
  real(rk8) , dimension(nzsr) :: awpc , awms , awmn , awmc

!      common /zsr/ awas,awss,awac,awsc,awan,awsn,awsb,awab,
!     &             awsa,awlc,awcs,awcn,awcc,awps,awpb,awpn,
!     &             awpc,awms,awmn,awmc

  !=================================================================
  ! Deliquescence relative humidities
  !=================================================================

  integer(ik4) wftyp
  real(rk8) :: drh2so4 ,  drnh42s4 , drnahso4 , drnacl ,   drnano3 
  real(rk8) :: drna2so4 , drnh4hs4 , drlc ,     drnh4no3 , drnh4cl
  real(rk8) :: drcaso4 ,  drcano32 , drcacl2 ,  drk2so4 ,  drkhso4
  real(rk8) :: drkno3 , drkcl , drmgso4 , drmgno32 , drmgcl2

!      common /drh / drh2so4,drnh42s4,drnahso4,drnacl,drnano3, 
!     &              drna2so4,drnh4hs4,drlc,drnh4no3,drnh4cl,
!     &              drcaso4,drcano32,drcacl2,drk2so4,drkhso4,
!     &              drkno3,drkcl,drmgso4,drmgno32,drmgcl2

  real(rk8) ::  drmlcab , drmlcas , drmasan , drmg1 , drmg2
  real(rk8) ::  drmg3 , drmh1 , drmh2 , drmi1 , drmi2
  real(rk8) ::  drmi3 , drmq1 , drmr1 , drmr2 , drmr3
  real(rk8) ::  drmr4 , drmr5 , drmr6 , drmr7 , drmr8
  real(rk8) ::  drmr9 , drmr10 , drmr11 , drmr12 , drmr13

!      common /mdrh/ drmlcab,drmlcas,drmasan,drmg1,drmg2,
!     &              drmg3,drmh1,drmh2,drmi1,drmi2,
!     &              drmi3,drmq1,drmr1,drmr2,drmr3,
!     &              drmr4,drmr5,drmr6,drmr7,drmr8,
!     &              drmr9,drmr10,drmr11,drmr12,drmr13,
!     &              wftyp

  real(rk8) :: drmo1 , drmo2 , drmo3 , drml1 , drml2
  real(rk8) :: drml3 , drmm1 , drmm2 , drmp1 , drmp2
  real(rk8) :: drmp3 , drmp4 , drmp5 , drmv1

!      common /mdrh2/ drmo1,drmo2,drmo3,drml1,drml2,
!     &               drml3,drmm1,drmm2,drmp1,drmp2,
!     &               drmp3,drmp4,drmp5,drmv1

  !=================================================================
  ! Variables for liquid aerosol phase
  !=================================================================

  real(rk8) , dimension(nions) :: molal
  real(rk8) , dimension(npair) :: molalr , m0
  real(rk8) , dimension(nions) :: z
  real(rk8) , dimension(npair) :: zz
  real(rk8) :: epsact
  real(rk8) , dimension(npair) :: gamou , gamin
  real(rk8) , dimension(ngasaq) :: gasaq
  real(rk8) , dimension(npair) :: gama
  real(rk8) :: coh , chno3 , chcl , water
  integer(ik4) :: iacalc
  real(rk8) :: ionic
  logical :: calaou , calain , frst , dryf

!      common /ions/ molal,molalr,gama,zz,
!     &              z,gamou,gamin,m0, 
!     &              gasaq,epsact,coh,chno3,
!     &              chcl,water,ionic,iacalc,
!     &              frst,calain,calaou,dryf
      
  !=================================================================
  ! Variables for solid aerosol phase
  !=================================================================

  real(rk8) :: ch2so4 ,  cnh42s4 , cnh4hs4 , cnacl ,   cna2so4 
  real(rk8) :: cnano3 ,  cnh4no3 , cnh4cl ,  cnahso4 , clc , ccaso4
  real(rk8) :: ccano32 , ccacl2 ,  ck2so4 ,  ckhso4 ,  ckno3 , ckcl
  real(rk8) :: cmgso4 ,  cmgno32 , cmgcl2

!      common /salt/ ch2so4,cnh42s4,cnh4hs4,cnacl,cna2so4, 
!     &              cnano3,cnh4no3,cnh4cl,cnahso4,clc,ccaso4,
!     &              ccano32,ccacl2,ck2so4,ckhso4,ckno3,ckcl,
!     &              cmgso4,cmgno32,cmgcl2

  !=================================================================
  ! Variables for gas phase
  !=================================================================

  real(rk8) :: gnh3 , ghno3 , ghcl 

!      common /gas / gnh3,ghno3,ghcl 

  !=================================================================
  ! Equilibrium constants
  !=================================================================

  real(rk8) :: xk1 , xk2 , xk3 , xk4 , xk5 , xk6 , xk7 , xk8 , xk9 , xk10
  real(rk8) :: xk11 , xk12 , xk13 , xk14 , xkw , xk21 , xk22
  real(rk8) :: xk31 , xk32 , xk41 , xk42
  real(rk8) :: xk15 , xk16 , xk17 , xk18 , xk19 , xk20 , xk23
  real(rk8) :: xk24 , xk25

!      common /equk/ xk1,xk2,xk3,xk4,xk5,xk6,xk7,xk8,xk9,xk10,
!     &              xk11,xk12,xk13,xk14,xkw,xk21,xk22,xk31,xk32,xk41,
!     &              xk42,xk15,xk16,xk17,xk18,xk19,xk20,xk23,
!     &              xk24,xk25,xk26,xk27

  !=================================================================
  ! Molecular Weights
  !=================================================================

  real(rk8) :: r
  real(rk8) , dimension(nions) :: imw
  real(rk8) , dimension(ncomp) :: wmw
  real(rk8) , dimension(npair) :: smw

!      common /othr/ r,imw,wmw,smw

  !=================================================================
  ! Solution/info variables
  !=================================================================

  character(len=15) :: scase
  real(rk8) :: sulratw , sulrat ,  sodrat , so4rat ,  crnarat , crrat

!      common /case/ sulratw,sulrat,sodrat,so4rat,crnarat,crrat,scase

  real(rk8) :: eps
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
  real(rk8) :: great , tiny1 , tiny2 , zero , one

!      common /cgen/ great,tiny,tiny2,zero,one,version

end module mod_che_common_isorropia

module mod_che_common_solut
  use mod_realkinds
  use mod_intkinds
  real(rk8) :: chi1 , chi2 , chi3 , chi4 , chi5 , chi6 , chi7 , chi8
  real(rk8) :: chi9 , chi10 , chi11 , chi12 , chi13 , chi14 , chi15
  real(rk8) :: chi16 , chi17 , psi1 , psi2 , psi3 , psi4 , psi5 , psi6
  real(rk8) :: psi7 , psi8 , psi9 , psi10 , psi11 , psi12 , psi13
  real(rk8) :: psi14 , psi15 , psi16 , psi17 , a1 , a2 , a3 , a4 , a5 , a6
  real(rk8) :: a7 , a8 , a9 , a10 , a11 , a12 , a13 , a14 , a15 , a16 , a17
end module mod_che_common_solut

module mod_che_common_caseg
  use mod_realkinds
  use mod_intkinds
  real(rk8) :: chi1 , chi2 , chi3 , chi4 , chi5 , chi6 , lamda
  real(rk8) :: psi1 , psi2 , psi3 , psi4 , psi5 , psi6 , psi7
  real(rk8) :: a1 , a2 , a3 , a4 , a5 , a6 , a7
end module mod_che_common_caseg

module mod_che_common_casej
  use mod_realkinds
  use mod_intkinds
  real(rk8) :: chi1 , chi2 , chi3 , lamda , kapa , psi1 , psi2 , psi3
  real(rk8) :: a1 , a2 , a3
end module mod_che_common_casej

module mod_che_common_casek
  use mod_realkinds
  use mod_intkinds
  real(rk8) :: chi1 , chi2 , chi3 , chi4 , lamda , kapa , psi1 , psi2 , psi3
  real(rk8) :: a1 , a2 , a3 , a4
end module mod_che_common_casek

module mod_che_common_caseo
  use mod_realkinds
  use mod_intkinds
  real(rk8) :: chi1 , chi2 , chi3 , chi4 , chi5 , chi6 , chi7 , chi8
  real(rk8) :: chi9 , lamda , psi1 , psi2 , psi3 , psi4 , psi5
  real(rk8) :: psi6 , psi7 , psi8 , psi9 , a1 , a2 , a3 , a4
  real(rk8) :: a5 , a6 , a7 , a8 , a9
end module mod_che_common_caseo
