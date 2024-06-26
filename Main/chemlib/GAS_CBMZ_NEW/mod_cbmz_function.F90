! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! The ODE Function of Chemical Model File
!
! Generated by KPP-2.2.3 symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
!
! File                 : mod_cbmz_Function.f90
! Time                 : Mon Nov 25 13:41:15 2013
! Working directory    : /scratch/ashalaby/kpp-2.2.3/compare/CBMZ
! Equation file        : mod_cbmz.kpp
! Output root filename : mod_cbmz
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

module mod_cbmz_function

  use mod_cbmz_parameters
  implicit none

  ! a - rate for each equation
  real(kind=dp) :: a(nreact)

  contains

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !
  ! fun - time derivatives of variables - agregate form
  !   arguments :
  !      v         - concentrations of variable species (local)
  !      f         - concentrations of fixed species (local)
  !      rct       - rate constants (local)
  !      vdot      - time derivative of variable species concentrations
  !
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine fun ( v, f, rct, vdot )

    ! v - concentrations of variable species (local)
    real(kind=dp) :: v(nvar)
    ! f - concentrations of fixed species (local)
    real(kind=dp) :: f(nfix)
    ! rct - rate constants (local)
    real(kind=dp) :: rct(nreact)
    ! vdot - time derivative of variable species concentrations
    real(kind=dp) :: vdot(nvar)

    ! computation of equation rates
    a(1) = rct(1)*v(56)
    a(2) = rct(2)*v(57)
    a(3) = rct(3)*v(21)
    a(4) = rct(4)*v(29)
    a(5) = rct(5)*v(17)
    a(6) = rct(6)*v(53)
    a(7) = rct(7)*v(53)
    a(8) = rct(8)*v(10)
    a(9) = rct(9)*v(9)*f(1)
    a(10) = rct(10)*v(9)*f(2)
    a(11) = 2.2e-10_dp*v(9)*v(35)
    a(12) = rct(12)*v(27)*f(1)
    a(13) = rct(13)*v(27)*v(53)
    a(14) = rct(14)*v(27)*v(56)
    a(15) = rct(15)*v(27)*v(56)
    a(16) = rct(16)*v(27)*v(55)
    a(17) = rct(17)*v(53)*v(55)
    a(18) = rct(18)*v(53)*v(56)
    a(19) = rct(19)*v(53)*v(54)
    a(20) = rct(20)*v(52)*v(53)
    a(21) = rct(21)*v(18)*v(54)
    a(22) = rct(22)*v(54)*v(55)
    a(23) = rct(23)*v(54)*v(56)
    a(24) = 2.2e-11_dp*v(54)*v(57)
    a(25) = rct(25)*v(21)*v(54)
    a(26) = rct(26)*v(29)*v(54)
    a(27) = rct(27)*v(17)*v(54)
    a(28) = rct(28)*v(52)*v(54)
    a(29) = rct(29)*v(10)*v(54)
    a(30) = rct(30)*v(52)*v(52)
    a(31) = rct(31)*v(35)*v(52)*v(52)
    a(32) = rct(32)*v(52)*v(55)
    a(33) = rct(33)*v(52)*v(56)
    a(34) = 5e-16_dp*v(52)*v(56)
    a(35) = rct(35)*v(17)
    a(36) = rct(36)*v(55)*v(57)
    a(37) = rct(37)*v(56)*v(57)
    a(38) = rct(38)*v(56)*v(57)
    a(39) = rct(39)*v(57)*v(57)
    a(40) = 3.5e-12_dp*v(52)*v(57)
    a(41) = 2e-21_dp*v(24)*v(35)
    a(42) = rct(42)*v(24)
    a(43) = rct(43)*v(55)*v(55)*f(1)
    a(44) = rct(44)*v(28)*v(54)
    a(45) = rct(45)*v(11)*v(54)
    a(46) = rct(46)*v(13)*v(54)
    a(47) = rct(47)*v(14)*v(54)
    a(48) = 8.1e-13_dp*v(30)*v(54)
    a(49) = rct(49)*v(20)*v(54)
    a(50) = rct(50)*v(39)
    a(51) = rct(51)*v(39)
    a(52) = 1e-11_dp*v(39)*v(54)
    a(53) = rct(53)*v(39)*v(57)
    a(54) = rct(54)*v(45)
    a(55) = rct(55)*v(45)*v(54)
    a(56) = rct(56)*v(45)*v(57)
    a(57) = rct(57)*v(36)
    a(58) = rct(58)*v(36)*v(54)
    a(59) = rct(59)*v(42)
    a(60) = 1.7e-11_dp*v(42)*v(54)
    a(61) = rct(61)*v(42)*v(57)
    a(62) = rct(62)*v(25)*v(53)
    a(63) = rct(63)*v(25)*v(54)
    a(64) = rct(64)*v(40)*v(53)
    a(65) = rct(65)*v(37)*v(53)
    a(66) = rct(66)*v(40)*v(54)
    a(67) = rct(67)*v(37)*v(54)
    a(68) = rct(68)*v(40)*v(57)
    a(69) = 2.5e-12_dp*v(37)*v(57)
    a(70) = rct(70)*v(8)*v(54)
    a(71) = rct(71)*v(12)*v(54)
    a(72) = 8.1e-12_dp*v(19)*v(55)
    a(73) = 4.1e-11_dp*v(26)*v(54)
    a(74) = 2.2e-11_dp*v(26)*v(57)
    a(75) = 1.4e-11_dp*v(15)*v(56)
    a(76) = 3e-11_dp*v(31)*v(54)
    a(77) = rct(77)*v(31)
    a(78) = rct(78)*v(31)*v(53)
    a(79) = rct(79)*v(38)*v(54)
    a(80) = rct(80)*v(38)*v(53)
    a(81) = rct(81)*v(38)*v(57)
    a(82) = 3.3e-11_dp*v(47)*v(54)
    a(83) = 7e-18_dp*v(47)*v(53)
    a(84) = rct(84)*v(47)
    a(85) = 1e-15_dp*v(47)*v(57)
    a(86) = rct(86)*v(22)
    a(87) = rct(87)*v(23)
    a(88) = rct(88)*v(49)
    a(89) = rct(89)*v(22)*v(54)
    a(90) = rct(90)*v(23)*v(54)
    a(91) = rct(91)*v(49)*v(54)
    a(92) = rct(92)*v(51)*v(54)
    a(93) = rct(93)*v(51)
    a(94) = rct(94)*v(56)*v(58)
    a(95) = rct(95)*v(7)
    a(96) = rct(96)*v(46)*v(55)
    a(97) = rct(97)*v(43)*v(55)
    a(98) = 4e-12_dp*v(50)*v(55)
    a(99) = rct(99)*v(55)*v(58)
    a(100) = 4e-12_dp*v(48)*v(55)
    a(101) = 4e-12_dp*v(44)*v(55)
    a(102) = 4e-12_dp*v(33)*v(55)
    a(103) = 4e-12_dp*v(32)*v(55)
    a(104) = 4e-12_dp*v(34)*v(55)
    a(105) = 4e-12_dp*v(41)*v(55)
    a(106) = 1.1e-12_dp*v(46)*v(57)
    a(107) = 2.5e-12_dp*v(43)*v(57)
    a(108) = 2.5e-12_dp*v(50)*v(57)
    a(109) = 4e-12_dp*v(57)*v(58)
    a(110) = 1.2e-12_dp*v(48)*v(57)
    a(111) = 4e-12_dp*v(44)*v(57)
    a(112) = 2.5e-12_dp*v(41)*v(57)
    a(113) = rct(113)*v(46)*v(52)
    a(114) = rct(114)*v(43)*v(52)
    a(115) = rct(115)*v(50)*v(52)
    a(116) = rct(116)*v(52)*v(58)
    a(117) = rct(117)*v(48)*v(52)
    a(118) = rct(118)*v(44)*v(52)
    a(119) = rct(119)*v(33)*v(52)
    a(120) = rct(120)*v(32)*v(52)
    a(121) = rct(121)*v(34)*v(52)
    a(122) = rct(122)*v(41)*v(52)
    a(123) = rct(123)*v(16)*v(54)
    a(124) = rct(124)*v(16)*v(57)

    ! aggregate function
    vdot(1) = 0.24_dp*a(62)+0.22_dp*a(64)+0.18_dp*a(65)+a(99)
    vdot(2) = a(45)
    vdot(3) = 0.52_dp*a(62)+0.22_dp*a(64)
    vdot(4) = 0.09_dp*a(64)+0.16_dp*a(65)+0.39_dp*a(80)+0.46_dp*a(83)+0.4_dp*a(116)
    vdot(5) = 0.6_dp*a(123)
    vdot(6) = a(122)
    vdot(7) = a(94)-a(95)
    vdot(8) = -a(70)
    vdot(9) = a(7)-a(9)-a(10)-a(11)
    vdot(10) = -a(8)-a(29)+a(30)+a(31)
    vdot(11) = -a(45)+0.4_dp*a(123)+a(124)
    vdot(12) = -a(71)
    vdot(13) = -a(46)+0.06_dp*a(64)+0.08_dp*a(65)
    vdot(14) = -a(47)+0.01_dp*a(64)+0.01_dp*a(65)
    vdot(15) = 0.4_dp*a(73)+a(74)-a(75)
    vdot(16) = -a(123)-a(124)
    vdot(17) = -a(5)-a(27)+a(33)-a(35)
    vdot(18) = -a(21)+0.08_dp*a(64)
    vdot(19) = 0.8_dp*a(70)+0.45_dp*a(71)-a(72)
    vdot(20) = -a(49)+0.03_dp*a(64)+0.04_dp*a(65)
    vdot(21) = -a(3)+a(22)-a(25)+a(34)
    vdot(22) = -a(86)-a(89)+a(113)
    vdot(23) = -a(87)-a(90)+a(114)
    vdot(24) = a(38)-a(41)-a(42)
    vdot(25) = -a(62)-a(63)
    vdot(26) = 0.12_dp*a(70)+0.05_dp*a(71)-a(73)-a(74)
    vdot(27) = a(1)+0.89_dp*a(2)+a(6)+a(9)+a(10)-a(12)-a(13)-a(14)-a(15)-a(16)
    vdot(28) = -a(44)+a(50)+a(51)+a(52)+a(53)+a(54)+a(59)+a(61)+0.24_dp*a(62) + &
               0.31_dp*a(64)+0.3_dp*a(65)+2_dp*a(76)+a(77)+0.69_dp*a(78) + &
               0.07_dp*a(80)+0.1_dp*a(83)+0.33_dp*a(84)+0.64_dp*a(85)+0.59_dp*a(104)
    vdot(29) = -a(4)+a(23)-a(26)+0.3_dp*a(40)+2_dp*a(41)+a(53)+a(56)+a(61)+a(74) + &
               0.07_dp*a(85)+a(124)
    vdot(30) = -a(48)-1.06_dp*a(64)-2.26_dp*a(65)-a(66)-2.23_dp*a(67)+1.1_dp*a(71) + &
               1.86_dp*a(85)-1.98_dp*a(88)+0.42_dp*a(91)-1.98_dp*a(93)-1.68_dp*a(98) - &
               a(101)+0.18_dp*a(102)+1.6_dp*a(103)-1.98_dp*a(108)-a(111)+2_dp*a(120)
    vdot(31) = 0.95_dp*a(72)+0.3_dp*a(73)-a(76)-a(77)-a(78)
    vdot(32) = a(81)-a(103)-a(120)
    vdot(33) = a(79)-a(102)-a(119)
    vdot(34) = 0.5_dp*a(82)-a(104)-a(121)
    vdot(35) = -a(11)+a(21)+a(28)-a(31)-a(41)
    vdot(36) = -a(57)-a(58)+0.07_dp*a(65)+0.23_dp*a(67)+0.09_dp*a(83)+0.03_dp*a(84) + &
               0.74_dp*a(88)+0.74_dp*a(93)+0.62_dp*a(98)+0.63_dp*a(104)+0.74_dp*a(108)
    vdot(37) = -a(65)-a(67)-a(69)
    vdot(38) = -a(79)-a(80)-a(81)
    vdot(39) = a(49)-a(50)-a(51)-a(52)-a(53)+a(62)+1.56_dp*a(63)+0.57_dp*a(64) + &
               a(66)+a(76)+0.7_dp*a(78)+0.6_dp*a(80)+0.15_dp*a(83)+0.2_dp*a(84) + &
               0.28_dp*a(85)+a(86)+0.3_dp*a(89)+a(96)+a(100)+0.5_dp*a(101) + &
               0.63_dp*a(102)+0.25_dp*a(104)+a(106)+a(110)+0.5_dp*a(111)
    vdot(40) = -a(64)-a(66)-a(68)
    vdot(41) = a(60)+a(63)+a(66)+a(67)+0.08_dp*a(70)+0.5_dp*a(71)+0.6_dp*a(73) + &
               a(76)+0.03_dp*a(78)+0.08_dp*a(79)+0.2_dp*a(80)+0.2_dp*a(82) + &
               0.07_dp*a(83)+0.93_dp*a(85)+0.4_dp*a(88)+0.41_dp*a(93)+0.34_dp*a(98) - &
               a(105)+0.4_dp*a(108)-a(112)-a(122)
    vdot(42) = -a(59)-a(60)-a(61)+0.04_dp*a(64)+0.07_dp*a(65)+0.8_dp*a(71) + &
               0.2_dp*a(78)+0.85_dp*a(83)+0.19_dp*a(91)+0.34_dp*a(104)
    vdot(43) = a(47)+0.06_dp*a(64)+0.05_dp*a(65)+0.1_dp*a(88)+0.7_dp*a(90) + &
               0.1_dp*a(93)-a(97)+0.08_dp*a(98)-a(107)+0.1_dp*a(108)-a(114)
    vdot(44) = a(68)+a(69)+a(92)-a(101)-a(111)-a(118)
    vdot(45) = -a(54)-a(55)-a(56)+0.22_dp*a(63)+0.47_dp*a(64)+1.03_dp*a(65)+a(66) + &
               1.77_dp*a(67)+0.03_dp*a(78)+0.15_dp*a(80)+0.02_dp*a(83)+0.07_dp*a(84) + &
               0.28_dp*a(85)+a(87)+0.3_dp*a(88)+0.3_dp*a(90)+0.04_dp*a(91) + &
               0.3_dp*a(93)+a(97)+0.25_dp*a(98)+0.5_dp*a(101)+0.8_dp*a(103) + &
               0.55_dp*a(104)+a(107)+0.3_dp*a(108)+0.5_dp*a(111)
    vdot(46) = a(46)+a(54)+a(57)+0.07_dp*a(64)+0.1_dp*a(65)+0.05_dp*a(83) + &
               0.7_dp*a(84)+0.7_dp*a(89)-a(96)+a(99)-a(106)+a(109)-a(113)
    vdot(47) = 0.65_dp*a(80)-a(82)-a(83)-a(84)-a(85)+0.91_dp*a(102)+0.2_dp*a(103)
    vdot(48) = a(58)+0.11_dp*a(65)-a(100)-a(110)-a(117)
    vdot(49) = -a(88)-a(91)+a(115)+a(117)+a(119)+a(121)
    vdot(50) = a(48)+0.03_dp*a(64)+0.09_dp*a(65)+0.77_dp*a(91)-a(98)-a(108)-a(115)
    vdot(51) = 0.05_dp*a(72)+a(75)+0.93_dp*a(85)-a(92)-a(93)+0.16_dp*a(98) + &
               0.5_dp*a(101)+0.09_dp*a(102)+0.8_dp*a(103)+0.5_dp*a(111)+a(118)+a(120)
    vdot(52) = a(5)+a(19)-a(20)+a(21)+a(24)-a(28)+a(29)-2_dp*a(30)-2_dp*a(31)-a(32) - &
               a(33)-a(34)+a(35)-a(40)+a(44)+a(45)+a(49)+2_dp*a(50)+a(52)+a(53) + &
               a(54)+a(59)+0.2_dp*a(62)+a(63)+0.26_dp*a(64)+0.22_dp*a(65)+a(66) + &
               a(67)+0.2_dp*a(70)+0.55_dp*a(71)+0.95_dp*a(72)+0.6_dp*a(73) + &
               2_dp*a(76)+a(77)+0.76_dp*a(78)+0.07_dp*a(80)+0.1_dp*a(83) + &
               0.33_dp*a(84)+0.93_dp*a(85)+a(86)+a(87)+0.9_dp*a(88)+0.9_dp*a(93) + &
               a(96)+a(97)+0.76_dp*a(98)+0.5_dp*a(101)+0.91_dp*a(102)+0.8_dp*a(103) + &
               a(104)+a(106)+a(107)+0.9_dp*a(108)+0.5_dp*a(111)-a(113)-a(114)-a(115) - &
               a(116)-a(117)-a(118)-a(119)-a(120)-a(121)-a(122)
    vdot(53) = -a(6)-a(7)+a(12)-a(13)-a(17)-a(18)-a(19)-a(20)-a(62)-a(64)-a(65) - &
               a(78)-a(80)-a(83)+0.4_dp*a(116)
    vdot(54) = a(3)+a(4)+2_dp*a(8)+2_dp*a(11)-a(19)+a(20)-a(21)-a(22)-a(23)-a(24) - &
               a(25)-a(26)-a(27)-a(28)-a(29)+a(32)+0.7_dp*a(40)-a(44)-a(45)-a(46) - &
               a(47)-a(48)-a(49)-a(52)-a(55)-a(58)-a(60)+0.12_dp*a(62)-a(63) + &
               0.33_dp*a(64)+0.6_dp*a(65)-a(66)-a(67)-a(70)-a(71)-a(73)-a(76) + &
               0.08_dp*a(78)-a(79)+0.27_dp*a(80)-a(82)+0.27_dp*a(83)+a(86)+a(87) + &
               a(88)-0.7_dp*a(89)-0.7_dp*a(90)-0.77_dp*a(91)-a(92)-a(123)
    vdot(55) = a(1)+0.11_dp*a(2)+a(3)+a(14)-a(16)-a(17)-a(22)-a(32)-a(36)+a(37) - &
               2_dp*a(43)-a(72)-a(96)-a(97)-a(98)-a(99)-a(100)-a(101)-a(102) - &
               a(103)-a(104)-a(105)
    vdot(56) = -a(1)+0.89_dp*a(2)+a(4)+a(5)-a(14)-a(15)+a(16)+a(17)-a(18)-a(23) + &
               a(24)+a(25)+a(27)+a(32)-a(33)-a(34)+a(35)+2_dp*a(36)-a(38)+2_dp*a(39) + &
               0.7_dp*a(40)+a(42)+2_dp*a(43)+0.95_dp*a(72)-a(75)+a(93)-a(94)+a(95) + &
               a(96)+a(97)+0.84_dp*a(98)+a(99)+a(100)+1.5_dp*a(101)+0.91_dp*a(102) + &
               1.2_dp*a(103)+a(104)+a(105)+a(106)+a(107)+a(108)+a(109)+a(110) + &
               1.5_dp*a(111)+a(112)
    vdot(57) = -a(2)+a(15)+a(18)-a(24)+a(26)-a(36)-a(37)-a(38)-2_dp*a(39)-a(40) + &
               a(42)-a(53)-a(56)-a(61)-a(68)-a(69)-a(74)-a(81)-a(85)-a(106)-a(107) - &
               a(108)-a(109)-a(110)-a(111)-a(112)-a(124)
    vdot(58) = a(55)+a(56)+a(57)+a(59)+a(60)+a(61)+0.13_dp*a(64)+0.19_dp*a(65)+a(76) + &
               a(77)+0.62_dp*a(78)+0.2_dp*a(80)+0.5_dp*a(82)+0.11_dp*a(83) + &
               0.97_dp*a(84)+0.07_dp*a(85)-a(94)+a(95)-a(99)+a(100)-a(109)+a(110)-a(116)

  end subroutine fun

end module mod_cbmz_function

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
