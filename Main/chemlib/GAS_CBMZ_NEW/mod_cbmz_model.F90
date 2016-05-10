module mod_cbmz_model

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !  Completely defines the model mod_cbmz
  !    by using all the associated modules
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  use mod_cbmz_precision
  use mod_cbmz_parameters
  use mod_cbmz_global
  use mod_cbmz_function
  use mod_cbmz_integrator
  use mod_cbmz_rates
  use mod_cbmz_jacobian
  use mod_cbmz_linearalgebra
  use mod_cbmz_monitor

end module mod_cbmz_model

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
