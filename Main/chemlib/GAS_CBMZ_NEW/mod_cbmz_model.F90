module mod_cbmz_model

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !  Completely defines the model mod_cbmz
  !    by using all the associated modules
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  use mod_cbmz_precision, only : dp, ik4, rk8
  use mod_cbmz_parameters, only : nvar
  use mod_cbmz_global, only : altabove, altbelow, altmid, atol, c, c_jval, &
    c_m, cfactor, dens, deptha, depthb, dt, fix, jval_c2h5cho, jval_ch2oa, &
    jval_ch2ob, jval_ch3choa, jval_ch3coch3, jval_chocho, jval_h2o2, &
    jval_hno2, jval_hno3, jval_hno4, jval_n2o5b, jval_no2, jval_no3a, &
    jval_no3b, jval_o31d, jval_o33p, jval_pan, rconst, rtol, stepmax, &
    stepmin, temp, tend, time, tstart, var, xr, xrin, xrout, zenith
  use mod_cbmz_integrator, only : integrate, reset_integrate
  use mod_cbmz_rates, only : loadperoxyparameters, update_rconst
  implicit none

end module mod_cbmz_model

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
