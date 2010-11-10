* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* script_math_demo.gs
*
* Illustrates the use of the scripting language math funcitons
* Written by M. Fiorino
* Modified by J.M. Adams
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

doformat = 1
dotrig = 1
doint = 1
domisc = 1
dovalnum = 1

* FORMATTING NUMBERS
if (doformat)
  print ' '
  print 'Formatting numbers:'
  v = 3.1456
  fmt = '%-6.1f'
  rc = math_format(fmt,v)
  print fmt' of 'v' = 'rc

  fmt = '%-6.2f'
  rc = math_format(fmt,v)
  print fmt' of 'v' = 'rc

  fmt = '%6.2f'
  rc = math_format(fmt,v)
  print fmt'  of 'v' = 'rc

  fmt = '%06.3f'
  rc = math_format(fmt,v)
  print fmt' of 'v' = 'rc

  fmt = '%6.0g'
  rc = math_format(fmt,v)
  print fmt'  of 'v' = 'rc

  fmt = '%12.3e'
  rc = math_format(fmt,v)
  print fmt' of 'v' = 'rc
  print ' '
endif

* TRIG FUNCTIONS
if (dotrig)
  print 'Trigonometric functions:'
  pi  = 3.1415926
  d2r = pi/180
  r2d = 1/d2r
  dang = 45
  ang = dang*d2r
  fmt = '%-6.4f'

  cos  = math_cos(ang)
  cosf = math_format(fmt,cos)
  print 'cos of 'dang' = 'cosf

  sin  = math_sin(ang)
  sinf = math_format(fmt,sin)
  print 'sin of 'dang' = 'sinf

  tan  = math_tan(ang)
  tanf = math_format(fmt,tan)
  print 'tan of 'dang' = 'tanf

  ang = math_acos(cos)
  dang = ang * r2d
  dangf = math_format(fmt,dang)
  print 'acos of 'cosf' = 'dangf

  ang = math_asin(sin)
  dang = ang * r2d
  dangf = math_format(fmt,dang)
  print 'asin of 'sinf' = 'dangf

  ang = math_atan(tan)
  dang = ang * r2d
  dangf = math_format(fmt,dang)
  print 'atan of 'tanf' = 'dangf

  u =  1
  v = -1
  rc = math_atan2(u,v)
  rc2 = rc*r2d
  rc2f = math_format(fmt,rc2)
  print 'atan2 of ('u','v') = 'rc2f

  dang = 30
  ang = dang*d2r

  cosh = math_cosh(ang)
  coshf = math_format(fmt,cosh)
  print 'cosh of 'dang' = 'coshf

  sinh = math_sinh(ang)
  sinhf = math_format(fmt,sinh)
  print 'sinh of 'dang' = 'sinhf

  tanh = math_tanh(ang)
  tanhf = math_format(fmt,tanh)
  print 'tanh of 'dang' = 'tanhf

  demo = dang
  ang = math_acosh(cosh)
  dang = ang * r2d
  dangf = math_format(fmt,dang)
  print 'acosh of 'coshf' = 'dangf' <-- math_acosh() may have a bug ... this value should be 'demo

  ang = math_asinh(sinh)
  dang = ang * r2d
  dangf = math_format(fmt,dang)
  print 'asinh of 'sinhf' = 'dangf

  ang = math_atanh(tanh)
  dang = ang * r2d
  dangf = math_format(fmt,dang)
  print 'atanh of 'tanhf' = 'dangf
  print ' '
endif

* VALNUM
if (dovalnum)
  print 'Evaluating strings to see if they are numbers:'
  print '0 = not a number'
  print '1 = integer'
  print '2 = not an integer'
  num = '3.1455'
  rc = valnum(num)
  print 'valnum of 'num' = 'rc

  num = '3'
  rc = valnum(num)
  print 'valnum of 'num' = 'rc

  num = '3e10'
  rc = valnum(num)
  print 'valnum of 'num' = 'rc

  num = '390210'
  rc = valnum(num)
  print 'valnum of 'num' = 'rc

  num = 'a string'
  rc = valnum(num)
  print 'valnum of 'num' = 'rc
  print ' '
endif

* INTEGER CONVERSION
if (doint) 
  print 'Real-to-integer conversion:'
  v = 3.0
  while(v < 4.0) 
    rc1 = math_nint(v)
    rc2 = math_int(v)
    print 'nint of 'v' = 'rc1'    int of 'v' = 'rc2
    v = v + 0.1
  endwhile
  print ' '
endif

* MISCELLANEOUS FUNCTIONS
if (domisc)
  print 'Exponents:'
  pow = math_pow(2,0.5);
  print '2 raised to the power 0.5 = 'pow
  print ' '

  print 'Exponential function:'
  num = math_exp(1)
  print 'exp(1) = 'num
  print ' '

  print 'Modulo operator:'
  fmod = math_fmod(5,2);
  print '5 modulo 2 (the remainder when 5 is divided by 2) = 'fmod
  print ' '

  print 'String operations:'
  s = 'this is a test'
  rc = math_strlen(s)
  print 'length of the string "'s'" = 'rc

  p = 2
  rc = wrdpos(s,p)
  print 'word 'p' of the string "'s'" starts at character 'rc
  print ' '
endif


return