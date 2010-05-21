*****************************************************************
*         HORIZONTAL FORWARD AND BACKWARD TRAJECTORIES          *
*           Department of Astronomy and Meteorology             *
*                    University of Barcelona                    *
*                                                               *
*               (C) Bernat Codina - spring 2001                 * 
*                      bcodina@am.ub.es                         *
*****************************************************************

  function main(arg)
    _At=6371229
    _PI=3.141592654
    _D2R=_PI/180
    _R2D=180/_PI

* Wind components in the ctl file. Must be changed if other
* notation is employed.
    _U='U'
    _V='V'

    if(inici())
      return
    endif
    if (arg='')
      say 'Click on the map to start a trajectory.'
      say 'Alternatively, you can also issue:'
      say '        traj lon_ini lat_ini'
      say ' '
      'q pos'
      x=subwrd(result,3)
      y=subwrd(result,4)
      'q xy2w 'x' 'y
      LON0=subwrd(result,3)
      LAT0=subwrd(result,6)
    else
      LON0=subwrd(arg,1)
      LAT0=subwrd(arg,2)
    endif  
    volta=1
    while (volta<=2)
      if (volta=1)
        say 'Forward trajectory (in red):'
        signe=1
        _color=2
        Tm=-999999
        TM=_Tmax
      else
        say 'Backward trajectory (in blue):'
        signe=-1
        _color=4
        TM=999999
        Tm=1
      endif
      LONt=LON0
      LATt=LAT0
      t=_Tini
      say t' 'LONt' 'LATt
      'set T ' t
      Ut=interp(_U,LONt,LATt)  
      Vt=interp(_V,LONt,LATt)
      Utm1=Ut
      Vtm1=Vt
      while (t<TM & t>Tm)
        Ut=Utm1
        Vt=Vtm1
        'set T ' t+signe
        LON0tm1=LONt
        LAT0tm1=LATt
        DR=1000000
        iter=0
        while(DR>0.1 & iter<15)
          Umig=signe*(Utm1+Ut)/2
          Vmig=signe*(Vtm1+Vt)/2
          'd sqrt('Umig'*'Umig'+'Vmig'*'Vmig')'
          V=subwrd(result,4)
          D=V*_DT
          'd atan2('Vmig','Umig')'
          alpha=subwrd(result,4)*_R2D
          rec=mou(LONt,LATt,D,alpha)
          LON1tm1=subwrd(rec,1)
          LAT1tm1=subwrd(rec,2)
          DR=dist(LON1tm1,LAT1tm1,LON0tm1,LAT0tm1)
          Utm1=interp(_U,LON1tm1,LAT1tm1)
          Vtm1=interp(_V,LON1tm1,LAT1tm1)
          LON0tm1=LON1tm1
          LAT0tm1=LAT1tm1
          iter=iter+1
        endwhile
        if (iter=15)
          say '--> Lack of convergence: DR='DR' m.'
        endif
        say t+signe' 'LON1tm1' 'LAT1tm1
        linia(LONt,LATt,LON1tm1,LAT1tm1)
        LONt=LON1tm1
        LATt=LAT1tm1
        t=t+signe
      endwhile
      volta=volta+1
    endwhile
    fi(0)
  return

  function mou(lon0,lat0,dist,alpha)
* Coordinates of a point located a distance dist (in meters) and angle
* alpha (in degrees and mathematical convention) from the point (lon0,lat0).
    if (dist<1);return(lon0' 'lat0);endif
    lon0=lon0*_D2R
    lat0=lat0*_D2R
    alpha=90-alpha
    if (alpha<0)
      alpha=360+alpha
    endif
    if (alpha>=360)
      alpha=alpha-360
    endif
    A=alpha*_D2R
    b=dist/_At
    c=_PI/2-lat0
    'd acos(cos('b')*cos('c')+sin('b')*sin('c')*cos('A'))'
    a=subwrd(result,4)
    lat1=(_PI/2-a)*_R2D
    'd asin(sin('b')*sin('A')/sin('a'))'
    B=subwrd(result,4)
    lon1=(lon0+B)*_R2D
  return(lon1' 'lat1)

  function interp(camp,lon,lat)
* Interpolates variable camp in (lon,lat).
    'q w2gr 'lon' 'lat
    x=subwrd(result,3);y=subwrd(result,6)
    x.2=int(x);y.3=int(y)
    x.3=x.2;y.4=y.3
    x.1=x.2+1;y.1=y.3+1
    x.4=x.1;y.2=y.1
    num=0;den=0
    i=1
    while(i<=4)
      if(x.i<_Xmin|x.i>_Xmax|y.i<_Ymin|y.i>_Ymax)
        fi(1)
      endif
      'set x 'x.i
      'set y 'y.i
      'q gr2w 'x.i' 'y.i
      rec=sublin(result,1)
      lon.i=subwrd(rec,3)
      lat.i=subwrd(rec,6)
      'd 'camp
      rec=sublin(result,1)
      if(subwrd(rec,1)!='Result')
        rec=sublin(result,2)
      endif
      z.i=subwrd(rec,4)
      d.i=dist(lon,lat,lon.i,lat.i)
      if (d.i=0)
        return (z.i)
        break
      endif
      num=num+1/d.i*z.i
      den=den+1/d.i
      i=i+1
    endwhile
  return (num/den)

  function dist(lon1,lat1,lon2,lat2)
* Distance between two points on the Earth surface
    phi=lon1*_D2R;theta=(90-lat1)*_D2R
    'd sin('theta')*cos('phi')'
    x1=subwrd(result,4)
    'd sin('theta')*sin('phi')'
    y1=subwrd(result,4)
    'd cos('theta')'
    z1=subwrd(result,4)
    phi=lon2*_D2R;theta=(90-lat2)*_D2R
    'd sin('theta')*cos('phi')'
    x2=subwrd(result,4)
    'd sin('theta')*sin('phi')'
    y2=subwrd(result,4)
    'd cos('theta')'
    z2=subwrd(result,4)
    x=y1*z2-y2*z1
    y=x2*z1-x1*z2
    z=x1*y2-x2*y1
    d2=x*x+y*y+z*z
    'd asin(sqrt('d2'))'
  return (subwrd(result,4)*_At)

  function int(x)
* Integer part of x
    i=1
    while(i<17)
      if(substr(x,i,1)=".");break;endif
      i=i+1
    endwhile
    if (i=17); return(x);endif
    if (i=1); return('0');endif
  return (substr(x,1,i-1))

  function linia(lon0,lat0,lon1,lat1)
* Draws a line limited by two circles in _color
    'set line '_color
    'q w2xy 'lon0' 'lat0
    x0=subwrd(result,3)
    y0=subwrd(result,6)
    'q w2xy 'lon1' 'lat1
    x1=subwrd(result,3)
    y1=subwrd(result,6)
    'draw line 'x0' 'y0' 'x1' 'y1
    'draw mark 3 'x0' 'y0' .07'
    'draw mark 3 'x1' 'y1' .07'
  return
  
  function inici
* Some verifications and variable initialization
  'q w2gr 1 1'
  if(subwrd(result,3)='environment')
    say 'Please open first a file and draw a chart.'
    return(1)
  endif
  'q dims'
  rec1=sublin(result,2)
  rec2=sublin(result,3)
  rec3=sublin(result,4)
  rec4=sublin(result,5)
  if(subwrd(rec1,3)='fixed' | subwrd(rec2,3)='fixed')
    say 'X and Y dimensions must be opened:'
    say rec1
    say rec2
    return(1)
  endif
  if(subwrd(rec3,3)!='fixed')
    say 'Vertical dimension must be fixed:'
    say rec3
    return(1)
  endif
  if(subwrd(rec4,3)!='fixed')
    say 'Time dimension must be fixed'
    say rec4
    return(1)
  endif
  _Xmin=subwrd(rec1,11)
  _Xmax=subwrd(rec1,13)
  _Ymin=subwrd(rec2,11)
  _Ymax=subwrd(rec2,13)
  _Tini=subwrd(rec4,9)
  'q file'
  rec=sublin(result,5)
  _Tmax=subwrd(rec,12)
  rec=sublin(result,2)
  fitxer=subwrd(rec,2)
  rc=0
  while(rc=0)
    rec=read(fitxer)
    rc=sublin(rec,1)
    lin=sublin(rec,2)
    if (subwrd(lin,1)='tdef' | subwrd(lin,1)='TDEF')
      interval=subwrd(lin,5)
      break
    endif
  endwhile
  if (rc!=0)
    say 'TDEF not found in 'fitxer
    return(1)
  endif
  rc=close(fitxer)
  i=1
  factor=0
  while(i<10)
    a=substr(interval,i,2)
    if(a='mn'|a='MN'|a='Mn')
      factor=60
      break
    endif
    if(a='hr'|a='HR'|a='Hr')
      factor=3600
      break
    endif
    i=i+1
  endwhile
  if(factor>0)
    _DT=factor*substr(interval,1,i-1)
  else
    say 'HR or MN missing in line TDEF of file 'fitxer
    return(1)
  endif
  return (0)

  function fi(i)
* Things are left as they were before execution of traj
    'set x ' _Xmin' '_Xmax
    'set y ' _Ymin' '_Ymax
    'set t ' _Tini
* Does anybody know how to nicely quit a GrADS script?
    if (i=1)
      say 'The trajectory can't be completed.'
      say 'Press Ctrl+C and Enter to quit.'
      pull nothing
      fi(1)
    endif
  return
