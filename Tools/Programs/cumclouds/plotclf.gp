
set key left top
set terminal postscript color
set output 'Cloud_fraction.ps'
set xlabel 'Depth of cloud [mb]'
set ylabel 'Total cloud fraction [%]'
set title  'RegCM cumulus cloud model'
plot 'dat.dat' u ($2*10):($3*100)   title '10 km grid', \
     'dat.dat' u ($5*10):($6*100)   title '25 km grid', \
     'dat.dat' u ($8*10):($9*100)   title '50 km grid', \
     'dat.dat' u ($11*10):($12*100) title '75 km grid', \
     'dat.dat' u ($14*10):($15*100) title '100km grid', \
     'dat.dat' u ($17*10):($18*100) title '125km grid', \
     'dat.dat' u ($20*10):($21*100) title '150km grid', \
     'dat.dat' u ($23*10):($24*100) title '175km grid', \
     'dat.dat' u ($26*10):($27*100) title '200km grid'
