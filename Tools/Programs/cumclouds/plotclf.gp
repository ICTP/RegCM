
set key left top
set terminal postscript color
set output 'Cloud_fraction.ps'
set xlabel 'Depth of cloud [mb]'
set ylabel 'Total cloud fraction [%]'
set title  'RegCM cumulus cloud model'
plot 'dat.dat' u ($2*10):($3*100)   title '10 km grid', \
     'dat.dat' u ($6*10):($7*100)   title '25 km grid', \
     'dat.dat' u ($10*10):($11*100)   title '50 km grid', \
     'dat.dat' u ($14*10):($15*100) title '75 km grid', \
     'dat.dat' u ($18*10):($19*100) title '100km grid', \
     'dat.dat' u ($22*10):($23*100) title '125km grid', \
     'dat.dat' u ($26*10):($27*100) title '150km grid', \
     'dat.dat' u ($30*10):($31*100) title '175km grid', \
     'dat.dat' u ($34*10):($35*100) title '200km grid', \
     'dat.dat' u ($10*10):($12*100) title 'old 50 km grid'
