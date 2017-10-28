#!/bin/tcsh
set VISCODIR = /Users/fpollitz/fred4/CIDER2015/VISCO2.5D
set EXDIR = $VISCODIR/MAINPROG
set WORKDIR = $VISCODIR/Work
set WORKDIR1 = $VISCODIR/Work1
set WORKDIR2 = $VISCODIR/Work2
#
\cp elastic.paramEX1-Maxwell $WORKDIR/elastic.param
\cp simulation-spherical.infoEX1 $WORKDIR/simulation-spherical.info
\cp source-spherical.paramEX1 $WORKDIR/source-spherical.param
\cp receivers-latlondepEX1.txt $WORKDIR/receivers-latlondep.txt
\cp $EXDIR/visco2pt5d $WORKDIR/.

cd $WORKDIR
\rm sporder.out
./visco2pt5d << ! > /dev/null
1
!

date > $WORKDIR1/visco2pt5d-progress1.txt
date > $WORKDIR2/visco2pt5d-progress2.txt

\cp $WORKDIR/elastic.param $WORKDIR1/.
\cp $WORKDIR/simulation-spherical.info $WORKDIR1/.
\cp $WORKDIR/source-spherical.param $WORKDIR1/.
\cp $EXDIR/visco2pt5d $WORKDIR1/.
\cp $WORKDIR/elastic.param $WORKDIR2/.
\cp $WORKDIR/simulation-spherical.info $WORKDIR2/.
\cp $WORKDIR/source-spherical.param $WORKDIR2/.
\cp $EXDIR/visco2pt5d $WORKDIR2/.

cd $WORKDIR1
ln -sf $WORKDIR/sporder.out sporder.out
\cp $WORKDIR/receivers-latlondep.txt receivers-latlondep.txt
./visco2pt5d << ! > /dev/null &
0
rec-j-m1
displ-j-m1
visco2pt5d-progress1.txt
0. 0.50
!

cd $WORKDIR2
ln -sf $WORKDIR/sporder.out sporder.out
\cp $WORKDIR/receivers-latlondep.txt receivers-latlondep.txt
./visco2pt5d << ! > /dev/null &
0
rec-j-m2
displ-j-m2
visco2pt5d-progress2.txt
0.50 1.00
!

wait

cd $WORKDIR
\cat $WORKDIR1/vertp-j $WORKDIR2/vertp-j > vertp-j
\cat $WORKDIR1/rec-j-m1 $WORKDIR2/rec-j-m2 > rec-j-m
\cat $WORKDIR1/displ-j-m1 $WORKDIR2/displ-j-m2 > displ-j-m

visco2pt5d << ! > /dev/null 
2
not-needed
not-needed
not-needed
0. 1.00
!
#
\mv visco2pt5d-stat_vertp.gmt visco2pt5d-stat_vertp-EX1.gmt
\mv visco2pt5d-statstrains_vertp.gmt visco2pt5d-statstrains_vertp-EX1.gmt
\mv visco2pt5d-stat-rec.gmt visco2pt5d-stat-rec-EX1.gmt
\mv visco2pt5d-statstrains-rec.gmt visco2pt5d-statstrains-rec-EX1.gmt
\mv visco2pt5d-post_vertp.gmt visco2pt5d-post_vertp-EX1.gmt
\mv visco2pt5d-poststrains_vertp.gmt visco2pt5d-poststrains_vertp-EX1.gmt
\mv visco2pt5d-post-rec.gmt visco2pt5d-post-rec-EX1.gmt
\mv visco2pt5d-poststrains-rec.gmt visco2pt5d-poststrains-rec-EX1.gmt
\mv visco2pt5d-stat-phival.gmt visco2pt5d-stat-phival-EX1.gmt
\mv visco2pt5d-post-phival.gmt visco2pt5d-post-phival-EX1.gmt
