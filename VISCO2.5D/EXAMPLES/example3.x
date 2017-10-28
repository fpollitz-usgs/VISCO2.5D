#!/bin/bash
#
export VISCODIR=/home/cig/VISCO2.5D
export EXDIR=$VISCODIR/MAINPROG
export WORKDIR=$VISCODIR/Work
export WORKDIR1=$VISCODIR/Work1
#
\cp elastic.paramEX3 $WORKDIR/elastic.param
\cp simulation-spherical.infoEX3 $WORKDIR/simulation-spherical.info
\cp source-spherical.paramEX3 $WORKDIR/source-spherical.param
\cp receivers-latlondepEX3.txt $WORKDIR/receivers-latlondep.txt
\cp $EXDIR/visco2pt5d $WORKDIR/.

cd $WORKDIR
\rm sporder.out
./visco2pt5d << ! > /dev/null
1
!

date > $WORKDIR1/visco2pt5d-progress1.txt

\cp $WORKDIR/elastic.param $WORKDIR1/.
\cp $WORKDIR/simulation-spherical.info $WORKDIR1/.
\cp $WORKDIR/source-spherical.param $WORKDIR1/.
\cp $EXDIR/visco2pt5d $WORKDIR1/.

cd $WORKDIR1
ln -sf $WORKDIR/sporder.out sporder.out
\cp $WORKDIR/receivers-latlondep.txt receivers-latlondep.txt
./visco2pt5d << ! > /dev/null
0
rec-j-m1
displ-j-m1
visco2pt5d-progress1.txt
0. 1.00
!

cd $WORKDIR
\cp $WORKDIR1/vertp-j vertp-j
\cp $WORKDIR1/rec-j-m1 rec-j-m
\cp $WORKDIR1/displ-j-m1  displ-j-m

./visco2pt5d << ! > /dev/null 
2
not-needed
not-needed
not-needed
0. 1.00
!
#
\mv visco2pt5d-stat_vertp.gmt visco2pt5d-stat_vertp-EX3.gmt
\mv visco2pt5d-stat-rec.gmt visco2pt5d-stat-rec-EX3.gmt
\mv visco2pt5d-post_vertp.gmt visco2pt5d-post_vertp-EX3.gmt
\mv visco2pt5d-post-rec.gmt visco2pt5d-post-rec-EX3.gmt
\mv visco2pt5d-stat-phival.gmt visco2pt5d-stat-phival-EX3.gmt
\mv visco2pt5d-post-phival.gmt visco2pt5d-post-phival-EX3.gmt
