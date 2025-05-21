#!/bin/sh

ln -sf ../EXAMPLES-VISCO3D/sporder.out
\cp ../EXAMPLES-VISCO3D/receivers-latlondep.txt receivers-latlondep.txt
./visco3d << ! > look
0
rec-j-m1
displ-j-m1
visco3d-progress1.txt
0. 0.14286
!
