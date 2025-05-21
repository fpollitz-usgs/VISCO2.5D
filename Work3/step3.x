#!/bin/sh

ln -sf ../EXAMPLES-VISCO3D/sporder.out
\cp ../EXAMPLES-VISCO3D/receivers-latlondep.txt receivers-latlondep.txt
./visco3d << ! > look
0
rec-j-m3
displ-j-m3
visco3d-progress3.txt
0.28572 0.42858
!
