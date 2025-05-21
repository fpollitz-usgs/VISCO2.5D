#!/bin/sh

ln -sf ../EXAMPLES-VISCO3D/sporder.out
\cp ../EXAMPLES-VISCO3D/receivers-latlondep.txt receivers-latlondep.txt
./visco3d << ! > look
0
rec-j-m5
displ-j-m5
visco3d-progress5.txt
0.57143 0.71429
!
