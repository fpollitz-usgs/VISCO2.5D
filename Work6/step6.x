#!/bin/sh

ln -sf ../EXAMPLES-VISCO3D/sporder.out
\cp ../EXAMPLES-VISCO3D/receivers-latlondep.txt receivers-latlondep.txt
./visco3d << ! > look
0
rec-j-m6
displ-j-m6
visco3d-progress6.txt
0.71429 0.85715
!
