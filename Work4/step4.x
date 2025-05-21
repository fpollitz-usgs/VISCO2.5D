#!/bin/sh

ln -sf ../EXAMPLES-VISCO3D/sporder.out
\cp ../EXAMPLES-VISCO3D/receivers-latlondep.txt receivers-latlondep.txt
./visco3d << ! > look
0
rec-j-m4
displ-j-m4
visco3d-progress4.txt
0.42858 0.57143
!
