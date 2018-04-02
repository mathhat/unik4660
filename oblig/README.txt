You might need to install hdf5 libs first, but writing "make" will compile the stuff for you. I'll leave the executable "read" so you won't be forced to do that.
The python script for visualizing streamlines is currently dealing with too dense of a data file to visualize. This might be fixed.

Interpolation is done at the top of the read.cpp script.
Changing line lenght is done by changing the L variable at line 106
Changing output filename for LIC images is done in the pgm function :) 
