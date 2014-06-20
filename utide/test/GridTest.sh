#!/bin/bash

#~/MATLAB/bin/matlab -nosplash -nodisplay -nodesktop -r "run grid_test.m; exit;"
~/MATLAB/bin/matlab -nosplash -nodisplay -nodesktop -r "run UTideCurrentVersion/grid_test.m; exit;"
#cp UTideCurrentVersion/matlabcoef.mat .
python grid_test.py
python comparison.py
python reconCompare.py
