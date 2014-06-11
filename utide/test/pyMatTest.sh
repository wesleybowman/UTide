#!/bin/bash

~/MATLAB/bin/matlab -nosplash -nodisplay -nodesktop -r "run UTideCurrentVersion/Utide_test.m; exit;"
cp UTideCurrentVersion/matlabcoef.mat .
python Utide_test.py
python comparison.py
