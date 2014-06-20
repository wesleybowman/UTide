#!/bin/bash

~/MATLAB/bin/matlab -nosplash -nodisplay -nodesktop -r "run UTideCurrentVersion/Utide_test.m; exit;"
#cp UTideCurrentVersion/matlabcoef.mat .
#cp UTideCurrentVersion/speedmatlabcoef.mat .
#cp UTideCurrentVersion/matlabrecon .
#cp UTideCurrentVersion/speedmatlabrecon .
python Utide_test.py
python comparison.py
python reconCompare.py
