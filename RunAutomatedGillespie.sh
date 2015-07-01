#!/usr/bin/env bash

F1="\"N\""
F2="\"O\""
F3="\"C=O\""
F4="\"C#N\""
F5="\"C#CC#N\"" 

echo `python Stmatrix_parser.py`

(echo "$F1"; echo "$F2"; echo "$F3"; echo "$F4"; echo "$F5" ; echo "") |  python Find_Food.py

echo `python Find_NonFood.py`

echo `python Find_Catalysts.py`

echo `../../../../Matlab/MATLAB_R2012b.app/bin/./matlab -r AutomatedGillespie`
