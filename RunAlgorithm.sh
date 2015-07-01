#!/usr/bin/env bash

NUM_MOLS=9
NUM_REACTIONS=4
CAT_PROB=0.05
F1="\"N\""
F2="\"O\""
F3="\"C=O\""
F4="\"C#N\""
F5="\"C#CC#N\""

(echo "$NUM_MOLS"; echo "$NUM_REACTIONS"; echo "$CAT_PROB"; echo "$F1"; echo "$F2"; echo "$F3"; echo "$F4"; echo "$F5" ; echo "") |  python Compute_maxRAF.py

echo `python Resave.py`

echo `dot -Teps -O RAFhighlight.txt`
