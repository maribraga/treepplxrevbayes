#!/bin/bash
tmux new-session -d -s "long1" ./run_inf.sh hostrep_inf_long.Rev 1;sleep 0.1
tmux new-session -d -s "long2" ./run_inf.sh hostrep_inf_long.Rev 2;sleep 0.1
