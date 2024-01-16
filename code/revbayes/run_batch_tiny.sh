#!/bin/bash
tmux new-session -d -s "tiny1" ./run_inf.sh hostrep_inf_tiny.Rev 1;sleep 0.1
tmux new-session -d -s "tiny2" ./run_inf.sh hostrep_inf_tiny.Rev 2;sleep 0.1
