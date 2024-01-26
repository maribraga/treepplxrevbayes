#!/bin/bash
tmux new-session -d -s "tiny_mu1" ./run_inf.sh hostrep_inf_tiny_mu.Rev 1;sleep 0.1
tmux new-session -d -s "tiny_mu2" ./run_inf.sh hostrep_inf_tiny_mu.Rev 2;sleep 0.1
