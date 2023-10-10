# single-molecule-clock-simulator
Kinetic Monte Carlo simulation of single-molecule clocks as described in Johnson-Buck and Shih, Nano Lett. 2017, 17, 12, 7940–7944.

# Detailed Description

Simulates single-molecule clocks comprising a linear reaction pathway
containing an arbitrary number of steps (N).  The steps can be irreversible
(k2 = 0) or reversible (k2 > 0).
As described in Johnson-Buck and Shih, Nano Lett. 2017, 17, 12, 7940–7944.

The Gillespie algorithm for kinetic Monte Carlo simulation is used.

Input arguments (all optional, but must be specified in the indicated order):

results = simulate_sm_clock(N, trials, k1, k2)

N = Number of steps in pathway (default = 50)
trials = Number of trajectories to simulate (default = 100)
k1 = Forward rate constant, in seconds (default = 1)
k2 = Reverse rate constant, in seconds (default = 2)

# Author Information
Copyright (C) 2016 Alex Johnson-Buck
