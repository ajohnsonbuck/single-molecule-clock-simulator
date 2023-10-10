function results = simulate_sm_clock(varargin)

% Copyright (C) 2016 Alex Johnson-Buck
%
% Simulates single-molecule clocks comprising a linear reaction pathway
% containing an arbitrary number of steps (N).  The steps can be irreversible
% (k2 = 0) or reversible (k2 > 0).
% As described in Johnson-Buck and Shih, Nano Lett. 2017, 17, 12, 7940–7944.

% The Gillespie algorithm for kinetic Monte Carlo simulation is used.

% Input arguments (all optional, but must be specified in the indicated order):
%
% results = simulate_sm_clock(N, trials, k1, k2)
%
% N = Number of steps in pathway (default = 50)
% trials = Number of trajectories to simulate (default = 100)
% k1 = Forward rate constant, in seconds (default = 1)
% k2 = Reverse rate constant, in seconds (default = 2)

args = varargin

[N,trials,k1,k2] = parse_inputs(args);

display_parameters(N,trials,k1,k2);

##kmat = cat(1,ones(1,N)*k1,ones(1,N)*k2); %Matrix of rate constants for all states; row 1 is forward, row 2 is reverse; columns represent different states
##kmat(2,1) = 0; % Backward rate constant for first state (state 0) is 0.  There is no state (-1).

ksum = k1+k2; % sum of rate constants for all states other than A = 0

DT = zeros(trials,1); % Initialize matrix of dwell times for reaching final state

for trial = 1:trials
    dt1 = single(random('exp',1/k1,10000,1)); % Initialize randomized dwell times
    dt2 = single(random('exp',1/ksum,10000,1));
    rxnp = random('unif',0,1,10000,1); % Initialize random reaction variable
    traj = zeros(10000,1); % initialize state vector
    time = zeros(10000,1); % initialize time vector
    pi = cat(2,k1/ksum, k2/ksum); % transition probability matrix for state A > 0
    pmat = cumsum(pi);
    A = 0; % Start in state 0
    t = 0; % initialize time
    s = 0; % initialize step

    if mod(trial,50)==0
       disp(num2str(trial));
    end

    while A ~= N
    s = s+1; % Update step
    if mod(s,10000)==0
       dt1 = cat(1,dt1,single(random('exp',1/k1,10000,1))); % add to existing random variable matrices if needed
       dt2 = cat(1,dt2,single(random('exp',1/ksum,10000,1)));
       rxnp = cat(1,rxnp,random('unif',0,1,10000,1));
       traj = cat(1,traj,zeros(10000,1));
       time = cat(1,time,zeros(10000,1));
    end
    if A > 0
      dt = dt2(s); % Time interval
      if pmat(1) > rxnp(s)
        A = A+1;
      else
        A = A-1;
      end
    else % if A = 0, first state occupied, transition must be forward
      A = A+1;
      dt = dt1(s);
    end
    t = t+dt; % Update time
    traj(s) = A; % Update state and time vectors
    time(s) = t;
    end
    DT(trial) = t; % Time to reach final state
    traj = traj(1:s); % Truncate time and trajectory vectors
    time = time(1:s);
end

if trials > 1
% Plot histogram of time to reach state N across all trials
figure(1)
[histy, histx] = hist(DT,0:1:max(DT)+1);
bar(histx,histy,'BarWidth', 1);
xlabel(strcat({'Time to reach state N = '},num2str(N)));
ylabel('Number of molecules');
end

% Plot trajectory of last simulated molecule
figure(2)
plot(time,traj);
title('Representative Individual Molecular Trajectory');
xlabel('Time');
ylabel('State');

% Store results in results struct for output
results.histx = histx;
results.histy = histy;
results.time = time;
results.traj = traj;
results.N = N;
results.trials = trials;
results.k1 = k1;
results.k2 = k2;

% histy_c = cumsum(histy);
end

function [N, trials, k1, k2] = parse_inputs(args)
  % Default arguments
  N = 50;
  trials = 100;
  k1 = 1;
  k2 = 0;

  % Parse arguments and provide defaults
  if length(args) > 3
    N = args{1,1}; %Number of steps in pathway (number of states = N + 1) i.e., (State 0) <--> (State 1) <--> (State 2) <-->...<--> (State N)
    trials = args{2}; % Number of molecules to simulate
    k1 = args{3}; % Rate constant of forward transition (0-->1)
    k2 = args{4}; % Rate constant of reverse transition (1-->0)
  elseif length(args) > 2
    k1 = args{3};
    trials = args{2};
    N = args{1,1};
  elseif length(args) > 1
    trials = args{2};
    N = args{1,1};
  elseif length(args) > 0
    N = args{1,1};
  else
    disp('No argument supplied; using defaults.');
  end
end

function display_parameters(N,trials,k1,k2)
  disp('-------------------------------------');
  disp('Parameters:');
  disp(strcat('N=',num2str(N)));
  disp(strcat('trials=',num2str(trials)));
  disp(strcat('k1=',num2str(k1)));
  disp(strcat('k2=',num2str(k2)));
  disp('-------------------------------------');
end
