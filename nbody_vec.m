close all;
clear all;
clc;

%% Create Your Own N-body Simulation (With Matlab/Octave)
% Philip Mocz (2020) Princeton Univeristy, @PMocz

% Simulate orbits of stars interacting due to gravity
% Code calculates pairwise forces according to Newton's Law of Gravity

% uses VECTORIZATION

%% Simulation parameters
N         = 100;   % Number of particles
t         = 0;     % current time of the simulation
tEnd      = 10;    % time at which simulation ends
dt        = 0.01;  % timestep
softening = 0.1;   % softening length
G         = 1;     % Newton's Gravitational Constant


%% Generate Initial Conditions
rng(42);                % set the random number generator seed

mass = 20*ones(N,1)/N;  % total mass of particles is 20
pos = randn(N,3);       % randomly selected positions and velocities
vel = randn(N,3);

% calculate initial gravitational accelerations
acc = getAcc_vec( pos, mass, G, softening );

% number of timesteps
Nt = ceil(tEnd/dt);


%% Simulation Main Loop
tic
for i = 1:Nt
    
    % (1/2) kick
    vel = vel + acc * dt/2;
    
    % drift
    pos = pos + vel * dt;
    
    % update accelerations
    acc = getAcc_vec( pos, mass, G, softening );
    
    % (1/2) kick
    vel = vel + acc * dt/2;
    
    % update time
    t = t + dt;
 
end
toc
