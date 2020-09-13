"""
Create Your Own N-body Simulation (With Julia)
Philip Mocz (2020) Princeton Univeristy, @PMocz

Simulate orbits of stars interacting due to gravity
Code calculates pairwise forces according to Newton's Law of Gravity

uses VECTORIZATION

"""

using Random
using Statistics

"""
Calculate the acceleration on each particle due to Newton's Law
pos  is an N x 3 matrix of positions
mass is an N x 1 vector of masses
G is Newton's Gravitational constant
softening is the softening length
a is N x 3 matrix of accelerations
"""
function getAcc( pos, mass, G, softening )

	# positions r = [x,y,z] for all particles
	x = pos[:,1];
	y = pos[:,2];
	z = pos[:,3];

	# matrix that stores all pairwise particle separations: r_j - r_i
	dx = x' .- x;
	dy = y' .- y;
	dz = z' .- z;

	# matrix that stores 1/r^3 for all particle pairwise particle separations
	inv_r3 = (dx.^2 .+ dy.^2 .+ dz.^2 .+ softening.^2).^(-1.5);

	ax = G .* (dx .* inv_r3) * mass;
	ay = G .* (dy .* inv_r3) * mass;
	az = G .* (dz .* inv_r3) * mass;

	# pack together the acceleration components
	a = [ax ay az];

	return a;
end



""" N-body simulation """
function main()
	# Simulation parameters
	N         = 100;    # Number of particles
	t         = 0;      # current time of the simulation
	tEnd      = 10;     # time at which simulation ends
	dt        = 0.01;   # timestep
	softening = 0.1;    # softening length
	G         = 1;      # Newton's Gravitational Constant

	# Generate Initial Conditions
	rng = MersenneTwister(42);    # set the random number generator seed

	mass = 20*ones(N,1)/N;        # total mass of particles is 20
	pos  = randn(rng,N,3);        # randomly selected positions and velocities
	vel  = randn(rng,N,3);

	# calculate initial gravitational accelerations
	acc = getAcc( pos, mass, G, softening );
	
	# number of timesteps
	Nt = convert(Int64,ceil(tEnd/dt));


	# Simulation Main Loop
	@time for i in 1:Nt
		# (1/2) kick
		vel .+= acc .* dt/2;

		# drift
		pos .+= vel .* dt;

		# update accelerations
		acc = getAcc( pos, mass, G, softening );

		# (1/2) kick
		vel .+= acc .* dt/2;

		# update time
		t += dt;
		
	end


end

main()
