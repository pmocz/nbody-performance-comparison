import numpy as np
import matplotlib.pyplot as plt
import time
from numba import jit

"""
Create Your Own N-body Simulation (With Python)
Philip Mocz (2020) Princeton Univeristy, @PMocz

Simulate orbits of stars interacting due to gravity
Code calculates pairwise forces according to Newton's Law of Gravity

uses FOR LOOPS with NUMBA

"""

@jit(nopython=True)
def getAcc( pos, mass, G, softening ):
	"""
    Calculate the acceleration on each particle due to Newton's Law 
	pos  is an N x 3 matrix of positions
	mass is an N x 1 vector of masses
	G is Newton's Gravitational constant
	softening is the softening length
	a is N x 3 matrix of accelerations
	"""
	
	N = pos.shape[0];
	a = np.zeros((N,3));
	
	for i in range(N):
		for j in range(N):
			dx = pos[j,0] - pos[i,0];
			dy = pos[j,1] - pos[i,1];
			dz = pos[j,2] - pos[i,2];
			inv_r3 = (dx**2 + dy**2 + dz**2 + softening**2)**(-1.5);
			a[i,0] +=  G * (dx * inv_r3) * mass[j];
			a[i,1] +=  G * (dy * inv_r3) * mass[j];
			a[i,2] +=  G * (dz * inv_r3) * mass[j];

	return a
	

def main():
	""" N-body simulation """
	
	# Simulation parameters
	N         = 100    # Number of particles
	t         = 0      # current time of the simulation
	tEnd      = 10.0   # time at which simulation ends
	dt        = 0.01   # timestep
	softening = 0.1    # softening length
	G         = 1.0    # Newton's Gravitational Constant
	
	# Generate Initial Conditions
	np.random.seed(17)            # set the random number generator seed
	
	mass = 20.0*np.ones(N)/N  # total mass of particles is 20
	pos  = np.random.randn(N,3)   # randomly selected positions and velocities
	vel  = np.random.randn(N,3)
	
	# calculate initial gravitational accelerations
	acc = getAcc( pos, mass, G, softening )
	
	# number of timesteps
	Nt = int(np.ceil(tEnd/dt))
	
	tStart = time.time()
	# Simulation Main Loop
	for i in range(Nt):
		# (1/2) kick
		vel += acc * dt/2.0
		
		# drift
		pos += vel * dt
		
		# update accelerations
		acc = getAcc( pos, mass, G, softening )
		
		# (1/2) kick
		vel += acc * dt/2.0
		
		# update time
		t += dt
	elapsed = time.time() - tStart
	print('Elapsed: %s' % elapsed)

	return 0
	


  
if __name__== "__main__":
  main()
