# QuantumWave

This C program implements a numerical solution to the time-dependent Schrödinger equation for a quantum particle in a potential field. Using the finite difference method, it calculates the spatial and temporal evolution of the particle's wavefunction. By modeling systems like finite potential wells, quantum barriers, or periodic lattices, the code can simulate the dynamics of quantum particles under realistic conditions. The output is the probability distribution of the particle's position over time, granting insight into its quantized energy states and wave-like behavior. Running simulations for different potentials allows investigation of quintessential quantum phenomena like tunneling and interference. The program provides an accessible tool for modeling fundamental quantum mechanics principles. It enables understanding these non-intuitive microscopic processes by transforming abstract equations into tangible wavefunction solutions. With versatile potential field configurations and visualization of the probability densities, this program opens up the nanoscale world of quantum particles to explore and better comprehend.

# how to run

C COMPILER ---> gcc -o QuantumWave QuantumWave.c

RUN ---> ./QuantumWave file timesteps

example --> ./QuantumWave test.dat 300 > answer.txt  |
-----------------------------------------------------|


----------------------------------------------
### tutorial https://youtu.be/1-FmUKs6yxY ###
----------------------------------------------

This program models the behavior of a single quantum particle, such as an electron, placed in a one-dimensional potential field. By numerically solving the Schrödinger equation, it allows calculating the wave function describing the quantum states of the particle.

The test presented examined the quantum state of an electron trapped in a finite potential well ranging from x=0 to x=0.5. The wave function exhibits a maximum probability of finding the particle at the center of the well, while its amplitude decreases with distance from it. This is the expected form of the wave function for the ground state of a particle in such a potential.

The test result confirms the correctness of the numerical implementation and allows verifying the model's consistency with the predictions of quantum mechanics for a simple, well-defined system. The program can be used to study more complex 1D systems, such as particles in periodic potentials, tunneling through barriers, or collisions. This will allow modeling and analyzing real quantum phenomena with applications in solid-state physics, among others.


# Simple explanation:
The correct answer is where the wave function reaches a value of 1 at a single point (x = 0.5), while being equal to 0 at all other points.

This stems from the fact that the wave function describes the probability of finding the quantum particle at a given point.

A value of 1 signifies certainty - a 100% probability of the particle occurring at that location.

Whereas 0 indicates the probability of finding the particle there is zero.

Such a form of the wave function, with a single global maximum, corresponds to a single quantum state of a particle trapped in the region around the point x = 0.5.

This is consistent with expectations for a simple model of a quantum potential well.

Any other wave function values (other than 1 or 0) would indicate an error in calculations or incorrect model assumption.

Therefore, the answer with a single maximum of 1 and surrounding zeros is correct and confirms the proper operation of the program.
