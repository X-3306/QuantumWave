# QuantumWave

that program uses the fourth-order Runge-Kutta method to solve the Schrödinger equation for a given potential. The program can take either a fixed potential or a file of potential values as input. The output is the wave function at each step, which can be used to calculate expectation values such as energy and momentum. The program uses arrays to store the real and imaginary parts of the wave function, computes the derivatives of the wave function at each step, and calculates the next step using the Runge-Kutta method. The program outputs the final wave function after a specified number of timesteps.

# how to run

C COMPILER - gcc -o QuantumWave QuantumWave.c

RUN - ./QuantumWave file timesteps

file is the name of the file that contains the potential data and timesteps is the number of timesteps for the simulation.
