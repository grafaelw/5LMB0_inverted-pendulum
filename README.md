# 5LMB0_inverted-pendulum
Assigments from 5LMB0 - Model Predictive Control course in TU EIndhoven for 2021/2022 academic year. The engineering problem is to implement a robust MPC, from linear MPC up to LPV-MPC to control the non-linear dynamics. The Electrically-driven Inverted Pendulum is based on Section 8.5.2 page 171 from Hanema J., Anticipative MPC for LPV, 2018. 

1. Graded Assignment 1

This assignment is to linearise the plant so that one can apply a linear MPC to control the pendulum's trajectory.
2. Graded Assignment 2

To obtain stability from MPC on the first assignment, one can use the Lyapunov's stability so that the stage cost can be maximised for stability (in [Part 1]:./Graded Project 2/part1.m). Other way around, one can also use LMI (Linear Matrix Inequalities) method to actually compute all of the components based on the dARE equations in [part 2]: ./Graded Project 2/part2.m 
3. Graded Assignment 3

This is the part where one can control the nonlinear dynamics using the LPV formulations. Here, one can run the first the [initialDataGenerator] 5LMB0_inverted-pendulum/Graded Project 3/initialDataGenerator.m to obtain the **Uk-en-Rho.mat** file. Afterwards, runing the [Simulink Simulation]: 5LMB0_inverted-pendulum/Graded Project 3/assignment3_simulink.m or the [Closed-Loop Simulation]: 5LMB0_inverted-pendulum/Graded Project 3/assignment3_ClosedLoop.m for Simulink and closed-loop simulations respectively
