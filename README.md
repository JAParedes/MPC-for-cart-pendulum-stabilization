# MPC for Inverted Pendulum on a Cart Stabilization
 Model Predictive Control (MPC) for stabilization of a cart-pendulum system. The files in this repository evaluate LQR, linear MPC, nonlinear MPC and Reference Governor (RG) algorithms. More details about the implementation are given in the **Inverted Pendulum on a Cart system stabilization via MPC control schemes** and **Inverted Pendulum on a Cart system stabilization via MPC control schemes_Presentation** .pdf files.

## Setup for CasADi in Matlab for Windows
* Download CasADi binaries from [here](https://web.casadi.org/get/), unzip the files and place them in a directory of your choice. Then, in Matlab, run the following command depending on the file you downloaded.
```
addpath('<yourpath>/casadi-3.6.6)
```
In order to run this command each time Matlab starts, run
```
open startup.m
```
and paste the `addpath` command in the *startup.m* file.

## Demonstration Video

The Youtube video in [this link](https://youtu.be/m_wtxx0UWF0) shows the results of each of the main files.

## Main Files
* **CartPendulum_Linear_MPC.m** : Implements linear MPC and LQR, and compares the performance of both. The objective is to stabilize the pendulum at the upwards position from an initial angle close to said upwards position.
* **CartPendulum_LinearMPC_HorizontalPosition.m** : Implements linear MPC and LQR, and compares the performance of both. The objective is to stabilize the pendulum at the upwards position while the cart moves horizontally.
* **CartPendulum_ReferenceGovernor_HorizontalPosition.m** : Implements LQR and Reference Governor (RG) aided LQR, and compares the performance of both. The objective is to stabilize the pendulum at the upwards position while the cart moves horizontally.
* **CartPendulum_NonlinearMPC_SwingUp.m** : Implements Nonlinear MPC (NMPC). The objective is to swing up the pendulum from a downwards position.
* **CartPendulum_SwitchingMPC_SwingUp.m** : Implements Switching NMPC for extremely-constrained wind-up maneuvers.The objective is to swing up the pendulum from a downwards position under small input and position constraints.

## Support Files
* **cartpend.m** : Differential equations that describe cart-pendulum dynamics.
* **cartpendCas.m** : Version of **cartpend.m** compatible with CasADi programming scheme.
* **Create_Inverted_Pendulum_on_a_Cart_animation.m** : Function that creates an animation from the .mat files created at the end of each of the main programs. Uncomment the instructions at the end of each program to create these files.
* **formQPMatrices.m** : Function that creates matrices to set MPC as a Quadratic Programming (QP) problem.
* **myQP.m** : Simple dual projected optimization algorithm for solving QP problems, used for implementing MPC at each step.
* **nmpc.m** : Implementation of Nonlinear MPC by applying Sequential Quadratic Programming (SQP) by solving QP subsystems and adding up the results.
* **switching_nmpc.m** : Implementation of Nonlinear MPC with similar implementation to **nmpc.m**. This scheme switches Q and Qf matrices depending on proximity to target bearing.
