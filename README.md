# MPC for cart pendulum stabilization
 Model Predictive Control (MPC) for stabilization of a cart-pendulum system. The files in this repository evaluate LQR, linear MPC, nonlinear MPC and Reference Governor (RG) algorithms. More details about the implementation are given in the **Cart-Pendulum system stabilization via MPC control schemes** and **Cart-Pendulum system stabilization via MPC control schemes_Presentation** .pdf files.

## Setup for CasADi in Matlab
* Download CasADi binaries from [here](https://web.casadi.org/get/), unzip the files and place them in a directory of your choice. Then, in Matlab, run the following command depending on the file you downloaded.
```
addpath('<yourpath>/casadi-3.6.6)
```
In order to run this command each time Matlab starts, run
```
open startup.m
```
and paste the `addpath` command in the *startup.m* file.

## Main Files
* **CartPendulum_Linear_MPC.m** : Implements linear MPC and LQR, and compares the performance of both. The objective is to stabilize the pendulum at the upwards position from an initial angle close to said upwards position.
* **CartPendulum_LinearMPC_HorizontalPosition.m** : Implements linear MPC and LQR, and compares the performance of both. The objective is to stabilize the pendulum at the upwards position while the cart moves horizontally.
* **CartPendulum_ReferenceGovernor_HorizontalPosition.m** : Implements LQR and Reference Governor (RG) aided LQR, and compares the performance of both. The objective is to stabilize the pendulum at the upwards position while the cart moves horizontally.
* **CartPendulum_NonlinearMPC_SwingUp.m** : Implements Nonlinear MPC (NMPC). The objective is to swing up the pendulum from a downwards position.
* **CartPendulum_SwitchingMPC_SwingUp.m** : Implements Switching NMPC for extremely-constrained wind-up maneuvers.The objective is to swing up the pendulum from a downwards position under small input and position constraints.

## Support Files
* **cartpend.m** : Differential equations that describe cart-pendulum dynamics.
* **cartpendCas.m** : Version of **cartpend.m** compatible with CasADi programming scheme.
* **createCartPendulumAnimation.m** : Function that creates an animation from the X matrix obtained from **CartPendulum_ReferenceGovernor_HorizontalPosition.m**.
* **formQPMatrices.m** : Function that creates matrices to set MPC as a Quadratic Programming (QP) problem.
* **myQP.m** : Simple dual projected optimization algorithm for solving QP problems, used for implementing MPC at each step.
* **nmpc.m** : Implementation of Nonlinear MPC by applying Sequential Quadratic Programming (SQP) by solving QP subsystems and adding up the results.
* **switching_nmpc.m** : Implementation of Nonlinear MPC with similar implementation to **nmpc.m**. This scheme switches Q and Qf matrices depending on proximity to target bearing.

## Media
* **windup_Switching_MPC.avi**: Video showing switching MPC being used to wind up pendulum from downwards position with extreme input constraints.
