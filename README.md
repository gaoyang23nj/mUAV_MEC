# mUAV_MEC
This project/repository provides the experiments about the joint optimization scheme of UAV trajectory, transmit power and GU-MEC association with multiple UAV and multiple GUs.
This joint optimization scheme tries to minimize the maximum AoT among all the GUs in the multi-UAV multi-GU MEC scenario,
where the UAV flying constraints, the co-channel inference and the energy constraints exist at the same time.
It also includes some other schemes (e.g., Peak Power, Static UAV, Local Comp.).

The whole project is implemented with Matlab.

## Implementation of the Algorithm (AoT OPT Alg.)
The main algorithm (i.e., AoT OPT Alg.) is implemented in [SolveP1bcdsca3subs.m](https://github.com/gaoyang23nj/mUAV_MEC/blob/main/SolveP1bcdsca3subs.m).

And this implementation is based on the general convex optimization problem solver (here we adopt CVX + Mosek).
This SolveP1bcdsca3subs.m can also adopt the other solvers, e.g., SDPT3 of CVX.

## Comparative Experiments
[VaryingKallAlgs.m](https://github.com/gaoyang23nj/mUAV_MEC/blob/main/VaryingK/VaryingKallAlgs.m) ealuates the performance of various schemes with varying K (i.e., the number of GUs).

[VaryingLkallAlgs.m](https://github.com/gaoyang23nj/mUAV_MEC/blob/main/VaryingLk/VaryingLkallAlgs.m) ealuates the performance of various schemes with varying L_k (i.e., the size of task-input data).

[CaseStudy](https://github.com/gaoyang23nj/mUAV_MEC/tree/main/CaseStudy) provides a case study about the optimization variables change with iteration.

