# Unscented Kalman Filter - SDCND Term 2 P2

In this project, I used an *Unscented Kalman Filter (UKF)* to track a vehicle moving in a figure-eight pattern. The difference between this model and the EKF is that the UKF is better at non-linear tracking. The non-linear model used here is the *Continuous Turn Rate and Velocity (CTRV)* model which helps us handle the non-linearity of turns better than the EKF from P1.The example two datasets start the vehicle in opposite directions.

To compile and run the C++ code, follow the steps below:

1. cd build
2. cmake ..
3. make
4. ./ExtendedKF

Then, run the simulator and my code will track the simluated vehicles path.

## Results
My implementation results in RMSE values as follows for Dataset 1:

`
X  : 0.0663
Y  : 0.0835
Vx : 0.3367
Vy : 0.2229 .
`
As reference, this is an improvement upon the EKF from P1:

`
X  : 0.0973
Y  : 0.0854
Vx : 0.4512
Vy : 0.4396 .
`