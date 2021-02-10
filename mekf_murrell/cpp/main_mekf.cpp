/*
 Numerical validation with MATLAB implementation.
 (Check the file ../matlab/main_numerical_val)
 
 References:
   [1] Markley, Crassidis - Fundamentals of Spacecraft Attitude Determination and 
        Control (2014)
   [2] Crassidis, Junkins - Optimal Estimation of Dynamic Systems (2nd ed.) (2011)

 Rishav (2021/2/5)
*/

#include "mekf.h"
#include <iostream>

int main()
{
    Quaternion q_init({0,0,0,1}); // Initial quaternion
    Vector<3> b_init({1.324,2.4532,3.4532}); // Initial bias
    Matrix<6,6> P_init; // Initial state noise covariance
    P_init = {{1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0},
            {0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1}};
    float dt = 0.134; // Time step

    // Noise standard deviation
    float sigma_r[2] = {1.345, 2.654}; // Vector measurements noise
    float sigma_g[3] = {4.453, 2.534, 3.234}; // Gyro noise 
    float sigma_b[3] = {2.435, 1.2345, 4.55}; // Gyro-bias noise
    
    // Vector measurements
    Vector<3> mb[2];
    mb[0] = {1.32,2.65,3.63}; // Eg. Magnetometer
    mb[1] = {4.6534,5.634,6.634}; // Eg. Accelerometer
    
    // Reference vectors 
    Vector<3> mr[2]; 
    mr[0] = {7.53,8.23,9.97}; // Magnetometer vector in inertial frame
    mr[1] = {1.42,2.432,3.546}; // Acc. due to gravity in inertial frame

    // Gyro reading 
    Vector<3> w({0.223,2.34,3.45});
    
    MEKF mekf(q_init, b_init, P_init, sigma_r, sigma_g, sigma_b, dt); // Initialize MEKF
    mekf.estimate(w - mekf.b_hat, mb, mr, dt); // Run MEKF

    // Display
    std::cout << "Estimated quaternion:" << mekf.q_hat; 
    std::cout << std::endl;
    std::cout << "Estimated bias:" << mekf.b_hat;
    
    return 0;
}