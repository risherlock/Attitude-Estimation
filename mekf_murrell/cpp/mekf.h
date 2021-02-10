/* 
 Multiplicative Extended Kalman Filter with Murrell's algorithm. ~~

 References:
   [1] Markley, Crassidis - Fundamentals of Spacecraft Attitude
       Determination and Control (2014)
   [2] Crassidis, Junkins - Optimal Estimation of Dynamic Systems
        (2nd ed.) (2011)

 Note:
    1. The equations mentioned in the comments references to [2].
 
 Rishav (2020/2/5)
*/

#ifndef _MEKF_H_
#define _MEKF_H_

#include "lanczos/lanczos.h"

class MEKF
{
private:
    Matrix<6,6> Q;      // Process noise covariance
    Vector<6> dx_hat;   // Updated state estimate
    Matrix<6,6> P_hat;  // Updated estimate covariance
    float var_r[2];     // Vector sensors' noise variance [1xn]

    void computeProcessCovariance(const float sigma_g[3], const float sigma_b[3], float dt);
    void propQuaternion(const Quaternion &q, const Vector<3> &w, const float dt);
    void upQuaternion(const Quaternion &q, const Vector<3> &del_alpha);
    void propStateCovariance(const Matrix<6,6> &P, const Matrix<6,6> &Q, 
                            const Vector<3> &w, float dt);

public:
    Quaternion  q_hat;
    Vector<3>   b_hat;
    
    MEKF(const Quaternion q_init, const Vector<3> &b_init, const Matrix<6,6> &P_init, 
        const float sigma_r[2], const float sigma_g[3], const float sigma_b[3], float dt);

    void estimate(const Vector<3> w, const Vector<3> mb[2], const Vector<3> mr[2], float dt);
};

MEKF::MEKF(const Quaternion q_init, const Vector<3> &b_init, const Matrix<6,6> &P_init, 
            const float sigma_r[2], const float sigma_g[3], const float sigma_b[3], float dt)
{
    q_hat = q_init;
    b_hat = b_init;
    P_hat = P_init;
    computeProcessCovariance(sigma_g, sigma_b, dt); // Eq.(7.46)
    
    // Compute variances from standard deviations
    for (uint8_t i_iters = 0; i_iters < 2; i_iters++)
    {
        var_r[i_iters] = sigma_r[i_iters] * sigma_r[i_iters]; 
    }
}

/*
 Inputs:
    mb[n] = n sensor measurement vectors in body frame [3xn]
    mr[n] = n measurement vectors in inertial frame [3xn]
    dt    = Sample time (sec)

 Outputs:
    Updates quaternion and bias
*/
void MEKF::estimate
(const Vector<3> w, const Vector<3> mb[2], const Vector<3> mr[2], float dt)
{   
    // Propagate state, gyro-bias, quaternion and estimate covariance
    // Note: Gyro bias propagation is constant. Eq.(7.42b) (i.e. b = b)
    propQuaternion(q_hat, w, dt); // Eq.(7.39)
    propStateCovariance(P_hat, Q, w, dt); // Eq.(7.43)
    dx_hat = {0,0,0,0,0,0};
    Matrix<3,3> A = q_hat.getDCM(); // Attitude matrix

    // Murrell's loop
    for (uint8_t i_iters = 0; i_iters < 2; i_iters++)
    {   
        Matrix<3,3> I3({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
        Matrix<6,6> I6({{1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0},
                        {0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1}}); 
        
        // Sensitivity matrix (H[3x6]) and residal (e[3x1])
        Vector<3> Ar = A * mr[i_iters];
        Matrix<3,6> H({{0, -Ar(2), Ar(1)}, {Ar(2), 0, -Ar(0)}, {-Ar(1), Ar(0), 0}});
        Vector<3> e = mb[i_iters] - Ar;
        
        // MEKF update
        Matrix<6,3> K = P_hat * trans(H) * inv(H * P_hat * trans(H) + var_r[i_iters] * I3);
        dx_hat = dx_hat + K * (e - H * dx_hat); // Eq.(7.30)
        P_hat = (I6 - K * H) * P_hat;
    }
    // Update q_hat and b_hat
    Vector<3> dq({dx_hat(0), dx_hat(1), dx_hat(2)});
    Vector<3> db({dx_hat(3), dx_hat(4), dx_hat(5)});
    upQuaternion(q_hat, dq); // Eq.(7.34)
    b_hat = b_hat + db; // Eq.(7.32)
}

// Discrete Q using gyro and gyro-bias noise standard deviations
void MEKF::computeProcessCovariance
(const float sigma_g[3], const float sigma_b[3], float dt)
{
    // Eq.(7.46)
    float sq_gx = sigma_g[0] * sigma_g[0];
    float sq_gy = sigma_g[1] * sigma_g[1];
    float sq_gz = sigma_g[2] * sigma_g[2];
    float sq_bx = sigma_b[0] * sigma_b[0];
    float sq_by = sigma_b[1] * sigma_b[1];
    float sq_bz = sigma_b[2] * sigma_b[2];

    float a = dt * dt * dt / 3.0;
    float b = 0.5 * dt * dt;
 
    // First quadrant
    Q.set(0, 0, sq_gx * dt + sq_bx * a);
    Q.set(1, 1, sq_gy * dt + sq_by * a);
    Q.set(2, 2, sq_gz * dt + sq_bz * a);

    // Second quadrant
    Q.set(0, 3, sq_bx * b);
    Q.set(1, 4, sq_by * b);
    Q.set(2, 5, sq_bz * b);
    
    // Third quadrant
    Q.set(3, 0, Q(0,3));
    Q.set(4, 1, Q(1,4));
    Q.set(5, 2, Q(2,5));

    // Fourth quadrant
    Q.set(3, 3, sq_bx * dt);
    Q.set(4, 4, sq_by * dt);
    Q.set(5, 5, sq_bz * dt);

    Matrix<6,6> Y;
    Y.set(0,0,-1);  Y.set(1,1,-1);  Y.set(2,2,-1);
    Y.set(3,3,1);   Y.set(4,4,1);   Y.set(5,5,1);

    Q = Y*Q*trans(Y); // YQY' is used in Eq.(7.43)
}

// Discrete-time quaternion propagation
void MEKF::propQuaternion
(const Quaternion &q, const Vector<3> &w, const float dt)
{
    float omega_tol = 1e-5;
    float n = norm(w);
    if (n > omega_tol)
    {   
        // Eq.(7.40)
        Matrix<4,4> Omega;
        float c = cos(0.5 * n * dt);
        float s = sin(0.5 * n * dt) / n;
        float x = w(0) * s;
        float y = w(1) * s;
        float z = w(2) * s;
        Omega = {{c, z, -y, x}, {-z, c, x, y}, {y, -x, c, z}, {-x, -y, -z, c}};
        q_hat = Omega * q_hat; // Eq.(7.39)
    }
}

// Update quaternion using error-quaternion
void MEKF::upQuaternion
(const Quaternion &q, const Vector<3> &del_alpha)
{   
    Matrix<4,3> Xi;
    Xi = {{q(3), -q(2), q(1)}, 
        {q(2), q(3), -q(0)}, 
        {-q(1), q(0), q(3)}, 
        {-q(0), -q(1), -q(2)}};
    
    q_hat = q + 0.5*Xi*del_alpha; // Eq.(7.34)
    q_hat = q_hat/norm(q_hat); // Quaternion normalization
}

// Discrete propagation of the state error covariance
void MEKF::propStateCovariance
(const Matrix<6,6> &P, const Matrix<6,6> &YQY_, const Vector<3> &w, float dt)
{
    float omega_tol = 1e-5;
    float n = norm(w);
    Matrix<6,6> Phi;
    
    if (n > omega_tol)
    {
        // Eq.(7.45)
        float s = sin(n*dt);
        float c = cos(n*dt);
        float p = s/n;
        float q = (1 - c)/(n*n);
        float r = (n*dt - s)/(n*n*n);

        Matrix<3,3> wX({{0, -w(2), w(1)} , {w(2), 0, -w(0)} , {-w(1), w(0), 0}});
        Matrix<3,3> I({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
        Matrix<3,3> wXwX = wX*wX;
        
        Matrix<3,3> Phi_00 = I - wX*p + wXwX*q; 
        Matrix<3,3> Phi_01 = wX*q - I*dt - wXwX*r;
        
        // Phi_00                   // Phi_01 
        Phi.set(0,0,Phi_00(0,0));   Phi.set(0,3,Phi_01(0,0));   
        Phi.set(0,1,Phi_00(0,1));   Phi.set(0,4,Phi_01(0,1));   
        Phi.set(0,2,Phi_00(0,2));   Phi.set(0,5,Phi_01(0,2));   
        Phi.set(1,0,Phi_00(1,0));   Phi.set(1,3,Phi_01(1,0));
        Phi.set(1,1,Phi_00(1,1));   Phi.set(1,4,Phi_01(1,1));
        Phi.set(1,2,Phi_00(1,2));   Phi.set(1,5,Phi_01(1,2));
        Phi.set(2,0,Phi_00(2,0));   Phi.set(2,3,Phi_01(2,0));
        Phi.set(2,1,Phi_00(2,1));   Phi.set(2,4,Phi_01(2,1));
        Phi.set(2,2,Phi_00(2,2));   Phi.set(2,5,Phi_01(2,2));
    }
    else
    {
        // Steady-state Phi: Eq.(7.45) and Eq.(7.51b)
        // Phi_00           // Phi_01 
        Phi.set(0,0,1);     Phi.set(0,3,-dt);   
        Phi.set(1,1,1);     Phi.set(1,4,-dt);   
        Phi.set(2,2,1);     Phi.set(2,5,-dt);   
    }

    // Phi_11
    Phi.set(3,3,1);
    Phi.set(4,4,1);
    Phi.set(5,5,1);
    P_hat = Phi*P*trans(Phi) + YQY_;
}

#endif // mekf.h
