#ifndef _GINAV_PREINTEGRATION_H_
#define _GINAV_PREINTEGRATION_H_


#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "types.h"

class Preintegration
{
private:
    double delta_t_;
    Eigen::Vector3d delta_p_;
    Eigen::Vector3d delta_v_;
    Eigen::Quaterniond delta_q_;
    Eigen::Vector3d bg_;
    Eigen::Vector3d ba_;
    Eigen::Matrix<double, 9, 9> cov_;

    Eigen::Vector3d g_;
    Eigen::Vector3d gyr_noise_;
    Eigen::Vector3d acc_noise_;
    Imu last_imu_;
    bool last_imu_flag_;

public:
    inline double delta_t() const { return delta_t_; }
    inline Eigen::Vector3d delta_p() const { return delta_p_; }
    inline Eigen::Vector3d delta_v() const { return delta_v_; }
    inline Eigen::Quaterniond delta_q() const { return delta_q_; }
    inline Eigen::Vector3d bg() const { return bg_; }
    inline Eigen::Vector3d ba() const { return ba_; }
    inline Eigen::Matrix<double, 9, 9> cov() const { return cov_; }

    void Integrate(const Imu& imu);
    void Reset();

    void evaluate(const State& state0, const State& state1, 
                  Eigen::Matrix<double, 9, 1>& residuals, Eigen::Matrix<double, 9, 15>& jacobians);


public:
    Preintegration(/* args */);
    ~Preintegration();
};




#endif // _GINAV_PREINTEGRATION_H_