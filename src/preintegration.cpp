#include "preintegration.h"


Preintegration::Preintegration()
{
    Reset();
}

void Preintegration::Reset()
{
    delta_t_ = 0.0;
    delta_p_ = Eigen::Vector3d::Zero();
    delta_v_ = Eigen::Vector3d::Zero();
    delta_q_ = Eigen::Quaterniond::Identity();
    bg_ = Eigen::Vector3d::Zero();
    ba_ = Eigen::Vector3d::Zero();
    cov_ = Eigen::Matrix<double, 9, 9>::Zero();
    last_imu_ = Imu();
    last_imu_flag_ = false;
    g_ << 0, 0, -GRAVITY;
}

void Preintegration::Integrate(const Imu& imu)
{
    /** state sequence: delta_q, delta_v, delta_p, bg, ba */
    double dt = 0.0;
    Eigen::Vector3d acc = Eigen::Vector3d::Zero();
    Eigen::Vector3d gyr = Eigen::Vector3d::Zero();

    if (!last_imu_flag_)
    {
        last_imu_ = imu;
        last_imu_flag_ = true;
        return;
    }

    dt = imu.timestamp - last_imu_.timestamp;
    acc = 0.5 * (imu.acc + last_imu_.acc) - ba_;
    gyr = 0.5 * (imu.gyr + last_imu_.gyr) - bg_;

    /** Update delta_p, delta_v, delta_q */
    delta_p_ = delta_p_ + delta_v_ * dt + 0.5 * delta_q_.toRotationMatrix() * acc * dt * dt;
    delta_v_ = delta_v_ + delta_q_.toRotationMatrix() * acc * dt;
    delta_q_ = delta_q_ * Eigen::Quaterniond(Eigen::AngleAxisd(dt * 0.5 * gyr.norm(), gyr.normalized()));
    delta_q_.normalize();

    /** Updata covariance matrix */
    Eigen::Matrix<double, 9, 9> A = Eigen::Matrix<double, 9, 9>::Zero();
    Eigen::Matrix<double, 9, 6> B = Eigen::Matrix<double, 9, 6>::Zero();
    Eigen::Matrix<double, 6, 6> Q = Eigen::Matrix<double, 6, 6>::Zero();

    A.block<3, 3>(0, 3) = Eigen::Matrix3d::Identity();
    A.block<3, 3>(3, 6) = -delta_q_.toRotationMatrix();
    A.block<3, 3>(6, 0) = -delta_q_.toRotationMatrix();
    A.block<3, 3>(6, 3) = Eigen::Matrix3d::Identity();

    B.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
    B.block<3, 3>(3, 3) = Eigen::Matrix3d::Identity();
    B.block<3, 3>(6, 4) = -delta_q_.toRotationMatrix();

    Q.block<3, 3>(0, 0) = acc_noise_.cwiseProduct(acc_noise_).asDiagonal();
    Q.block<3, 3>(3, 3) = gyr_noise_.cwiseProduct(gyr_noise_).asDiagonal();

    cov_ = A * cov_ * A.transpose() + B * Q * B.transpose();
}

void Preintegration::evaluate(const State& state0, const State& state1, 
                              Eigen::Matrix<double, 9, 1>& residuals, 
                              Eigen::Matrix<double, 9, 15>& jacobians)
{
    /** state sequence: delta_q, delta_v, delta_p, bg, ba */
    double dt = state1.t - state0.t;
    Eigen::Vector3d corrected_p_ = Eigen::Vector3d::Zero();
    Eigen::Vector3d corrected_v_ = Eigen::Vector3d::Zero();
    Eigen::Quaterniond corrected_q_ = Eigen::Quaterniond::Identity();
    
    /** 暂时不考虑偏置，即一个预积分窗口内的偏置都相同 */
    corrected_p_ = delta_p_;    
    corrected_v_ = delta_v_;
    corrected_q_ = delta_q_;

    // Residuals
    residuals.block<3, 1>(0, 0) = 2 * (corrected_q_.inverse() * state0.q.inverse() * state1.q).vec();
    residuals.block<3, 1>(3, 0)  = state0.q.inverse() * (state1.v - state0.v - g_ * delta_t_) - corrected_v_;
    residuals.block<3, 1>(6, 0)  = state0.q.inverse() * (state1.p - state0.p - state0.v * delta_t_ -
                                                       0.5 * g_ * delta_t_ * delta_t_) - corrected_p_;

    // Jacobians
    jacobians.block<3, 3>(0, 0) = -2 * (corrected_q_.inverse() * state0.q.inverse() * state1.q).toRotationMatrix();
    jacobians.block<3, 3>(0, 3) = -2 * (corrected_q_.inverse() * state0.q.inverse() * state1.q).toRotationMatrix() * state0.q.inverse().vec();
    jacobians.block<3, 3>(0, 6) = -2 * (corrected_q_.inverse() * state0.q.inverse() * state1.q).toRotationMatrix() * state1.q.inverse().vec();


}

