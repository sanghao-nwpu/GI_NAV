
#include "rotation.h"
#include "earth.h"
#include "ginav.h"

Status Ginav::Predict() {
    double dt = 0.0;
    double timestamp = 0.0;
    Eigen::Vector3d gyr, acc, acc_n;
    Eigen::Matrix3d R, Rvi;
    Eigen::Vector3d gyr_noise, acc_noise;
    Eigen::Vector3d v, p;
    Eigen::Quaterniond q;
    Eigen::Vector3d delta;
    Eigen::Matrix<double, STATE_DIM, STATE_DIM> Phi = Eigen::Matrix<double, STATE_DIM, STATE_DIM>::Zero();
    Eigen::Matrix<double, STATE_DIM, 12> G = Eigen::Matrix<double, STATE_DIM, 12>::Zero();
    Eigen::Matrix<double, 12, 12> Q = Eigen::Matrix<double, 12, 12>::Zero();

    if (last_imu_flag_)
    {
        dt = imu_array_.back().timestamp - imu_array_[imu_array_.size() - 2].timestamp;
    }

    R = state_.quat.toRotationMatrix();
    v = state_.vel;
    p = state_.pos;
    q = state_.quat;


    timestamp = imu_array_.back().timestamp;
    Rvi = config_.Rvi;
    gyr = Rvi * (imu_array_.back().gyr - config_.gyr_bias);
    acc = Rvi * (imu_array_.back().acc - config_.acc_bias);

    // 状态预测
    state_.timestamp = timestamp;
    acc_n = state_.quat.toRotationMatrix() * (acc - state_.ab) + g_;
    state_.quat *= Rotation::RotvecToQuat((gyr - state_.gb) * dt);
    state_.vel = v + acc_n * dt;
    state_.pos = p + v * dt + 0.5 * acc_n * dt * dt;

    // 设置状态转移矩阵
    Phi.setZero();
    Phi.block<3, 3>(0, 0) = Rotation::RotvecToMatrix(-(gyr - state_.gb) * dt); /* Phi_rr */
    Phi.block<3, 3>(3, 0) = R * Rotation::SkewSymmetric(acc - state_.ab) * dt; /* Phi_vr */
    Phi.block<3, 3>(3, 3) = Eigen::Matrix3d::Identity();
    Phi.block<3, 3>(6, 3) = Eigen::Matrix3d::Identity() * dt;
    Phi.block<3, 3>(6, 6) = Eigen::Matrix3d::Identity();

    if (config_.estimate_gyr_bias)
    {
        Phi.block<3, 3>(0, 9) = -Eigen::Matrix3d::Identity() * dt;
        Phi.block<3, 3>(9, 9) = Eigen::Matrix3d::Identity();
    }
    if (config_.estimate_acc_bias)
    {
        Phi.block<3, 3>(3, 12) = -R * dt;
        Phi.block<3, 3>(12, 12) = Eigen::Matrix3d::Identity();
    }

    // 设置噪声驱动矩阵
    G.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity() * dt;  /* G_r_gyr */
    G.block<3, 3>(3, 3) = -R * dt;   /* G_v_acc */
    if (config_.estimate_gyr_bias)
    {
        G.block<3, 3>(9, 6) = Eigen::Matrix<double, 3, 3>::Identity() * dt;   
    }
    if (config_.estimate_acc_bias)
    {
        G.block<3, 3>(12, 9) = Eigen::Matrix<double, 3, 3>::Identity() * dt;
    }
    
    // 设置噪声矩阵
    gyr_noise = Rvi * config_.gyr_noise;
    acc_noise = Rvi * config_.acc_noise;
    Q(0, 0) = gyr_noise[0] * gyr_noise[0];
    Q(1, 1) = gyr_noise[1] * gyr_noise[1];
    Q(2, 2) = gyr_noise[2] * gyr_noise[2];
    Q(3, 3) = acc_noise[0] * acc_noise[0];
    Q(4, 4) = acc_noise[1] * acc_noise[1];
    Q(5, 5) = acc_noise[2] * acc_noise[2];

    if (config_.estimate_gyr_bias)
    {
        // TODO: 
        Q(6, 6) = 1e-6 * 1e-6;
        Q(7, 7) = 1e-6 * 1e-6;
        Q(8, 8) = 1e-6 * 1e-6;
    }

    if (config_.estimate_acc_bias)
    {
        // TODO：
        Q(9, 9) = 1e-6 * 1e-6;
        Q(10, 10) = 1e-6 * 1e-6;
        Q(11, 11) = 1e-6 * 1e-6;
    }

    // covariance_ = Phi * (covariance_ + G * Q * (G.transpose())) * (Phi.transpose());
    state_.cov = Phi * state_.cov * Phi.transpose() + G * Q * G.transpose();
    
    return Status::OK;
}


Status Ginav::Correct(const Eigen::MatrixXd& H, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& N)
{
    Eigen::MatrixXd P = state_.cov;
    Eigen::MatrixXd PHT = P * H.transpose();
    Eigen::MatrixXd S = H * PHT + N;

    Eigen::MatrixXd K = PHT * S.inverse();

    //计算更新量
    Eigen::VectorXd delta = K * Z;
    Eigen::Vector3d delta_rot = delta.segment(0, 3);

    //状态更新
    state_.quat = state_.quat * Rotation::RotvecToQuat(delta_rot);
    state_.vel = state_.vel + delta.segment(3, 3);
    state_.pos = state_.pos + delta.segment(6, 3);

    state_.gb = state_.gb + delta.segment(9, 3);
    state_.ab = state_.ab + delta.segment(12, 3);
    
    //协方差更新
    Eigen::MatrixXd IKH = Eigen::MatrixXd::Identity(STATE_DIM, STATE_DIM) - K * H;
    Eigen::MatrixXd P_new = IKH * P * (IKH.transpose()) + K * N * (K.transpose()); 
    state_.cov = 0.5 * (P_new + P_new.transpose());

    return Status::OK;
}


Status Ginav::GnssUpdate()
{
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(5, STATE_DIM);
    Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(5, 1);
    Eigen::MatrixXd N = Eigen::MatrixXd::Zero(5, 5);
    Eigen::Matrix3d Hv = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d Hp = Eigen::Matrix3d::Identity();
    Eigen::Vector3d vel_dev = Eigen::Vector3d::Zero();
    Eigen::Vector3d pos_dev = Eigen::Vector3d::Zero();
    Eigen::Vector3d pos = Eigen::Vector3d::Zero();
    Eigen::Vector3d vel = Eigen::Vector3d::Zero();

    pos = Earth::BlhToEnu(gnss_ori_, gnss_.blh);
    pos_dev = gnss_.dev;
    vel(0) = gnss_.vel_forward * sin(gnss_.vel_track * DEG2RAD);
    vel(1) = gnss_.vel_forward * cos(gnss_.vel_track * DEG2RAD);
    vel(2) = gnss_.vel_upward;
    vel_dev(0) = 0.1;
    vel_dev(1) = 0.1;
    vel_dev(2) = 0.1;

    N(0, 0) = vel_dev(0) * vel_dev(0);
    N(1, 1) = vel_dev(1) * vel_dev(1);
    N(2, 2) = pos_dev(0) * pos_dev(0);
    N(3, 3) = pos_dev(1) * pos_dev(1);
    N(4, 4) = pos_dev(2) * pos_dev(2);
    
    Z(0, 0) = vel(0) - state_.vel(0);
    Z(1, 0) = vel(1) - state_.vel(1); 
    Z(2, 0) = pos(0) - state_.pos(0);
    Z(3, 0) = pos(1) - state_.pos(1);
    Z(4, 0) = pos(2) - state_.pos(2);

    H.block<2, 3>(0, 3) = Hv.block<2, 3>(0, 0);

    H.block<3, 3>(2, 6) = Hp;

    Correct(H, Z, N);

    return Status::OK;
}


Status Ginav::ProcessImu(const Imu& imu)
{
    Status status = Status::OK;
    imu_array_.push_back(imu);
    if (200 < imu_array_.size())
    {
        imu_array_.pop_front();
    }

    status = CalculateVehicleStatus();
    if (Status::OK != status) return status;
    
    // 未初始化则直接返回
    if (NavStatus::UNINITED == state_.status)
    {
        return Status::STATUS_OK;
    }

    //IMU预测
    status = Predict();
    if (Status::STATUS_OK != status) return status;
    
    //Nhc更新
    // status = NhcUpdate();
    // if (Status::STATUS_OK != status) return status;

    // 更新上一帧IMU有效标识
    last_imu_flag_ = true;
 
    return Status::STATUS_OK;
}


Status Ginav::ProcessGnss(const Gnss& gnss)
{
    Status status = Status::STATUS_OK;
    
    gnss_array_.push_back(gnss);
    if (5 < gnss_array_.size())
    {
        gnss_array_.pop_front();
    }

    // 陀螺偏置未初始化
    if (!stat_info_.flag)
    {
        return Status::STATUS_OK;
    }

    gnss_info_.timestamp = gnss.timestamp;
    gnss_info_.sol_type = gnss.sol_type;

    // GNSS 无效，直接返回
    if (0 == gnss.sol_type)
    {
        return Status::STATUS_OK;
    }

    gnss_info_.pos = Earth::BlhToEnu(config_.gnss_ori, gnss.wgs);
    gnss_info_.pos_dev = gnss.pos_dev;

    // 测试：固定噪声
    gnss_info_.pos_dev(0) = 0.05;
    gnss_info_.pos_dev(1) = 0.05;
    gnss_info_.pos_dev(2) = 0.05;
    
    gnss_info_.vel_forward = gnss.vel_forward;
    gnss_info_.vel_track = gnss.vel_track * DEG2RAD;
    gnss_info_.vel_upward = gnss.vel_upward;
    gnss_info_.vel(0) = gnss.vel_forward * sin(gnss_info_.vel_track);
    gnss_info_.vel(1) = gnss.vel_forward * cos(gnss_info_.vel_track);
    gnss_info_.vel(2) = gnss.vel_upward;

    gnss_info_.vel_dev = {0.1, 0.1, 0.1};
    
    if (NavStatus::UNINITED == state_.status)
    {
        if (gnss.vel_forward > 1.0)
        {
            std::cout << "vel forward: " << gnss_info_.vel_forward << ";\n";
            std::cout << "vel track(rad): " << gnss_info_.vel_track << ";\n";
            std::cout << gnss_info_.vel << std::endl; 
            status = Initilize();
        }
        return status;
    }

    status = GnssUpdate();
    
    return status;
}


Status Ginav::Initilize()
{
    Eigen::Vector3d euler = Eigen::Vector3d::Zero();
    euler(2) = 0.5 * M_PI - gnss_info_.vel_track;
    state_.pos = gnss_info_.pos;
    state_.vel = gnss_info_.vel;
    state_.dcm = Rotation::EulerToDcm(euler);
    state_.gb = config_.gyr_bias;
    state_.ab = Eigen::Vector3d::Zero();
    state_.status = NavStatus::LOW_PRECISION;
    
    covariance_.setZero();
    covariance_.block<3, 3>(0, 0) = pow(1.0e+0 * DEG2RAD, 2) * Eigen::Matrix3d::Identity();  // Attitude
    covariance_.block<3, 3>(3, 3) = gnss_info_.vel_dev.cwiseProduct(gnss_info_.vel_dev).asDiagonal(); // velocity
    covariance_.block<3, 3>(6, 6) = gnss_info_.pos_dev.cwiseProduct(gnss_info_.pos_dev).asDiagonal(); // Position
    if (config_.gyrbias_flag)
    {
        covariance_.block<3, 3>(9, 9) = pow(3.0e-3 * DEG2RAD, 2) * Eigen::Matrix3d::Identity(); // gyr bias
    }
    if (config_.accbias_flag)
    {
        covariance_.block<3, 3>(12, 12) = pow(1.0e-1, 2) * Eigen::Matrix3d::Identity(); // Acc bias cov
    }
   
    std::cout << "state initilized! gnss time: " << gnss_info_.timestamp << "; pos: " << state_.pos(0) << ", " << state_.pos(1) << std::endl;

    return Status::STATUS_OK;
}


void Ginav::GetResult(Result& result)
{
    result.timestamp = state_.timestamp;
    result.pos = state_.pos;
    result.vel_enu = state_.vel;
    // result.vel_body = state_.quat.toRotationMatrix().transpose() * state_.vel;
    result.vel_body = state_.dcm.transpose() * state_.vel;
    
    result.euler = Rotation::DcmToEuler(state_.dcm) * RAD2DEG;
    result.gb = state_.gb;
    result.ab = state_.ab;
    result.status = state_.status;
}