#ifndef GINAV_H
#define GINAV_H

#include <deque>
#include <Eigen/Dense>

#include "types.h"

const int STATE_DIM = 15;

class Ginav {
public:
    
    Ginav();
    ~Ginav();

    Status ProcessGnss(const Gnss& gnss);
    Status ProcessImu(const Imu& imu);
    Result GetResult();

private:
    bool imu_flag_;
    bool gnss_flag_;
    Imu imu_;
    Gnss gnss_;
    double g_;
    Eigen::Vector3d gnss_ori_;
    State state_;
    VehicleStatus vehicle_status_;

    std::deque<Imu> imu_array_;
    std::deque<Imu> imu_stat_array_;
    std::deque<Gnss> gnss_array_;

    Status Initialize();
    Status CalculateVehicleStatus();
    Status Predict();
    Status NhcUpdate();
    Status GnssUpdate();
    Status Correct(const Eigen::MatrixXd& H, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& N);
};
#endif // GINAV_H