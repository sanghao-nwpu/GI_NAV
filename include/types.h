#ifndef _GINAV_TYPES_H_
#define _GINAV_TYPES_H_

#include <Eigen/Dense>

// 1. 基本类型定义

const double DEG2RAD = 0.017453292519943295769236907684886;
const double RAD2DEG = 57.295779513082320876798154814105;

const double GRAVITY = 9.80665;

enum class Status { 
    OK = 0, 
    ERROR = 1, 
};

enum class MoveStatus { 
    UNKNOWN = 0, 
    MOVING = 1, 
    STOPED = 2, 
};

enum class NavStatus {
    UNINITIALIZED = 0, 
    LOW_PRECISION = 1, 
    HIGH_PRECISION = 2, 
};

struct VehicleStatus 
{
    /* data */
    double timestamp;
    MoveStatus status;
};


struct Imu { 
    double timestamp;
    Eigen::Vector3d gyr;
    Eigen::Vector3d acc;    
};

struct Gnss { 
    double timestamp;
    Eigen::Vector3d blh;
    Eigen::Vector3d dev;
    double vel_forward;
    double vel_track;   // 北偏东为正，deg
    double vel_upward;
    double hdop;
    int sat_num;
    int type;   // 0: None, 1: Single, 4: RTK fixed, 5: RTK float
};


struct State { 
    double t;
    Eigen::Vector3d p;
    Eigen::Quaterniond q;
    Eigen::Vector3d v;
    Eigen::Vector3d bg;
    Eigen::Vector3d ba;
    NavStatus status;
};

struct Result 
{
    double timestamp;
    Eigen::Vector3d pos;
    Eigen::Vector3d euler;
    Eigen::Vector3d vel_body;
    Eigen::Vector3d acc_body;
    Eigen::Vector3d gyr_bias;
    Eigen::Vector3d acc_bias;
    NavStatus status;
};



#endif // _GINAV_TYPES_H_