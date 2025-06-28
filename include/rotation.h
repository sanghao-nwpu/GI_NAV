#ifndef ROTATION_H
#define ROTATION_H

#include <Eigen/Geometry>

namespace Rotation {


inline Eigen::Matrix3d RotvecToMatrix(const Eigen::Vector3d &rotvec) {
    double angle = rotvec.norm();
    Eigen::Vector3d vec = rotvec.normalized();
    return Eigen::AngleAxisd(angle, vec).toRotationMatrix();
}

inline Eigen::Quaterniond MatrixToQuat(const Eigen::Matrix3d &matrix) {
    return Eigen::Quaterniond(matrix);
}

inline Eigen::Matrix3d QuatToMatrix(const Eigen::Quaterniond &quaternion) {
    return quaternion.toRotationMatrix();
}

// ZYX旋转顺序, 前右下的IMU, 输出RPY
inline Eigen::Vector3d MatrixToEuler(const Eigen::Matrix3d &dcm) {
    Eigen::Vector3d euler;

    euler[1] = atan(-dcm(2, 0) / sqrt(dcm(2, 1) * dcm(2, 1) + dcm(2, 2) * dcm(2, 2)));

    if (dcm(2, 0) <= -0.999) {
        euler[0] = atan2(dcm(2, 1), dcm(2, 2));
        euler[2] = atan2((dcm(1, 2) - dcm(0, 1)), (dcm(0, 2) + dcm(1, 1)));
    } else if (dcm(2, 0) >= 0.999) {
        euler[0] = atan2(dcm(2, 1), dcm(2, 2));
        euler[2] = M_PI + atan2((dcm(1, 2) + dcm(0, 1)), (dcm(0, 2) - dcm(1, 1)));
    } else {
        euler[0] = atan2(dcm(2, 1), dcm(2, 2));
        euler[2] = atan2(dcm(1, 0), dcm(0, 0));
    }

    // heading 0~2PI
    if (euler[2] < 0) {
        euler[2] = M_PI * 2 + euler[2];
    }

    return euler;
}

inline Eigen::Vector3d QuatToEuler(const Eigen::Quaterniond &quaternion) {
    return MatrixToEuler(quaternion.toRotationMatrix());
}

inline Eigen::Quaterniond RotvecToQuat(const Eigen::Vector3d &rotvec) {
    double angle = rotvec.norm();
    Eigen::Vector3d vec = rotvec.normalized();
    return Eigen::Quaterniond(Eigen::AngleAxisd(angle, vec));
}

inline Eigen::Vector3d QuatToRotvec(const Eigen::Quaterniond &quaternion) {
    Eigen::AngleAxisd axisd(quaternion);
    return axisd.angle() * axisd.axis();
}

// RPY --> C_b^n, ZYX顺序
inline Eigen::Matrix3d EulerToMatrix(const Eigen::Vector3d &euler) {
    return Eigen::Matrix3d(Eigen::AngleAxisd(euler[2], Eigen::Vector3d::UnitZ()) *
                    Eigen::AngleAxisd(euler[1], Eigen::Vector3d::UnitY()) *
                    Eigen::AngleAxisd(euler[0], Eigen::Vector3d::UnitX()));
}

inline Eigen::Quaterniond EulerToQuat(const Eigen::Vector3d &euler) {
    return Eigen::Quaterniond(Eigen::AngleAxisd(euler[2], Eigen::Vector3d::UnitZ()) *
                        Eigen::AngleAxisd(euler[1], Eigen::Vector3d::UnitY()) *
                        Eigen::AngleAxisd(euler[0], Eigen::Vector3d::UnitX()));
}

inline Eigen::Matrix3d SkewSymmetric(const Eigen::Vector3d &vector) {
    Eigen::Matrix3d mat;
    mat << 0, -vector(2), vector(1), vector(2), 0, -vector(0), -vector(1), vector(0), 0;
    return mat;
}

inline Eigen::Matrix4d QuaternionLeft(const Eigen::Quaterniond &q) {
    Eigen::Matrix4d ans;
    ans(0, 0)             = q.w();
    ans.block<1, 3>(0, 1) = -q.vec().transpose();
    ans.block<3, 1>(1, 0) = q.vec();
    ans.block<3, 3>(1, 1) = q.w() * Eigen::Matrix3d::Identity() + SkewSymmetric(q.vec());
    return ans;
}

inline Eigen::Matrix4d QuaternionRight(const Eigen::Quaterniond &p) {
    Eigen::Matrix4d ans;
    ans(0, 0)             = p.w();
    ans.block<1, 3>(0, 1) = -p.vec().transpose();
    ans.block<3, 1>(1, 0) = p.vec();
    ans.block<3, 3>(1, 1) = p.w() * Eigen::Matrix3d::Identity() - SkewSymmetric(p.vec());
    return ans;
}

}

#endif // ROTATION_H