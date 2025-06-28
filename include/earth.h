#ifndef EARTH_H
#define EARTH_H

#include <Eigen/Dense>

#define WGS84_WIE (7.2921151467E-5)            /** 地球自转角速度*/
#define WGS84_F (0.0033528106647474805)        /** 地球扁率 */
#define WGS84_RA (6378137.0000000000)          /** 长半轴a */
#define WGS84_RB (6356752.3142451793)          /** 短半轴b */
#define WGS84_GM0 (398600441800000.00)         /** 地球引力常数 */
#define WGS84_SQURE_E1 (0.0066943799901413156) /** 第一偏心率平方 */
#define WGS84_SQURE_E2 (0.0067394967422764341) /** 第二偏心率平方 */

namespace Earth {

/**
 * @brief 计算地球上某点的重力加速度
 * @param blh 纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 * @return 重力加速度(m/s^2)
 */
inline double CalculateGravity(const Eigen::Vector3d& blh)
{
    double squre_sin_lat = 0.0;
    double gravity = 0.0;
    squre_sin_lat = sin(blh.x());
    squre_sin_lat *= squre_sin_lat;
    /** g = 9.7803267715 * (1 + 0.0052790414 * sin^2(lat) + 0.0000232718 * sin^2(lat)) +
     * blh.z * (0.0000000043977311 * sin^2(lat) - 0.0000030876910891)
     * + 0.0000000000007211 * z^2 */
    gravity = 1 + 0.0052790414 * squre_sin_lat;
    gravity += 0.0000232718 * squre_sin_lat * squre_sin_lat;
    gravity = 9.7803267715 * gravity;
    gravity += blh.z() * (0.0000000043977311 * squre_sin_lat - 0.0000030876910891);
    gravity += 0.0000000000007211 * blh.z() * blh.z();
    return gravity;
}

/**
 * @brief 计算子午圈半径
 * @param lat 纬度坐标(latitude{rad})
 * @return 子午圈半径(m)
 */
inline double CalculateMeridianRadius(const double lat)
{
    double meridian_radius = 0.0;
    double squre_sin_lat = 0.0;
    squre_sin_lat = sin(lat);
    squre_sin_lat *= squre_sin_lat;
    meridian_radius = WGS84_RA * (1.0 - WGS84_SQURE_E1) / pow(1 - WGS84_SQURE_E1 * squre_sin_lat, 1.5);
    return meridian_radius;
}

/**
 * @brief 计算卯酉圈半径
 * @param lat 纬度坐标(latitude{rad})
 * @return 卯酉圈半径(m)
 */
inline double CalculatePrimeVerticalRadius(const double lat)
{
    double prime_vertical_radius = 0.0;
    double sinlat = sin(lat);
    prime_vertical_radius = WGS84_RA / sqrt(1.0 - WGS84_SQURE_E1 * sinlat * sinlat);
    return prime_vertical_radius;
}

/**
 * @brief 计算东北天系(enu)到地心地固系(ecef)的转换矩阵
 * @param blh 纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 * @return 转换矩阵
 */
inline Eigen::Matrix3d CalculateEnuToEcefMatrix(const Eigen::Vector3d& blh)
{
    Eigen::Matrix3d enu_to_ecef_matrix = Eigen::Matrix3d::Zero();
    double sinlat = sin(blh.x());
    double coslat = cos(blh.x());
    double sinlon = sin(blh.y());
    double coslon = cos(blh.y());

    enu_to_ecef_matrix(0, 0) = -sinlon;
    enu_to_ecef_matrix(0, 1) = -coslon * sinlat;
    enu_to_ecef_matrix(0, 2) = coslat * coslon;

    enu_to_ecef_matrix(1, 0) = coslon;
    enu_to_ecef_matrix(1, 1) = -sinlon * sinlat;
    enu_to_ecef_matrix(1, 2) = coslat * sinlon;

    enu_to_ecef_matrix(2, 0) = 0.0;
    enu_to_ecef_matrix(2, 1) = coslat;
    enu_to_ecef_matrix(2, 2) = sinlat;

    return enu_to_ecef_matrix;
}

/**
 * @brief 地心地固系(ecef)坐标转换为纬经高坐标(blh)
 * @param ecef 地心地固系坐标(x{m}, y{m}, z{m})
 * @return 纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 */
inline Eigen::Vector3d EcefToBlh(const Eigen::Vector3d& ecef)
{
    Eigen::Vector3d blh = Eigen::Vector3d::Zero();
    double p = sqrt(ecef.x() * ecef.x() + ecef.y() * ecef.y());
    double rn;
    double lat, lon;
    double h = 0.0, h2 = 0.0;

    // 初始状态
    lat = atan2(ecef.z(), (p * (1.0 - WGS84_SQURE_E1)));
    lon = atan2(ecef.y(), ecef.x());

    while (true)
    {
        h2 = h;
        rn = CalculatePrimeVerticalRadius(lat);
        h = p / cos(lat) - rn;
        // lat = atan(ecef.z() / (p * (1.0 - WGS84_SQURE_E1 * rn / (rn + h))));
        lat = atan2(ecef.z() + rn * WGS84_SQURE_E1 * sin(lat), p);
        if (fabs(h - h2) <= 1.0e-4)
        {
            break; // 当条件满足时退出循环
        }
    }

    blh.x() = lat;
    blh.y() = lon;
    blh.z() = h;

    return blh;
}

/**
 * @brief 纬经高坐标(blh)转换为地心地固系(ecef)坐标
 * @param blh 纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 * @return 地心地固系坐标(x{m}, y{m}, z{m})
 */
inline Eigen::Vector3d BlhToEcef(const Eigen::Vector3d& blh)
{
    Eigen::Vector3d ecef = Eigen::Vector3d::Zero();
    double coslat, sinlat, coslon, sinlon;
    double rnh, rn;

    coslat = cos(blh.x());
    sinlat = sin(blh.x());
    coslon = cos(blh.y());
    sinlon = sin(blh.y());

    rn = CalculatePrimeVerticalRadius(blh.x());
    rnh = rn + blh.z();
    ecef.x() = rnh * coslat * coslon;
    ecef.y() = rnh * coslat * sinlon;
    ecef.z() = (rnh - rn * WGS84_SQURE_E1) * sinlat;
    return ecef;
}

/**
 * @brief 纬经高坐标(blh)转换为东北天系(enu)坐标
 * @param ref_blh 参考点的纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 * @param blh 待转换点的纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 * @return 东北天系坐标(east{m}, north{m}, up{m})
 */
inline Eigen::Vector3d BlhToEnu(const Eigen::Vector3d& ref_blh, const Eigen::Vector3d& blh)
{
    Eigen::Vector3d enu = Eigen::Vector3d::Zero();
    Eigen::Vector3d ref_ecef = Eigen::Vector3d::Zero();
    Eigen::Vector3d ecef = Eigen::Vector3d::Zero();
    Eigen::Matrix3d enu_to_ecef_matrix = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d ecef_to_enu_matrix = Eigen::Matrix3d::Zero();
    Eigen::Vector3d delta_ecef = Eigen::Vector3d::Zero();

    ref_ecef = BlhToEcef(ref_blh);
    ecef = BlhToEcef(blh);
    enu_to_ecef_matrix = CalculateEnuToEcefMatrix(ref_blh);
    ecef_to_enu_matrix = enu_to_ecef_matrix.transpose();

    delta_ecef = ecef - ref_ecef;
    enu = ecef_to_enu_matrix * delta_ecef;

    return enu;
}

/**
 * @brief 东北天系(enu)坐标转换为纬经高坐标(blh)
 * @param ref_blh 参考点的纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 * @param enu 待转换点的东北天系坐标(east{m}, north{m}, up{m})
 * @return 纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 */
inline Eigen::Vector3d EnuToBlh(const Eigen::Vector3d& ref_blh, const Eigen::Vector3d& enu)
{
    Eigen::Vector3d blh = Eigen::Vector3d::Zero();
    Eigen::Vector3d ref_ecef = Eigen::Vector3d::Zero();
    Eigen::Vector3d ecef = Eigen::Vector3d::Zero();
    Eigen::Matrix3d enu_to_ecef_matrix = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d ecef_to_enu_matrix = Eigen::Matrix3d::Zero();
    Eigen::Vector3d delta_ecef = Eigen::Vector3d::Zero();

    ref_ecef = BlhToEcef(ref_blh);
    enu_to_ecef_matrix = CalculateEnuToEcefMatrix(ref_blh);

    delta_ecef = enu_to_ecef_matrix * enu;
    ecef = ref_ecef + delta_ecef;

    blh = EcefToBlh(ecef);

    return blh;
}

}

#endif // EARTH_H