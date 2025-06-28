#ifndef _GINAV_OB_GINS_H_
#define _GINAV_OB_GINS_H_

#include <ceres/ceres.h>
#include <Eigen/Dense>
#include <preintegration.h>
#include "types.h"

class ImuFactor : public ceres::CostFunction {
public:
    explicit ImuFactor(std::shared_ptr<Preintegration> preintegration)
        : preint_(std::move(preintegration)) {

    }

    bool Evaluate(const double *const *parameters, double *residuals, double **jacobians) const override {

        // parameters: vel[3], bg[3], ba[3]

        preint_->imuErrorEvaluate(parameters, residuals);

        if (jacobians) {
            if (jacobians[0]) {
                preintegration_->imuErrorJacobian(jacobians[0]);
            }
        }

        return true;
    }

private:
    std::shared_ptr<Preintegration> preint_;
};

class GnssFactor {

};

class NhcFactor {

};


class ObGins {
public:
    ObGins() {}


};

#endif // _GINAV_OB_GINS_H_