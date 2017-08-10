#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  ekf_.P_ = MatrixXd(4,4);
  ekf_.P_ <<  1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1000, 0,
              0, 0, 0, 1000;

  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ <<  1, 0, 1, 0,
              0, 1, 0, 1,
              0, 0, 1, 0,
              0, 0, 0, 1;

  const double uncertain = 1000;
  const double certain = 0.01;
  VectorXd x_in(4);
  MatrixXd P_in(4, 4);
  P_in << certain, 0, 0, 0,
      0, certain, 0, 0,
      0, 0, uncertain, 0,
      0, 0, 0, uncertain;
  MatrixXd F_in(4, 4);
  MatrixXd Q_in(4, 4);

  ekf_.Init(x_in, P_in, F_in, H_laser_, R_laser_, Q_in);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  cout << "FusionEKF::ProcessMeasurement " << endl;


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    double px, py, vx, vy;
    // RADAR case
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      cout << "first measurement RADAR: " << endl;
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      double rho_dot = measurement_pack.raw_measurements_[2];

      px = rho * std::cos(phi);
      py = rho * std::sin(phi);
      vx = rho_dot * std::cos(phi);
      vy = rho_dot * std::sin(phi);
    }
    // LIDAR sensor
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      cout << "first measurement LASER: " << endl;
      /**
      Initialize state.
      */
      px = measurement_pack.raw_measurements_[0];
      py = measurement_pack.raw_measurements_[1];
      vx = 0;
      vy = 0;
    }

    if (fabs(px) < 1e-4){
      px = 1e-4;
    }
    if (fabs(py) < 1e-4){
      py = 1e-4;
    }

    ekf_.x_ << px, py, vx, vy;
    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    cout << "first measurement: " << endl;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

   cout << "Prediction: " << endl;
  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
   double noise_ax = 9.;
   double noise_ay = 9.;
   double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
   ekf_.F_ << 1, 0, dt, 0,
              0, 1, 0, dt,
              0, 0, 1, 0,
              0, 0, 0, 1;
  // update process noise matrix
  double dt_4 = std::pow(dt, 4) / 4.0;
  double dt_3 = std::pow(dt, 3) / 2.0;
  double dt_2 = std::pow(dt, 2);

  cout << "Update process noise matrix " << endl;
  ekf_.Q_ << dt_4 * noise_ax, 0, dt_3 * noise_ax, 0,
              0, dt_4 * noise_ay, 0, dt_3 * noise_ay,
              dt_3 * noise_ax, 0, dt_2 * noise_ax, 0,
              0, dt_3 * noise_ay, 0, dt_2 * noise_ay;
  cout << "B4 ekf_.Predict();" << endl;
  ekf_.Predict();
  previous_timestamp_ = measurement_pack.timestamp_;
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  cout << "Update: " << endl;
  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;

    auto square = [](double value) { return std::pow(value, 2); };
    VectorXd H_x(3);
    double px = ekf_.x_(0);
    double py = ekf_.x_(1);
    double vx = ekf_.x_(2);
    double vy = ekf_.x_(3);

    if (fabs(px) < 1e-4 || std::sqrt(square(px) + square(py)) < 1e-4) {
      // No update when zero division
      return;
    }

    H_x << std::sqrt(square(px) + square(py)),
        std::atan2(py, px),
        (px * vx + py * vy) / std::sqrt(square(px) + square(py));

    VectorXd y = measurement_pack.raw_measurements_ - H_x;

    // because range is between -pi and pi
    while (y(1) > M_PI) {
      y(1) -= M_PI;
    }

    while (y(1) < -M_PI) {
      y(1) += M_PI;
    }

    ekf_.Update(y);
    // ekf_.UpdateRADAR(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_ - ekf_.H_ * ekf_.x_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
