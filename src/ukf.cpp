#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 25.0;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3.0 - n_aug_;
  is_initialized_ = false;
  weights_ = VectorXd(2*n_aug_+1);
  weights_.fill(0.5/(lambda_ + n_aug_));
  weights_(0) = lambda_/(lambda_ + n_aug_);
  Xsig_pred_ = MatrixXd(n_x_,2*n_aug_ + 1);
  Xsig_pred_.fill(0.0);
 
  lidarNisStream_.open("lidarNis.csv",fstream::out);
  radarNisStream_.open("radarNis.csv",fstream::out);

}

UKF::~UKF() {
lidarNisStream_.close();
radarNisStream_.close();
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if (!is_initialized_) {
    time_us_ = meas_package.timestamp_;
    double init_vel = 5.0;
    double init_psi = 0.0;
    double init_psi_dot = 5.0;

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

      x_<<meas_package.raw_measurements_[0],
          meas_package.raw_measurements_[1],
          init_vel,
          init_psi,
          init_psi_dot;
      P_.setIdentity();
      P_(0,0) = 0.15;
      P_(1,1) = 0.15;

    }
    else {

      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];

      x_<< rho*cos(phi),
           rho*sin(phi),
           init_vel,
           init_psi,
           init_psi_dot;
      P_.setIdentity();
    }

    is_initialized_ = true;
    cout<<"UKF initialized"<<endl;
  }
  double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);

  if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
    UpdateRadar(meas_package);
  }
  else {
    UpdateLidar(meas_package);
  }

  cout<<"x_ : \n"<<x_<<endl;
  cout<<"P_: \n"<<P_<<endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Complete this function! Estimate the object's location. Modify the state
  TODO:
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  cout<<"In prediction, delta_t: "<<delta_t<<endl;
  int n_aug = 7;
  VectorXd x_aug = VectorXd(n_aug);
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;
  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);
  P_aug.setZero();
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented sigma points
   Xsig_aug.setZero();
   Xsig_aug.col(0) = x_aug;
   for (int i =0; i< n_aug_; ++i)
   {
     Xsig_aug.col(i+1) = x_aug + sqrt(lambda_+n_aug_)*A.col(i);
     Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_)*A.col(i);
   }	
    //predict sigma points
    //avoid division by zero
    //write predicted sigma points into right column
    Xsig_pred_.setZero();
    int num_cols = Xsig_aug.cols();

    for (int i=0; i < num_cols; ++i) {
        VectorXd sig_point = Xsig_aug.col(i);
        VectorXd sig_pt_predict(5), sig_pt_delta1(5), sig_pt_delta2(5);
        sig_pt_delta2 << (0.5*delta_t*delta_t*cos(sig_point(3))*sig_point(5)),
                         (0.5*delta_t*delta_t*sin(sig_point(3))*sig_point(5)),
                         (delta_t*sig_point(5)),
                         (0.5*delta_t*delta_t*sig_point(6)),
                         (delta_t*sig_point(6));
        // if psi_dot is zero
        if(sig_point(n_x_-1) < 0.001) {
          sig_pt_delta1<<sig_point(2)*cos(sig_point(3))*delta_t,
                         sig_point(2)*sin(sig_point(3))*delta_t,
                         0.0,
                         0.0,
                         0.0;
        }
        else {
          double temp_ratio = sig_point(2)/sig_point(4);
          sig_pt_delta1<<  temp_ratio*(sin(sig_point(3)+(sig_point(4)*delta_t)) - sin(sig_point(3))),
                           temp_ratio*(-cos(sig_point(3)+(sig_point(4)*delta_t)) + cos(sig_point(3))),
                           0.0,
                           sig_point(4)*delta_t,
                           0.0;
        }
        sig_pt_predict = sig_point.head(n_x_)+sig_pt_delta1+sig_pt_delta2;
        Xsig_pred_.col(i) = sig_pt_predict;
    }
    

  //predict state mean
    VectorXd x = VectorXd(n_x_);
    x.fill(0.0);
    for (int i = 0; i < Xsig_pred_.cols();++i) {
        x+=(weights_(i)*Xsig_pred_.col(i));
    }
    cout<<"x predict: \n"<<x<<endl;
    x_ = x;

    //predict state covariance matrix
    MatrixXd P = MatrixXd(5,5);
    P.fill(0.0);
    for(int i =0; i< Xsig_pred_.cols(); ++i)
    {
        VectorXd diff_sig = Xsig_pred_.col(i)- x_;
        //angle normalization
        while (diff_sig(3)> M_PI) diff_sig(3)-=2.*M_PI;
        while (diff_sig(3)<-M_PI) diff_sig(3)+=2.*M_PI;

        P+=weights_(i)*diff_sig*(diff_sig.transpose());
    }
    P_ = P;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  cout<<"In Lidar update"<<endl;

  VectorXd z = meas_package.raw_measurements_;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(2, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(2);
  MatrixXd S = MatrixXd(2,2);
  MatrixXd Tc = MatrixXd(n_x_,2);

  Zsig.fill(0.0);
 for (int i = 0; i< Xsig_pred_.cols(); ++i)
 {
     VectorXd sp = Xsig_pred_.col(i);
     VectorXd sp_lidar_space(2);
     sp_lidar_space.fill(0.0);

     sp_lidar_space(0) = sp(0);
     sp_lidar_space(1) = sp(1);

     Zsig.col(i) = sp_lidar_space;
 }


 //calculate mean predicted measurement
 z_pred.fill(0.0);
 for (int i = 0; i< Zsig.cols(); ++i) {
   z_pred+=weights_(i)*Zsig.col(i);
 }
 //calculate innovation covariance matrix S
 S.fill(0.0);
 MatrixXd R(2,2);
 R.fill(0.0);
 R(0,0) = std_laspx_* std_laspx_;
 R(1,1) = std_laspy_* std_laspy_;

 for(int i = 0; i < Zsig.cols(); ++i) {
     VectorXd z_diff = Zsig.col(i) - z_pred;
     S+= (weights_(i)*z_diff*z_diff.transpose());
 }

 S+=R;

  //calculate cross correlation matrix
  Tc.fill(0.0);

  for(int i =0; i < Xsig_pred_.cols(); ++i ) {

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3) < M_PI) x_diff(3)+= 2.*M_PI;
    while (x_diff(3) > M_PI) x_diff(3)-= 2.*M_PI;

    VectorXd z_diff = Zsig.col(i) - z_pred;
    Tc += (weights_(i)*x_diff*z_diff.transpose());

  }
 //calculate Kalman gain K;
 MatrixXd K = Tc*S.inverse();

 //update state mean and covariance matrix
 x_+=(K*(z-z_pred));
 P_-=(K*S*K.transpose());
 double lidar_nis = (z-z_pred).transpose()*S.inverse()*(z-z_pred);
 cout<<"Lidar NIS: "<<lidar_nis<<endl;
 tools.writeValueToStream(lidar_nis,lidarNisStream_);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  cout<<"Update radar "<<endl;
  VectorXd z = meas_package.raw_measurements_;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(3, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(3);
  MatrixXd S = MatrixXd(3,3);
  MatrixXd Tc = MatrixXd(n_x_,3);

  Zsig.fill(0.0);
 for (int i = 0; i< Xsig_pred_.cols(); ++i)
 {
     VectorXd sp = Xsig_pred_.col(i);
     VectorXd sp_radar_space(3);
     sp_radar_space.fill(0.0);
     double px_2 = pow(sp(0),2), py_2= pow(sp(1),2);
     sp_radar_space(0) = sqrt(px_2 +py_2);
     sp_radar_space(1) = atan2(sp(1),sp(0));
     if(px_2 + py_2 > 0.01)
       sp_radar_space(2) = ((sp(0)*cos(sp(3))*sp(2)) + (sp(1)*sin(sp(3))*sp(2)))/ sp_radar_space(0);

     Zsig.col(i) = sp_radar_space;
 }

 //calculate mean predicted measurement
 z_pred.fill(0.0);
 for (int i = 0; i< Zsig.cols(); ++i) {
   z_pred+=weights_(i)*Zsig.col(i);
 }

 //calculate innovation covariance matrix S
 S.fill(0.0);
 MatrixXd R(3,3);
 R.fill(0.0);
 R(0,0) = std_radr_* std_radr_;
 R(1,1) = std_radphi_* std_radphi_;
 R(2,2) = std_radrd_* std_radrd_;

 for(int i = 0; i < Zsig.cols(); ++i) {
     VectorXd x_diff = Zsig.col(i) - z_pred;

     while(x_diff(1) < M_PI) x_diff(1)+= 2.*M_PI;
     while(x_diff(1) > M_PI) x_diff(1)-= 2.*M_PI;

     S+= (weights_(i)*x_diff*x_diff.transpose());
 }

 S+=R;

  //calculate cross correlation matrix
  Tc.fill(0.0);

  for(int i =0; i < Xsig_pred_.cols(); ++i ) {

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3) < M_PI) x_diff(3)+= 2.*M_PI;
    while (x_diff(3) > M_PI) x_diff(3)-= 2.*M_PI;

    VectorXd z_diff = Zsig.col(i) - z_pred;
    while (z_diff(1) < M_PI) z_diff(1)+= 2.*M_PI;
    while (z_diff(1) > M_PI) z_diff(1)-= 2.*M_PI;

    Tc += (weights_(i)*x_diff*z_diff.transpose());

  }

 //calculate Kalman gain K;
 MatrixXd K = Tc*S.inverse();

 //update state mean and covariance matrix
 x_+=(K*(z-z_pred));
 P_-=(K*S*K.transpose());
 double radar_nis = (z-z_pred).transpose()*S.inverse()*(z-z_pred);
 cout<<"Radar NIS: "<<radar_nis<<endl;
 tools.writeValueToStream(radar_nis,radarNisStream_);
}
