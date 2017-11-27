#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
{
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd::Zero(5);

  // initial covariance matrix
  P_ = MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.25;// 1.8; //mean is 0.231 m/s^2

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0;//0.7;

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

  //
  n_x_ = 5;
  n_aug_ = n_x_ + 2;
  n_col_aug_ = 2 * n_aug_ + 1;
  lambda_ = 3 - n_aug_;

  P_aug_ = MatrixXd::Zero(n_aug_, n_aug_);
  P_aug_.topLeftCorner(5,5) = P_;
  Q_ = MatrixXd::Zero(2,2);
  Q_ << std_a_*std_a_, 0, 0, std_yawdd_*std_yawdd_;
  P_aug_.bottomRightCorner(2,2) = Q_;
  //cerr << "P_aug:\n" << P_aug_ << endl;

  // weights
  weights_ = VectorXd(n_col_aug_);
  weights_(0) = lambda_ /(lambda_ + n_aug_);
  for (int i=1; i<n_col_aug_; ++i)
	  weights_(i) = 0.5 / (lambda_ + n_aug_);
  //cerr << "weights_:\n" << weights_ << endl;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
	static long previous_timestamp;
	const VectorXd& raw_meas = meas_package.raw_measurements_;
	if (!is_initialized_)
	{
		x_ << 1, 1, 0, 0, 0;

		if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
		{
			float rho = raw_meas(0);
			float phi = raw_meas(1);
			x_(0) = rho * cos(phi);
			x_(1) = rho * sin(phi);
		}
		else // LASER
		{
			x_(0) = raw_meas(0);
			x_(1) = raw_meas(1);
		}
		previous_timestamp = meas_package.timestamp_;
		is_initialized_ = true;
		return;
	}

	double dt =  (meas_package.timestamp_ - previous_timestamp) / 1000000.0;
	cerr << "dt = " << dt << endl;
	previous_timestamp = meas_package.timestamp_;
	Prediction(dt);
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
		UpdateRadar(meas_package);
	else
	if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
		UpdateLidar(meas_package);
	cerr << "Leaving Process\n";
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t)
{cerr << "Prediction()\n";
	double dt = delta_t;
	// P1 - Generate sigma points
	VectorXd x_aug(n_aug_);
	x_aug << x_, 0, 0;

	P_aug_.topLeftCorner(5,5) = P_;
	/*
	MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
	P_aug.topLeftCorner(5,5) = P_;
	MatrixXd Q(2,2);
	Q << std_a_*std_a_, 0, 0, std_yawdd_*std_yawdd_;
	P_aug.bottomRightCorner(2,2) = Q;
	cerr << "P_aug:\n" << P_aug << endl;
	*/

	//create square root matrix
	MatrixXd P_aug_sqrt = P_aug_.llt().matrixL();
	cerr << "P_aug_sqrt:\n" << P_aug_sqrt << endl;
	//create augmented sigma points
	MatrixXd P_pos = sqrt(lambda_ + n_aug_) * P_aug_sqrt;
	cerr << "P_pos:\n" << P_pos << endl;
	P_pos.colwise() += x_aug;
	MatrixXd P_neg = -1.0 * sqrt(lambda_ + n_aug_) * P_aug_sqrt;
	cerr << "P_neg:\n" << P_neg << endl;
	P_neg.colwise() += x_aug;
	MatrixXd Xsig_aug(n_aug_, n_col_aug_);
	Xsig_aug << x_aug, P_pos, P_neg;
	cerr << "Xsig_aug:\n" << Xsig_aug << endl;
	// P2 - Predict sigma points
	Xsig_pred_ = MatrixXd(n_x_, n_col_aug_);
	for (int k=0; k < n_col_aug_; ++k)
	{
		VectorXd c = Xsig_aug.col(k);
		VectorXd x = c.head(5);
		double px  = c[0];
		double py  = c[1];
		double v   = c[2];
		double psi = c[3];
		double psi_dot = c[4];
		double nu_a = c[5];
		double nu_psidd = c[6];
		VectorXd a(5), b(5);
		if (fabs(psi_dot) < 0.000001)
		{cerr << "psi_dot is 0\n";
			a << v*cos(psi)*dt, v*sin(psi)*dt, 0, psi_dot*dt, 0;
			b << 0.5*dt*dt*cos(psi)*nu_a, 0.5*dt*dt*sin(psi)*nu_a, dt*nu_a, 0.5*dt*dt*nu_psidd, dt*nu_psidd;
		}
		else
		{
			a << v*(sin(psi + psi_dot*dt) - sin(psi))/psi_dot,
					v*(-cos(psi + psi_dot*dt) + cos(psi))/psi_dot,
					0,
					psi_dot*dt,
					0;
			b << 0.5*dt*dt*cos(psi)*nu_a, 0.5*dt*dt*sin(psi)*nu_a, dt*nu_a, 0.5*dt*dt*nu_psidd,
					dt*nu_psidd;
		}
		Xsig_pred_.col(k) = x + a + b;
		//cerr << "Xsig_pred[:," << k << "]:\n" << Xsig_pred_.col(k) << endl;
	}

	// P3 - Predict mean and covariance
	// predict state mean
	x_ = weights_(0) * Xsig_pred_.col(0);
	for (int k=1; k < n_col_aug_; ++k)
	    x_ += weights_(k) * Xsig_pred_.col(k);
	cerr << "[x_]:\n" << x_ << endl;
	// predict state covariance matrix
	P_.fill(0.0);
	for (int k=0; k<n_col_aug_; ++k)
	{
		VectorXd a = Xsig_pred_.col(k) - x_;
		P_ += weights_(k) * (a * a.transpose());
	}
	cerr << "weights_:\n" << weights_ << endl;
	cerr << "[x_aug]:\n" << x_aug << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package)
{cerr << "UpdateLidar\n";
		// U1 - Predict Measurement
		//transform sigma points into measurement space
		static int n_z = 2; //px, py
		MatrixXd Zsig = MatrixXd(n_z, n_col_aug_);
		VectorXd z_pred = VectorXd::Zero(n_z);
		MatrixXd S = MatrixXd::Zero(n_z, n_z);
		cerr << "Xsig_pred_:\n" << Xsig_pred_ << endl;
		for (int k=0; k<n_col_aug_; ++k)
		{
			VectorXd x = Xsig_pred_.col(k);
			double px = x(0);
			double py = x(1);
			Zsig.col(k) << px, py;
		}
		cerr << "Zsig:\n" << Zsig << endl;
		//calculate mean predicted measurement
		for (int k=0; k<n_col_aug_; ++k)
			z_pred += weights_(k) * Zsig.col(k);
		//calculate measurement covariance matrix S
		for (int k=0; k<n_col_aug_; ++k)
		{
			VectorXd z_diff = Zsig.col(k) - z_pred;
			S += weights_(k) * (z_diff * z_diff.transpose());
		}
		MatrixXd R = MatrixXd::Zero(n_z, n_z);
		R(0,0) = std_laspx_ * std_laspx_;
		R(1,1) = std_laspy_ * std_laspy_;
		S += R;
		cerr << "-= 1 =-\n";
		// U2 - Update State
		// calculate cross correlation matrix
		VectorXd& z = meas_package.raw_measurements_;
		MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
		for (int k = 0; k < n_col_aug_; ++k)
		{
			VectorXd x_diff = Xsig_pred_.col(k) - x_;
			VectorXd z_diff = Zsig.col(k) - z_pred;
			Tc += weights_(k) * (x_diff * z_diff.transpose());
		}
		cerr << "-= 2 =-\n";
		//calculate Kalman gain K;
		MatrixXd K = Tc * S.inverse();
		cerr << "-= 3 =-\n";
		//update state mean and covariance matrix
		cerr << "z:\n" << z << endl;
		cerr << "z_pred:\n" << z_pred << endl;
		VectorXd z_err = z - z_pred;
		cerr << "-= 4 =-\n";
		x_ += K * z_err;
		P_ += -1.0 * (K * S * K.transpose());

		// NIS
		double e = z_err.transpose() * S.inverse() * z_err;
		cerr << "LASER NIS - " << e << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package)
{cerr << "UpdateRadar\n";
	// U1 - Predict Measurement
	//transform sigma points into measurement space
	static int n_z = 3;
	MatrixXd Zsig = MatrixXd(n_z, n_col_aug_);
	VectorXd z_pred = VectorXd::Zero(n_z);
	MatrixXd S = MatrixXd::Zero(n_z, n_z);
	cerr << "Xsig_pred_:\n" << Xsig_pred_ << endl;
	for (int k=0; k<n_col_aug_; ++k)
	{
		VectorXd x = Xsig_pred_.col(k);
		double px = x(0);
		double py = x(1);
		double v  = x(2);
		double psi = x(3);
		double rho = sqrt(px*px + py*py);
		if (rho < 0.00001) rho = 0.00001;
		double phi = atan2(py, px);
		double rho_dot = v * (px*cos(psi) + py*sin(psi)) / rho;
		Zsig.col(k) << rho, phi, rho_dot;
	}
	cerr << "Zsig:\n" << Zsig << endl;
	//calculate mean predicted measurement
	for (int k=0; k<n_col_aug_; ++k)
		z_pred += weights_(k) * Zsig.col(k);
	//calculate measurement covariance matrix S
	for (int k=0; k<n_col_aug_; ++k)
	{
		VectorXd z_diff = Zsig.col(k) - z_pred;
		S += weights_(k) * (z_diff * z_diff.transpose());
	}
	MatrixXd R = MatrixXd::Zero(n_z, n_z);
	R(0,0) = std_radr_ * std_radr_;
	R(1,1) = std_radphi_ * std_radphi_;
	R(2,2) = std_radrd_ * std_radrd_;
	S += R;

	// U2 - Update State
	// calculate cross correlation matrix
	VectorXd& z = meas_package.raw_measurements_;
	MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
	for (int k = 0; k < n_col_aug_; ++k)
	{
		VectorXd x_diff = Xsig_pred_.col(k) - x_;
		VectorXd z_diff = Zsig.col(k) - z_pred;
		Tc += weights_(k) * (x_diff * z_diff.transpose());
	}
	//calculate Kalman gain K;
	MatrixXd K = Tc * S.inverse();
	//update state mean and covariance matrix
	VectorXd z_err = z - z_pred;
	x_ += K * z_err;
	P_ += -1.0 * (K * S * K.transpose());

	// NIS
	double e = z_err.transpose() * S.inverse() * z_err;
	cerr << "RADAR NIS - " << e << endl;
}
