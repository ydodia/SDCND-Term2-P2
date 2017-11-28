#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
		const vector<VectorXd> &ground_truth)
{
	VectorXd rmse = VectorXd::Zero(4);
	VectorXd diff = VectorXd::Zero(4);
	unsigned int n = estimations.size();
	for (auto k = 0; k < n; ++k)
	{
		diff = estimations[k] - ground_truth[k];
		diff = diff.array().pow(2.0);
		rmse += VectorXd(diff);
	}
	rmse /= n;
	rmse = rmse.array().sqrt();
	return rmse;
}
