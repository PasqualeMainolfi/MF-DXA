#include "armadillo"
#include "mfdxa.h"
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <vector>

data_vector get_profile(data_vector& data)
{
	double average = 0.0;
	for (const auto& x : data) {
		average += x;
	}
	size_t n = data.size();
	average /= n;

	data_vector profile;
	double cumulative_sum = 0.0;
	for (const auto& value : data) {
		cumulative_sum += value;
		profile.push_back(cumulative_sum - average);
	}

	return profile;
}

data_matrix split_data(data_vector& data, size_t winsize)
{
	data_matrix m;
	size_t n = data.size();
	size_t n_frames = std::ceil((float) n / (float) winsize);

	size_t hop = 0;
	for (size_t i = 0; i < n_frames; ++i) {
		std::vector<double> frame(winsize, 0.0);
		for (size_t j = 0; j < winsize; ++j) {
			if (j + hop >= n ) break;
			frame[j] = data[j + hop];
		}
		m.push_back(frame);
		hop += winsize; // NO OVERLAP!
	}

	return m;
}

arma::vec fit_polynomials(data_vector& x, arma::mat& poly_matrix)
{
	arma::vec frame(x);
	arma::vec coeffs = arma::solve(poly_matrix, frame);
	return coeffs;
}


double calculate_local_coovariance(data_vector& x, arma::vec& coeffs, arma::mat& poly_matrix)
{
	arma::vec frame(x);
	arma::vec trend = poly_matrix * coeffs;
	arma::vec detrend = frame - trend;
	double value = arma::sum(arma::square(detrend));
	return value;
}

double calculate_fluctuation(data_vector& local_coovariance, double q)
{
	arma::vec c(local_coovariance);
	arma::vec c_sign = arma::sign(c);
	arma::vec c_abs = arma::pow(arma::abs(c), q / 2.0);
	double fluct = arma::sum(c_sign % c_abs) / (double) local_coovariance.size();
	fluct = std::pow(fluct, q);
	return fluct;
}

arma::mat get_poly_matrix(arma::vec& time, size_t frame_length, int order)
{
	arma::mat A(frame_length, order + 1);
	for (int i = 0; i <= order; ++i) {
		A.col(i) = arma::pow(time, i);
	}
	return A;
}

// MF-DXA CLASS

MFDXA::MFDXA(std::vector<size_t> win_sizes, int poly_order)
: win_sizes(win_sizes), poly_order(poly_order)
{ }

MFDXA::~MFDXA()
{ }

void MFDXA::get_profiles(data_vector& x, data_vector& y)
{
	if (x.size() != y.size()) {
		std::cout << "[ERROR] x and y must have same length!" << std::endl;
		std::exit(1);
	}

	this->X = get_profile(x);
	this->Y = get_profile(y);

}

void MFDXA::split_series(size_t w_size)
{
	this->x_frames = split_data(this->X, w_size);
	this->y_frames = split_data(this->Y, w_size);
}

Coeffs MFDXA::fit_poly(data_vector& x_frame, data_vector& y_frame, arma::mat& poly_matrix)
{
	Coeffs coeffs;
	coeffs.x_coeffs = fit_polynomials(x_frame, poly_matrix);
	coeffs.y_coeffs = fit_polynomials(y_frame, poly_matrix);
	return coeffs;
}

FluctuationData MFDXA::calculate_local_coovariance_and_fluctuation(data_vector& x, data_vector& y, double q, bool verbose)
{
	this->get_profiles(x, y);
	data_vector fluct;
	data_vector sizes;
	data_vector q_vector;

	for (const auto& frame_size : this->win_sizes) {
		this->split_series(frame_size);
		size_t n_frames = this->x_frames.size();
		arma::vec t = arma::regspace(0, frame_size - 1);
		arma::mat pm = get_poly_matrix(t, frame_size, this->poly_order);
		data_vector local_coovariance;
		for (size_t i = 0; i < n_frames; i++) {
			data_vector x_frame = this->x_frames[i];
			data_vector y_frame = this->y_frames[i];

			Coeffs coeffs = this->fit_poly(x_frame, y_frame, pm);
			double x_coov = calculate_local_coovariance(x_frame, coeffs.x_coeffs, pm);
			double y_coov = calculate_local_coovariance(y_frame, coeffs.y_coeffs, pm);
			double lc = (x_coov * y_coov) / (double) frame_size;
			local_coovariance.push_back(lc);
		}
		double f = calculate_fluctuation(local_coovariance, q);
		fluct.push_back(std::log(f));
		sizes.push_back(std::log(frame_size));
		if (verbose) {
			std::cout << "[INFO] Frame size: [" << frame_size << "] q: [" << q << "]" << std::endl;
			std::cout << "Fluctuation: [" << f << "]" << std::endl;
		}
	}

	FluctuationData data = {
		.fluctuation_vector = fluct,
		.sizes = sizes,
	};

	return data;
}
