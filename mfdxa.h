#ifndef MFMXA_H
#define MFMXA_H

#include <cstddef>
#include <cstdlib>
#include <vector>
#include "armadillo"

using data_vector = std::vector<double>;
using data_matrix = std::vector<std::vector<double>>;

struct FluctuationData
{
	data_vector fluctuation_vector;
	data_vector sizes;
};

struct Coeffs
{
	arma::vec x_coeffs;
	arma::vec y_coeffs;
};

class MFDXA
{
public:
	std::vector<size_t> win_sizes;
	int poly_order;

	MFDXA(std::vector<size_t> win_sizes, int poly_order);

	~MFDXA();

	void get_profiles(data_vector& x, data_vector& y);
	void split_series(size_t w_size);
	Coeffs fit_poly(data_vector& x_frame, data_vector& y_frame, arma::mat& poly_matrix);
	FluctuationData calculate_local_coovariance_and_fluctuation(data_vector& x, data_vector& y, double q, bool verbose);

private:
	data_vector X;
	data_vector Y;
	data_matrix x_frames;
	data_matrix y_frames;

};


#endif
