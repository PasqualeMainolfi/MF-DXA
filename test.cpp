#include "mfdxa.h"
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>
#include <random>

const std::vector<size_t> WIN_SIZES = { 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8820 };

const std::vector<double> Q = { -2.0, -1.0, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0 };
const int POLI_ORDER = 2;
const size_t N = 16384;


int main() {

	std::vector<double> x_serie;
	std::vector<double> y_serie;

	std::random_device rd;
	std::mt19937 gen;
	std::uniform_real_distribution<double> dis(-10.0, 10.0);

	for (size_t i = 0; i < N; ++i) {
		double sine = std::sin(2.0 * M_PI * (double) i / 44100.0);
		double r1 = dis(gen) / 2.0;
		double r2 = dis(gen) / 3.0;
		x_serie.push_back(sine + r1);
		y_serie.push_back(sine + r2);
	}


	MFDXA mfdxa(WIN_SIZES, POLI_ORDER);

	double q = Q[5];
	FluctuationData fluctuation = mfdxa.calculate_local_coovariance_and_fluctuation(x_serie, y_serie, q, false);

	std::cout << "FLUCTUATION VECTOR" << std::endl;
	for (size_t i = 0; i < fluctuation.fluctuation_vector.size(); ++i) {
		std::cout << fluctuation.fluctuation_vector[i] << " ";
	}
	std::cout << std::endl;

	data_vector f = fluctuation.fluctuation_vector;
	data_vector s = fluctuation.sizes;

	std::ofstream data_file("temp_data.txt");
	for (size_t i = 0; i < f.size(); ++i) {
		data_file << s[i] << " " << f[i] << std::endl;
	}
	data_file.close();

	std::ofstream data_plot("plot_data.gp");
	data_plot << "set title 'MF-DXA'\n";
	data_plot << "set xlabel 'log(s)'\n";
	data_plot << "set ylabel 'log(Fs(s))'\n";
	data_plot << "plot 'temp_data.txt' with linespoints lc rgb 'black'\n";
	data_plot << "pause -1\n";
	data_plot.close();

	return 0;
}
