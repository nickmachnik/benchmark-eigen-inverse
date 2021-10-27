#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

auto rnd_mat_inv(size_t n) -> double
{
    using std::chrono::duration;
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;

    Eigen::MatrixXd m = Eigen::MatrixXd::Random(n, n);

    auto t1 = high_resolution_clock::now();
    // Eigen::MatrixXd sol = (m.transpose() * m).llt().solve(Eigen::MatrixXd::Identity(n, n));
    Eigen::MatrixXd sol = (m.transpose() * m).inverse();
    auto t2 = high_resolution_clock::now();

    duration<double, std::milli> ms_double = t2 - t1;

    return ms_double.count();
}

auto main(int argc, char *argv[]) -> int
{
    Eigen::setNbThreads(4);

    std::vector<size_t> sizes{10, 50, 100, 500, 1000, 5000};
    size_t sample_size{10};

    for (size_t i = 0; i < sizes.size(); ++i) {
        double time_sum{0.0};
        std::vector<double> times(100);

        for (size_t j = 0; j < sample_size; ++j) {
            times[j] = rnd_mat_inv(sizes[i]);
            time_sum += times[j];
        }

        double time_mean = time_sum / sample_size;

        double sum_squares{0.0};
        for (auto &&i : times) {
            double dev = i - time_mean;
            sum_squares += dev * dev;
        }
        double std_dev = std::sqrt(sum_squares / sample_size);

        std::cout << sizes[i] << "," << time_mean << "," << std_dev << std::endl;
    }

    return 0;
}