#include <vector>
#include <cmath>
#include <complex>
#include <iostream>
#include <omp.h>

struct DOAConfig {
    int n_target = 0;
    int n_samples = 512;
    double amplitude = 2;
    // double snr_db = 30;

    int antenna[8][3] = {{0, -9, 9}, {0, -1, 9}, {0, 9, 9},
                         {0, -5, 1}, {0, 5, -1},
                         {0, -9, -9}, {0, 1, -9}, {0, 9, -9}};

    int n_ant = sizeof(antenna) / sizeof(antenna[0]);

    long c = 3*1e10;     // cm
    long fc = 2.4*1e9;   // 2.4G
    
    double lambda = (double)c/fc; // meter
    double d = lambda/2;
    
    double r = 1;
    std::vector<double> doaAzimuth = {};
    std::vector<double> doaElevation = {};

};

std::vector<std::complex<float>> getSteering_v(const DOAConfig &cfg, const float &azimuth, const float &elevation) {
    std::vector<std::complex<float>> steering_v(8);

    for (int i = 0; i < 8; ++i) {

        float phaseShift_pi = (2 * cfg.d / cfg.lambda) * 
                                (cfg.antenna[i][0] * cos(azimuth) * sin(elevation) + 
                                cfg.antenna[i][1] * sin(azimuth) * sin(elevation) + 
                                cfg.antenna[i][2] * cos(elevation));

        steering_v[i] = cos(phaseShift_pi) + sin(phaseShift_pi);
    }

    return steering_v;
}

int main() {

    int num_procs = omp_get_num_procs();
    omp_set_num_threads(num_procs);
    std::cout << "Using " << num_procs << " threads.\n";

    int azi_offset = -60;
    int ele_offset = 60;

    #pragma omp parallel for collapse(2)
    for (int azi = 0; azi < 120; ++azi) {
        for (int ele = 0; ele < 36; ++ele) {
            float curr_azi = M_PI / 180.0 * (azi + azi_offset);
            float curr_ele = M_PI / 180.0 * (ele + ele_offset);
        }
    }

    return 0;
}