#include <vector>
#include <string>
#include <cstdint>
#include <cmath>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py=pybind11;

/* Compute the maximum value of two integers. */
int max(int a, int b) {
    if (a > b) {
        return a;
    }
    return b;
}

/* Compute the minimum value of two integers. */
int min(int a, int b) {
    if (a < b) {
        return a;
    }
    return b;
}

/* Compute the correlation coefficient, ignoring NaNs. */
float correlation_coefficient(const std::vector<float> x, const std::vector<float> y) {
    float sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, square_sum_x = 0.0, square_sum_y = 0.0;
    size_t i, n = 0;
    size_t N = x.size();
    for (i = 0; i < N; i++) {
        if (!std::isnan(x[i]) && !std::isnan(y[i])) {
            sum_x += x[i];
            sum_y += y[i];
            sum_xy += x[i] * y[i];
            square_sum_x += x[i] * x[i];
            square_sum_y += y[i] * y[i];
            n += 1;
        }
    }
    if (n > 0) {
        return (n * sum_xy - sum_x*sum_y)/sqrt((n*square_sum_x - sum_x*sum_x) * (n*square_sum_y - sum_y*sum_y));
    } else {
        return 0.0;
    }
}
// float correlation_coefficient(const std::vector<float> x, const std::vector<float> y) {
//     float sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, square_sum_x = 0.0, square_sum_y = 0.0;
//     size_t i, n = 0;
//     size_t N = x.size();
//     for (i = 0; i < N; i++) {
//         if (!isnan(x[i]) && !isnan(y[i])) {
//             sum_x += x[i];
//             sum_y += y[i];
//             sum_xy += x[i] * y[i];
//             square_sum_x += x[i] * x[i];
//             square_sum_y += y[i] * y[i];
//             n += 1;
//         }
//     }
//     if (n > 0) {
//         return (n * sum_xy - sum_x*sum_y)/sqrt((n*square_sum_x - sum_x*sum_x) * (n*square_sum_y - sum_y*sum_y));
//     } else {
//         return 0.0;
//     }
// }

/* Compute the cross-correlogram between two spike trains. Calculations and corrections are done using Equation 5.10
from "Analysis of Parallel Spike Trains" Sonja Grun & Stefan Rotter [eds] 2010 Volume 7, ISBN 978-1-4419-5674-3 */
std::vector<float> crosscorrelogram(const std::vector<int> spikes_1, const std::vector<int> spikes_2, float bin_width, float max_lag, std::string normalization) {
    float T, mean_val, denom;
    size_t i;
    size_t N = int(2 * max_lag / bin_width + 1);
    std::vector<float> result(N);
    size_t num_1 = spikes_1.size();
    size_t num_2 = spikes_2.size();
    size_t j;
    size_t j0 = 0;
    for (i = 0; i < num_1; i++) {
        j = j0;
        while ((j < num_2) && (float(spikes_2[j] - spikes_1[i]) < -max_lag - bin_width/2.0)) {
            j += 1;
        }
        j0 = j;
        while ((j < num_2) && (float(spikes_2[j] - spikes_1[i]) < max_lag + bin_width/2.0)) {
            result[int((float(spikes_2[j] - spikes_1[i]) + max_lag + 0.5*bin_width)/bin_width)] += 1.0;
            j += 1;
        }
    }
    if (normalization == "correlation") {
        int min_time = min(spikes_1[0], spikes_2[0]);
        int max_time = max(spikes_1[num_1-1], spikes_2[num_2-1]);
        T = float(max_time - min_time + 1);
        mean_val = float(num_1) * float(num_2)/T;
        denom = sqrt((float(num_1)-(float(num_1) * float(num_1)/T)) * (float(num_2)-(float(num_2) * float(num_2)/T)));
        /* Normalize the cross correlation */
        for (size_t i = 0; i < N; i++) {
            result[i] = (result[i] - mean_val) / denom;
        }
    } else if (normalization == "probability") {
        /* Normalize by the spike count. */
        float spike_count = 0.0;
        for (size_t i = 0; i < N; i++) {
            spike_count += result[i];
        }
        for (size_t i = 0; i < N; i++) {
            result[i] /= spike_count;
        }
    }
    
    return result;
}

/* Compute the cross-correlogram between two spike trains. Calculations and corrections are done using Equation 5.10
from "Analysis of Parallel Spike Trains" Sonja Grun & Stefan Rotter [eds] 2010 Volume 7, ISBN 978-1-4419-5674-3 */
// std::vector<float> crosscorrelogram_corr(const std::vector<int> spikes_1, const std::vector<int> spikes_2, float bin_width, float max_lag) {
//     int min_time, max_time;
//     float T, mean_val, denom;
//     size_t i;
//     size_t N = int(2 * max_lag / bin_width + 1);
//     std::vector<float> result(N);
//     size_t num_1 = spikes_1.size();
//     size_t num_2 = spikes_2.size();
//     min_time = min(spikes_1[0], spikes_2[0]);
//     max_time = max(spikes_1[num_1-1], spikes_2[num_2-1]);
//     T = float(max_time - min_time + 1);
//     mean_val = float(num_1) * float(num_2)/T;
//     denom = sqrt((float(num_1)-(float(num_1) * float(num_1)/T)) * (float(num_2)-(float(num_2) * float(num_2)/T)));
//     size_t j;
//     size_t j0 = 0;
//     for (i = 0; i < num_1; i++) {
//         j = j0;
//         while ((j < num_2) && (float(spikes_2[j] - spikes_1[i]) < -max_lag - bin_width/2.0)) {
//             j += 1;
//         }
//         j0 = j;
//         while ((j < num_2) && (float(spikes_2[j] - spikes_1[i]) < max_lag + bin_width/2.0)) {
//             result[int((float(spikes_2[j] - spikes_1[i]) + max_lag + 0.5*bin_width)/bin_width)] += 1.0;
//             j += 1;
//         }
//     }
//     /* Normalize the cross correlation */
//     for (size_t i = 0; i < N; i++) {
//         result[i] = (result[i] - mean_val) / denom;
//     }
//     return result;
// }

std::vector<float> autocorrelogram(const std::vector<int> spike_times, float bin_width, float max_lag) {
    float T, mean_val, denom, diff_time;
    size_t i;
    size_t N = int(max_lag / bin_width + 1);
    std::vector<float> result(N);
    size_t num_spikes = spike_times.size() - 1;
    T = float(spike_times[num_spikes-1] - spike_times[0] + 1);
    mean_val = float(num_spikes) * float(num_spikes)/T;
    denom = float(num_spikes) - mean_val;
    for (i = 0; i < num_spikes; i++) {
        diff_time = float(spike_times[i+1] - spike_times[i]);
        if (diff_time <= max_lag) {
            result[int(diff_time/bin_width)] += 1.0;
        }
    }
    /* Normalize the auto-correlation */
    for (size_t i = 0; i < N; i++) {
        result[i] = (result[i] - mean_val) / denom;
    }
    return result;
}