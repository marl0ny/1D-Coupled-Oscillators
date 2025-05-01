/* See pg. 429 to 430 of "An Introduction Computer Simulation Methods"
by Harvey Gould, Jan Tobochnik, and Wolfgang Christian for a general 
outline of the Metropolis algorithm. */
#include <vector>


struct MetropolisResultInfo {
    int accepted_count, rejection_count;
};

MetropolisResultInfo metropolis(
    std::vector<double> &configs, 
    const std::vector<double> &x0,
    const std::vector<double> &delta,
    double (* dist_func)(const std::vector<double> &x, void *params),
    int steps, void *params);
