#ifndef _CONFIGURATIONS_
#define _CONFIGURATIONS_

#include "initial_normal_mode_wave_function.hpp"

#include <vector>

namespace generate_configurations {

    void stationary_state(
        std::vector<double> &configurations, float &acceptance_rate, 
        double relative_delta, int number_of_mc_steps,
        int number_of_oscillators, const std::vector<double> &omega,
        const InitialNormalModeWaveFunction wave_func,
        double t, double m, double hbar);
    
    void coherent_state(
        std::vector<double> &configurations, float &acceptance_rate,
        double relative_delta, int number_of_mc_steps,
        int number_of_oscillators, const std::vector<double> &omega,
        const InitialNormalModeWaveFunction wave_func,
        double t, double m, double hbar);
    
    
    

    
};


#endif