#include "gl_wrappers.hpp"
#include "parameters.hpp"
#include "histogram.hpp"
#include "initial_normal_mode_wave_function.hpp"

#ifndef _SIM_2D_
#define _SIM_2D_

namespace sim_2d {

struct Frames {
    TextureParams view_tex_params;
    TextureParams configs_view_tex_params;
    TextureParams hist_tex_params;
    TextureParams initial_values_tex_params;
    Quad tmp;
    Quad initial_values;
    Quad hist;
    RenderTarget configs_view;
    RenderTarget view;
    WireFrame quad_wire_frame;
    Frames(
        const SimParams &sim_params,
        int view_width, int view_height);
};

struct GLSLPrograms {
    uint32_t scale;
    uint32_t copy;
    uint32_t height_map;
    uint32_t configs_view;
    uint32_t modes;
    GLSLPrograms();
};


class Simulation {
    GLSLPrograms m_programs;
    Frames m_frames;
    histogram::Histogram2D m_hist;
    std::vector<double> m_configs;  // Stores the Monte Carlo samples
    std::vector<double> m_slow_dst;
    // The above stores the transform matrix that takes the the position
    // of the oscillator in position coordinates and gives it in normal
    // coordinates.
    InitialNormalModeWaveFunction m_initial_wave_func;
    void compute_coherent_state_configurations(SimParams &sim_params);
    void compute_squeezed_state_configurations(SimParams &sim_params);
    void compute_stationary_state_configurations(SimParams &sim_params);
    void compute_single_excitations_configurations(SimParams &sim_params);
    public:
    Simulation(const SimParams &sim_params,
        int view_width, int view_height);
    void compute_configurations(SimParams &sim_params);
    const RenderTarget &render_view(const SimParams &sim_params);
    void reset_oscillator_count(int number_of_oscillators);
    void cursor_set_initial_wave_function(
        SimParams &sim_params, Vec2 cursor_pos);
    void set_relative_standard_deviation(float val);
};

};

#endif