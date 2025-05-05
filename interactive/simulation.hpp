#include "gl_wrappers.hpp"
#include "parameters.hpp"

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
    uint32_t copy_r;
    uint32_t height_map;
    uint32_t configs_view;
    uint32_t modes;
    GLSLPrograms();
};

struct Histogram2D {
    IVec2 dimensions;
    Vec2 min_val, range;
    std::vector<float> arr;
};

class Simulation {
    GLSLPrograms m_programs;
    Frames m_frames;
    Histogram2D m_hist;
    std::vector<double> m_configs;  // Stores the Monte Carlo samples
    std::vector<double> m_slow_dst;
    // The above stores the transform matrix from going from position
    // to normal coordinates.
    std::vector<double> m_initial_wf;  // The initial wave function
    std::vector<double> sigma;
    public:
    Simulation(const SimParams &sim_params,
        int view_width, int view_height);
    void compute_configurations(SimParams &sim_params);
    const RenderTarget &render_view(const SimParams &sim_params);
    void reset_oscillator_count(int number_of_oscillators);
    void cursor_set_initial_wf(SimParams &sim_params, Vec2 cursor_pos);
};

};

#endif