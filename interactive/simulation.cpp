#include "simulation.hpp"
#include "harmonic.hpp"
#include "configs_view.hpp"
#include "metropolis.hpp"
#include <complex>
#include <vector>

using namespace::sim_2d;

using std::complex;
using std::vector;

typedef std::vector<double> Arr1D;
typedef std::vector<int> ArrI1D;

#define PI 3.141592653589793


static void add_data_point(Histogram2D &hist, double x, double y, double val) {
    // x = min_val.x + range*(i/dimensions.x);
    // y = min_val.y + range*(j/dimensions.y);
    int i = double(hist.dimensions.x)*(x - hist.min_val.x)/hist.range.x;
    int j = double(hist.dimensions.y)*(y - hist.min_val.y)/hist.range.y;
    int ind = j*hist.dimensions[0] + i;
    if (ind >= 0 && ind < hist.arr.size())
        hist.arr[ind] += val;
}

static double coherent_state_prod(
    const Arr1D &x, double t, 
    const Arr1D &x0, const Arr1D &p0, double m,
    const Arr1D &omega, double hbar
) {
    std::complex<double> prod = 1.0;
    for (int i = 0; i < x.size(); i++)
        prod *= coherent_state(x[i], t, x0[i], p0[i], m, omega[i], hbar);
    return abs(prod)*abs(prod);
}

struct CoherentStateProdData {
    double t, m, hbar;
    Arr1D x0, p0, omega;
};

static double coherent_state_prod_dist_func(
    const Arr1D &x, void *data_ptr
) {
    CoherentStateProdData *data = (CoherentStateProdData *)data_ptr;
    return coherent_state_prod(
        x, data->t, data->x0, data->p0, data->m, data->omega, data->hbar);
}

static double stationary_states_prod(
    const Arr1D &x, double t,
    const ArrI1D &excitations,
    double m, const Arr1D &omega, double hbar
) {
    std::complex<double> prod = 1.0;
    for (int i = 0; i < excitations.size(); i++)
        prod *= 
            stationary_state(excitations[i], x[i], t, m, omega[i], hbar);
    return abs(prod)*abs(prod);
}

struct StationaryStatesProdData {
    double t, m, hbar;
    ArrI1D excitations;
    Arr1D omega;
};

static double stationary_states_prod_dist_func(
    const Arr1D &x, void *data_ptr
) {
    StationaryStatesProdData *data = (StationaryStatesProdData *)data_ptr;
    return stationary_states_prod(x, data->t, data->excitations,
         data->m, data->omega, data->hbar);
}

static double squeezed_state_prod(
    const Arr1D &x, double t, 
    const Arr1D &x0, const Arr1D &p0, const Arr1D &sigma0,
    double m, const Arr1D &omega, double hbar
) {
    std::complex<double> prod = 1.0;
    for (int i = 0; i < x.size(); i++)
        prod *= squeezed_state(
            x[i], t, x0[i], p0[i], 
            sigma0[i], m, omega[i], hbar);
    return abs(prod)*abs(prod);
}

struct SqueezedStateProdData {
    double t, m, hbar;
    Arr1D x0, p0, sigma0, omega;
};

static double squeezed_state_prod_dist_func(
    const Arr1D &x, void *data_ptr
) {
    SqueezedStateProdData *data = (SqueezedStateProdData *)data_ptr;
    return squeezed_state_prod(
        x, data->t, data->x0, data->p0, 
        data->sigma0, data->m, data->omega, data->hbar);
}

// struct squeezed_state_prod_dist()

static void c_sq_matrix_mul(
    double *dst, const double *m, const double *v, int n) {
    for (int i = 0; i < n; i++) {
        dst[i] = 0.0;
        for (int j = 0; j < n; j++) {
            dst[i] += m[i*n + j]*v[j];
        }
    }
}

static void c_copy(double *dst, const double *v, int n) {
    for (int i = 0; i < n; i++)
        dst[i] = v[i];
}

static WireFrame get_quad_wire_frame() {
    return WireFrame(
        {{"position", Attribute{
            3, GL_FLOAT, false,
            0, 0}}},
        {-1.0, -1.0, 0.0, -1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, -1.0, 0.0},
        {0, 1, 2, 0, 2, 3},
        WireFrame::TRIANGLES
    );
}

Frames::Frames(const SimParams &sim_params, 
    int view_width, int view_height):
    view_tex_params(
        {
            .format=GL_RGBA32F,
            .width=(uint32_t)view_width,
            .height=(uint32_t)view_height,
            .wrap_s=GL_CLAMP_TO_EDGE,
            .wrap_t=GL_CLAMP_TO_EDGE,
            .min_filter=GL_NEAREST,
            .mag_filter=GL_NEAREST,
        }
    ),
    configs_view_tex_params(
        {
            .format=GL_RGBA32F,
            .width=(uint32_t)view_width,
            .height=(uint32_t)view_height,
            // .width=(uint32_t)sim_params.numberOfOscillators,
            // .height=(uint32_t)sim_params.numberOfOscillators,
            .wrap_s=GL_REPEAT,
            .wrap_t=GL_REPEAT,
            .min_filter=GL_NEAREST,
            .mag_filter=GL_NEAREST,
        }
    ),
    hist_tex_params(
        {
            .format=GL_R32F,
            .width=(uint32_t)sim_params.numberOfOscillators,
            .height=(uint32_t)view_height/4,
            .wrap_s=GL_REPEAT,
            .wrap_t=GL_REPEAT,
            .min_filter=GL_NEAREST,
            .mag_filter=GL_NEAREST,
        }
    ),
    initial_values_tex_params(
        {
            .format=GL_RGBA32F,
            .width=(uint32_t)sim_params.numberOfOscillators,
            .height=(uint32_t)1,
            .wrap_s=GL_REPEAT,
            .wrap_t=GL_REPEAT,
            .min_filter=GL_NEAREST,
            .mag_filter=GL_NEAREST,
        }
    ),
    tmp(configs_view_tex_params),
    initial_values(initial_values_tex_params),
    hist(hist_tex_params),
    configs_view(configs_view_tex_params),
    view(view_tex_params),
    quad_wire_frame(get_quad_wire_frame())
    {
}

GLSLPrograms::GLSLPrograms() {
    this->copy = Quad::make_program_from_path("./shaders/util/copy.frag");
    this->copy_r = Quad::make_program_from_path("./shaders/util/copy-r.frag");
    this->scale = Quad::make_program_from_path("./shaders/util/scale.frag");
    this->height_map = Quad::make_program_from_path(
        "./shaders/util/height-map.frag");
    this->modes = Quad::make_program_from_path("./shaders/modes.frag");
    this->configs_view = make_program_from_paths(
        "./shaders/configs-view.vert", "./shaders/util/uniform-color.frag");
};

Simulation::Simulation(
    const SimParams &sim_params,
    int view_width, int view_height):
    m_programs(),
    m_frames(sim_params, view_width, view_height) {
    int n = sim_params.numberOfOscillators;
    m_configs = Arr1D(
        sim_params.numberOfMCSteps*n);
    m_slow_dst = Arr1D(n*n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            m_slow_dst[i*n + j] 
                = 2.0*sin(PI*(i + 1)*(j + 1)/(n + 1))/sqrt(2.0*(n + 1));
    m_initial_wf = Arr1D(n);
    Arr1D tmp (n);
    m_hist = {
        .dimensions=IVec2{.ind{n, (int)m_frames.hist_tex_params.height}},
        .min_val={.x=0.0, -20.0}, .range={.x=float(n), .y=40.0},
        .arr=std::vector<float>(n*m_frames.hist_tex_params.height)
    };
    for (int i = 0; i < n; i++)
        tmp[i] = 10.0*exp(-0.5*pow((double(i) - n/2.0)/(n*0.05), 2.0));
        // tmp[i] = 10.0*sin(4.0*PI*(i + 1)/(n + 1));
    c_sq_matrix_mul(&m_initial_wf[0], &m_slow_dst[0], &tmp[0], n);
}

void Simulation::compute_configurations(SimParams &sim_params) {
    int n = sim_params.numberOfOscillators;
    CoherentStateProdData data = {
        .t=sim_params.t, .hbar=1.0, .m=1.0,
        .x0=Arr1D(n), .p0=Arr1D(n), .omega=Arr1D(n)
    };
    auto delta = Arr1D(n);
    auto x = Arr1D(n);
    std::vector<float> initial_values_pixels {};

    for (int i = 0; i < n; i++) {
        data.x0[i] = m_initial_wf[i];
        // data.x0[i] = 5.0*exp(-0.5*pow((i - n/2)/(n/30), 2.0));
        // data.x0[i] = (i == 20)? 20.0: 0.0;
        data.p0[i] = 0.0;
        data.omega[i] = 2.0*sin(0.5*PI*(i + 1)/(n + 1));
        delta[i] = sim_params.relativeDelta/sqrt(data.omega[i]);
        if (sim_params.t < sim_params.dt) {
            x[i] = data.x0[i];
        } else {
            // x[i] = m_configs[
            //     sim_params.numberOfOscillators
            //     *(sim_params.numberOfMCSteps - 2) + i];
            x[i] = data.x0[i]*cos(data.omega[i]*sim_params.t);
        }
        initial_values_pixels.push_back(data.x0[i]);
        initial_values_pixels.push_back(data.p0[i]);
        initial_values_pixels.push_back(data.omega[i]);
        initial_values_pixels.push_back(0.0);
    }
    m_frames.initial_values.set_pixels(initial_values_pixels);
    // for (int i = 0; i < sim_params.numberOfMCSteps; i++)
    //     for (int k = 0; k < n; k++)
    //         m_configs[i*n + k] = data.x0[k];
    // for (int i = 0; i < sim_params.numberOfMCSteps; i++)
    //     for (int k = 0; k < n; k++)
    //         if (m_configs[i*n + k] != 0.0)
    //             printf("(%d, %d): %g\n", i, k, m_configs[i*n + k]);
    // printf("%g\n", sim_params.t);
    auto info = metropolis(
        m_configs, x, delta, 
        coherent_state_prod_dist_func,
        sim_params.numberOfMCSteps, (void *)&data);
    sim_params.acceptanceRate
         = float(info.accepted_count)/float(info.accepted_count + info.rejection_count);
}

const RenderTarget &Simulation::render_view(
    const SimParams &sim_params) {
    m_frames.configs_view.clear();
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    if (sim_params.displayType.selected == 2) {
        for (int i = 0; i < m_hist.arr.size(); i++)
            m_hist.arr[i] = 0.0;
        for (int k = 0; k < sim_params.numberOfMCSteps; k++) {
            int n = sim_params.numberOfOscillators;
            m_hist.min_val.y = -40.0;
            m_hist.range.y = 80.0;
            for (int j = 0; j < n; j++)
                add_data_point(
                    m_hist, 
                    double(j), m_configs[n*k + j] - 20.0,
                    sim_params.alphaBrightness);
            m_hist.min_val.y = -20.0;
            m_hist.range.y = 40.0;
            Arr1D res = Arr1D(n);
            c_sq_matrix_mul(
                &res[0], &m_slow_dst[0], 
                &m_configs[k*n], n);
            c_copy(&m_configs[k*n], &res[0], n);
            for (int j = 0; j < n; j++)
                add_data_point(
                    m_hist, 
                    double(j), m_configs[n*k + j] + 10.0,
                    sim_params.alphaBrightness);
        }
        m_frames.hist.set_pixels(&m_hist.arr[0]);
        m_frames.configs_view.draw(
            m_programs.height_map, 
            {{"tex", &m_frames.hist}}, 
            m_frames.quad_wire_frame
        );
    } else {
        WireFrame wire_frame = configs_view::get_configs_view_wire_frame(
            m_configs, 
            sim_params.numberOfMCSteps, 
            (sim_params.displayType.selected == 0)?
                configs_view::LINES_NO_ENDPOINTS:
                configs_view::DISCONNECTED_LINES);
        Vec3 c = sim_params.colorOfSamples2;
        m_frames.configs_view.draw(
            m_programs.configs_view,
            {
                {"scaleY", float(1.0F/40.0F)},
                {"color", Vec4{.r=c.r, c.g, c.b, sim_params.alphaBrightness}},
                {"yOffset", float(-0.5)},
            },
            wire_frame
        );
        for (int k = 0; k < sim_params.numberOfMCSteps; k++) {
            int n = sim_params.numberOfOscillators;
            Arr1D res = Arr1D(n);
            c_sq_matrix_mul(
                &res[0], &m_slow_dst[0], 
                &m_configs[k*n], n);
            c_copy(&m_configs[k*n], &res[0], n);
        }
    }
    m_frames.configs_view.draw(
        m_programs.modes,
        {
            {"t", sim_params.t},
            {"m", 1.0F},
            {"omega", 1.0F},
            {"hbar", 1.0F},
            {"scale", 40.0F},
            {"offset", float(-0.25F)},
            {"initialValuesTex", &m_frames.initial_values},
            {"brightness", sim_params.modesBrightness},
            {"colorPhase", int(sim_params.colorPhase)}
        },
        m_frames.quad_wire_frame
    );
    if (sim_params.displayType.selected == 0 ||
        sim_params.displayType.selected == 1) {
        WireFrame wire_frame = configs_view::get_configs_view_wire_frame(
            m_configs, sim_params.numberOfMCSteps,
            (sim_params.displayType.selected == 0)?
                configs_view::LINES_WITH_ZERO_ENDPOINTS:
                configs_view::DISCONNECTED_LINES);
        Vec3 c = sim_params.colorOfSamples1;
        m_frames.configs_view.draw(
            m_programs.configs_view,
            {
                {"scaleY", float(1.0F/20.0F)},
                {"color", Vec4{.r=c.r, c.g, c.b, sim_params.alphaBrightness}},
                {"yOffset", float(0.5)},
            },
            wire_frame
        );
    }
    glDisable(GL_BLEND);
    m_frames.view.draw(
        m_programs.scale,
        {
            {"tex", &m_frames.configs_view},
                {"scale", 1.0F}},
        m_frames.quad_wire_frame
    );
    return m_frames.view;
}

void Simulation::reset_oscillator_count(int number_of_oscillators) {
    m_frames.initial_values_tex_params = {
        .format=GL_RGBA32F,
        .width=(uint32_t)number_of_oscillators,
        .height=(uint32_t)1,
        .wrap_s=GL_REPEAT,
        .wrap_t=GL_REPEAT,
        .min_filter=GL_NEAREST,
        .mag_filter=GL_NEAREST,
    };
    m_frames.initial_values.reset(m_frames.initial_values_tex_params);
    int n = number_of_oscillators;
    m_slow_dst = Arr1D(n*n);
    m_initial_wf.resize(n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            m_slow_dst[i*n + j] 
                = 2.0*sin(PI*(i + 1)*(j + 1)/(n + 1))/sqrt(2.0*(n + 1));
    int view_height = m_frames.configs_view.texture_dimensions()[1];
    m_frames.hist_tex_params = {
        .format=GL_R32F,
        .width=(uint32_t)number_of_oscillators,
        .height=(uint32_t)view_height/4,
        .wrap_s=GL_REPEAT,
        .wrap_t=GL_REPEAT,
        .min_filter=GL_NEAREST,
        .mag_filter=GL_NEAREST,
    };
    m_frames.hist.reset(m_frames.hist_tex_params);
    m_hist = {
        .dimensions=IVec2{.ind{n, (int)m_frames.hist_tex_params.height}},
        .min_val={.x=0.0, -20.0}, .range={.x=float(n), .y=40.0},
        .arr=std::vector<float>(n*m_frames.hist_tex_params.height)
    };
}

void Simulation::cursor_set_initial_wf(SimParams &sim_params, Vec2 cursor_pos) {
    sim_params.t = 0.0;
    printf("%g, %g\n", cursor_pos.x, cursor_pos.y);
    int oscillator_pos = int(cursor_pos.x * sim_params.numberOfOscillators);
    if (cursor_pos.y < 0.5) {
        for (int i = 0; i < sim_params.numberOfOscillators; i++) {
            if (i == oscillator_pos)
                m_initial_wf[i] = 80.0*(cursor_pos.y - 0.25);
            else
                m_initial_wf[i] = 0.0;
        }
    } else {
        Arr1D tmp(sim_params.numberOfOscillators);
        for (int i = 0; i < sim_params.numberOfOscillators; i++) {
            int n = sim_params.numberOfOscillators;
            tmp[i] = 40.0*(cursor_pos.y - 0.75)*
                exp(-0.5*pow((double(i) - oscillator_pos)/(n*0.05), 2.0));
            // if (i == oscillator_pos)
            //     tmp[i] = 40.0*(cursor_pos.y - 0.75);
            // else
            //     tmp[i] = 0.0;
        }
        c_sq_matrix_mul(
            &m_initial_wf[0], &m_slow_dst[0], &tmp[0],
            sim_params.numberOfOscillators);
    }

}