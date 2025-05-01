#include "glfw_window.hpp"
#include "interactor.hpp"
#include "parameters.hpp"
#include "simulation.hpp"

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#include <emscripten/bind.h>
#endif
#include <functional>

#include "wasm_wrappers.hpp"

static std::function <void()> s_loop;
#ifdef __EMSCRIPTEN__
static void s_main_loop() {
    s_loop();
}
#endif

void oscillators(MainGLFWQuad main_render,
    int window_width, int window_height,
    sim_2d::SimParams &params,
    Interactor interactor) {
    sim_2d::Simulation sim(params, window_width, window_height);
    {
        s_sim_params_set = [&params, &sim](int c, Uniform u) {
            if (c == params.NUMBER_OF_OSCILLATORS) {
                sim.reset_oscillator_count(u.i32);
            }
            params.set(c, u);
        };
        s_sim_params_get = [&params](int c) -> Uniform {
            return params.get(c);
        };
        s_selection_set = [&params](
            int c, int val) {
            if (c == params.DISPLAY_TYPE)
                params.displayType.selected = val;
        };
    }
    size_t step_count = 0;
    s_loop = [&] {
        sim.compute_configurations(params);
        main_render.draw(sim.render_view(params));
        step_count++;
        params.t += params.dt;
        Vec2 pos = interactor.get_mouse_position();
        if (step_count % 3 == 0) {
            std::string text_content = "Acceptance rate: "
                + std::to_string(100.0*params.acceptanceRate) 
                + "% (33-50% ideal)";
            edit_label_display(params.ACCEPTANCE_RATE_LABEL, text_content);
        }
        auto poll_events = [&] {
            glfwPollEvents();
            interactor.click_update(main_render.get_window());
            if (pos.x > 0.0 && pos.x < 1.0 && 
                pos.y > 0.0 && pos.y < 1.0 && interactor.left_pressed()) {
                sim.cursor_set_initial_wf(params, pos);
            }
        };
        poll_events();
        glfwSwapBuffers(main_render.get_window());
    };
    #ifdef __EMSCRIPTEN__
    emscripten_set_main_loop(s_main_loop, 0, true);
    #else
    while (!glfwWindowShouldClose(main_render.get_window()))
        s_loop();
    #endif

}

int main(int argc, char *argv[]) {
    int window_width = 1500, window_height = 1500;
    // Construct the main window quad
    if (argc >= 3) {
        window_width = std::atoi(argv[1]);
        window_height = std::atoi(argv[2]);
    }
    auto main_quad = MainGLFWQuad(window_width, window_height);
    // Initialize Interactor instance
    Interactor interactor(main_quad.get_window());
    sim_2d::SimParams sim_params;
    oscillators(main_quad, 
        window_width, window_height, sim_params, interactor);
    return 1;
}