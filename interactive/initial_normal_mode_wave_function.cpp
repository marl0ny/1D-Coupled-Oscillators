#include "initial_normal_mode_wave_function.hpp"


InitialNormalModeWaveFunction
::InitialNormalModeWaveFunction(size_t size) {
    this->size = size;
    for (int i = 0; i < MAX_SIZE; i++)
        s[i] = 1.0;
}

void InitialNormalModeWaveFunction::resize(size_t new_size) {
    this->size = new_size;
} 