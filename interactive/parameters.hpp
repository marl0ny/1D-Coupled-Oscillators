#include "gl_wrappers.hpp"

namespace sim_2d {

#ifndef _PARAMETERS_
#define _PARAMETERS_

struct Button {};

typedef std::string Label;

typedef std::vector<std::string> EntryBoxes;

struct SelectionList {
    int selected;
    std::vector<std::string> options;
};

struct SimParams {

struct LineDivider {};
    int stepsPerFrame = (int)(1);
    float dt = (float)(0.1F);
    float t = (float)(0.0F);
    int numberOfOscillators = (int)(64);
    LineDivider lineDivMonteCarlo = LineDivider{};
    Label labelMonteCarlo = Label{};
    float relativeDelta = (float)(0.5F);
    Label acceptanceRateLabel = Label{};
    float acceptanceRate = (float)(0.0F);
    int numberOfMCSteps = (int)(1000);
    LineDivider lineDivSampleColor = LineDivider{};
    Label labelSamples = Label{};
    float alphaBrightness = (float)(0.01F);
    Vec3 colorOfSamples1 = (Vec3)(Vec3 {.ind={0.0, 0.85, 1.0}});
    Vec3 colorOfSamples2 = (Vec3)(Vec3 {.ind={0.0, 1.0, 0.0}});
    SelectionList displayType = SelectionList{0, {"Lines", "Scatter"}};
    LineDivider lineDivNormalModeWaveFunc = LineDivider{};
    Label labelNormalModeWaveFunc = Label{};
    bool colorPhase = (bool)(false);
    float modesBrightness = (float)(1.25F);
    enum {
        STEPS_PER_FRAME=0,
        DT=1,
        T=2,
        NUMBER_OF_OSCILLATORS=3,
        LINE_DIV_MONTE_CARLO=4,
        LABEL_MONTE_CARLO=5,
        RELATIVE_DELTA=6,
        ACCEPTANCE_RATE_LABEL=7,
        ACCEPTANCE_RATE=8,
        NUMBER_OF_M_C_STEPS=9,
        LINE_DIV_SAMPLE_COLOR=10,
        LABEL_SAMPLES=11,
        ALPHA_BRIGHTNESS=12,
        COLOR_OF_SAMPLES1=13,
        COLOR_OF_SAMPLES2=14,
        DISPLAY_TYPE=15,
        LINE_DIV_NORMAL_MODE_WAVE_FUNC=16,
        LABEL_NORMAL_MODE_WAVE_FUNC=17,
        COLOR_PHASE=18,
        MODES_BRIGHTNESS=19,
    };
    void set(int enum_val, Uniform val) {
        switch(enum_val) {
            case STEPS_PER_FRAME:
            stepsPerFrame = val.i32;
            break;
            case DT:
            dt = val.f32;
            break;
            case T:
            t = val.f32;
            break;
            case NUMBER_OF_OSCILLATORS:
            numberOfOscillators = val.i32;
            break;
            case RELATIVE_DELTA:
            relativeDelta = val.f32;
            break;
            case ACCEPTANCE_RATE:
            acceptanceRate = val.f32;
            break;
            case NUMBER_OF_M_C_STEPS:
            numberOfMCSteps = val.i32;
            break;
            case ALPHA_BRIGHTNESS:
            alphaBrightness = val.f32;
            break;
            case COLOR_OF_SAMPLES1:
            colorOfSamples1 = val.vec3;
            break;
            case COLOR_OF_SAMPLES2:
            colorOfSamples2 = val.vec3;
            break;
            case COLOR_PHASE:
            colorPhase = val.b32;
            break;
            case MODES_BRIGHTNESS:
            modesBrightness = val.f32;
            break;
        }
    }
    Uniform get(int enum_val) const {
        switch(enum_val) {
            case STEPS_PER_FRAME:
            return {(int)stepsPerFrame};
            case DT:
            return {(float)dt};
            case T:
            return {(float)t};
            case NUMBER_OF_OSCILLATORS:
            return {(int)numberOfOscillators};
            case RELATIVE_DELTA:
            return {(float)relativeDelta};
            case ACCEPTANCE_RATE:
            return {(float)acceptanceRate};
            case NUMBER_OF_M_C_STEPS:
            return {(int)numberOfMCSteps};
            case ALPHA_BRIGHTNESS:
            return {(float)alphaBrightness};
            case COLOR_OF_SAMPLES1:
            return {(Vec3)colorOfSamples1};
            case COLOR_OF_SAMPLES2:
            return {(Vec3)colorOfSamples2};
            case COLOR_PHASE:
            return {(bool)colorPhase};
            case MODES_BRIGHTNESS:
            return {(float)modesBrightness};
        }
        return Uniform(0);
    }
    void set(int enum_val, int index, std::string val) {
        switch(enum_val) {
        }
    }
};
#endif
}
