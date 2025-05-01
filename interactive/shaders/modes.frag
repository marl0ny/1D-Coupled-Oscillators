#if (__VERSION__ >= 330) || (defined(GL_ES) && __VERSION__ >= 300)
#define texture2D texture
#else
#define texture texture2D
#endif

#if (__VERSION__ > 120) || defined(GL_ES)
precision highp float;
#endif

#if __VERSION__ <= 120
varying vec2 UV;
#define fragColor gl_FragColor
#else
in vec2 UV;
out vec4 fragColor;
#endif

uniform int numberOfOscillators;
uniform bool colorPhase;
uniform float t;
uniform float m;
uniform float omega;
uniform float hbar;
uniform float scale;
uniform float offset;
uniform float brightness;
uniform sampler2D initialValuesTex;


#define complex vec2

const float PI = 3.141592653589793;


complex expI(float angle) {
    return complex(cos(angle), sin(angle));
}

complex cExp(complex z) {
    return exp(z.r)*complex(cos(z.y), sin(z.y));
}

complex mul(complex z, complex w) {
    return complex(z.x*w.x - z.y*w.y, z.x*w.y + z.y*w.x);
}

complex conj(complex z) {
    return complex(z[0], -z[1]);
}

/* See Shankar, pg. 610 */
complex coherentState(
    float x, float t, 
    float x0, float p0,
    float m, float omega, float hbar
) {
    complex z = complex(sqrt((m*omega)/(2.0*hbar))*x0,
                        sqrt(1.0/(2.0*m*omega*hbar))*p0);
    z = mul(z, expI(-omega*t));
    return pow(m*omega/(PI*hbar), 0.25)*(
        cExp(- z.r*z
             - complex(m*omega/(2.0*hbar)*x*x, 0.0)
             + sqrt(2.0*m*omega/hbar)*x*z));
}


vec3 argumentToColor(float argVal) {
    float maxCol = 1.0;
    float minCol = 50.0/255.0;
    float colRange = maxCol - minCol;
    if (argVal <= PI/3.0 && argVal >= 0.0) {
        return vec3(maxCol,
                    minCol + colRange*argVal/(PI/3.0), minCol);
    } else if (argVal > PI/3.0 && argVal <= 2.0*PI/3.0){
        return vec3(maxCol - colRange*(argVal - PI/3.0)/(PI/3.0),
                    maxCol, minCol);
    } else if (argVal > 2.0*PI/3.0 && argVal <= PI){
        return vec3(minCol, maxCol,
                    minCol + colRange*(argVal - 2.0*PI/3.0)/(PI/3.0));
    } else if (argVal < 0.0 && argVal > -PI/3.0){
        return vec3(maxCol, minCol,
                    minCol - colRange*argVal/(PI/3.0));
    } else if (argVal <= -PI/3.0 && argVal > -2.0*PI/3.0){
        return vec3(maxCol + (colRange*(argVal + PI/3.0)/(PI/3.0)),
                    minCol, maxCol);
    } else if (argVal <= -2.0*PI/3.0 && argVal >= -PI){
        return vec3(minCol,
                    minCol - (colRange*(argVal + 2.0*PI/3.0)/(PI/3.0)), 
                    maxCol);
    }
    else {
        return vec3(minCol, maxCol, maxCol);
    }
}

void main() {
    vec4 initialXPOmega = texture2D(initialValuesTex, UV);
    float x = (2.0*(UV[1] - 0.5 - offset))*scale;
    float x0 = initialXPOmega[0];
    float p0 = initialXPOmega[1];
    float omega = initialXPOmega[2];
    complex amplitude = coherentState(x, t, x0, p0, m, omega, hbar);
    if (colorPhase)
        fragColor = brightness*vec4(
            argumentToColor(atan(amplitude.y, amplitude.x)),
            dot(amplitude, amplitude)
        );
    else
        fragColor = brightness*vec4(vec3(1.0), dot(amplitude, amplitude));
}