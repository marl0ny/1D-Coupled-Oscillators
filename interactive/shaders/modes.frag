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
uniform int waveFunctionType;


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

/* Formula for the harmonic oscillator coherent state from
Shankar, pg. 610, equation 21.1.132 in exercise 21.2.18. */
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

/* Recursive Hermite polynomials.

See Shankar, pg. 195, 7.3.21 to obtain the base cases,
then 7.3.35 for the recursive relation itself.
*/
float hermite(int n, float x) {
    float prev1, prev2, res;
    for (int i = 0; i <= n; i++) {
        if (i == 0) {
            res = 1.0;
            prev2 = 1.0;
        } else if (i == 1) {
            res = 2.0*x;
            prev1 = 2.0*x;
        } else if (i > 1) {
            res = 2.0*(x*prev1 - float(i-1)*prev2);
            prev1 = prev2, prev1 = res;
        }
    }
    return res;
}

int factorial(int n) {
    int res = 1, prev = 1;
    for (int i = 0; i <= n; i++) {
        if (i > 1) {
            res = i*prev;
            prev = res;
        }
    }
    return res;
}

/* Energy eigenstates for the harmonic oscillator referenced from
Shankar, pg 195, eq 7.3.22.*/
complex stationary_state(
    int n, float x, float t,
    float m, float omega, float hbar
) {
    complex i = complex(0.0, 1.0);
    float norm_factor = pow(
        m*omega/(PI*hbar*pow(2.0, 2.0*float(n))
        *pow(float(factorial(n)), 2.0)), 0.25);
    return norm_factor*(
        exp(-0.5*m*omega*x*x/hbar - i*omega*(float(n) + 0.5)*t)
        *hermite(n, x*sqrt(m*omega/hbar)));
}

complex squeezedStatePhaseFactor(
    float x, float t, 
    float x0, float p0, float sigma0,
    float m, float omega, float hbar
) {
    float c = cos(omega*t);
    float s = sin(omega*t);
    float c2 = c*c;
    float c3 = c*c*c;
    float s2 = s*s;
    float hbar2 = hbar*hbar;
    float x2 = x*x;
    float x02 = x0*x0;
    float m2 = m*m;
    float m3 = m*m*m;
    float sigma4 = sigma0*sigma0*sigma0*sigma0;
    float p02 = p0*p0;
    float omega2 = omega*omega;
    float omega3 = omega*omega*omega;
    float omega_t = omega*t;
    float phase = (
        2.0*c3*m3*omega3*sigma4*x2 
        + c*m*omega*(hbar2*s2*x2 + hbar2*s2*x02 
                     - 4.0*m2*omega2*sigma4*x2 + 8.0*m*omega*p0*s*sigma4*x 
                     - 4.0*p02*s2*sigma4)/2.0
        + hbar2*s2*x0*(-m*omega*x + p0*s)
        )/(hbar*s*(4.0*c2*m2*omega2*sigma4 + hbar2*s2));
    float eps = 1e-30;
    if (abs(mod(omega_t, (2.0*PI))) < eps)
        return mul(cExp(complex(0.0, p0*x/hbar)), 
                   complex(1.0, -1.0)/sqrt(2.0));
    else if (abs(mod(omega_t, (PI))) < eps)
        mul(cExp(complex(0.0, -p0*x/hbar)), 
            complex(1.0, -1.0)/sqrt(2.0));
    return mul(cExp(complex(0.0, phase)),
               complex(1.0, -1.0)/sqrt(2.0));
}

complex squeezedState(
    float x, float t, 
    float x0, float p0, float sigma0,
    float m, float omega, float hbar, bool omit_phase
) {
    float c = cos(t*omega), s = sin(t*omega);
    float hbar2 = hbar*hbar;
    float m2 = m*m;
    float omega2 = omega*omega;
    float sigma4 = sigma0*sigma0*sigma0*sigma0;
    float xT = x0*c + p0*s/(m*omega);
    float sT = sqrt(hbar2*s*s + 4.0*m2*omega2*sigma4*c*c)
        /(2.0*m*omega*sigma0);
    complex phaseFactor = (colorPhase)?
        complex(1.0, 0.0):
        squeezedStatePhaseFactor(
            x, t, x0, p0, sigma0, m, omega, hbar);
    return exp(-0.25*pow((x - xT)/sT, 2.0))/sqrt(sT*sqrt(2.0*PI))
        *phaseFactor;
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