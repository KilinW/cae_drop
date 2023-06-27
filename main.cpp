/**
 * Copyright (c) 2022 Raspberry Pi (Trading) Ltd.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

#include <stdio.h>
#include <math.h>
#include <iostream>
#include "pico/stdlib.h"
#include "pico/cyw43_arch.h"

using namespace std;

class piezo_film{
public:
    piezo_film(); //init
    bool debug;
    // 1 # cantilever film dimensions
    float M1;           //tip mass in kg (referring to raindrop)
    float J1;           //tip mass moment of inertia
    long int Ys;        //substrate young's modulus
    int rho_s;          //substarte density
    unsigned int Yp;        //piezo young's moudulus
    int rho_p;          //piezo density
    float h1s;          //substrate thickness(m)
    float h1p;          //piezo thickness(m)
    float L1;           //ubstrate length(m)
    float Lp1;          //0 means the piezo cover from the fixed end
    float Lp2;          //piezo length
    float b1s;          //substrate width
    float b1p;          //piezo width
    float b1pp;         //electrode width
    // 2 # piezoelectric parameters
    float d31;          //piezoelectric constant(m/V)
    float vtheta;       //piezoelectric coupling coefficient, vtheta = Yp*d31*b1pp*h1pc
    float epsilon;      //absolute permittivity 3500*8.854*1e-12 # F/m
    float Cp;           //Capaticance of the piezo layer
    // 4 # external circuit parameters
    long int R;         //external circuit load resistance
    float zeta;         //damping ratio (Including internal damping and air damping)
    float force;
    float time_end;
    int time_step;
    // 3 # cantilever mechanical parameters
    float m1;           //unit mass of the cantilever
    float n;            //stiffness ratio of substrate to piezo modulus
    float h1pa;         //piezo layer surface distance to neutral axis
    float h1sa;         //substrate layer surface distance to neutral axis
    float h1pc;         
    float h1a;          //distance from neutral axis to bottom of substrate
    float h1b;          //distance from neutral axis to top of substrate
    float h1c;          //distance from neutral axis to top of piezo layer
    float EI1;          //bending stiffness of the cantilever

    float alpha;
    float det_value;
    float cos;
    float sin;
    float cosh;
    float sinh;
    float cos_L;
    float sin_L;
    float cosh_L;
    float sinh_L;
    float cos_2L;
    float sin_2L;
    float cosh_2L;
    float sinh_2L;
    float sigma;
    float A1;    

    void cal_properties();
    float get_A1();
    float det_A();
    float phi(float x);
    float d_phi(float x);
};

piezo_film::piezo_film(){
    debug = false;
    M1 = 11.0 * 1e-6;
    J1 = M1 * 25.0 * 1e-6;
    Ys = 0.5 * 1e9;
    rho_s = 1390;
    Yp = 2.45 * 1e9;
    rho_p = 1780;
    h1s = 0.5 * 1e-3;
    h1p = 0.028 * 1e-3;
    L1 = 20 * 1e-3;
    Lp1 = 0 * L1;
    Lp2 = 1 * L1;
    b1s = 13 * 1e-3;
    b1p = b1s;
    b1pp = b1p;
    d31 = -190 * 1e-12;
    vtheta = 0.85*1e-7;
    epsilon = 15.93 * 1e-9;
    Cp = 0.5 * 1e-9;
    R = 1 * 1e6;
    zeta = 0.043;
    force = 0.08*24;
    time_end = 0.2;
    time_step = 1000;
    cal_properties();
}

void piezo_film::cal_properties(){
    m1 = rho_s*b1s*h1s + rho_p*b1p*h1p;
    n = Ys * b1s / ( Yp *b1p );
    h1pa = ( h1p*h1p + 2*n*h1p*h1s + n * h1s*h1s ) / 2 / ( h1p + n*h1s );
    h1sa = ( h1p*h1p + 2*h1p*h1s + n * h1s*h1s ) / 2 / ( h1p + n*h1s );
    h1pc = n * h1s * ( h1p+h1s ) / 2 / ( h1p + n*h1s );
    h1a = -1 * h1sa;
    h1b = h1pa - h1p;
    h1c = h1pa;
    EI1 = b1s/3*Ys*(powf(h1b,3)-powf(h1a,3))+b1p/3*Yp*(powf(h1c,3)-powf(h1b,3));
    get_A1();
}

float piezo_film::get_A1(){
    alpha = 1.0;                                            // ??
    float get_A1_var = alpha*L1;
    cos_L = std::cos(get_A1_var);
    sin_L = std::sin(get_A1_var);
    cosh_L = std::cosh(get_A1_var);
    sinh_L = std::sinh(get_A1_var);
    cos_2L = std::cos(2*get_A1_var);
    sin_2L = std::sin(2*get_A1_var);
    cosh_2L = std::cosh(2*get_A1_var);
    sinh_2L = std::sinh(2*get_A1_var);
    sigma = (sin_L - sinh_L + (alpha*M1/m1)*(cos_L - cosh_L))/(cos_L + cosh_L - (alpha*M1/m1)*(sin_L - sinh_L));
    A1 = sqrt((4*alpha)/((-sigma*sigma*sin_2L
                            +sigma*sigma*sinh_2L
                            -4*sigma*sigma*cos_L*sinh_L
                            +4*(sigma*sigma+1)*sin_L*cosh_L
                            -2*sigma*cos_2L
                            +2*sigma*cosh_2L
                            +8*sigma*sin_L*sinh_L
                            +4*alpha*L1
                            +sin_2L
                            +sinh_2L
                            +4*cos_L*sinh_L)*m1));
    if (debug ==true)
        cout << "Natural frequency omega: {self.alpha**2*math.sqrt(self.EI1/self.m1)}\n";
    return A1;
}

float piezo_film::det_A(){
                                                            // ??
    return det_value;
}

float piezo_film::phi(float x){
    return A1*(std::cos(alpha*x)
            -std::cosh(alpha*x)
            +sigma*std::sin(alpha*x)
            -sigma*std::sinh(alpha*x));
};

float piezo_film::d_phi(float x){
    return A1*alpha*(-std::sin(alpha*x)
            -std::sinh(alpha*x)
            +sigma*std::cos(alpha*x) 
            -sigma*std::cosh(alpha*x));
};

int main() {
    stdio_init_all();
    if (cyw43_arch_init()) {
        printf("Wi-Fi init failed");
        //return -1;
    }
    piezo_film A;
    while (true) {
        printf("Blinking!\r\n");
        cyw43_arch_gpio_put(CYW43_WL_GPIO_LED_PIN, 1);
        sleep_ms(250);
        cyw43_arch_gpio_put(CYW43_WL_GPIO_LED_PIN, 0);
        sleep_ms(250);
    }
    return 0;
}

