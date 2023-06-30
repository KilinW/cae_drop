/**
 * Copyright (c) 2022 Raspberry Pi (Trading) Ltd.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <tuple>
#include "pico/stdlib.h"
#include "pico/cyw43_arch.h"
#include "scicpp/core.hpp"

using namespace std;

double solution(double t)
{

    int C = 1;
    return (C * std::exp(t) - t - 1);
}

std::tuple<std::vector<double>, std::vector<double>> methodRuneKutt(float (*function)(float t, float para), double runSize)
{

    std::vector<double> diffEq;
    std::vector<double> time;
    double xi, k1, k2, k3, k4, k5, k6;
    float dt = 0.5;
    float t = 0;
    double x0 = solution(0);

    time.push_back(0);
    diffEq.push_back(x0);

    for (int i = 1; i < runSize; i++)
    {
        k1 = dt * function(t, diffEq[i - 1]);
        k2 = dt * function(t + 0.5 * dt, diffEq[i - 1] + 0.5 * k1);
        k3 = dt * function(t + 0.5 * dt, diffEq[i - 1] + 0.5 * k2);
        k4 = dt * function(t + dt, diffEq[i - 1] + k3);
        xi = diffEq[i-1] + 1.0/6.0*(k1 + 2*k2 + 2*k3 + k4 );

        diffEq.push_back(xi);
        time.push_back(i);
        t += dt;
    }

    return std::make_tuple (time, diffEq);
}


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
    float cos, sin, cosh, sinh;
    float cos_L, sin_L, cosh_L, sinh_L;
    float cos_2L, sin_2L, cosh_2L, sinh_2L;
    float sigma;
    float A1;
    float eq_solver_arr[5]; // [ Cp, R, vphi, zeta, omega ]
    float eta1, eta1dot, v2, eta1ddot,v2dot;
    float eq_solver_arr2[3];

    void cal_properties();
    float get_A1();
    float det_A();
    float phi(float x);
    float d_phi(float x);
    float* cantilever_actuator_eq_solver(float t,float in_[],float para[]);
    void voltage(float para[]);
    void time_span(float time_span,float step);
    void tip_mass(float mass);
    void set_force(float force_func);
    void total_damping(float damping);
    void load_resistance(float resistance);
    void piezo_coupling(float coupling);
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
    alpha = 1.0;//scicpp::units::quantity::root(det_A());
        // ??? root(self.det_A, 0).x[0]
        // https://github.com/tvanderbruggen/Scicpp/tree/master
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
    if (debug == true)
        cout << "Natural frequency omega: {self.alpha**2*math.sqrt(self.EI1/self.m1)}\n";
    return A1;
};

float piezo_film::det_A(){
    cos = std::cos(alpha*L1);
    sin = std::sin(alpha*L1);
    cosh = std::cosh(alpha*L1);
    sinh = std::sinh(alpha*L1);
    det_value = (1 + cos*cosh
                +(M1*alpha/m1)*(cos*sinh-sin*cosh)
                -(J1*pow(alpha,3)/m1)*(cosh*sin+sinh*cos)
                -(M1*J1*pow(alpha,4)/m1*m1)*(cos*cosh-1));
    return det_value;
};

float piezo_film::phi(float x){
    return A1*(std::cos(alpha*x)
            -std::cosh(alpha*x)
            +sigma*std::sin(alpha*x)
            -sigma*std::sinh(alpha*x));
};

float piezo_film::d_phi(float x){
    return A1*alpha*(-std::sin(alpha*x)
                    -std::sinh(alpha*x)
                    +sigma*(std::cos(alpha*x) 
                            -std::cosh(alpha*x)));
};

float* piezo_film::cantilever_actuator_eq_solver(float t,float in_[],float para[]){
    eq_solver_arr[0] = *para; // [ Cp, R, vphi, zeta, omega ]
    eta1 = in_[ 0 ];
    eta1dot = in_[ 1 ];
    v2 = in_[ 2 ];
    eta1ddot = -2 * eq_solver_arr[3] * eq_solver_arr[4] * eta1dot - eq_solver_arr[4]*eq_solver_arr[4] * eta1 - eq_solver_arr[2]*v2 - force;
    v2dot = -1 / ( eq_solver_arr[0]*eq_solver_arr[1] ) * v2 + eq_solver_arr[2]/eq_solver_arr[0]*eta1dot;
    eq_solver_arr2[0] = eta1dot;
    eq_solver_arr2[1] = eta1ddot;
    eq_solver_arr2[2] = v2dot;
    return eq_solver_arr2;
};

void piezo_film::voltage(float para[]){
    eq_solver_arr[2] = vtheta*(d_phi(Lp2) - d_phi(Lp1));
    para[0] = eq_solver_arr[0];
    para[1] = eq_solver_arr[1];
    para[2] = eq_solver_arr[2];
    para[3] = eq_solver_arr[3];
    para[4] = alpha*alpha*sqrt(EI1/m1);
    std::tuple<std::vector<double>, std::vector<double>>  diffRuneKutt = methodRuneKutt(cantilever_actuator_eq_solver, 5);
    std::vector<double> time = std::get<0>(diffRuneKutt);
    std::vector<double> solutionRuneKutt = std::get<1>(diffRuneKutt);
    return std::make_tuple (time, solutionRuneKutt);
    // return solve_ivp ????
    // https://medium.com/geekculture/numerical-computation-in-c-part-2-33dcdbe6c726
};

void piezo_film::time_span(float time_span,float step){
    time_end = time_span;
    time_step = step;
};
void piezo_film::tip_mass(float mass){
    M1 = mass;
};
void piezo_film::set_force(float force_func){
    force = force_func;
};
void piezo_film::total_damping(float damping){
    zeta = damping;
};
void piezo_film::load_resistance(float resistance){
    R = resistance;
};
void piezo_film::piezo_coupling(float coupling){
    vtheta = coupling;
};

int main() {
    stdio_init_all();
    if (cyw43_arch_init()) {
        printf("Wi-Fi init failed");
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

