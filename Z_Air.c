/**
 * @file    Z_Air.c
 * @brief   Compute air refractivity
 *
 *
 */

#include <math.h>
#include <stdio.h>

//double rhocoeff = 1.0;
//double C_ls = 2.686777447e25; // Loschmidt constant
//double C_Na = 6.0221413e23; // Avogadro number

double Z_Air(double P, double T, double RH)
{
    double Z, Z0; //rho0;
    //double rho;
    //double tc; // C
    double rhocoeff1;
    //double An = 4.446e-6;
    //double Bn = 6.4e-13;
    //double Cn = -1.07e-16;
    double P0 = 101325.0;

    double TK;
    double f;   // enhancement factor
    double Psv; // water vapor saturation
    double xv;

    double A = 1.2378847e-5;
    double B = -1.9121316e-2;
    double C = 33.93711047;
    double D = -6.3431645e3;

    double alpha = 1.00062;
    double beta  = 3.14e-8;
    double gamma = 5.6e-7;

    double a0 = 1.58123e-6;
    double a1 = -2.9331e-8;
    double a2 = 1.1043e-10;
    //double b0 = 5.707e-6;
    double b1 = -2.051e-8;
    double c0 = 1.9898e-4;
    double c1 = -2.376e-6;
    double d  = 1.83e-11;
    double e  = -0.765e-8;

    double RH0 = 0.0;
    double T0, TK0, f0, Psv0, xv0;

    TK = T + 273.15;

    Z = 1.00001 - 5.8057e-9 * P + 2.6402e-16 * P * P - 3.3297e-7 * T +
        1.2420e-10 * P * T - 2.0158e-18 * P * P * T + 2.4925e-9 * T * T -
        6.2873e-13 * P * T * T + 5.4174e-21 * P * P * T * T - 3.5e-7 * RH -
        5.0e-9 * RH * RH;

    Psv = exp(A * TK * TK + B * TK + C + D / TK);
    printf("Water vapor saturation pressure Psv0 = %f Pa\n", Psv);
    f = alpha + beta * P + gamma * T * T;
    printf("enhancement factor f = %f\n", f);
    xv = RH * 0.01 * f * Psv / P;
    Z  = 1.0 -
        P / TK *
            (a0 + a1 * T + a2 * T * T + (c0 + b1 * T) * xv +
             (c0 + c1 * T) * xv * xv) +
        P * P / TK / TK * (d + e * xv * xv);

    Z0 = 1.00001 - 5.8057e-9 * P0 + 2.6402e-16 * P0 * P0;

    T0   = 0.0;
    TK0  = 273.15;
    Psv0 = exp(A * TK0 * TK0 + B * TK0 + C + D / TK0);
    printf("Water vapor saturation pressure Psv0 = %f Pa\n", Psv0);
    f0 = alpha + beta * P0 + gamma * T0 * T0;
    printf("enhancement factor f0 = %f\n", f0);
    xv0 = RH0 * 0.01 * f0 * Psv0 / P0;
    Z0  = 1.0 -
         P0 / TK0 *
             (a0 + a1 * T0 + a2 * T0 * T0 + (c0 + b1 * T0) * xv0 +
              (c0 + c1 * T0) * xv0 * xv0) +
         P0 * P0 / TK0 / TK0 * (d + e * xv0 * xv0);

    rhocoeff1 = 101325.0 / P;
    rhocoeff1 *= (T + 273.15) /
                 273.15; // particle count ratio : more particles at lower temp
    rhocoeff1 *= Z / Z0;

    //rho0 = C_ls / C_Na; // [mol.m^-3]

    //rhocoeff = rhocoeff1;

    printf("Z = %.8f   Z0 = %.8f\n", Z, Z0);

    //rho = rho0 / rhocoeff1;
    //rhocoeff = rhocoeff1;
    // * (1.0 + Bn/An*rho0 + Cn/An*rho0*rho0) / (1.0 + Bn/An*rho + Cn/An*rho*rho) / (1.0 + Bn/An*rho + Cn/An*rho*rho);

    return (Z);
}
