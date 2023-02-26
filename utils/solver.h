#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <chrono>
#include "constants.h"
#include "COO.h"


#define ind(x, y, z) ((x)+(y)*Nx+(z)*Nx*Ny)


/*-----------------------------------------------------------------------------
                    S   C   H   E   M   E
         *
         |
         |
        [2]
         |
         |
*-[3]---[0]---[1]-*
         |
         |
        [4]
         |
         |
         *



           |                   |                  |
           |                   |                  |
           |       A1          |       A3         |
           |                   |                  |
           |                   |                  |
           |                   |                  |
 MATRIX =  |-------------------|------------------|
           |                   |                  |
           |                   |                  |
           |         A2        |     A4           |
           |                   |                  |
           |                   |                  |



 now on zero time layer all p = 0


 ---------------------
 func S(s[i],s[j], p[i], p[j])
 if p[i]>p[j] return s[i]
 else return s[j]
 ----------------------
 k(s[i]) = 1

 ---------------------
check = true
for i = 0 : n * m
    i = 6
    [1]: 2 * Kx[i] * Kx[i+1] / ((hx^2) * (Kx[i] + Kx[i+1])) = Tau[1]
    [3]: 2 * Kx[i] * Kx[i-1] / ((hx^2) * (Kx[i] + Kx[i-1])) = Tau[3]
    [2]: 2 * Ky[i] * Ky[i-Nx] / ((hy^2) * (Ky[i] + Ky[i-Nx])) = Tau[2]
    [4]: 2 * Ky[i] * Ky[i+Nx] / ((hy^2) * (Ky[i] + Ky[i+Nx])) = Tau[4]
    [0]: -(Tau[1] + Tau[2] + Tau[3] + Tau[4]) = Tau[0]

    A = (Nx * Ny) x (Nx * Ny)

    A[i, i + 1] = Tau[1]
    A[i, i] = Tau[0]
    A[i, i - 1] = Tau[3]
    A[i, i - Nx] = Tau[2]
    A[i, i + Nx] = Tau[4]
    b = 0
    -*------------------*---------------------*-----------------*-------------------
    CALCULATING RESIDUAL:
    r_0[i] = phi[i] * (s[i]-s_prev[i])/dt + { Tau[1]*(p[i]-p[i+1])*S(i,i+1) + Tau[2]*(p[i]-p[i-nx])*S(i,i-nx)+
        +(Tau[3]*(p[i]-p[i-1])*S(i,i-1)+(Tau[4]*(p[i]-p[i+nx])*S(i,i+nx) }
    if i==well_index {r_0[i] -= k(s[i])*WI*(p_well - p[i]) }


    r_w[i] =  phi[i] * (-s[i]+s_initial[i])/dt - (+) { Tau[1]*(p[i]-p[i+1])*(1-S(i,i+1)) + Tau[2]*(p[i]-p[i-nx])*(1-S(i,i-nx))+
        +(Tau[3]*(p[i]-p[i-1])*(1-S(i,i-1))+(Tau[4]*(p[i]-p[i+nx])*(1-S(i,i+nx)) }
    if i==well_index {r_0[i] -= k(s[i])*WI*(p_well - p[i]) }

    if (abs(r_o[i])+abs(r_w[i]))>eps{
        check = false;
    }
    ---------------------------------------------------------------------
    MAKING J-MATRIX







 ------------------------------------------------------------------------
    if i < Nx:
        [2]: b = -2Ky[i] / (hy**2) * (dirichlet_up)
            A[i,i] += 2Ky[i] / (hy**2)
        [3]: b += -2Kx[i] / (hx**2) * (dirichlet_left)
            A[i,i] += 2Kx[i] / (hx**2)

                    S   C   H   E   M   E
--------------------------------------------------------------------------






*/



double k_o(double p_well, double p_i, double s_i);

double k_w(double p_well, double p_i, double s_i);

double S(double s_i, double s_j, double p_i, double p_j);

/* Calculate well index for cell */
double computeWellIndex(
        double kx, double ky
);


COO get_SLAE(
        COO A,
        double b[Nx * Ny],
        std::vector<double> kx,
        std::vector<double> ky,
        std::vector<double> kz,
        std::vector<double> s,
        std::vector<double> phi,
        std::vector<double> p
);

