#include <iostream>
#include <iostream>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>

#include "simulation.h"

Simulation::Simulation(double L_x, int N_t, int N_x, double d_t, int max_i_t, double to_l, double n_u)
{
    this->Lx = L_x;
    this->Nt = N_t;
    this->Nx = N_x;
    this->dt = d_t;
    this->max_it = max_i_t;
    this->tol = to_l;
    this->nu = n_u;
    this->dx = Lx/(Nx-1);


    u = new double[Nx+1]();
    un = new double[Nx+1]();
    p = new double[Nx]();
    
}

Simulation::~Simulation()
{
    delete[] u;
    delete[] un;
    delete[] p;
}

void Simulation::vICs(double* &u, const double &dx, const int &Nx) {
    for (int i = 0; i < Nx+1 ; ++i) {
        u[i] = sin(M_PI*i*dx);
    }
}

void Simulation::advdiff(double* &u, double* &un, const double &dt, const double &dx, const double &nu, const int &Nx)
{
    for(int i = 1; i < Nx; ++i)
    {
        u[i] += dt*nu*(un[i+1]-2.0*un[i]+un[i-1])/(dx*dx);
        u[i] -= dt*(un[i]-un[i-1])/(2*dx);
    }
}

void Simulation::pcalc(double* &p, double* &u, const double &dx, const int &Nx, const int &max_it, const double &tol)
{
    double* pn = new double[Nx];
    for(int iter = 0; iter < max_it; ++iter)
    {
        copyarr(p,pn,Nx);
        for(int i = 1; i < Nx-1; ++i)
        {
            double b = -0.5*dx*(u[i+1]-u[i-1]);
            p[i] = 0.5*(p[i+1]+p[i-1])+b;
        }
        p[0] = p[Nx-2];
        p[Nx-1] = p[1];

        double err = 0.0;
        for(int i = 0; i < Nx ; ++i)
        {
            err = std::max(err, std::abs(pn[i]-p[i]));
        }
        if(err < tol)
        {
            break;
        }
    }
    delete[] pn;
}

void Simulation::pcorr(double* &u, double* &p, const double &dt, const double &dx, const int &Nx)
{
    for(int i = 1; i < Nx -1; ++i)
    {
        u[i] -= 2.0*dt*(p[i+1]-p[i])/dx;
    }
}

void Simulation::vbcs(double* &u, const int &Nx)
{
    u[0] = 1.0;
    u[Nx-1] = 1.0;

}

void Simulation::filewrite(const double* arr, const int &size, const std::string& filename) {
    std::ofstream outfile(filename);
    if (outfile.is_open()) {
        for (int i = 0; i < size; ++i) {
            outfile << arr[i] << "\n";
        }
        outfile.close();
    } else {
        std::cerr << "Unable to open file " << filename << "\n";
    }
}

void Simulation::copyarr(double* &ori, double* &cpy, const int &size)
{
    for(int i = 0; i < size; ++i)
    {
        cpy[i] = ori[i];
    }
}

void Simulation::run()
{
    vICs(u,dx,Nx);
    for(int n = 0; n < Nt; ++n)
    {
        copyarr(u,un,Nx);
        advdiff(u,un,dt,dx,nu,Nx);
        pcalc(p,u,dx,Nx,max_it,tol);
        pcorr(u,p,dt,dx,Nx);
        vbcs(u,Nx);
    }
    filewrite(u,Nx,"velocity.txt");
    filewrite(p,Nx,"pressure.txt");

}