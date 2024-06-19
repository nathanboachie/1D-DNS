#include <iostream>
#include <iostream>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>

class Simulation{
public:

    Simulation(double L_x, int N_t, int N_x, double d_t, int max_i_t, double to_l, double n_u);
    ~Simulation();

    void vICs(double* &u, const double &dx, const int &Nx);
    void advdiff(double* &u, double* &un, const double &dt, const double &dx, const double &nu, const int &Nx);
    void pcalc(double* &p, double* &u, const double &dx, const int &Nx, const int &max_it, const double &tol);
    void pcorr(double* &u, double* &p, const double &dt, const double &dx, const int &Nx);
    void vbcs(double* &u, const int &Nx);
    void filewrite(const double* arr, const int &size, const std::string& filename);
    void run();

private:
    //Variables
    double Lx;
    double dt;
    double tol;
    double nu;
    double dx;
    int Nt;
    int Nx;
    int max_it;
    double* u = nullptr;
    double* un = nullptr;
    double* p = nullptr;

    //Helper methods
    void copyarr(double* &ori, double* &cpy, const int &size);
};