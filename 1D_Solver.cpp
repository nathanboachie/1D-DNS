#include <iostream>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include "simulation.h"

int main() {
    // Create an instance of the Simulation class and run the simulation
    Simulation sim(2.0, 400, 44, 0.001, 100, 0.0004, 0.01);
    sim.run();
    return 0;
}