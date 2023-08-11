#include "TBTK/Solver/Diagonalizer.h"
#include "TBTK/TBTK.h"
#include <complex>

using namespace std;
using namespace TBTK;

complex<double> i(0, 1);

int main(){

    // Initialize TBTK.
    Initialize();

    // Parameters.
    const unsigned int SIZE = 500;
    const double t = -1;
    const double mu = -1;

    // Set up the Model.
    Model model;
    for(unsigned int x = 0; x < SIZE; x++) {
        model << HoppingAmplitude(t, {x+1}, {x}) + HC;
    }
    model.construct();
    model.setChemicalPotential(mu);

    // Set up the Solver.
    Solver::Diagonalizer solver;
    solver.setModel(model);
    solver.run();
    
    cout << "Program finished normally\n";
    return 0;

}

