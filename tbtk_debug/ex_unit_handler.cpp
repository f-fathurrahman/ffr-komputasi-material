#include "TBTK/TBTK.h"
#include "TBTK/UnitHandler.h"
#include <complex>
using namespace std;
using namespace TBTK;
int main(int argc, char **argv){
    //Initialize TBTK.
    Initialize();
    //Initialize the UnitHandler.
    UnitHandler::setScales(
        {"1 rad", "1 C", "1 pcs", "1 meV", "1 m", "1 K", "1 s"}
    );
    return 0;
}
