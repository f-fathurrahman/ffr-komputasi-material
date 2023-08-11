#include "TBTK/Streams.h"
#include "TBTK/TBTK.h"
#include <complex>

using namespace std;
using namespace TBTK;

int main(int argc, char **argv){
    //Initialize TBTK.
    Initialize();
    Streams::out << "Hello quantum world!\n";
    return 0;
}


