using MyGaussianBasis
using MyMolecules
using TensorOperations

function RCCSDgrad(x...)
    RCCSDgrad(MyMolecule(), x...)
end

function RCCSDgrad(mol::MyMolecule, x...)
    dtype = Options.get("deriv_type")
    if dtype == "analytic"
        throw(MyFermiException("Invalid or unsupported derivative type for RCCSD: \"$dtype\""))
    elseif dtype == "findif"
        MyFermi.gradient_findif(MyFermi.CoupledCluster.RCCSD, mol, x...)
    else
        throw(MyFermiException("Invalid or unsupported derivative type: \"$dtype\""))
    end
end

### Analytic graidents go here ###