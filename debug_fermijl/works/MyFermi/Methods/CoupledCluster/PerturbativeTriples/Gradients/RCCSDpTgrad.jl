using MyGaussianBasis
using MyMolecules
using TensorOperations

function RCCSDpTgrad(x...)
    RCCSDpTgrad(MyMolecule(), x...)
end

function RCCSDpTgrad(mol::MyMolecule, x...)
    dtype = Options.get("deriv_type")
    if dtype == "analytic"
        throw(MyFermiException("Invalid or unsupported derivative type for RCCSDpT: \"$dtype\""))
    elseif dtype == "findif"
        MyFermi.gradient_findif(MyFermi.CoupledCluster.RCCSDpT, mol, x...)
    else
        throw(MyFermiException("Invalid or unsupported derivative type: \"$dtype\""))
    end
end

### Analytic graidents go here ###