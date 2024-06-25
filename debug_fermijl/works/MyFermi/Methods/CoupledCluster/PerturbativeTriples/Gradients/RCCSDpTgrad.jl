using GaussianBasis
using Molecules
using TensorOperations

function RCCSDpTgrad(x...)
    RCCSDpTgrad(Molecule(), x...)
end

function RCCSDpTgrad(mol::Molecule, x...)
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