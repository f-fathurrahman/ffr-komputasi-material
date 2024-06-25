using GaussianBasis
using Molecules
using TensorOperations

function UHFgrad(x...)
    UHFgrad(Molecule(), x...)
end

function UHFgrad(mol::Molecule, x...)
    dtype = Options.get("deriv_type")
    if dtype == "analytic"
        throw(MyFermiException("Invalid or unsupported derivative type for UHF: \"$dtype\""))
    elseif dtype == "findif"
        MyFermi.gradient_findif(MyFermi.HartreeFock.UHF, mol, x...)
    else
        throw(MyFermiException("Invalid or unsupported derivative type: \"$dtype\""))
    end
end

### Analytic graidents go here ###