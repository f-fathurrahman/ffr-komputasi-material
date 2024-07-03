using MyGaussianBasis
using MyMolecules
using TensorOperations

function UHFgrad(x...)
    UHFgrad(MyMolecule(), x...)
end

function UHFgrad(mol::MyMolecule, x...)
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