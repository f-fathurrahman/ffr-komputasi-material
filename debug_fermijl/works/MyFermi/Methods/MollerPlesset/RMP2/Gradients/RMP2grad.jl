using GaussianBasis
using Molecules
using TensorOperations

function RMP2grad(x...)
    RMP2grad(Molecule(), x...)
end

function RMP2grad(mol::Molecule, x...)
    dtype = Options.get("deriv_type")
    if dtype == "analytic"
        throw(MyFermiException("Invalid or unsupported derivative type for RMP2: \"$dtype\""))
    elseif dtype == "findif"
        MyFermi.gradient_findif(MyFermi.MollerPlesset.RMP2, mol, x...)
    else
        throw(MyFermiException("Invalid or unsupported derivative type: \"$dtype\""))
    end
end

### Analytic graidents go here ###