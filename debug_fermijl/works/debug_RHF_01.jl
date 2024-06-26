using MyFermi

@molecule {
    O        1.2091536548      1.7664118189     -0.0171613972
    H        2.1984800075      1.7977100627      0.0121161719
    H        0.9197881882      2.4580185570      0.6297938830
}

@set {
    basis sto-3g
    df false
    diis true
}

println("Options.Current:")
for k in keys(MyFermi.Options.Current)
    println("key = $k, val = $(MyFermi.Options.Current[k])")
end