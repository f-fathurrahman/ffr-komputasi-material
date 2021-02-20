# Atomic charge for a good starting guess
def atom_cmpb(x,y):
    "Comparison function for sorting of bound states"
    if abs(x[2]-y[2])>1e-4:
        return cmp(x[2],y[2])
    else:
        return cmp(x[1],y[1])
