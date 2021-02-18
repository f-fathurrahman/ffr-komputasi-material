import numpy as np
from math import pi, sqrt

##############################################################
###     Everything connected with fcc crystal structure  #####
###     Could be generalized to an arbitrary crystal    #####
##############################################################
class FccLattice:
    "Handles fcc crystal. Could be generalized to arbitrary crystal."
    def __init__(self, LatConst):
        self.a0 = np.array([0.5*LatConst,0.5*LatConst,0])
        self.a1 = np.array([0.5*LatConst,0,0.5*LatConst])
        self.a2 = np.array([0,0.5*LatConst,0.5*LatConst])
        Vc = np.dot(np.cross(self.a0,self.a1),self.a2) # Volume
        self.Volume = abs(Vc)
        print("Volume is ", self.Volume)
        self.b0 = (2*pi/Vc)*np.cross(self.a1,self.a2)
        self.b1 = (2*pi/Vc)*np.cross(self.a2,self.a0)
        self.b2 = (2*pi/Vc)*np.cross(self.a0,self.a1)
        # Special points in Brillouin zone
        brs = 2*pi/LatConst
        self.GPoint = [0,0,0]
        self.LPoint = np.array([0.5,  0.5,  0.5])*brs
        self.KPoint = np.array([0.75, 0.75, 0])*brs
        self.XPoint = np.array([1.0,  0.0,  0])*brs
        self.WPoint = np.array([1.0,  0.5,  0])*brs
        
    def RMuffinTin(self):
        return 0.5*sqrt(np.dot(self.a0,self.a0)) # Spheres just touch

    def GenerateReciprocalVectors(self, q, CutOffK):
        # Many reciprocal vectors are generated and later only the shortest are used
        Kmesh0=[]
        for n in range(-q,q+1):
            for l in range(-q,q+1):
                for m in range(-q, q+1):
                    vec = n*self.b0+l*self.b1+m*self.b2
                    if np.dot(vec,vec) <= CutOffK**2:
                        Kmesh0.append(vec)
                    
        Kmesh0.sort(lambda x,y: cmp(np.dot(x,x), np.dot(y,y)))
        self.Km = np.array(Kmesh0)
        print("K-mesh size=", len(self.Km))
        
    def ChoosePointsInFBZ(self, nkp, type=0): # Chooses the path in the 1BZ we will use
        
        def kv0(iq, q):
            return (iq-int((q+1.5)/2)+1)/(q+0.0)
        
        if type==0: # Choose mesh in the 1BZ to cover the whole space - for SC calculation
            kp=[]
            for i0 in range(nkp):
                r0 = kv0(i0,nkp)
                print('r0=', r0)
                for i1 in range(nkp):
                    r1 = kv0(i1,nkp)
                    for i2 in range(nkp):
                        r2 = kv0(i2,nkp)
                        k = self.b0*r0+self.b1*r1+self.b2*r2
                        kp.append(k)
            print("Number of all k-points=", len(kp))

            kpc = []
            for k in kp:
                kpc.append(np.sort(k))
                
            # ChooseIrreducible k-points only
            # The function performs all symmetry operations of a cubic point-group to each k-point and
            # keeps only thos k-points which can not be obtained from another k-point by group operation.
            # These k-points are obviously irreducible.
            irkp = []       # temporary list where irreducible k points will be stored
            wkp  = []       # temporary weights
            while len(kpc)>0: # continues until all k-points are grouped into irreducible classes 
                tk = kpc[0]               # we concentrate on the k-point which is the first in the list
                irkp.append(tk)          # the first can be stored as irreducible
                wkp.append(0)            # and the weights for this irreducible k-point is set to zero
                # We go over 48 symmetry operations of cubic system:
                # Each wector component can change sign: 2^3=8 possibilities
                # All permutations of components: 3!=6
                # Since the operations are independent, we have 3!*2^3=48 operations == number of cubic point group operations
                for ix in [-1,1]:  # three loops for all possible sign changes 
                    for iy in [-1,1]:
                        for iz in [-1,1]:
                            nk = np.sort([ix*tk[0], iy*tk[1], iz*tk[2]]) # sorted so that we do not need to try all permutations
                            ii=0
                            while ii<len(kpc): # This permutation and sign change leads to some element still in the list of k-points?
                                diff = sum(abs(nk - kpc[ii]))
                                if diff<1e-6:
                                    del kpc[ii] # These two k-points are the same
                                    wkp[-1] += 1.
                                else:
                                    ii+=1

            # irreducible k-points are stored in the output vectors
            self.wkp = np.array(wkp)/np.sum(wkp)
            self.kp = np.array(irkp)

            print("Number of irreducible k points is: ", len(self.kp))
            #for ik,k in enumerate(self.kmesh):
            #    print "%10.6f"*3 % tuple(k), '  ', self.wkp[ik]
            
        else:        # Choose one particular path in the 1BZ - for plotting purposes
            nkp = 4*int(nkp/4.)+1
            print("nkp=", nkp)
            self.kp = np.zeros((nkp,3), dtype=float)
            N0=nkp/4

            self.Points = [('$\Gamma$', 0), ('$X$', N0), ('$L$', 2*N0), ('$\Gamma$', 3*N0), ('$K$', 4*N0)]
            for i in range(N0): self.kp[i,:]      = self.GPoint + (self.XPoint-self.GPoint)*i/(N0-0.)
            for i in range(N0): self.kp[N0+i,:]   = self.XPoint + (self.LPoint-self.XPoint)*i/(N0-0.)
            for i in range(N0): self.kp[N0*2+i,:] = self.LPoint + (self.GPoint-self.LPoint)*i/(N0-0.)
            for i in range(N0): self.kp[N0*3+i,:] = self.GPoint + (self.KPoint-self.GPoint)*i/(N0-0.)
            self.kp[4*N0] = self.KPoint
