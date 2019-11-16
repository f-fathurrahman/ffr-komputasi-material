from ase.io import Trajectory
from ase.units import Ry
import sys

filnam = sys.argv[1]
basnam = filnam.replace(".traj","")
filout = open(basnam + ".xsf", "w")

traj = Trajectory(sys.argv[1])

Nimages = len(traj)
cell = traj[0].cell

filout.write("ANIMSTEPS {}\n".format(Nimages))
filout.write("CRYSTAL\n")
filout.write("PRIMVEC\n")
filout.write("{:18.10f} {:18.10f} {:18.10f}\n".format(cell[0,0], cell[0,1], cell[0,2]))
filout.write("{:18.10f} {:18.10f} {:18.10f}\n".format(cell[1,0], cell[1,1], cell[1,2]))
filout.write("{:18.10f} {:18.10f} {:18.10f}\n".format(cell[2,0], cell[2,1], cell[2,2]))

for idx_img in range(Nimages):

    atoms = traj[idx_img]
    forces_Ry = atoms.get_forces()/Ry # convert forces to Ry

    filout.write("PRIMCOORD {}\n".format(idx_img+1))
    filout.write("{} 1\n".format(len(atoms)))

    for i,a in enumerate(atoms):
        filout.write("{} {:18.10f} {:18.10f} {:18.10f}".format(a.symbol, a.position[0], a.position[1], a.position[2]))
        filout.write("{:18.10f} {:18.10f} {:18.10f}\n".format(forces_Ry[i][0], forces_Ry[i][1], forces_Ry[i][2]))


#for i,atoms in enumerate(traj):
#    print("i = ", i + 1)
#    filtmp = "TEMP_COORD_" + str(i+1) + ".xyz"
#    atoms.write(filtmp)
#    if i == 0:
#        os.system("cat " + filtmp + " > " + filout)
#    else:
#        os.system("cat " + filtmp + " >> " + filout )
#    os.system("rm " + filtmp)
#
#print("Done")

filout.close()