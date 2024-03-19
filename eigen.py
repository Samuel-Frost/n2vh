from ase.io import read
import numpy as np
from numpy.linalg import norm

initial = read('initial.xyz')
displacement = read('displacement.xyz')
c = 299792458 
prob = 2.95e-05

int_pos = (initial.get_positions())
disp_pos = np.array(displacement.get_positions())
dif = int_pos - disp_pos
dif /= 7.14
wavenumbers = np.loadtxt('wavenumbers')
k = 0 # for atom
atoms = np.loadtxt('phonon')
# 192 modes each with 64 ions :)
index = [3, 5, 7]
tot = np.array([0.0, 0.0, 0.0])
ions = [[0]]*64
modes = []
tot = 0
mode = np.array([0.0, 0.0, 0.0])
#for i in range(len(atoms)):
total = 0
for k in range(0, 64):
    modes = []
    dif[k] = dif[k] / norm(dif[k])
    for j in range(192):
        clean = []
        for i in range(64):
            clean.append(np.delete(atoms[i + (j*64)], index))

        clean = np.array(clean)
        x, y, z = clean.T[2], clean.T[3], clean.T[4]
        eigenvector = np.array([x[k], y[k], z[k]]) 
        mode += eigenvector
        
        modes.append(np.dot(eigenvector, dif[0]))
    total += sum(modes)
    print(sum(modes), total)
quit()
# normalising the displacement, and normalising the eigenvalues with respect
# to the displacement yields the same results
# however normalising the eigenvectors such that norm(sum(x), sum(y), sum(z) = 1
# gives a lower and different number

#norm_factor = (norm(mode))
#print(norm(mode))

dif[k] = dif[k] / norm(dif[k])

norm_factor = 1
#print(dif[k])
#norm_factor = norm(dif[k])
#print(norm_factor)


modes = []
for j in range(192):
    clean = []
    for i in range(64):
        clean.append(np.delete(atoms[i + (j*64)], index))

    clean = np.array(clean)
    x, y, z = clean.T[2], clean.T[3], clean.T[4]
    eigenvector = np.array([x[k], y[k], z[k]]) / norm_factor
    
    modes.append(np.dot(eigenvector, dif[k]))

# when you do NO normalising at all, the norm of the dispalcement, is equal to the norm of 
# all of the np.dot(eigenvector, dif[k])
# does that imply that you can represent the displacement fully by a sum of the 
# modal eigenvectors?

#print(len(modes))
#modes = modes/sum(modes) # ???
#print(sum(modes))


# interestingly if you normalising to the sum of the modes, then the previous normalisation of
# the displacement doesn't seem to make a difference anymore
# nor does normalising the eigenvectors with respect to themselves either
#modes = modes/sum(modes)

#print(np.array(modes))
print(sum(modes))
#modes = modes/norm(dif[0])
#print(modes)
#print(sum(modes))
print(sum(wavenumbers.T[1] * modes))
print(sum(wavenumbers.T[1] * 100 * c * modes / 1e12))
#print(prob * sum(wavenumbers.T[1]*100*c * modes) / 1e9)
#print(prob * 100 * 1e12 / 1e9)
#print(np.linalg.norm(tot))
#print(np.dot(x, dif.T[0]))
#print(np.dot(y, dif.T[1]))  
#print(np.dot(z, dif.T[2]))  