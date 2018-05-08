import numpy as np
import cmath
import math
import time
#import hashlib

print("   ")
file_name="inputs.in"
with open(file_name, 'r') as fn:
    basis_fn = fn.readline().strip()
    psit_fn = fn.readline().strip()
    out_fn = fn.readline().strip()
    Tsteps = int(fn.readline())
    pnum = int(fn.readline())
    plen = int(fn.readline())

start_i = (2*pnum-1)-1
end_i = ((2*pnum-1)+ (2*plen-1))-1
indicies = range(start_i,end_i+1)
####################################
#Import Basis Sets
with open(basis_fn, 'r') as fn:
    basis = np.loadtxt(fn,dtype=int)
basis_Len = basis.shape[0]
nMax = int(basis.shape[1]/2)
####################################

####################################
#Import Psi(t)
psi_t = np.zeros(shape=(Tsteps,basis_Len),dtype=complex)
with open(psit_fn, 'r') as fn:
    for t in range(0, Tsteps):
        for n in range(0, basis_Len):
            re, im = fn.readline().split()
            psi_t[t, n] = complex(float(re), float(im))
####################################

print("  ")
print("=============================================")
print('    Number of rotors:        '+str(nMax))
print('    Number of time steps:    '+str(Tsteps))
print('    Subsystem indicies:      '+str(pnum)+' - '+str(pnum+(plen-1)))
print('    Psi Basis size:          '+str(basis_Len))
print('')
print('   Output filename: '+out_fn)
print("=============================================")
print("  ")

# -------------------------------------------------------------------------------
#                           Compute the sub-basis
new = np.zeros(shape=(basis_Len,plen*2),dtype=int)
not_basis = np.zeros(shape=(basis_Len,2*nMax-(2*plen)),dtype=int)
for n in range(0, basis_Len):
    new[n,:] = np.take(basis[n,:],indicies)
    dummy1 = np.take(basis[n,:], range(0, start_i))
    dummy2 = np.take(basis[n,:], range(end_i + 1, 2 * nMax))
    not_basis[n,:] = np.concatenate((dummy1, dummy2))
sub_basis_long = new
sub_basis = np.unique(new,axis=0)
sub_basis_Len = sub_basis.shape[0]

basis_dict = dict()
for i in range(0,sub_basis_Len):
   basis_dict[tuple(sub_basis[i,:])] = i

# -------------------------------------------------------------------------------

t1=time.time()
# -------------------------------------------------------------------------------
#                     Compute the reduced density matrix

vN = np.zeros(shape=(Tsteps), dtype=float)

for t in range(0,Tsteps):
    vec = psi_t[t,:]
    rho = np.zeros(shape=(sub_basis_Len,sub_basis_Len),dtype=complex)

    tI = time.time()
    #tot1 = 0
    #tot4 = 0
    for i in range(0,basis_Len):

        #t1S = time.time()
        v1p = sub_basis_long[i,:]
        #EYE = np.argmax(np.all(v1p == sub_basis, axis=1))
        EYE = basis_dict[tuple(v1p)]
        #t1F = time.time()
        #tot1 = tot1 + (t1F - t1S)

        v2p = not_basis[i,:]

        for j in range(0,basis_Len):

            v2 = not_basis[j,:]

            if np.all(v2 == v2p):
                #t4S = time.time()
                v1 = sub_basis_long[j,:]
                #JAY = np.argmax(np.all(v1 == sub_basis, axis=1))
                JAY = basis_dict[tuple(v1)]
                rho[EYE,JAY] = rho[EYE,JAY] + vec[i]*np.conjugate(vec[j])
                #t4F = time.time()
                #tot4 = tot4 + (t4F - t4S)


    #for i in range(0, sub_basis_Len):
    #    for j in range(0, i):
    #        rho[i,j] = np.conjugate(rho[j,i])


    #te1 = time.time()
    eigenvals = np.linalg.eigvals(rho)
    #re_eigenvals = np.zeros(shape=(eigenvals.shape[0]),dtype=float)
    tol = 10 ** (-18)
    tot = 0.0
    for n in range(0,eigenvals.shape[0]):
        val = eigenvals[n]
        if val.imag < tol:
            if val.real > tol:
                tot = tot - val.real*math.log(val.real)
            elif val.real < -tol:
                print('WARNING:NEGATIVE EIGEVALUE  ' + str(n))
                print(val)
                #print(val)
                #print('')

        else:
            print('WARNING: COMPLEX EIGEVALUE  '+str(n))

    vN[t] = tot

    #te2 = time.time()
    tF = time.time()
    #print(tot1)
    #print(tot4)
    #print(te2-te1)
    #print('')
    #print(tot1  + tot4)
    if divmod(t, 100)==0:
        print(tF)
# -------------------------------------------------------------------------------
    #print(tot)
t2=time.time()

#print(vN)

print("=============================================")
print('     time = '+str(t2-t1)+' sec.')
print("=============================================")
print("  ")
print("  ")


with open(out_fn, 'w') as fn:
    for s in vN:
       fn.write("%19.16f \n" % s)