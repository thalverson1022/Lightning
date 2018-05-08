import numpy as np
import time

print("   ")
file_name="inputs.in"
with open(file_name, 'r') as fn:
    H0_basis_fn = fn.readline().strip()
    Initial_state_fn = fn.readline().strip()
    H1_basis_fn = fn.readline().strip()
    H1_eigen_fn = fn.readline().strip()
    out_fn = fn.readline().strip()

####################################
#Import Basis Sets
with open(H0_basis_fn, 'r') as fn:
    basis0 = np.loadtxt(fn,dtype=int)
basis0_Len = basis0.shape[0]

with open(H1_basis_fn, 'r') as fn:
    basis1 = np.loadtxt(fn,dtype=int)
basis1_Len = basis1.shape[0]
####################################


####################################
#Import State vector
psi0 = np.zeros(shape=(basis0_Len),dtype=complex)
with open(Initial_state_fn, 'r') as fn:
    for n in range(0,basis0_Len):
        re, im = fn.readline().split()
        psi0[n] = complex(float(re), float(im))
####################################


####################################
#Import H(t) info
with open(H1_eigen_fn, 'r') as fn:
    nMax = int(fn.readline())
    LMax = int(fn.readline())
    gamma = float(fn.readline())
    delta_t = float(fn.readline())
    Tsteps = int(fn.readline())
    vals_path = fn.readline().strip()
    vecs_path = fn.readline().strip()
####################################

eigenvals = np.zeros(shape=(Tsteps,basis1_Len),dtype=float)
for t in range(0,Tsteps):
    dummy_fn=vals_path+'out-t'+str(t+1)+'.dat'
    with open(dummy_fn, 'r') as fn:
        for j, line in enumerate(fn):
            eigenvals[t,j] = float(line)

print("  ")
print("=============================================")
print('    Number of rotors:        '+str(nMax))
print('    L parameter:             '+str(LMax))
print('    Gamma:                   '+str(gamma))
print('    Time Step:               '+str(delta_t))
print('    Number of time steps:    '+str(Tsteps))
print('    Psi0 Basis size:         '+str(basis0_Len))
print('    Psi(t) Basis size:       '+str(basis1_Len))
print('')
print('   Output filename: '+out_fn)
print("=============================================")
print("  ")


t1=time.time()
# -------------------------------------------------------------------------------
#                              Compute Psi(t)
psit = np.zeros(shape=(Tsteps,basis1_Len),dtype=complex)
psit[0] = psi0
for t in range(1,Tsteps):

    dummy_fn = vecs_path + 'out-t' + str(t + 1) + '.dat'
    #t_f1 = time.time()
    VECS = np.zeros(shape=(basis1_Len,basis1_Len),dtype=complex)
    with open(dummy_fn, 'r') as fn:
        for i in range(0,basis1_Len):
            for j in range(0, basis1_Len):
                re,im = fn.readline().split()
                VECS[i,j]=complex(float(re),float(im))
                #print[VECS[i,j]]
    #t_f2 = time.time()

    #t_r1 = time.time()
    rhot = np.zeros(shape=(basis1_Len, basis1_Len, basis1_Len), dtype=complex)
    for n in range(0, basis1_Len):
        vec = VECS[n, :]
        rhot[n]=np.outer(vec,np.conjugate(vec))
    #t_r2 = time.time()

    #t_s1 = time.time()
    phat = np.zeros(shape=(basis1_Len, basis1_Len), dtype=complex)
    for n in range(0, basis1_Len):
        eps = eigenvals[t,n]#-eigenvals[0,0]
        phat = phat + np.exp(-complex(0.0,1.0)*eps*delta_t)*rhot[n,:,:]
    #t_s2 = time.time()

    #t_m1 = time.time()
    psit[t] = np.matmul(phat,psit[t-1])
    #t_m2 = time.time()
    # print(t_f2-t_f1)
    # print(t_r2-t_r1)
    # print(t_s2 - t_s1)
    # print(t_m2 - t_m1)
    # print('')
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

t2=time.time()

print("=============================================")
print('     time = '+str(t2-t1)+' sec.')
print("=============================================")
print("  ")
print("  ")



# np.set_printoptions(precision=5)
# np.set_printoptions(suppress=True)
# print(psit)

# #Write output to file
with open(out_fn, 'w') as fn:
    for t in range(Tsteps):
        for n in range(basis1_Len):
            re = psit[t,n].real
            im = psit[t,n].imag
            fn.write("%19.16f    %19.16f \n" % (re,im))
           # fn.write(str(re) + "  " + str(im) + '\n')



