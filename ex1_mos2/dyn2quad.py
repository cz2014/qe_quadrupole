import numpy as np
from numpy.core.numeric import tensordot 
from scipy.optimize import leastsq
import matplotlib.pyplot as plt

def reader_info(FILENAME='./info', ndim=2):
    fo = open(FILENAME, 'r')
    line = fo.readline()
    while line:
        line = line.split()
        if not line:
            line = fo.readline()
            continue

        if line[1] == 'symmetry':
            nsym = int(line[0])
            sym_matrix = np.zeros((nsym, 3, 3))
            for i in range(nsym//6):
                line = fo.readline()
                line = fo.readline().split()
                sym_matrix[i*6:(i+1)*6, 0] = np.reshape([int(s) for s in line], (6,3))
                line = fo.readline().split()
                sym_matrix[i*6:(i+1)*6, 1] = np.reshape([int(s) for s in line], (6,3))
                line = fo.readline().split()
                sym_matrix[i*6:(i+1)*6, 2] = np.reshape([int(s) for s in line], (6,3))

            if nsym%6 != 0:
                tmp = nsym%6
                line = fo.readline()
                line = fo.readline().split()
                sym_matrix[i*6:i*6+tmp, 0] = np.reshape([int(s) for s in line], (tmp,3))
                line = fo.readline().split()
                sym_matrix[i*6:i*6+tmp, 1] = np.reshape([int(s) for s in line], (tmp,3))
                line = fo.readline().split()
                sym_matrix[i*6:i*6+tmp, 2] = np.reshape([int(s) for s in line], (tmp,3))

        line = fo.readline()

    fo.close()

    if ndim == 2:
        sym_matrix = sym_matrix[:,0:2,0:2]

    return sym_matrix 

def sym_f2c(rec_vec, sym_f):
    rinv = np.linalg.inv(rec_vec)
    nsym = len(sym_f)
    sym_c = np.zeros((nsym, 3, 3))

    for i in range(nsym):
        sym_c[i] = np.dot(np.dot(rinv, sym_f[i].T), rec_vec)

    return sym_c 

def trunc_sym(spq, sym):
    small_sym = []
    for i in range(len(sym)):
        tmpq = np.dot(sym[i],spq)
        if np.linalg.norm(tmpq-spq)<1e-5:
            small_sym.append(np.copy(sym[i]))

    return np.array(small_sym)

def symtensor(ts, sym):
    sumts = np.zeros(np.shape(ts), dtype=ts.dtype) 

    nsym = len(sym)
    for i in range(nsym):
        sumts += np.dot(np.dot(np.linalg.inv(sym[i]), ts), sym[i])

    sumts = sumts / nsym

    return sumts 

def reader_zue(file="./dyn"):
    fo = open(file, 'r')

    s2c = lambda line: [complex(float(line.split()[0]), float(line.split()[1])), \
                        complex(float(line.split()[2]), float(line.split()[3])), \
                        complex(float(line.split()[4]), float(line.split()[5]))]

    lines = fo.readlines()
    natm = int(lines[2].split()[1])
    zue = np.zeros((natm,3,3), dtype=complex)
    for i in range(len(lines)):
        if "Effective Charges U-E: Z_{s,alpha}{beta}" in lines[i]:
            i += 2
            for iatm in range(natm):
                i += 1
                for ipol in range(3):
                    zue[iatm,ipol,:] = s2c(lines[i])
                    i += 1

            break 
        else:
            i += 1

    return zue 

def reader_dynheader(file):
    def s2a(type, s):
        return [type(et) for et in s.split()]

    with open(file, 'r') as fo:
        l = fo.readlines()
    
    ntyp, nat, ibrav = s2a(int, l[2][:14])
    if ibrav == 4:
        celldm = s2a(float, l[2][14:])
        alat2bohrz = celldm[0]#*celldm[2]
        Lz = celldm[0]*celldm[2]
    else:
        raise Exception("not implemented for this ibrav")

    ntyp = int(l[2].split()[0])
    nat = int(l[2].split()[1])
    pos = np.zeros((nat, 3))
    for iat in range(nat):
        ind0 = 3+ntyp
        pos[iat] = s2a(float, l[ind0+iat])[2:5]

    posz = pos[:,2] * alat2bohrz
    posz = posz - np.sum(posz)/nat

    return nat, posz, Lz


def reader_epsilon(file):
    fo = open(file, 'r')

    s2f = lambda line: [float(line.split()[i]) for i in [0,1,2]] 

    lines = fo.readlines()
    epsilon = np.zeros((3,3))
    for i in range(len(lines)):
        if "Dielectric Tensor" in lines[i]:
            i += 2
            for ipol in range(3):
                epsilon[ipol,:] = s2f(lines[i])
                i += 1
            break 
        else:
            i += 1

    return epsilon 

def reader_q(file="./dyn"):
    fo = open(file, 'r')

    lines = fo.readlines()
    alat = float(lines[2].split()[3])
    ntype = int(lines[2].split()[0])
    natm = int(lines[2].split()[1])

    pos_iline = [3+ntype+i for i in range(natm)] 
    pos_cart = np.array([ [ float(xs) for xs in lines[il].split()[2:5] ] for il in pos_iline])

    for i in range(len(lines)):
        if "Dynamical  Matrix in cartesian axes" in lines[i]:
            i += 2
            line = lines[i].split()
            q = [float(line[ichunk]) for ichunk in [3,4,5]]
            return np.array(q), alat, pos_cart, natm 
    else: 
        i += 1
    
    Exception("No q points found!")
    return 0 

def zue2Q(izuelist, pattern, alat):
    """ izuelist: list of izue;
        pattern: (index of zue2, index of zue1, q length)"""
    
    Qtmp = np.zeros((3,3,3), dtype=complex) # [gamma, alpha, beta]
    for ipol, pa in enumerate(pattern):
        # Qtmp[ipol] = np.real( 1j*(izuelist[pa[0]] - izuelist[pa[1]])/(pa[2]*2*np.pi/alat) )
        Qtmp[ipol] = 1j*(izuelist[pa[0]] - izuelist[pa[1]])/(pa[2]*2*np.pi/alat)
        # print('\n'.join([''.join(['{:50}'.format(item) for item in row]) for row in Qtmp[ipol]])) 


    Q = np.zeros((3,3,3), dtype=complex) # beta, alpha, gamma
    for ipol in [0,1,2]:
        Q[ipol] = Qtmp[:,:,ipol] #+ Qtmp[:,:,ipol].T  

    return Q 

def reader_zue_series(prefix, nfile):
    rho_all = []
    q_all = []
    for ifile in range(nfile):
        dynfile = prefix+str(ifile+1)
        qpoint, alat, _, natm = reader_q(dynfile)
        
        tmpzue = reader_zue(dynfile)
        rho_all.append( np.copy(tmpzue[:,:,0]) )
        q_all.append( np.copy(qpoint) ) 
    rho_all = np.array(rho_all) # [iq, iatm, 3]
    q_all = np.array(q_all) 

    fo = open(prefix+"_sum", 'w')
    for iatm in range(natm):
        fo.write("atom no. "+str(iatm)+"\n")
        for ifile in range(nfile):
            for x in q_all[ifile,:]: fo.write("%6.2f "% (x))
            fo.write("   ")
            for x in rho_all[ifile,iatm,:]: fo.write("%+6.3e %+6.3e  "% (np.real(x), np.imag(x)))
            fo.write("\n")
        fo.write("\n\n")
    
    fo.close()
    return q_all*2*np.pi/alat, rho_all # [ifile, ipol], [ifile,iatm,ipol]

def corr_rho(q_all, rho_all, epsilon, cmode="3d", a=1):
    for ifile in range(len(q_all)):
        iq = q_all[ifile]
        if np.linalg.norm(iq) > 1e-6:
            if cmode in ["3d", "2d"]:
                corr = np.dot(iq, np.dot(epsilon, iq)) / np.linalg.norm(iq)**2 
            elif cmode == "2d_d": # wrong: not used 
                # corr = 1 + 70*np.linalg.norm(iq) #as
                corr = 1 + 51*np.linalg.norm(iq) #mos2
                print(iq, corr)
            rho_all[ifile] = rho_all[ifile] * corr 

    return 


def list2symarray(l):
    return np.array([[l[0],l[5],l[4]],[l[5],l[1],l[3]],[l[4],l[3],l[2]]])

def fit_Z(q_all, rho_all):
    rho_Z = np.imag(rho_all)
    err_for_z = lambda zue, xq, rho: [-1*np.dot(zue, iq) - irho for iq, irho in zip(xq, rho)]

    natm = len(rho_Z[0]) 
    """ fit zue """
    zue_fitted = np.zeros((natm,3,3)) 
    zue_err = []
    for iatm in range(natm):
        for ipol in range(3):
            res_z = leastsq(err_for_z, [0,0,0], args=(q_all, rho_Z[:,iatm,ipol])) 
            zue_fitted[iatm,ipol,:] = res_z[0]
            zue_err.extend( [err, val] for err, val in \
                zip(err_for_z(res_z[0], q_all, rho_Z[:,iatm,ipol]), rho_Z[:,iatm,ipol]) )
    zue_err = np.array(zue_err) 

    return zue_fitted, zue_err

def fit_Q(q_all, rho_all):
    rho_Q = np.real(rho_all)
    err_for_q = lambda quad, xq, rho: [-0.5*np.dot(iq, np.dot(list2symarray(quad), iq)) - irho \
        for iq, irho in zip(xq, rho)] 

    natm = len(rho_Q[0]) 
    """ fit quad """
    quad_fitted = np.zeros((natm,3,6)) 
    quad_err = []
    for iatm in range(natm):
        for ipol in range(3):
            res_q = leastsq(err_for_q, [0,0,0,0,0,0],args=(q_all, rho_Q[:,iatm,ipol]))
            quad_fitted[iatm,ipol,:] = res_q[0]
            quad_err.extend( [err, val] for err, val in \
                zip(err_for_q(res_q[0], q_all, rho_Q[:,iatm,ipol]), rho_Q[:,iatm,ipol]) ) 
    quad_err = np.array(quad_err)

    return quad_fitted, quad_err 

def rot(tensor, rm, axes):
    perm = np.arange(len(np.shape(tensor)))
    perm[axes[0]] = perm[-1]
    perm[-1] = axes[0]

    res = np.tensordot(tensor, rm, axes=axes)
    res = np.transpose(res, axes=perm)

    return res


def corr_zue(zue0, bcc, epsil=None, rotz=0):
    zue = np.copy(zue0)

    if bcc in ['2d']:
        ezz = epsil[2,2]
        zue[:,2,:] = zue[:,2,:]/ezz
        zue[:,:2,2] = zue[:,:2,2]/ezz

    th = rotz/180*np.pi
    rm = np.array([[np.cos(th), -np.sin(th), 0], 
                [np.sin(th), np.cos(th), 0], [0,0,1]])
    
    zue = rot(zue, rm, axes=(1,0))
    zue = rot(zue, rm, axes=(2,0))

    return zue

def corr_quad(quad0, zue0, bcc, posz, epsil=None, rotz=0, sym='r3'):

    def rotm(theta):
        th = theta/180*np.pi
        return np.array([[np.cos(th), -np.sin(th), 0], 
                [np.sin(th), np.cos(th), 0], [0,0,1]])

    def symarray2list(symm):
        return np.array([symm[0,0], symm[1,1], symm[2,2], symm[1,2], symm[0,2], symm[0,1]])

    def cz_tdot(syml, rm):
        symmat = list2symarray(syml)
        symmat = rot(symmat, rm, axes=(0,0))
        symmat = rot(symmat, rm, axes=(1,0))
        return symarray2list(symmat)

    def rotquad(quad, rm):
        quad = rot(quad, rm, axes=(1,0))
        for iat in range(len(quad)):
            for b in [0,1,2]:
                quad[iat, b, :] = cz_tdot(quad[iat, b, :], rm)
        return quad

    quad = np.copy(quad0)

    if bcc in ['2d']:
        ezz = epsil[2,2]
        chi = (epsil[:2,:2] - np.array([[1,0],[0,1]]))
        for iat in range(len(quad)):
            for b in [0,1,2]:
                quad[iat,b,2] = (quad[iat,b,2] + 2*posz[iat]*zue0[iat,b,2])/ezz
                quad[iat,b,3] = (quad[iat,b,3] + posz[iat]*zue0[iat,b,1])/ezz
                quad[iat,b,4] = (quad[iat,b,4] + posz[iat]*zue0[iat,b,0])/ezz

                qzzb = quad[iat,b,2]
                quad[iat,b,0] = quad[iat,b,0] - chi[0,0]*qzzb
                quad[iat,b,1] = quad[iat,b,1] - chi[1,1]*qzzb
                quad[iat,b,5] = quad[iat,b,5] - chi[0,1]*qzzb


    rm = rotm(rotz)
    quad = rotquad(quad, rm)

    if sym in ['r3']:
        print("=======================================================")
        print("Caution: 3-fold symmetry assumed!")
        rm1 = rotm(120)
        rm2 = rotm(-120)
        quadunr = np.copy(quad)
        quad = quadunr + rotquad(quadunr, rm1) + rotquad(quadunr, rm2)
        quad = quad/3.0

    return quad

def writer_forepw(dir, quad, epszz):
    fo = open(dir+'quadrupole.fmt', 'w')
    fo.write("atom   dir       Qxx         Qyy         Qzz         Qyz         Qxz         Qxy\n")
    for iatm in range(len(quad)):
        for idir in [0,1,2]:
            fo.write( "%6i%6i" %(iatm+1, idir+1))
            fo.write( "%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f\n" % tuple(quad[iatm,idir]))
    fo.write("%16.12f" % (epszz))
    fo.close()


if __name__ == "__main__":
    interact = True

    ############# now test a series of q 
    dirprefix = "./ds/r17_inse/"
    dynprefix = dirprefix+'dyn'
    nfile = 16
    
    if interact:
        st = input('Enter the folder name of dyn files:\n').strip()
        if st[-1] != '/': st = st + '/'
        dirprefix = st
        dynprefix = dirprefix+'dyn'
        st = input('Enter the number of dyn files:\n')
        nfile = int(st)


    bcc = "2d"
    # bcc = False 
    # boundary condition correction
    natm, posz, Lz = reader_dynheader(dynprefix+'0')

    q_all, rho_all = reader_zue_series(dynprefix, nfile)
    if bcc:
        epsilon = reader_epsilon(dynprefix+'0')
        # print(np.linalg.norm(epsilon)) 
        # print("eps matrix:")
        # print(epsilon) 
        corr_rho(q_all, rho_all, epsilon, bcc) 
        # print("atomic posz:")
        # print(posz)
        print("=======================================================")
        print("L_min:", Lz*(1-1/epsilon[2,2]))
        print("Caution: L_min should smaller than 15!")


    # fit zue
    zue0, err = fit_Z(q_all, rho_all)
    zue = corr_zue(zue0, bcc, epsilon, rotz=0)
    fo = open(dynprefix+"_zue", 'w')
    for iatm in range(natm):
        fo.write("atom no. "+str(iatm)+"\n")
        for ipol in range(3):
            for x in zue[iatm,ipol,:]: fo.write("%8.4f "% (x))
            fo.write("\n")
        fo.write("\n\n")
    fo.close()

    plt.plot(err[:,1]+err[:,0], marker='o', label='raw')
    plt.plot(err[:,1], marker='o', label='fitted')
    plt.legend()
    # plt.show() 
    plt.savefig(dirprefix+'fit_zue.png')
    plt.close()

    # fit quad
    quad0, err2 = fit_Q(q_all, rho_all)
    quad = corr_quad(quad0, zue0, bcc, posz, epsilon, rotz=0, sym='r3')
    fo = open(dynprefix+"_quad", 'w')
    for iatm in range(natm):
        fo.write("atom no. "+str(iatm)+"\n")
        for ipol in range(3):
            for x in quad[iatm,ipol,:]: fo.write("%8.4f "% (x))
            fo.write("\n")
        fo.write("\n\n")
    fo.close()

    writer_forepw(dirprefix, quad, epsilon[2,2])


    plt.plot(err2[:,1]+err2[:,0], marker='o', label='raw')
    plt.plot(err2[:,1], marker='o', label='fitted')
    plt.legend()
    # plt.show() 
    plt.savefig(dirprefix+'fit_quad.png')
    plt.close()
