import numpy as np

import utils as h
import reconstruction as rc
import maths as mth

def estimate_aff_hom(cams, vps):
    # your code here
    C1 = cams[-2]
    C2 = cams[-1]

    vps1 = vps[-2]
    vps2 = vps[-1]

    es = rc.estimate_3d_points(C1, C2, vps1.T, vps2.T)
    es = es/(es[3])
    A = [es.T]
    U, S, V = np.linalg.svd(A)

    Ha = V.T[:, 3];
    Ha = Ha/Ha[3]

    aff_hom = np.identity(4)
    aff_hom[3,:] = Ha.T

    return aff_hom 

def estimate_euc_hom(cams_aff, vps):

    v1 = vps[0]
    v2 = vps[1]
    v3 = vps[2]

    row1 = [v1[0]*v2[0], v1[0]*v2[1]+v1[1]*v2[0], v1[0]*1+1*v2[0], v1[1]*v2[1], v1[1]*1+1*v2[1], 1*1]
    row2 = [v1[0]*v3[0], v1[0]*v3[1]+v1[1]*v3[0], v1[0]*1+1*v3[0], v1[1]*v3[1], v1[1]*1+1*v3[1], 1*1]
    row3 = [v2[0]*v3[0], v2[0]*v3[1]+v2[1]*v3[0], v2[0]*1+1*v3[0], v2[1]*v3[1], v2[1]*1+1*v3[1], 1*1]
    row4 = [0,1,0,0,0,0]
    row5 = [1,0,0,-1,0,0]

    A = np.vstack((row1,row2,row3,row4,row5))

    U, S, V = np.linalg.svd(A)

    W = V.T[:,-1]


    row1 = [W[0], W[1], W[2]]
    row2 = [W[1], W[3], W[4]]
    row3 = [W[2], W[4], W[5]]
    M = cams_aff[:, 0:2]
    Wp  =  np.vstack((row1,row2,row3))
    r = M.T @ Wp @ M
    ri = np.linalg.inv(r)
    A = np.linalg.cholesky(ri)

    Ha = np.eye(4, 4)
    Ha[0: 2, 0: 2] = np.linalg.inv(A)
    return Ha