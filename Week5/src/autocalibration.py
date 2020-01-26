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

