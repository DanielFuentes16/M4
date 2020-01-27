import cv2
import numpy as np

import utils as h
import maths as mth
import scipy
from scipy import linalg, matrix
from scipy.stats import kurtosis, skew


def compute_proj_camera(F, i):
    # Result 9.15 of MVG (v = 0, lambda = 1). It assumes P1 = [I|0]
    # your code here

    P = np.array([[1.,0.,0.,0.], [0.,1.,0.,0.], [0.,0.,1.,0.]])
    u, s, vh = scipy.linalg.svd(F.conj().T)
    e12 = np.array(vh.conj().T[:, 2])
    skw = np.array([[0.,-e12[2],e12[1]], [e12[2],0.,-e12[0]], [-e12[1],e12[0],0.]])
    M1= skw @ F


    Pprime = np.hstack((M1,e12.reshape(3,1)))
    #print(Pprime)
    return Pprime

def estimate_3d_points(P1, P2, xr1, xr2):
    # Triangulate 3D points from camera matrices
    Xprj = cv2.triangulatePoints(P1, P2, xr1, xr2) 

    # Divide by the last column 
    Xprj = Xprj / Xprj[3, :]

    if h.debug >2:
        print("  X estimated:\n", Xprj)

    return Xprj

def compute_reproj_error(Xx, cam, x1, x2):
    P1 =  cam[-2]
    P2 =  cam[-1]
    x1p = P1 @ Xx
    x2p = P2 @ Xx



    x1p = x1p/x1p[2,:]
    x2p = x2p/x2p[2,:]

    x1p = np.delete(x1p.T,2, axis = 1)
    x2p = np.delete(x2p.T, 2, axis=1)

    d1 = (x1p[:,0]- x1[:,0])**2 + (x1p[:,1]- x1[:,1])**2
    d2 = (x2p[:,0]- x2[:,0])**2 + (x2p[:,1]- x2[:,1])**2
    error = d1 + d2
    #d = np.linalg.norm(x1p[:2,:]-x1.T, axis = 0) + np.linalg.norm(x2p[:2,:]-x2.T, axis = 0)
    return np.sum(error)

def transform(aff_hom, Xprj, cams_pr):
    # your code here

    H = np.linalg.inv(aff_hom)
    Xaff = aff_hom @ Xprj

    cams_aff1 = (cams_pr[0] @ aff_hom )
    cams_aff2 = (cams_pr[1] @ H )
    cams_aff = [cams_aff1, cams_aff2]
    return Xaff, cams_aff


def resection(tracks, img):
    # your code here

    return P
