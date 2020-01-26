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

def compute_reproj_error(X, cam): 
    # your code here 

    return error

def transform(aff_hom, Xprj, cams_pr):
    # your code here

    return Xaff, cams_aff

def resection(tracks, img):
    # your code here

    return P
