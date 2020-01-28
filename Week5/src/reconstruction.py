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
    # TODO take into account more cameras. Check that this code works.
    H = np.linalg.inv(aff_hom)
    Xaff = aff_hom @ Xprj
    cams_aff1 = (cams_pr[-2] @ H )
    cams_aff2 = (cams_pr[-1] @ H )
    cams_aff = [cams_aff1, cams_aff2]
    return Xaff, cams_aff


def resection(tracks, img):
    # your code here
    n = 0
    print("--------------- RESECTION, image", img)
    for i in range(len(tracks)):
        if img in tracks[i].views:
            if abs(tracks[i].pt[3]) > 10**-15:
                n += 1
    points_2d = np.empty((n, 3))
    points_3d = np.empty((n, 4))

    count = 0
    for i in range(len(tracks)):
        if img in tracks[i].views:
            if abs(tracks[i].pt[3]) > 10 ** -15:
                points_2d[count] = np.hstack((tracks[i].views[img], np.array([1])))
                points_3d[count] = tracks[i].pt
                count += 1

    print("---------- CHECKING points 3d --------------")
    for i in range(n):
        print(points_3d[i])

    #################### i) a) Normalization ####################################
    t = - np.mean(points_2d, axis=0)
    print("t", t)
    s = np.sqrt(2) / np.std(points_2d)
    print("s", s)
    T = np.array([[s, 0, t[0]],
                     [0, s, t[1]],
                     [0, 0, 1]])

    t = - np.mean(points_3d, axis=0)
    print("t3d", t)
    s = np.sqrt(3) / np.std(points_3d)
    print("s3d", s)
    U = np.array([[s, 0, 0, t[0]],
                     [0, s, 0, t[1]],
                     [0, 0, s, t[2]],
                     [0, 0, 0, 1]])

    points_2d = T.dot(points_2d.T).T
    points_3d = U.dot(points_3d.T).T

    print("Normalization done")

    ################### i) b) DLT ##################################################

    A = np.empty((2*n, 12))
    for i in range(n):
        x = points_2d[i][0]
        y = points_2d[i][1]
        w = points_2d[i][2]
        X = points_3d[i]
        A[2 * i] = np.array([0, 0, 0, 0,
                             -w * X[0], -w * X[1], -w * X[2], -w * X[3],
                             y * X[0], y * X[1], y * X[2], y * X[3]])
        A[2 * i + 1] = np.array([w * X[0], w * X[1], w * X[2], w * X[3],
                             0, 0, 0, 0,
                             -x * X[0], -x * X[1], -x * X[2], -x * X[3]])

    u, s, v = np.linalg.svd(A)
    p = v.T[:, -1]
    row1 = [p[0], p[1], p[2], p[3]]
    row2 = [p[4], p[5], p[6], p[7]]
    row3 = [p[8], p[9], p[10], p[11]]
    P = np.vstack((row1, row2, row3))

    print("DLT DONE")
    print("--------- P ----------- ")
    print(P)
    ####################### ii) Minimize geometric error ##############################
    iteration = 10 # TODO iterations not working, check derivatives.
    rate = 0.001
    while (iteration < 5):
        error = 0
        dP = np.zeros((3, 4))
        for i in range(n):
            w = points_2d[i][2]
            x = points_2d[i][0] / w
            y = points_2d[i][1] / w
            p3d = points_3d[i]
            #print("oints_3d[i]", points_3d[i])
            p3d /= p3d[3]
            #print("p3d", p3d)
            p_3d_2d = P.dot(points_3d[i])
            #print("p_3d_2d", p_3d_2d)
            w2 = p_3d_2d[2]
            x2 = p_3d_2d[0] / w2
            y2 = p_3d_2d[1] / w2
            #print("x2, y2", x2, y2)

            error += (x-x2)**2 + (y-y2)**2

            x3d = p3d[0]
            y3d = p3d[1]
            z3d = p3d[2]
            w3d = p3d[3]

            dP[0][0] += 2 * (x - x2) * x3d / w2
            dP[0][1] += 2 * (x - x2) * y3d / w2
            dP[0][2] += 2 * (x - x2) * z3d / w2
            dP[0][3] += 2 * (x - x2) * w3d / w2

            dP[1][0] += 2 * (y - y2) * x3d / w2
            dP[1][1] += 2 * (y - y2) * y3d / w2
            dP[1][2] += 2 * (y - y2) * z3d / w2
            dP[1][3] += 2 * (y - y2) * w3d / w2

            dP[2][0] += 2 * (x - x2) * x2*w2 * (-P[2][0])/w2 + 2 * (y - y2) * y2*w2 * (-P[2][0])/w2
            dP[2][1] += 2 * (x - x2) * x2 * w2 * (-P[2][1]) / w2 + 2 * (y - y2) * y2 * w2 * (-P[2][1]) / w2
            dP[2][2] += 2 * (x - x2) * x2 * w2 * (-P[2][2]) / w2 + 2 * (y - y2) * y2 * w2 * (-P[2][2]) / w2
            dP[2][3] += 2 * (x - x2) * x2 * w2 * (-P[2][3]) / w2 + 2 * (y - y2) * y2 * w2 * (-P[2][3]) / w2

        P -= dP*rate
        print("--------- ITERATION ", iteration)
        print(".... Error: ", error)
        print(".... dP: ", dP)
        print(".... P: ", P)

        iteration +=1

        ######################## Denormalization ###########################################
    print ("np.linalg.inv(T)", np.linalg.inv(T))
    print("P", P)
    print("U", U)
    print("shape np.linalg.inv(T)", np.linalg.inv(T).shape)
    print("shape P", P.shape)
    print("shape U", U.shape)
    P = np.linalg.inv(T) @ P @ U
    return P