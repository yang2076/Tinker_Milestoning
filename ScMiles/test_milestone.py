import sys, os
from scipy.spatial import Voronoi
import numpy as np
AnchorNum = 12
points = np.array([[-70.,90.],[-70.,30.],[-70.,60.],[-70.,0.],[-70.,-30.],[-70.,-70.],[-30.,-70.],[0.,0.],[0.,-70.],[30.,-70.],[60.,-70.],[90.,-70.]])
def check_dihedral_pbc(point):
    point_ = point
    if(point > 180.0):
        point_ -= 360.0
    elif(point < -180.0):
        point_ += 360.0
    return point_

def cross_point(p1, p2, p3, p4):
    x = y = 0
    x1,y1,x2,y2 = p1[0], p1[1], p2[0], p2[1] 
    x3,y3,x4,y4 = p3[0], p3[1], p4[0], p4[1]
    p31 = np.array(p1) - np.array(p3)
    p34 = np.array(p4) - np.array(p3)
    p32 = np.array(p2) - np.array(p3)
    p23 = - p32
    p21 = np.array(p1) - np.array(p2)
    p24 = np.array(p4) - np.array(p2)

    if(max(x1,x2)>=min(x3,x4)
    and max(x3,x4)>=min(x1,x2)
    and max(y1,y2)>=min(y3,y4)
    and max(y3,y4)>=min(y1,y2)
    and np.dot(np.cross(p31,p34),np.cross(p32,p34))<=0
    and np.dot(np.cross(p23,p21),np.cross(p24,p21))<=0):
        if (x2 - x1 == 0.):
            k1 = None
            b1 = 0
        else:
            k1 = (y2 - y1) * 1.0 / (x2 - x1)
            b1 = y1 * 1.0 - x1 * k1 * 1.0
    
        if (x4 - x3) == 0:
            k2 = None
            b2 = 0
        else:
            k2 = (y4 - y3) * 1.0 / (x4 - x3)
            b2 = y3 * 1.0 - x3 * k2 * 1.0
    
        if k1 is None:
            if not k2 is None:
                x = x1
                y = k2 * x1 + b2
        elif k2 is None:
            x = x3
            y = k1 * x3 + b1
        elif not k2 == k1:
            x = (b2 - b1) * 1.0 / (k1 - k2)
            y = k1 * x * 1.0 + b1 * 1.0
    return np.array([x, y])

def Voronoi_tessellation(points):
    # For 2D dihedral voronoi tessellation only
    right_ext = np.hstack((np.ones((AnchorNum,1))*360.0,np.zeros((AnchorNum,1))))
    upper_ext = np.hstack((np.zeros((AnchorNum,1)),np.ones((AnchorNum,1))*360.0))
    points_copy1 = right_ext + points
    points_copy2 = upper_ext + points
    points_copy3 = right_ext + upper_ext + points
    points_final = np.vstack((points, points_copy1, points_copy2, points_copy3))
    possible_MS = []
    possible_MS_property = {}
    vor = Voronoi(points_final)
    vertices = list(vor.vertices)
    ridge_points = list(vor.ridge_points)
    ridge_vertices = list(vor.ridge_vertices)
    for k, i in enumerate(ridge_points):
        pt_pbc_1 = [0., 0.]
        pt_pbc_2 = [0., 0.]
        anchor_1 = min(i[0]%AnchorNum, i[1]%AnchorNum)+1
        anchor_2 = max(i[0]%AnchorNum, i[1]%AnchorNum)+1
        index_1 = ridge_vertices[k][0]
        index_2 = ridge_vertices[k][1]
        if(index_1 == -1 or index_2 == -1):
            continue
        elif((anchor_1, anchor_2) not in possible_MS):
            possible_MS.append((anchor_1,anchor_2))
            ver_1 = vertices[index_1]
            ver_2 = vertices[index_2]
            pt_pbc_1 = cross_point(ver_1, ver_2, [180.,-180.],[180.,180.])
            pt_pbc_2 = cross_point(ver_1, ver_2, [-180.,180.],[180.,180.])
            if(ver_1[0] < ver_2[0]):
                ver_1[0] = check_dihedral_pbc(ver_1[0])
                ver_1[1] = check_dihedral_pbc(ver_1[1])
                ver_2[0] = check_dihedral_pbc(ver_2[0])
                ver_2[1] = check_dihedral_pbc(ver_2[1])
                possible_MS_property[(anchor_1,anchor_2)] = [ver_1, ver_2, pt_pbc_1, pt_pbc_2]
            elif(ver_1[0] > ver_2[0]):
                ver_1[0] = check_dihedral_pbc(ver_1[0])
                ver_1[1] = check_dihedral_pbc(ver_1[1])
                ver_2[0] = check_dihedral_pbc(ver_2[0])
                ver_2[1] = check_dihedral_pbc(ver_2[1])
                possible_MS_property[(anchor_1,anchor_2)] = [ver_2, ver_1, pt_pbc_1, pt_pbc_2]
            elif(ver_1[0] == ver_2[0] and ver_1[1] <= ver_2[1]):
                ver_1[0] = check_dihedral_pbc(ver_1[0])
                ver_1[1] = check_dihedral_pbc(ver_1[1])
                ver_2[0] = check_dihedral_pbc(ver_2[0])
                ver_2[1] = check_dihedral_pbc(ver_2[1])
                possible_MS_property[(anchor_1,anchor_2)] = [ver_1, ver_2, pt_pbc_1, pt_pbc_2]
            elif(ver_1[0] == ver_2[0] and ver_1[1] > ver_2[1]):
                ver_1[0] = check_dihedral_pbc(ver_1[0])
                ver_1[1] = check_dihedral_pbc(ver_1[1])
                ver_2[0] = check_dihedral_pbc(ver_2[0])
                ver_2[1] = check_dihedral_pbc(ver_2[1])
                possible_MS_property[(anchor_1,anchor_2)] = [ver_2, ver_1, pt_pbc_1, pt_pbc_2]
    print(possible_MS)
    return possible_MS_property

t = Voronoi_tessellation(points)
for key in t.keys():
    print(key, t[key])
