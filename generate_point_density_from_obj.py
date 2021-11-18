import numpy as np
import re
import pyvista as pv


def ray_triangle_intersection(ray_start, ray_vec, v1,v2,v3):
    """Moeller–Trumbore intersection algorithm.

    Parameters
    ----------
    ray_start : np.ndarray
        Length three numpy array representing start of point.

    ray_vec : np.ndarray
        Direction of the ray.

    triangle : np.ndarray
        ``3 x 3`` numpy array containing the three vertices of a
        triangle.

    Returns
    -------
    bool
        ``True`` when there is an intersection.

    tuple
        Length three tuple containing the distance ``t``, and the
        intersection in unit triangle ``u``, ``v`` coordinates.  When
        there is no intersection, these values will be:
        ``[np.nan, np.nan, np.nan]``

    """
    # define a null intersection
    null_inter = np.array([np.nan, np.nan, np.nan])

    # break down triangle into the individual points
    eps = 0.000001

    # compute edges
    edge1 = v2 - v1
    edge2 = v3 - v1
    pvec = np.cross(ray_vec, edge2)
    det = edge1.dot(pvec)

    if abs(det) < eps:  # no intersection
        return False, null_inter
    inv_det = 1. / det
    tvec = ray_start - v1
    u = tvec.dot(pvec) * inv_det

    if u < 0. or u > 1.:  # if not intersection
        return False, null_inter

    qvec = np.cross(tvec, edge1)
    v = ray_vec.dot(qvec) * inv_det
    if v < 0. or u + v > 1.:  # if not intersection
        return False, null_inter

    t = edge2.dot(qvec) * inv_det
    if t < eps:
        return False, null_inter

    return True, np.array([t, u, v])

reComp = re.compile("(?<=^)(v |vn |vt |f )(.*)(?=$)", re.MULTILINE)
with open('horse.obj') as f:
    data = [txt.group() for txt in reComp.finditer(f.read())]
v_arr, vn_arr, vt_arr, f_arr = [], [], [], []
for line in data:
    tokens = line.split(' ')
    if tokens[0] == 'v':
        v_arr.append([float(c) for c in tokens[1:]])
    elif tokens[0] == 'vn':
        vn_arr.append([float(c) for c in tokens[1:]])
    elif tokens[0] == 'vt':
        vn_arr.append([float(c) for c in tokens[1:]])
    elif tokens[0] == 'f':
        f_arr.append([[int(i) if len(i) else 0 for i in c.split('/')] for c in tokens[1:]])
vertices = []
for face in f_arr:
    #print(face)
    for tp in face:
        vertices += v_arr[tp[0]-1]
#print(vertices)

valid_points = []
while len(valid_points)<100000:
    
    count1 = 0
    count2 = 0
    q0 = [np.random.random(),np.random.random(),np.random.random()]
    w  = [1.0,1.0,1.0]
    for face in f_arr:
        p1 = v_arr[(face[0])[0]-1]
        p2 = v_arr[(face[1])[0]-1]
        p3 = v_arr[(face[2])[0]-1]
        #print(p3)

        inter, tuv = ray_triangle_intersection(np.array(q0),-1*np.array(w),np.array(p1),np.array(p2),np.array(p3))
        if (inter):
            count1 += 1
        inter, tuv = ray_triangle_intersection(np.array(q0),np.array(w),np.array(p1),np.array(p2),np.array(p3))
        if (inter):
            count2 += 1
    if(((count1 % 2) == 1) and ((count2 % 2) == 1)):
        #print("instide")
        #print(q0)
        valid_points += [q0]
        if(len(valid_points)%10 == 0):
            print(len(valid_points))
    #else:
        #print("outside")

txt =""
for s in range(0,len(valid_points)):
    txt += str(valid_points[s][0])+" "+str(valid_points[s][1])+" "+str(valid_points[s][2])+"\n"
file = open("./horse.dat", "w")
file.write(txt)
file.close()