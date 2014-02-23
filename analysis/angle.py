import math
import MD.base.points as points


#project vector from center to end ontk z vector
def face_check(VW,Z,k,A,B,side1,side2,L,kappa=1/10.):
    angle_cut=kappa*math.pi
    v1 = points.vector1d(VW[k][A],Z[k][A*6+side1],L[k])
    v2 = points.vector1d(VW[k][B],Z[k][B*6+side2],L[k])
    angle = points.angle_between(v1,v2)
    if ((angle%(math.pi/2) < angle_cut) or
        (angle%(math.pi/2) > math.pi/2- angle_cut)):
        return 1
    return 0
