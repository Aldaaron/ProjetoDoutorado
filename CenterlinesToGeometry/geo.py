import cubit
from settings import *
import math

def createVertex(x,y,z):
    cubit.cmd("create vertex %f %f %f" % (x,y,z))
    vertexID = cubit.get_last_id("vertex")
    return vertexID

def moveVertex(v1, v2):
    cubit.cmd("move Vertex %d location vertex %f" % (v1,v2))
    
def moveVertexL(v,x,y,z):
    cubit.cmd("move Vertex %d location %f %f %f" % (v,x,y,z))

def createLine(v1,v2):
    cubit.cmd("create curve vertex %d %d" % (v1,v2))
    lineID = cubit.get_last_id("curve")
    return lineID
    
def createSpline(v1,v2,v3):
    cubit.cmd("create curve spline vertex %d %d %d" % (v1,v2,v3))
    splineID = cubit.get_last_id("curve")
    return splineID

def splitCurve(c,v):
    cubit.cmd("split curve %d at vertex %d" % (c,v))
    splineID = cubit.get_last_id("curve")
    return splineID

def createSpline2(pts):
    points = ""
    for p in pts:
        points += str(p) + " "
    cubit.cmd("create curve spline vertex %s" % (points))
    splineID = cubit.get_last_id("curve")
    return splineID

def createEllipse(v1,v2,v3):
    cubit.cmd("create curve vertex %d vertex %d vertex %d  ellipse" % (v1,v2,v3))
    ellipseID = cubit.get_last_id("curve")
    return ellipseID

def createEllipseFull(v1,v2,v3):
    cubit.cmd("create curve vertex %d vertex %d vertex %d  ellipse start angle 0 stop angle 360" % (v1,v2,v3))
    ellipseID = cubit.get_last_id("curve")
    return ellipseID
    
def createVertexOnCurveFraction(c,f):
    cubit.cmd("create vertex on curve %d fraction %f from start" % (c,f))
    vertexID = cubit.get_last_id("vertex")
    return vertexID
    
def createVertexOnCurveDistance(c,d):
    cubit.cmd("create vertex on curve %d distance %f from start" % (c,d))
    vertexID = cubit.get_last_id("vertex")
    return vertexID

def createCircleNormal(v,r,c):
    cubit.cmd("create curve arc center vertex %d radius %f normal curve %d" % (v,r,c))
    circleID = cubit.get_last_id("curve")
    return circleID

def createCircleNormal2(cv,iv,fv,r,c):
    cubit.cmd("create curve arc center vertex %d %d %d radius %f normal curve %d" % (cv,iv,fv,r,c))
    circleID = cubit.get_last_id("curve")
    return circleID
    
def createArc(r,c,n,sa,ea):
    cubit.cmd("create curve arc radius %f center location at vertex %d normal curve %d start angle %f stop angle %f" % (r,c,n,sa,ea))
    arcID = cubit.get_last_id("curve")
    return arcID

def createArcNormal(r,c,n,sa,ea):
    cubit.cmd("create curve arc radius %f center location at vertex %d normal %f %f %f start angle %f stop angle %f" % (r,c,n[0],n[1],n[2],sa,ea))
    arcID = cubit.get_last_id("curve")
    return arcID

def createCombineCurve(curves):
    cs = ""
    for c in curves:
        cs += str(c) + " "
    cubit.cmd("create curve combine curve %s" % (cs))
    curveID = cubit.get_last_id("curve")
    return curveID 
    
def createSweep(s,c):
    cubit.cmd("sweep surface %d along curve %d" % (s,c))
    volID = cubit.get_last_id("volume")
    return volID 

def createSkinCurve2(curves):
    circles = ""
    for c in curves:
        circles += str(c) + " "
    cubit.cmd("create surface skin curve %s" % (circles))
    surfID = cubit.get_last_id("surface")
    return surfID

def createSkinCurve(c1,c2):
    cubit.cmd("create surface skin curve %d %d" % (c1,c2)) 
    surfID = cubit.get_last_id("surface")
    return surfID
    
def createSurfaceCurve(c):
    cubit.cmd("create surface curve %d" % (c))
    surfID = cubit.get_last_id("surface")
    return surfID

def createSurfaceCurve2(curves):
    cc = ""
    for c in curves:
        cc += str(c) + " "
    cubit.cmd("create surface curve %s" % (cc))
    surfID = cubit.get_last_id("surface")
    return surfID
    
def createVolume(surfs):
    volID = cubit.get_owning_volume("surface", surfs[-1])
    ss = ""
    for s in surfs:
        ss += str(s) + " "
    cubit.cmd("create volume surface %s noheal" % (ss))
    #volID = cubit.get_last_id("volume")
    return volID

def createVolumeFromAllSurfaces():
    cubit.cmd("create volume surface all noheal")
    volID = cubit.get_last_id("volume")
    return volID

def mergeAllVertexCurves():
    cubit.cmd("merge vertex all")
    cubit.cmd("merge curve all")
    cubit.cmd("merge surface all")

def imprintMergeAll():
    cubit.cmd("imprint all")
    cubit.cmd("merge all")

def imprintBody(bodies):
    bs = ""
    for b in bodies:
        bs += str(b) + " "
    cubit.cmd("imprint body %s" % (bs))
    
def mergeBody(bodies):
    bs = ""
    for b in bodies:
        bs += str(b) + " "
    cubit.cmd("merge body %s" % (bs))
    
def imprintVolumes(v1,v2):
    cubit.cmd("imprint volume %d %d" % (v1,v2))
    
def mergeVolumes(v1,v2):
    cubit.cmd("merge volume %d %d" % (v1,v2))
    
def unionAll():
    cubit.cmd("unite volume all")

def unionVolumes(v1, v2):
    cubit.cmd("unite volume %d %d"% (v1, v2))
    return v1

def stepUnion(vols):
    volID = vols[0]
    for x in range(0, len(vols)-1):
        cubit.cmd("unite volume %d %d" % (volID,vols[x+1]))
    return volID
    
def createCone(h,r,t):
    cubit.cmd("create frustum height %f radius %f top %f" % (h,r,t))
    volID = cubit.get_last_id('volume')
    return volID
    
def alignVolume(v,s1,s2): 
    cubit.cmd("align Volume %d surface %d with surface %d " % (v,s1,s2))

def colorVolume(v,c):
    cubit.cmd("color volume %d %s" % (v,c))

def colorVolume2(c):
    cubit.cmd("color volume all %s" % (c))

def deleteVolume(v):
    cubit.cmd("delete Volume %d" % (v))

def deleteSurface(s):
    cubit.cmd("delete Surface %d" % (s))

def deleteCurve(c):
    cubit.cmd("delete Curve %d" % (c))

def deleteCurve2(c):
    cubit.cmd("delete Curve %d keep" % (c))

def deleteVertex(v):
    cubit.cmd("delete Vertex %d" % (v))

def spheCoords(coords):
    radius = math.sqrt(coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2]);
    theta  = math.acos(coords[2] / radius);
    phi    = math.atan2(coords[1], coords[0]);
    coords[0] = radius
    coords[1] = theta
    coords[2] = phi
    
def cartCoords(coords):
    x = coords[0] * math.sin(coords[1]) * math.cos(coords[2]);
    y = coords[0] * math.sin(coords[1]) * math.sin(coords[2]);
    z = coords[0] * math.cos(coords[1]);
    coords[0] = x
    coords[1] = y
    coords[2] = z
    
def degreeToRadians(degree):
    return ((degree*math.pi)/180)

def radiansToDegree(radians):
    return ((radians*180)/math.pi)
    
def distance(p1, p2):
    return math.sqrt(math.pow(p1[0]-p2[0],2) + math.pow(p1[1]-p2[1],2) + math.pow(p1[2]-p2[2],2))

def getAngles(p1, p2):
    deltaX = p2[0]-p1[0];
    deltaY = p2[1]-p1[1];
    deltaZ = p2[2]-p1[2];
    angles = [math.atan2(deltaY, deltaX), math.atan2(deltaZ, deltaX)]
    return angles

def rotate(point, origin, angle):
    deltaY = point[1] - origin[1];
    deltaX = point[0] - origin[0];
    deltaZ = point[2] - origin[2];
    ay = math.atan2(deltaY, deltaX)
    az = math.atan2(deltaZ, deltaX)
    dist = distance(point,origin)
    coords = [0,0,0]
    angleR = degreeToRadians(angle)
    coords[0] += dist
    coords[1] += degreeToRadians(90)
    coords[2] += ay+angleR
    cartCoords(coords)
    point[0] = coords[0] + origin[0]
    point[1] = coords[1] + origin[1]
    point[2] = coords[2] + origin[2]

def rotateAndMove(point, origin, angleY, angleZ, dist):
    deltaY = point[1] - origin[1];
    deltaX = point[0] - origin[0];
    deltaZ = point[2] - origin[2];
    ay = math.atan2(deltaY, deltaX)
    az = math.atan2(deltaZ, deltaX)
    dist += distance(point,origin)
    coords = [0,0,0]
    angleYR = degreeToRadians(angleY)
    angleZR = degreeToRadians(angleZ)
    coords[0] += dist
    coords[1] += az+angleZR+degreeToRadians(90)
    coords[2] += ay+angleYR
    cartCoords(coords)
    point[0] = coords[0] + origin[0]
    point[1] = coords[1] + origin[1]
    point[2] = coords[2] + origin[2]
    
def rotate2(point, origin, angle, dist):
    deltaY = point[1] - origin[1];
    deltaX = point[0] - origin[0];
    deltaZ = point[2] - origin[2];
    ay = math.atan2(deltaY, deltaX)
    ay = ay+degreeToRadians(angle)
    az = math.atan2(deltaZ, deltaX)
    az = az+degreeToRadians(90)
    if deltaX < 0:
        ay = ay
        az = az+degreeToRadians(-180)
    coords = [0,0,0]
    coords[0] += dist
    coords[1] += az
    coords[2] += ay
    cartCoords(coords)
    point[0] = coords[0] + origin[0]
    point[1] = coords[1] + origin[1]
    point[2] = coords[2] + origin[2]
        