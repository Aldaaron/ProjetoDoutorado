import cubit
from settings import *
import math

def createVertex(x,y,z):
    printM("create vertex %f %f %f" % (x,y,z))
    cubit.cmd("create vertex %f %f %f \n" % (x,y,z))
    vertexID = cubit.get_last_id("vertex")
    printM("new Vertex " + str(vertexID) + " created")
    return vertexID

def moveVertex(v1, v2):
    printM("move Vertex %d location vertex %f" % (v1,v2))
    cubit.cmd("move Vertex %d location vertex %f  \n" % (v1,v2))
    
def moveVertexL(v,x,y,z):
    printM("move Vertex %d location %f %f %f" % (v,x,y,z))
    cubit.cmd("move Vertex %d location %f %f %f" % (v,x,y,z))

def createLine(v1,v2):
    printM("create curve vertex %d %d" % (v1,v2))
    cubit.cmd("create curve vertex %d %d" % (v1,v2))
    lineID = cubit.get_last_id("curve")
    printM("new Curve " + str(lineID) + " created")
    return lineID
    
def createSpline(v1,v2,v3):
    printM("create curve spline vertex %d %d %d " % (v1,v2,v3))
    cubit.cmd("create curve spline vertex %d %d %d" % (v1,v2,v3))
    splineID = cubit.get_last_id("curve")
    printM("new Curve " + str(splineID) + " created")
    return splineID

def createSpline2(pts):
    points = ""
    for p in pts:
        points += str(p) + " "
    printM("create curve spline vertex %s " % (points))
    cubit.cmd("create curve spline vertex %s " % (points))
    splineID = cubit.get_last_id("curve")
    printM("new Curve " + str(splineID) + " created")
    return splineID

def createEllipse(v1,v2,v3):
    printM("create curve vertex %d vertex %d vertex %d  ellipse" % (v1,v2,v3))
    cubit.cmd("create curve vertex %d vertex %d vertex %d  ellipse" % (v1,v2,v3))
    ellipseID = cubit.get_last_id("curve")
    printM("new Curve " + str(ellipseID) + " created")
    return ellipseID

def createEllipse2(v1,v2,v3):
    printM("create curve vertex %d %d %d ellipse " % (v1,v2,v3))
    cubit.cmd("create surface ellipse vertex %d %d %d " % (v1,v2,v3))
    ellipseS = cubit.get_last_id("surface")
    ellipseID = cubit.get_last_id("curve")
    cubit.cmd("curve %d copy " % (ellipseID))
    ellipseID = cubit.get_last_id("curve")
    deleteSurface(ellipseS)
    deleteVertex(v2)
    printM("new Curve " + str(ellipseID) + " created")
    return ellipseID
    
def createVertexOnCurveFraction(c,f):
    printM("create vertex on curve %d fraction %f from start" % (c,f))
    cubit.cmd("create vertex on curve %d fraction %f from start \n" % (c,f))
    vertexID = cubit.get_last_id("vertex")
    printM("new Vertex " + str(vertexID) + " created")
    return vertexID
    
def createVertexOnCurveDistance(c,d):
    printM("create vertex on curve %d distance %f from start" % (c,d))
    cubit.cmd("create vertex on curve %d distance %f from start \n" % (c,d))
    vertexID = cubit.get_last_id("vertex")
    printM("new Vertex " + str(vertexID) + " created")
    return vertexID

def createCircleNormal(v,r,c):
    printM("create curve arc center vertex %d radius %f normal curve %d" % (v,r,c))
    cubit.cmd("create curve arc center vertex %d radius %f normal curve %d \n" % (v,r,c))
    circleID = cubit.get_last_id("curve")
    printM("new Curve " + str(circleID) + " created")
    return circleID

def createCircleNormal2(cv,iv,fv,r,c):
    printM("create curve arc center vertex %d %d %d radius %f normal curve %d" % (cv,iv,fv,r,c))
    cubit.cmd("create curve arc center vertex %d %d %d radius %f normal curve %d \n" % (cv,iv,fv,r,c))
    circleID = cubit.get_last_id("curve")
    printM("new Curve " + str(circleID) + " created")
    return circleID
    
def createArc(r,c,n,sa,ea):
    printM("create curve arc radius %f center location at vertex %d normal curve %d start angle %f stop angle %f" % (r,c,n,sa,ea))
    cubit.cmd("create curve arc radius %f center location at vertex %d normal curve %d start angle %f stop angle %f" % (r,c,n,sa,ea))
    arcID = cubit.get_last_id("curve")
    printM("new Curve " + str(arcID) + " created")
    return arcID
    
def createSweep(s,c):
    printM("sweep surface %d along curve %d" % (s,c))
    cubit.set_modified()
    cubit.cmd("sweep surface %d along curve %d" % (s,c))
    if not cubit.is_modified():
        print "erro sweep surface %d" % (c)
    volID = cubit.get_last_id("volume")
    printM("new Volume " + str(volID) + " created")
    return volID 

def createSkinCurve2(curves):
    circles = ""
    for c in curves:
        circles += str(c) + " "
    printM("create surface skin curve %s" % (circles))
    cubit.set_modified()
    cubit.cmd("create surface skin curve %s\n" % (circles))
    if not cubit.is_modified():
        print "Erro curvas: " + circles
    surfID = cubit.get_last_id("surface")
    printM("new Surface " + str(surfID) + " created")
    return surfID

def createSkinCurve(c1,c2):
    printM("create surface skin curve %d %d" % (c1,c2))
    cubit.set_modified()
    cubit.cmd("create surface skin curve %d %d\n" % (c1,c2)) 
    if not cubit.is_modified():
        surfID = -1
    else:
        surfID = cubit.get_last_id("surface")
    printM("new Surface " + str(surfID) + " created")
    return surfID
    
def createSurfaceCurve(c):
    printM("create surface curve %d" % (c))
    cubit.cmd("create surface curve %d\n" % (c))
    surfID = cubit.get_last_id("surface")
    printM("new Surface " + str(surfID) + " created")
    return surfID

def createSurfaceCurve2(curves):
    cc = ""
    for c in curves:
        cc += str(c) + " "
    printM("create surface curve %s" % (cc))
    cubit.cmd("create surface curve %s\n" % (cc))
    surfID = cubit.get_last_id("surface")
    printM("new Surface " + str(surfID) + " created")
    return surfID
    
def createVolume(surfs):
    ss = ""
    for s in surfs:
        ss += str(s) + " "
    printM("create volume surface %s heal keep" % (ss))
    cubit.cmd("create volume surface %s heal keep \n" % (ss))
    volID = cubit.get_last_id("volume")
    printM("new Volume " + str(volID) + " created")
    return volID
    
def union(vols):
    vs = ""
    for v in vols:
        vs += str(v) + " "
    printM("unite volume %s keep" % (vs))
    cubit.cmd("unite volume %s keep \n" % (vs))

def stepUnion(vols):
    volID = vols[0]
    for x in range(0, len(vols)-1):
        cubit.set_modified()
        printM("unite volume %d %d" % (volID,vols[x+1]))
        cubit.cmd("unite volume %d %d" % (volID,vols[x+1]))
        if not cubit.is_modified():
            print "erro union volumes %d %d" % (volID,vols[x+1])
            return -1
            #printM("move Volume %d x %f" % (vols[x+1],0.01))
            #cubit.cmd("move Volume %d x %f y %f z %f\n" % (vols[x+1],0.0,0.0,0.01))
            #printM("unite volume %d %d" % (vols[0],vols[x+1]))
            #cubit.cmd("unite volume %d %d \n" % (vols[0],vols[x+1]))
    return volID
    
def createCone(h,r,t):
    printM("create frustum height %f radius %f top %f" % (h,r,t))
    cubit.cmd("create frustum height %f radius %f top %f \n" % (h,r,t))
    volID = cubit.get_last_id('volume')
    printM("new Volume " + str(volID) + " created")
    return volID
    
def alignVolume(v,s1,s2): 
    printM("align Volume %d surface %d with surface %d " % (v,s1,s2))
    cubit.cmd("align Volume %d surface %d with surface %d  \n" % (v,s1,s2))

def deleteVolume(v):
    printM("delete Volume %d" % (v))
    cubit.cmd("delete Volume %d \n" % (v))

def deleteSurface(s):
    printM("delete Surface %d" % (s))
    cubit.cmd("delete Surface %d \n" % (s))

def deleteCurve(c):
    printM("delete Curve %d" % (c))
    cubit.cmd("delete Curve %d \n" % (c))

def deleteVertex(v):
    printM("delete Vertex %d" % (v))
    cubit.cmd("delete Vertex %d \n" % (v))

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

def printM(s):
    if manualPrint:
        print s
        