import geo
from vessel import Vessel
import cubit
import math
import random
from settings import *

class Tree(object):

    def __init__(self, points, links):
        self.points = points
        self.links = links
        self.vertices = []
        self.vessels = []
        self.linksFilhos=[[] for i in range(len(self.links))]            
        
    def preProcess(self):
        for f in range(0, len(self.links)):
            for p in range(0, len(self.links)):
                if int(self.links[f][0]) == int(self.links[p][1]):
                    self.linksFilhos[p].append(self.links[f])
                    pip = self.points[int(self.links[p][0])]
                    pi = self.points[int(self.links[f][0])]
                    pf = self.points[int(self.links[f][1])]
                    d = geo.distance(pi, pf)
                    angles = geo.getAngles(pip, pi)
                    ay = geo.radiansToDegree(angles[0])
                    az = geo.radiansToDegree(angles[1])
                       
                    deltaY = pi[0]-pip[0];
                    deltaX = pi[1]-pip[1];
                    #deltaZ = pi[2]-pip[2];
                    ay = math.atan2(deltaY, deltaX)
                    #az = math.atan2(deltaZ, deltaX)
                    ay = geo.radiansToDegree(ay)
                    deltaY = pf[0]-pi[0];
                    deltaX = pf[1]-pi[1];
                    #deltaZ = pf[2]-pi[2];
                    ay2 = math.atan2(deltaY, deltaX)
                    #az2 = math.atan2(deltaZ, deltaX)
                    ay2 = geo.radiansToDegree(ay2)
                    a = ay - ay2
                    if a > 180:
                        a -= 360
                    if a < -180:
                        a += 360
                    positive = True
                    if a < 0:
                        positive = False
                    if abs(a) > 50:
                        dif = abs(a) - 50
                        dif2 = abs(a)
                        if positive:
                            dif = -dif
                            dif2= -dif2
                        geo.rotate(pf,pi,dif)   
                        if (d*scale) < (self.links[f][2]*6):
                            geo.rotateAndMove(pf, pip, dif2, 0, (0.5*6)/scale)                        
        for p in range(0, len(self.links)):
            if len(self.linksFilhos[p]) == 2:
                f1 = self.linksFilhos[p][0]
                f2 = self.linksFilhos[p][1]
                f1i = self.points[int(f1[0])]
                f1f = self.points[int(f1[1])]
                f2i = self.points[int(f2[0])]
                f2f = self.points[int(f2[1])]
                pi = self.points[int(self.links[p][0])]
                pf = self.points[int(self.links[p][1])]
                deltaY = f1f[0]-f1i[0];
                deltaX = f1f[1]-f1i[1];
                #deltaZ = f1f[2]-f1i[2];
                ay = math.atan2(deltaY, deltaX)
                #az = math.atan2(deltaZ, deltaX)
                ay = geo.radiansToDegree(ay)
                deltaY = f2f[0]-f2i[0];
                deltaX = f2f[1]-f2i[1];
                #deltaZ = f2f[2]-f2i[2];
                ay2 = math.atan2(deltaY, deltaX)
                #az2 = math.atan2(deltaZ, deltaX)
                ay2 = geo.radiansToDegree(ay2)
                deltaY = f1i[0]-pi[0];
                deltaX = f1i[1]-pi[1];
                #deltaZ = f1i[2]-pi[2];
                ay3 = math.atan2(deltaY, deltaX)
                #az3 = math.atan2(deltaZ, deltaX)
                ay3 = geo.radiansToDegree(ay)
                a = ay - ay2
                if a > 180:
                    a -= 360
                if a < -180:
                    a += 360
                a2 = ay3 - ay
                if a2 > 180:
                    a2 -= 360
                if a2 < -180:
                    a2 += 360
                positive = True
                if a2 < 0:
                    positive = False
                if abs(a) < 20:
                    dif = 25 - abs(a)
                    if positive:
                        dif = -dif  
                    print "aaaaaaaaaa"
                    geo.rotate(f1f,f1i,dif)   
                  
    def genVertices(self, scale):
        for p in self.points:
            vertexID = geo.createVertex(p[0]*scale,p[1]*scale,p[2]*scale)
            self.vertices.append(vertexID)
            
    def genCenterlines(self):
        for l in self.links:
            root = None
            for x in range(0, len(self.vessels)):
                if int(l[0]+1) == self.vessels[x].vF:
                    root = self.vessels[x]
            lineID = geo.createLine(int(l[0]+1),int(l[1]+1))
            vessel = Vessel(lineID, int(l[0]+1),int(l[1]+1), l[2], root)
            self.vessels.append(vessel)
            if root is not None:
                if root.son1 is None:
                    root.son1 = vessel
                    root.nSons += 1
                else:
                    root.son2 = vessel 
                    root.nSons += 1
            
    def genSurfaces(self):
        for v in self.vessels:
            v.genSurfaces2()
            #v.genVolume()
            
    def union(self):
        geo.imprintMergeAll()
        self.vol = self.vessels[0].volID
        geo.unionAll()
        geo.colorVolume(self.vol, "red")
    
    def clean(self):
        for v in self.vessels:
            v.clean()
    
    def genPoints(self):
        vertexID = geo.createVertex(self.points[0][0],self.points[0][1],self.points[0][2])
        VI = vertexID
        points = []
        points.append(VI)
        for x in range(0, len(self.points)-1):
            vertexID = geo.createVertex(self.points[x+1][0],self.points[x+1][1],self.points[x+1][2])
            VF = vertexID
            lineID = geo.createLine(VI,VF)
            points.append(vertexID)
            dist = cubit.get_distance_between(VI, VF)
            #np = 4
            #np += int(dist/(self.raio*2))
            np = 10
            inc = 1.0/np
            for x in range(1, np):
                vertexID = geo.createVertexOnCurveFraction(lineID, inc*x)
                points.append(vertexID)
            points.append(VF)
            VI = VF
        pts = []
        for x in range(0, len(points)-1):
            pt1 = cubit.vertex(points[x]).coordinates()
            pt2 = cubit.vertex(points[x+1]).coordinates()
            for k in range(0, 100):
                rx = random.random()*0.5
                ry = random.random()*0.5
                rz = random.random()*0.5
                pt = [pt1[0]+rx, pt1[1]+ry, pt1[2]+rz]
                pts.append(pt)
        target = open("nodes.node", 'w')
        target.truncate()
        target.write("# Node count, 3 dim, no attribute, no boundary marker\n")
        target.write("%d  3  0  0\n" % (len(pts)))
        target.write("# Node index, node coordinates\n")
        for x in range(0, len(pts)):
            target.write("%d  %f %f %f\n" % (x+1,pts[x][0],pts[x][1],pts[x][2]))
        target.close()
        
    def mesh(self):
        cubit.cmd("volume %d scheme Tetmesh \n" % (self.vol))
        #cubit.cmd("set tetmesher interior points on \n")
        #cubit.cmd("set tetmesher optimize level 3 optimize overconstrained  off sliver  off  \n")
        #cubit.cmd("set tetmesher boundary recovery  off  \n")
        #cubit.cmd("volume %d  tetmesh growth_factor 1.0 \n" % (self.vol))
        #cubit.silent_cmd("mesh volume %d \n" % (self.vol))