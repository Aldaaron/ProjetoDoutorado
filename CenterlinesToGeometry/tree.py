import geo
from vessel import Vessel
import cubit
import math

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
                        if positive:
                            dif = -dif
                        geo.rotate(pf,pi,dif)                   
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
                if abs(a) < 10:
                    dif = 10 - abs(a)
                    if positive:
                        dif = -dif  
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
            
    def genVolume(self):
        for v in self.vessels:
            cod = v.genVolume2()
            if cod == -1:
                return
            
    def union(self):
        vols = []
        for v in self.vessels:
            vols.append(v.volID)
        #geo.union(vols)
        geo.stepUnion(vols)
        self.vol = vols[0]
    
    def clean(self):
        for v in self.vessels:
            v.clean()
            
    
    def mesh(self):
        cubit.cmd("volume %d scheme Tetmesh \n" % (self.vol))
        #cubit.cmd("set tetmesher interior points on \n")
        #cubit.cmd("set tetmesher optimize level 3 optimize overconstrained  off sliver  off  \n")
        #cubit.cmd("set tetmesher boundary recovery  off  \n")
        #cubit.cmd("volume %d  tetmesh growth_factor 1.0 \n" % (self.vol))
        #cubit.silent_cmd("mesh volume %d \n" % (self.vol))