import geo
import geo2
from vessel import Vessel
import cubit
import math
import random
from settings import *
from branch import Branch
from geo import distance, getAngles
from tube import Tube
import tube

class Arvore(object):
    
    def __init__(self, points, links):
        for p in points:
            p[0] = p[0]*scale
            p[1] = p[1]*scale
            p[2] = p[2]*scale
        for l in links:
            l[2] = l[2]*scale
        self.points = points
        self.root = Branch(int(links[0][0]),int(links[0][1]), links[0][2])
        heap = [self.root]
        self.points.append(points[self.root.initialPoint])
        while not (len(heap) == 0):
            atual = heap.pop()
            self.points.append(points[atual.finalPoint])
            for l in links:
                ip = int(l[0])
                if ip == atual.finalPoint:
                    son = Branch(int(l[0]),int(l[1]), l[2])
                    son.root = atual
                    if atual.son1 == None:
                        atual.son1 = son
                        heap.append(son)
                    else:
                        atual.son2 = son
                        heap.append(son)
        self.sort()
        self.vol = None
                
    def draw(self):
        for p in self.points:
            geo.createVertex(p[0],p[1],p[2])
        heap = [self.root]
        while not (len(heap) == 0):
            atual = heap.pop()
            lineID = geo.createLine(atual.initialPoint+1, atual.finalPoint+1)
            atual.index = lineID
            if atual.son1 is not None:
                heap.append(atual.son1)
            if atual.son2 is not None:
                heap.append(atual.son2)   
    
    def sort(self):
        heap = [self.root]
        sortedPoints = []
        sortedPoints.append(self.points[self.root.initialPoint])
        self.root.initialPoint = 0
        while not (len(heap) == 0):
            atual = heap.pop()
            sortedPoints.append(self.points[atual.finalPoint])
            atual.finalPoint = len(sortedPoints)-1
            if atual.son1 is not None:
                atual.son1.initialPoint = atual.finalPoint
                heap.append(atual.son1)
            if atual.son2 is not None:
                atual.son2.initialPoint = atual.finalPoint
                heap.append(atual.son2)  
        self.points = sortedPoints
                           
    def imprime(self):
        heap = [self.root]
        while not (len(heap) == 0):
            atual = heap.pop()
            pi = self.points[atual.initialPoint]
            pf = self.points[atual.finalPoint]
            lenght = geo.distance(pi, pf)
            if lenght < atual.radius*3:
                print "Segmento muito pequeno "
                print atual.initialPoint+1
            #print str(atual.initialPoint) + str(self.points[atual.initialPoint])+ " " + str(atual.finalPoint) + str(self.points[atual.finalPoint])
            #print str(atual.initialPoint) + " " + str(atual.finalPoint)
            if atual.son1 is not None:
                heap.append(atual.son1)
            if atual.son2 is not None:
                heap.append(atual.son2)
                
    def fixSizes(self):
        heap = [self.root]
        while not (len(heap) == 0):
            atual = heap.pop()
            ip = self.points[atual.initialPoint]
            fp = self.points[atual.finalPoint]
            lenght = geo.distance(ip, fp)
            if lenght < atual.radius*4:
                t=(atual.radius*4.5)/lenght
                npf = [(1-t)*ip[0]+t*fp[0],((1-t)*ip[1]+t*fp[1]),((1-t)*ip[2]+t*fp[2])]
                #self.points[atual.finalPoint] = npf
                xyz = [npf[0]-fp[0], npf[1]-fp[1], npf[2]-fp[2]]
                if atual.son1 is not None:
                    self.move(atual.son1,xyz)
                if atual.son2 is not None:
                    self.move(atual.son2,xyz)
            if atual.son1 is not None:
                heap.append(atual.son1)
            if atual.son2 is not None:
                heap.append(atual.son2)
                
    def move(self, branch, xyz):
        heap = [branch]
        while not (len(heap) == 0):
            atual = heap.pop()
            if atual.root.son1 == atual:
                self.points[atual.initialPoint][0] += xyz[0]
                self.points[atual.initialPoint][1] += xyz[1]
                self.points[atual.initialPoint][2] += xyz[2]
            self.points[atual.finalPoint][0] += xyz[0]
            self.points[atual.finalPoint][1] += xyz[1]
            self.points[atual.finalPoint][2] += xyz[2]
            if atual.son1 is not None:
                heap.append(atual.son1)
            if atual.son2 is not None:
                heap.append(atual.son2)
    
    def checkAngles(self):
        heap = []
        if self.root.son1 is not None:
            heap.append(self.root.son1)
        if self.root.son2 is not None:
            heap.append(self.root.son2)
        while not (len(heap) == 0):
            atual = heap.pop()
            initialRoot = self.points[atual.root.initialPoint]
            finalRoot = self.points[atual.root.finalPoint]
            anglesRoot = geo.getAngles(initialRoot, finalRoot)
            ay = geo.radiansToDegree(anglesRoot[0])
            initial = self.points[atual.initialPoint]
            final = self.points[atual.finalPoint]
            angles = geo.getAngles(initial, final)
            ay2 = geo.radiansToDegree(angles[0])
            a = ay2-ay
            dif = abs(a)
            if atual.son1 is not None:
                heap.append(atual.son1)
            if atual.son2 is not None:
                heap.append(atual.son2)  
        print "FIM"   
        
    def split(self):
        heap = [self.root]
        while not (len(heap) == 0):
            atual = heap.pop()
            if atual.son1 is not None and atual.son2 is not None:
                self.bif(atual) 
#                 sonVertex = geo.createVertex(coords[0],coords[1],coords[2])
            if atual.son1 is not None:
                heap.append(atual.son1)
            if atual.son2 is not None:
                heap.append(atual.son2) 
        self.sort() 
        
    def split2(self):
        heap = [self.root]
        while not (len(heap) == 0):
            atual = heap.pop()
            if atual.son1 is not None and atual.son2 is not None:
                #self.bif(atual) 
                self.bif2(atual)
                
            if atual.son1 is not None:
                heap.append(atual.son1)
            if atual.son2 is not None:
                heap.append(atual.son2) 
        self.sort() 
    
    def bif2(self, branch):
        ip = self.points[branch.initialPoint]
        fp = self.points[branch.finalPoint]
        ips = self.points[branch.son1.initialPoint]
        fps = self.points[branch.son1.finalPoint]
        ips2 = self.points[branch.son2.initialPoint]
        fps2 = self.points[branch.son2.finalPoint]
        lenght = geo.distance(ip, fp)
        t=(lenght-branch.radius)/lenght
        ipaux = [(1-t)*ip[0]+t*fp[0],((1-t)*ip[1]+t*fp[1]),((1-t)*ip[2]+t*fp[2])]
        fp = ipaux
        self.points[branch.finalPoint] = ipaux
        lenght = geo.distance(ipaux, fps)
        t=(branch.radius*2)/lenght
        fpaux = [(1-t)*ipaux[0]+t*fps[0],((1-t)*ipaux[1]+t*fps[1]),((1-t)*ipaux[2]+t*fps[2])]
        ips = fpaux
        
        ipaux2 = [fp[0],fp[1],fp[2]]
        lenght = geo.distance(fp, fps2)
        t=(branch.radius*2)/lenght
        fpaux2 = [(1-t)*fp[0]+t*fps2[0],((1-t)*fp[1]+t*fps2[1]),((1-t)*fp[2]+t*fps2[2])]
        ips2 = fpaux2
        
#         vecR = geo2.getVector(ip, fp)
#         vecAux = geo2.getVector(ipaux, fpaux)
#         vecS = geo2.getVector(ips, fps)
#         vecX = [1,0,0]
#         vecY = [0,1,0]
#         vecZ = [0,0,1]
#         if branch == self.root:
#             print "Root"
#             print geo2.getAngle(vecAux, vecX)
#             print geo2.getAngle(vecAux, vecY)
#             print geo2.getAngle(vecAux, vecZ)
#             print ""
#             print "Son"
#             print geo2.getAngle(vecS, vecX)
#             print geo2.getAngle(vecS, vecY)
#             print geo2.getAngle(vecS, vecZ)

        self.points.append(fpaux)
        geo.createVertex(fpaux[0], fpaux[1], fpaux[2])
        new1 = Branch(branch.finalPoint, len(self.points)-1, branch.radius*0.75)
        new1.root = branch
        new1.son1 = branch.son1
        new1.son1.initialPoint = len(self.points)-1
        new1.son1.root = new1
        branch.son1 = new1
           
        self.points.append(fpaux2)
        geo.createVertex(fpaux2[0], fpaux2[1], fpaux2[2])
        new2 = Branch(branch.finalPoint, len(self.points)-1, branch.radius*0.75)
        new2.root = branch
        new2.son2 = branch.son2
        new2.son2.initialPoint = len(self.points)-1
        new2.son2.root = new2
        branch.son2 = new2
        
    def smoth(self):
        change = True
        count = 0
        while change:
            change = False
            heap = []
            if self.root.son1 is not None:
                heap.append(self.root.son1)
            if self.root.son2 is not None:
                heap.append(self.root.son2)
            while not (len(heap) == 0):
                atual = heap.pop()
                initialRoot = self.points[atual.root.initialPoint]
                finalRoot = self.points[atual.root.finalPoint]
                anglesRoot = geo.getAngles(initialRoot, finalRoot)
                ay = geo.radiansToDegree(anglesRoot[0])
                initial = self.points[atual.initialPoint]
                final = self.points[atual.finalPoint]
                angles = geo.getAngles(initial, final)
                ay2 = geo.radiansToDegree(angles[0])
                a = ay2-ay
                if a > 180:
                    a -= 360
                if a < -180:
                    a += 360
                dif = abs(a)
                rot = False
                bothSons = False
                #if atual.root.son1 is not None and atual.root.son2 is not None:
                #    bothSons = True
                if dif > 25: #and not bothSons:
                    print dif
                    change = True
                    positive = True
                    if a < 0:
                        positive = False
                    if dif > 60:
                        angleRot = 20
                    else:
                        angleRot = 10
                    if positive:
                        angleRot = -angleRot
                    rot = True
#                 if atual.son1 is not None and atual.son2 is not None:
#                     if dif < 10:
#                         print " DIFFFF"
#                         change = True
#                         positive = True
#                         if a < 0:
#                             positive = False
#                         angleRot = 20
#                         if not positive:
#                             angleRot = -angleRot
#                         rot = True
                if rot:  
                    count = count+1
                    print count  
                    vertexCoords = [(initial[0]+2*final[0])/3,(initial[1]+2*final[1])/3,(initial[2]+2*final[2])/3]
                    vertexCoords[0] = vertexCoords[0]
                    vertexCoords[1] = vertexCoords[1]
                    vertexCoords[2] = vertexCoords[2]
                    geo.rotate(vertexCoords,initial,angleRot) 
                    self.points.append(vertexCoords)
                    rootRadius = atual.root.radius
                    if atual.root.son1 is not None and atual.root.son2 is not None:
                        rootRadius = rootRadius*0.75
                    new = Branch(len(self.points)-1, atual.finalPoint, (2*atual.radius+rootRadius)/3)
                    new.root = atual
                    new.son1 = atual.son1
                    new.son2 = atual.son2
                    if atual.son1 is not None:
                        heap.append(atual.son1)
                        new.son1.root = new
                    if atual.son2 is not None:
                        heap.append(atual.son2)
                        new.son2.root = new
                        atual.son2 = None 
                    atual.finalPoint = new.initialPoint
                    atual.son1 = new
                else:
                    if atual.son1 is not None:
                        heap.append(atual.son1)
                    if atual.son2 is not None:
                        heap.append(atual.son2)
        self.sort()  
        print "FIM"  
    
    def smoth2(self):
        heap = [self.root]
        while not (len(heap) == 0):
            atual = heap.pop()
            if atual.son1 is not None and atual.son2 is not None:
                rSon1 = atual.son1.son1  
                r1i = self.points[atual.son1.initialPoint] 
                r1f = self.points[atual.son1.finalPoint]
                f1i = self.points[rSon1.initialPoint]
                f1f = self.points[rSon1.finalPoint]
                t = 0.9
                ref = [(1-t)*f1i[0]+t*f1f[0],((1-t)*f1i[1]+t*f1f[1]),((1-t)*f1i[2]+t*f1f[2])]
                self.splitBranch(rSon1, ref)
                f1f = self.points[rSon1.finalPoint]
                v1 = geo.createVertex(r1i[0], r1i[1], r1i[2])
                v2 = geo.createVertex(f1i[0], f1i[1], f1i[2])
                v3 = geo.createVertex(f1f[0], f1f[1], f1f[2])
                spAux = geo.createSpline(v1, v2, v3)
                sp = geo.splitCurve(spAux, v2)
                newV = geo.createVertexOnCurveFraction(sp, 0.66)
                newPoint = cubit.vertex(newV).coordinates()
                self.points.append(newPoint)
                self.splitBranch(rSon1, newPoint)
                newV = geo.createVertexOnCurveFraction(sp, 0.33)
                newPoint = cubit.vertex(newV).coordinates()
                self.points.append(newPoint)
                self.splitBranch(rSon1, newPoint)
                
                rSon2 = atual.son2.son1
                r2i = self.points[atual.son2.initialPoint] 
                r2f = self.points[atual.son2.finalPoint]
                f2i = self.points[rSon2.initialPoint]
                f2f = self.points[rSon2.finalPoint]
                t = 0.9
                ref = [(1-t)*f2i[0]+t*f2f[0],((1-t)*f2i[1]+t*f2f[1]),((1-t)*f2i[2]+t*f2f[2])]
                self.splitBranch(rSon2, ref)
                f2f = self.points[rSon2.finalPoint]
                v1 = geo.createVertex(r2i[0], r2i[1], r2i[2])
                v2 = geo.createVertex(f2i[0], f2i[1], f2i[2])
                v3 = geo.createVertex(f2f[0], f2f[1], f2f[2])
                spAux = geo.createSpline(v1, v2, v3)
                sp = geo.splitCurve(spAux, v2)
                newV = geo.createVertexOnCurveFraction(sp, 0.66)
                newPoint = cubit.vertex(newV).coordinates()
                self.points.append(newPoint)
                self.splitBranch(rSon2, newPoint)
                newV = geo.createVertexOnCurveFraction(sp, 0.33)
                newPoint = cubit.vertex(newV).coordinates()
                self.points.append(newPoint)
                self.splitBranch(rSon2, newPoint)
                
                heap.append(rSon1)
                heap.append(rSon2)
            else:
                if atual.son1 is not None:
                    heap.append(atual.son1)
                if atual.son2 is not None:
                    heap.append(atual.son2) 
                    
    def save(self):
        pointsFile = open(folder+"pontos2.txt", 'w')
        linesFile = open(folder+"linhas2.txt", 'w')
        for p in self.points:
            pointsFile.write("%f\t%f\t%f\n" % (p[0],p[1],p[2]))
        pointsFile.close()
        heap = [self.root]
        while not (len(heap) == 0):
            atual = heap.pop()
            linesFile.write("%d\t%d\t%f\n" % (atual.initialPoint,atual.finalPoint,atual.radius))
            if atual.son1 is not None:
                heap.append(atual.son1)
            if atual.son2 is not None:
                heap.append(atual.son2) 
        linesFile.close()
        
    def save2(self):
        pointsFile = open(folder+"pontos3.txt", 'w')
        linesFile = open(folder+"linhas3.txt", 'w')
        for p in self.points:
            pointsFile.write("%f\t%f\t%f\n" % (p[0],p[1],p[2]))
        pointsFile.close()
        heap = [self.root]
        while not (len(heap) == 0):
            atual = heap.pop()
            linesFile.write("%d\t%d\t%f\n" % (atual.initialPoint,atual.finalPoint,atual.radius))
            if atual.son1 is not None:
                heap.append(atual.son1)
            if atual.son2 is not None:
                heap.append(atual.son2) 
        linesFile.close()
        
    def makeGeometry(self):
        for p in self.points:
            geo.createVertex(p[0],p[1],p[2])
        heap = [self.root]
        while not (len(heap) == 0):
            atual = heap.pop()
            lineID = geo.createLine(atual.initialPoint+1, atual.finalPoint+1)
            atual.index = lineID
            if atual.son1 is not None:
                heap.append(atual.son1)
            if atual.son2 is not None:
                heap.append(atual.son2) 
        
        heap = [self.root]
        count = 0
        while not (len(heap) == 0):
            count += 1
            if count < 2:
                atual = heap.pop()
            else:
                break
            if atual.root is None:
                if atual.son1 is not None and atual.son2 is not None:
                    ipr = self.points[atual.initialPoint]
                    fpr = self.points[atual.finalPoint]
                    fp = self.points[atual.son1.finalPoint]
                    lenght = geo.distance(ipr, fpr)
                    t=(-atual.radius)/lenght
                    ref = [(1-t)*ipr[0]+t*fpr[0],((1-t)*ipr[1]+t*fpr[1]),((1-t)*ipr[2]+t*fpr[2])]
                    u0 = [fpr[0]-ref[0],fpr[1]-ref[1],fpr[2]-ref[2]]
                    u1 = [ref[0]-fp[0],ref[1]-fp[1],ref[2]-fp[2]]
                    u = geo2.crossProduct(u0, u1)
                    u = geo2.getVector4_3(u)
                    geo2.normalizeVector(u)
                    arcsi = self.genArc(atual, atual.initialPoint+1,ref,u,atual.radius)
                    t=(lenght-atual.radius)/lenght
                    ref = [(1-t)*ipr[0]+t*fpr[0],((1-t)*ipr[1]+t*fpr[1]),((1-t)*ipr[2]+t*fpr[2])]
                    arcsf = self.genArc(atual, atual.finalPoint+1,ref,u,atual.radius)
                    atual.tube = Tube(arcsi,arcsf,atual.initialPoint+1,atual.finalPoint+1)
                    t=(lenght-atual.radius*0.5)/lenght
                    ref = [(1-t)*ipr[0]+t*fpr[0],((1-t)*ipr[1]+t*fpr[1]),((1-t)*ipr[2]+t*fpr[2])]
                    coords = geo2.rotateByAxis(geo2.getVector4_3(ref), geo2.getVector4_3(fpr), u, -90)
                    v1 = geo.createVertex(coords[0],coords[1],coords[2])
                    coords = geo2.rotateByAxis(geo2.getVector4_3(ref), geo2.getVector4_3(fpr), u, 90)
                    v2 = geo.createVertex(coords[0],coords[1],coords[2])
                    self.genSons(atual, arcsf,v1,v2,u) 
                    if atual.son1.son1 is not None:
                        heap.append(atual.son1.son1)
                    if atual.son1.son2 is not None:
                        heap.append(atual.son1.son2) 
                    if atual.son2.son1 is not None:
                        heap.append(atual.son2.son1)
                    if atual.son2.son2 is not None:
                        heap.append(atual.son2.son2)
                else:
                    if atual.son1 is not None:
                        heap.append(atual.son1)     
            else:
                if atual.son1 is not None and atual.son2 is not None:
                    ipr = self.points[atual.initialPoint]
                    fpr = self.points[atual.finalPoint]
                    fp = self.points[atual.son1.finalPoint]
                    lenght = geo.distance(ipr, fpr)
                    t=(-atual.radius)/lenght
                    ref = [(1-t)*ipr[0]+t*fpr[0],((1-t)*ipr[1]+t*fpr[1]),((1-t)*ipr[2]+t*fpr[2])]
                    u0 = [fpr[0]-ref[0],fpr[1]-ref[1],fpr[2]-ref[2]]
                    u1 = [ref[0]-fp[0],ref[1]-fp[1],ref[2]-fp[2]]
                    u = geo2.crossProduct(u0, u1)
                    u = geo2.getVector4_3(u)
                    geo2.normalizeVector(u)
                    arcsi = []
                    for i in range(0, 4):
                        arcsi.append([atual.root.tube.finalVertexs[i], atual.root.tube.finalArcs[i]])
                    t=(lenght-atual.radius)/lenght
                    ref = [(1-t)*ipr[0]+t*fpr[0],((1-t)*ipr[1]+t*fpr[1]),((1-t)*ipr[2]+t*fpr[2])]
                    arcsf = self.genArc(atual, atual.finalPoint+1,ref,u,atual.radius)
                    atual.tube = Tube(arcsi,arcsf,atual.initialPoint+1,atual.finalPoint+1)  
                    t=(lenght-atual.radius*0.5)/lenght
                    ref = [(1-t)*ipr[0]+t*fpr[0],((1-t)*ipr[1]+t*fpr[1]),((1-t)*ipr[2]+t*fpr[2])]
                    coords = geo2.rotateByAxis(geo2.getVector4_3(ref), geo2.getVector4_3(fpr), u, -90)
                    v1 = geo.createVertex(coords[0],coords[1],coords[2])
                    coords = geo2.rotateByAxis(geo2.getVector4_3(ref), geo2.getVector4_3(fpr), u, 90)
                    v2 = geo.createVertex(coords[0],coords[1],coords[2])
                    self.genSons(atual, arcsf,v1,v2,u)
                    if atual.son1.son1 is not None:
                        heap.append(atual.son1.son1)
                    if atual.son1.son2 is not None:
                        heap.append(atual.son1.son2) 
                    if atual.son2.son1 is not None:
                        heap.append(atual.son2.son1)
                    if atual.son2.son2 is not None:
                        heap.append(atual.son2.son2)
                else:
                    ipr = self.points[atual.initialPoint]
                    fpr = self.points[atual.finalPoint]
                    if atual.son1 is not None:
                        fp = self.points[atual.son1.finalPoint]
                        lenght = geo.distance(ipr, fpr)
                        t=(-atual.radius)/lenght
                        ref = [(1-t)*ipr[0]+t*fpr[0],((1-t)*ipr[1]+t*fpr[1]),((1-t)*ipr[2]+t*fpr[2])]
                        u0 = [fpr[0]-ref[0],fpr[1]-ref[1],fpr[2]-ref[2]]
                        u1 = [ref[0]-fp[0],ref[1]-fp[1],ref[2]-fp[2]]
                        u = geo2.crossProduct(u0, u1)
                        u = geo2.getVector4_3(u)
                        geo2.normalizeVector(u)
                    else:
                        fp = self.points[atual.root.initialPoint]
                        lenght = geo.distance(ipr, fpr)
                        t=(-atual.radius)/lenght
                        ref = [(1-t)*ipr[0]+t*fpr[0],((1-t)*ipr[1]+t*fpr[1]),((1-t)*ipr[2]+t*fpr[2])]
                        fpr = self.points[atual.finalPoint]
                        u0 = [fpr[0]-ref[0],fpr[1]-ref[1],fpr[2]-ref[2]]
                        u1 = [ref[0]-fp[0],ref[1]-fp[1],ref[2]-fp[2]]
                        angle = self.getAngle(atual)
                        if angle < 0:
                            u = geo2.crossProduct(u0, u1)
                        else:
                            u = geo2.crossProduct(u1, u0)
                        u = geo2.getVector4_3(u)
                        geo2.normalizeVector(u)
                    arcsi = []
                    for i in range(0, 4):
                        arcsi.append([atual.root.tube.finalVertexs[i], atual.root.tube.finalArcs[i]])
                    t=(lenght-atual.radius)/lenght
                    ref = [(1-t)*ipr[0]+t*fpr[0],((1-t)*ipr[1]+t*fpr[1]),((1-t)*ipr[2]+t*fpr[2])]
                    arcsf = self.genArc(atual, atual.finalPoint+1,ref,u,atual.radius)
                    
                    atual.tube = Tube(arcsi,arcsf,atual.initialPoint+1,atual.finalPoint+1)   
                    if atual.son1 is not None:
                        #atual.son1.arcsi = arcsf
                        heap.append(atual.son1)
        #self.genSurfaces() 
        #self.genVolumes()

    def genVolumes(self):
        heap = [self.root]
        while not (len(heap) == 0):
            atual = heap.pop()
            atual.tube.genVolume()
            if atual == self.root:
                self.vol = atual.tube.vol
            else:
                geo.imprintVolumes(self.vol, atual.tube.vol)
                geo.mergeVolumes(self.vol, atual.tube.vol)
                self.vol = geo.unionVolumes(self.vol, atual.tube.vol)
            if atual.son1 is not None:
                heap.append(atual.son1)
            if atual.son2 is not None:
                heap.append(atual.son2) 
#         geo.imprintMergeAll() 
#         geo.unionAll()
#         cubit.cmd('Color Define "%s" RGB %f %f %f' % ("darkred", 0.75,0,0))
#         geo.colorVolume(self.vol,'user "darkred"')
    
    def genSurfaces(self):
        heap = [self.root]
        count = 0
        while not (len(heap) == 0):
            count += 1
            atual = heap.pop()
            print atual.initialPoint+1
            #if atual.son1 is not None:
            atual.tube.genSurfaces()
            if atual.son1 is not None:
                heap.append(atual.son1)
            if atual.son2 is not None:
                heap.append(atual.son2) 
    
    def genSons(self, atual, arcsi,v1,v2,u2):
        eli1 = geo.createEllipse(arcsi[3][0],v1,atual.finalPoint+1)
        eli2 = geo.createEllipse(arcsi[1][0],v1,atual.finalPoint+1)
        ip = self.points[atual.son1.initialPoint]
        fp = self.points[atual.son1.finalPoint]
        fps = self.points[atual.son1.son1.finalPoint]
        lenght = geo.distance(ip, fp)
        t=(lenght-atual.son1.radius)/lenght
        ref = [(1-t)*ip[0]+t*fp[0],((1-t)*ip[1]+t*fp[1]),((1-t)*ip[2]+t*fp[2])]
        u0 = [fp[0]-ip[0],fp[1]-ip[1],fp[2]-ip[2]]
        u1 = [fp[0]-fps[0],fp[1]-fps[1],fp[2]-fps[2]]
        u = geo2.crossProduct(u0, u1)
        u = geo2.getVector4_3(u)
        geo2.normalizeVector(u)
        arcsf1 = self.genArc2(atual.son1, atual.son1.finalPoint+1,ref,u2,atual.son1.radius,1)
        atual.son1.tube = Tube([[arcsi[3][0],eli1],[v1,eli2],arcsi[1],arcsi[2]], arcsf1, atual.son1.initialPoint+1, atual.son1.finalPoint+1)
         
        eli3 = geo.createEllipse(arcsi[3][0],v2,atual.finalPoint+1)
        eli4 = geo.createEllipse(arcsi[1][0],v2,atual.finalPoint+1)
        ip = self.points[atual.son2.initialPoint]
        fp = self.points[atual.son2.finalPoint]
        fps = self.points[atual.son2.son1.finalPoint]
        lenght = geo.distance(ip, fp)
        t=(lenght-atual.son2.radius)/lenght
        ref = [(1-t)*ip[0]+t*fp[0],((1-t)*ip[1]+t*fp[1]),((1-t)*ip[2]+t*fp[2])]
        u0 = [fp[0]-ip[0],fp[1]-ip[1],fp[2]-ip[2]]
        u1 = [fp[0]-fps[0],fp[1]-fps[1],fp[2]-fps[2]]
        u = geo2.crossProduct(u0, u1)
        u = geo2.getVector4_3(u)
        geo2.normalizeVector(u)
        arcsf2 = self.genArc2(atual.son2, atual.son2.finalPoint+1,ref,u2,atual.son2.radius,2)
        atual.son2.tube = Tube([[arcsi[1][0],eli3],[v2,eli4],arcsi[3],arcsi[0]], arcsf2, atual.son2.initialPoint+1, atual.son2.finalPoint+1)

    def getAngle(self, son):
        initialRoot = self.points[son.root.initialPoint]
        finalRoot = self.points[son.root.finalPoint]
        anglesRoot = geo.getAngles(initialRoot, finalRoot)
        ay = geo.radiansToDegree(anglesRoot[0])
        initial = self.points[son.initialPoint]
        final = self.points[son.finalPoint]
        angles = geo.getAngles(initial, final)
        ay2 = geo.radiansToDegree(angles[0])
        a = ay2-ay
        if a > 180:
            a -= 360
        if a < -180:
            a += 360
            
        az = geo.radiansToDegree(anglesRoot[1])
        initial = self.points[son.initialPoint]
        final = self.points[son.finalPoint]
        angles = geo.getAngles(initial, final)
        az2 = geo.radiansToDegree(angles[1])
        b = az2-az
        if b > 180:
            b -= 360
        if b < -180:
            b += 360
        return a+b
    
    def bif(self, branch):
        pi = self.points[branch.initialPoint]
        pf = self.points[branch.finalPoint]
        pi1 = self.points[branch.son1.initialPoint]
        pf1 = self.points[branch.son1.finalPoint]
        pi2 = self.points[branch.son2.initialPoint]
        pf2 = self.points[branch.son2.finalPoint]
        lenght = geo.distance(pi, pf)
            
        t=(lenght-1*branch.radius)/lenght
        newPf = [(1-t)*pi[0]+t*pf[0],((1-t)*pi[1]+t*pf[1]),((1-t)*pi[2]+t*pf[2])]
        self.points[branch.finalPoint] = newPf
        #geo.createVertex(newPf[0],newPf[1],newPf[2])
        t=(lenght+0*branch.radius)/lenght
        newPoint1 = [(1-t)*pi[0]+t*pf[0],((1-t)*pi[1]+t*pf[1]),((1-t)*pi[2]+t*pf[2])]
#         u0 = [newPf[0]-newPoint1[0],newPf[1]-newPoint1[1],newPf[2]-newPoint1[2]]
#         u1 = [newPoint1[0]-pf1[0],newPoint1[1]-pf1[1],newPoint1[2]-pf1[2]]
#         u0 = [newPoint1[0]-newPf[0],newPoint1[1]-newPf[1],newPoint1[2]-newPf[2]]
#         u1 = [pf1[0]-newPf[0],pf1[1]-newPf[1],pf1[2]-newPf[2]]
        u0 = [pf1[0]-newPf[0],pf1[1]-newPf[1],pf1[2]-newPf[2]]
        u1 = [pf2[0]-newPf[0],pf2[1]-newPf[1],pf2[2]-newPf[2]]
        u = geo2.crossProduct(u0, u1)
        u = geo2.getVector4_3(u)
        geo2.normalizeVector(u)
        newPoint1 = geo2.rotateByAxis(geo2.getVector4_3(newPoint1), geo2.getVector4_3(newPf), u, 40)
        #geo.createVertex(newPoint1[0],newPoint1[1],newPoint1[2])
        t=(lenght+0*branch.radius)/lenght
        newPoint2 = [(1-t)*pi[0]+t*pf[0],((1-t)*pi[1]+t*pf[1]),((1-t)*pi[2]+t*pf[2])]
#         u0 = [newPoint2[0]-newPf[0],newPoint2[1]-newPf[1],newPoint2[2]-newPf[2]]
#         u1 = [pf2[0]-newPf[0],pf2[1]-newPf[1],pf2[2]-newPf[2]]
#         u = geo2.crossProduct(u0, u1)
#         u = geo2.getVector4_3(u)
#         geo2.normalizeVector(u)
        newPoint2 = geo2.rotateByAxis(geo2.getVector4_3(newPoint2), geo2.getVector4_3(newPf), u, -40)
        #geo.createVertex(newPoint2[0],newPoint2[1],newPoint2[2])
        d1 = geo.distance(pf1, newPoint1)+geo.distance(pf2, newPoint2)
        d2 = geo.distance(pf1, newPoint2)+geo.distance(pf2, newPoint1)
        if d1>d2:
            aux = newPoint1
            newPoint1 = newPoint2
            newPoint2 = aux
        self.points.append(newPoint1)
        radius1 = branch.radius
        radius2 = branch.radius
#         radius1 = branch.son1.radius
#         radius2 = branch.son2.radius
#         sumRadius = radius1 + radius2
#         while sumRadius > 1.2:
#             radius1 = radius1*0.95
#             radius2 = radius2*0.95
#             sumRadius = radius1 + radius2
        new1 = Branch(branch.finalPoint, len(self.points)-1,radius1)
        new1.root = branch
        new1.son1 = branch.son1
        new1.son1.initialPoint = len(self.points)-1
        new1.son1.root = new1
        branch.son1 = new1
           
        self.points.append(newPoint2)
        new2 = Branch(branch.finalPoint, len(self.points)-1, radius2)
        new2.root = branch
        new2.son2 = branch.son2
        new2.son2.initialPoint = len(self.points)-1
        new2.son2.root = new2
        branch.son2 = new2
    
    def genArc(self, branch, point, ref, u, radius):
        ipr = self.points[point-1]
        coords = geo2.rotateByAxis(geo2.getVector4_3(ref), geo2.getVector4_3(ipr), u, -90)
        vAuxI1 = geo.createVertex(coords[0],coords[1],coords[2])
        circleI = geo.createCircleNormal2(point,vAuxI1,vAuxI1,radius,branch.index)
        vAuxI2 = geo.createVertexOnCurveFraction(circleI, 0.25)
        arc = geo.createCircleNormal2(point, vAuxI1, vAuxI2, radius, branch.index)
        arcI1 = [vAuxI1, arc]
        vAuxI3 = geo.createVertexOnCurveFraction(circleI, 0.5)
        arc = geo.createCircleNormal2(point, vAuxI2, vAuxI3, radius, branch.index)
        arcI2 = [vAuxI2, arc]
        vAuxI4 = geo.createVertexOnCurveFraction(circleI, 0.75)
        arc = geo.createCircleNormal2(point, vAuxI3, vAuxI4, radius, branch.index)
        arcI3 = [vAuxI3, arc]
        arc = geo.createCircleNormal2(point, vAuxI4, vAuxI1, radius, branch.index)
        arcI4 = [vAuxI4, arc]
        geo.deleteCurve2(circleI)

        arcs = [arcI1,arcI2,arcI3,arcI4]
        return arcs
    
    def genArc2(self, branch, point, ref, u, radius, son):
        ipr = self.points[branch.initialPoint]
        fpr = self.points[branch.finalPoint]
        lenght = geo.distance(ipr, fpr)
        if son == 1:
            off = 40
        else:
            off = -40
         
        t=(lenght-branch.radius*0.6)/lenght
        ref2 = [(1-t)*ipr[0]+t*fpr[0],((1-t)*ipr[1]+t*fpr[1]),((1-t)*ipr[2]+t*fpr[2])]
        ipr = self.points[point-1]
        coords = geo2.rotateByAxis(geo2.getVector4_3(ref2), geo2.getVector4_3(ipr), u, -90+off)
        vAux1 = geo.createVertex(coords[0],coords[1],coords[2])
#         coords = geo2.rotateByAxis(geo2.getVector4_3(ref2), geo2.getVector4_3(ipr), u, 90+off)
#         v2 = geo.createVertex(coords[0],coords[1],coords[2])
            
        ipr = self.points[point-1]
        coords = geo2.rotateByAxis(geo2.getVector4_3(ref), geo2.getVector4_3(ipr), u, -90+off)
        vAux2 = geo.createVertex(coords[0],coords[1],coords[2])
        circleI = geo.createCircleNormal2(point,vAux2,vAux2,radius,branch.root.index)
        
        vAux3 = geo.createVertexOnCurveFraction(circleI, 0.25)
        eliFull = geo.createEllipseFull(vAux3, vAux1, point)
        
        geo.deleteCurve2(circleI)
        geo.deleteVertex(vAux1)
        geo.deleteVertex(vAux2)
        
        v1 = vAux3
        v2 = geo.createVertexOnCurveFraction(eliFull, 0.25)
        arc = geo.createEllipse(v1,v2,point)
        arcI1 = [v1, arc]
        v3 = geo.createVertexOnCurveFraction(eliFull, 0.5)
        arc = geo.createEllipse(v3,v2,point)
        arcI2 = [v2, arc]
        v4 = geo.createVertexOnCurveFraction(eliFull, 0.75)
        arc = geo.createEllipse(v3,v4,point)
        arcI3 = [v3, arc]
        arc = geo.createEllipse(v1,v4,point)
        arcI4 = [v4, arc]
        

#         arc = geo.createEllipse(vAuxI2,v1,point)
#         arcI1 = [v1, arc]
#         
#         v2 = geo.createVertexOnCurveFraction(eliFull, 0.5)
#         arc = geo.createEllipse(vAuxI2,v2,point)
#         arcI2 = [vAuxI2, arc]
#         
#         vAuxI4 = geo.createVertexOnCurveFraction(circleI, 0.75)
#         arc = geo.createEllipse(vAuxI4,v2,point)
#         arcI3 = [v2, arc]
#         
#         arc = geo.createEllipse(vAuxI4,v1,point)
#         arcI4 = [vAuxI4, arc]
#         
#         geo.deleteCurve2(circleI)
#         geo.deleteCurve2(eliFull)
#         geo.deleteVertex(vAuxI1)
# 
        geo.deleteCurve2(eliFull)
        arcs = [arcI1,arcI2,arcI3,arcI4]
        return arcs

    def splitBranch(self,branch, point):
        self.points.append(point)
        radius = branch.radius
        new = Branch(len(self.points)-1, branch.finalPoint, radius)
        new.root = branch
        new.son1 = branch.son1
        new.son2 = branch.son2
        if branch.son1 is not None:
            new.son1.root = new
        if branch.son2 is not None:
            new.son2.root = new
            branch.son2 = None 
        branch.finalPoint = new.initialPoint
        branch.son1 = new
    
    def check(self):
        heap = [self.root]
        while not (len(heap) == 0):
            atual = heap.pop()
            print atual.initialPoint
            if atual.son1 is not None and atual.son2 is not None:
                print (atual.son1.radius/atual.radius+atual.son2.radius/atual.radius)/2
                heap.append(atual.son1)
                heap.append(atual.son2)  
            elif atual.son1 is not None:
                print atual.son1.radius/atual.radius
                heap.append(atual.son1) 