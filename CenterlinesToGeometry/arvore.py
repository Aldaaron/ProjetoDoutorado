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
            #print str(atual.initialPoint) + str(self.points[atual.initialPoint])+ " " + str(atual.finalPoint) + str(self.points[atual.finalPoint])
            print str(atual.initialPoint) + " " + str(atual.finalPoint)
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
            print "Check " + str(atual.root.index) + " " + str(atual.index) + " Angles = " + str(ay) + " " + str(ay2) + " dif " +  str(dif)
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
                pi = self.points[atual.initialPoint]
                pf = self.points[atual.finalPoint]
                pi1 = self.points[atual.son1.initialPoint]
                pf1 = self.points[atual.son1.finalPoint]
                pi2 = self.points[atual.son2.initialPoint]
                pf2 = self.points[atual.son2.finalPoint]
                anglesP = geo.getAngles(pi, pf)
                ayp = geo.radiansToDegree(anglesP[0])
                anglesF1 = geo.getAngles(pi1, pf1)
                ayf1 = geo.radiansToDegree(anglesF1[0])
                anglesF2 = geo.getAngles(pi1, pf2)
                ayf2 = geo.radiansToDegree(anglesF2[0])
                apf1 = ayf1-ayp
                if apf1 > 180:
                    apf1 -= 360
                if apf1 < -180:
                    apf1 += 360
                apf2 = ayf2-ayp
                if apf2 > 180:
                    apf2 -= 360
                if apf2 < -180:
                    apf2 += 360
                dif = apf1-apf2
                dif = abs(dif)
                if dif < 30:
                    lf1 = geo.distance(pi1, pf1)
                    lf2 = geo.distance(pi2, pf2)
                    if lf1 < lf2:
                        positive = True
                        if apf1 < 0:
                            positive = False
                        angleRot = 15
                        if not positive:
                            angleRot = -angleRot
                        geo.rotate(pf1,pi1,angleRot) 
                    else:
                        positive = True
                        if apf2 < 0:
                            positive = False
                        angleRot = 15
                        if not positive:
                            angleRot = -angleRot
                        geo.rotate(pf2,pi2,angleRot) 
                lenght = geo.distance(pi, pf)
                if lenght < atual.radius:
                    print "Erro segmento mto pequeno para bifurcar"
                else:
                    t=(lenght-atual.radius)/lenght
                    pivo = [(1-t)*pi[0]+t*pf[0],((1-t)*pi[1]+t*pf[1]),((1-t)*pi[2]+t*pf[2])]
                    self.points[atual.finalPoint] = pivo
                    geo.createVertex(pivo[0],pivo[1],pivo[2])
                    newPoint1 = list(pf)
                    geo.rotate2(newPoint1, pivo, -25, 2*atual.radius)
                    geo.createVertex(newPoint1[0],newPoint1[1],newPoint1[2])
                    newPoint2 = list(pf)
                    geo.rotate2(newPoint2, pivo, 25, 2*atual.radius)
                    geo.createVertex(newPoint2[0],newPoint2[1],newPoint2[2])
                    d1 = geo.distance(pf1, newPoint1)+geo.distance(pf2, newPoint2)
                    d2 = geo.distance(pf1, newPoint2)+geo.distance(pf2, newPoint1)
                    if d1>d2:
                        aux = newPoint1
                        newPoint1 = newPoint2
                        newPoint2 = aux
                          
                    self.points.append(newPoint1)
                    new1 = Branch(atual.finalPoint, len(self.points)-1, atual.radius*0.75)
                    new1.root = atual
                    new1.son1 = atual.son1
                    new1.son1.initialPoint = len(self.points)-1
                    new1.son1.root = new1
                    atual.son1 = new1
                      
                    self.points.append(newPoint2)
                    new2 = Branch(atual.finalPoint, len(self.points)-1, atual.radius*0.75)
                    new2.root = atual
                    new2.son2 = atual.son2
                    new2.son2.initialPoint = len(self.points)-1
                    new2.son2.root = new2
                    atual.son2 = new2
                              
#                 sonVertex = geo.createVertex(coords[0],coords[1],coords[2])
            if atual.son1 is not None:
                heap.append(atual.son1)
            if atual.son2 is not None:
                heap.append(atual.son2) 
        self.sort() 
        
    def split3(self):
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
                pi = self.points[atual.initialPoint]
                pf = self.points[atual.finalPoint]
                pi1 = self.points[atual.son1.initialPoint]
                pf1 = self.points[atual.son1.finalPoint]
                pi2 = self.points[atual.son2.initialPoint]
                pf2 = self.points[atual.son2.finalPoint]
                anglesP = geo.getAngles(pi, pf)
                ayp = geo.radiansToDegree(anglesP[0])
                anglesF1 = geo.getAngles(pi1, pf1)
                ayf1 = geo.radiansToDegree(anglesF1[0])
                anglesF2 = geo.getAngles(pi1, pf2)
                ayf2 = geo.radiansToDegree(anglesF2[0])
                apf1 = ayf1-ayp
                if apf1 > 180:
                    apf1 -= 360
                if apf1 < -180:
                    apf1 += 360
                apf2 = ayf2-ayp
                if apf2 > 180:
                    apf2 -= 360
                if apf2 < -180:
                    apf2 += 360
                dif = apf1-apf2
                dif = abs(dif)
                if dif < 30:
                    lf1 = geo.distance(pi1, pf1)
                    lf2 = geo.distance(pi2, pf2)
                    if lf1 < lf2:
                        positive = True
                        if apf1 < 0:
                            positive = False
                        angleRot = 15
                        if not positive:
                            angleRot = -angleRot
                        geo.rotate(pf1,pi1,angleRot) 
                    else:
                        positive = True
                        if apf2 < 0:
                            positive = False
                        angleRot = 15
                        if not positive:
                            angleRot = -angleRot
                        geo.rotate(pf2,pi2,angleRot) 
                lenght = geo.distance(pi, pf)
                if lenght < atual.radius:
                    print "Erro segmento mto pequeno para bifurcar"
                else:
                    self.bif(atual) 
#                 sonVertex = geo.createVertex(coords[0],coords[1],coords[2])
            if atual.son1 is not None:
                heap.append(atual.son1)
            if atual.son2 is not None:
                heap.append(atual.son2) 
        self.sort() 
                
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
        while not (len(heap) == 0):
            atual = heap.pop()
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
#                     bodies = []
#                     for s in surfs:
#                         bodies.append(cubit.get_owning_body("surface", s))
#                     geo.imprintBody(bodies) 
#                     geo.mergeBody(bodies) 
#                     vol = geo.createVolume(surfs) 
                    t=(lenght-atual.radius*0.25)/lenght
                    ref = [(1-t)*ipr[0]+t*fpr[0],((1-t)*ipr[1]+t*fpr[1]),((1-t)*ipr[2]+t*fpr[2])]
                    coords = geo2.rotateByAxis(geo2.getVector4_3(ref), geo2.getVector4_3(fpr), u, -90)
                    v1 = geo.createVertex(coords[0],coords[1],coords[2])
                    coords = geo2.rotateByAxis(geo2.getVector4_3(ref), geo2.getVector4_3(fpr), u, 90)
                    v2 = geo.createVertex(coords[0],coords[1],coords[2])
                    self.genSons(atual, arcsf,v1,v2)
#                     geo.imprintVolumes(volsF[0], volsF[1])
#                     geo.mergeVolumes(volsF[0], volsF[1])
#                     volFilhos = geo.unionVolumes(volsF[0], volsF[1])
#                     geo.imprintVolumes(vol, volFilhos)
#                     geo.mergeVolumes(vol, volFilhos)
#                     self.vol = geo.unionVolumes(vol, volFilhos) 
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
#                     bodies = []
#                     for s in surfs:
#                         bodies.append(cubit.get_owning_body("surface", s))
#                     geo.imprintBody(bodies) 
#                     geo.mergeBody(bodies) 
#                     vol = geo.createVolume(surfs) 
                    t=(lenght-atual.radius*0.25)/lenght
                    ref = [(1-t)*ipr[0]+t*fpr[0],((1-t)*ipr[1]+t*fpr[1]),((1-t)*ipr[2]+t*fpr[2])]
                    coords = geo2.rotateByAxis(geo2.getVector4_3(ref), geo2.getVector4_3(fpr), u, -90)
                    v1 = geo.createVertex(coords[0],coords[1],coords[2])
                    coords = geo2.rotateByAxis(geo2.getVector4_3(ref), geo2.getVector4_3(fpr), u, 90)
                    v2 = geo.createVertex(coords[0],coords[1],coords[2])
                    self.genSons(atual, arcsf,v1,v2)
#                     geo.imprintVolumes(volsF[0], volsF[1])
#                     geo.mergeVolumes(volsF[0], volsF[1])
#                     volFilhos = geo.unionVolumes(volsF[0], volsF[1])
#                     geo.imprintVolumes(vol, volFilhos)
#                     geo.mergeVolumes(vol, volFilhos)
#                     vol = geo.unionVolumes(vol, volFilhos) 
#                     geo.imprintVolumes(self.vol, vol)
#                     geo.mergeVolumes(self.vol, vol)
#                     vol = geo.unionVolumes(self.vol, vol)
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
#                     bodies = []
#                     for s in surfs:
#                         bodies.append(cubit.get_owning_body("surface", s))
#                     geo.imprintBody(bodies) 
#                     geo.mergeBody(bodies) 
#                     vol = geo.createVolume(surfs)
#                     geo.imprintVolumes(self.vol, vol)
#                     geo.mergeVolumes(self.vol, vol)
#                     vol = geo.unionVolumes(self.vol, vol) 
                    if atual.son1 is not None:
                        #atual.son1.arcsi = arcsf
                        heap.append(atual.son1) 
#         cubit.cmd('Color Define "%s" RGB %f %f %f' % ("darkred", 0.75,0,0))
#         geo.colorVolume(self.vol, 'user "darkred"')
        self.genVolumes()

    def genVolumes(self):
        heap = [self.root]
        while not (len(heap) == 0):
            atual = heap.pop()
            atual.tube.genVolume()
            if atual.son1 is not None:
                heap.append(atual.son1)
            if atual.son2 is not None:
                heap.append(atual.son2)  
        
    def makeGeometry2(self):
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
        while not (len(heap) == 0):
            atual = heap.pop()
            if atual.root is None:
                circleI = geo.createCircleNormal(atual.initialPoint+1,atual.radius,atual.index)
                circleF = geo.createCircleNormal(atual.finalPoint+1,atual.radius,atual.index)
                surfs = []
                surf = geo.createSurfaceCurve(circleI)
                surfs.append(surf)
                surf = geo.createSkinCurve(circleI,circleF)
                surfs.append(surf)
                surf = geo.createSurfaceCurve(circleF)
                surfs.append(surf)
                geo.mergeAllVertexCurves()
#                 bodies = []
#                 for s in surfs:
#                     bodies.append(cubit.get_owning_body("surface", s))
#                 geo.imprintBody(bodies)  
#                 geo.mergeBody(bodies) 
                vol = geo.createVolume(surfs)
                if atual.son1 is not None and atual.son2 is not None:
                    volf1 = self.genSon(atual.son1)
                    volf2 = self.genSon(atual.son2)
                    geo.imprintVolumes(volf1, volf2)
                    geo.mergeVolumes(volf1, volf2)
                    geo.mergeAllVertexCurves()
                    volFilhos = geo.unionVolumes(volf1, volf2)
                    geo.imprintVolumes(volFilhos, vol)
                    geo.mergeVolumes(volFilhos, vol)
                    geo.mergeAllVertexCurves()
                    self.vol = geo.unionVolumes(volFilhos, vol) 
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
                circleI = geo.createCircleNormal(atual.initialPoint+1,atual.root.radius,atual.root.index)
                circleF = geo.createCircleNormal(atual.finalPoint+1,atual.radius,atual.index)
                surfs = []
                surf = geo.createSurfaceCurve(circleI)
                surfs.append(surf)
                surf = geo.createSkinCurve(circleI,circleF)
                surfs.append(surf)
                surf = geo.createSurfaceCurve(circleF)
                surfs.append(surf)
                geo.mergeAllVertexCurves()
#                 bodies = []
#                 for s in surfs:
#                     bodies.append(cubit.get_owning_body("surface", s))
#                 #geo.imprintBody(bodies)  
#                 #geo.mergeBody(bodies) 
                vol = geo.createVolume(surfs)
                if atual.son1 is not None and atual.son2 is not None:
                    volf1 = self.genSon(atual.son1)
                    volf2 = self.genSon(atual.son2)
                    geo.imprintVolumes(volf1, volf2)
                    geo.mergeVolumes(volf1, volf2)
                    geo.mergeAllVertexCurves()
                    volFilhos = geo.unionVolumes(volf1, volf2)
                    geo.imprintVolumes(volFilhos, vol)
                    geo.mergeVolumes(volFilhos, vol)
                    geo.mergeAllVertexCurves()
                    vol = geo.unionVolumes(volFilhos, vol)
                    geo.imprintVolumes(self.vol, vol)
                    geo.mergeVolumes(self.vol, vol)
                    self.vol = geo.unionVolumes(self.vol, vol)
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
                    geo.imprintVolumes(self.vol, vol)
                    geo.mergeVolumes(self.vol, vol)
                    geo.mergeAllVertexCurves()
                    self.vol = geo.unionVolumes(self.vol, vol)
        cubit.cmd('Color Define "%s" RGB %f %f %f' % ("darkred", 0.75,0,0))
        geo.colorVolume(self.vol, 'user "darkred"')

    
    def genSons(self, atual, arcsi,v1,v2):
        eli1 = geo.createEllipse(arcsi[3][0],v1,atual.finalPoint+1)
        eli2 = geo.createEllipse(arcsi[1][0],v1,atual.finalPoint+1)
        ip = self.points[atual.son1.initialPoint]
        fp = self.points[atual.son1.finalPoint]
        fps = self.points[atual.son1.son1.finalPoint]
        lenght = geo.distance(ip, fp)
        t=(lenght-atual.radius*0.85)/lenght
        ref = [(1-t)*ip[0]+t*fp[0],((1-t)*ip[1]+t*fp[1]),((1-t)*ip[2]+t*fp[2])]
        u0 = [fp[0]-ip[0],fp[1]-ip[1],fp[2]-ip[2]]
        u1 = [fp[0]-fps[0],fp[1]-fps[1],fp[2]-fps[2]]
        u = geo2.crossProduct(u0, u1)
        u = geo2.getVector4_3(u)
        geo2.normalizeVector(u)
        arcsf1 = self.genArc(atual.son1, atual.son1.finalPoint+1,ref,u,atual.radius*0.85)
        atual.son1.tube = Tube([[arcsi[3][0],eli1],[v1,eli2],arcsi[1],arcsi[2]], arcsf1, atual.son1.initialPoint+1, atual.son1.finalPoint+1)
         
        eli3 = geo.createEllipse(arcsi[3][0],v2,atual.finalPoint+1)
        eli4 = geo.createEllipse(arcsi[1][0],v2,atual.finalPoint+1)
        ip = self.points[atual.son2.initialPoint]
        fp = self.points[atual.son2.finalPoint]
        fps = self.points[atual.son2.son1.finalPoint]
        lenght = geo.distance(ip, fp)
        t=(lenght-atual.radius*0.75)/lenght
        ref = [(1-t)*ip[0]+t*fp[0],((1-t)*ip[1]+t*fp[1]),((1-t)*ip[2]+t*fp[2])]
        u0 = [fp[0]-ip[0],fp[1]-ip[1],fp[2]-ip[2]]
        u1 = [fp[0]-fps[0],fp[1]-fps[1],fp[2]-fps[2]]
        u = geo2.crossProduct(u0, u1)
        u = geo2.getVector4_3(u)
        geo2.normalizeVector(u)
        arcsf2 = self.genArc(atual.son2, atual.son2.finalPoint+1,ref,u,atual.radius*0.75)
        atual.son2.tube = Tube([[arcsi[1][0],eli3],[v2,eli4],arcsi[3],arcsi[0]], arcsf2, atual.son2.initialPoint+1, atual.son2.finalPoint+1)
#         bodies = []
#         for s in surfs:
#             bodies.append(cubit.get_owning_body("surface", s))
#         geo.imprintBody(bodies)
#         geo.mergeBody(bodies)
#         volf2 = geo.createVolume(surfs)
#         arcsf2 = [arcsf2[2],arcsf2[3],arcsf2[0],arcsf2[1]]
#         atual.son2.son1.arcsi = arcsf2
#         return [volf1,volf2]
    
    def genSon(self, son):
        initialPoint = self.points[son.root.initialPoint]
        finalPoint = self.points[son.root.finalPoint]
        angles = geo.getAngles(initialPoint, finalPoint)
        ay = geo.radiansToDegree(angles[0]+geo.degreeToRadians(-90))
        az = geo.radiansToDegree(angles[1])
        coords = list(finalPoint)
        geo.rotateAndMove(coords, finalPoint, ay, -90, son.root.radius)
        vAuxRoot1 = geo.createVertex(coords[0],coords[1],coords[2])
        coords = list(finalPoint)
        geo.rotateAndMove(coords, finalPoint, ay, 0, son.root.radius)
        vAuxRoot2 = geo.createVertex(coords[0],coords[1],coords[2])
        coords = list(finalPoint)
        geo.rotateAndMove(coords, finalPoint, ay, 90, son.root.radius)
        vAuxRoot3 = geo.createVertex(coords[0],coords[1],coords[2])
        coords = list(finalPoint)
        geo.rotateAndMove(coords, finalPoint, ay, 180, son.root.radius)
        vAuxRoot4 = geo.createVertex(coords[0],coords[1],coords[2])  
        angle = self.getAngle(son)   
        if angle > 0:
            initialArc1 = geo.createCircleNormal2(son.initialPoint+1,vAuxRoot4,vAuxRoot1,son.root.radius,son.root.index)
            initialArc2 = geo.createCircleNormal2(son.initialPoint+1,vAuxRoot3,vAuxRoot4,son.root.radius,son.root.index)
            coords = list(finalPoint)
            ay2 = geo.radiansToDegree(angles[0]+geo.degreeToRadians(-90))
            geo.rotateAndMove(coords, finalPoint, ay2, 0, son.root.radius/2)
            pivo = geo.createVertex(coords[0],coords[1],coords[2])
            initialLine1 = geo.createEllipse(vAuxRoot3,pivo,son.initialPoint+1)
            initialLine2 = geo.createEllipse(vAuxRoot1,pivo,son.initialPoint+1)
        else:
            initialArc1 = geo.createCircleNormal2(son.initialPoint+1,vAuxRoot1,vAuxRoot2,son.root.radius,son.root.index)
            initialArc2 = geo.createCircleNormal2(son.initialPoint+1,vAuxRoot2,vAuxRoot3,son.root.radius,son.root.index)
            coords = list(finalPoint)
            ay2 = geo.radiansToDegree(angles[0]+geo.degreeToRadians(90))
            geo.rotateAndMove(coords, finalPoint, ay2, 0, son.root.radius/2)
            pivo = geo.createVertex(coords[0],coords[1],coords[2])
            initialLine1 = geo.createEllipse(vAuxRoot3,pivo,son.initialPoint+1)
            initialLine2 = geo.createEllipse(vAuxRoot1,pivo,son.initialPoint+1) 
          
        initialPoint = self.points[son.initialPoint]
        finalPoint = self.points[son.finalPoint]
        angles = geo.getAngles(initialPoint, finalPoint)
        ay = geo.radiansToDegree(angles[0]+geo.degreeToRadians(-90))
        az = geo.radiansToDegree(angles[1])
        coords = list(finalPoint)
        geo.rotateAndMove(coords, finalPoint, ay, -90, son.radius)
        vAuxSon1 = geo.createVertex(coords[0],coords[1],coords[2])
        coords = list(finalPoint)
        geo.rotateAndMove(coords, finalPoint, ay, 0, son.radius)
        vAuxSon2 = geo.createVertex(coords[0],coords[1],coords[2])
        coords = list(finalPoint)
        geo.rotateAndMove(coords, finalPoint, ay, 90, son.radius)
        vAuxSon3 = geo.createVertex(coords[0],coords[1],coords[2])
        coords = list(finalPoint)
        geo.rotateAndMove(coords, finalPoint, ay, 180, son.radius)
        vAuxSon4 = geo.createVertex(coords[0],coords[1],coords[2])
          
        if angle > 0:
            finalArc1 = geo.createCircleNormal2(son.finalPoint+1,vAuxSon4,vAuxSon1,son.radius,son.index)
            finalArc2 = geo.createCircleNormal2(son.finalPoint+1,vAuxSon3,vAuxSon4,son.radius,son.index)
            finalArc3 = geo.createCircleNormal2(son.finalPoint+1,vAuxSon2,vAuxSon3,son.radius,son.index)
            finalArc4 = geo.createCircleNormal2(son.finalPoint+1,vAuxSon1,vAuxSon2,son.radius,son.index)
            horizontalLine1 = geo.createLine(vAuxRoot1,vAuxSon1)  
            horizontalLine2 = geo.createLine(vAuxRoot4,vAuxSon4) 
            horizontalLine3 = geo.createLine(vAuxRoot3,vAuxSon3)
            horizontalLine4 = geo.createLine(pivo,vAuxSon2)
        else:
            finalArc1 = geo.createCircleNormal2(son.finalPoint+1,vAuxSon1,vAuxSon2,son.radius,son.index)
            finalArc2 = geo.createCircleNormal2(son.finalPoint+1,vAuxSon2,vAuxSon3,son.radius,son.index)
            finalArc3 = geo.createCircleNormal2(son.finalPoint+1,vAuxSon3,vAuxSon4,son.radius,son.index)
            finalArc4 = geo.createCircleNormal2(son.finalPoint+1,vAuxSon4,vAuxSon1,son.radius,son.index)
            horizontalLine1 = geo.createLine(vAuxRoot1,vAuxSon1)
            horizontalLine2 = geo.createLine(vAuxRoot2,vAuxSon2)
            horizontalLine3 = geo.createLine(vAuxRoot3,vAuxSon3) 
            horizontalLine4 = geo.createLine(pivo,vAuxSon4) 
 
        curves_s1 = [initialArc1,horizontalLine1,horizontalLine2,finalArc1]
        curves_s2 = [initialArc2,horizontalLine2,horizontalLine3,finalArc2]
        curves_s3 = [initialLine1,horizontalLine3,horizontalLine4,finalArc3]
        curves_s4 = [initialLine2,horizontalLine4,horizontalLine1,finalArc4]
        curves_s5 = [finalArc1,finalArc2,finalArc3,finalArc4]
        curves_s6 = [initialLine1, initialLine2, initialArc1,initialArc2]
        surfs = []
        surfs.append(geo.createSurfaceCurve2(curves_s1))
        surfs.append(geo.createSurfaceCurve2(curves_s2))
        surfs.append(geo.createSurfaceCurve2(curves_s3))
        surfs.append(geo.createSurfaceCurve2(curves_s4))
        surfs.append(geo.createSurfaceCurve2(curves_s5))
        surfs.append(geo.createSurfaceCurve2(curves_s6))
        geo.mergeAllVertexCurves()
#         bodies = []
#         for s in surfs:
#             bodies.append(cubit.get_owning_body("surface", s))
#         geo.imprintBody(bodies)
#         geo.mergeBody(bodies)
        volID = geo.createVolume(surfs)
        return volID

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
        t=(lenght-branch.radius)/lenght
        newPf = [(1-t)*pi[0]+t*pf[0],((1-t)*pi[1]+t*pf[1]),((1-t)*pi[2]+t*pf[2])]
        self.points[branch.finalPoint] = newPf
        geo.createVertex(newPf[0],newPf[1],newPf[2])
        t=(lenght+branch.radius)/lenght
        newPoint1 = [(1-t)*pi[0]+t*pf[0],((1-t)*pi[1]+t*pf[1]),((1-t)*pi[2]+t*pf[2])]
        u0 = [newPf[0]-newPoint1[0],newPf[1]-newPoint1[1],newPf[2]-newPoint1[2]]
        u1 = [newPoint1[0]-pf1[0],newPoint1[1]-pf1[1],newPoint1[2]-pf1[2]]
        u = geo2.crossProduct(u0, u1)
        u = geo2.getVector4_3(u)
        geo2.normalizeVector(u)
        newPoint1 = geo2.rotateByAxis(geo2.getVector4_3(newPoint1), geo2.getVector4_3(newPf), u, 25)
        #geo.rotate2(newPoint1, pivo, -25, 2*atual.radius)
        geo.createVertex(newPoint1[0],newPoint1[1],newPoint1[2])
        t=(lenght+branch.radius)/lenght
        newPoint2 = [(1-t)*pi[0]+t*pf[0],((1-t)*pi[1]+t*pf[1]),((1-t)*pi[2]+t*pf[2])]
        #geo.rotate2(newPoint2, pivo, 25, 2*atual.radius)
        newPoint2 = geo2.rotateByAxis(geo2.getVector4_3(newPoint2), geo2.getVector4_3(newPf), u, -25)
        geo.createVertex(newPoint2[0],newPoint2[1],newPoint2[2])
        d1 = geo.distance(pf1, newPoint1)+geo.distance(pf2, newPoint2)
        d2 = geo.distance(pf1, newPoint2)+geo.distance(pf2, newPoint1)
        if d1>d2:
            aux = newPoint1
            newPoint1 = newPoint2
            newPoint2 = aux
        self.points.append(newPoint1)
        new1 = Branch(branch.finalPoint, len(self.points)-1, branch.radius*0.75)
        new1.root = branch
        new1.son1 = branch.son1
        new1.son1.initialPoint = len(self.points)-1
        new1.son1.root = new1
        branch.son1 = new1
           
        self.points.append(newPoint2)
        new2 = Branch(branch.finalPoint, len(self.points)-1, branch.radius*0.75)
        new2.root = branch
        new2.son2 = branch.son2
        new2.son2.initialPoint = len(self.points)-1
        new2.son2.root = new2
        branch.son2 = new2
        
    def genArcs(self, branch):
        ipr = self.points[branch.initialPoint]
        fpr = self.points[branch.finalPoint]
        ip = self.points[branch.son1.initialPoint]
        fp = self.points[branch.son1.finalPoint]
        lenght = geo.distance(ipr, fpr)
        t=(-branch.radius)/lenght
        ref = [(1-t)*ipr[0]+t*fpr[0],((1-t)*ipr[1]+t*fpr[1]),((1-t)*ipr[2]+t*fpr[2])]
        
        u0 = [fpr[0]-ref[0],fpr[1]-ref[1],fpr[2]-ref[2]]
        u1 = [ref[0]-fp[0],ref[1]-fp[1],ref[2]-fp[2]]
        u = geo2.crossProduct(u0, u1)
        u = geo2.getVector4_3(u)
        geo2.normalizeVector(u)
        
        coords = geo2.rotateByAxis(geo2.getVector4_3(ref), geo2.getVector4_3(ipr), u, -90)
        vAuxI1 = geo.createVertex(coords[0],coords[1],coords[2])
        circleI = geo.createCircleNormal2(branch.initialPoint+1,vAuxI1,vAuxI1,branch.radius,branch.index)
        vAuxI2 = geo.createVertexOnCurveFraction(circleI, 0.25)
        arc = geo.createCircleNormal2(branch.initialPoint+1, vAuxI1, vAuxI2, branch.radius, branch.index)
        arcI1 = [vAuxI1, arc]
        vAuxI3 = geo.createVertexOnCurveFraction(circleI, 0.5)
        arc = geo.createCircleNormal2(branch.initialPoint+1, vAuxI2, vAuxI3, branch.radius, branch.index)
        arcI2 = [vAuxI2, arc]
        vAuxI4 = geo.createVertexOnCurveFraction(circleI, 0.75)
        arc = geo.createCircleNormal2(branch.initialPoint+1, vAuxI3, vAuxI4, branch.radius, branch.index)
        arcI3 = [vAuxI3, arc]
        arc = geo.createCircleNormal2(branch.initialPoint+1, vAuxI4, vAuxI1, branch.radius, branch.index)
        arcI4 = [vAuxI4, arc]
        geo.deleteCurve2(circleI)

        t=(lenght-branch.radius)/lenght
        ref = [(1-t)*ipr[0]+t*fpr[0],((1-t)*ipr[1]+t*fpr[1]),((1-t)*ipr[2]+t*fpr[2])]
        
        u0 = [fpr[0]-ref[0],fpr[1]-ref[1],fpr[2]-ref[2]]
        u1 = [ref[0]-fp[0],ref[1]-fp[1],ref[2]-fp[2]]
        u = geo2.crossProduct(u0, u1)
        u = geo2.getVector4_3(u)
        geo2.normalizeVector(u)
        
        coords = geo2.rotateByAxis(geo2.getVector4_3(ref), geo2.getVector4_3(fpr), u, -90)
        vAuxF1 = geo.createVertex(coords[0],coords[1],coords[2])
        circleF = geo.createCircleNormal2(branch.finalPoint+1,vAuxF1,vAuxF1,branch.radius,branch.index)
        vAuxF2 = geo.createVertexOnCurveFraction(circleF, 0.25)
        arc = geo.createCircleNormal2(branch.finalPoint+1, vAuxF1, vAuxF2, branch.radius, branch.index)
        arcF1 = [vAuxF1, arc]
        vAuxF3 = geo.createVertexOnCurveFraction(circleF, 0.5)
        arc = geo.createCircleNormal2(branch.finalPoint+1, vAuxF2, vAuxF3, branch.radius, branch.index)
        arcF2 = [vAuxF2, arc]
        vAuxF4 = geo.createVertexOnCurveFraction(circleF, 0.75)
        arc = geo.createCircleNormal2(branch.finalPoint+1, vAuxF3, vAuxF4, branch.radius, branch.index)
        arcF3 = [vAuxF3, arc]
        arc = geo.createCircleNormal2(branch.finalPoint+1, vAuxF4, vAuxF1, branch.radius, branch.index)
        arcF4 = [vAuxF4, arc]
        geo.deleteCurve2(circleF)
        arcs = [arcI1,arcI2,arcI3,arcI4,arcF1,arcF2,arcF3,arcF4]
        return arcs
    
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