import geo
from vessel import Vessel
import cubit
import math
import random
from settings import *
from branch import Branch
from geo import distance

class Arvore(object):
    
    def __init__(self, points, links):
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
                
    def draw(self):
        for p in self.points:
            geo.createVertex(p[0]*scale,p[1]*scale,p[2]*scale)
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
                lenght = geo.distance(pi, pf)
                if lenght < atual.radius:
                    print "Erro segmento mto pequeno para bifurcar"
                else:
                    t=(lenght-atual.radius)/lenght
                    print t
                    pivo = [(1-t)*pi[0]+t*pf[0],((1-t)*pi[1]+t*pf[1]),((1-t)*pi[2]+t*pf[2])]
                    self.points[atual.finalPoint] = pivo
                    geo.createVertex(pivo[0]*scale,pivo[1]*scale,pivo[2]*scale)
                    newPoint1 = list(pf)
                    geo.rotate2(newPoint1, pivo, -25, 2*atual.radius)
                    geo.createVertex(newPoint1[0]*scale,newPoint1[1]*scale,newPoint1[2]*scale)
                    newPoint2 = list(pf)
                    geo.rotate2(newPoint2, pivo, 25, 2*atual.radius)
                    geo.createVertex(newPoint2[0]*scale,newPoint2[1]*scale,newPoint2[2]*scale)
                    d1 = geo.distance(pf1, newPoint1)+geo.distance(pf2, newPoint2)
                    d2 = geo.distance(pf1, newPoint2)+geo.distance(pf2, newPoint1)
                    if d1>d2:
                        aux = newPoint1
                        newPoint1 = newPoint2
                        newPoint2 = aux
                          
                    self.points.append(newPoint1)
                    new1 = Branch(atual.finalPoint, len(self.points)-1, atual.radius)
                    new1.root = atual
                    new1.son1 = atual.son1
                    new1.son1.initialPoint = len(self.points)-1
                    new1.son1.root = new1
                    atual.son1 = new1
                      
                    self.points.append(newPoint2)
                    new2 = Branch(atual.finalPoint, len(self.points)-1, atual.radius)
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
                
    def smoth(self):
        change = True
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
                if dif > 40:
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
                    vertexCoords = [(initial[0]+final[0])/2,(initial[1]+final[1])/2,(initial[2]+final[2])/2]
                    vertexCoords[0] = vertexCoords[0]
                    vertexCoords[1] = vertexCoords[1]
                    vertexCoords[2] = vertexCoords[2]
                    geo.rotate(vertexCoords,initial,angleRot) 
                    self.points.append(vertexCoords)
                    new = Branch(len(self.points)-1, atual.finalPoint, atual.radius)
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