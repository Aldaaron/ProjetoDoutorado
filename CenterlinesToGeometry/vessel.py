import geo
import cubit
import math
from settings import *

class Vessel(object):

    def __init__(self, index, vI, vF, raio, root):
        self.index = index
        self.vI = vI
        self.vF = vF
        #self.raio = raio
        self.raio = 0.5
        self.nSons = 0
        self.points = []
        self.circles = []
        self.surfs = []
        self.vols = []
        self.root = root
        self.volID = -1
        self.son = False
        self.cf1 = -1
        self.cf2 = -1
        self.son1 = None
        self.son2 = None
        
    def genIntervals(self):
        self.points.append(self.vI)
        dist = cubit.get_distance_between(self.vI, self.vF)
        np = 4
        np += int(dist/(self.raio*2))
        #np = 4
        inc = 1.0/np
        for x in range(1, np):
            vertexID = geo.createVertexOnCurveFraction(self.index, inc*x)
            self.points.append(vertexID)
        self.points.append(self.vF)
    
    def genVolume(self):
        if self.root is not None:
            self.linkToRoot()
        self.genIntervals()
        
        if self.nSons == 2:
            if self.root is None:
                circleID = geo.createCircleNormal(self.vI,self.raio,self.index)
                self.circles.append(circleID)
            circleID = geo.createCircleNormal(self.points[-2],self.raio,self.index)
            self.circles.append(circleID)
            body = []
            surfID = geo.createSkinCurve(self.circles[0],self.circles[-1])
            body.append(surfID)    
            surfID = geo.createSurfaceCurve(self.circles[0])
            body.append(surfID)
            surfID = geo.createSurfaceCurve(self.circles[len(self.circles)-1])
            body.append(surfID)
            self.vols.append(geo.createVolume(body))
            self.vols.append(self.genSon2(-30))
            self.vols.append(self.genSon2(30))
            self.volID = geo.stepUnion(self.vols)
        else: 
            if self.root is None:
                circleID = geo.createCircleNormal(self.vI,self.raio,self.index)
                self.circles.append(circleID)
            circleID = geo.createCircleNormal(self.points[-1],self.raio,self.index)
            self.circles.append(circleID)
            body = []
            surfID = geo.createSkinCurve(self.circles[0],self.circles[-1])
            body.append(surfID)  
            surfID = geo.createSurfaceCurve(self.circles[0])
            body.append(surfID)
            surfID = geo.createSurfaceCurve(self.circles[len(self.circles)-1])
            body.append(surfID)
            self.volID = geo.createVolume(body)
            if self.nSons == 1:
                self.vf1 == self.vF
                self.cf1 == surfID
    
    def linkToRoot(self):
        if self.root.nSons == 2:
            otherSon = self.root.son1
            if otherSon == self:
                otherSon = self.root.son2
            d1 = cubit.get_distance_between(self.vF, self.root.vf1)+cubit.get_distance_between(otherSon.vF, self.root.vf2)
            d2 = cubit.get_distance_between(self.vF, self.root.vf2)+cubit.get_distance_between(otherSon.vF, self.root.vf1)
            if d1<d2:
                self.vI = self.root.vf1
                lineID = geo.createLine(self.root.vf1,self.vF)
                self.index = lineID
                self.circles.append(self.root.cf1)
            else:
                self.vI = self.root.vf2 
                lineID = geo.createLine(self.root.vf2,self.vF)
                self.index = lineID
                self.circles.append(self.root.cf2)
        else:
            self.vI = self.root.vf1
            lineID = geo.createLine(self.root.vf1,self.vF)
            self.index = lineID
            self.circles.append(self.root.cf1)   
    
    def genSon(self, angle):
        initialPoint = cubit.vertex(self.points[-2]).coordinates()
        finalPoint = cubit.vertex(self.points[-1]).coordinates() 
#
        coords = list(initialPoint)
        geo.rotateAndMove(coords, initialPoint, 0, -90, self.raio)
        vAuxRoot = geo.createVertex(coords[0],coords[1],coords[2])
        coords = list(initialPoint)
        geo.rotateAndMove(coords, initialPoint, 0, 90, self.raio)
        vAuxRoot2 = geo.createVertex(coords[0],coords[1],coords[2])
        initialLine = geo.createLine(vAuxRoot,vAuxRoot2)     
        if angle > 0:
            initialArc = geo.createCircleNormal2(self.points[-2],vAuxRoot2,vAuxRoot,self.raio,self.index)
        else:
            initialArc = geo.createCircleNormal2(self.points[-2],vAuxRoot,vAuxRoot2,self.raio,self.index)
#        
        coords = list(finalPoint)
        geo.rotate2(coords, initialPoint, angle, 2*self.raio)
        sonVertex = geo.createVertex(coords[0],coords[1],coords[2])   
        alongLine = geo.createLine(sonVertex,self.points[-2])
       
        sonPoint = cubit.vertex(sonVertex).coordinates()
        coords = list(sonPoint)
        geo.rotateAndMove(coords, sonPoint, 0, -90, 0.75*self.raio)
        vAuxSon = geo.createVertex(coords[0],coords[1],coords[2])
        coords = list(sonPoint)
        geo.rotateAndMove(coords, sonPoint, 0, 90, 0.75*self.raio)
        vAuxSon2 = geo.createVertex(coords[0],coords[1],coords[2])
        if angle > 0:
            finalArc1 = geo.createCircleNormal2(sonVertex,vAuxSon2,vAuxSon,0.75*self.raio,self.index)
            finalArc2 = geo.createCircleNormal2(sonVertex,vAuxSon,vAuxSon2,0.75*self.raio,self.index)
        else:
            finalArc1 = geo.createCircleNormal2(sonVertex,vAuxSon,vAuxSon2,0.75*self.raio,self.index)
            finalArc2 = geo.createCircleNormal2(sonVertex,vAuxSon2,vAuxSon,0.75*self.raio,self.index)
#
        horizontalLine1 = geo.createLine(vAuxRoot,vAuxSon)
        horizontalLine2 = geo.createLine(vAuxRoot2,vAuxSon2)
        if angle < 0:
            self.vf1 = sonVertex
            self.cf1 = geo.createCircleNormal(sonVertex,0.75*self.raio,self.index)
        else:
            self.vf2 = sonVertex
            self.cf2 = geo.createCircleNormal(sonVertex,0.75*self.raio,self.index)      
#
        curves_s1 = [initialArc,horizontalLine1,horizontalLine2,finalArc1]
        curves_s2 = [initialLine,horizontalLine1,horizontalLine2,finalArc2]
        curves_s3 = [finalArc1,finalArc2]
        curves_s4 = [initialLine,initialArc]
        surfs = []
        surfs.append(geo.createSurfaceCurve2(curves_s1))
        surfs.append(geo.createSurfaceCurve2(curves_s2))
        surfs.append(geo.createSurfaceCurve2(curves_s3))
        surfs.append(geo.createSurfaceCurve2(curves_s4))
        return geo.createVolume(surfs)     
    
    def genSon2(self, angle):
        initialPoint = cubit.vertex(self.points[-2]).coordinates()
        finalPoint = cubit.vertex(self.points[-1]).coordinates() 
        angles = geo.getAngles(initialPoint, finalPoint)
        ay = geo.radiansToDegree(angles[0]+geo.degreeToRadians(-90))
        az = geo.radiansToDegree(angles[1])
        coords = list(initialPoint)
        geo.rotateAndMove(coords, initialPoint, ay, -90, self.raio)
        vAuxRoot1 = geo.createVertex(coords[0],coords[1],coords[2])
        coords = list(initialPoint)
        geo.rotateAndMove(coords, initialPoint, ay, 0, self.raio)
        vAuxRoot2 = geo.createVertex(coords[0],coords[1],coords[2])
        coords = list(initialPoint)
        geo.rotateAndMove(coords, initialPoint, ay, 90, self.raio)
        vAuxRoot3 = geo.createVertex(coords[0],coords[1],coords[2])
        coords = list(initialPoint)
        geo.rotateAndMove(coords, initialPoint, ay, 180, self.raio)
        vAuxRoot4 = geo.createVertex(coords[0],coords[1],coords[2])     
        if angle > 0:
            initialArc1 = geo.createCircleNormal2(self.points[-2],vAuxRoot4,vAuxRoot1,self.raio,self.index)
            initialArc2 = geo.createCircleNormal2(self.points[-2],vAuxRoot3,vAuxRoot4,self.raio,self.index)
            coords = list(initialPoint)
            ay2 = geo.radiansToDegree(angles[0]+geo.degreeToRadians(-90))
            geo.rotateAndMove(coords, initialPoint, ay2, 0, self.raio/2)
            pivo = geo.createVertex(coords[0],coords[1],coords[2])
            initialLine1 = geo.createEllipse(vAuxRoot3,pivo,self.points[-2])
            initialLine2 = geo.createEllipse(vAuxRoot1,pivo,self.points[-2])
            #initialLine1 = geo.createLine(vAuxRoot3,self.points[-2])
            #initialLine2 = geo.createLine(self.points[-2],vAuxRoot1)   
        else:
            initialArc1 = geo.createCircleNormal2(self.points[-2],vAuxRoot1,vAuxRoot2,self.raio,self.index)
            initialArc2 = geo.createCircleNormal2(self.points[-2],vAuxRoot2,vAuxRoot3,self.raio,self.index)
            coords = list(initialPoint)
            ay2 = geo.radiansToDegree(angles[0]+geo.degreeToRadians(90))
            geo.rotateAndMove(coords, initialPoint, ay2, 0, self.raio/2)
            pivo = geo.createVertex(coords[0],coords[1],coords[2])
            initialLine1 = geo.createEllipse(vAuxRoot3,pivo,self.points[-2])
            initialLine2 = geo.createEllipse(vAuxRoot1,pivo,self.points[-2])
            #initialLine1 = geo.createLine(self.points[-2],vAuxRoot3)
            #initialLine2 = geo.createLine(vAuxRoot1,self.points[-2])   
        
        coords = list(finalPoint)
        geo.rotate2(coords, initialPoint, angle, 2*self.raio)
        sonVertex = geo.createVertex(coords[0],coords[1],coords[2])   
        alongLine = geo.createLine(sonVertex,self.points[-2])
       
        sonPoint = cubit.vertex(sonVertex).coordinates()
        coords = list(sonPoint)
        geo.rotateAndMove(coords, sonPoint, ay, -90, 0.75*self.raio)
        vAuxSon1 = geo.createVertex(coords[0],coords[1],coords[2])
        coords = list(sonPoint)
        geo.rotateAndMove(coords, sonPoint, ay, 0, 0.75*self.raio)
        vAuxSon2 = geo.createVertex(coords[0],coords[1],coords[2])
        coords = list(sonPoint)
        geo.rotateAndMove(coords, sonPoint, ay, 90, 0.75*self.raio)
        vAuxSon3 = geo.createVertex(coords[0],coords[1],coords[2])
        coords = list(sonPoint)
        geo.rotateAndMove(coords, sonPoint, ay, 180, 0.75*self.raio)
        vAuxSon4 = geo.createVertex(coords[0],coords[1],coords[2])
        
        if angle > 0:
            finalArc1 = geo.createCircleNormal2(sonVertex,vAuxSon4,vAuxSon1,0.75*self.raio,self.index)
            finalArc2 = geo.createCircleNormal2(sonVertex,vAuxSon3,vAuxSon4,0.75*self.raio,self.index)
            finalArc3 = geo.createCircleNormal2(sonVertex,vAuxSon2,vAuxSon3,0.75*self.raio,self.index)
            finalArc4 = geo.createCircleNormal2(sonVertex,vAuxSon1,vAuxSon2,0.75*self.raio,self.index)
            horizontalLine1 = geo.createLine(vAuxRoot1,vAuxSon1)  
            horizontalLine2 = geo.createLine(vAuxRoot4,vAuxSon4) 
            horizontalLine3 = geo.createLine(vAuxRoot3,vAuxSon3)
            horizontalLine4 = geo.createLine(pivo,vAuxSon2)
        else:
            finalArc1 = geo.createCircleNormal2(sonVertex,vAuxSon1,vAuxSon2,0.75*self.raio,self.index)
            finalArc2 = geo.createCircleNormal2(sonVertex,vAuxSon2,vAuxSon3,0.75*self.raio,self.index)
            finalArc3 = geo.createCircleNormal2(sonVertex,vAuxSon3,vAuxSon4,0.75*self.raio,self.index)
            finalArc4 = geo.createCircleNormal2(sonVertex,vAuxSon4,vAuxSon1,0.75*self.raio,self.index)
            horizontalLine1 = geo.createLine(vAuxRoot1,vAuxSon1)
            horizontalLine2 = geo.createLine(vAuxRoot2,vAuxSon2)
            horizontalLine3 = geo.createLine(vAuxRoot3,vAuxSon3)
            horizontalLine4 = geo.createLine(pivo,vAuxSon4) 
        if angle < 0:
            self.vf1 = sonVertex
            self.cf1 = geo.createCircleNormal(sonVertex,0.75*self.raio,self.index)
        else:
            self.vf2 = sonVertex
            self.cf2 = geo.createCircleNormal(sonVertex,0.75*self.raio,self.index)      

        curves_s1 = [initialArc1,horizontalLine1,horizontalLine2,finalArc1]
        curves_s2 = [initialArc2,horizontalLine2,horizontalLine3,finalArc2]
        curves_s3 = [initialLine1,horizontalLine3,horizontalLine4,finalArc3]
        curves_s4 = [initialLine2,horizontalLine4,horizontalLine1,finalArc4]
        curves_s5 = [finalArc1,finalArc2,finalArc3,finalArc4]
        curves_s6 = [initialLine1,initialLine2,initialArc1,initialArc2]
        surfs = []
        surfs.append(geo.createSurfaceCurve2(curves_s1))
        surfs.append(geo.createSurfaceCurve2(curves_s2))
        surfs.append(geo.createSurfaceCurve2(curves_s3))
        surfs.append(geo.createSurfaceCurve2(curves_s4))
        surfs.append(geo.createSurfaceCurve2(curves_s5))
        surfs.append(geo.createSurfaceCurve2(curves_s6))
        return geo.createVolume(surfs)
        
    def clean(self):
        for p in self.points:
            if p != self.vI and p != self.vF:
                geo.deleteVertex(p)
        geo.deleteCurve(self.index)  
        for c in self.circles:
            geo.deleteCurve(c)  
        for s in self.surfs:
            geo.deleteSurface(s)  
        #geo.deleteVolume(self.volID)
            