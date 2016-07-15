import geo
import cubit
import math
from settings import *
from __builtin__ import range

class Vessel(object):

    def __init__(self, index, vI, vF, raio, root):
        self.index = index
        self.vI = vI
        self.vF = vF
        self.raio = raio
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
        self.vf1 = -1
        self.vf2 = -1
        self.son1 = None
        self.son2 = None
        self.splitLines = None
        self.splitPoints = None
        self.linkCirclesf1 = None
        self.linkCirclesf2 = None
        self.volf1 = -1
        self.volf2 = -1
        
    def genIntervals(self):
        self.points.append(self.vI)
        dist = cubit.get_distance_between(self.vI, self.vF)
        np = 4
        np += int(dist/(self.raio*10))
        inc = 1.0/np
        for x in range(1, np):
            vertexID = geo.createVertexOnCurveFraction(self.index, inc*x)
            self.points.append(vertexID)
        self.points.append(self.vF)
            
    def linkToRoot2(self):
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
            
    def genSurfaces(self):
        a = 1
        
    
    def genSurfaces2(self):
        if self.root is not None:
            self.linkToRoot2()
        self.genIntervals()
        if self.nSons == 2:
            if self.root is None:
                for x in range(0, len(self.points)-2):
                    circleID = geo.createCircleNormal(self.points[x],self.raio,self.index)
                    self.circles.append(circleID)
            else:
                for x in range(1, len(self.points)-2):
                    raio = (x*self.raio + (len(self.points)-2-x) * self.root.raio*0.75) / (len(self.points)-2)
                    circleID = geo.createCircleNormal(self.points[x],raio,self.index)
                    self.circles.append(circleID)
            linksF1 = self.genSon2(-25)
            linksF2 = self.genSon2(25)
            linkCircles = []
            linkCircles.extend(linksF1)
            linkCircles.extend(linksF2)
            circleID = geo.createCombineCurve(linkCircles)
            self.circles.append(circleID)
            geo.imprintVolumes(self.volf1, self.volf2)
            geo.mergeVolumes(self.volf1, self.volf2)
            volFilhos = geo.unionVolumes(self.volf1, self.volf2)

            surfs = []
            surfID = geo.createSurfaceCurve(self.circles[0])
            surfs.append(surfID)
            for x in range(0, len(self.circles)-1):
                surfID = geo.createSkinCurve(self.circles[x],self.circles[x+1])
                surfs.append(surfID)  
            surfID = geo.createSurfaceCurve(self.circles[-1])
            surfs.append(surfID)
            bodies = []
            for s in surfs:
                bodies.append(cubit.get_owning_body("surface", s))
            geo.imprintBody(bodies)  
            geo.mergeBody(bodies) 
            volID = geo.createVolume(surfs)
            geo.imprintVolumes(volFilhos, volID)
            geo.mergeVolumes(volFilhos, volID)
            self.volID = geo.unionVolumes(volFilhos, volID) 
        else: 
            if self.root is None:
                for x in range(0, len(self.points)):
                    circleID = geo.createCircleNormal(self.points[x],self.raio,self.index)
                    self.circles.append(circleID)
            else:
                for x in range(1, len(self.points)):
                    raio = (x*self.raio + (len(self.points)-1-x) * self.root.raio*0.75) / (len(self.points)-1)
                    circleID = geo.createCircleNormal(self.points[x],raio,self.index)
                    self.circles.append(circleID)
            surfs = []
            surfID = geo.createSurfaceCurve(self.circles[0])
            surfs.append(surfID)
            for x in range(0, len(self.circles)-1):
                surfID = geo.createSkinCurve(self.circles[x],self.circles[x+1])
                surfs.append(surfID)  
            surfID = geo.createSurfaceCurve(self.circles[-1])
            surfs.append(surfID)
            bodies = []
            for s in surfs:
                bodies.append(cubit.get_owning_body("surface", s))
            geo.imprintBody(bodies)  
            geo.mergeBody(bodies) 
            volID = geo.createVolume(surfs)
            if self.nSons == 1:
                self.vf1 = self.vF
                self.cf1 = self.circles[-1]
    
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
            #initialArc1 = geo.createEllipse(vAuxRoot4,vAuxRoot1,self.points[-2])
            #initialArc2 = geo.createEllipse(vAuxRoot3,vAuxRoot4,self.points[-2])
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
            #initialArc1 = geo.createEllipse(vAuxRoot1,vAuxRoot2,self.points[-2])
            #initialArc2 = geo.createEllipse(vAuxRoot2,vAuxRoot3,self.points[-2])
            print "Arc = " + str(initialArc1)  + " Points " + str(vAuxRoot1) + " - "+str(vAuxRoot2)
            print "Arc = " + str(initialArc2)  + " Points " + str(vAuxRoot2) + " - "+str(vAuxRoot3)
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
        finalPoint = cubit.vertex(sonVertex).coordinates() 
        angles = geo.getAngles(initialPoint, finalPoint)
        ay = geo.radiansToDegree(angles[0]+geo.degreeToRadians(-90))
        az = geo.radiansToDegree(angles[1])
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
            finalArc1 = geo.createCircleNormal2(sonVertex,vAuxSon1,vAuxSon4,0.75*self.raio,alongLine)
            finalArc2 = geo.createCircleNormal2(sonVertex,vAuxSon4,vAuxSon3,0.75*self.raio,alongLine)
            finalArc3 = geo.createCircleNormal2(sonVertex,vAuxSon3,vAuxSon2,0.75*self.raio,alongLine)
            finalArc4 = geo.createCircleNormal2(sonVertex,vAuxSon2,vAuxSon1,0.75*self.raio,alongLine)
            horizontalLine1 = geo.createLine(vAuxRoot1,vAuxSon1)  
            horizontalLine2 = geo.createLine(vAuxRoot4,vAuxSon4) 
            horizontalLine3 = geo.createLine(vAuxRoot3,vAuxSon3)
#             horizontalLine4 = geo.createLine(self.points[-2],vAuxSon2)
            horizontalLine4 = geo.createLine(pivo,vAuxSon2)
        else:
            finalArc1 = geo.createCircleNormal2(sonVertex,vAuxSon2,vAuxSon1,0.75*self.raio,alongLine)
            finalArc2 = geo.createCircleNormal2(sonVertex,vAuxSon3,vAuxSon2,0.75*self.raio,alongLine)
            finalArc3 = geo.createCircleNormal2(sonVertex,vAuxSon4,vAuxSon3,0.75*self.raio,alongLine)
            finalArc4 = geo.createCircleNormal2(sonVertex,vAuxSon1,vAuxSon4,0.75*self.raio,alongLine)
            horizontalLine1 = geo.createLine(vAuxRoot1,vAuxSon1)
            horizontalLine2 = geo.createLine(vAuxRoot2,vAuxSon2)
            horizontalLine3 = geo.createLine(vAuxRoot3,vAuxSon3)
#             horizontalLine4 = geo.createLine(self.points[-2],vAuxSon4) 
            horizontalLine4 = geo.createLine(pivo,vAuxSon4) 

        #curves_s1 = [initialArc2,initialLine1,horizontalLine2,horizontalLine4,finalArc2, finalArc3]
        #curves_s2 = [initialArc1,initialLine2,horizontalLine2,horizontalLine4,finalArc1, finalArc4]
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
        bodies = []
        for s in surfs:
            bodies.append(cubit.get_owning_body("surface", s))
        geo.imprintBody(bodies)
        geo.mergeBody(bodies)
        volIDf = geo.createVolume(surfs)
        if angle < 0:
            self.vf1 = sonVertex
            self.cf1 = geo.createCircleNormal(sonVertex,0.75*self.raio,alongLine)
            self.volf1 = volIDf
        else:
            self.vf2 = sonVertex
            self.cf2 = geo.createCircleNormal(sonVertex,0.75*self.raio,alongLine)      
            self.volf2 = volIDf
        link = [initialArc1, initialArc2]
        return link
    
    def genVolume(self):
        geo.imprintBody(self.surfs)
        geo.mergeBody(self.surfs)
        self.volID = geo.createVolume(self.surfs)
        
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
            