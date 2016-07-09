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
            for x in range(0, len(self.points)-1):
                if x != 0 or self.root is None:
                    raio = self.raio
                    if self.root is not None:
                        raio = (x*self.raio + (len(self.points)-1-x) * self.root.raio*0.75) / (len(self.points)-1)
                    circleID = geo.createCircleNormal(self.points[x],raio,self.index)
                    self.circles.append(circleID)
            body1 = []
            for x in range(0, len(self.circles)-1):
                surfID = geo.createSkinCurve(self.circles[x],self.circles[x+1])
                body1.append(surfID)    
            surfID = geo.createSurfaceCurve(self.circles[0])
            body1.append(surfID)
            surfID = geo.createSurfaceCurve(self.circles[len(self.circles)-1])
            body1.append(surfID)
            self.vols.append(geo.createVolume(body1))
            
            v1 = cubit.vertex(self.points[0]).coordinates()
            v2 = cubit.vertex(self.points[-1]).coordinates()
            x = v1[0]
            y = v1[1]
            z = v1[2]
            x2 = v2[0]
            y2 = v2[1]
            z2 = v2[2]
            d = cubit.get_distance_between(self.points[0], self.points[-2])
            dt = d + self.raio*0.1
            t = dt/d
            
            vAux = geo.createVertex(x,y,z)
            x3 = x - t*x +t*x2;
            y3 = y - t*y +t*y2;
            z3 = z - t*z +t*z2;
            geo.moveVertexL(vAux,x3,y3,z3)
            
            self.branch(vAux)
            
            circleID = geo.createCircleNormal(vAux,0.25*self.raio,self.index)  
            body2 = []
            surfID = geo.createSkinCurve(self.circles[-1],circleID)
            body2.append(surfID)  
            surfID = geo.createSurfaceCurve(self.circles[-1])
            body2.append(surfID)
            surfID = geo.createSurfaceCurve(circleID)
            body2.append(surfID)
            self.vols.append(geo.createVolume(body2))
            self.volID = geo.stepUnion(self.vols)
        else: 
            for x in range(0, len(self.points)):
                if x != 0 or self.root is None:
                    raio = self.raio
                    if self.root is not None:
                        raio = (x*self.raio + (len(self.points)-1-x) * self.root.raio*0.75) / (len(self.points)-1)
                    circleID = geo.createCircleNormal(self.points[x],raio,self.index)
                    self.circles.append(circleID)
            body = []
            for x in range(0, len(self.circles)-1):
                surfID = geo.createSkinCurve(self.circles[x],self.circles[x+1])
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
    
    def branch(self, vID):
        v1 = cubit.vertex(self.points[-2]).coordinates()
        v2 = cubit.vertex(vID).coordinates()
        x = v1[0]
        y = v1[1]
        z = v1[2]
        x2 = v2[0]
        y2 = v2[1]
        z2 = v2[2]
        v1ID = geo.createVertex(x,y,z)
        deltaY = y2-y;
        deltaX = x2-x;
        deltaZ = z2-z;
        ay = math.atan2(deltaY, deltaX)
        az = math.atan2(deltaZ, deltaX)
        angle = (90) * (math.pi/180)
        xT = math.cos(angle) * (x2 - x) - math.sin(angle) * (y2 - y) + x;
        yT = math.sin(angle) * (x2 - x) + math.cos(angle) * (y2 - y) + y;
        vAux = geo.createVertex(x,y,z)
        x3 = x + self.raio * math.cos(0*0.017453292519+az)*math.cos(-90*0.017453292519+ay);
        y3 = y + self.raio * math.cos(0*0.017453292519+az)*math.sin(-90*0.017453292519+ay);
        z3 = z + self.raio * math.sin(0*0.017453292519+az);
        geo.moveVertexL(vAux,x3,y3,z3)
        angle = -angle
        xT = math.cos(angle) * (x2 - x) - math.sin(angle) * (y2 - y) + x;
        yT = math.sin(angle) * (x2 - x) + math.cos(angle) * (y2 - y) + y;
        vAux2 = geo.createVertex(x,y,z)
        x3 = x + self.raio * math.cos(0*0.017453292519+az)*math.cos(90*0.017453292519+ay);
        y3 = y + self.raio * math.cos(0*0.017453292519+az)*math.sin(90*0.017453292519+ay);
        z3 = z + self.raio * math.sin(0*0.017453292519+az);
        geo.moveVertexL(vAux2,x3,y3,z3)
        l1 = geo.createLine(vAux,self.points[-2])
        l2 = geo.createLine(vAux2,self.points[-2])
        vt1 = geo.createVertexOnCurveDistance(l1,0.75*self.raio)
        vt2 = geo.createVertexOnCurveDistance(l2,0.75*self.raio)
        circleIDf1i = geo.createCircleNormal(vt1,0.75*self.raio,self.index)
        circleIDf2i = geo.createCircleNormal(vt2,0.75*self.raio,self.index)      
        angle = (-30) * (math.pi/180)
        self.vf1 = geo.createVertex(0,0,0)
        coordsf1 = []
        coordsf1.append(0)
        coordsf1.append(0)
        coordsf1.append(0)
        origin = []
        origin.append(v1[0])
        origin.append(v1[1])
        origin.append(v1[2])
        geo.spheCoords(origin)
        coordsf1[0] += 2*(self.raio)
        angleR = geo.degreeToRadians(30)
        coordsf1[1] += geo.degreeToRadians(90)
        coordsf1[2] += ay-angleR
        geo.cartCoords(coordsf1)
        coordsf1[0] += x
        coordsf1[1] += y
        coordsf1[2] += z
        geo.moveVertexL(self.vf1,coordsf1[0],coordsf1[1],coordsf1[2])
        self.vf2 = geo.createVertex(0,0,0)
        coordsf2 = []
        coordsf2.append(0)
        coordsf2.append(0)
        coordsf2.append(0)
        coordsf2[0] += 2*(self.raio)
        coordsf2[1] += geo.degreeToRadians(90)
        coordsf2[2] += ay+angleR
        geo.cartCoords(coordsf2)
        coordsf2[0] += x
        coordsf2[1] += y
        coordsf2[2] += z
        geo.moveVertexL(self.vf2,coordsf2[0],coordsf2[1],coordsf2[2])
        l1 = geo.createLine(self.vf1,self.points[-2])
        l2 = geo.createLine(self.vf2,self.points[-2])
        circleIDf1f = geo.createCircleNormal(self.vf1,0.75*self.raio,self.index)
        circleIDf2f = geo.createCircleNormal(self.vf2,0.75*self.raio,self.index)
        self.cf1 = circleIDf1f
        self.cf2 = circleIDf2f   
        bodyf1 = []
        surfID = geo.createSkinCurve(circleIDf1i,circleIDf1f)
        bodyf1.append(surfID)  
        surfID = geo.createSurfaceCurve(circleIDf1i)
        bodyf1.append(surfID)
        surfID = geo.createSurfaceCurve(circleIDf1f)
        bodyf1.append(surfID)
        bodyf2 = []
        surfID = geo.createSkinCurve(circleIDf2i,circleIDf2f)
        bodyf2.append(surfID)  
        surfID = geo.createSurfaceCurve(circleIDf2i)
        bodyf2.append(surfID)
        surfID = geo.createSurfaceCurve(circleIDf2f)
        bodyf2.append(surfID)   
        self.vols.append(geo.createVolume(bodyf1))
        self.vols.append(geo.createVolume(bodyf2))
    
    def genVolume2(self):
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
            self.branch2()
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
    
    def branch2(self):
        v1 = cubit.vertex(self.points[-2]).coordinates()
        v2 = cubit.vertex(self.points[-1]).coordinates()
        x = v1[0]
        y = v1[1]
        z = v1[2]
        x2 = v2[0]
        y2 = v2[1]
        z2 = v2[2]
        v1ID = geo.createVertex(x,y,z)
        deltaY = y2-y;
        deltaX = x2-x;
        deltaZ = z2-z;
        ay = math.atan2(deltaY, deltaX)
        az = math.atan2(deltaZ, deltaX)
        angle = (-90) * (math.pi/180)
        vAux = geo.createVertex(x,y,z)
        x3 = x + self.raio * math.cos(angle+az)*math.cos(ay);
        y3 = y + self.raio * math.cos(angle+az)*math.sin(ay);
        z3 = z + self.raio * math.sin(angle+az);
        geo.moveVertexL(vAux,x3,y3,z3)
        angle = -angle
        vAux2 = geo.createVertex(x,y,z)
        x3 = x + self.raio * math.cos(angle+az)*math.cos(ay);
        y3 = y + self.raio * math.cos(angle+az)*math.sin(ay);
        z3 = z + self.raio * math.sin(angle+az);
        geo.moveVertexL(vAux2,x3,y3,z3)
        lp = geo.createLine(vAux,vAux2)
        cp1 = geo.createCircleNormal2(self.points[-2],vAux2,vAux,self.raio,self.index)
        cp2 = geo.createCircleNormal2(self.points[-2],vAux,vAux2,self.raio,self.index)
        #vt1 = geo.createVertexOnCurveDistance(l1,0.75*self.raio)
        #vt2 = geo.createVertexOnCurveDistance(l2,0.75*self.raio)
        #circleIDf1i = geo.createCircleNormal(vt1,0.75*self.raio,self.index)
        #circleIDf2i = geo.createCircleNormal(vt2,0.75*self.raio,self.index)      
        angle = (-30) * (math.pi/180)
        self.vf1 = geo.createVertex(0,0,0)
        coordsf1 = []
        coordsf1.append(0)
        coordsf1.append(0)
        coordsf1.append(0)
        origin = []
        origin.append(v1[0])
        origin.append(v1[1])
        origin.append(v1[2])
        geo.spheCoords(origin)
        coordsf1[0] += 2*(self.raio)
        angleR = geo.degreeToRadians(30)
        coordsf1[1] += geo.degreeToRadians(90)
        coordsf1[2] += ay-angleR
        geo.cartCoords(coordsf1)
        coordsf1[0] += x
        coordsf1[1] += y
        coordsf1[2] += z
        geo.moveVertexL(self.vf1,coordsf1[0],coordsf1[1],coordsf1[2])
        self.vf2 = geo.createVertex(0,0,0)
        coordsf2 = []
        coordsf2.append(0)
        coordsf2.append(0)
        coordsf2.append(0)
        coordsf2[0] += 2*(self.raio)
        coordsf2[1] += geo.degreeToRadians(90)
        coordsf2[2] += ay+angleR
        geo.cartCoords(coordsf2)
        coordsf2[0] += x
        coordsf2[1] += y
        coordsf2[2] += z
        geo.moveVertexL(self.vf2,coordsf2[0],coordsf2[1],coordsf2[2])
        l1 = geo.createLine(self.vf1,self.points[-2])
        l2 = geo.createLine(self.vf2,self.points[-2])
        circleIDf1f = geo.createCircleNormal(self.vf1,0.75*self.raio,self.index)
        circleIDf2f = geo.createCircleNormal(self.vf2,0.75*self.raio,self.index)
        self.cf1 = circleIDf1f
        self.cf2 = circleIDf2f
        angle = (-90) * (math.pi/180)
        vAuxf11 = geo.createVertex(0,0,0)
        f1c = cubit.vertex(self.vf1).coordinates()
        x3 = f1c[0] + 0.75*self.raio * math.cos(angle+az)*math.cos(ay);
        y3 = f1c[1] + 0.75*self.raio * math.cos(angle+az)*math.sin(ay);
        z3 = f1c[2] + 0.75*self.raio * math.sin(angle+az);
        geo.moveVertexL(vAuxf11,x3,y3,z3)
        angle = -angle
        vAuxf12 = geo.createVertex(0,0,0)
        x3 = f1c[0] + 0.75*self.raio * math.cos(angle+az)*math.cos(ay);
        y3 = f1c[1] + 0.75*self.raio * math.cos(angle+az)*math.sin(ay);
        z3 = f1c[2] + 0.75*self.raio * math.sin(angle+az);
        geo.moveVertexL(vAuxf12,x3,y3,z3)
        lf1 = geo.createLine(vAuxf11,vAuxf12)
        cf1_1 = geo.createCircleNormal2(self.vf1,vAuxf12,vAuxf11,0.75*self.raio,self.index)
        cf1_2 = geo.createCircleNormal2(self.vf1,vAuxf11,vAuxf12,0.75*self.raio,self.index)
        ldirf1_1 = geo.createLine(vAux,vAuxf11)
        ldirf1_2 = geo.createLine(vAux2,vAuxf12)
        angle = (-90) * (math.pi/180)
        vAuxf21 = geo.createVertex(0,0,0)
        f2c = cubit.vertex(self.vf2).coordinates()
        x3 = f2c[0] + 0.75*self.raio * math.cos(angle+az)*math.cos(ay);
        y3 = f2c[1] + 0.75*self.raio * math.cos(angle+az)*math.sin(ay);
        z3 = f2c[2] + 0.75*self.raio * math.sin(angle+az);
        geo.moveVertexL(vAuxf21,x3,y3,z3)
        angle = -angle
        vAuxf22 = geo.createVertex(0,0,0)
        x3 = f2c[0] + 0.75*self.raio * math.cos(angle+az)*math.cos(ay);
        y3 = f2c[1] + 0.75*self.raio * math.cos(angle+az)*math.sin(ay);
        z3 = f2c[2] + 0.75*self.raio * math.sin(angle+az);
        geo.moveVertexL(vAuxf22,x3,y3,z3)
        lf2 = geo.createLine(vAuxf21,vAuxf22)
        cf2_1 = geo.createCircleNormal2(self.vf2,vAuxf22,vAuxf21,0.75*self.raio,self.index)
        cf2_2 = geo.createCircleNormal2(self.vf2,vAuxf21,vAuxf22,0.75*self.raio,self.index)
        ldirf2_1 = geo.createLine(vAux,vAuxf21)
        ldirf2_2 = geo.createLine(vAux2,vAuxf22)
        #building filho1 -----------------------------------------------------------------
        curvesf1_s1 = [cp1,ldirf1_1,ldirf1_2,cf1_1]
        curvesf1_s2 = [lp,ldirf1_1,ldirf1_2,cf1_2]
        curvesf1_s3 = [cf1_1,cf1_2]
        curvesf1_s4 = [lp,cp1]
        surfsf1 = []
        surfsf1.append(geo.createSurfaceCurve2(curvesf1_s1))
        surfsf1.append(geo.createSurfaceCurve2(curvesf1_s2))
        surfsf1.append(geo.createSurfaceCurve2(curvesf1_s3))
        surfsf1.append(geo.createSurfaceCurve2(curvesf1_s4))
        volIDf1 = geo.createVolume(surfsf1)
        #building filho2 -----------------------------------------------------------------
        curvesf2_s1 = [lp,ldirf2_1,ldirf2_2,cf2_1]
        curvesf2_s2 = [cp2,ldirf2_1,ldirf2_2,cf2_2]
        curvesf2_s3 = [cf2_1,cf2_2]
        curvesf2_s4 = [lp,cp2]
        surfsf2 = []
        surfsf2.append(geo.createSurfaceCurve2(curvesf2_s1))
        surfsf2.append(geo.createSurfaceCurve2(curvesf2_s2))
        surfsf2.append(geo.createSurfaceCurve2(curvesf2_s3))
        surfsf2.append(geo.createSurfaceCurve2(curvesf2_s4))
        volIDf2 = geo.createVolume(surfsf2)
        vols = [volIDf1,volIDf2]
        self.vols.append(geo.stepUnion(vols))
    
    def genSon(self, angle):
        initialPoint = cubit.vertex(self.points[-2]).coordinates()
        finalPoint = cubit.vertex(self.points[-1]).coordinates()
        angles = geo.getAngles(initialPoint,finalPoint)  
#        
        coords = [initialPoint[0], initialPoint[1], initialPoint[2]]
        geo.rotateAndMove(coords, initialPoint, -90, self.raio)
        vAux = geo.createVertex(coords[0],coords[1],coords[2])
        coords = [initialPoint[0], initialPoint[1], initialPoint[2]]
        geo.rotateAndMove(coords, initialPoint, 90, self.raio)
        vAux2 = geo.createVertex(coords[0],coords[1],coords[2])
        lp = geo.createLine(vAux,vAux2)
        cp1 = geo.createCircleNormal2(self.points[-2],vAux2,vAux,self.raio,self.index)
        cp2 = geo.createCircleNormal2(self.points[-2],vAux,vAux2,self.raio,self.index)
#   
#        angle = (-30) * (math.pi/180)
#        self.vf1 = geo.createVertex(0,0,0)
#        coordsf1 = []
#        coordsf1.append(0)
#        coordsf1.append(0)
#        coordsf1.append(0)
#        origin = []
#        origin.append(v1[0])
#        origin.append(v1[1])
#        origin.append(v1[2])
#        geo.spheCoords(origin)
#        coordsf1[0] += 2*(self.raio)
#        angleR = geo.degreeToRadians(30)
#       
        
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
            