import geo
import cubit
import math
from settings import *
from __builtin__ import range

class Vessel(object):

    def __init__(self, index, vI, vF, raio, root):
        self.index = index
        self.index2 = index
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
        np += int(dist/(self.raio*5))
        inc = 1.0/np
        for x in range(1, np):
            if self.root is not None and x >= np-2:
                vertexID = geo.createVertexOnCurveFraction(self.index, inc*x)
            else:
                vertexID = geo.createVertexOnCurveFraction(self.index2, inc*x)
            self.points.append(vertexID)
        self.points.append(self.vF)
        
    def genSegment(self, vi, vf, normal, raio, linkCurves):
        initialPoint = cubit.vertex(vi).coordinates()
        finalPoint = cubit.vertex(vf).coordinates() 
        angles = geo.getAngles(initialPoint, finalPoint)
        ay = geo.radiansToDegree(angles[0]+geo.degreeToRadians(-90))
        inc = 360/nSplit
        pi = []
        arcsi = []
        if linkCurves is None:
            for i in xrange(-180,180, inc):
                coords = list(initialPoint)
                geo.rotateAndMove(coords, initialPoint, ay, i, raio)
                p = geo.createVertex(coords[0],coords[1],coords[2])
                pi.append(p)
            for i in range(0,len(pi)-1):
                arc = geo.createCircleNormal2(vi,pi[i],pi[i+1],raio,normal)
                arcsi.append(arc)
            arc = geo.createCircleNormal2(vi,pi[-1],pi[0],raio,normal)
            arcsi.append(arc)
        else:
            for i in range(0,len(linkCurves)/2):
                pi.append(linkCurves[(len(linkCurves)/2)+i])
                arcsi.append(linkCurves[i])
        pf = []
        arcsf = []
        for i in xrange(-180,180, inc):
            coords = list(finalPoint)
            geo.rotateAndMove(coords, finalPoint, ay, i, raio)
            p = geo.createVertex(coords[0],coords[1],coords[2])
            pf.append(p)
        for i in range(0,len(pf)-1):
            arc = geo.createCircleNormal2(vf,pf[i],pf[i+1],raio,normal)
            arcsf.append(arc)
        arc = geo.createCircleNormal2(vf,pf[-1],pf[0],raio,normal)
        arcsf.append(arc)
        horizontals = []
        for i in range(0,len(pf)):
            horizontals.append(geo.createLine(pi[i],pf[i]))
        curves = []
        for i in range(0,nSplit-1):
            curves = [arcsi[i],horizontals[i],horizontals[i+1],arcsf[i]]
            surf = geo.createSurfaceCurve2(curves)
            self.surfs.append(surf)
        curves = [arcsi[-1],horizontals[-1],horizontals[0],arcsf[-1]]
        surf = geo.createSurfaceCurve2(curves)
        self.surfs.append(surf)
        if vi == self.vI:
            self.surfs.append(geo.createSurfaceCurve2(arcsi))
        nextLinkElements = []
        nextLinkElements.extend(arcsf)
        nextLinkElements.extend(pf)
        return nextLinkElements
    
    def linkToRoot(self):
        if self.root.nSons == 2:
            otherSon = self.root.son1
            if otherSon == self:
                otherSon = self.root.son2
            d1 = cubit.get_distance_between(self.vF, self.root.vf1)+cubit.get_distance_between(otherSon.vF, self.root.vf2)
            d2 = cubit.get_distance_between(self.vF, self.root.vf2)+cubit.get_distance_between(otherSon.vF, self.root.vf1)
            if d1<d2:
                self.vI = self.root.vf1
                lineID = geo.createLine(self.vI,self.vF)
                self.index = lineID
                linkCircles = self.root.linkCirclesf1
            else:
                self.vI = self.root.vf2 
                lineID = geo.createLine(self.vI,self.vF)
                self.index = lineID
                linkCircles = self.root.linkCirclesf2
        else:
            self.vI = self.root.vf1
            lineID = geo.createLine(self.vI,self.vF)
            self.index = lineID
            linkCircles = self.root.linkCirclesf1 
        return linkCircles

    def genSurfaces(self):
        nextLinkElements = None
        if self.root is not None:
            nextLinkElements = self.linkToRoot()
        self.genIntervals()
        if self.nSons == 2:
            for i in range(0, len(self.points)-2):
                raio = self.raio
                if self.root is not None:
                    if i == 0:
                        raio = 0.75*self.root.raio
                    else:
                        raio = (i*self.raio + (len(self.points)-2-i) * self.root.raio*0.75) / (len(self.points)-2)
                nextLinkElements = self.genSegment(self.points[i], self.points[i+1], self.index, raio, nextLinkElements)
            self.genSon(nextLinkElements, 30)
            self.genSon(nextLinkElements, -30)
        else:
            for i in range(0, len(self.points)-1):
                raio = self.raio
                if self.root is not None:
                    if i == 0:
                        raio = 0.75*self.root.raio
                    else:
                        raio = (i*self.raio + (len(self.points)-1-i) * self.root.raio*0.75) / (len(self.points)-1)
                nextLinkElements = self.genSegment(self.points[i], self.points[i+1], self.index, raio, nextLinkElements)
            self.vf1 = self.vF
            self.linkCirclesf1 = nextLinkElements
            close = []
            for i in range(0, (len(nextLinkElements)/2)):
                close.append(nextLinkElements[i])
            self.surfs.append(geo.createSurfaceCurve2(close))
            
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
                lineIDAux = geo.createSpline(self.root.points[-2], self.root.vf1,self.vF)
                print [self.root.points[-2], self.root.vf1,self.vF]
                lineID2 = geo.splitCurve(lineIDAux, self.root.vf1)
                self.index = lineID
                self.index2 = lineID2
                self.circles.append(self.root.cf1)
            else:
                self.vI = self.root.vf2 
                lineID = geo.createLine(self.root.vf2,self.vF)
                lineIDAux = geo.createSpline(self.root.points[-2],self.root.vf2,self.vF)
                print [self.root.points[-2], self.root.vf1,self.vF]
                lineID2 = geo.splitCurve(lineIDAux, self.root.vf2)
                self.index = lineID
                self.index2 = lineID2
                self.circles.append(self.root.cf2)
        else:
            self.vI = self.root.vf1
            lineID = geo.createLine(self.root.vf1,self.vF)
            lineIDAux = geo.createSpline(self.root.points[-2], self.root.vf1,self.vF)
            lineID2 = geo.splitCurve(lineIDAux, self.root.vf1)
            self.index = lineID
            self.index2 = lineID2
            self.circles.append(self.root.cf1)     
            
    def genSurfaces2(self):
        if self.root is not None:
            self.linkToRoot2()
        self.genIntervals()
        if self.nSons == 2:
            if self.root is None:
                for x in range(0, len(self.points)-2):
                    circleID = geo.createCircleNormal(self.points[x],self.raio,self.index2)
                    self.circles.append(circleID)
            else:
                for x in range(1, len(self.points)-2):
                    raio = (x*self.raio + (len(self.points)-2-x) * self.root.raio*0.75) / (len(self.points)-2)
                    circleID = geo.createCircleNormal(self.points[x],raio,self.index2)
                    self.circles.append(circleID)
            linksF1 = self.genSon2(-30)
            linksF2 = self.genSon2(30)
            linkCircles = []
            linkCircles.extend(linksF1)
            linkCircles.extend(linksF2)
#             circleID = geo.createCombineCurve(linkCircles)
#             self.circles.append(circleID)
            geo.imprintVolumes(self.volf1, self.volf2)
            geo.mergeVolumes(self.volf1, self.volf2)
            volFilhos = geo.unionVolumes(self.volf1, self.volf2)
            initialPoint = cubit.vertex(self.points[-2]).coordinates()
            finalPoint = cubit.vertex(self.points[-1]).coordinates() 
            angles = geo.getAngles(initialPoint, finalPoint)
            ay = geo.radiansToDegree(angles[0]+geo.degreeToRadians(-90))
            az = geo.radiansToDegree(angles[1])
            coords = list(initialPoint)
            geo.rotateAndMove(coords, initialPoint, ay, -90, self.raio)
            vEli2 = geo.createVertex(coords[0],coords[1],coords[2])
            coords = list(initialPoint)
            geo.rotateAndMove(coords, initialPoint, ay, 0, self.raio)
            vEli1 = geo.createVertex(coords[0],coords[1],coords[2])
            circleID = geo.createEllipseFull(vEli1, vEli2, self.points[-2])
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
                    if x < len(self.points)-2:
                        circleID = geo.createCircleNormal(self.points[x],raio,self.index2)
                    else:
                        circleID = geo.createCircleNormal(self.points[x],raio,self.index2)
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
    
    def genSon(self, nextLinkElements, angle):
        coords = list(cubit.vertex(self.points[-1]).coordinates())
        geo.rotate2(coords, cubit.vertex(self.points[-2]).coordinates() , angle, 2*self.raio)
        sonVertex = geo.createVertex(coords[0],coords[1],coords[2])   
        alongLine = geo.createLine(self.points[-2],sonVertex)
        initialPoint = cubit.vertex(self.points[-2]).coordinates()
        finalPoint = cubit.vertex(sonVertex).coordinates() 
        angles = geo.getAngles(initialPoint, finalPoint)
        ay = geo.radiansToDegree(angles[0]+geo.degreeToRadians(-90)) 
        piAux = [] 
        arcsiAux = [] 
        for i in range(0, len(nextLinkElements)/2):
            piAux.append(nextLinkElements[len(nextLinkElements)/2+i]) 
            arcsiAux.append(nextLinkElements[i])
        inc = 360/nSplit
        pi = []
        arcsi = []
        if angle > 0:
            self.splitPoints = []
            self.splitLines = []
            sl = geo.createLine(piAux[(len(piAux)/4)],piAux[((len(piAux)*3)/4)])
            for i in xrange(0,180, inc):
                pct = float(i)/180
                self.splitPoints.append(geo.createVertexOnCurveFraction(sl, pct))
            self.splitPoints.append(piAux[((len(piAux)*3)/4)])
            self.splitPoints[0] = piAux[(len(piAux)/4)]
            self.splitPoints[len(self.splitPoints)/2] = self.points[-2]
            for i in range(0,len(self.splitPoints)-1):
                self.splitLines.append(geo.createLine(self.splitPoints[i],self.splitPoints[i+1])) 
            for i in range(0,(len(piAux)/4)):  
                pi.append(piAux[i])
                arcsi.append(arcsiAux[i])
            for i in range(0,len(self.splitPoints)-1):
                pi.append(self.splitPoints[i])
                arcsi.append(self.splitLines[i])
            for i in range((len(piAux)*3)/4,len(piAux)):  
                pi.append(piAux[i])
                arcsi.append(arcsiAux[i])
        else:
            for i in range(0,len(self.splitPoints)/2):
                pi.append(self.splitPoints[(len(self.splitPoints)/2)-i])
                arcsi.append(self.splitLines[(len(self.splitLines)/2)-1-i])
            for i in range((len(piAux)/4),(len(piAux)*3)/4):  
                pi.append(piAux[i])
                arcsi.append(arcsiAux[i])
            for i in range(0,len(self.splitPoints)/2):
                pi.append(self.splitPoints[-i-1])
                arcsi.append(self.splitLines[-i-1])
        sonPoint = cubit.vertex(sonVertex).coordinates()
        pf = []
        arcsf = []
        for i in xrange(-180,180, inc):
            coords = list(sonPoint)
            geo.rotateAndMove(coords, sonPoint, ay, i, 0.75*self.raio)
            pf.append(geo.createVertex(coords[0],coords[1],coords[2]))
        for i in range(0,len(pf)-1):
            arc = geo.createCircleNormal2(sonVertex,pf[i],pf[i+1],0.75*self.raio,alongLine)
            arcsf.append(arc)
        arc = geo.createCircleNormal2(sonVertex,pf[-1],pf[0],0.75*self.raio,alongLine)
        arcsf.append(arc)
        horizontals = []
        for i in range(0,nSplit):
            horizontals.append(geo.createLine(pi[i],pf[i]))
        curves = []
        for i in range(0,nSplit-1):
            curves = [arcsi[i],horizontals[i],horizontals[i+1],arcsf[i]]
            self.surfs.append(geo.createSurfaceCurve2(curves))
        curves = [arcsi[-1],horizontals[-1],horizontals[0],arcsf[-1]]
        self.surfs.append(geo.createSurfaceCurve2(curves))
        if angle > 0:
            self.vf1 = sonVertex
            self.linkCirclesf1 = list(arcsf)
            self.linkCirclesf1.extend(pf)
        else:
            self.vf2 = sonVertex  
            self.linkCirclesf2 = list(arcsf)
            self.linkCirclesf2.extend(pf)  
        self.surfs.append(geo.createSurfaceCurve2(arcsf))
    
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
            