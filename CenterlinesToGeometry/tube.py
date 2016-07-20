import cubit
import geo

class Tube(object):

    def __init__(self, initialArcs, finalArcs, v1, v2):
        self.initialArcs = []
        self.finalArcs = []
        self.lines = []
        self.surfs = []
        self.initialVertexs = []
        self.finalVertexs = []
        
        self.alignPoints(initialArcs,finalArcs, v1,v2)
        for i in range(0, len(self.initialVertexs)):
            self.lines.append(geo.createLine(self.initialVertexs[i], self.finalVertexs[i]))
        self.alignCurves(initialArcs,finalArcs)
        self.genSurfaces()
    
    def genSurfaces(self):
        for i in range(0, len(self.initialArcs)-1):
            cvs = [self.initialArcs[i],self.finalArcs[i], self.lines[i],self.lines[i+1]]
            self.surfs.append(geo.createSurfaceCurve2(cvs))
        cvs = [self.initialArcs[-1],self.finalArcs[-1], self.lines[-1],self.lines[0]]
        self.surfs.append(geo.createSurfaceCurve2(cvs))
        self.surfs.append(geo.createSurfaceCurve2(self.initialArcs))
        self.surfs.append(geo.createSurfaceCurve2(self.finalArcs))
        
    def alignPoints(self,initialArcs,finalArcs, v1, v2):
        for i in range(0, len(initialArcs)):
            iVt = initialArcs[i][0]
            fIndex = 0
            minAngle = 360
            for f in range(0, len(finalArcs)):
                vt = finalArcs[f][0]
                angle = self.getAngle(v1, v2, iVt, vt)
                if angle < minAngle:
                    minAngle = angle
                    fIndex = f
            self.initialVertexs.append(initialArcs[i][0])
            self.finalVertexs.append(finalArcs[fIndex][0])
        
    def alignCurves(self,initialArcs,finalArcs):
        for i in range(0, len(self.initialVertexs)-1):
            for j in range(0, len(initialArcs)):
                arc = initialArcs[j][1]
                vts = cubit.get_relatives("curve", arc, "vertex")
                d00 = cubit.get_distance_between(vts[0], self.initialVertexs[i])
                d11 = cubit.get_distance_between(vts[1], self.initialVertexs[i+1])
                d1 = d00 + d11
                d01 = cubit.get_distance_between(vts[0], self.initialVertexs[i+1])
                d10 = cubit.get_distance_between(vts[1], self.initialVertexs[i])
                d2 = d01+d10
                if d1 < 0.00001 or d2 < 0.00001:
                    self.initialArcs.append(arc)
                    break
        for j in range(0, len(initialArcs)):
            arc = initialArcs[j][1]
            vts = cubit.get_relatives("curve", arc, "vertex")
            d00 = cubit.get_distance_between(vts[0], self.initialVertexs[-1])
            d11 = cubit.get_distance_between(vts[1], self.initialVertexs[0])
            d1 = d00 + d11
            d01 = cubit.get_distance_between(vts[0], self.initialVertexs[-1])
            d10 = cubit.get_distance_between(vts[1], self.initialVertexs[0])
            d2 = d01+d10
            if d1 < 0.00001 or d2 < 0.00001:
                self.initialArcs.append(arc)
                break
        
        for i in range(0, len(self.finalVertexs)-1):
            for j in range(0, len(finalArcs)):
                arc = finalArcs[j][1]
                vts = cubit.get_relatives("curve", arc, "vertex")
                d00 = cubit.get_distance_between(vts[0], self.finalVertexs[i])
                d11 = cubit.get_distance_between(vts[1], self.finalVertexs[i+1])
                d1 = d00 + d11
                d01 = cubit.get_distance_between(vts[0], self.finalVertexs[i+1])
                d10 = cubit.get_distance_between(vts[1], self.finalVertexs[i])
                d2 = d01+d10
                if d1 < 0.00001 or d2 < 0.00001:
                    self.finalArcs.append(arc)
                    break
        for j in range(0, len(finalArcs)):
            arc = finalArcs[j][1]
            vts = cubit.get_relatives("curve", arc, "vertex")
            d00 = cubit.get_distance_between(vts[0], self.finalVertexs[-1])
            d11 = cubit.get_distance_between(vts[1], self.finalVertexs[0])
            d1 = d00 + d11
            d01 = cubit.get_distance_between(vts[0], self.finalVertexs[-1])
            d10 = cubit.get_distance_between(vts[1], self.finalVertexs[0])
            d2 = d01+d10
            if d1 < 0.00001 or d2 < 0.00001:
                self.finalArcs.append(arc)
                break
            
    def genVolume(self):
        bodies = []
        for s in self.surfs:
            bodies.append(cubit.get_owning_body("surface", s))
        geo.imprintBody(bodies) 
        geo.mergeBody(bodies) 
        self.vol = geo.createVolume(self.surfs)
            
    def getAngle(self, v1, v2, v3, v4):
        initialPoint1 = cubit.vertex(v1).coordinates()
        finalPoint1 = cubit.vertex(v2).coordinates()
        anglesRoot = geo.getAngles(initialPoint1, finalPoint1)
        ay = geo.radiansToDegree(anglesRoot[0])
        initialPoint2 = cubit.vertex(v3).coordinates()
        finalPoint2 = cubit.vertex(v4).coordinates()
        angles = geo.getAngles(initialPoint2, finalPoint2)
        ay2 = geo.radiansToDegree(angles[0])
        a = ay2-ay
        if a > 180:
            a -= 360
        if a < -180:
            a += 360
            
        az = geo.radiansToDegree(anglesRoot[1])
        az2 = geo.radiansToDegree(angles[1])
        b = az2-az
        if b > 180:
            b -= 360
        if b < -180:
            b += 360
        return abs(a)+abs(b)
        