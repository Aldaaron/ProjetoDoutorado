from __future__ import print_function
import cubit
import geo
import geo2
import sys
import math

def eprint(*args, **kwargs):
        print(*args, file=sys.stderr, **kwargs)

class Tube(object):

        
    def __init__(self, arcsI, arcsF):
        self.initialArcs = []
        self.finalArcs = []
        self.surfs = []
        if len(arcsI) == 1:
            self.initialArcs.append(arcsI[0])
            self.finalArcs.append(arcsF[0])
        else:
#             arcI = geo.createCombineCurve(arcsI)
#             arcF = geo.createCombineCurve(arcsF)
#             geo.deleteCurve(arcsI[0])
#             geo.deleteCurve(arcsI[1])
#             geo.deleteCurve(arcsF[0])
#             geo.deleteCurve(arcsF[1])
#             self.initialArcs.append(arcI)
#             self.finalArcs.append(arcF)
            self.alignCurves(arcsI, arcsF)
            
    def checkPoints(self):
        if len(self.finalVertexs) != len(set(self.finalVertexs)):
            eprint("Erro")
            eprint(self.finalVertexs)
            
    def genSurfaces(self):
#         if len(self.initialArcs) == 2:
#             self.surfs.append(geo.createSurfaceCurve2([self.initialArcs[0], self.line]))
#             self.surfs.append(geo.createSurfaceCurve2([self.initialArcs[1], self.line]))
#         else:
#             self.surfs.append(geo.createSurfaceCurve2(self.initialArcs))
        for i in range(0, len(self.initialArcs)):
            self.surfs.append(geo.createSkinCurve(self.initialArcs[i], self.finalArcs[i]))
#         self.surfs.append(geo.createSurfaceCurve2(self.finalArcs))

    def genSurfacesI(self):
        self.surfs.append(geo.createSurfaceCurve2(self.initialArcs))
        
    def genSurfacesF(self):
        self.surfs.append(geo.createSurfaceCurve2(self.finalArcs))
        
    def checkVertex(self, list, v1,v2,v):
        fIndex = 0
        minAngle = 360
        for f in range(0, len(list)):
            vt = list[f]
            angle = self.getAngle(v1, v2, v, vt)
            if angle < minAngle:
                minAngle = angle
                fIndex = f
        return fIndex
            
    def alignCurves(self,initialArcs,finalArcs):
        self.initialArcs.append(initialArcs[0])
        self.initialArcs.append(initialArcs[1])
        mid1i = geo.createVertexOnCurveFraction(initialArcs[0],0.5)
        mid2i = geo.createVertexOnCurveFraction(initialArcs[1],0.5)
        mid1f = geo.createVertexOnCurveFraction(finalArcs[0],0.5)
        mid2f = geo.createVertexOnCurveFraction(finalArcs[1],0.5)
        d1 = cubit.get_distance_between(mid1i, mid1f)
        d2 = cubit.get_distance_between(mid1i, mid2f)
        if d1 < d2:
            self.finalArcs.append(finalArcs[0])
            self.finalArcs.append(finalArcs[1])
        else:
            self.finalArcs.append(finalArcs[1])
            self.finalArcs.append(finalArcs[0])
        geo.deleteVertex(mid1i)
        geo.deleteVertex(mid2i)
        geo.deleteVertex(mid1f)
        geo.deleteVertex(mid2f)
            
            
    def genVolume(self):
#         bodies = []
#         for s in self.surfs:
#             bodies.append(cubit.get_owning_body("surface", s))
#         geo.imprintBody(bodies) 
#         geo.mergeBody(bodies) 
        self.vol = geo.createVolume(self.surfs)
    
    def getAngle(self, v1, v2, v3, v4):
        initialPoint1 = cubit.vertex(v1).coordinates()
        finalPoint1 = cubit.vertex(v2).coordinates()
        initialPoint2 = cubit.vertex(v3).coordinates()
        finalPoint2 = cubit.vertex(v4).coordinates()
        vec1 = geo2.getVector(initialPoint1,finalPoint1)
        vec2 = geo2.getVector(initialPoint2,finalPoint2)
        dot = geo2.dotProduct(vec1, vec2)
        n1 = geo2.normVector(vec1)
        n2 = geo2.normVector(vec2)
        cosA = dot/(n1*n2)
        if cosA > 1.0:
            cosA = 1.0
        if cosA < -1.0:
            cosA = -1.0
        angle = math.acos(cosA)
        return abs(geo.radiansToDegree(angle))
        