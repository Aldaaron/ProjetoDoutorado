from __future__ import print_function
import cubit
import geo
import geo2
import sys
import math

def eprint(*args, **kwargs):
        print(*args, file=sys.stderr, **kwargs)

class Tube(object):

#     def __init__(self, initialArcs, finalArcs, v1, v2):
#         self.initialArcs = []
#         self.finalArcs = []
#         self.lines = []
#         self.surfs = []
#         self.initialVertexs = []
#         self.finalVertexs = []
# 
#         self.initialVertexs = [initialArcs[0][0],initialArcs[1][0]]
#         self.finalVertexs = [0,0]
#         iList = [initialArcs[0][0],initialArcs[1][0]]
#         fList = [finalArcs[0][0],finalArcs[1][0]]
#         self.alignPoints2(0,iList,fList, v1,v2, None)
#         self.checkPoints()
#         for i in range(0, len(self.initialVertexs)):
#             self.lines.append(geo.createLine(self.initialVertexs[i], self.finalVertexs[i]))  
#         self.alignCurves2(initialArcs,finalArcs) 
        
    def __init__(self, arcsI, arcsF):
        self.initialArcs = []
        self.finalArcs = []
        self.surfs = []
        if len(arcsI) == 1:
            self.initialArcs.append(arcsI[0])
            self.finalArcs.append(arcsF[0])
        else:
            self.alignCurves2(arcsI, arcsF)
            
    def checkPoints(self):
        if len(self.finalVertexs) != len(set(self.finalVertexs)):
            eprint("Erro")
            eprint(self.finalVertexs)
            
    def genSurfaces(self):
        for i in range(0, len(self.initialArcs)):
            self.surfs.append(geo.createSkinCurve(self.initialArcs[i], self.finalArcs[i]))
        self.surfs.append(geo.createSurfaceCurve2(self.initialArcs))
        self.surfs.append(geo.createSurfaceCurve2(self.finalArcs))

    def genSurfacesI(self):
        self.surfs.append(geo.createSurfaceCurve2(self.initialArcs))
        
    def genSurfacesF(self):
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
    
    def alignPoints2(self,i, iList, fList, v1, v2, add):
        finalVertexs = []
        for j in range(i, len(iList)): 
            fIndex = self.checkVertex(fList, v1, v2, iList[j])
            if j == i:
                ele = fList[fIndex]
            if fList[fIndex] == ele:
                finalVertexs.append(fList[fIndex])
        if len(finalVertexs) == 1:
            self.finalVertexs[i] = ele
            fList.remove(ele)
            if add is not None:
                fList.append(add)
            if i < 1:
                self.alignPoints2(i+1, iList, fList, v1, v2, None)
        else:
            angles = []
            for k in range(0, len(finalVertexs)):
                angles.append(self.getAngle(v1, v2, iList[i], finalVertexs[k]))
            if angles.index(min(angles)) == 0:
                self.finalVertexs[i] = ele
                fList.remove(ele)
                if add is not None:
                    fList.append(add)
                if i < 1:
                    self.alignPoints2(i+1, iList, fList, v1, v2, None)
            else:
                fList.remove(ele)
                self.alignPoints2(i, iList, fList, v1, v2, ele)
        
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
            d01 = cubit.get_distance_between(vts[0], self.initialVertexs[0])
            d10 = cubit.get_distance_between(vts[1], self.initialVertexs[-1])
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
            d01 = cubit.get_distance_between(vts[0], self.finalVertexs[0])
            d10 = cubit.get_distance_between(vts[1], self.finalVertexs[-1])
            d2 = d01+d10
            if d1 < 0.00001 or d2 < 0.00001:
                self.finalArcs.append(arc)
                break
            
    def alignCurves2(self,initialArcs,finalArcs):
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
        