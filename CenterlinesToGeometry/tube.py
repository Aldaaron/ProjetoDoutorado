from __future__ import print_function
import cubit
import geo
import geo2
import sys
import math

def eprint(*args, **kwargs):
        print(*args, file=sys.stderr, **kwargs)

class Tube(object):

    def __init__(self, initialArcs, finalArcs, v1, v2):
        self.initialArcs = []
        self.finalArcs = []
        self.lines = []
        self.surfs = []
        self.initialVertexs = []
        self.finalVertexs = []

        self.initialVertexs = [initialArcs[0][0],initialArcs[1][0],initialArcs[2][0],initialArcs[3][0]]
        self.finalVertexs = [0,0,0,0]
        iList = [initialArcs[0][0],initialArcs[1][0],initialArcs[2][0],initialArcs[3][0]]
        fList = [finalArcs[0][0],finalArcs[1][0],finalArcs[2][0],finalArcs[3][0]]
        self.alignPoints2(0,iList,fList, v1,v2, None)
        self.checkPoints()
        for i in range(0, len(self.initialVertexs)):
            self.lines.append(geo.createLine(self.initialVertexs[i], self.finalVertexs[i]))  
        self.alignCurves(initialArcs,finalArcs)
        
#         if len(self.initialArcs) < 4:
#             eprint("Faltam arcos iniciais: " + str(self.initialArcs) + str(initialArcs))
#         if len(self.finalArcs) < 4:
#             eprint("Faltam arcos finais: " + str(self.finalArcs)+ str(finalArcs))
#         if len(self.lines) < 4:
#             eprint("Faltam linhas: " + str(self.lines))

    def checkPoints(self):
        if len(self.finalVertexs) != len(set(self.finalVertexs)):
            eprint("Erro")
            eprint(self.finalVertexs)
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
            if i < 3:
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
                if i < 3:
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
        