import geo
import geo2
import cubit
import math
import random
from settings import *
from branch import Branch
from geo import distance, getAngles
from tube import Tube
import tube

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
        self.surfVols = []
    
    def preProcess(self):
        self.fixSizes()
        self.split()
        self.smoth()
        self.clean()
        cubit.reset()
        
    def clean(self):
        for x in range(0, len(self.points)):
            self.points[x][0] = round(self.points[x][0], 5)
            self.points[x][1] = round(self.points[x][1], 5)
            self.points[x][2] = round(self.points[x][2], 5)
            
    def smothRadius(self):
        heap = [self.root]
        while not (len(heap) == 0):
            atual = heap.pop()
            if atual.son1 is not None and atual.son2 is not None:
                atual.son1.son1.nSplit = 1
                atual.son2.son1.nSplit = 1
                atual.son1.son1.rootRadius = atual.radius
                atual.son2.son1.rootRadius = atual.radius
                fpr = self.points[atual.finalPoint]
                fp = self.points[atual.son1.finalPoint]
                fp2 = self.points[atual.son2.finalPoint]
                u_0 = [fp[0]-fpr[0],fp[1]-fpr[1],fp[2]-fpr[2]]
                u_1 = [fp2[0]-fpr[0],fp2[1]-fpr[1],fp2[2]-fpr[2]]
                u = geo2.crossProduct(u_0, u_1)
                u = geo2.getVector4_3(u)
                geo2.normalizeVector(u)
                atual.son1.son1.u = u
                atual.son2.son1.u = u
                
                if atual.son1.son1 is not None:
                    heap.append(atual.son1.son1)
                if atual.son2.son1 is not None:
                    heap.append(atual.son2.son1)
            else:
                if atual.son1 is not None:
                    atual.son1.nSplit =  atual.nSplit + 1
                    atual.son1.rootRadius = atual.rootRadius
                    atual.son1.u = atual.u
                    heap.append(atual.son1)
                
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
            pi = self.points[atual.initialPoint]
            pf = self.points[atual.finalPoint]
            lenght = geo.distance(pi, pf)
            if lenght < atual.radius*3:
                print "Segmento muito pequeno "
                print atual.initialPoint+1
            #print str(atual.initialPoint) + str(self.points[atual.initialPoint])+ " " + str(atual.finalPoint) + str(self.points[atual.finalPoint])
            #print str(atual.initialPoint) + " " + str(atual.finalPoint)
            if atual.son1 is not None:
                heap.append(atual.son1)
            if atual.son2 is not None:
                heap.append(atual.son2)
                
    def fixSizes(self):
        heap = [self.root]
        while not (len(heap) == 0):
            atual = heap.pop()
            ip = self.points[atual.initialPoint]
            fp = self.points[atual.finalPoint]
            lenght = geo.distance(ip, fp)
            if lenght < atual.radius*5:
                t=(atual.radius*6)/lenght
                npf = [(1-t)*ip[0]+t*fp[0],((1-t)*ip[1]+t*fp[1]),((1-t)*ip[2]+t*fp[2])]
                #self.points[atual.finalPoint] = npf
                xyz = [npf[0]-fp[0], npf[1]-fp[1], npf[2]-fp[2]]
                if atual.son1 is not None:
                    self.move(atual.son1,xyz)
                if atual.son2 is not None:
                    self.move(atual.son2,xyz)
            if atual.son1 is not None:
                heap.append(atual.son1)
            if atual.son2 is not None:
                heap.append(atual.son2)
                
    def move(self, branch, xyz):
        heap = [branch]
        while not (len(heap) == 0):
            atual = heap.pop()
            if branch.root.son1 == atual:
                self.points[atual.initialPoint][0] += xyz[0]
                self.points[atual.initialPoint][1] += xyz[1]
                self.points[atual.initialPoint][2] += xyz[2]
            self.points[atual.finalPoint][0] += xyz[0]
            self.points[atual.finalPoint][1] += xyz[1]
            self.points[atual.finalPoint][2] += xyz[2]
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
                self.bif(atual) 
            if atual.son1 is not None:
                heap.append(atual.son1)
            if atual.son2 is not None:
                heap.append(atual.son2) 
        self.sort() 
    
    def smoth(self):
        heap = [self.root]
        while not (len(heap) == 0):
            atual = heap.pop()
            if atual.son1 is not None and atual.son2 is not None:
                rSon1 = atual.son1.son1  
                r1i = self.points[atual.son1.initialPoint] 
                r1f = self.points[atual.son1.finalPoint]
                f1i = self.points[rSon1.initialPoint]
                f1f = self.points[rSon1.finalPoint]
                v1 = geo.createVertex(r1i[0], r1i[1], r1i[2])
                v2 = geo.createVertex(f1i[0], f1i[1], f1i[2])
                v3 = geo.createVertex(f1f[0], f1f[1], f1f[2])
                spAux = geo.createSpline(v1, v2, v3)
                sp = geo.splitCurve(spAux, v2)
                newV = geo.createVertexOnCurveFraction(sp, 0.9)
                newPoint = cubit.vertex(newV).coordinates()
                newPoint = [newPoint[0],newPoint[1],newPoint[2]]
                self.points.append(newPoint)
                self.splitBranch(rSon1, newPoint)
                newV = geo.createVertexOnCurveFraction(sp, 0.75)
                newPoint = cubit.vertex(newV).coordinates()
                newPoint = [newPoint[0],newPoint[1],newPoint[2]]
                self.points.append(newPoint)
                self.splitBranch(rSon1, newPoint)
                newV = geo.createVertexOnCurveFraction(sp, 0.6)
                newPoint = cubit.vertex(newV).coordinates()
                newPoint = [newPoint[0],newPoint[1],newPoint[2]]
                self.points.append(newPoint)
                self.splitBranch(rSon1, newPoint)
                newV = geo.createVertexOnCurveFraction(sp, 0.45)
                newPoint = cubit.vertex(newV).coordinates()
                newPoint = [newPoint[0],newPoint[1],newPoint[2]]
                self.points.append(newPoint)
                self.splitBranch(rSon1, newPoint)
                newV = geo.createVertexOnCurveFraction(sp, 0.3)
                newPoint = cubit.vertex(newV).coordinates()
                newPoint = [newPoint[0],newPoint[1],newPoint[2]]
                self.points.append(newPoint)
                self.splitBranch(rSon1, newPoint)
                newV = geo.createVertexOnCurveFraction(sp, 0.15)
                newPoint = cubit.vertex(newV).coordinates()
                newPoint = [newPoint[0],newPoint[1],newPoint[2]]
                self.points.append(newPoint)
                self.splitBranch(rSon1, newPoint)
                
                rSon2 = atual.son2.son1
                r2i = self.points[atual.son2.initialPoint] 
                r2f = self.points[atual.son2.finalPoint]
                f2i = self.points[rSon2.initialPoint]
                f2f = self.points[rSon2.finalPoint]
#                 t = 0.9
#                 ref = [(1-t)*f2i[0]+t*f2f[0],((1-t)*f2i[1]+t*f2f[1]),((1-t)*f2i[2]+t*f2f[2])]
#                 self.splitBranch(rSon2, ref)
#                 f2f = self.points[rSon2.finalPoint]
                v1 = geo.createVertex(r2i[0], r2i[1], r2i[2])
                v2 = geo.createVertex(f2i[0], f2i[1], f2i[2])
                v3 = geo.createVertex(f2f[0], f2f[1], f2f[2])
                spAux = geo.createSpline(v1, v2, v3)
                sp = geo.splitCurve(spAux, v2)
                newV = geo.createVertexOnCurveFraction(sp, 0.9)
                newPoint = cubit.vertex(newV).coordinates()
                newPoint = [newPoint[0],newPoint[1],newPoint[2]]
                self.points.append(newPoint)
                self.splitBranch(rSon2, newPoint)
                newV = geo.createVertexOnCurveFraction(sp, 0.75)
                newPoint = cubit.vertex(newV).coordinates()
                newPoint = [newPoint[0],newPoint[1],newPoint[2]]
                self.points.append(newPoint)
                self.splitBranch(rSon2, newPoint)
                newV = geo.createVertexOnCurveFraction(sp, 0.6)
                newPoint = cubit.vertex(newV).coordinates()
                newPoint = [newPoint[0],newPoint[1],newPoint[2]]
                self.points.append(newPoint)
                self.splitBranch(rSon2, newPoint)
                newV = geo.createVertexOnCurveFraction(sp, 0.45)
                newPoint = cubit.vertex(newV).coordinates()
                newPoint = [newPoint[0],newPoint[1],newPoint[2]]
                self.points.append(newPoint)
                self.splitBranch(rSon2, newPoint)
                newV = geo.createVertexOnCurveFraction(sp, 0.3)
                newPoint = cubit.vertex(newV).coordinates()
                newPoint = [newPoint[0],newPoint[1],newPoint[2]]
                self.points.append(newPoint)
                self.splitBranch(rSon2, newPoint)
                newV = geo.createVertexOnCurveFraction(sp, 0.15)
                newPoint = cubit.vertex(newV).coordinates()
                newPoint = [newPoint[0],newPoint[1],newPoint[2]]
                self.points.append(newPoint)
                self.splitBranch(rSon2, newPoint)
                
                heap.append(rSon1)
                heap.append(rSon2)
            else:
                if atual.son1 is not None:
                    heap.append(atual.son1)
                if atual.son2 is not None:
                    heap.append(atual.son2) 
        self.sort() 
                    
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
        
    def save2(self):
        pointsFile = open(folder+"pontos3.txt", 'w')
        linesFile = open(folder+"linhas3.txt", 'w')
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
        count = 0
        while not (len(heap) == 0):
            count += 1
            if count < 1500:
                atual = heap.pop()
            else:
                break
            if atual.root is None:
                circleI = geo.createCircleNormal(atual.initialPoint+1, atual.radius, atual.index)
            else:
                circleI = atual.root.link
            if atual.nSplit > 0 and atual.nSplit < 6:
                ip = self.points[atual.initialPoint]
                fp = self.points[atual.finalPoint]
                lenght = geo.distance(ip, fp)
                t=(lenght-atual.getRadius()[1])/lenght
                ref = [(1-t)*ip[0]+t*fp[0],((1-t)*ip[1]+t*fp[1]),((1-t)*ip[2]+t*fp[2])]
                circleF = self.genArcFull(atual, atual.finalPoint+1, ref)
                if circleF == -1:
                    circleF = geo.createCircleNormal(atual.finalPoint+1, atual.radius, atual.index)
            else:
                circleF = geo.createCircleNormal(atual.finalPoint+1, atual.radius, atual.index)
            atual.circleI = circleI
            atual.circleF = circleF
            if atual.son1 is not None and atual.son2 is not None:
                self.genSons(atual)
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
                    atual.link = circleF#geo.copyCurve(circleF)
                    heap.append(atual.son1)
        self.genSurfaces()
        self.genVolumes()
        self.vol = cubit.get_entities("volume")[0]
        bb = geo.createBoudingBox(self.vol, self.root.radius/4)
        self.vol = geo.subtractVolumes(bb, self.vol)
        cubit.cmd('Color Define "%s" RGB %f %f %f' % ("darkred", 0.75,0,0))
        geo.colorVolume(self.vol,'user "darkred"')
#         fix = True
#         while fix:
#             fix = self.fixCurves()
            
    def boundary(self):
        boundary = []
        surfs = cubit.get_relatives("volume", self.vol, "surface")
        for x in range(6, len(surfs)):
            tris = cubit.get_surface_tris(surfs[x])
            for y in range(0, len(tris)):
                nodes = cubit.get_connectivity("tri", tris[y])
                p1 = cubit.get_nodal_coordinates(nodes[0])
                p2 = cubit.get_nodal_coordinates(nodes[1])
                p3 = cubit.get_nodal_coordinates(nodes[2])
                u_0 = [p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2]]
                u_1 = [p3[0]-p1[0],p3[1]-p1[1],p3[2]-p1[2]]
                u1 = geo2.crossProduct(u_0, u_1)
                u1 = geo2.getVector4_3(u1)
                geo2.normalizeVector(u1)
#                 if y == 0:
#                     center_point = cubit.get_center_point("tri", tris[y])
#                     cubit.cmd("create curve location %f %f %f direction %f %f %f length %f" %(center_point[0],center_point[1],center_point[2],u1[0],u1[1],u1[2],self.root.radius))
                element = [nodes[0]-1, nodes[1]-1, nodes[2]-1, u1[0], u1[1], u1[2]]
                boundary.append(element)
        return boundary
        
    def mesh(self):
#         surfs = cubit.get_relatives("volume", self.vol, "surface")
#         for x in range(0, len(surfs)):
#             id = surfs[x]
#             cubit.cmd("surface %d sizing function type skeleton min_size auto max_size auto max_gradient 1.5 min_num_layers_2d 1 min_num_layers_1d 1" % (id))
#             cubit.cmd("surface %d scheme TriMesh geometry approximation angle 15" % (id)) 
#             cubit.cmd("Trimesher surface gradation 1.3") 
#             cubit.cmd("mesh surface %d" % (id)) 
#         surfs = cubit.get_relatives("volume", self.vol, "surface")
#         cubit.cmd("surface %d %d %d %d %d %d size auto factor 5" % (surfs[0],surfs[1],surfs[2],surfs[3],surfs[4],surfs[5]))
#         cubit.cmd("surface %d %d %d %d %d %d scheme TriMesh geometry approximation angle 15" % (surfs[0],surfs[1],surfs[2],surfs[3],surfs[4],surfs[5]))
#         cubit.cmd("Trimesher surface gradation 1.3") 
#         cubit.cmd("mesh surface %d %d %d %d %d %d" % (surfs[0],surfs[1],surfs[2],surfs[3],surfs[4],surfs[5]))
        cubit.cmd("volume %d size auto factor 4" % (self.vol)) 
        cubit.cmd("volume %d scheme Tetmesh proximity layers off geometry approximation angle 15" % (self.vol)) 
        cubit.cmd("volume %d tetmesh growth_factor 1.3" % (self.vol)) 
        cubit.cmd("Trimesher surface gradation 1.3") 
        cubit.cmd("Trimesher volume gradation 1.3") 
        cubit.cmd("Set Tetmesher Optimize Level 6") 
        cubit.cmd("mesh volume %d" % (self.vol)) 
        
    def fixCurves(self):
        small = self.root.radius/100
        cvs = cubit.get_relatives("volume", self.vol, "curve")
        for x in range(0, len(cvs)):
            type1 = cubit.get_curve_type(cvs[x])
            virtual1 = cubit.is_virtual("curve", cvs[x])
            if type1 == "Ellipse curve" or type1 == "Arc curve" or virtual1 or type1 == "Spline curve":
                vts1 = cubit.get_relatives("curve", cvs[x], "vertex")
                for y in range(0, len(cvs)):
                    type2 = cubit.get_curve_type(cvs[y])
                    virtual2 = cubit.is_virtual("curve", cvs[y])
                    vts2 = cubit.get_relatives("curve", cvs[y], "vertex")
                    equal = False
                    arcEli = type2 == "Ellipse curve" or type2 == "Arc curve" or type2 == "Spline curve"
                    equal = type1 == type2 or virtual2 or (virtual1 and arcEli)
                    if equal and x != y and len(vts1) > 1 and len(vts2) > 1:
                        if vts1[0] == vts2[0] or vts1[0] == vts2[1] or vts1[1] == vts2[0] or vts1[1] == vts2[1]:
                            cubit.set_modified()
                            geo.compositeCurve(cvs[x], cvs[y])
                            if cubit.is_modified():
                                return True
                    
        return False                            
        
    def saveMesh(self):
        nNodes = cubit.get_node_count()
        nTets = cubit.get_tet_count()
        meshFile = open(folder+"mesh3D_noboxbound.xml", 'w')
        meshFile.write('<mesh celltype="tetrahedron" dim="3">\n')
        meshFile.write('  <nodes size="%d">\n' % (nNodes))
        for x in range(0, nNodes):
            coords = cubit.get_nodal_coordinates(x+1)
            meshFile.write('    <node id="%d" x="%f" y="%f" z="%f"/>\n' % (x,coords[0],coords[1],coords[2]))
        meshFile.write('  </nodes>\n')
        
        meshFile.write('  <elements size="%d">\n' % (nTets))
        for x in range(0, nTets):
            nodes = cubit.get_connectivity("tet", x+1)
            meshFile.write('    <element id="%d" v0="%d" v1="%d" v2="%d" v3="%d"/>\n' % (x,nodes[0]-1,nodes[1]-1,nodes[2]-1,nodes[3]-1))
        meshFile.write('  </elements>\n')
        meshFile.write('  <element_data type="fiber_transversely_isotropic">\n')
        for x in range(0, nTets):
            meshFile.write('    <element id="%d">\n' %(x))
            meshFile.write('      <fiber>1.000000,0.000000,0.000000</fiber>\n')
            meshFile.write('    </element>\n')
        meshFile.write('  </element_data>\n')
        
        meshFile.write('  <boundary celltype="triangle" dim="2">\n')
        bsurfs = cubit.get_relatives("volume", self.vol, "surface")
        ec = 0
        for x in range(6, len(bsurfs)):
            #if x is not 16 and x is not 14:
            tris = cubit.get_surface_tris(bsurfs[x])
            surf = cubit.surface(bsurfs[x])
            for y in range(0, len(tris)):
                cp = cubit.get_center_point("tri", tris[y])
                norm = surf.normal_at([cp[0],cp[1],cp[2]])
                #cubit.cmd("create curve location %f %f %f direction %f %f %f length %f" % (cp[0],cp[1],cp[2],norm[0],norm[1],norm[2],200))

                nodes = cubit.get_connectivity("tri", tris[y])
                element = [nodes[0]-1, nodes[1]-1, nodes[2]-1]
                meshFile.write('    <element id="%d" marker="%d" v0="%d" v1="%d" v2="%d" nx="%f" ny="%f" nz="%f"/>\n' % (ec,1.0,element[0],element[1],element[2],norm[0],norm[1],norm[2]))
                ec = ec+1
        meshFile.write('  </boundary>\n')

        
        
        #meshFile.write('  <boundary celltype="triangle" dim="2">\n')
        #surfs = cubit.get_relatives("volume", self.vol, "surface")
        #ec = 0
        #for x in range(0, 2):
        #    tris = cubit.get_surface_tris(surfs[x])
        #    for y in range(0, len(tris)):
        #        nodes = cubit.get_connectivity("tri", tris[y])
        #        element = [nodes[0]-1, nodes[1]-1, nodes[2]-1]
        #        ec = ec+1
        #        meshFile.write('    <element id="%d" marker="%d" v0="%d" v1="%d" v2="%d"/>\n' % (ec,x,element[0],element[1],element[2]))
        #meshFile.write('    <element id="%d" marker="%d" v0="%d" v1="%d" v2="%d"/>\n' % (ec+1,2,element[0],element[1],element[2]))
        #meshFile.write('  </boundary>\n')
        
        meshFile.write('</mesh>\n')
        #meshFile.write('<poisson>\n')
        #meshFile.write('  <neumann>\n')
        #meshFile.write('    <node id="0" marker="0" value="-1.0" />\n')
        #meshFile.write('    <node id="1" marker="1" value="1.0" />\n')
        #meshFile.write('  </neumann>\n')
        #meshFile.write('  <dirichlet>\n')
        #meshFile.write('    <node id="2" marker="2" value="0.0" />\n')
        #meshFile.write('  </dirichlet>\n')
        #meshFile.write('</poisson>\n')
        
        #meshFile.write('<electrophysiology>\n')
        #meshFile.write('  <stimuli number="1">\n')
        #bb = cubit.get_bounding_box("volume", self.vol)
        #x0 = bb[0]
        #x1 = 0.8*x0 + 0.2*bb[1]
        #y0 = bb[3]
        #y1 = 0.8*y0 + 0.2*bb[4]
        #z0 = bb[6]
        #z1 = 0.8*z0 + 0.2*bb[7]
        #meshFile.write('    <stim start="0.00" duration="4.00" value="-35.7140" x0="%f" x1="%f" y0="%f" y1="%f" z0="%f" z1="%f" />\n' % (x0,x1,y0,y1,z0,z1))
        #meshFile.write('  </stimuli>\n')
        #meshFile.write('</electrophysiology>\n')
        meshFile.close()
    
    def saveMesh3D(self):
        cubit.cmd("brick x 5000 y 5000 z 5000")
        cubit.cmd("create cylinder height 4000 radius 500")
        cubit.cmd("subtract body 2 from body 1")
        volID = cubit.get_last_id("volume")
        self.vol = volID
        self.mesh()
        
        nNodes = cubit.get_node_count()
        meshFile = open(folder+"mesh3D2.xml", 'w')
        meshFile.write('<mesh celltype="tetrahedron" dim="3">\n')
        meshFile.write('  <nodes size="%d">\n' % (nNodes))
        for x in range(0, nNodes):
            coords = cubit.get_nodal_coordinates(x+1)
            meshFile.write('    <node id="%d" x="%f" y="%f" z="%f"/>\n' % (x,coords[0],coords[1],coords[2]))
        meshFile.write('  </nodes>\n')
        
        nTets = cubit.get_tet_count()
        meshFile.write('  <elements size="%d">\n' % (nTets))
        for x in range(0, nTets):
            nodes = cubit.get_connectivity("tet", x+1)
            meshFile.write('    <element id="%d" v0="%d" v1="%d" v2="%d" v3="%d"/>\n' % (x,nodes[0]-1,nodes[1]-1,nodes[2]-1,nodes[3]-1))
        meshFile.write('  </elements>\n')
        meshFile.write('  <element_data type="fiber_transversely_isotropic">\n')
        for x in range(0, nTets):
            meshFile.write('    <element id="%d">\n' %(x))
            meshFile.write('      <fiber>1.000000,0.000000,0.000000</fiber>\n')
            meshFile.write('    </element>\n')
        meshFile.write('  </element_data>\n')
        meshFile.write('  <boundary celltype="triangle" dim="2">\n')
        bsurfs = [10, 11, 12]
        ec = 0
        for x in range(0, len(bsurfs)):
            tris = cubit.get_surface_tris(bsurfs[x])
            surf = cubit.surface(bsurfs[x])
            for y in range(0, len(tris)):
                cp = cubit.get_center_point("tri", tris[y])
                norm = surf.normal_at([cp[0],cp[1],cp[2]])
                #cubit.cmd("create curve location %f %f %f direction %f %f %f length %f" % (cp[0],cp[1],cp[2],norm[0],norm[1],norm[2],200))

                nodes = cubit.get_connectivity("tri", tris[y])
                element = [nodes[0]-1, nodes[1]-1, nodes[2]-1]
                meshFile.write('    <element id="%d" marker="%d" v0="%d" v1="%d" v2="%d" nx="%f" ny="%f" nz="%f"/>\n' % (ec,1.0,element[0],element[1],element[2],norm[0],norm[1],norm[2]))
                ec = ec+1
        meshFile.write('  </boundary>\n')
        meshFile.write('</mesh>\n')
        #meshFile.write('<poisson>\n')
        #meshFile.write('  <neumann>\n')
        #meshFile.write('    <node id="1" marker="1" value="1.0" />\n')
        #meshFile.write('  </neumann>\n')
        #meshFile.write('</poisson>\n')
        #meshFile.write('<electrophysiology>\n')
        #meshFile.write('  <stimuli number="1">\n')
        #bb = cubit.get_bounding_box("volume", self.vol)
        #x0 = bb[0]
        #x1 = 0.8*x0 + 0.2*bb[1]
        #y0 = bb[3]
        #y1 = 0.8*y0 + 0.2*bb[4]
        #z0 = bb[6]
        #z1 = 0.8*z0 + 0.2*bb[7]
        #meshFile.write('    <stim start="0.00" duration="4.00" value="-35.7140" x0="%f" x1="%f" y0="%f" y1="%f" z0="%f" z1="%f" />\n' % (x0,x1,y0,y1,z0,z1))
        #meshFile.write('  </stimuli>\n')
        #meshFile.write('</electrophysiology>\n')
        meshFile.close()
        
    def saveMesh2D2(self):
        cubit.cmd("create surface rectangle width 5000 height 5000 zplane")
        cubit.cmd("create surface rectangle width 1000 height 1000 zplane")
        cubit.cmd("subtract body 2 from body 1")
        surfID = cubit.get_last_id("surface")
        cubit.cmd("surface all size auto factor 4")
        cubit.cmd("surface all scheme TriMesh")
        cubit.cmd("mesh surface all") 
        
        nNodes = cubit.get_node_count()
        meshFile = open(folder+"mesh2D2.xml", 'w')
        meshFile.write('<mesh celltype="triangle" dim="2">\n')
        meshFile.write('  <nodes size="%d">\n' % (nNodes))
        for x in range(0, nNodes):
            coords = cubit.get_nodal_coordinates(x+1)
            meshFile.write('    <node id="%d" x="%f" y="%f" z="%f"/>\n' % (x,coords[0],coords[1],coords[2]))
        meshFile.write('  </nodes>\n')
        

        tris = cubit.get_surface_tris(surfID)
        meshFile.write('  <elements size="%d">\n' % (len(tris)))
        for x in range(0, len(tris)):
            nd = cubit.get_connectivity("tri", tris[x])
            meshFile.write('    <element id="%d" v0="%d" v1="%d" v2="%d"/>\n' % (x,nd[0]-1,nd[1]-1,nd[2]-1))
        meshFile.write('  </elements>\n')
        meshFile.write('  <element_data type="fiber_transversely_isotropic">\n')
        for x in range(0, len(tris)):
            meshFile.write('    <element id="%d">\n' %(x))
            meshFile.write('      <fiber>1.000000,0.000000,0.000000</fiber>\n')
            meshFile.write('    </element>\n')
        meshFile.write('  </element_data>\n')
        
        meshFile.write('  <boundary celltype="line" dim="1">\n')
        
        eds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

        ec = 1
        surf = cubit.surface(surfID)
        
        for y in range(0, len(eds)):
            nodes = cubit.get_connectivity("edge", eds[y])
            element = [nodes[0]-1, nodes[1]-1]
            #cp = cubit.get_center_point("edge", eds[y])
            #norm = surf.normal_at([cp[0],cp[1],cp[2]])
            #cubit.cmd("create curve location %f %f %f direction %f %f %f length %f" % (cp[0],cp[1],cp[2], norm[0], norm[1], norm[2], 200))
            #ec = ec+1
            #meshFile.write('    <element id="%d" marker="%d" v0="%d" v1="%d" nx="%f" ny="%f" nz="%f"/>\n' % (ec,1,element[0],element[1], norm[0], norm[1], norm[2]))
        #meshFile.write('    <element id="%d" marker="%d" v0="%d" v1="%d" v2="%d"/>\n' % (ec+1,2,element[0],element[1],element[2]))
        meshFile.write('  </boundary>\n')
        meshFile.write('</mesh>\n')
        meshFile.write('<poisson>\n')
        meshFile.write('  <neumann>\n')
        meshFile.write('    <node id="1" marker="1" value="1.0" />\n')
        meshFile.write('  </neumann>\n')
        #meshFile.write('  <dirichlet>\n')
        #meshFile.write('    <node id="2" marker="2" value="0.0" />\n')
        #meshFile.write('  </dirichlet>\n')
        meshFile.write('</poisson>\n')
        
        meshFile.write('<electrophysiology>\n')
        meshFile.write('  <stimuli number="1">\n')
        bb = cubit.get_bounding_box("surface", surfID)
        x0 = bb[0]
        x1 = 0.8*x0 + 0.2*bb[1]
        y0 = bb[3]
        y1 = 0.8*y0 + 0.2*bb[4]
        #z0 = bb[6]
        #z1 = 0.8*z0 + 0.2*bb[7]
        meshFile.write('    <stim start="0.00" duration="4.00" value="-35.7140" x0="%f" x1="%f" y0="%f" y1="%f" />\n' % (x0,x1,y0,y1))
        meshFile.write('  </stimuli>\n')
        meshFile.write('</electrophysiology>\n')
        meshFile.close()
        
    def saveMeshPoisson(self):
        nNodes = cubit.get_node_count()
        nTets = cubit.get_tet_count()
        meshFile = open(folder+"meshPoisson.xml", 'w')
        meshFile.write('<mesh celltype="tetrahedron" dim="3">\n')
        meshFile.write('  <nodes size="%d">\n' % (nNodes))
        for x in range(0, nNodes):
            coords = cubit.get_nodal_coordinates(x+1)
            meshFile.write('    <node id="%d" x="%f" y="%f" z="%f"/>\n' % (x,coords[0],coords[1],coords[2]))
        meshFile.write('  </nodes>\n')
        meshFile.write('  <elements size="%d">\n' % (nTets))
        for x in range(0, nTets):
            nodes = cubit.get_connectivity("tet", x+1)
            meshFile.write('    <element id="%d" v0="%d" v1="%d" v2="%d" v3="%d"/>\n' % (x,nodes[0]-1,nodes[1]-1,nodes[2]-1,nodes[3]-1))
        meshFile.write('  </elements>\n')
        meshFile.write('  <boundary celltype="triangle" dim="2">\n')
        
        surfs = cubit.get_relatives("volume", self.vol, "surface")
        ec = 0
        for x in range(0, 2):
            tris = cubit.get_surface_tris(surfs[x])
            for y in range(0, len(tris)):
                nodes = cubit.get_connectivity("tri", tris[y])
                element = [nodes[0]-1, nodes[1]-1, nodes[2]-1]
                ec = ec+1
                meshFile.write('    <element id="%d" marker="%d" v0="%d" v1="%d" v2="%d"/>\n' % (ec,x,element[0],element[1],element[2]))
        meshFile.write('    <element id="%d" marker="%d" v0="%d" v1="%d" v2="%d"/>\n' % (ec+1,2,element[0],element[1],element[2]))
        meshFile.write('  </boundary>\n')
        
        meshFile.write('</mesh>\n')
        meshFile.write('<poisson>\n')
        meshFile.write('  <neumann>\n')
        meshFile.write('    <node id="0" marker="0" value="-1.0" />\n')
        meshFile.write('    <node id="1" marker="1" value="1.0" />\n')
        meshFile.write('  </neumann>\n')
        meshFile.write('  <dirichlet>\n')
        meshFile.write('    <node id="2" marker="2" value="0.0" />\n')
        meshFile.write('  </dirichlet>\n')
        meshFile.write('</poisson>\n')
        meshFile.close()
        
    def genSurfaces(self):
        heap = [self.root]
        count = 0
        #self.root.tube.genSurfacesI()
        while not (len(heap) == 0):
            count += 1
            if count < 1500:
                atual = heap.pop()
            else:
                break
            circles = []
            circles.append(atual.circleI)
            end = False
            surfs = []
            while not end:
                if atual.son1 is not None and atual.son2 is not None:
                    end = True
                    circles.append(atual.circleF)
                elif atual.son1 is None and atual.son2 is None:
                    end = True
                    circles.append(atual.circleF)
                else:
                    if atual.root.son1 is not None and atual.root.son2 is not None:
                        atual.tube.genSurfaces()
                        surfs.extend(atual.tube.surfs)
                    circles.append(atual.circleF)
                    atual = atual.son1
            if len(circles) == 2:
                surfs.append(geo.createSurfaceCurve(circles[0]))
                surfs.append(geo.createSkinCurve(circles[0], circles[1]))
            if len(circles) > 2:
                circles.remove(circles[0])
                surfs.append(geo.createSkinCurve2(circles))
            if atual.son1 is None and atual.son2 is None:
                surfs.append(geo.createSurfaceCurve(circles[-1]))
            self.surfVols.extend(surfs)
            
            if atual.son1 is not None:
                heap.append(atual.son1)
            if atual.son2 is not None:
                heap.append(atual.son2) 
                
    def genVolumes(self):
#         for x in range(0, len(self.surfVols)):
# #             bodies = []
# #             for s in self.surfVols[x]:
# #                 bodies.append(cubit.get_owning_body("surface", s))
# #             geo.imprintBody(bodies) 
# #             geo.mergeBody(bodies) 
#             vol = geo.createVolume(self.surfVols[x])
#             if self.vol is None:
#                 self.vol = vol
#             else:
# #                 geo.imprintVolumes(self.vol, vol)
# #                 geo.mergeVolumes(self.vol, vol)
#                 self.vol = geo.unionVolumes(self.vol, vol)
        self.vol = geo.createVolumeFromAllSurfaces()
        free_curve_id_list = cubit.get_list_of_free_ref_entities("curve")
        for id in free_curve_id_list:
            geo.deleteCurve(id)
        free_vertex_id_list = cubit.get_list_of_free_ref_entities("vertex")
        for id in free_vertex_id_list:
            geo.deleteVertex(id)
    
    def genSons(self, branch):
        ipr = self.points[branch.initialPoint]
        fpr = self.points[branch.finalPoint]
        fp = self.points[branch.son1.finalPoint]
        fp2 = self.points[branch.son2.finalPoint]
        lenght = geo.distance(ipr, fpr)
        t=(lenght-branch.radius)/lenght
        ref = [(1-t)*ipr[0]+t*fpr[0],((1-t)*ipr[1]+t*fpr[1]),((1-t)*ipr[2]+t*fpr[2])]
        t2=(lenght+branch.radius)/lenght
        ref2 = [(1-t2)*ipr[0]+t2*fpr[0],((1-t2)*ipr[1]+t2*fpr[1]),((1-t2)*ipr[2]+t2*fpr[2])]
        u_0 = [fp[0]-fpr[0],fp[1]-fpr[1],fp[2]-fpr[2]]
        u_1 = [fp2[0]-fpr[0],fp2[1]-fpr[1],fp2[2]-fpr[2]]
        u_2 = [ref[0]-fp2[0],ref[1]-fp2[1],ref[2]-fp2[2]]
        u1 = geo2.crossProduct(u_0, u_1)
        u1 = geo2.getVector4_3(u1)
        geo2.normalizeVector(u1)
        u2 = geo2.crossProduct(u_0, u_2)
        u2 = geo2.getVector4_3(u2)
        geo2.normalizeVector(u2)
        coords = geo2.rotateByAxis(geo2.getVector4_3(ref), geo2.getVector4_3(fpr), u1, -90)
        vAux1 = geo.createVertex(coords[0],coords[1],coords[2])
        vAP = geo.createVertex(ref2[0],ref2[1],ref2[2])
#         
        circleAux = geo.createCircleNormal2(branch.finalPoint+1,vAux1,vAux1,branch.radius,branch.index)
        v1 = geo.createVertexOnCurveFraction(circleAux, 0.25)
        v2 = geo.createVertexOnCurveFraction(circleAux, 0.75)
        geo.deleteCurve(circleAux)
        arc1 = geo.createCircleNormal2(branch.finalPoint+1, v1, v2, branch.radius, branch.index)
        arc2 = geo.createCircleNormal2(branch.finalPoint+1, v2, v1, branch.radius, branch.index)
        vts1 = cubit.get_relatives("curve", arc1, "vertex")
        vts2 = cubit.get_relatives("curve", arc2, "vertex")
        l1 = geo.createLine(vts1[0], vts1[1])
        #l2 = geo.createLine(vts2[0], vts2[1])
        t=(lenght-branch.radius*1)/lenght
        u_3 = [cubit.vertex(vts1[0]).coordinates()[0]-cubit.vertex(vts1[1]).coordinates()[0],cubit.vertex(vts1[0]).coordinates()[1]-cubit.vertex(vts1[1]).coordinates()[1],cubit.vertex(vts1[0]).coordinates()[2]-cubit.vertex(vts1[1]).coordinates()[2]]
        u3 = geo2.getVector4_3(u_3)
        geo2.normalizeVector(u3)
        refE1 = [(1-t)*ipr[0]+t*fpr[0],((1-t)*ipr[1]+t*fpr[1]),((1-t)*ipr[2]+t*fpr[2])]
        coords = geo2.rotateByAxis(geo2.getVector4_3(refE1), geo2.getVector4_3(fpr), u3, 90)
        vE = geo.createVertex(coords[0],coords[1],coords[2])
        mid = geo.createVertexOnCurveFraction(l1, 0.5)
        eli = geo.createEllipseFull(vts1[0], vE, mid)
        eli1 = geo.splitCurve2(eli, vts1[0], vts1[1])
        eli2 = eli1-1
        
        e1 = geo.createEllipse(vts1[0], vAP, branch.finalPoint+1)
        e2 = geo.createEllipse(vts1[1], vAP, branch.finalPoint+1)
        e3 = geo.createCombineCurve([e1,e2])
        #e3 = geo.createSpline(vts1[0], vAP, vts1[1])
        
        circleI = branch.circleF
        
        uaux = cubit.vertex(v1).coordinates() 
        uaux2 = cubit.vertex(v2).coordinates() 
        u = [uaux2[0]-uaux[0],uaux2[1]-uaux[1],uaux2[2]-uaux[2]]
        u = geo2.getVector4_3(u)
        geo2.normalizeVector(u)
        ip = self.points[branch.son1.initialPoint]
        fp = self.points[branch.son1.finalPoint]
        fps = self.points[branch.son1.son1.finalPoint]
        lenght = geo.distance(ip, fp)
        t=(lenght-branch.son1.radius)/lenght
        ref = [(1-t)*ip[0]+t*fp[0],((1-t)*ip[1]+t*fp[1]),((1-t)*ip[2]+t*fp[2])]
        arcsf1 = self.genArc2(branch.son1, branch.son1.finalPoint+1,ref,u1,branch.son1.radius,1)
        branch.son1.circleI = circleI
        branch.son1.circleF = arcsf1[2]
        branch.son1.link = arcsf1[2]
        branch.son1.tube = Tube([arc1,e3], arcsf1)
        branch.son1.tube.line = l1
        
        ip = self.points[branch.son2.initialPoint]
        fp = self.points[branch.son2.finalPoint]
        fps = self.points[branch.son2.son1.finalPoint]
        lenght = geo.distance(ip, fp)
        t=(lenght-branch.son2.radius)/lenght
        ref = [(1-t)*ip[0]+t*fp[0],((1-t)*ip[1]+t*fp[1]),((1-t)*ip[2]+t*fp[2])]
        arcsf2 = self.genArc2(branch.son2, branch.son2.finalPoint+1,ref,u1,branch.son2.radius,2)
        branch.son2.circleI = circleI
        branch.son2.circleF = arcsf2[2]
        branch.son2.link = arcsf2[2]
        branch.son2.tube = Tube([arc2,e3], arcsf2)
        branch.son2.tube.line = l1
                
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
            
        t=(lenght-1.1*branch.radius)/lenght
        newPf = [(1-t)*pi[0]+t*pf[0],((1-t)*pi[1]+t*pf[1]),((1-t)*pi[2]+t*pf[2])]
        self.points[branch.finalPoint] = newPf
        #geo.createVertex(newPf[0],newPf[1],newPf[2])
        t=(lenght+0*branch.radius)/lenght
        newPoint1 = [(1-t)*pi[0]+t*pf[0],((1-t)*pi[1]+t*pf[1]),((1-t)*pi[2]+t*pf[2])]
#         u0 = [newPf[0]-newPoint1[0],newPf[1]-newPoint1[1],newPf[2]-newPoint1[2]]
#         u1 = [newPoint1[0]-pf1[0],newPoint1[1]-pf1[1],newPoint1[2]-pf1[2]]
#         u0 = [newPoint1[0]-newPf[0],newPoint1[1]-newPf[1],newPoint1[2]-newPf[2]]
#         u1 = [pf1[0]-newPf[0],pf1[1]-newPf[1],pf1[2]-newPf[2]]
        u0 = [pf1[0]-newPf[0],pf1[1]-newPf[1],pf1[2]-newPf[2]]
        u1 = [pf2[0]-newPf[0],pf2[1]-newPf[1],pf2[2]-newPf[2]]
        u = geo2.crossProduct(u0, u1)
        u = geo2.getVector4_3(u)
        geo2.normalizeVector(u)
        newPoint1 = geo2.rotateByAxis(geo2.getVector4_3(newPoint1), geo2.getVector4_3(newPf), u, 40)
        vecR = geo2.getVector(pi, newPf)
        vecF1 = geo2.getVector(newPf, newPoint1)
        #geo.createVertex(newPoint1[0],newPoint1[1],newPoint1[2])
        t=(lenght+0*branch.radius)/lenght
        newPoint2 = [(1-t)*pi[0]+t*pf[0],((1-t)*pi[1]+t*pf[1]),((1-t)*pi[2]+t*pf[2])]
#         u0 = [newPoint2[0]-newPf[0],newPoint2[1]-newPf[1],newPoint2[2]-newPf[2]]
#         u1 = [pf2[0]-newPf[0],pf2[1]-newPf[1],pf2[2]-newPf[2]]
#         u = geo2.crossProduct(u0, u1)
#         u = geo2.getVector4_3(u)
#         geo2.normalizeVector(u)
        newPoint2 = geo2.rotateByAxis(geo2.getVector4_3(newPoint2), geo2.getVector4_3(newPf), u, -40)
        vecF2 = geo2.getVector(newPf, newPoint2)
        ang3 = geo2.getAngle(vecF1, vecF2)
        if ang3 < 60:
            newPoint1 = geo2.rotateByAxis(geo2.getVector4_3(newPoint1), geo2.getVector4_3(newPf), u, 20)
            newPoint2 = geo2.rotateByAxis(geo2.getVector4_3(newPoint2), geo2.getVector4_3(newPf), u, -20)
        #geo.createVertex(newPoint2[0],newPoint2[1],newPoint2[2])
        d1 = geo.distance(pf1, newPoint1)+geo.distance(pf2, newPoint2)
        d2 = geo.distance(pf1, newPoint2)+geo.distance(pf2, newPoint1)
        if d1>d2:
            aux = newPoint1
            newPoint1 = newPoint2
            newPoint2 = aux
        self.points.append(newPoint1)
        radius1 = branch.radius
        radius2 = branch.radius
#         radius1 = branch.son1.radius
#         radius2 = branch.son2.radius
#         sumRadius = radius1 + radius2
#         while sumRadius > 1.2:
#             radius1 = radius1*0.95
#             radius2 = radius2*0.95
#             sumRadius = radius1 + radius2
        new1 = Branch(branch.finalPoint, len(self.points)-1,radius1)
        new1.root = branch
        new1.son1 = branch.son1
        new1.son1.initialPoint = len(self.points)-1
        new1.son1.root = new1
        branch.son1 = new1
           
        self.points.append(newPoint2)
        new2 = Branch(branch.finalPoint, len(self.points)-1, radius2)
        new2.root = branch
        new2.son1 = branch.son2
        new2.son1.initialPoint = len(self.points)-1
        new2.son1.root = new2
        branch.son2 = new2
    
    def genArc(self, branch, point, ref, u, radius):
        ipr = self.points[point-1]
        coords = geo2.rotateByAxis(geo2.getVector4_3(ref), geo2.getVector4_3(ipr), u, -90)
        vAuxI1 = geo.createVertex(coords[0],coords[1],coords[2])
        circleI = geo.createCircleNormal2(point,vAuxI1,vAuxI1,radius,branch.index)
        vAuxI2 = geo.createVertexOnCurveFraction(circleI, 0.25)
        vAuxI3 = geo.createVertexOnCurveFraction(circleI, 0.5)
        vAuxI4 = geo.createVertexOnCurveFraction(circleI, 0.75)
        #arc = geo.createCircleNormal2(point, vAuxI1, vAuxI2, radius, branch.index)
        #arcI1 = [vAuxI1, arc]
        arc = geo.createCircleNormal2(point, vAuxI2, vAuxI4, radius, branch.index)
        arcI2 = [vAuxI2, arc]
        #arc = geo.createCircleNormal2(point, vAuxI3, vAuxI4, radius, branch.index)
        #arcI3 = [vAuxI3, arc]
        arc = geo.createCircleNormal2(point, vAuxI4, vAuxI2, radius, branch.index)
        arcI4 = [vAuxI4, arc]
        geo.deleteCurve2(circleI)

        #arcs = [arcI1,arcI2,arcI3,arcI4]
        arcs = [arcI2,arcI4]
        return arcs
    
    def genArc2(self, branch, point, ref, u, radius, son):
        ipr = self.points[branch.initialPoint]
        fpr = self.points[branch.finalPoint]
        lenght = geo.distance(ipr, fpr)
        if son == 1:
            off = 40
        else:
            off = -40
        off = 0
         
        t=(lenght-branch.radius*0.6)/lenght
        ref2 = [(1-t)*ipr[0]+t*fpr[0],((1-t)*ipr[1]+t*fpr[1]),((1-t)*ipr[2]+t*fpr[2])]
        ipr = self.points[point-1]
        coords = geo2.rotateByAxis(geo2.getVector4_3(ref2), geo2.getVector4_3(ipr), u, -90+off)
        vAux1 = geo.createVertex(coords[0],coords[1],coords[2])
            
        ipr = self.points[point-1]
        coords = geo2.rotateByAxis(geo2.getVector4_3(ref), geo2.getVector4_3(ipr), u, -90+off)
        vAux2 = geo.createVertex(coords[0],coords[1],coords[2])
        circleI = geo.createCircleNormal2(point,vAux2,vAux2,radius,branch.index)
        
        vAux3 = geo.createVertexOnCurveFraction(circleI, 0.25)
        eliFullAux = geo.createEllipseFull(vAux3, vAux1, point)
        eliFull = geo.createEllipseFull(vAux3, vAux1, point)
        
        geo.deleteCurve2(circleI)
        geo.deleteVertex(vAux1)
        geo.deleteVertex(vAux2)
         
        v1 = geo.createVertexOnCurveFraction(eliFullAux, 0)
        v2 = geo.createVertexOnCurveFraction(eliFullAux, 0.5)
        eli1 = geo.splitCurve2(eliFullAux, v1, v2)
        eli2 = eli1-1
        arcs = [eli1,eli2,eliFull]
        return arcs
    
    def genArcFull(self, branch, point, ref):
        radius = branch.getRadius()
        u = branch.getU(self.points)
        ipr = self.points[branch.initialPoint]
        fpr = self.points[branch.finalPoint]
        lenght = geo.distance(ipr, fpr)
         
        t=(lenght-radius[0]*radius[1])/lenght
        ref2 = [(1-t)*ipr[0]+t*fpr[0],((1-t)*ipr[1]+t*fpr[1]),((1-t)*ipr[2]+t*fpr[2])]
        ipr = self.points[point-1]
        coords = geo2.rotateByAxis(geo2.getVector4_3(ref2), geo2.getVector4_3(ipr), u, -90)
        vAux1 = geo.createVertex(coords[0],coords[1],coords[2])
            
        
        coords = geo2.rotateByAxis(geo2.getVector4_3(ref), geo2.getVector4_3(ipr), u, -90)
        vAux2 = geo.createVertex(coords[0],coords[1],coords[2])
        circleI = geo.createCircleNormal2(point,vAux2,vAux2,radius[1],branch.index)
        
        vAux3 = geo.createVertexOnCurveFraction(circleI, 0.25)
        eliFull = geo.createEllipseFull(vAux3, vAux1, point)
        geo.deleteCurve(circleI)
        #geo.deleteVertex(vAux2)
        geo.deleteVertex(vAux1)
        #geo.deleteVertex(vAux3)
        if eliFull == circleI:
            tube.eprint("Erro Ellipse")
            return -1
        else:
            return eliFull

    def splitBranch(self,branch, point):
        self.points.append(point)
        radius = branch.radius
        new = Branch(len(self.points)-1, branch.finalPoint, radius)
        new.root = branch
        new.son1 = branch.son1
        new.son2 = branch.son2
        if branch.son1 is not None:
            new.son1.root = new
        if branch.son2 is not None:
            new.son2.root = new
            branch.son2 = None 
        branch.finalPoint = new.initialPoint
        branch.son1 = new
    
    def check(self):
        heap = [self.root]
        while not (len(heap) == 0):
            atual = heap.pop()
            print atual.initialPoint
            if atual.son1 is not None and atual.son2 is not None:
                print (atual.son1.radius/atual.radius+atual.son2.radius/atual.radius)/2
                heap.append(atual.son1)
                heap.append(atual.son2)  
            elif atual.son1 is not None:
                print atual.son1.radius/atual.radius
                heap.append(atual.son1) 