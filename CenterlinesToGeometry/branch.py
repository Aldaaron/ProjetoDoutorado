import geo
import cubit
import math
import random
import geo2
from settings import *

class Branch(object):

    def __init__(self, initialPoint, finalPoint, radius):
        self.initialPoint = initialPoint
        self.finalPoint = finalPoint
        self.radius = radius
        self.son1 = None
        self.son2 = None
        self.root = None
        self.tube = None
        self.circleI = 0
        self.circleF = 0
        self.nSplit = 0
        self.rootRadius = radius
        self.u = None
        self.surfs = []
        self.delete = False
        
    def getRadius(self):
        smallR = 0.6
        if self.nSplit == 1:
            smallR = 0.7
        if self.nSplit == 2:
            smallR = 0.8 
        if self.nSplit == 3:
            smallR = 0.85 
        if self.nSplit == 4:
            smallR = 0.9 
        if self.nSplit == 5:
            smallR = 0.95 
        if self.nSplit == 6:
            smallR = 1 
        if self.nSplit == 7:
            smallR = 1
        r = (self.radius*self.nSplit + self.rootRadius*(7-self.nSplit))/7 
        return [smallR, r]
    
    def getU(self, points):
        fpr = points[self.root.initialPoint]
        fpi = points[self.initialPoint]
        fpf = points[self.finalPoint]
        u_0 = [fpr[0]-fpi[0],fpr[1]-fpi[1],fpr[2]-fpi[2]]
        u_1 = [fpf[0]-fpi[0],fpf[1]-fpi[1],fpf[2]-fpi[2]]
        u = geo2.crossProduct(u_0, u_1)
        u = geo2.getVector4_3(u)
        geo2.normalizeVector(u)
        return u
        
    
