import geo
from vessel import Vessel
import cubit
import math
import random
from settings import *

class Branch(object):

    def __init__(self, initialPoint, finalPoint, radius):
        self.initialPoint = initialPoint
        self.finalPoint = finalPoint
        self.radius = radius
        self.son1 = None
        self.son2 = None
        self.root = None
    
