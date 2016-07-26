import cubit
import sys
sys.path.append('/home/bruno/git/local/CenterlinesToGeometry')
from tree import Tree
from settings import *
from arvore import Arvore
import geo
import geo2

def main():
    cubit.cmd("set warning off")
    fp = open(folder+"pontos3.txt")
    pontos = [ map(float,line.split('\t')) for line in fp ]
    fl = open(folder+"linhas3.txt")
    linhas = [ map(float,line.split('\t')) for line in fl ]
#     for l in linhas:
#         l[2] = 0.5
#     arvore = Tree(pontos, linhas)
#     #arvore.genPoints()
#     arvore.preProcess()
#     arvore.genVertices(scale)
#     arvore.genCenterlines()
#     arvore.genSurfaces()
#     arvore.union()     
    #arvore.clean() 
    arvore = Arvore(pontos, linhas)
#    arvore.fixSizes()
#    arvore.draw()
#    arvore.imprime()
#    arvore.split2()
#    arvore.split()
#    arvore.smoth()
#    arvore.smoth2()
#    arvore.save2()
    arvore.makeGeometry2()
#    arvore.check()
    fp.close()
    fl.close()
    
main()