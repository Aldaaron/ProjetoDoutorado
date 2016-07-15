import sys
sys.path.append('/home/bruno/git/local/CenterlinesToGeometry')
import cubit
from tree import Tree
from settings import *
from arvore import Arvore
    
def main():
    cubit.cmd('set geometry engine acis')
    cubit.cmd('set developer commands on')
    if not cubitPrint:
        cubit.cmd('set echo off') # mostra informacoes
    fp = open(folder+"pontos2.txt")
    pontos = [ map(float,line.split('\t')) for line in fp ]
    fl = open(folder+"linhas2.txt")
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
    #arvore.draw()
    #arvore.split()
    #arvore.smoth()
    #arvore.save()
    arvore.makeGeometry()
    
main()