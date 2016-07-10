import sys
sys.path.append('/home/bruno/git/local/CenterlinesToGeometry')
import cubit
from tree import Tree
from settings import *
    
def main():
    cubit.cmd('set geometry engine acis')
    cubit.cmd('set developer commands on')
    if not cubitPrint:
        cubit.cmd('set echo off') # mostra informacoes
    fp = open(folder+"pontos.txt")
    pontos = [ map(float,line.split('\t')) for line in fp ]
    fl = open(folder+"linhas.txt")
    linhas = [ map(float,line.split('\t')) for line in fl ]
    arvore = Tree(pontos, linhas)
    arvore.preProcess()
    arvore.genVertices(scale)
    arvore.genCenterlines()
    arvore.genVolume()
    arvore.union()
    arvore.clean()

main()