import cubit
import sys
sys.path.append('/home/bruno/git/local/CenterlinesToGeometry')
from settings import *
from arvore import Arvore
import geo
import geo2

def main():
    cubit.cmd("set warning off")
    cubit.cmd("set developer commands on")
    fp = open(folder+"pontos.txt")
    pontos = [ map(float,line.split('\t')) for line in fp ]
    fl = open(folder+"linhas.txt")
    linhas = [ map(float,line.split('\t')) for line in fl ]
    arvore = Arvore(pontos, linhas)
    arvore.preProcess()
    #arvore.draw()
    arvore.smothRadius()
    arvore.makeGeometry()
    arvore.mesh()
    arvore.saveMesh()
    #arvore.saveMesh3D()
    fp.close()
    fl.close()
    
main()