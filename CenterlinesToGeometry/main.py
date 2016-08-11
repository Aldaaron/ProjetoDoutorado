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
    fp = open(folder+"pontos3.txt")
    pontos = [ map(float,line.split('\t')) for line in fp ]
    fl = open(folder+"linhas3.txt")
    linhas = [ map(float,line.split('\t')) for line in fl ]
    arvore = Arvore(pontos, linhas)
#     arvore.fixSizes()
#     arvore.draw()
#    arvore.imprime()
#     arvore.split()
#     arvore.smoth()
#     arvore.save2()
    arvore.smothRadius()
    arvore.makeGeometry()
    arvore.mesh()
    arvore.saveMesh()
#    arvore.check()
    fp.close()
    fl.close()
    
main()