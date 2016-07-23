import math
import sys
import geo

#Carrega a matriz identidade em "matrix"
def getIdentity():
    mat = [[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,1.0]]
    return mat

#Copia o conteudo do vetor3d "src" para o vetor4d "dest",
#atribuindo 1 para o valor da ultima coordenada.
def getVector4_3(src):
    dest = [src[0],src[1],src[2],1.0]
    return dest

#Multiplica duas matrizes de 3 dimensoes
def multMat(src1, src2):
    dest = getIdentity()
    dest[0][0] = src1[0][0] * src2[0][0] + src1[0][1] * src2[1][0] + src1[0][2] * src2[2][0] + src1[0][3] * src2[3][0]
    dest[0][1] = src1[0][0] * src2[0][1] + src1[0][1] * src2[1][1] + src1[0][2] * src2[2][1] + src1[0][3] * src2[3][1]
    dest[0][2] = src1[0][0] * src2[0][2] + src1[0][1] * src2[1][2] + src1[0][2] * src2[2][2] + src1[0][3] * src2[3][2]
    dest[0][3] = src1[0][0] * src2[0][3] + src1[0][1] * src2[1][3] + src1[0][2] * src2[2][3] + src1[0][3] * src2[3][3]
 
    dest[1][0] = src1[1][0] * src2[0][0] + src1[1][1] * src2[1][0] + src1[1][2] * src2[2][0] + src1[1][3] * src2[3][0]
    dest[1][1] = src1[1][0] * src2[0][1] + src1[1][1] * src2[1][1] + src1[1][2] * src2[2][1] + src1[1][3] * src2[3][1]
    dest[1][2] = src1[1][0] * src2[0][2] + src1[1][1] * src2[1][2] + src1[1][2] * src2[2][2] + src1[1][3] * src2[3][2]
    dest[1][3] = src1[1][0] * src2[0][3] + src1[1][1] * src2[1][3] + src1[1][2] * src2[2][3] + src1[1][3] * src2[3][3]
 
    dest[2][0] = src1[2][0] * src2[0][0] + src1[2][1] * src2[1][0] + src1[2][2] * src2[2][0] + src1[2][3] * src2[3][0]
    dest[2][1] = src1[2][0] * src2[0][1] + src1[2][1] * src2[1][1] + src1[2][2] * src2[2][1] + src1[2][3] * src2[3][1]
    dest[2][2] = src1[2][0] * src2[0][2] + src1[2][1] * src2[1][2] + src1[2][2] * src2[2][2] + src1[2][3] * src2[3][2]
    dest[2][3] = src1[2][0] * src2[0][3] + src1[2][1] * src2[1][3] + src1[2][2] * src2[2][3] + src1[2][3] * src2[3][3]
 
    dest[3][0] = src1[3][0] * src2[0][0] + src1[3][1] * src2[1][0] + src1[3][2] * src2[2][0] + src1[3][3] * src2[3][0]
    dest[3][1] = src1[3][0] * src2[0][1] + src1[3][1] * src2[1][1] + src1[3][2] * src2[2][1] + src1[3][3] * src2[3][1]
    dest[3][2] = src1[3][0] * src2[0][2] + src1[3][1] * src2[1][2] + src1[3][2] * src2[2][2] + src1[3][3] * src2[3][2]
    dest[3][3] = src1[3][0] * src2[0][3] + src1[3][1] * src2[1][3] + src1[3][2] * src2[2][3] + src1[3][3] * src2[3][3]
    
    return dest 
#Agrega uma matriz de translacao com os parametros tx(translacao em x),
#ty(translacao em y), tz(translacao em z) na matriz do indice passado.
def translate(mat, tx, ty, tz):
    t = [[1.0,0.0,0.0,tx],[0.0,1.0,0.0,ty],[0.0,0.0,1.0,tz],[0.0,0.0,0.0,1.0]]
    return multMat(mat, t)
 
#Efetua a multiplicacao de uma matrix4d por um vetor4d
def multMatVec(mat, vec):
    dest = [0.0,0.0,0.0,0.0]
    dest[0] = mat[0][0] * vec[0] + mat[0][1] * vec[1] + mat[0][2] * vec[2] + mat[0][3] * vec[3]
    dest[1] = mat[1][0] * vec[0] + mat[1][1] * vec[1] + mat[1][2] * vec[2] + mat[1][3] * vec[3]
    dest[2] = mat[2][0] * vec[0] + mat[2][1] * vec[1] + mat[2][2] * vec[2] + mat[2][3] * vec[3]
    dest[3] = mat[3][0] * vec[0] + mat[3][1] * vec[1] + mat[3][2] * vec[2] + mat[3][3] * vec[3]
    return dest
    
#Agrega uma matriz de rotacao em torno de X com angulo theta passado em graus na matriz do indice passado.
def rotateX(mat, theta):
    theta = (theta*math.pi)/180.0
    rotateX = [[1.0,0.0,0.0,0.0],
               [0.0,math.cos(theta),-math.sin(theta),0.0],
               [0.0,math.sin(theta),math.cos(theta),0.0],
               [0.0,0.0,0.0,1.0]]
    return multMat(mat, rotateX)
    
#Agrega uma matriz de rotacao em torno de Y com angulo theta passado em graus na matriz do indice passado.
def rotateY(mat, theta):
    theta = (theta*math.pi)/180.0
    rotateY = [[math.cos(theta),0.0,math.sin(theta),0.0],
               [0.0,1.0,0.0,0.0],
               [-math.sin(theta),0.0,math.cos(theta),0.0],
               [0.0,0.0,0.0,1.0]]
    return multMat(mat, rotateY)
    
#Agrega uma matriz de rotacao em torno de X com angulo theta passado em graus na matriz do indice passado.
def rotateZ(mat, theta):
    theta = (theta*math.pi)/180.0
    rotateZ = [[math.cos(theta),-math.sin(theta),0.0,0.0],
               [math.sin(theta),math.cos(theta),0.0,0.0],
               [0.0,0.0,1.0,0.0],
               [0.0,0.0,0.0,1.0]]
    return multMat(mat, rotateZ)
    
def rotateByAxis(point, origin, direc, theta):
    a = direc[0]
    b = direc[1]
    c = direc[2]
    d = math.sqrt(b*b+c*c)
    rotX = [[1.0,0.0,0.0,0.0],
               [0.0,c/d,-b/d,0.0],
               [0.0,b/d,c/d,0.0],
               [0.0,0.0,0.0,1.0]]
    rotXi = [[1.0,0.0,0.0,0.0],
               [0.0,c/d,b/d,0.0],
               [0.0,-b/d,c/d,0.0],
               [0.0,0.0,0.0,1.0]]
    rotY = [[d,0.0,-a,0.0],
               [0.0,1.0,0.0,0.0],
               [a,0.0,d,0.0],
               [0.0,0.0,0.0,1.0]]
    rotYi = [[d,0.0,a,0.0],
               [0.0,1.0,0.0,0.0],
               [-a,0.0,d,0.0],
               [0.0,0.0,0.0,1.0]]
    
    mat = translate(getIdentity(), origin[0], origin[1], origin[2])
    mat = multMat(mat, rotXi)
    mat = multMat(mat, rotYi)
    mat = rotateZ(mat, theta)
    mat = multMat(mat, rotY)
    mat = multMat(mat, rotX)
    mat = translate(mat, -origin[0], -origin[1], -origin[2])
    dest = multMatVec(mat, point)
    return [dest[0],dest[1],dest[2]]
    
#Calcula a normalizacao do vetor4d "src" em "dest"
def normalizeVector(vec):
    norm = normVector(vec)
    vec[0] = vec[0]/norm
    vec[1] = vec[1]/norm
    vec[2] = vec[2]/norm
    
def normVector(vect):
    return math.sqrt((vect[0]*vect[0])+(vect[1]*vect[1])+(vect[2]*vect[2]))

def getVector(v1,v2):
    dest = [v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]]
    return dest

def dotProduct(v1,v2):
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

def crossProduct(src1, src2):
    dest = [0.0,0.0,0.0]
    dest[0] = (src1[1]*src2[2])-(src2[1]*src1[2])
    dest[1] = (src1[2]*src2[0])-(src2[2]*src1[0])
    dest[2] = (src1[0]*src2[1])-(src2[0]*src1[1])
    return dest

def getAngle(vec1, vec2):
    dot = dotProduct(vec1, vec2)
    n1 = normVector(vec1)
    n2 = normVector(vec2)
    cosA = dot/(n1*n2)
    if cosA > 1.0:
        cosA = 1.0
    if cosA < -1.0:
        cosA = -1.0
    angle = math.acos(cosA)
    return abs(geo.radiansToDegree(angle))

def getAngleUns(vec1, vec2):
    dot = dotProduct(vec1, vec2)
    n1 = normVector(vec1)
    n2 = normVector(vec2)
    cosA = dot/(n1*n2)
    if cosA > 1.0:
        cosA = 1.0
    if cosA < -1.0:
        cosA = -1.0
    angle = math.acos(cosA)
    return geo.radiansToDegree(angle)