import math
import sys

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

# def invertMatrix(dest, src):
#     b = getIdentity()
#     index = [0.0,0.0,0.0,0.0]
#     for i in range(0, 4):
#         index[i] = 0
#         for j in range(0, 4):
#             if (i==j):
#                 b[i][j] = 1
#             else:
#                 b[i][j] = 0
#         dest[i][j] = 0
# #     gaussian(src, index);
# 
# #     for i in range(0, 3):
# #         for j in range(i+1, 4):
# #             for k in range(0, 4):
# #                 b[index[j]][k] -= src[index[j]][i]*b[index[i]][k]
# # 
# #     for i in range(0, 4):
# #         dest[4-1][i] = b[index[4-1]][i]/src[index[4-1]][4-1]
# #         for j in range(2, -1,-1):
# #             dest[j][i] = b[index[j]][i]
# #             for k in range(j+1, 4):
# #                 dest[j][i] -= src[index[j]][k]*dest[k][i]
# #             dest[j][i] /= src[index[j]][j]
# #     return dest
# 
# def gaussian(a,index):
#     c = [0.0,0.0,0.0,0.0]
#     for i in range(0, 4):
#         index[i] = i
#         c[i] = 0
# 
#     for i in range(0, 4):
#         c1 = 0
#         for j in range(0, 4):
#             if (a[i][j] < 0):
#                 c0 = a[i][j]*-1.0
#             else:
#                 c0 = a[i][j]
#             if (c0 > c1):
#                 c1 = c0
#         c[i] = c1
#     k = 0
#     for j in range(0, 3):
#         pi1 = 0
#         for i in range(j, 4):
#             if(a[index[i]][j] < 0):
#                 pi0 = a[index[i]][j]*-1.0
#             else:
#                 pi0 = a[index[i]][j]
#             pi0 /= c[index[i]]
#             if (pi0 > pi1):
#                 pi1 = pi0
#                 k = i
#         itmp = index[j]
#         index[j] = index[k]
#         index[k] = itmp
#         for i in range(j+1, 4):
#             pj = a[index[i]][j]/a[index[j]][j]
#             a[index[i]][j] = pj
#             for l in range(j+1, 4):
#                 a[index[i]][l] -= pj*a[index[j]][l]
#                 
def inv(M):
    """
    return the inv of the matrix M
    """
    #clone the matrix and append the identity matrix
    # [int(i==j) for j in range_M] is nothing but the i(th row of the identity matrix
    m2 = [row[:]+[int(i==j) for j in range(len(M) )] for i,row in enumerate(M) ]
    # extract the appended matrix (kind of m2[m:,...]
    return [row[len(M[0]):] for row in m2] if gauss_jordan(m2) else None
  
def gauss_jordan(m, eps = 1.0/(10**10)):
    """Puts given matrix (2D array) into the Reduced Row Echelon Form.
       Returns True if successful, False if 'm' is singular.
       NOTE: make sure all the matrix items support fractions! Int matrix will NOT work!
       Written by Jarno Elonen in April 2005, released into Public Domain"""
    (h, w) = (len(m), len(m[0]))
    for y in range(0,h):
        maxrow = y
        for y2 in range(y+1, h):    # Find max pivot
            if abs(m[y2][y]) > abs(m[maxrow][y]):
                maxrow = y2
        (m[y], m[maxrow]) = (m[maxrow], m[y])
        if abs(m[y][y]) <= eps:     # Singular?
            return False
        for y2 in range(y+1, h):    # Eliminate column y
            c = m[y2][y] / m[y][y]
            for x in range(y, w):
                m[y2][x] -= m[y][x] * c
    for y in range(h-1, 0-1, -1): # Backsubstitute
        c  = m[y][y]
        for y2 in range(0,y):
            for x in range(w-1, y-1, -1):
                m[y2][x] -=  m[y][x] * m[y2][y] / c
        m[y][y] /= c
        for x in range(h, w):       # Normalize row y
            m[y][x] /= c
    return True

def rot(x,y,z,a,b,c,u,v,w,t):
    t = (t*math.pi)/180.0
    vet = [0.0,0.0,0.0]
    vet[0] = (a*(v*v+w*w) - u*(b*v + c*w - u*x - v*y - w*z))*(1-math.cos(t))+x*math.cos(t)+(-c*v+b*w-w*y+v*z)*math.sin(t)
    vet[1] = (b*(u*u+w*w) - v*(a*u + c*w - u*x - v*y - w*z))*(1-math.cos(t))+y*math.cos(t)+(c*u-a*w+w*x-u*z)*math.sin(t)
    vet[2] = (c*(u*u+v*v) - w*(a*u + b*v - u*x - v*y - w*z))*(1-math.cos(t))+z*math.cos(t)+(-b*u+a*v-v*x+u*y)*math.sin(t)
    return vet
    
def crossProduct(src1, src2):
    dest = [0.0,0.0,0.0]
    dest[0] = (src1[1]*src2[2])-(src2[1]*src1[2])
    dest[1] = (src1[2]*src2[0])-(src2[2]*src1[0])
    dest[2] = (src1[0]*src2[1])-(src2[0]*src1[1])
    return dest