#!/usr/bin/python
import math

def printMatrix (matrix):
   """print matrix. self-explanatory"""

   tot = matrix[0][0]

   print round(matrix[0][0],5),'\t',round(matrix[0][1],5),'\t',round(matrix[0][2],5)
   print round(matrix[1][0],5),'\t',round(matrix[1][1],5),'\t',round(matrix[1][2],5)
   print round(matrix[2][0],5),'\t',round(matrix[2][1],5),'\t',round(matrix[2][2],5)
   print '\n\n sqrt(s00)=', round(math.sqrt(matrix[0][0]),5), 'sqrt(sLL)=', round(math.sqrt(matrix[1][1]),5), 'sqrt(sRR)=', round(math.sqrt(matrix[2][2]),5)

infile = open('ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad/SystematicOutput_el_mu_lephad_bTag.script2.txt','r')

m = 3
corrSystMatrix = [[0 for x in range(m)] for y in range(m)]
corrStatMatrix = [[0 for x in range(m)] for y in range(m)]
corrMatrix = [[0 for x in range(m)] for y in range(m)]

for line in infile:
    name= line.split('&')[0].split('\t')[0]
    d0= float(line.split('&')[1].split('\t')[0])
    dL= float(line.split('&')[2].split('\t')[0])
    dR= float(line.split('&')[3].split('\\')[0])
    
    s00 = sLL = sRR = s0L = s0R = sLR = 0

    #if 'BTAG_bTag' in name:
    print 'name:', name, 'd0:', d0, 'dL:', dL, 'dR:',dR
    s00 = d0*d0
    sLL = dL*dL
    sRR = dR*dR
    s0L = d0*dL
    s0R = d0*dR
    sLR = dL*dR
        
    corrSystMatrix[0][0] = corrSystMatrix[0][0] + s00
    corrSystMatrix[1][1] = corrSystMatrix[1][1] + sLL
    corrSystMatrix[2][2] = corrSystMatrix[2][2] + sRR
    corrSystMatrix[0][1] = corrSystMatrix[0][1] + s0L
    corrSystMatrix[1][0] = corrSystMatrix[1][0] + s0L
    corrSystMatrix[0][2] = corrSystMatrix[0][2] + s0R
    corrSystMatrix[2][0] = corrSystMatrix[2][0] + s0R
    corrSystMatrix[1][2] = corrSystMatrix[1][2] + sLR
    corrSystMatrix[2][1] = corrSystMatrix[2][1] + sLR

stat0 = 0.007
statL = 0.005
statR = 0.004
corrSystMatrix[0][0] = corrSystMatrix[0][0] + stat0*stat0
corrSystMatrix[1][1] = corrSystMatrix[1][1] + statL*statL
corrSystMatrix[2][2] = corrSystMatrix[2][2] + statR*statR
corrSystMatrix[0][1] = corrSystMatrix[0][1] + stat0*statL
corrSystMatrix[1][0] = corrSystMatrix[1][0] + statL*stat0
corrSystMatrix[0][2] = corrSystMatrix[0][2] + stat0*statR
corrSystMatrix[2][0] = corrSystMatrix[2][0] + statR*stat0
corrSystMatrix[1][2] = corrSystMatrix[1][2] + statL*statR
corrSystMatrix[2][1] = corrSystMatrix[2][1] + statR*statL


printMatrix(corrSystMatrix)
