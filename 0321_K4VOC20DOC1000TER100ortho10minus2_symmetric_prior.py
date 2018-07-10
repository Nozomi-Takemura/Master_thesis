import os
os.getcwd()
import numpy as np
import scipy
import pandas as pd
import re
 ##Retrieved from http://www.hongliangjie.com/2010/09/30/generate-synthetic-data-for-lda/
#X. Yi, L. Hong,
## This program is to generated synthetic data for LDA
import math
import random
import numpy
import numpy.random
import sys
from scipy import io
import os
## define some constant
lp = [1000,100,10]
ortho = [0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001]
TOPIC_N = 4
VOCABULARY_SIZE = 10
#FILE_NAME = sys.argv[1]
#FILE_NAME = "test"
#for dn in lp:
#for tm in lp:
dn = 1000
tm = 100 
#TOPIC_N = 2
#VOCABULARY_SIZE = 20
DOC_NUM = dn
TERM_PER_DOC = tm
expe_count = 0
for i in range(1,101):
    expe_count = expe_count + 1
    ortho_count = 0
    #beta = [0.01 for i in range(VOCABULARY_SIZE)]
    pri = []
    prior = 0.05
    for i in range(VOCABULARY_SIZE):
        prior = prior 
        pri.append(prior)	
    beta = pri
    #alpha = [0.9 for i in range(TOPIC_N)]
    #alpha = [0.5 for i in range(TOPIC_N)]
    alpha = [0.25, 0.25, 0.25, 0.25]
    FILE_NAME = sys.argv[1]
    #for p in ortho:
    p = 0.01
    #ortho_count = ortho_count + 1
    ortho_count = 2
    ## generate multinomial distribution over words for each topic
    a = 0
    zero = p
    zero_out = zero*10000000
    os.system ("mkdir " + "K2ND" + str(dn) + "TE" + str(tm) + "VOC" + str(VOCABULARY_SIZE) + "ortho_minus" + str(ortho_count) + "index_expe" + str(expe_count))
    path = "./K2ND" + str(dn) + "TE" + str(tm) + "VOC" + str(VOCABULARY_SIZE) + "ortho_minus" + str(ortho_count) + "index_expe" + str(expe_count) + "/"
    while (a == 0):
        phi = []
        for i in range(TOPIC_N):
            topic = numpy.random.mtrand.dirichlet(beta, size = 1)
            phi.append(topic)
        p1 = np.squeeze(np.asarray(phi[0]))
        p2 = np.squeeze(np.asarray(phi[1]))
        #a = 1
        dp = np.dot(p1,p2)
        if dp < zero:
            a = 1
    ## generate words for each document
    output_f = open(str(path) + FILE_NAME+'.doc','w')
    z_f = open(str(path) + FILE_NAME+'.z','w')
    theta_f = open(str(path) + FILE_NAME+'.theta','w')
    for i in range(DOC_NUM):
        buffer = {}
        z_buffer = {} ## keep track the true z
        ## first sample theta
        theta = numpy.random.mtrand.dirichlet(alpha,size = 1)
        for j in range(TERM_PER_DOC):
            ## first sample z
            z = numpy.random.multinomial(1,theta[0],size = 1)
            z_assignment = 0
            for k in range(TOPIC_N):
                if z[0][k] == 1:
                    break
                z_assignment += 1
            if not z_assignment in z_buffer:
                z_buffer[z_assignment] = 0
            z_buffer[z_assignment] = z_buffer[z_assignment] + 1
            ## sample a word from topic z
            w = numpy.random.multinomial(1,phi[z_assignment][0],size = 1)
            w_assignment = 0
            for k in range(VOCABULARY_SIZE):
                if w[0][k] == 1:
                    break
                w_assignment += 1
            if not w_assignment in buffer:
                buffer[w_assignment] = 0
            buffer[w_assignment] = buffer[w_assignment] + 1
        ## output
        output_f.write(str(i)+'\t'+str(TERM_PER_DOC)+'\t')
        for word_id, word_count in buffer.items():
            output_f.write(str(word_id)+':'+str(word_count)+' ')
        output_f.write('\n')
        z_f.write(str(i)+'\t'+str(TERM_PER_DOC)+'\t')
        for z_id, z_count in z_buffer.items():
            z_f.write(str(z_id)+':'+str(z_count)+' ')
        z_f.write('\n')
        theta_f.write(str(i)+'\t')
        for k in range(TOPIC_N):
            theta_f.write(str(k)+':'+str(theta[0][k])+' ')
        theta_f.write('\n')
    z_f.close()
    theta_f.close()
    output_f.close()

    ## output phi
    output_f = open(str(path) + FILE_NAME+'.phi','w')
    for i in range(TOPIC_N):
        output_f.write(str(i)+'\t')
        for j in range(VOCABULARY_SIZE):
            output_f.write(str(j)+':'+str(phi[i][0][j])+' ')
        output_f.write('\n')
    output_f.close()

    ## output hyper-parameters
    output_f = open(str(path) + FILE_NAME+'.hyper','w')
    output_f.write('TOPIC_N:'+str(TOPIC_N)+'\n')
    output_f.write('VOCABULARY_SIZE:'+str(VOCABULARY_SIZE)+'\n')
    output_f.write('DOC_NUM:'+str(DOC_NUM)+'\n')
    output_f.write('TERM_PER_DOC:'+str(TERM_PER_DOC)+'\n')
    output_f.write('alpha:'+str(alpha[0])+'\n')
    output_f.write('beta:'+str(beta[0])+'\n')
    output_f.close()
    #observed word matrix of standard unit vector is created.
    #wc = pd.read_csv('testfi.data.doc', header=None, sep='\t')
    #hy = pd.read_csv('testfi.data.hyper', header=None, sep='\t')
    wc = pd.read_csv(str(path) + FILE_NAME + '.doc', header=None, sep='\t')
    hy = pd.read_csv(str(path) + FILE_NAME + '.hyper', header=None, sep='\t')
    nvoc = int(hy[0][1].split(':')[1])#total vocaburary in the corpus
    ndoc = int(hy[0][2].split(':')[1])# # of documets
    nter = int(hy[0][3].split(':')[1])# # of words in a document <-- assume same
    matrix = numpy.zeros((ndoc*nter, nvoc))
    count = 0
    count_total = 0
    a = wc[0]
    for k in range(len(a)):
        #print('documet number:'+ k)
        b = []
        b = wc[2][k]
        b = b.strip(' \t\n\r')
        bb = re.split(':|\s+', b)
        for i in range(0,len(bb),2):
            #word index
            index = int(bb[i])
            #count
            count = int(bb[i+1])
            for j in range(count): #take j is from 0 to count-1
                matrix[count_total+j][index] = 1
                #matrix[k*]
                #print (j)t (k)
            count_total += int(bb[i+1]) 
    scipy.io.savemat('observed_word_matrix_DN' + str(dn) + 'NT' + str(tm) + "VOC" + str(VOCABULARY_SIZE) + "ortho_minus" + str(ortho_count) + "index_expe" + str(expe_count) +'.mat', mdict={'matrix': matrix})

    #ORIGINAL PHI(WORD-TOPIC-DISTRIBUTION) MATRIX
    pp = pd.read_csv(str(path) + FILE_NAME + '.phi', header=None, sep='\t')
    phi_matrix = numpy.zeros((nvoc,len(pp[1])))
    for j in range(len(pp[1])):
        ab = pp[1][j]
        ab = ab.strip(' \t\n\r')
        ab = re.split(':|\s+', ab)
        for i in range(0,len(ab),2):
                #word index
                index = int(ab[i])
                #probability
                pro = float(ab[i+1])
                phi_matrix[index][j] = pro 
    scipy.io.savemat('parameter_original_DN' + str(dn) + 'NT' + str(tm) + "VOC" + str(VOCABULARY_SIZE) + "ortho_minus" + str(ortho_count) + "index_expe" + str(expe_count) +'.mat', mdict={'phi_matrix_original': phi_matrix})
    #brief check for matrix generation
# import os
# os.getcwd()
# import numpy as np
# import scipy
# import pandas as pd
# import re
#  ##Retrieved from http://www.hongliangjie.com/2010/09/30/generate-synthetic-data-for-lda/
# #X. Yi, L. Hong,
# ## This program is to generated synthetic data for LDA
# import math
# import random
# import numpy
# import numpy.random
# import sys
# from scipy import io
# import os
# ## define some constant
# lp = [1000,100,10]
# ortho = [0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001]
# TOPIC_N = 2
# VOCABULARY_SIZE = 10
# #FILE_NAME = sys.argv[1]
# #FILE_NAME = "test"
# #for dn in lp:
# #for tm in lp:
# dn = 1000
# tm = 100 
# #TOPIC_N = 2
# #VOCABULARY_SIZE = 10
# DOC_NUM = dn
# TERM_PER_DOC = tm
# ortho_count = 0
# #beta = [0.01 for i in range(VOCABULARY_SIZE)]
# pri = []
# prior = 0.01
# for i in range(VOCABULARY_SIZE):
#     prior = prior + 0.01
#     pri.append(prior)    
# beta = pri
# #alpha = [0.9 for i in range(TOPIC_N)]
# #alpha = [0.5 for i in range(TOPIC_N)]
# alpha = [0.2, 0.8]
# FILE_NAME = sys.argv[1]
# #for p in ortho:
# p = 0.01
# #ortho_count = ortho_count + 1
# ortho_count = 2
# ## generate multinomial distribution over words for each topic
# a = 0
# zero = p
# zero_out = zero*10000000
# os.system ("mkdir " + "K2ND" + str(dn) + "TE" + str(tm) + "VOC" + str(VOCABULARY_SIZE) + "ortho_minus" + str(ortho_count))
# path = "K2ND" + str(dn) + "TE" + str(tm) + "VOC" + str(VOCABULARY_SIZE) + "ortho_minus" + str(ortho_count)
# while (a == 0):
#     phi = []
#     for i in range(TOPIC_N):
#         topic = numpy.random.mtrand.dirichlet(beta, size = 1)
#         phi.append(topic)
#     p1 = np.squeeze(np.asarray(phi[0]))
#     p2 = np.squeeze(np.asarray(phi[1]))
#     #a = 1
#     dp = np.dot(p1,p2)
#     if dp < zero:
#         a = 1
# ## generate words for each document
# output_f = open(FILE_NAME+'.doc','w')
# z_f = open(FILE_NAME+'.z','w')
# theta_f = open(FILE_NAME+'.theta','w')
# for i in range(DOC_NUM):
#     buffer = {}
#     z_buffer = {} ## keep track the true z
#     ## first sample theta
#     theta = numpy.random.mtrand.dirichlet(alpha,size = 1)
#     for j in range(TERM_PER_DOC):
#         ## first sample z
#         z = numpy.random.multinomial(1,theta[0],size = 1)
#         z_assignment = 0
#         for k in range(TOPIC_N):
#             if z[0][k] == 1:
#                 break
#             z_assignment += 1
#         if not z_assignment in z_buffer:
#             z_buffer[z_assignment] = 0
#         z_buffer[z_assignment] = z_buffer[z_assignment] + 1
#         ## sample a word from topic z
#         w = numpy.random.multinomial(1,phi[z_assignment][0],size = 1)
#         w_assignment = 0
#         for k in range(VOCABULARY_SIZE):
#             if w[0][k] == 1:
#                 break
#             w_assignment += 1
#         if not w_assignment in buffer:
#             buffer[w_assignment] = 0
#         buffer[w_assignment] = buffer[w_assignment] + 1
#     ## output
#     output_f.write(str(i)+'\t'+str(TERM_PER_DOC)+'\t')
#     for word_id, word_count in buffer.items():
#         output_f.write(str(word_id)+':'+str(word_count)+' ')
#     output_f.write('\n')
#     z_f.write(str(i)+'\t'+str(TERM_PER_DOC)+'\t')
#     for z_id, z_count in z_buffer.items():
#         z_f.write(str(z_id)+':'+str(z_count)+' ')
#     z_f.write('\n')
#     theta_f.write(str(i)+'\t')
#     for k in range(TOPIC_N):
#         theta_f.write(str(k)+':'+str(theta[0][k])+' ')
#     theta_f.write('\n')
# z_f.close()
# theta_f.close()
# output_f.close()
 
# ## output phi
# output_f = open(FILE_NAME+'.phi','w')
# for i in range(TOPIC_N):
#     output_f.write(str(i)+'\t')
#     for j in range(VOCABULARY_SIZE):
#         output_f.write(str(j)+':'+str(phi[i][0][j])+' ')
#     output_f.write('\n')
# output_f.close()
 
# ## output hyper-parameters
# output_f = open(FILE_NAME+'.hyper','w')
# output_f.write('TOPIC_N:'+str(TOPIC_N)+'\n')
# output_f.write('VOCABULARY_SIZE:'+str(VOCABULARY_SIZE)+'\n')
# output_f.write('DOC_NUM:'+str(DOC_NUM)+'\n')
# output_f.write('TERM_PER_DOC:'+str(TERM_PER_DOC)+'\n')
# output_f.write('alpha:'+str(alpha[0])+'\n')
# output_f.write('beta:'+str(beta[0])+'\n')
# output_f.close()
# #observed word matrix of standard unit vector is created.
# #wc = pd.read_csv('testfi.data.doc', header=None, sep='\t')
# #hy = pd.read_csv('testfi.data.hyper', header=None, sep='\t')
# wc = pd.read_csv(FILE_NAME + '.doc', header=None, sep='\t')
# hy = pd.read_csv(FILE_NAME + '.hyper', header=None, sep='\t')
# nvoc = int(hy[0][1].split(':')[1])#total vocaburary in the corpus
# ndoc = int(hy[0][2].split(':')[1])# # of documets
# nter = int(hy[0][3].split(':')[1])# # of words in a document <-- assume same
# matrix = numpy.zeros((ndoc*nter, nvoc))
# count = 0
# count_total = 0
# a = wc[0]
# for k in range(len(a)):
#     #print('documet number:'+ k)
#     b = []
#     b = wc[2][k]
#     b = b.strip(' \t\n\r')
#     bb = re.split(':|\s+', b)
#     for i in range(0,len(bb),2):
#         #word index
#         index = int(bb[i])
#         #count
#         count = int(bb[i+1])
#         for j in range(count): #take j is from 0 to count-1
#             matrix[count_total+j][index] = 1
#             #matrix[k*]
#             #print (j)t (k)
#         count_total += int(bb[i+1]) 
# scipy.io.savemat('observed_word_matrix_DN' + str(dn) + 'NT' + str(tm) + "VOC" + str(VOCABULARY_SIZE) + "ortho_minus" + str(ortho_count) + '.mat', mdict={'matrix': matrix})

# #ORIGINAL PHI(WORD-TOPIC-DISTRIBUTION) MATRIX
# pp = pd.read_csv(FILE_NAME + '.phi', header=None, sep='\t')
# phi_matrix = numpy.zeros((nvoc,len(pp[1])))
# for j in range(len(pp[1])):
#     ab = pp[1][j]
#     ab = ab.strip(' \t\n\r')
#     ab = re.split(':|\s+', ab)
#     for i in range(0,len(ab),2):
#             #word index
#             index = int(ab[i])
#             #probability
#             pro = float(ab[i+1])
#             phi_matrix[index][j] = pro 
#scipy.io.savemat('parameter_original_DN' + str(dn) + 'NT' + str(tm) + "VOC" + str(VOCABULARY_SIZE) + "ortho_minus" + str(ortho_count) + '.mat', mdict={'phi_matrix_original': phi_matrix})
#brief check for matrix generation
