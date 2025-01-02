import numpy as np


#Função de elim de linhas LI fornecido pelo professor Cristiano Arbex
def makeMatrixFullRank(A):
    ''' Esta função recebe uma matriz, que pode ser em numpy,
        e retorna dois argumentos:
          - A matriz com linhas eliminadas
          - Uma lista indicando quais linhas foram eliminadas.
    '''
    if np.linalg.matrix_rank(A) == A.shape[0]: return A, []
    row = 1
    rowsEliminated = []
    counter = 0
    while 1:
        counter += 1
        B = A[0:(row+1), :]
        C = np.linalg.qr(B.T)[1]
        C[np.isclose(C, 0)] = 0
        if not np.any(C[row, :]):
            rowsEliminated.append(counter)
            A = np.delete(A, (row), axis=0)
        else:
            row += 1
        # end if
        if row >= A.shape[0]: break
    # end for
    return A, rowsEliminated
# end makeMatrixFullRank

#######################################################


def Abc(M):
  obj = M[0][-1]
  c = M[0][:-1]
  b = []
  A = []
  for r in M[0:]:
    b.append(r[-1])
    A.append(r[:-1])
  A = np.array(A)
  b = np.array(b)

  return A,b,c,obj


def Pivot(M, l1, l2, mul):
  print(M[l1])
  aux = M[l1] + M[l2]*mul
  M[l1] = aux
  print(M[l2]*mul)
  print(M[l1])
  return M

def MatrizInit(filename):
    with open(filename,'r') as f:
        #lendo arquivo
        nvars = int(f.readline())
        nres = int(f.readline())
        VarsRes = list(map(int, f.readline().split()))
        coefs = list(map(int, f.readline().split()))
        #restricoes
        res = []
        for l in f:
            res.append(l.split())
        
        #iterando pelas restricoes e extraindo info delas
        b = []
        nNewvars = 0
        isNewvar = [0]*nvars #check para ver se aquela restricoes recebe variavel de folga
        count1 = 0
        for r in res:
            b.append(r[-1]) #construindo b
            if(r[nvars] != '=='): #se a restriçao ja não é uma igualdade
                nNewvars += 1
                isNewvar[count1] = 1
            count1 += 1
        b = list(map(int, b)) #b é matriz de int

        c = coefs + [0]*(nNewvars + 1) #c ampliado no numero de novas vars

        #construindo A
        A = []
        count2 = 0
        for r in res:
            if(isNewvar[count2] == 0):
                newres = r[:nvars] + [0]*(nNewvars-1)
            else:
                newres = r[:nvars] + [0]*count2 + [1] + [0]*(nNewvars-(count2+1))
            newres = list(map(int, newres))
            if(r[nvars] == '>='):
                newres = [-x for x in newres[:nvars]]+newres[nvars:]
                b[count2] *= -1
            A.append(newres)
            count2 += 1


    M = [c]
    aux = 0
    for l in A:
        M.append(l+[b[aux]])
        aux += 1

 
    M = np.array(M,dtype=float)


    return(M)

    
M = MatrizInit('test.txt')
print(M)
(A,b,c,obj) = Abc(M)
print(A,b,c,obj)
    

