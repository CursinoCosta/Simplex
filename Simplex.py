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


def Abc(M): #recebe uma matriz M e retorna suas parte matiz A, vetores b,c e o valor objetivo obj
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


def Pivot(M, lalvo, l2, mul, id):
  aux = M[lalvo] + M[l2]*mul
  M[lalvo] = aux.round(decimals=2)

  auxi = id[lalvo] + id[l2]*mul
  id[lalvo] = auxi.round(decimals=2)
  return M

def CanonCheck(M): #checa se uma matriz esta na forma canonica
  if(np.all(M[0][:-1] != 0)):
    return False
  else:
    pos = [False]*(M.shape[0]-1)
    count = 0
    for col in M.T:
      if(np.all(np.isin(col,[0,1])) and col.sum()==1 and M[0][count]==0):
        for i in range(1,M.shape[0]):
          if(M[i][count]==1):
            pos[i-1] = True
            break
            
      count+=1
    if(np.all(pos)): 
      return True
  return False

def MatrizInit(filename): #le a entrada e retorna a FPI correspondente
    with open(filename,'r') as f:
        #lendo arquivo
        nvars = int(f.readline())
        nres = int(f.readline())
        VarsRes = list(map(float, f.readline().split()))
        coefs = list(map(float, f.readline().split()))
        #restricoes
        res = []
        for l in f:
            res.append(l.split())
    
    #eliminando restricoes linearmente dependentes
    AuxA = [] 
    for r in res:
        reduc = r[:nvars] + [r[-1]]
        reduc = list(map(float, reduc))
        AuxA.append(reduc)
    AuxA = np.array(AuxA)

    _ ,  elim = makeMatrixFullRank(AuxA)

    nres -= np.size(elim)
    res = np.delete(res,elim,axis=0).tolist()

    #iterando pelas restricoes e extraindo info delas
    b = []
    nNewvars = 0
    isNewvar = [0]*nres #check para ver se aquela restricoes recebe variavel de folga
    count1 = 0
    for r in res:
        b.append(r[-1]) #construindo b
        if(r[nvars] != '=='): #se a restriçao ja não é uma igualdade
            nNewvars += 1
            isNewvar[count1] = 1
            if(r[nvars] == '>='):
               VarsRes.append(-1)
            else:
               VarsRes.append(1)
        count1 += 1
    b = list(map(float, b)) #b é matriz de float

    c = coefs + [0]*(nNewvars + 1) #c ampliado no numero de novas vars

    #construindo A
    A = []
    count2 = 0
    for r in res:
        if(isNewvar[count2] == 0):
            newres = r[:nvars] + [0]*(nNewvars-1)
        else:
            newres = r[:nvars] + [0]*count2 + [1] + [0]*(nNewvars-(count2+1))

        newres = list(map(float, newres))
        A.append(newres)
        count2 += 1


    M = [c]
    aux = 0
    for l in A:
        M.append(l+[b[aux]])
        aux += 1


 
    M = np.array(M,dtype=float)

#tornando todas as variaveis restritas >=0
    count = 0
    for i in VarsRes:
        if(i != 1):
            if(i == -1): 
                for l in M:
                    l[count] *= -1
                VarsRes[count] = 1
            else:
                newM = []
                for l in M:
                    newl = list(l)
                    newl.insert(count+1,l[count]*-1)
                    newM.append(newl)
                M = np.array(newM)
                nvars +=1
                VarsRes[count] = 1
                VarsRes.insert(count+1, 1)
        count += 1


    return(M)

def protdiv(objs, col):
  out = []
  for i in range(len(col)):
    if(col[i] <= 0):
      out.append(float('inf'))
    else:
      out.append(round(objs[i]/col[i],2))
  print(out) #bom print pro passo a passo
  return out

def findIdent(M): #funcao que acha o que falta da identeidade e para a canonizacao da matriz
    Ipos = [False]*(M.shape[0]-1)
    Cpos = [False]*(M.shape[0]-1)
    count = 0
    incompletos = []
    for col in M.T:
        if(np.all(np.isin(col[1:],[0,1])) and col[1:].sum()==1):
            for i in range(1,M.shape[0]):
                if(M[i][count]==1):
                    Ipos[i-1] = True
                    if(M[0][count]==0):
                        Cpos[i-1] = True
                    else:
                        incompletos.append((count,i))
                break
        count += 1
    return Ipos,incompletos

def Simplex(M,id,exec): #simplex
    print('iniciei')
    for i in range(1,M.shape[0]):
       if(M[i][-1] < 0):
          M[i] *= -1
          id[i] *= -1
    M[0] *= -1
    cminus = M[0]
    print('************','\n',id,'\n',M,'\n---------\n')
    itr = 1
    otimo = False
    if(not CanonCheck(M)):
       print('not canon')
       print(M)      
       (Ipos,incompletos) = findIdent(M)

       if(np.any(Ipos)):
            for col,row in incompletos:
                M = Pivot(M,0,row,-(M[0][col]/M[row][col]),id)

       if(np.any(Ipos.count(False))):
            auxc = [0]*(M.shape[1]-1) + [1]*Ipos.count(False) + [0]
            auxM = M.copy()
            I = np.insert(np.identity(auxM.shape[0]-1),0,-1,axis=0)
            for i in range(len(Ipos)):
                if(Ipos[i]): continue
                else:
                  auxM = np.insert(auxM,auxM.shape[1]-1,I.T[i],axis=1)
            auxM[0] = auxc
            (auxIpos,auxincompletos) = findIdent(auxM)
            for col,row in auxincompletos:
                auxM = Pivot(auxM,0,row,-(auxM[0][col]/auxM[row][col]),id)
            auxM[0] = -auxM[0]
            print(id)
            exec +=1
            id = Simplex(auxM,id,exec)
            
            print('b',type(id))
            print(exec)
            return('parada')
            print('a',auxM)
            print(id)
            if(auxM[0][-1] != 0):
                return('Inviavel')
            else:
                return('canonize')
    while(not otimo):
        M[0][np.isclose(M[0], 0)] = 0
        print('iteração: ',itr)
        print('c: ',cminus)

        pivo = np.argmin(cminus[:-1])
        print('pivo: ',pivo)

        lchoice = (np.argmin(protdiv(M.T[-1][1:],M.T[pivo][1:]))+1) 
        print('lchoice: ',lchoice)  

        temp = (id[lchoice]/M[lchoice][pivo]).round(decimals=2)
        M[lchoice] = (M[lchoice]/M[lchoice][pivo]).round(decimals=2)
        id[lchoice] = temp

        for i in range(M.shape[0]):
            if(i == lchoice or M[i][pivo] == 0): continue
            M = Pivot(M,i,lchoice,-(round(M[i][pivo]/M[lchoice][pivo],2)),id)

        cminus = (M[0])
        if(np.all(cminus[:-1] >= 0)):
            otimo = True
        itr+=1
        print(M,'\n---------\n')
        print(id,'\n---------\n')
    out = id
    print('a',type(out))
    return out  


    
M = MatrizInit('test.txt')
id = np.insert(np.identity(M.shape[0]-1),0,0,axis=0)
#print(M)
Simplex(M,id,0)

# (A,b,c,obj) = Abc(M)
# print(A,b,c,obj)
    

