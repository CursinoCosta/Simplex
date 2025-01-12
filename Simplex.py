import numpy as np
import sys


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
    while row >= A.shape[0]:
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
        #if : break
    # end for
    return A, rowsEliminated
# end makeMatrixFullRank

#######################################################


def Pivot(M, lalvo, l2, mul, id, decimals):
  aux = M[lalvo] + M[l2]*mul
  M[lalvo] = np.round(aux,decimals)

  auxi = id[lalvo] + id[l2]*mul
  id[lalvo] = np.round(auxi,decimals)
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
        ogvars = nvars
        nres = int(f.readline())
        VarsRes = list(map(float, f.readline().split()))
        coefs = f.readline().split()
        minmax = coefs[0]
        coefs = list(map(float, coefs[1:]))
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
    isNewvar = [False]*nres #check para ver se aquela restricoes recebe variavel de folga
    count = 0
    for r in res:
        b.append(r[-1]) #construindo b
        if(r[nvars] != '=='): #se a restriçao ja não é uma igualdade
            nNewvars += 1
            isNewvar[count] = True
            if(r[nvars] == '>='):
               VarsRes.append(-1)
            else:
               VarsRes.append(1)
        count += 1
    b = list(map(float, b)) #b é matriz de float

    c = coefs + [0]*(nNewvars + 1) #c ampliado no numero de novas vars

    #construindo A
    A = []
    count = 0
    for r in res:
        if(not isNewvar[count]):
            newres = r[:nvars] + [0]*(nNewvars)
        else:
            newres = r[:nvars] + [0]*count + [1] + [0]*(nNewvars-(count+1))

        newres = list(map(float, newres))
        A.append(newres)
        count += 1
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

    if(minmax == 'min'):
        M[0] *= -1
    return(M, minmax, ogvars)


def Mprint(id,M):
    for i in range(M.shape[0]):
        print(id[i],' | ',M[i])
        

def protdiv(objs, col,decimals = 2):
  out = []
  for i in range(len(col)):
    if(col[i] <= 0):
      out.append(float('inf'))
    else:
      out.append(round(objs[i]/col[i],decimals))
  print('Escolhendo linha:')    
  print(out) #bom print pro passo a passo
  return out

def findIdent(M): #funcao que acha o que falta da identeidade e para a canonizacao da matriz
    Ipos = [False]*(M.shape[0]-1)
    Cpos = [-1]*(M.shape[1]-1)
    count = 0
    incompletos = []
    for col in M.T:
        if(np.all(np.isin(col[1:],[0,1])) and col[1:].sum()==1):
            for i in range(1,M.shape[0]):
                if(M[i][count]==1):
                    Ipos[i-1] = True
                    if(M[0][count]==0):
                        Cpos[count] = i
                    else:
                        incompletos.append((count,i))
                    break
        count += 1
    return Ipos,incompletos, Cpos

def Caminhador(M, id,policy, decimals = 2):
    cminus = M[0]
    print('Caminhando no espaço de busca:')
    itr = 1
    otimo = False
    if(np.all(cminus[:-1] >= 0)):
        otimo = True
    while(not otimo):
        M[0][np.isclose(M[0], 0)] = 0
        print('iteração: ',itr)
        print('c: ',cminus)
        Mprint(id,M)

        pivo = np.argmax((cminus[:-1])*-1)
        if(policy == 'smallest'):
            pivo = np.argmin(cminus[:-1])
        if(policy == 'bland'):
            counter = 0
            for item in cminus[:-1]:
                if(item < 0):
                    pivo = counter
                    break
                counter += 1
        print('pivo: ',pivo)

        divs = protdiv(M.T[-1][1:],M.T[pivo][1:],decimals)
        if(divs.count(float('inf')) == len(divs)):
            return('ilimitado')
        lchoice = (np.argmin(divs)+1) 
        print('linha escolhida: ',lchoice)  

        temp = np.round((id[lchoice]/M[lchoice][pivo]),decimals)
        M[lchoice] = np.round((M[lchoice]/M[lchoice][pivo]),decimals)
        id[lchoice] = temp

        for i in range(M.shape[0]):
            if(i == lchoice or M[i][pivo] == 0): continue
            M = Pivot(M,i,lchoice,-(round(M[i][pivo]/M[lchoice][pivo])),id,decimals)

        cminus = (M[0])
        if(np.all(cminus[:-1] >= 0)):
            otimo = True
        itr+=1
        Mprint(id,M)
        print('xxxxxxxxxxxxxxxx\n')
    print('----------------')
    return M, id



def Simplex(M,id,minmax = 'max', policy  = 'largest', decimals = 2): #simplex
    backup = M.copy()
    for i in range(1,M.shape[0]):
       if(M[i][-1] < 0):
          M[i] *= -1
          id[i] *= -1
    M[0] *= -1
    print('Tableau incial:')
    Mprint(id,M)
    print('\n')
    itr = 1
    otimo = False
    if(not CanonCheck(M)):
       print('----------------')    
       print('Processo de canonização')   
       (Ipos,incompletos,_) = findIdent(M)

       if(np.any(Ipos)):
            for col,row in incompletos:
                M = Pivot(M,0,row,-(M[0][col]/M[row][col]),id,decimals)
            print('Pivoteamento para c')
            Mprint(id,M)

       if(np.any(Ipos.count(False))):
            auxc = [0]*(M.shape[1]-1) + [1]*Ipos.count(False) + [0]
            auxM = M.copy()
            I = np.insert(np.identity(auxM.shape[0]-1),0,-1,axis=0)
            n_aux_r = 0
            for i in range(len(Ipos)):
                if(Ipos[i]): continue
                else:
                    auxM = np.insert(auxM,auxM.shape[1]-1,I.T[i],axis=1)
                    n_aux_r += 1
            auxM[0] = auxc
            (auxIpos,auxincompletos,_) = findIdent(auxM)
            for col,row in auxincompletos:
                auxM = Pivot(auxM,0,row,-(auxM[0][col]/auxM[row][col]),id,decimals)
            auxM[0] = -auxM[0]
            tup = Caminhador(auxM,id,policy,decimals)
            if(type(tup) == str):
                return("ilimitado")
            (auxM, id) = tup
            if(auxM[0][-1] != 0):
                return('inviavel')
            else:
                M = np.delete(auxM, [x for x in range(-(n_aux_r+1),-1)],axis=1)
                M[0] = backup[0]
                (fIpos,f_incompletos,_) = findIdent(M)
                if(np.all(fIpos)):
                    for col,row in f_incompletos:
                        M = Pivot(M,0,row,-(M[0][col]/M[row][col]),id,decimals)
                    M[0] *= -1
                else: return('Erro!')
            print('Canonizada:')
            Mprint(id,M)
    tup = Caminhador(M,id,policy,decimals)
    if(type(tup) == str):
        return('ilimitado')
    (M,id) = tup
    return M, id  



filename = sys.argv[1]
args = sys.argv[2:]
decimals = 2
digits = 1
policy = 'largest'

for param, val in zip(args[::2],args[1::2]):
    if(param == '--decimals'):
        decimals = int(val)
    if(param == '--digits'):
        digits = int(val)
    if(param == '--policy'):
        policy = val
    
np.set_printoptions(formatter={'float': lambda x: f"{x: .{digits}f}"})

(M, minmax, nvars) = MatrizInit(filename)
id = np.insert(np.identity(M.shape[0]-1),0,0,axis=0)

status = Simplex(M,id,minmax,policy,decimals)
if(type(status) == str):
    print('Status: ',status)
    exit(0)
else:
    (M,id) = status
l = M[0][:-1].copy().tolist() #numero de 0 em c
multi = l.count(0) > M.shape[0]

__, _, Cpos = findIdent(M)
x = []
for i in range(M.shape[1]-1):
    x.append(M[Cpos[i]][-1])


otimo = M[0][-1]
if(minmax == 'min'): 
    otimo *= -1



if(multi):
    M[0][np.isclose(M[0], 0)] = 0
    print('iteração: Otimo 2')
    print('c: ',M[0])

    for i in range(M.shape[1]):
        if(not (np.all(np.isin(M.T[i][1:],[0,1])) and M.T[i][1:].sum()==1) and M.T[i][0] == 0):
            pivo = i
            break
    print('pivo: ',pivo)

    divs = protdiv(M.T[-1][1:],M.T[pivo][1:],decimals)
    lchoice = (np.argmin(divs)+1) 
    print('lchoice: ',lchoice)  

    temp = np.round((id[lchoice]/M[lchoice][pivo]),decimals)
    M[lchoice] = np.round((M[lchoice]/M[lchoice][pivo]),decimals)
    id[lchoice] = temp

    for i in range(M.shape[0]):
        if(i == lchoice or M[i][pivo] == 0): continue
        M = Pivot(M,i,lchoice,-(round(M[i][pivo]/M[lchoice][pivo],decimals)),id,decimals)

    __, _, Cpos= findIdent(M)
    x2 = []
    for i in range(M.shape[1]-1):
        x2.append(M[Cpos[i]][-1])
    print(M,'\n---------\n')
    print(id,'\n---------\n')
    print('xxxxxxxx')

print('Status: otimo')
print('Objetivo: ',otimo)
if(multi):
    print('Solucoes:')
    print(x[:nvars])
    print(x2[:nvars])
else:
    print('Solucao:')
    print(x[:nvars])
print('Dual:')
print(id[0])

# (A,b,c,obj) = Abc(M)
# print(A,b,c,obj)
    

