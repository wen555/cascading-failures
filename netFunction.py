import networkx as nx

def edge2file(edge, fileName):
    file=open(fileName,'a')
    file.seek(0)
    file.truncate()
    edgeList=list(edge)
    for l in edgeList:
        edgeStr=str(1+l[0])+' '+str(1+l[1])+' '+'1'+'\n'
        file.write(edgeStr)
    file.close()
def degreefile(edge, fileName):
    file=open(fileName,'a')
    file.seek(0)
    file.truncate()
    edgeList=list(edge)
    for l in edgeList:
        edgeStr=str(l[1])+'\n'
        file.write(edgeStr)
    file.close()
def matrix_file(matrixlist, fileName):
    file=open(fileName,'a')
    file.seek(0)
    file.truncate()
    # matrixStr=str(matrixlist)
    # print(matrixStr)
    for l in matrixlist:
        for i in l:
            file.write(str(i)+' ')
        file.write('\n')
    file.close()
# score=nx.degree_centrality(G)
# score = sorted(score.items(), key=lambda item: item[1], reverse=True)
def Net_algorthm(file_place,fun,x):
    # fun=str(fun)
    score=eval(fun)(x)
    # print(fun, score)
    output = []
    for key,value in score.items():
        output.append(value)
    # print(output)
    fout = open(file_place+fun+'.txt', 'w')
    fout.seek(0)
    fout.truncate()
    for target in output:
        fout.write(str(target) + "\n")
def edge_betweenness_centrality2file(ebc, fileName):
    file=open(fileName,'a')
    file.seek(0)
    file.truncate()
    for key, values in ebc.items():
        edgeList = list(key)
        edgeStr=str(1+edgeList[0])+' '+str(1+edgeList[1])+' '+str(values)+'\n'
        file.write(edgeStr)
