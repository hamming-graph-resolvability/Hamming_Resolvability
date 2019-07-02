from scipy.io import loadmat
import pickle

def parseData(filename):
    mat = loadmat(filename)
    data = mat['tests'][0]
    ILP_times = {}
    for d in data:
        d = d[0]
        As = list(ILP_times.keys())
        k = d[0][0][0]
        a = d[1][0][0]
        R = d[2][0]
        R = R.split(',')
        isResolving = d[3][0][0]
        times = d[4][0]
        if a not in As:
            ILP_times[a] = {}     
        ks = list(ILP_times[a].keys())
        if k not in ks:
            ILP_times[a][k] = []
        for time in times:
            ILP_times[a][k].append((R,isResolving,time))
    return ILP_times

filename = str(input('Filename: '))
ILP_times = parseData(filename)
with open('ILP_times_all.dict','wb') as f:
    pickle.dump(ILP_times,f)
