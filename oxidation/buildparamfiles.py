import json, codecs
import numpy
import sys
import os
import random
import itertools
class jsonComm(object):
    """docstring for jsonComm."""
    def __init__(self):
        super(jsonComm, self).__init__()
        self.fn = "params.json"
        try:
            open(self.fn,'x')
            self.f = open(self.fn, 'r+')
            data = {"lastUID":0, "iterations":[]}
            json.dump(data, self.f, indent=2)
        except FileExistsError:
            self.f = open(self.fn, 'r')
            data = json.load(self.f)
        self.fs = data
        self.f = open(self.fn, 'w')

    def jsonoutput(self,params):
        UID = self.fs["lastUID"]+1
        self.fs["lastUID"] = UID
        data = {"UID":UID,"params":params}
        self.fs["iterations"].insert(0,data)
        json.dump(self.fs,self.f,indent=2)

    def jsonoutputgrid(self,params):
        ## Build grid layout
        ranges = []
        i=0
        for item in params:
            key = params[item]
            ranges[i] = numpy.linspace(key[1],key[2],2) #linspace num (3rd arg) can be increased.
            i=i+1
        list(itertools.product(*ranges))
        #https://stackoverflow.com/questions/798854/all-combinations-of-a-list-of-lists
        ## Write grid layout to file  
        UID = self.fs["lastUID"]+1
        self.fs["lastUID"] = UID
        data = {"UID":UID,"params":params}
        self.fs["iterations"].insert(0,data)
    def jsondump(self):
        json.dump(self.fs,self.f,indent=2)

def randparams(params):
    for item in params:
        key = params[item]
        key[0] = random.uniform(key[1],key[2])
    return (params)

def randparamsgrid(params, j):
    for item in params:
        key = params[item]
        gridVals = numpy.linspace(key[1],key[2],2)
        for val in gridVals:
            key[0] = val
            j.jsonoutputgrid(params)
    j.jsondump()


    return (params)

def main():
    params = {
        'Vco':[11,1,10],
        'Ftot':[5,1,10],
        'conversion':[.5,0,1],
        'tspan':[10,1,100],
        'mH': [70,60,100],
        'MWH':[120,100,200],
        'delta':[5,1,10],
        'Vr': [10,1,100],
        'tau':[10,1,10],
        'n':[1,1,10],
        'x':[1,1,10]
    }

    j = jsonComm()
    try:
        input = sys.argv[1]
        if(input=='-i' or input=='--init'):
            params = randparamsgrid(params,j)
    except IndexError:
        params = randparams(params)
        j.jsonoutput(params)


if __name__ == "__main__":
    main()
