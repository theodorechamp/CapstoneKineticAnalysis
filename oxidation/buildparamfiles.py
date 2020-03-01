import json, codecs
import math as m
import os
import random
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

def randparams(params):
    for item in params:
        key = params[item]
        key[0] = random.uniform(key[1],key[2])
    return (params)

def main():
    constants = {
        'k1':10**-1.55,
        'gamma1':2.5,
        'n1':2.0,
        'k3':10**-1,
        'n3':1.5,
        'tstep':.1
    }
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
    params = randparams(params)
    j.jsonoutput(params)

    #
    # jsonoutput(params,testCaseID,"testparams.json")

if __name__ == "__main__":
    main()
