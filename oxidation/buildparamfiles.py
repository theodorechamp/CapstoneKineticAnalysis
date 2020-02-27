import json, codecs
import math as m
import os
from numpy.random import randint
class jsonComm(object):
    """docstring for jsonComm."""
    def __init__(self):
        super(jsonComm, self).__init__()
        try:
            open("params.txt",'x')
            self.f = open("params.txt", 'r+')
            data = {"lastUID":0, "iterations":[]}
            json.dump(data, self.f, indent=2)
        except FileExistsError:
            self.f = open("params.txt", 'r')
            data = json.load(self.f)
        self.fs = data
        print(data)
        self.f = open("params.txt", 'w')


    def jsonoutput(self,params):
        print(self.fs)
        UID = self.fs["lastUID"]+1
        self.fs["lastUID"] = UID
        data = {"UID":UID,"params":params}
        self.fs["iterations"].insert(0,data)
        print(self.fs)
        json.dump(self.fs,self.f,indent=2)

def randparams(params):
    for key in params:
        test = 0
        params[key] = randint(20)
        print(params)
    return ([params, testCaseID])

def main():
    params = {
        'k1':10**-1.55,
        'gamma1':2.5,
        'n1':2.0,
        'k3':10**-1,
        'n3':1.5,
        'Vco':11,
        'Ftot':5,
        'conversion':.5,
        'tspan':10,
        'tstep':.1,
        'mH': 70,
        'MWH':120,
        'delta':5,
        'Vr': 10,
        'tau':1,
        'n':1,
        'x':1 }
    j = jsonComm()
    j.jsonoutput(params)

    # params = randparams(params)
    # jsonoutput(params,testCaseID,"testparams.json")

if __name__ == "__main__":
    main()
