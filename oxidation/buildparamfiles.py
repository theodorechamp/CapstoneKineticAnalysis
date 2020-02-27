import json, codecs
import math as m
import os
from numpy.random import randint
class jsonComm(object):
    """docstring for jsonComm."""
    def __init__(self):
        super(jsonComm, self).__init__()
        if(open("params.json","x"))
            self.f_obj = open("params.json","+")
            data = {"lastUID":0, "iterations":[]}
            json.dump(data, f_obj)
        else
            self.f_obj = open("params.json", "+")

<<<<<<< HEAD
def jsonoutput(params, fn):
    savelocation = os.getcwd() + "/paramfiles/" + fn
    with open(savelocation,"a") as f_obj:
        json.dump(params, f_obj)
=======
>>>>>>> b808caf61111adde11e3eede02cfa5e5a072d535

    def createfile():
        test = 1

    def jsonoutput(params, testCaseID, fn):
        #If file doesn't exist, initialize with [] for json --> add later
        # if not (f_obj[testCaseID]):
        #     f_obj[testCaseID] = []
        # f_obj[testCaseID].append(params)
        data = {}
        data[testCaseID] = []
        data[testCaseID].append(params)
        print(data)
        json.dump(data, f_obj)

    def getUID():
        f = json.load(f_obj)
        UID = f['lastUID']+1
        return(UID)

def randparams(params):
    for key in params:
        test = 0
        params[key] = randint(20)
        print(params)
    return ([params, testCaseID])

def main():
    params = {
<<<<<<< HEAD
        'k1':10**-1.55,
        'gamma1':2,
=======
        'k1':10**-1.55,1,10],
        'gamma1':2.5,
>>>>>>> b808caf61111adde11e3eede02cfa5e5a072d535
        'n1':2.0,
        'k3':10**-1,
        'n3':1.5,
        'Vco':10,
        'Ftot':5,
        'conversion':.5,
        'tspan':10,
        'tstep':.1,
        'mH': 70,
        'MWH':120,
        'alpha0':5,
        'Vr': 10,
        'tau':1,
        'n':1,
        'x':1 }
    j = jsonComm

    params = randparams(params)
    jsonoutput(params,testCaseID,"testparams.json")

if __name__ == "__main__":
    main()
