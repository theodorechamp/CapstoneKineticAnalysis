import json, codecs
import math as m
import os
from numpy.random import randint

def jsonoutput(params, testCaseID, fn):
    #If file doesn't exist, initialize with [] for json --> add later
    savelocation = os.getcwd() + "/paramfiles/" + fn
    with open(savelocation,"w") as f_obj:
        # if not (f_obj[testCaseID]):
        #     f_obj[testCaseID] = []
        # f_obj[testCaseID].append(params)
        data = {}
        data[testCaseID] = []
        data[testCaseID].append(params)
        print(data)
        json.dump(data, f_obj)

def randparams(params):
    testCaseID = "test"
    for key in params:
        test = 0
        params[key] = randint(20)
        print(params)
    return ([params, testCaseID])

def main():
    testCaseID = "YA"
    params = {
        'k1':10**-1.55,
        'gamma1':2.5,
        'n1':2.0,
        'k3':10**-1,
        'n3':1.5,
        'Vco':10,
        'Ftot':5,
        'conversion':.5,
        'tspan':10,
        'tstep':.1 }

    [params, testCaseID] = randparams(params)
    jsonoutput(params,testCaseID,"testparams.json")




if __name__ == "__main__":
    main()
