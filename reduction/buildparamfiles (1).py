
# coding: utf-8

# In[ ]:


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
        'R':0.008314472, # kJ/(K*mol)
        'T':1573 # K --> could be varied from 1553 K to 1693 K later on
        }
    
    params = {
        'A':[812.41,812.41,1808], # [value,min of range,max of range]; pre-exponential factor vol^2/(mol^2*s)
        'Ea':[256,256,262], # activation energy kJ/mol
        'conversion':[0,0,1] # vary conversion and look at how the time needed changes
    }
    
    j = jsonComm()
    params = randparams(params)
    j.jsonoutput(params)

    #
    # jsonoutput(params,testCaseID,"testparams.json")

if __name__ == "__main__":
    main()

