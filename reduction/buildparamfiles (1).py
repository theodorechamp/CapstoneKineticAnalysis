
# coding: utf-8

# In[ ]:


import json, codecs
import numpy
import sys
import os
import random
import itertools

class jsonComm(object):
    """docstring for jsonComm."""
    def __init__(self):
        ## params.json validation and json skeleton construction
        super(jsonComm, self).__init__()
        self.fn = "params.json"
        try: #check if file exists to read in pre-existing data, if not, create new file with json skeleton.
            open(self.fn,'x')
            self.f = open(self.fn, 'r+')
            data = {"lastUID":0, "iterations":[]}
        except FileExistsError:
            self.f = open(self.fn, 'r')
            # data = json.load(self.f)
            try:
                data = json.load(self.f)
            except json.decoder.JSONDecodeError: #if file is empty, default skeleton is used.
                print('here')
                data = {"lastUID":0, "iterations":[]}
        self.fs = data
        self.f = open(self.fn, 'w')

    def jsonoutput(self,params):
        ##prepend the passed params to params.json iterations list and increments lastUID.
        UID = self.fs["lastUID"]+1
        self.fs["lastUID"] = UID
        data = {"UID":UID,"params":params}
        self.fs["iterations"].insert(0,data)
        json.dump(self.fs,self.f,indent=2)

    def jsonoutputgrid(self,params):
        ##Builds grid of all possible variable combinations.
        # Build grid layout
        ranges = []
        i=0
        for item in params:
            key = params[item]
            ranges.append(numpy.linspace(key[1],key[2],2)) #linspace num (3rd arg) can be increased.
            i=i+1
        params = list(itertools.product(*ranges)) #The itertools.product() tool computes the cartesian product of input iterables. It is equivalent to nested for-loops. https://stackoverflow.com/questions/798854/all-combinations-of-a-list-of-lists

        # Write grid layout to file
        i=0
        for l in params:
            data={ #We can get rid of the dictionary and just print list of values to simplify json file, if desired.
                'A':l[0],
                'Ea':l[1],
                'conversion':l[2],
                }
            UID = self.fs["lastUID"]+1
            self.fs["lastUID"] = UID
            dataFormatted = {"UID":UID,"params": data}
            self.fs["iterations"].insert(0,dataFormatted)
            i=i+1
        json.dump(self.fs,self.f,indent=2)

    def randparams(self, params):
        ##randomizes the values of params in the range provided internal to the dicitionary.
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
    try: #checks if gridoutput is desired (called "initialization" here as it provides a basis for successive iteration).
        input = sys.argv[1]
        if(input=='-i' or input=='--init'):
            params = j.jsonoutputgrid(params)
    except IndexError:
        params = j.randparams(params)
        j.jsonoutput(params)

if __name__ == "__main__":
    main()
    

    
    
    

