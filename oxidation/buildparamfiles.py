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
                'Vco':l[0],
                'Ftot':l[1],
                'conversion':l[2],
                'tspan':l[3],
                'mH': l[4],
                'MWH':l[5],
                'delta':l[6],
                'Vr': l[7],
                'tau':l[8],
                'n':l[9],
                'x':l[10]
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
    try: #checks if gridoutput is desired (called "initialization" here as it provides a basis for successive iteration).
        input = sys.argv[1]
        if(input=='-i' or input=='--init'):
            params = j.jsonoutputgrid(params)
    except IndexError:
        params = j.randparams(params)
        j.jsonoutput(params)

if __name__ == "__main__":
    main()
