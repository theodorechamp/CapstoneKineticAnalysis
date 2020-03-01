from execODE import *
import json
import os

def parsejson(fn):
    with open(fn, "r") as f_obj:
        data = json.load(f_obj)
        for i in data['iterations']:
            #each i will be a different set of parameters
            params = i['params']
            execODE(params)



def main():
    with open("../global.val", 'r') as f_obj:
        globalvals = json.load(f_obj)
        fn_params = globalvals['oxiParamFN']
        parsejson(fn_params)

if __name__ == "__main__":
    main()
