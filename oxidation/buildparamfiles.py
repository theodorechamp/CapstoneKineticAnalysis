import json
import math as m
import os

def jsonoutput(params, fn):
    savelocation = os.getcwd() + "/paramfiles/" + fn
    with open(savelocation,"w") as f_obj:
        json.dump(params, f_obj)


def main():
    params = {
        'k1':10**-1.55,
        'gamma1':2.5,
        'n1':2.0,
        'k3':10**-1,
        'n3':1.5,
        'Vco':10,
        'Ftot':5 }
    jsonoutput(params,"testparams.json")




if __name__ == "__main__":
    main()
