
# coding: utf-8

# In[ ]:


import json
import math as m
import os

def jsonoutput(params, fn):
    savelocation = os.getcwd() + "/paramfiles/" + fn
    with open(savelocation,"w") as f_obj:
        json.dump(params, f_obj)

def main():
    params = {
        'A':1998.196, # vol^2/(mol^2*s) 
        'Ea':254, # kJ/mol
        'R':0.008314472, # kJ/(K*mol)
        'T':1573, # K
        'conversion':.5, # vary conversion and look at how the time needed changes
            }
    jsonoutput(params,"testparams.json")




if __name__ == "__main__":
    main()

