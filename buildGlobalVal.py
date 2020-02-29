import json

vals = {
    'oxiParamFN':"params.json"
}

def main():
    fn = "global.val"
    with open(fn, "w") as f_obj:
        json.dump(vals,f_obj)



if __name__ == "__main__":
    main()
