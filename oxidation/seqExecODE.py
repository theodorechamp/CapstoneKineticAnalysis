from execODE import *


def parsejson(fn):
    with open(fn, "r") as f_obj:
        data = json.load(f_obj)
def main():
    print("Hello")
    parsejson("testparams.json")

if __name__ == "__main__":
    main()
