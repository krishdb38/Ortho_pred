import argparse
parser = argparse.ArgumentParser(description="Process some integers")

parser.add_argument("-s",dest="accumulate",action="store_const",const=sum , default=max,help="sum the integers (default max")
parser.add_argument("integer",metavar="N",type = int , nargs="+",help="An integer for the accumulator")
args = parser.parse_args()
print(args.accumulate(args.integer))
print()