import argparse
parser = argparse.ArgumentParser()
#parser.add_argument("echo")
parser.add_argument("--verbose", help="increase output verbosity" , action="store_true")
args = parser.parse_args()
#print(args.echo)
print(args.verbose,type(args.verbose))
