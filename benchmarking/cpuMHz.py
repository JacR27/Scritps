import sys

def main():
    totalMHz = 0
    for line in sys.stdin.readlines():
        MHz = line.strip().split()[-1]
        totalMHz = totalMHz + float(MHz)
        
    sys.stdout.write(str(totalMHz)+"\n")

main()
