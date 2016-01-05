import sys
lines = []
totalsize = 0
for line in sys.stdin.readlines():
    info = line.strip().split()
    lines.append(info)
    totalsize = totalsize + int(info[4])
print (totalsize)
    
