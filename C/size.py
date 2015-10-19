import sys
total = 0
for line in sys.stdin.readlines():
    fileInfo = line.strip().split();
    size = fileInfo[4]
    total += int(size)
print(total)
