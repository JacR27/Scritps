import sys

def fillter_and_sort_tags():
    for line in sys.stdin:
        line = line.strip().split('\t')
        print("\t".join(line[0:11] + sorted([i for i in line[11:] if not i.startswith("NM") and not i.startswith('MD')])))

def parse_bam_diff():
    firstline = set()
    secondline = set()
    for line in sys.stdin:
        if line.startswith('H0CPQALXX'):
            readname =line.strip()
        elif line.strip().startswith('<'):
            firstline = filtertags(line)
        elif line.strip().startswith('>'):
            secondline = filtertags(line)
            diff = firstline.symmetric_difference(secondline)
            if diff:
                print(readname)
                print("\t".join(diff))
        else:
            print('error: {}'.format(line))
            


def filtertags(line):
    return set([i for i in line.strip().split()[1:]]) #if not i.strip('-').isdigit()])#and not i.startswith('MD') and not i.startswith('NM') and not i.startswith('AS') and not i.startswith('SM') 
                #and not i.startswith('OC') and not i.startswith('SA')])# and not i=='a3' and not i=='4a3' and not i=='a1' and not i=='b1' and not i=='4a1' and not i == '4b1' 
                #and not i=='8a1'and not i==8b1' #and not i== 'cb1' and not i == 'c61' and not len(i) == 3])

def main():
    parse_bam_diff()




if __name__ == "__main__":
    main()
