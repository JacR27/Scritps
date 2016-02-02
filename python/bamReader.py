import numpy as np
import struct as st

def main(bamfile):
    intSizeMap = {"c":1, "C":1, "s":2, "S":2, "i":4, "I":4 }
    intTypeMap = {"c":"b", "C":"B", "s":"h", "S":"H", "i":"i", "I":"I" }
    cigarstring = 'MIDNSHP=X'
    sequencestring = '=ACMGRSVTWYHKDBN'
    with open(bamfile, 'rb') as bh:
        magic = bh.read(4)
        assert magic == 'BAM\x01'
        l_text = st.unpack('<i',bh.read(4))[0]
        text = bh.read(l_text)
        n_ref = st.unpack('<i',bh.read(4))[0]
        total_name_len = 0
        for i in range(n_ref):
            l_name = st.unpack('<i',bh.read(4))[0]
            name = bh.read(l_name)
            l_ref = bh.read(4)
            total_name_len += l_name
        bytes_in_header = 4 + 4 + l_text + 4 + n_ref *(4 + 4) + total_name_len
        block_size = st.unpack('<i',bh.read(4))[0]
        refID = st.unpack('<i',bh.read(4))[0]
        pos = st.unpack('<i',bh.read(4))[0]
        bin_mq_nl =  st.unpack('<I',bh.read(4))[0]
        bin, mq_nl = divmod(bin_mq_nl , 2**16)
        mq, l_read_name = divmod(mq_nl , 2**8)
        flag_nc = st.unpack('<I',bh.read(4))[0]
        flag, n_cigar_op = divmod(flag_nc , 2**16)
        l_seq = st.unpack('<i',bh.read(4))[0]
        next_refID = st.unpack('<i',bh.read(4))[0]
        next_pos =  st.unpack('<i',bh.read(4))[0]
        tlen =  st.unpack('<i',bh.read(4))[0]
        read_name= bh.read(l_read_name)
        cigar = st.unpack('<{}i'.format(n_cigar_op),bh.read(4*n_cigar_op))
        seq = st.unpack('{}B'.format(((l_seq+1)/2)),bh.read(((l_seq+1)/2)))
        qual = st.unpack("{}B".format(l_seq),bh.read(l_seq))
        bytes_left_in_block = block_size - (8*4 + l_read_name + 4*n_cigar_op + ((l_seq+1)/2) + l_seq)
        tags = []
        while(bytes_left_in_block > 0):
            tag = bh.read(2)
            val_type = bh.read(1)
            if val_type == 'Z':
                char = bh.read(1)
                Z_string=[char]
                count = 1
                while char !="\0":
                    count +=1
                    char = bh.read(1)
                    Z_string.append(char)
                    if count > 100:
                        break
                tags.append("{}:{}:{}".format(tag,val_type,"".join(Z_string)))
                bytes_left_in_block -= 3+count
                
            else:
                val_size = intSizeMap[val_type]
                val_ptype = intTypeMap[val_type]
                value = st.unpack("{}{}".format(1,val_ptype), bh.read(val_size))[0]
                bytes_left_in_block -= 3+val_size
                tags.append("{}:{}:{}".format(tag,val_type,value))
        

        seqformated = []
        for i in seq:
            f,s = divmod(i,2**4)
            seqformated.extend([sequencestring[f],sequencestring[s]])
        
        
        formatedcigar = []
        for i in cigar:
            op_len = divmod(i, 2**4)
            op_len , op = divmod(i,2**4)
            op = cigarstring[op]
            formatedcigar.extend([str(op_len),op])
        
        
        
        
            record = [str(i) for i in [read_name, flag, refID ,pos, mq, "".join(formatedcigar), next_refID, next_pos, tlen, "".join(seqformated),
                                       "".join([chr(i+33) for i in qual])] + tags]
        print("\t".join(record))
        
    
def sortTags(recored):
    tagOrder = {"SM":1, "AS":2, "RG":3, "NM":4, "BC":5, "OC":6, "SA":7}
    record = recored.strip().split()
    return "\t".join(record[:11]+sorted(record[11:], key=lambda tag: tagOrder[tag[:2]]))

def printHeader(header):
    for i in sorted(header,key=lambda line: line[:3]):
        print(i)

def processRecords():
    for line in sys.stdin:
        print(sortTags(line))

def processHeaderAndFirstRecord():
    header = []
    for line in sys.stdin:
        if line.startswith("@"):
            if "scramble" not in line:
                header.append(line.strip())
        else:
            printHeader(header)
            print(sortTags(line))
            return


if __name__ == '__main__':
    main('bam4.bam')
