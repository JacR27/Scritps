#import numpy as np
import struct as st
import sys
import time


def readBamHeader(filehandle):
    magic = filehandle.read(4)
    l_text = filehandle.read(4)
    l_textU = st.unpack('<i',l_text)[0]
    text = filehandle.read(l_textU)
    n_ref = filehandle.read(4)
    n_refU = st.unpack('<i',n_ref)[0]
    refs = []
    for i in range(n_refU):
        l_name = filehandle.read(4)
        l_nameU = st.unpack('<i',l_name)[0]
        name = filehandle.read(l_nameU)
        l_ref = filehandle.read(4)
        refs.append((l_name,name,l_ref))
    return(magic,l_text,text,n_ref,refs)


def processHeader(text):
    header = []    
    for line in text.strip().split("\n"):
        if "scramble" not in line:
            header.append(line)
    header = sorted(header,key=lambda line: line[:3])+[""]
    header = "\n".join(header)
    return header, st.pack('<i',len(header))


def printheader(filehandle ,magic, l_text, text, n_ref, refs):
    filehandle.write(magic)
    filehandle.write(l_text)
    filehandle.write(text)
    filehandle.write(n_ref)
    for ref in refs:
        for i in ref:
            filehandle.write(i)


def processTags(tags):
    intSizeMap = {"c":1, "C":1, "s":2, "S":2, "I":4}
    intTypeMap = {"c":"b", "C":"B", "s":"h", "S":"H", "I":"I"}
    newtags =[]
    tsize = 0
    for tag in tags:
        ttype = intSizeMap.get(tag[1])
        if ttype:
            newValue = st.pack('<i',st.unpack("{}{}".format("<1",intTypeMap[tag[1]]), tag[2])[0])
            tag = (tag[0],"i",newValue,)
        newtags.append(tag)
        tsize += sum([len(i) for i in tag])
    return tsize, sortTags(newtags)


def sortTags(tags):
    tagOrder = {"SM":1, "AS":2, "RG":3, "NM":4, "BC":5, "OC":6, "SA":7}
    return sorted(tags, key=lambda tag: tagOrder[tag[0]])


def readnullterminatedstring(filehandle):
    char = filehandle.read(1)
    stringArray=[char]
    while char !="\0":
        char = filehandle.read(1)
        stringArray.append(char)
    return "".join(stringArray)


def readtags(filehandle,bytes_left_in_block):
    intSizeMap = {"c":1, "C":1, "s":2, "S":2, "i":4, "I":4 }
    intTypeMap = {"c":"b", "C":"B", "s":"h", "S":"H", "i":"i", "I":"I" }
    tags = []
    while(bytes_left_in_block):
        tag = filehandle.read(2)
        val_type = filehandle.read(1)
        if val_type == 'Z':
            value = readnullterminatedstring(filehandle)
            tags.append((tag,val_type,value,))
            bytes_left_in_block -= 3+len(value)
        else:
            val_size = intSizeMap[val_type]
            val_ptype = intTypeMap[val_type]
            value = filehandle.read(val_size)
            bytes_left_in_block -= 3+val_size
            tags.append((tag,val_type,value,))
    return tags


def readblock(filehandle):
        refID_pos = filehandle.read(4+4)
        bin_mq_nl =  filehandle.read(4)
        bin_mq_nlU = st.unpack('<I',bin_mq_nl)[0]
        l_read_name = bin_mq_nlU % 2**8
        flag_nc = filehandle.read(4)
        flag_ncU = st.unpack('<I',flag_nc)[0]
        n_cigar_op = flag_ncU % 2**16
        l_seq = filehandle.read(4)
        l_seqU = st.unpack('<i',l_seq)[0]
        next_refID_next_pos_tlen_read_name_cigar_seq_qual = filehandle.read(4+4+4+ l_read_name + 4*n_cigar_op + ((l_seqU+1)/2) + l_seqU)    
        block_size_no_tags = (8*4 + l_read_name + 4*n_cigar_op + ((l_seqU+1)/2) + l_seqU)
        return refID_pos, bin_mq_nl, flag_nc, l_seq, next_refID_next_pos_tlen_read_name_cigar_seq_qual, block_size_no_tags
    

def printRecord(filehandle, block_size, refID_pos, bin_mq_nl, flag_nc, l_seq, next_refID_next_pos_tlen_read_name_cigar_seq_qual, tags):
    for i in ([block_size, refID_pos, bin_mq_nl, flag_nc, l_seq, next_refID_next_pos_tlen_read_name_cigar_seq_qual] + [i for tag in tags for i in tag]):
        filehandle.write(i)


def main():
    #read, filter sort and print header
    magic, l_text, text, n_refs, refs = readBamHeader(sys.stdin)
    text, l_text = processHeader(text)
    #printheader(sys.stdout, magic, l_text, text,n_refs, refs)
    
    # main loop for bam records, return when run out of data.
    #while False:
    for i in range(1):
        block_size = sys.stdin.read(4)
        if not block_size:
            return
        block_sizeU = st.unpack('<i',block_size)[0]
        
        refID_pos, bin_mq_nl, flag_nc, l_seq, next_refID_next_pos_tlen_read_name_cigar_seq_qual, block_size_no_tags = readblock(sys.stdin)
        bytes_left_in_block = block_sizeU - block_size_no_tags

        tags = readtags(sys.stdin,bytes_left_in_block)
        tagSize, tags = processTags(tags)
        block_size = st.pack('<i',(block_size_no_tags)+tagSize)

        printRecord(sys.stdout, block_size, refID_pos, bin_mq_nl, flag_nc, l_seq, next_refID_next_pos_tlen_read_name_cigar_seq_qual, tags)



if __name__ == '__main__':
    
    start_time = time.time()
    main()
    sys.stderr.write("time elapsed in minutes: {}".format((time.time()-start_time)/60))
