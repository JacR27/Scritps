import numpy as np
import struct as st
import sys

def readBamHeader(filehandle):
    magic = bh.read(4)
    l_text = bh.read(4)
    l_textU = st.unpack('<i',l_text)[0]
    text = bh.read(l_textU)
    n_ref = bh.read(4)
    n_refU = st.unpack('<i',n_ref)[0]
    refs = []
    for i in range(n_refU):
        l_name = bh.read(4)
        l_nameU = st.unpack('<i',l_name)[0]
        name = bh.read(l_nameU)
        l_ref = bh.read(4)
        refs.append((l_name,name,l_ref))
    return(magic,l_text,n_ref,refs)
    

def main2():
    intSizeMap = {"c":1, "C":1, "s":2, "S":2, "i":4, "I":4 }
    intTypeMap = {"c":"b", "C":"B", "s":"h", "S":"H", "i":"i", "I":"I" }
    out = sys.stdout
    bh = sys.stdin
    
    
    magic = bh.read(4)
    l_text = bh.read(4)
    l_textU = st.unpack('<i',l_text)[0]
    text = bh.read(l_textU)
    n_ref = bh.read(4)
    n_refU = st.unpack('<i',n_ref)[0]
    refs = []
    for i in range(n_refU):
        l_name = bh.read(4)
        l_nameU = st.unpack('<i',l_name)[0]
        name = bh.read(l_nameU)
        l_ref = bh.read(4)
        refs.append((l_name,name,l_ref))
  
    out.write(magic)
    header , l_header = processHeader(text)
    out.write(st.pack('<i',l_header))
    out.write(header)
    out.write(n_ref)
    for ref in refs:
        for i in ref:
            out.write(i)


    while True:
    #for it in range(1):
        block_size = bh.read(4)
        if not block_size:
            return
        block_sizeU = st.unpack('<i',block_size)[0]
        refID_pos = bh.read(4+4)
        bin_mq_nl =  bh.read(4)
        bin_mq_nlU = st.unpack('<I',bin_mq_nl)[0]
        l_read_name = bin_mq_nlU % 2**8
        flag_nc = bh.read(4)
        flag_ncU = st.unpack('<I',flag_nc)[0]
        n_cigar_op = flag_ncU % 2**16
        l_seq = bh.read(4)
        l_seqU = st.unpack('<i',l_seq)[0]
        next_refID_next_pos_tlen_read_name_cigar_seq_qual = bh.read(4+4+4+ l_read_name + 4*n_cigar_op + ((l_seqU+1)/2) + l_seqU)
        block_size_no_tags = (8*4 + l_read_name + 4*n_cigar_op + ((l_seqU+1)/2) + l_seqU)
        bytes_left_in_block = block_sizeU - block_size_no_tags
        
        
        tags = []
        while(bytes_left_in_block):
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
                     assert count < 100, "string didn't end"
                 tags.append((tag,val_type,"".join(Z_string),))
                 bytes_left_in_block -= 3+count
             else:
                 val_size = intSizeMap[val_type]
                 val_ptype = intTypeMap[val_type]
                 value = bh.read(val_size)
                 bytes_left_in_block -= 3+val_size
                 tags.append((tag,val_type,value,))
        tagSize = block_size_no_tags
        tagSize, tags = processTags(tags)
        block_size = st.pack('<i',(block_size_no_tags)+tagSize)


        for i in ([block_size, refID_pos, bin_mq_nl, flag_nc, l_seq, next_refID_next_pos_tlen_read_name_cigar_seq_qual] + [i for tag in tags for i in tag]):
                    out.write(i)



def main(bamfile,outfile):
    intSizeMap = {"c":1, "C":1, "s":2, "S":2, "i":4, "I":4 }
    intTypeMap = {"c":"b", "C":"B", "s":"h", "S":"H", "i":"i", "I":"I" }
    cigarstring = 'MIDNSHP=X'
    sequencestring = '=ACMGRSVTWYHKDBN'
    #out = sys.stdout
    printHead =1
    with open(outfile, 'wb') as out:
        with open(bamfile, 'rb') as bh:
            magic = bh.read(4)
            assert magic == 'BAM\x01'
            l_text = bh.read(4)
            l_textU = st.unpack('<i',l_text)[0]
            text = bh.read(l_textU)
            n_ref = bh.read(4)
            n_refU = st.unpack('<i',n_ref)[0]
            total_name_len = 0
            refs = []
            for i in range(n_refU):
                l_name = bh.read(4)
                l_nameU = st.unpack('<i',l_name)[0]
                name = bh.read(l_nameU)
                l_ref = bh.read(4)
                refs.append((l_name,name,l_ref))
                total_name_len += l_nameU
                
            if printHead:
                out.write(magic)
                header , l_header = processHeader(text)
                out.write(st.pack('<i',l_header))
                out.write(header)
                #out.write(l_text)
                #out.write(text)
                out.write(n_ref)
                for ref in refs:
                    for i in ref:
                        out.write(i)
        
            bytes_in_header = 4 + 4 + l_textU + 4 + n_refU *(4 + 4) + total_name_len
            
            
            
            while True:
 #           for dlkfs in [1]:
                block_size = bh.read(4)
                if not block_size:
                    return
                block_sizeU = st.unpack('<i',block_size)[0]
                refID = bh.read(4)
                pos = bh.read(4)
                bin_mq_nl =  bh.read(4)
                bin_mq_nlU = st.unpack('<I',bin_mq_nl)[0]
                bin, mq_nl = divmod(bin_mq_nlU, 2**16)
                mq, l_read_name = divmod(mq_nl , 2**8)

                flag_nc = bh.read(4)
                flag_ncU = st.unpack('<I',flag_nc)[0]
                flag, n_cigar_op = divmod(flag_ncU , 2**16)
                
                l_seq = bh.read(4)
                l_seqU = st.unpack('<i',l_seq)[0]
                next_refID = bh.read(4)
                next_pos =  bh.read(4)
                tlen =  bh.read(4)
                read_name= bh.read(l_read_name)
                cigar = bh.read(4*n_cigar_op)
                cigarU = st.unpack('<{}i'.format(n_cigar_op),cigar)
                seq = bh.read(((l_seqU+1)/2))
                seqU = st.unpack('{}B'.format((l_seqU+1)/2),seq)
                qual = bh.read(l_seqU)
                qualU = st.unpack("{}B".format(l_seqU),qual)
                
                
                
                bytes_left_in_block = block_sizeU - (8*4 + l_read_name + 4*n_cigar_op + ((l_seqU+1)/2) + l_seqU)
                block_size_no_tags = bytes_left_in_block
                
                
            
                tags = []
                while(bytes_left_in_block):
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
                            assert count < 100, "string didn't end"
                        tags.append((tag,val_type,"".join(Z_string),))
                        bytes_left_in_block -= 3+count
                    else:
                        val_size = intSizeMap[val_type]
                        val_ptype = intTypeMap[val_type]
                        value = bh.read(val_size)
                        valueU = st.unpack("{}{}".format(1,val_ptype), value)[0]
                        bytes_left_in_block -= 3+val_size
                        tags.append((tag,val_type,value,))
                tagSize = block_size_no_tags
                tagSize, tags = processTags(tags)
                block_size = st.pack('<i',(block_sizeU - block_size_no_tags)+tagSize)
                
        
                        
                for i in ([block_size, refID, pos, bin_mq_nl, flag_nc, l_seq, next_refID, next_pos, tlen, read_name, cigar, seq, qual] + [i for tag in tags for i in tag]):
                    out.write(i)

def hellow():
                seqformated = []
                for i in seqU:
                    f,s = divmod(i,2**4)
                    seqformated.extend([sequencestring[f],sequencestring[s]])
        
        
                formatedcigar = []
                for i in cigarU:
                    op_len = divmod(i, 2**4)
                    op_len , op = divmod(i,2**4)
                    op = cigarstring[op]
                    formatedcigar.extend([str(op_len),op])
        
            
            
        
                record = [str(i) for i in [read_name, flag, st.unpack('<I', refID)[0] ,st.unpack('<I', pos)[0], mq, "".join(formatedcigar), 
                                           st.unpack('<I', next_refID)[0], st.unpack('<I', next_pos)[0], st.unpack('<I', tlen)[0], "".join(seqformated),
                                           "".join([chr(i+33) for i in qualU])] + [("{}:{}:".format(*tag)) for tag in tags]]
                print("\t".join(record))
        
    
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

def processHeader(text):
    header = []
    
    for line in text.strip().split("\n"):
        if "scramble" not in line:
            header.append(line)
    #print(header)    
    header = sorted(header,key=lambda line: line[:3])+[""]
    header = "\n".join(header)
    return header, len(header)
            


if __name__ == '__main__':
#    main("bam6.bam","out.bam")
    main2()
