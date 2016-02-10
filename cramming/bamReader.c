#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#define MAXNAMELEN 100
enum tags {SM, AS, RG, NM, BC, OC, SA}


main()
{
  int i;
  char magic[4];
  int32_t l_text;
  char *text;
  int32_t n_ref;
  int32_t l_name;
  char *name;
  int32_t l_ref;
  char *refs;

  fread(magic,sizeof(char),4,stdin);
  fread(&l_text, sizeof(int32_t),1, stdin);
  text = malloc(sizeof(char) * l_text);
  fread(text, sizeof(char),l_text,stdin);
  fread(&n_ref, sizeof(int32_t),1, stdin);
  //refs = malloc(n_ref * (sizeof(int32_t)*2 + MAXNAMELEN));
  
  fwrite(magic,1,4,stdout);
  fwrite(&l_text, sizeof(int32_t),1, stdout);
  fwrite(text, 1,l_text,stdout);
  fwrite(&n_ref, sizeof(int32_t),1, stdout);
  free(text);
  
  for (i=0; i<n_ref; ++i){
    fread(&l_name, sizeof(int32_t),1, stdin);
    name = malloc(sizeof(char) * l_name);
    fread(name, sizeof(char), l_name,stdin);
    fread(&l_ref, sizeof(int32_t), 1, stdin);
    
    fwrite(&l_name, sizeof(int32_t),1, stdout);
    fwrite(name, sizeof(char), l_name,stdout);
    fwrite(&l_ref, sizeof(int32_t), 1, stdout);
    free(name);
  }
  //free(refs);
  int32_t block_size;
  int32_t refID;
  int32_t pos;
  uint32_t bin_mq_nl;
  uint32_t flag_nc;
  int32_t l_seq;
  int32_t next_refID;
  int32_t next_pos;
  int32_t tlen;
  char *read_name;
  uint32_t *cigar;
  uint8_t *seq;
  char *qual;

  char l_read_name;
  int n_cigar_op;
  int block_size_no_tags;
  int bytes_left_in_block;
  int seq_size;

  char val_type;
  char *tags;
  int count;
  int c;

  char *new_tags;
  int n_tags;
  int tag_byte_diff;
  int i;
  while (1){
    fread(&block_size, sizeof(int32_t), 1, stdin);
    if(feof(stdin))
      break;
    fread(&refID, sizeof(int32_t), 1, stdin);
    fread(&pos, sizeof(int32_t), 1, stdin);
    fread(&bin_mq_nl, sizeof(uint32_t), 1, stdin);
    fread(&flag_nc, sizeof(uint32_t), 1, stdin);
    fread(&l_seq, sizeof(int32_t), 1, stdin);
    fread(&next_refID, sizeof(int32_t), 1, stdin);
    fread(&next_pos, sizeof(int32_t), 1, stdin);
    fread(&tlen, sizeof(int32_t), 1, stdin);
    
    l_read_name = bin_mq_nl%256;
    n_cigar_op = flag_nc%65536;
    seq_size = (l_seq+1)/2;
    
    read_name = malloc(sizeof(char) * l_read_name);
    fread(read_name, sizeof(char), l_read_name, stdin);
    
    cigar = malloc(sizeof(uint32_t) * n_cigar_op);
    fread(cigar, sizeof(uint32_t), n_cigar_op, stdin);
    
    seq = malloc(sizeof(uint8_t)*seq_size);
    fread(seq, sizeof(uint8_t), seq_size, stdin);
    
    qual = malloc(sizeof(char) * l_seq);
    fread(qual, sizeof(char), l_seq, stdin);
    
    block_size_no_tags = 32 + l_read_name + 4*n_cigar_op + seq_size + l_seq;
    bytes_left_in_block = block_size - block_size_no_tags;
    
    tags = malloc(bytes_left_in_block);
    
    count = 0;
    n_tags = 0;
    tag_byte_diff = 0;
    //fprintf(stderr,"%d\n",bytes_left_in_block);
    //fprintf(stderr, "l_seq: %d\tl_read_name:%d\tn_cigar_op: %d\t\n", l_seq, l_read_name, n_cigar_op);
    while (bytes_left_in_block != count){
      ++n_tags;
    //for (i = 0; i<10; i++){ 
      //fprintf(stderr,"enterloop\n");
      fread(&tags[count], sizeof(char), 2, stdin);
      //fprintf(stderr,"tags read\n");
      count +=2;
      fread(&tags[count], sizeof(char), 1, stdin);
      val_type = tags[count];
      //fprintf(stderr,"read valew %c\n",val_type);
      count += 1;
      switch (val_type){
      case 'Z':
	while (1){
	  c = getchar();
	  tags[count] = c;
	  count++;
	  if (c == '\0')
	    break;
	}
	break;
      case 'c':
	fread(&tags[count], sizeof(int8_t), 1, stdin);
	count += sizeof(int8_t);
	tag_byte_diff += 3;
	break;
      case 'C':
	fread(&tags[count], sizeof(uint8_t), 1, stdin);
	count += sizeof(uint8_t);
	tag_byte_diff += 3;
	break;
      case 's':
	fread(&tags[count], sizeof(int16_t), 1, stdin);
	count += sizeof(int16_t);
	tag_byte_diff += 2;
	break;
      case 'S':
	fread(&tags[count], sizeof(uint16_t), 1, stdin);
	count += sizeof(uint16_t);
	tag_byte_diff += 2;
	break;
      case 'i':
	fread(&tags[count], sizeof(int32_t), 1, stdin);
	count += sizeof(int32_t);
	break;
      case 'I':
	fread(&tags[count], sizeof(uint32_t), 1, stdin);
	count += sizeof(uint32_t);
	break;
      }
      //fprintf(stderr,"count: %d\n", count);
    }

    new_tags = malloc(bytes_left_in_block+tag_byte_diff);
    

    fwrite(&block_size, sizeof(int32_t), 1, stdout);
    fwrite(&refID, sizeof(int32_t), 1, stdout);
    fwrite(&pos, sizeof(int32_t), 1, stdout);
    fwrite(&bin_mq_nl, sizeof(uint32_t), 1, stdout);
    fwrite(&flag_nc, sizeof(uint32_t), 1, stdout);
    fwrite(&l_seq, sizeof(int32_t), 1, stdout);
    fwrite(&next_refID, sizeof(int32_t), 1, stdout);
    fwrite(&next_pos, sizeof(int32_t), 1, stdout);
    fwrite(&tlen, sizeof(int32_t), 1, stdout);
    fwrite(read_name, sizeof(char), l_read_name, stdout);
    fwrite(cigar, sizeof(uint32_t), n_cigar_op, stdout);
    fwrite(seq, sizeof(uint8_t), seq_size, stdout);
    fwrite(qual, sizeof(char), l_seq, stdout);    
    for (i = 0; i < n_tags; ++i){
      
      
    fwrite(tags, sizeof(char), bytes_left_in_block, stdout);
    free(read_name);
    free(cigar);
    free(seq);
    free(qual);
    free(tags);
  }
  
    
  
  
  

}
