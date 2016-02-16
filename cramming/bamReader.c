#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include <strings.h>



  
//#define MAXSTRINGLEN 100
#define MAXNAMELEN 10000
#define MAXLINES 10000
#define MAXSTRINGLEN 100
//#define MAXTAGS 100
#define MAXTAGS 100

struct tag readtag();
void getstring(char *);
int compare(const void *s1, const void *s2);

union Uvalue {
  int32_t ival;
  char sval[MAXSTRINGLEN];
};

struct tag {
  int len;
  int new_len;
  int order;
  char tag_name[3];
  char tag_value_type;
  char old_tag_value_type;
  union Uvalue value;
};

//struct block {

main()
{
  int ref;
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
  refs = malloc(n_ref * (sizeof(int32_t)*2 + MAXNAMELEN));
  
  int nlines;
  char *stor;
  char *lineptr[MAXLINES];

  stor = malloc(sizeof(char) * l_text + MAXLINES);
  nlines = readlines(lineptr, text, stor, l_text);
  insertion(lineptr,nlines);
  
  


  fwrite(magic,1,4,stdout);
  fwrite(&l_text, sizeof(int32_t),1, stdout);
  //fwrite(text, 1,l_text,stdout);
  writelines(lineptr,nlines);
  fwrite(&n_ref, sizeof(int32_t),1, stdout);
  free(text);
  free(stor);
  
  
  for (ref=0; ref<n_ref; ++ref){
    fread(&l_name, sizeof(int32_t),1, stdin);
    name = malloc(sizeof(char) * l_name);
    fread(name, sizeof(char), l_name,stdin);
    fread(&l_ref, sizeof(int32_t), 1, stdin);
    
    fwrite(&l_name, sizeof(int32_t),1, stdout);
    fwrite(name, sizeof(char), l_name,stdout);
    fwrite(&l_ref, sizeof(int32_t), 1, stdout);
    free(name);
  }
  

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

  int32_t new_block_size;
  int count;
  struct tag new1_tags[MAXTAGS];
  int new_count;
  int n_tags;

  while (1){
  //int i;
  //for (i=0; i<1; i++){
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
    
    count = 0;
    new_count = 0;
    n_tags = 0;
  
    while (bytes_left_in_block != count){
      new1_tags[n_tags] = readtag();
      count += new1_tags[n_tags].len;
      new_count += new1_tags[n_tags].new_len;
      ++n_tags;
  }
    
    new_block_size = block_size_no_tags + new_count; 
    //fprintf(stderr,"block size : %d\nblock size no tags: %d\n bytes_left in block: %d\n newblocksize: %d\ncount: %d\n new count: %d\n", block_size, 
    //	    block_size_no_tags, bytes_left_in_block, new_block_size, count, new_count);
    qsort(new1_tags, n_tags, sizeof(struct tag), compare);

    int i;
   

    fwrite(&new_block_size, sizeof(int32_t), 1, stdout);
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
  
    //for (i = 0; i<n_tags; i++){
    //printtag(new1_tags+i);
    //}
    for (i = 0; i<n_tags; i++){
      outputtag(new1_tags+i);
    }
  
    free(read_name);
    free(cigar);
    free(seq);
    free(qual);
  }
}

struct tag readtag()
{
  struct tag tg;
  int len = 3;
  tg.new_len = 3;
  char tag_name[3] = "\0\0";
  char tag_value_type;
  union Uvalue value;
  int tmp= 0;
  fread(tag_name, sizeof(char), 2, stdin);
  //"AS,SM,RG,NM,BC,OC,SA"
  if (strcmp(tag_name,"SM") == 0)
    tg.order =1;
  if (strcmp(tag_name,"AS") == 0)
    tg.order =2;
  if (strcmp(tag_name,"RG") == 0)
    tg.order =3;
  if (strcmp(tag_name,"NM") == 0)
    tg.order =4;
  if (strcmp(tag_name,"BC") == 0)
    tg.order =5;
  if (strcmp(tag_name,"OC") == 0)
    tg.order =6;
  if (strcmp(tag_name,"SA") == 0)
    tg.order =7;
  fread(&tag_value_type, sizeof(char), 1, stdin);
  switch (tag_value_type){
  case 'Z':
    getstring(value.sval);
//fprintf(stderr, "%s\n", value.sval);
    len += (strlen(value.sval)+1);
    tg.new_len += (strlen(value.sval)+1);
    tg.old_tag_value_type = tag_value_type;
    break;
  case 'c':
    fread(&tmp, sizeof(int8_t), 1, stdin);
    value.ival = tmp;
    len += sizeof(int8_t);
    tg.old_tag_value_type = tag_value_type;
    tag_value_type = 'i';
    tg.new_len += sizeof(int32_t);
    break;
  case 'C':
    fread(&tmp, sizeof(uint8_t), 1, stdin);
    value.ival = tmp;
    len += sizeof(uint8_t);
    tg.old_tag_value_type = tag_value_type;
    tag_value_type = 'i';
    tg.new_len += sizeof(int32_t);
    break;
  case 's':
    fread(&tmp, sizeof(int16_t), 1, stdin);
    value.ival = tmp;
    len += sizeof(int16_t);
    tg.old_tag_value_type = tag_value_type;
    tag_value_type = 'i';
    tg.new_len += sizeof(int32_t);
    break;
  case 'S':
    fread(&tmp, sizeof(uint16_t), 1, stdin);
    value.ival = tmp;
    len += sizeof(uint16_t);
    tg.old_tag_value_type = tag_value_type;
    tag_value_type = 'i';
    tg.new_len += sizeof(int32_t);
    break;
  case 'i':
    fread(&tmp, sizeof(int32_t), 1, stdin);
    value.ival = tmp;
    len += sizeof(int32_t);
    tg.old_tag_value_type = tag_value_type;
    tag_value_type = 'i';
    tg.new_len += sizeof(int32_t);
    break;
  case 'I':
    fread(&tmp, sizeof(uint32_t), 1, stdin);
    value.ival = tmp;
    len += sizeof(uint32_t);
    tg.old_tag_value_type = tag_value_type;
    tag_value_type = 'i';
    tg.new_len += sizeof(int32_t);
    break;
  }
  strcpy(tg.tag_name,tag_name);
  tg.tag_value_type = tag_value_type;
  tg.value = value;
  tg.len = len;
  //fprintf(stderr, "tag len: %d\n ",tg.len);
  //fprintf(stderr, "tag name: ");
  //fwrite(tg.tag_name, sizeof(char),2,stderr);
  //fprintf(stderr, "\n");
  //fprintf(stderr, "tag value type: %c\n ",tg.tag_value_type);
  //if ((tg.tag_value_type) == 'Z'){
  //fprintf(stderr, "tag value: %s\n ",value.sval);
    //fprintf(stderr, "%s\n", value.sval);
  //}
  //else
  //fprintf(stderr, "tag value: %d\n ",value.ival);
  return tg;
}  

printtag(struct tag *t1)
{
  fprintf(stderr, "tag len: %d\n ",t1->len);
  fprintf(stderr, "new tag len: %d\n ",t1->new_len);
  fprintf(stderr, "tag name: ");
  fwrite(t1->tag_name, sizeof(char),2,stderr);
  fprintf(stderr, "\n");
  fprintf(stderr, "tag order %d\n", t1->order);
 //fprintf(stderr, "tag name: %\n ",t1->tag_name);
  fprintf(stderr, "tag value type: %c\n ",t1->tag_value_type);
  fprintf(stderr, "old tag value type: %c\n ",t1->old_tag_value_type);
  
  if ((t1->tag_value_type) == 'Z')
    fprintf(stderr, "tag value: %s\n ",(t1->value).sval);
  else
    fprintf(stderr, "tag value: %d\n ",(t1->value).ival);
  fprintf(stderr, "\n");
}

outputtag(struct tag *t1)
{
  fwrite(t1->tag_name, sizeof(char),2,stdout);
  fprintf(stdout, "%c",t1->tag_value_type);
  if ((t1->tag_value_type) == 'Z')
    fwrite((t1->value).sval,sizeof(char),(t1->len)-3,stdout);
  else{
    //fprintf(stderr, "tag value: %d\n ",(t1->value).ival);
    int32_t i = (t1->value).ival + 0;
    //fprintf(stderr, "%d\n", pt);
    fwrite(&i,sizeof(int32_t),1,stdout);
  }
}


getstring(char *s)
{
  //char string[MAXSTRINGLEN];
  char c;
  //char *s = string;
  while (1){
    c = getchar();
    *s++ = c;
    if (c == '\0')
      break;
  }
  //return string;
}

sortheader()
{
  int32_t l_text;
  char *text;
  int nlines;
  char *stor;
  char *lineptr[MAXLINES];
  fread(&l_text, sizeof(int32_t),1, stdin);
  text = malloc(sizeof(char) * l_text);
  stor = malloc(sizeof(char) * l_text + MAXLINES);
  fread(text, sizeof(char),l_text,stdin);
  nlines = readlines(lineptr, text, stor, l_text);
  fprintf(stderr, "%d\n", nlines);
  //fwrite(stor, 1,l_text+ MAXLINES,stdout);
  insertion(lineptr,nlines);
  writelines(lineptr,nlines);
  //fwrite(text, 1,l_text,stdout);
}

int compare(const void *s1, const void *s2)
{
  char order1[] = "AS,SM,RG,NM,BC,OC,SA";
  char *order = order1;
  struct tag *t1 = (struct tag *)s1;
  struct tag *t2 = (struct tag *)s2;
  
  int tagcompare = t1->order - t2->order;
  return tagcompare;
}


readbam()
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
      fread(&tags[count], sizeof(char), 2, stdin);
      count +=2;
      fread(&tags[count], sizeof(char), 1, stdin);
      val_type = tags[count];
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
      
      
    fwrite(tags, sizeof(char), bytes_left_in_block, stdout);
    free(read_name);
    free(cigar);
    free(seq);
    free(qual);
    free(tags);
  }
}


void insertion(char *v[], int len)
{
  int i, j;
  char *temp;
  for (i = 1; i < len; i++){
    temp = v[i];
    j = i - 1;
    while ((strncmp(temp, v[i], 3) > 0) && (j >= 0)){
      v[j + 1] = v[j];
      j = j - 1;
    }
    v[j + 1] = temp;
  }
}

void qmsort(char *v[],char *s[], int left, int right)
{
  ;
}

void swap(char *v[], int i, int j)
{
  char *temp;

  temp = v[i];
  v[i] = v[j];
  v[j] = temp;
}


#define MAXLEN 1000

int readlines(char *lineptr[], char *text, char *stor, int32_t l_text)
{
  int i, nlines;
  int isnewline;
  nlines = 0;
  //fprintf(stderr,"l_text: %d\n",l_text);
  for (i = 0; i < l_text; i++){
    //fprintf(stderr, "char %d\n", i);
    //
    if (isnewline){
      //fprintf(stderr, "is new line: %d\n", nlines);
      //fprintf(stderr, "first letter: %c\n", text[i]);
      lineptr[nlines] = &stor[i+nlines];
      isnewline = 0;
    }
    //fprintf(stderr, "checkpoint2\n");
    stor[i+nlines] = text[i];
    //if new line assign a pointer to it in the next avalable lineptr pointer array
    if (text[i] == '\n'){
      //fprintf(stderr, "end of line: %d\n", nlines);
      nlines++;
      stor[i+nlines] = '\0';
      isnewline = 1;
    }
    //fprintf(stderr, "checkpoint3\n");
  }
  return nlines;
}


int mgetline(char *s, int lim)
{
  int c;
  char *t=s;
  
  while(--lim)
    ;
}


void writelines(char *lineptr[],int nlines)
{
  int i;
  for(i=0;i<nlines;i++)
    printf("%s",lineptr[i]);
}



