#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <zlib.h>

#define CHUNK 256000
#  define SET_BINARY_MODE(file)

int def(FILE *source, FILE *dest, int level)
{
  int ret, flush;
  unsigned have;
  z_stream strm;
  unsigned char in[CHUNK];
  unsigned char out[CHUNK];
  int windowBits = 15;
  int memLevel = 8;
  
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  ret = deflateInit2(&strm, level, Z_DEFLATED, windowBits, memLevel, Z_HUFFMAN_ONLY);
  if (ret != Z_OK)
    return ret;

  do {
    strm.avail_in = fread(in, 1, CHUNK, source);
    if (ferror(source)) {
      (void)deflateEnd(&strm);
      return Z_ERRNO;
    }
    flush = feof(source) ? Z_FINISH : Z_NO_FLUSH;
    strm.next_in = in;
    
    do{
      strm.avail_out = CHUNK;
      strm.next_out = out;
      ret = deflate(&strm, flush); 
      assert(ret != Z_STREAM_ERROR);
      have = CHUNK - strm.avail_out;
      if (fwrite(out, 1, have, dest) != have || ferror(dest)) {
	(void)deflateEnd(&strm);
	return Z_ERRNO;
      }
    } while (strm.avail_out == 0);
    assert(strm.avail_in == 0);
  } while (flush != Z_FINISH);
  assert (ret == Z_STREAM_END);
  (void)deflateEnd(&strm);
  return Z_OK;
    
}


/* report a zlib or i/o error */
void zerr(int ret)
{
  fputs("zpipe: ", stderr);
  switch (ret) {
  case Z_ERRNO:
    if (ferror(stdin))
      fputs("error reading stdin\n", stderr);
    if (ferror(stdout))
      fputs("error writing stdout\n", stderr);
    break;
  case Z_STREAM_ERROR:
    fputs("invalid compression level\n", stderr);
    break;
  case Z_DATA_ERROR:
    fputs("invalid or incomplete deflate data\n", stderr);
    break;
  case Z_MEM_ERROR:
    fputs("out of memory\n", stderr);
    break;
  case Z_VERSION_ERROR:
    fputs("zlib version mismatch!\n", stderr);
  }
}

int main()
{
  int ret;
  ret = Z_OK;
  ret = def(stdin, stdout, Z_DEFAULT_COMPRESSION);
  printf("compression done\n");
  if (ret != Z_OK)
    zerr(ret);
  return ret;
}