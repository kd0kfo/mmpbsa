#include "Zipper.h"

int mmpbsa_utils::Zipper::zip(const std::string& source_filename, const std::string& dest_filename)
{
	FILE* source = fopen(source_filename.c_str(),"r");
	if(source == NULL)
		throw mmpbsa::ZipperException("zip: could not open " + source_filename,mmpbsa::FILE_IO_ERROR);

	FILE* destination = fopen(dest_filename.c_str(),"w");
	if(destination == NULL)
		throw mmpbsa::ZipperException("zip: could not write to " + dest_filename,mmpbsa::FILE_IO_ERROR);

	int retval = fzip(source,destination,Z_DEFAULT_COMPRESSION);
	fclose(source);
	fclose(destination);
	return retval;
}

int mmpbsa_utils::Zipper::fzip(FILE *source, FILE *dest, int level)
{
    int ret, flush;
    unsigned have;
    z_stream strm;
    unsigned char in[CHUNK];
    unsigned char out[CHUNK];

    /* allocate deflate state */
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    ret = deflateInit2(&strm, level,Z_DEFLATED,ZIPPER_WINDOW_BITS,9,Z_DEFAULT_STRATEGY);
    if (ret != Z_OK)
        return ret;

    /* compress until end of file */
    do {
        strm.avail_in = fread(in, 1, CHUNK, source);
        if (ferror(source)) {
            (void)deflateEnd(&strm);
            return Z_ERRNO;
        }
        flush = feof(source) ? Z_FINISH : Z_NO_FLUSH;
        strm.next_in = in;

        /* run deflate() on input until output buffer not full, finish
           compression if all of source has been read in */
        do {
            strm.avail_out = CHUNK;
            strm.next_out = out;
            ret = deflate(&strm, flush);    /* no bad return value */
            assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
            have = CHUNK - strm.avail_out;
            if (fwrite(out, 1, have, dest) != have || ferror(dest)) {
                (void)deflateEnd(&strm);
                return Z_ERRNO;
            }
        } while (strm.avail_out == 0);
        assert(strm.avail_in == 0);     /* all input will be used */

        /* done when last data in file processed */
    } while (flush != Z_FINISH);
    assert(ret == Z_STREAM_END);        /* stream will be complete */

    /* clean up and return */
    (void)deflateEnd(&strm);
    return Z_OK;
}

int mmpbsa_utils::Zipper::unzip(const std::string& source_filename, const std::string& dest_filename)
{
	FILE* source = fopen(source_filename.c_str(),"r");
	if(source == NULL)
		throw mmpbsa::ZipperException("unzip: could not open " + source_filename,mmpbsa::FILE_IO_ERROR);

	FILE* destination = fopen(dest_filename.c_str(),"w");
	if(destination == NULL)
		throw mmpbsa::ZipperException("unzip: could not write to " + dest_filename,mmpbsa::FILE_IO_ERROR);

	int retval = funzip(source,destination);
	fclose(source);
	fclose(destination);
	return retval;
}

int mmpbsa_utils::Zipper::funzip(FILE *source, FILE *dest)
{
    int ret;
    unsigned have;
    z_stream strm;
    unsigned char in[CHUNK];
    unsigned char out[CHUNK];

    /* allocate inflate state */
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    ret = inflateInit2(&strm,ZIPPER_WINDOW_BITS);
    if (ret != Z_OK)
        return ret;

    /* decompress until deflate stream ends or end of file */
    do {
        strm.avail_in = fread(in, 1, CHUNK, source);
        if (ferror(source)) {
            (void)inflateEnd(&strm);
            return Z_ERRNO;
        }
        if (strm.avail_in == 0)
            break;
        strm.next_in = in;

        /* run inflate() on input until output buffer not full */
        do {
            strm.avail_out = CHUNK;
            strm.next_out = out;
            ret = inflate(&strm, Z_NO_FLUSH);
            assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
            switch (ret) {
            case Z_NEED_DICT:
                ret = Z_DATA_ERROR;     /* and fall through */
            case Z_DATA_ERROR:
            case Z_MEM_ERROR:
                (void)inflateEnd(&strm);
                return ret;
            }
            have = CHUNK - strm.avail_out;
            if (fwrite(out, 1, have, dest) != have || ferror(dest)) {
                (void)inflateEnd(&strm);
                return Z_ERRNO;
            }
        } while (strm.avail_out == 0);

        /* done when inflate() says it's done */
    } while (ret != Z_STREAM_END);

    /* clean up and return */
    (void)inflateEnd(&strm);
    return ret == Z_STREAM_END ? Z_OK : Z_DATA_ERROR;
}

/* report a zlib or i/o error */
void mmpbsa_utils::Zipper::zerr(int ret)
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


void mmpbsa_utils::Zipper::write_oct(char* buffer, const int& number, const int& width)
{
  std::ostringstream buf;
  buf.width(width);
  buf.fill('0');
  buf << std::oct << number;
  strcpy(buffer,buf.str().c_str());

}

void mmpbsa_utils::Zipper::write_dec(char* buffer, const int& number, const int& width)
{
  std::ostringstream buf;
  buf.width(width);
  buf.fill('0');
  buf << number;
  strcpy(buffer,buf.str().c_str());

}

int mmpbsa_utils::Zipper::checksum(const char* header,const size_t& header_size)
{
  int sum = 0;
  for(size_t i = 0;i<header_size;i++)
    {
      sum += header[i];
    }
  return sum;
}

std::stringstream* mmpbsa_utils::Zipper::funtar(FILE* in_file,const std::string& filename)
{
	std::stringstream* returnMe;
	std::istringstream int_buff;
	size_t filesize,amount_read;
	char* in_data;
	std::string block_filename;
	while(!feof(in_file))
	{
		in_data = new char[101];
		fgets(in_data,101,in_file);//extract filename from header
		block_filename = in_data;
		fseek(in_file,24,SEEK_CUR);//jump to file size entry
		delete[] in_data;in_data = new char[13];
		fgets(in_data,13,in_file);//get file size to help and
		int_buff.str(in_data);
		int_buff >> std::oct >> filesize;
		fseek(in_file,376,SEEK_CUR);//jump to file content (header is a 512 byte block)
		if(int_buff.fail())
			throw mmpbsa::ZipperException("mmpbsa_utils::Zipper::untar: Could not read file size from tar file.",mmpbsa::DATA_FORMAT_ERROR);
		if(block_filename == filename)//does it match what we're looking for. Note: prefix (i.e. subdirectory) is included in tar-header entry.
		{
			delete[] in_data;in_data = new char[filesize];
			amount_read = fread(in_data,sizeof(char),filesize,in_file);//actual file contents; thus, don't use fgets because it will append '\0'
			if(amount_read != filesize)
				throw mmpbsa::ZipperException("mmpbsa_utils::Zipper::untar: Could not read entire buffer from tar file.",mmpbsa::DATA_FORMAT_ERROR);
			returnMe = new std::stringstream;
			*returnMe << in_data;
			filesize = 512 - filesize % 512;
			if(filesize != 0)
				fseek(in_file,filesize,SEEK_CUR);//Advance to next tar entry
			delete[] in_data;
			return returnMe;
		}

		fseek(in_file,filesize + (512 - filesize % 512),SEEK_CUR);
	}
	return 0;
}

void mmpbsa_utils::Zipper::tar(const std::vector<std::string>& filenames,FILE* out_file,const std::string& dir,const std::string& prefix)
{

  if(filenames.size() == 0)
	  throw mmpbsa::ZipperException("tar: At least one input file name is needed.");

  if(out_file == 0)
	  throw mmpbsa::ZipperException("tar: the output file pointer is null.");

  using namespace std;
  fstream input_file;

  char header[512];
  std::string prefix_filename, full_filename;
  ostringstream copier;
  vector<string>::const_iterator file = filenames.begin();
  for(;file != filenames.end();file++)
  {
    input_file.open((dir + *file).c_str());
	  if(input_file.fail())
	  {
		  input_file.clear();
		  continue;//not all files actually exist. That's ok.
	  }
	  struct stat buffer;
	  full_filename = dir + *file;
	  int retval = stat(full_filename.c_str(),&buffer);
	  if(retval != 0)
	  {
		  std::ostringstream error;
		  error << "tar: stat error value for " << *file << errno << std::endl;
		  error << strerror(errno) << std::endl;
		  throw mmpbsa::ZipperException(error);
	  }

	  prefix_filename = prefix + *file;
	  copier.str("");
	  copier << input_file.rdbuf();
	  create_header(header,copier.str().c_str(),copier.str().size(),buffer,prefix_filename.c_str());


	  fwrite(header,1,sizeof(header),out_file);
	  fwrite(copier.str().c_str(),1,copier.str().size(),out_file);

	  pad_tarfile(copier.str().c_str(), out_file);
	  input_file.close();
  }
}

void mmpbsa_utils::Zipper::pad_tarfile(const char* data_buffer, FILE* out_file)
{
	size_t buffer_size = sizeof(data_buffer);
	if(buffer_size % 512 != 0)
	{
		char zeros[512 - buffer_size % 512];
		for(size_t i = 0;i<sizeof(zeros);i++)
			zeros[i] = 0;
		fwrite(zeros,1,sizeof(zeros),out_file);
	}
}

struct stat mmpbsa_utils::Zipper::default_stat(const size_t& file_size)
{
	struct stat file_stat;
	file_stat.st_uid = getuid();
	file_stat.st_gid = getgid();
	file_stat.st_size = file_size;
	file_stat.st_mtime = 394179000;//FIX THIS!!! Replace with real time of day stamp
	file_stat.st_mode = S_IRWXU | S_IRGRP | S_IROTH;
	return file_stat;
}

void mmpbsa_utils::Zipper::create_header(char* header, const char* data, const size_t& data_size, const struct stat& file_stat, const char* data_filename)
{
	size_t curr_offset = 0;

	for(size_t i = 0;i<512;i++)
		header[i] = 0;

	strcpy(header,data_filename);
	curr_offset += 100;

	write_oct(header+curr_offset,file_stat.st_mode,7);curr_offset += 8;
	write_oct(header+curr_offset,file_stat.st_uid,7);curr_offset += 8;
	write_oct(header+curr_offset,file_stat.st_gid,7);curr_offset += 8;
	write_oct(header+curr_offset,file_stat.st_size,11);curr_offset += 12;
	write_oct(header+curr_offset,file_stat.st_mtime,11);curr_offset += 12;
	strcpy(header+curr_offset,"        ");curr_offset += 8;
	strcpy(header+curr_offset,"0");curr_offset++;

	strcpy(header+curr_offset," ");

	write_oct(header+148,checksum(header,512),6);
	header[154] = 0;header[155] = 32;

}
