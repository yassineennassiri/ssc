#include "lib_pv_shade_loss_db.h"
#include <stdio.h>
//#include <wx/wx.h>

#include "lib_miniz.h"
// include in only one file in project due to single fileness without header of miniz
//#define MINIZ_HEADER_FILE_ONLY
//#include "../../../Storage/Compression/miniz_v115_r4/tinfl.c"

typedef unsigned char uint8;
typedef unsigned short uint16;
typedef unsigned int uint;

short DB8::get_vmpp(size_t i)
{
	if (i >= 0 && i<6045840) // uint16 check
		return (short)((p_vmpp[2 * i + 1] << 8) | p_vmpp[2 * i]); 
	else 
		return -1;
};

short DB8::get_impp(size_t i)
{ 
	if (i >= 0 && i<6045840) // uint16 check
		return (short)((p_impp[2 * i + 1] << 8) | p_impp[2 * i]); 
	else 
		return -1; 
};

short DB8::get_vs(size_t i)
{ 
	if (i >= 0 && i<30229200)// uint16 check
		return (short)((p_vs[2 * i + 1] << 8) | p_vs[2 * i]); 
	else 
		return -1; 
};

short DB8::get_is(size_t i)
{ 
	if (i >= 0 && i<30229200)// uint16 check
		return (short)((p_is[2 * i + 1] << 8) | p_is[2 * i]); 
	else 
		return -1; 
};

int DB8::get_index(const size_t &N, const size_t &d, const  size_t &t, const size_t &S, const  db_type &DB_TYPE)
{
	size_t ret_ndx=-1;
	size_t length=0, offset=0;
	size_t length_t =10, length_d=10;
	size_t iN = 0, id = 0, it = 0;

	// ret_ndx==0 is an error condition.
	// check N
	if ((N < 1) || (N>8)) return ret_ndx;
	// check d
	if ((d < 1) || (d>10)) return ret_ndx;
	// check t
	if ((t < 1) || (t>10)) return ret_ndx;

	// check S value for validity
	// find number of s vectors
	size_t size_s = n_choose_k(t + N - 1, t);
	if ((S < 1) || (S>size_s)) return ret_ndx;



	switch (DB_TYPE)
	{
		case VMPP:
			length = 8;
			offset = 0;
			break;
		case IMPP:
			length = 8;
			offset = p_vmpp_uint8_size / 2; // short offset
			break;
		case VS:
			length = 40;
			offset = p_vmpp_uint8_size / 2 + p_impp_uint8_size / 2; // short offset
			break;
		case IS:
			length = 40;
			offset = p_vmpp_uint8_size / 2 + p_impp_uint8_size / 2 + p_vs_uint8_size / 2; // short offset
			break;
	}
	if (length == 0) return ret_ndx;
	ret_ndx = 0; // independent vectors for vmpp,impp,vs and is so offset=0

	size_t t_ub = 11; // upper bound of t index for iteration
	size_t d_ub = 10; // upper bound of d index for iteration
	do
	{
		iN++;
		d_ub = ((iN == N) ? d : 10);
		id = 0;
		do
		{
			id++;
			t_ub = (((iN==N) && (id==d)) ? t : 11);
			for (it = 1; it < t_ub; it++)
			{
				// find number of s vectors
				size_s = n_choose_k(it + iN - 1, it);
				// multiply by length of each S vector
				ret_ndx += size_s*length;
			}
		} while (id < d_ub);
	} while (iN < N);
	ret_ndx += (S - 1)*length;
	return ret_ndx;
}

size_t DB8::n_choose_k(size_t n, size_t k)
{
	if (k > n) return 0;
	if (k * 2 > n) k = n - k;
	if (k == 0) return 1;

	size_t result = n;
	for (size_t i = 2; i <= k; ++i) {
		result *= (n - i + 1);
		result /= i;
	}
	return result;
}

std::vector<double> DB8::get_vector(const size_t &N, const size_t &d, const size_t &t, const size_t &S, const db_type &DB_TYPE)
{
	std::vector<double> ret_vec;
	size_t length = 0;
	switch (DB_TYPE)
	{
	case VMPP:
		length = 8;
		break;
	case IMPP:
		length = 8;
		break;
	case VS:
		length = 40;
		break;
	case IS:
		length = 40;
		break;
	}
	if (length == 0) return ret_vec;
	size_t ndx = get_index(N, d, t, S, DB_TYPE);
	for (size_t i = 0; i < length; i++)
	{ // could replace with single get!
		if (DB_TYPE == VMPP)
			ret_vec.push_back((double)get_vmpp(ndx + i) / 1000.0);
		else if (DB_TYPE == IMPP)
			ret_vec.push_back((double)get_impp(ndx + i) / 1000.0);
		else if (DB_TYPE == VS)
			ret_vec.push_back((double)get_vs(ndx + i) / 1000.0);
		else if (DB_TYPE == IS)
			ret_vec.push_back((double)get_is(ndx + i) / 1000.0);
	}
	return ret_vec;
}

void DB8::init()
{
	p_vmpp_uint8_size = 12091680; // uint8 size from matlab
	p_impp_uint8_size = 12091680; // uint8 size from matlab
	p_vmpp = (unsigned char *)malloc(p_vmpp_uint8_size);//malloc(12091680); uint8 size
	p_impp = (unsigned char *)malloc(p_impp_uint8_size);//malloc(12091680); uint8 size
	p_vs_uint8_size = 60458400; // uint8 size from matlab
	p_is_uint8_size = 60458400; // uint8 size from matlab
	p_vs = (unsigned char *)malloc(p_vs_uint8_size);
	p_is = (unsigned char *)malloc(p_is_uint8_size);
	decompress_file_to_uint8();
}

DB8::~DB8()
{
	if (p_vmpp)
		free(p_vmpp);
	if (p_impp)
		free(p_impp);
	if (p_vs)
		free(p_vs);
	if (p_is)
		free(p_is);
}



bool DB8::decompress_file_to_uint8()
{
	size_t status;
	FILE *pInfile;
	size_t infile_size;
	long file_loc;
	uint8 *pCmp_data, *pTmp_data;

	// Open input file.
	pInfile = fopen("DB8_uint8_bin.mgz", "rb");
//	pInfile = fopen("DB8_vmpp_impp_uint8_bin.mgz", "rb");
	if (!pInfile)
	{
		printf("Failed opening input file!\n");
		return EXIT_FAILURE;
	}
	// Determine input file's size.
	fseek(pInfile, 0, SEEK_END);
	file_loc = ftell(pInfile);
	fseek(pInfile, 0, SEEK_SET);

	if ((file_loc < 0) || (file_loc > INT_MAX))
	{
		// This is not a limitation of miniz or tinfl, but this example.
		printf("File is too large to be processed by this example.\n");
		return EXIT_FAILURE;
	}

	infile_size = (uint)file_loc;

	pCmp_data = (uint8 *)malloc(infile_size);
	if (!pCmp_data)
	{
		printf("Out of memory!\n");
		free(pCmp_data);
		return EXIT_FAILURE;
	}
	if (fread(pCmp_data, 1, infile_size, pInfile) != infile_size)
	{
		printf("Failed reading input file!\n");
		free(pCmp_data);
		return EXIT_FAILURE;
	}

	size_t mem_size = p_vmpp_uint8_size + p_impp_uint8_size + p_vs_uint8_size + p_is_uint8_size;

	pTmp_data = (uint8 *)malloc(mem_size);

	status = tinfl_decompress_mem_to_mem((void *)pTmp_data, mem_size, pCmp_data, infile_size, TINFL_FLAG_PARSE_ZLIB_HEADER);

	free(pCmp_data);

	memcpy(p_vmpp, pTmp_data, p_vmpp_uint8_size);
	memcpy(p_impp, pTmp_data + p_vmpp_uint8_size, p_impp_uint8_size);
	memcpy(p_vs, pTmp_data + p_vmpp_uint8_size + p_impp_uint8_size, p_vs_uint8_size);
	memcpy(p_is, pTmp_data + p_vmpp_uint8_size + p_impp_uint8_size + p_vs_uint8_size, p_is_uint8_size);

	free(pTmp_data);

	if (status == TINFL_DECOMPRESS_MEM_TO_MEM_FAILED)
	{
		printf("tinfl_decompress_mem_to_mem() failed with status %i!\n", status);
		return EXIT_FAILURE;
	}

	return true;
};