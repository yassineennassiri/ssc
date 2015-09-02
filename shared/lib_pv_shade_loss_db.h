//extern const unsigned char vmpp[12091680];
//extern const unsigned char impp[12091680];
//extern const unsigned char vs[60458400];
//extern const unsigned char is[60458400];
#include <vector>

class DB8
{
public:
	enum db_type{VMPP, IMPP, VS, IS};
	//	DB8() {
	//		init();
	//	};
	~DB8();
	void init();
	short vmpp(size_t ndx){
		return get_vmpp(ndx);
	};
	short impp(size_t ndx){
		return get_impp(ndx);
	};
	short vs(size_t ndx){
		return get_vs(ndx);
	};
	short is(size_t ndx){
		return get_is(ndx);
	};
	std::vector<double> get_vector(const size_t &N, const size_t &d, const size_t &t, const size_t &S, const db_type &DB_TYPE);
	size_t n_choose_k(size_t n, size_t k);
	bool get_index(const size_t &N, const size_t &d, const size_t &t, const size_t &S, const db_type &DB_TYPE, size_t* ret_ndx);
private:
	unsigned char *p_vmpp;
	unsigned char *p_impp;
	unsigned char *p_vs;
	unsigned char *p_is;
	short get_vmpp(size_t i);
	short get_impp(size_t i);
	short get_vs(size_t i);
	short get_is(size_t i);
	bool decompress_file_to_uint8();
	size_t p_vmpp_uint8_size;
	size_t p_impp_uint8_size;
	size_t p_vs_uint8_size;
	size_t p_is_uint8_size;
};
