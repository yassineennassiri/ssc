#include "core.h"
#include "lib_util.h"
#include "lib_pv_shade_loss_db.h"

static var_info _cm_vtab_test_pv_shade_loss_db[] = {
/*   VARTYPE           DATATYPE         NAME                           LABEL                                UNITS     META                      GROUP                      REQUIRED_IF                 CONSTRAINTS                      UI_HINTS*/
	{ SSC_INPUT, SSC_NUMBER, "number_strings", "Number of strings (1-8)", "", "", "PV Shade Loss DB", "*", "INTEGER,MIN=1,MAX=8", "" },
	{ SSC_INPUT, SSC_NUMBER, "diffuse_frac", "Diffuse fraction (1-10)", "", "", "PV Shade Loss DB", "*", "INTEGER,MIN=1,MAX=10", "" },
	/* TODO calculated from str_shade_fracs below */
	{ SSC_INPUT, SSC_NUMBER, "max_str_shade", "Maximum shade fraction (1-10)", "", "", "PV Shade Loss DB", "*", "INTEGER,MIN=1,MAX=10", "" },

	/* TODO will be array inputs of length number_strings */
	{ SSC_INPUT, SSC_NUMBER, "str_shade_frac_index", "Index in database based on str_shade_frac (TBD)", "", "", "PV Shade Loss DB", "*", "INTEGER,MIN=1", "" },
//	{ SSC_INPUT, SSC_ARRAY, "str_shade_fracs", "Shading fractions for strings (1-10)", "", "", "PV Shade Loss DB", "*", "length=number_strings", "" },



	{ SSC_OUTPUT, SSC_ARRAY, "vmpp", "Maximum power point voltage", "V", "", "PV Shade Loss DB", "*", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "impp", "Maximum power point current", "A", "", "PV Shade Loss DB", "*", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "vs", "Off power point voltage", "V", "", "PV Shade Loss DB", "*", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "is", "Off power point current", "A", "", "PV Shade Loss DB", "*", "", "" },


var_info_invalid };

class cm_test_pv_shade_loss_db : public compute_module
{
public:

	cm_test_pv_shade_loss_db()
	{
		add_var_info( _cm_vtab_test_pv_shade_loss_db );
	}

	void exec() throw(general_error)
	{
		int N = as_integer("number_strings");
		int d = as_integer("diffuse_frac");
		int t = as_integer("max_str_shade");
		int S = as_integer("str_shade_frac_index");

		
		DB8 db8;
		db8.init();
		std::vector<double>test_vmpp = db8.get_vector(N, d, t, S, DB8::VMPP);
		std::vector<double>test_impp = db8.get_vector(N, d, t, S, DB8::IMPP);
		std::vector<double>test_vs = db8.get_vector(N, d, t, S, DB8::VS);
		std::vector<double>test_is = db8.get_vector(N, d, t, S, DB8::IS);


		ssc_number_t *p_vmpp = allocate("vmpp", test_vmpp.size());
		ssc_number_t *p_impp = allocate("impp", test_impp.size());
		ssc_number_t *p_vs = allocate("vs", test_vs.size());
		ssc_number_t *p_is = allocate("is", test_is.size());

		for (size_t i = 0; i < test_vmpp.size(); i++)
			p_vmpp[i] = test_vmpp[i];
		for (size_t i = 0; i < test_impp.size(); i++)
			p_impp[i] = test_impp[i];
		for (size_t i = 0; i < test_vs.size(); i++)
			p_vs[i] = test_vs[i];
		for (size_t i = 0; i < test_is.size(); i++)
			p_is[i] = test_is[i];


	}
};

DEFINE_MODULE_ENTRY( test_pv_shade_loss_db, "Testing PV Shade Loss DB for strings from Sara MacAlpine", 1 )
