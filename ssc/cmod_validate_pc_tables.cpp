/*******************************************************************************************************
*  Copyright 2017 Alliance for Sustainable Energy, LLC
*
*  NOTICE: This software was developed at least in part by Alliance for Sustainable Energy, LLC
*  (“Alliance”) under Contract No. DE-AC36-08GO28308 with the U.S. Department of Energy and the U.S.
*  The Government retains for itself and others acting on its behalf a nonexclusive, paid-up,
*  irrevocable worldwide license in the software to reproduce, prepare derivative works, distribute
*  copies to the public, perform publicly and display publicly, and to permit others to do so.
*
*  Redistribution and use in source and binary forms, with or without modification, are permitted
*  provided that the following conditions are met:
*
*  1. Redistributions of source code must retain the above copyright notice, the above government
*  rights notice, this list of conditions and the following disclaimer.
*
*  2. Redistributions in binary form must reproduce the above copyright notice, the above government
*  rights notice, this list of conditions and the following disclaimer in the documentation and/or
*  other materials provided with the distribution.
*
*  3. The entire corresponding source code of any redistribution, with or without modification, by a
*  research entity, including but not limited to any contracting manager/operator of a United States
*  National Laboratory, any institution of higher learning, and any non-profit organization, must be
*  made publicly available under this license for as long as the redistribution is made available by
*  the research entity.
*
*  4. Redistribution of this software, without modification, must refer to the software by the same
*  designation. Redistribution of a modified version of this software (i) may not refer to the modified
*  version by the same designation, or by any confusingly similar designation, and (ii) must refer to
*  the underlying software originally provided by Alliance as “System Advisor Model” or “SAM”. Except
*  to comply with the foregoing, the terms “System Advisor Model”, “SAM”, or any confusingly similar
*  designation may not be used to refer to any modified version of this software or any modified
*  version of the underlying software originally provided by Alliance without the prior written consent
*  of Alliance.
*
*  5. The name of the copyright holder, contributors, the United States Government, the United States
*  Department of Energy, or any of their employees may not be used to endorse or promote products
*  derived from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
*  IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
*  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER,
*  CONTRIBUTORS, UNITED STATES GOVERNMENT OR UNITED STATES DEPARTMENT OF ENERGY, NOR ANY OF THEIR
*  EMPLOYEES, BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
*  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
*  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
*  IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
*  THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*******************************************************************************************************/

#include "core.h"
#include "common.h"
//#include <algorithm>
#include "csp_solver_util.h"
#include "sco2_pc_csp_int.h"

static var_info _cm_vtab_validate_pc_tables[] = {
//   VARTYPE        DATATYPE        NAME                    LABEL                                                   UNITS          META  GROUP   REQUIRED_IF  CONSTRAINTS  UI_HINTS*/
    { SSC_INPUT,    SSC_NUMBER,     "htf",                  "Integer code for HTF used in PHX",                       "",           "",    "",      "*",     "",       "" },
    { SSC_INPUT,    SSC_MATRIX,     "htf_props",            "User defined HTF property data",                         "", "7 columns (T,Cp,dens,visc,kvisc,cond,h), at least 3 rows", "", "?=[[0]]", "", "" },
    { SSC_INPUT,    SSC_NUMBER,     "T_htf_hot_des",        "HTF design hot temperature (PHX inlet)",                 "C",          "",    "",      "*",     "",       "" },
    { SSC_INPUT,    SSC_NUMBER,     "dT_PHX_hot_approach",  "Temp diff btw hot HTF and turbine inlet",                "C",          "",    "",      "*",     "",       "" },
    { SSC_INPUT,    SSC_NUMBER,     "T_amb_des",            "Ambient temperature",                                    "C",          "",    "",      "*",     "",       "" },
    { SSC_INPUT,    SSC_NUMBER,     "dT_mc_approach",       "Temp diff btw ambient air and main compressor inlet",    "C",          "",    "",      "*",     "",       "" },
    { SSC_INPUT,    SSC_NUMBER,     "site_elevation",       "Site elevation",                                         "m",          "",    "",      "?=300.0","",      "" },
    { SSC_INPUT,    SSC_NUMBER,     "W_dot_net_des",        "Design cycle power output (no cooling parasitics)",      "MWe",        "",    "",      "*",     "",       "" },
    { SSC_INPUT,    SSC_NUMBER,     "design_method",        "1 = Specify efficiency, 2 = Specify total recup UA",     "",           "",    "",      "?=1",   "",       "" },
    { SSC_INPUT,    SSC_NUMBER,     "eta_thermal_des",      "Power cycle thermal efficiency",                         "",           "",    "",      "?=-1.0","",       "" },
    { SSC_INPUT,    SSC_NUMBER,     "UA_recup_tot_des",     "Total recuperator conductance",                          "kW/K",       "",    "",      "?=-1.0","",       "" },
    { SSC_INPUT,    SSC_NUMBER,     "is_recomp_ok",         "1 = Yes, 0 = simple cycle only",                         "",           "",    "",      "?=1",   "",       "" },
    { SSC_INPUT,    SSC_NUMBER,     "is_PR_fixed",          "0 = No, >0 = fixed pressure ratio",                      "",           "",    "",      "?=0",   "",       "" },
    // Cycle Design
    { SSC_INPUT,    SSC_NUMBER,     "eta_isen_mc",          "Design main compressor isentropic efficiency",           "-",          "",    "",      "*",     "",       "" },
    { SSC_INPUT,    SSC_NUMBER,     "eta_isen_rc",          "Design re-compressor isentropic efficiency",             "-",          "",    "",      "*",     "",       "" },
    { SSC_INPUT,    SSC_NUMBER,     "eta_isen_t",           "Design turbine isentropic efficiency",                   "-",          "",    "",      "*",     "",       "" },
    { SSC_INPUT,    SSC_NUMBER,     "LT_recup_eff_max",     "Maximum allowable effectiveness in LT recuperator",      "-",          "",    "",      "?=1.0", "",       "" },
    { SSC_INPUT,    SSC_NUMBER,     "HT_recup_eff_max",     "Maximum allowable effectiveness in HT recuperator",      "-",          "",    "",      "?=1.0", "",       "" },
    { SSC_INPUT,    SSC_NUMBER,     "P_high_limit",         "High pressure limit in cycle",                           "MPa",        "",    "",      "*",     "",       "" },
    // PHX Design
    { SSC_INPUT,    SSC_NUMBER,     "dT_PHX_cold_approach", "Temp diff btw cold HTF and cold CO2",                    "C",          "",    "",      "?=-1.0","",       "" },
    // Air Cooler Design
    { SSC_INPUT,    SSC_NUMBER,     "fan_power_frac",       "Fraction of net cycle power consumed by air cooler fan", "",           "",    "",      "?=0.01", "",      "" },
    { SSC_INPUT,    SSC_NUMBER,     "deltaP_cooler_frac",   "Fraction of cycle high pres. that is design point cooler CO2 pres. drop", "", "", "", "?=0.002", "",      "" },
    // Off Design UDPC Options
    { SSC_INPUT,    SSC_NUMBER,     "is_generate_udpc",     "1 = generate udpc tables, 0 = only calculate design point cyle", "",   "",    "",      "?=1",    "",      "" },
    { SSC_INPUT,    SSC_NUMBER,     "is_apply_default_htf_mins", "1 = yes (0.5 rc, 0.7 simple), 0 = no, only use 'm_dot_htf_ND_low'", "", "", "",   "?=1",    "",      "" },
    // User Defined Power Cycle Table Inputs
    { SSC_INOUT,    SSC_NUMBER,     "T_htf_hot_low",        "Lower level of HTF hot temperature",					  "C",          "",    "",      "",       "",      "" },
    { SSC_INOUT,    SSC_NUMBER,     "T_htf_hot_high",	    "Upper level of HTF hot temperature",					  "C",		    "",    "",      "",       "",      "" },
    { SSC_INOUT,    SSC_NUMBER,     "n_T_htf_hot",		    "Number of HTF hot temperature parametric runs",		   "",		    "",    "",      "",       "",      "" },
    { SSC_INOUT,    SSC_NUMBER,     "T_amb_low",		    "Lower level of ambient temperature",					  "C",		    "",    "",      "",       "",      "" },
    { SSC_INOUT,    SSC_NUMBER,     "T_amb_high",		    "Upper level of ambient temperature",					  "C",		    "",    "",      "",       "",      "" },
    { SSC_INOUT,    SSC_NUMBER,     "n_T_amb",			    "Number of ambient temperature parametric runs",		   "",		    "",    "",      "",       "",      "" },
    { SSC_INOUT,    SSC_NUMBER,     "m_dot_htf_ND_low",	    "Lower level of normalized HTF mass flow rate",			   "",		    "",    "",      "",       "",      "" },
    { SSC_INOUT,    SSC_NUMBER,     "m_dot_htf_ND_high",    "Upper level of normalized HTF mass flow rate",			   "",		    "",    "",      "",       "",      "" },
    { SSC_INOUT,    SSC_NUMBER,     "n_m_dot_htf_ND",	    "Number of normalized HTF mass flow rate parametric runs", "",		    "",    "",      "",       "",      "" },

    // ** Design OUTPUTS **
    // System Design Solution
    { SSC_OUTPUT,   SSC_NUMBER,     "T_htf_cold_des",       "HTF design cold temperature (PHX outlet)",               "C",          "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "m_dot_htf_des",        "HTF mass flow rate",                                     "kg/s",       "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "eta_thermal_calc",     "Calculated cycle thermal efficiency",                    "-",          "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "m_dot_co2_full",       "CO2 mass flow rate through HTR, PHX, turbine",           "kg/s",       "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "recomp_frac",          "Recompression fraction",                                 "-",          "",    "",      "?=1.2345",     "",       "" },
    // Compressor																															
    { SSC_OUTPUT,   SSC_NUMBER,     "P_comp_in",            "Compressor inlet pressure",                              "MPa",        "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "P_comp_out",           "Compressor outlet pressure",                             "MPa",        "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "mc_phi_des",           "Compressor design flow coefficient",					  "",           "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "mc_tip_ratio_des",     "Compressor design tip speed ratio",                      "",           "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "mc_n_stages",          "Compressor stages",                                      "",           "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "mc_N_des",             "Compressor design shaft speed",                          "rpm",        "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "mc_D",                 "Compressor diameter",                                    "m",          "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "mc_phi_surge",         "Compressor flow coefficient where surge occurs",         "",           "",    "",      "?=1.2345",     "",       "" },
    // Recompressor																															
    { SSC_OUTPUT,   SSC_NUMBER,     "rc_phi_des",           "Recompressor design flow coefficient",                   "",           "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "rc_tip_ratio_des",     "Recompressor 1st stage design tip speed ratio",          "",           "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "rc_n_stages",          "Recompressor stages",                                    "",           "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "rc_N_des",             "Recompressor design shaft speed",                        "rpm",        "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "rc_D",                 "Recompressor first stage diameter",                      "m",          "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "rc_phi_surge",         "Compressor flow coefficient where surge occurs",         "",           "",    "",      "?=1.2345",     "",       "" },
    // Turbine																																
    { SSC_OUTPUT,   SSC_NUMBER,     "t_nu_des",             "Turbine design velocity ratio",                          "",           "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "t_tip_ratio_des",	    "Turbine design tip speed ratio",                         "",           "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "t_N_des",			    "Turbine design shaft speed",	                          "rpm",        "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "t_D",                  "Turbine diameter",                                       "m",          "",    "",      "?=1.2345",     "",       "" },
    // Recuperators																				 											
    { SSC_OUTPUT,   SSC_NUMBER,     "UA_recup_total",       "Total recuperator UA",                                   "kW/K",       "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "UA_LTR",               "Low temp recuperator UA",                                "kW/K",       "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "eff_LTR",              "Low temp recuperator effectiveness",                     "",           "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "NTU_LTR",              "Low temp recuperator NTU",                               "",           "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "UA_HTR",               "High temp recuperator UA",                               "kW/K",       "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "eff_HTR",              "High temp recuperator effectiveness",                    "",           "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "NTU_HTR",              "High temp recuperator NTRU",                             "",           "",    "",      "?=1.2345",     "",       "" },
    // PHX Design Solution																													
    { SSC_OUTPUT,   SSC_NUMBER,     "UA_PHX",               "PHX Conductance",                                        "kW/K",       "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "eff_PHX",              "PHX effectiveness",                                      "",           "",    "",      "?=1.2345",     "",       "" },
    { SSC_OUTPUT,   SSC_NUMBER,     "NTU_PHX",              "PHX NTU",                                                "",           "",    "",      "?=1.2345",     "",       "" },
    // Air Cooler Design
    // ?????
    // State Points
    { SSC_OUTPUT,   SSC_ARRAY,      "T_co2_des",            "Array of cycle CO2 state point temps",	                 "C",          "",    "",      "?=[1.2,2.3,3,4]", "",     "" },
    { SSC_OUTPUT,   SSC_ARRAY,      "P_co2_des",            "Array of cycle CO2 state point pressures",              "MPa",        "",    "",      "?=[1.2,2.3,3,4]", "",     "" },

    // Power Cycle Tables
    { SSC_OUTPUT,   SSC_MATRIX,     "T_htf_ind",            "Parametric of HTF temperature w/ ND HTF mass flow rate levels", "",   "",    "",      "?=[[0,1,2,3,4,5,6,7,8,9,10,11,12][0,1,2,3,4,5,6,7,8,9,10,11,12]]",     "",       "" },
    { SSC_OUTPUT,   SSC_MATRIX,     "T_amb_ind",            "Parametric of ambient temp w/ HTF temp levels",         "",     "",   "",             "?=[[0,1,2,3,4,5,6,7,8,9,10,11,12][0,1,2,3,4,5,6,7,8,9,10,11,12]]",     "",       "" },
    { SSC_OUTPUT,   SSC_MATRIX,     "m_dot_htf_ND_ind",     "Parametric of ND HTF mass flow rate w/ ambient temp levels",    "",   "",    "",      "?=[[0,1,2,3,4,5,6,7,8,9,10,11,12][0,1,2,3,4,5,6,7,8,9,10,11,12]]",     "",       "" },

var_info_invalid };


class cm_validate_pc_tables : public compute_module
{
public:
    cm_validate_pc_tables()
    {
        add_var_info(_cm_vtab_validate_pc_tables);
    }

    void exec() throw(general_error)
    {
        // Set model parameters with SSC inputs
        C_sco2_rc_csp_template::S_des_par mut_par;      // structure holding parameters for model under test
        set_model_params(mut_par);
        
        C_sco2_rc_csp_template *mut;            // model under test, TODO - generalize the type of model
        int out_type = -1;
        std::string out_msg = "";


        util::matrix_t<double> indep_1_w_3, indep_2_w_1, indep_3_w_2;
        
        // Get SSC inputs


        // Setup analytical model
        C_sco2_rc_csp_template::S_des_par mut_des_par;       // TODO - populate structure of parameters (or get it as an input)

        // Run design simulation
        try
        {
            mut->design(mut_des_par);
        }
        catch (C_csp_exception &csp_exception)
        {
            // Report warning before exiting with error
            while (mut->mc_messages.get_message(&out_type, &out_msg)) {
                log(out_msg);
            }
            throw exec_error("model_under_test", csp_exception.m_error_message);
        }

        // Generate regression models
        try
        {
            mut->generate_ud_pc_tables(
                indep_1_lo, indep_1_hi, indep_1_n,
                indep_2_lo, indep_2_hi, indep_2_n,
                indep_3_lo, indep_3_hi, indep_3_n,
                // Outputs:
                indep_1_w_3, indep_2_w_1, indep_3_w_2);
        }
        catch (C_csp_exception &csp_exception)
        {
            while (mut->mc_messages.get_message(&out_type, &out_msg)) {
                log(out_msg);
            }
            throw exec_error("model_under_test", csp_exception.m_error_message);
        }

        // Set SSC outputs
        int ncols = (int)indep_1_w_3.ncols();
        int n_indep_1_w_3 = (int)indep_1_w_3.nrows();
        int n_indep_2_w_1 = (int)indep_2_w_1.nrows();
        int n_indep_3_w_2 = (int)indep_3_w_2.nrows();

        ssc_number_t *p_indep_1_w_3 = allocate("indep_1_w_3", n_indep_1_w_3, ncols);
        std::copy(indep_1_w_3.data(), indep_1_w_3.data() + n_indep_1_w_3*ncols, p_indep_1_w_3);
        ssc_number_t *p_indep_2_w_1 = allocate("indep_2_w_1", n_indep_2_w_1, ncols);
        std::copy(indep_2_w_1.data(), indep_2_w_1.data() + n_indep_2_w_1 * ncols, p_indep_2_w_1);
        ssc_number_t *p_indep_3_w_2 = allocate("indep_3_w_2", n_indep_3_w_2, ncols);
        std::copy(indep_3_w_2.data(), indep_3_w_2.data() + n_indep_3_w_2 * ncols, p_indep_3_w_2);

        // If all calls were successful, log to SSC any messages from sco2_recomp_csp
        while (mut->mc_messages.get_message(&out_type, &out_msg)) {
            log(out_msg);
        }
        log("\n Tested model tables complete");
    }

    void set_model_params(C_sco2_rc_csp_template::S_des_par mut_params) {
        mut_params.m_hot_fl_code = as_integer("htf");							//[-] Integer code for HTF
        mut_params.mc_hot_fl_props = as_matrix("htf_props");					//[-] Custom HTF properties
        mut_params.m_T_htf_hot_in = as_double("T_htf_hot_des") + 273.15;			//[K] Convert from C
        mut_params.m_phx_dt_hot_approach = as_double("dT_PHX_hot_approach");	//[K/C] Temperature difference between hot HTF and turbine CO2 inlet
        mut_params.m_T_amb_des = as_double("T_amb_des") + 273.15;				//[K] Convert from C
        mut_params.m_dt_mc_approach = as_double("dT_mc_approach");				//[K/C] Temperature difference between ambient air and main compressor inlet
        mut_params.m_elevation = as_double("site_elevation");					//[m] Site elevation
        mut_params.m_W_dot_net = as_double("W_dot_net_des")*1000.0;			//[kWe] Convert from MWe, cycle power output w/o cooling parasitics
        mut_params.m_eta_thermal = as_double("eta_thermal_des");				//[-] Cycle thermal efficiency

        mut_params.m_design_method = as_integer("design_method");			//[-] 1 = Specify efficiency, 2 = Specify total recup UA
        if (mut_params.m_design_method == 1)
        {
            mut_params.m_eta_thermal = as_double("eta_thermal_des");				//[-] Cycle thermal efficiency
            if (mut_params.m_eta_thermal < 0.0)
            {
                log("For cycle design method = 1, the input cycle thermal efficiency must be greater than 0", SSC_ERROR, -1.0);
                return;
            }
            mut_params.m_UA_recup_tot_des = std::numeric_limits<double>::quiet_NaN();
        }
        else if (mut_params.m_design_method == 2)
        {
            mut_params.m_UA_recup_tot_des = as_double("UA_recup_tot_des");		//[kW/K] Total recuperator conductance
            if (mut_params.m_UA_recup_tot_des < 0.0)
            {
                log("For cycle design method = 2, the input total recuperator conductance must be greater than 0", SSC_ERROR, -1.0);
                return;
            }
            mut_params.m_eta_thermal = std::numeric_limits<double>::quiet_NaN();
        }
        else
        {
            std::string err_msg = util::format("The input cycle design method, %d, is invalid. It must be 1 or 2.", mut_params.m_design_method);
            log(err_msg, SSC_ERROR, -1.0);
        }

        mut_params.m_is_recomp_ok = as_integer("is_recomp_ok");

        double mc_PR_in = as_double("is_PR_fixed");		//[-]
        if (mc_PR_in != 0.0)
        {
            if (mc_PR_in < 0.0)
            {
                mut_params.m_PR_mc_guess = mc_PR_in * 1.E3;		//[kPa] convert from MPa
            }
            else
            {
                mut_params.m_PR_mc_guess = mc_PR_in;			//[-] Pressure Ratio!
            }
            mut_params.m_fixed_PR_mc = true;
        }
        else
        {
            mut_params.m_PR_mc_guess = std::numeric_limits<double>::quiet_NaN();
            mut_params.m_fixed_PR_mc = false;
        }

        // Cycle design parameters: hardcode pressure drops, for now
        // Define hardcoded sco2 design point parameters
        std::vector<double> DP_LT(2);
        /*(cold, hot) positive values are absolute [kPa], negative values are relative (-)*/
        DP_LT[0] = 0;
        DP_LT[1] = 0;
        /*(cold, hot) positive values are absolute [kPa], negative values are relative (-)*/
        std::vector<double> DP_HT(2);
        DP_HT[0] = 0;
        DP_HT[1] = 0;
        /*(cold, hot) positive values are absolute [kPa], negative values are relative (-)*/
        std::vector<double> DP_PC(2);
        DP_PC[0] = 0;
        DP_PC[1] = 0;
        /*(cold, hot) positive values are absolute [kPa], negative values are relative (-)*/
        std::vector<double> DP_PHX(2);
        DP_PHX[0] = 0;
        DP_PHX[1] = 0;
        mut_params.m_DP_LT = DP_LT;
        mut_params.m_DP_HT = DP_HT;
        mut_params.m_DP_PC = DP_PC;
        mut_params.m_DP_PHX = DP_PHX;
        mut_params.m_N_sub_hxrs = 10;
        mut_params.m_N_turbine = 3600.0;
        mut_params.m_tol = 1.E-3;
        mut_params.m_opt_tol = 1.E-3;

        // Remaining cycle design parameters
        mut_params.m_LT_eff_max = as_double("LT_recup_eff_max");
        mut_params.m_HT_eff_max = as_double("HT_recup_eff_max");
        mut_params.m_eta_mc = as_double("eta_isen_mc");
        mut_params.m_eta_rc = as_double("eta_isen_rc");
        mut_params.m_eta_t = as_double("eta_isen_t");
        mut_params.m_P_high_limit = as_double("P_high_limit")*1000.0;		//[kPa], convert from MPa		

                                                                                // PHX design parameters
        mut_params.m_phx_dt_cold_approach = as_double("dT_PHX_cold_approach");
        if (mut_params.m_phx_dt_cold_approach < 0.0)
        {
            mut_params.m_phx_dt_cold_approach = mut_params.m_phx_dt_hot_approach;
        }

        // Air cooler parameters
        mut_params.m_frac_fan_power = as_double("fan_power_frac");
        mut_params.m_deltaP_cooler_frac = as_double("deltaP_cooler_frac");
    }
};

DEFINE_MODULE_ENTRY(validate_pc_tables, "Regression model validation framework", 0)
