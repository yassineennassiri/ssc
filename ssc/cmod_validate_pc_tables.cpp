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
#include <map>
#include <chrono>
#include <random>
//#include <algorithm>
#include "csp_solver_util.h"
#include "sco2_pc_csp_int.h"

static var_info _cm_vtab_validate_pc_tables[] = {
//   VARTYPE        DATATYPE        NAME                    LABEL                                                   UNITS          META  GROUP   REQUIRED_IF  CONSTRAINTS  UI_HINTS*/
    { SSC_INPUT,    SSC_STRING,     "model_name",           "Name of model to test (e.g., 'sco2_recomp_csp_scale')"   "",           "",    "",      "*",     "",       "" },
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

    // Regression vs. basis model comparison metrics
    { SSC_OUTPUT,   SSC_ARRAY,      "dQ_dot_regr_mns_basis",   "Regression model - basis model cycle input heat",      "MWt",       "",    "",              "",     "",       "" },
    { SSC_OUTPUT,   SSC_ARRAY,      "dW_dot_regr_mns_basis",   "Regression model - basis model cycle output power",    "MWe",       "",    "",              "",     "",       "" },
    { SSC_OUTPUT,   SSC_ARRAY,      "pcdQ_dot_regr_mns_basis", "Perc. diff. regr. vs. basis model cycle input heat",   "MWt",       "",    "",              "",     "",       "" },
    { SSC_OUTPUT,   SSC_ARRAY,      "pcdW_dot_regr_mns_basis", "Perc. diff. regr. vs. basis model cycle output power", "MWe",       "",    "",              "",     "",       "" },

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
        // Select and initialize model
        string model_name = as_string("model_name");
        C_sco2_rc_csp_template *mut;
        C_sco2_recomp_csp sco2_recomp_csp_direct;
        C_sco2_recomp_csp_10MWe_scale sco2_recomp_csp_scale;
        if (model_name == "sco2_recomp_csp_direct") {
            mut = &sco2_recomp_csp_direct;
        }
        else if (model_name == "sco2_recomp_csp_scale") {
            mut = &sco2_recomp_csp_scale;
        }
        else {
            throw exec_error("model_under_test", "model name not found");
        }
        
        // Compile model parameters from SSC inputs
        C_sco2_rc_csp_template::S_des_par mut_par;
        if (compile_params(mut_par)) {
            throw exec_error("model_under_test", "error in model parameters");
        };

        // Pass through callback function (with update percent) and pointer
        mut->mf_callback_update = ssc_cmod_update;
        mut->mp_mf_update = (void*)(this);

        // Run design simulation
        int out_type = -1;
        std::string out_msg = "";
        try
        {
            mut->design(mut_par);
        }
        catch (C_csp_exception &csp_exception)
        {
            // Report warning before exiting with error
            while (mut->mc_messages.get_message(&out_type, &out_msg)) {
                log(out_msg);
            }

            throw exec_error("sco2_csp_system", csp_exception.m_error_message);
        }

        // Set design outputs and state points
        output_design_vals(mut);

        if (as_integer("is_generate_udpc") == 0) {
            log("\n Design calculations complete; no off-design cases requested");
            return;
        }

        // Generate regression models, get and output main effect tables and parameters
        util::matrix_t<double> T_htf_parametrics, T_amb_parametrics, m_dot_htf_ND_parametrics;
        if (create_regressions(mut, mut_par, T_htf_parametrics, T_amb_parametrics, m_dot_htf_ND_parametrics)) {
            throw exec_error("model_under_test", "regressions failed");
        }
        output_regressions(T_htf_parametrics, T_amb_parametrics, m_dot_htf_ND_parametrics);  // output main effect tables

        // Generate interaction effect tables
        C_ud_power_cycle custom_pc;
        custom_pc.init(                 // the temperature parameters are in C
            T_htf_parametrics, as_double("T_htf_hot_des"), as_double("T_htf_hot_low"), as_double("T_htf_hot_high"),
            T_amb_parametrics, as_double("T_amb_des"), as_double("T_amb_low"), as_double("T_amb_high"),
            m_dot_htf_ND_parametrics, as_double("m_dot_htf_ND_parametrics"), as_double("m_dot_htf_ND_low"), as_double("m_dot_htf_ND_high"));

        // Generate sample of independent parameters
        // TODO - create orthogonal samples (e.g., Latin Hypercubes)
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator(seed);
        std::uniform_real_distribution<double> distr_T_htf_hot(as_double("T_htf_hot_low"), as_double("T_htf_hot_high"));
        std::uniform_real_distribution<double> distr_T_amb(as_double("T_amb_low"), as_double("T_amb_high"));
        std::uniform_real_distribution<double> distr_m_dot_htf_ND(as_double("m_dot_htf_ND_low"), as_double("m_dot_htf_ND_high"));
        const std::vector<int>::size_type nSamples = 20;
        std::vector<double> T_htf_hot(nSamples), T_amb(nSamples), m_dot_htf_ND(nSamples);
        for (int i = 0; i < nSamples; i++) {
            T_htf_hot.push_back(distr_T_htf_hot(generator));
            T_amb.push_back(distr_T_amb(generator));
            m_dot_htf_ND.push_back(distr_m_dot_htf_ND(generator));
        }
        // Expand samples into a full factorial
        double n_ff = pow(nSamples, 3);
        std::vector<double> T_htf_hot_ff(n_ff), T_amb_ff(n_ff), m_dot_htf_ND_ff(n_ff);
        for (std::vector<int>::size_type i = 0; i != T_htf_hot.size(); i++) {
            for (std::vector<int>::size_type j = 0; j != T_amb.size(); j++) {
                for (std::vector<int>::size_type k = 0; k != m_dot_htf_ND.size(); k++) {
                    T_htf_hot_ff.push_back(T_htf_hot.at(i));
                    T_amb_ff.push_back(T_amb.at(i));
                    m_dot_htf_ND_ff.push_back(m_dot_htf_ND.at(i));
                }
            }
        }

        // Calculate regression model heat and power from sample set
        double Q_dot_des, W_dot_des;
        W_dot_des = as_double("W_dot_net_des");
        Q_dot_des = W_dot_des / as_double("eta_thermal_des");
        std::vector<double> Q_dot_regr_ff(n_ff), W_dot_regr_ff(n_ff);
        for (std::vector<int>::size_type i = 0; i != T_htf_hot_ff.size(); i++) {
            Q_dot_regr_ff.push_back(custom_pc.get_Q_dot_HTF_ND(T_htf_hot_ff.at(i), T_amb_ff.at(i), m_dot_htf_ND_ff.at(i)) * Q_dot_des);
            W_dot_regr_ff.push_back(custom_pc.get_W_dot_gross_ND(T_htf_hot_ff.at(i), T_amb_ff.at(i), m_dot_htf_ND_ff.at(i)) * W_dot_des);
        }

        // Calculate basis model heat and power from sample set
        C_sco2_rc_csp_template::S_od_par mut_od_par;
        const double od_opt_tol = 0.1;                   // doesn't appear to be used
        C_sco2_rc_csp_template::E_off_design_strategies od_strategy = C_sco2_rc_csp_template::E_off_design_strategies::E_MAX_ETA;       // arbitrarily chosen
        const C_sco2_rc_csp_template::S_od_solved *mut_od_sol;
        std::vector<double> Q_dot_basis_ff(n_ff), W_dot_basis_ff(n_ff);
        for (std::vector<int>::size_type i = 0; i != T_htf_hot_ff.size(); i++) {
            mut_od_par.m_T_htf_hot = T_htf_hot_ff.at(i);
            mut_od_par.m_T_amb = T_amb_ff.at(i);
            mut_od_par.m_m_dot_htf = m_dot_htf_ND_ff.at(i) * as_double("m_dot_htf_des");
            mut->off_design_nested_opt(mut_od_par, od_strategy, od_opt_tol);
            mut_od_sol = mut->get_od_solved();
            Q_dot_basis_ff.push_back(mut_od_sol->ms_rc_cycle_od_solved.m_Q_dot);
            W_dot_basis_ff.push_back(mut_od_sol->ms_rc_cycle_od_solved.m_W_dot_net);
        }
        
        // Compare regression and basis models and output metrics
        std::vector<double> dQ_dot_regr_mns_basis(n_ff), dW_dot_regr_mns_basis(n_ff);
        std::vector<double> pcdQ_dot_regr_mns_basis(n_ff), pcdW_dot_regr_mns_basis(n_ff);
        for (std::vector<int>::size_type i = 0; i != Q_dot_regr_ff.size(); i++) {
            dQ_dot_regr_mns_basis.at(i) = Q_dot_regr_ff.at(i) - Q_dot_basis_ff.at(i);
            dW_dot_regr_mns_basis.at(i) = W_dot_regr_ff.at(i) - W_dot_basis_ff.at(i);
            pcdQ_dot_regr_mns_basis.at(i) = dQ_dot_regr_mns_basis.at(i) / Q_dot_basis_ff.at(i) * 100;
            pcdW_dot_regr_mns_basis.at(i) = dW_dot_regr_mns_basis.at(i) / W_dot_basis_ff.at(i) * 100;
        }
        ssc_number_t *dQ_dot_regr_mns_basis_cm = allocate("dQ_dot_regr_mns_basis", n_ff);
        std::copy(dQ_dot_regr_mns_basis.begin(), dQ_dot_regr_mns_basis.end(), dQ_dot_regr_mns_basis_cm);
        ssc_number_t *dW_dot_regr_mns_basis_cm = allocate("dW_dot_regr_mns_basis", n_ff);
        std::copy(dW_dot_regr_mns_basis.begin(), dW_dot_regr_mns_basis.end(), dW_dot_regr_mns_basis_cm);
        ssc_number_t *pcdQ_dot_regr_mns_basis_cm = allocate("pcdQ_dot_regr_mns_basis", n_ff);
        std::copy(pcdQ_dot_regr_mns_basis.begin(), pcdQ_dot_regr_mns_basis.end(), pcdQ_dot_regr_mns_basis_cm);
        ssc_number_t *pcdW_dot_regr_mns_basis_cm = allocate("pcdW_dot_regr_mns_basis", n_ff);
        std::copy(pcdW_dot_regr_mns_basis.begin(), pcdW_dot_regr_mns_basis.end(), pcdW_dot_regr_mns_basis_cm);
        // TODO - output T_htf_hot_ff, T_amb_ff and m_dot_htf_ND_ff
    }

    int compile_params(C_sco2_rc_csp_template::S_des_par &mut_params) {
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
                return 1;
            }
            mut_params.m_UA_recup_tot_des = std::numeric_limits<double>::quiet_NaN();
        }
        else if (mut_params.m_design_method == 2)
        {
            mut_params.m_UA_recup_tot_des = as_double("UA_recup_tot_des");		//[kW/K] Total recuperator conductance
            if (mut_params.m_UA_recup_tot_des < 0.0)
            {
                log("For cycle design method = 2, the input total recuperator conductance must be greater than 0", SSC_ERROR, -1.0);
                return 1;
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

        return 0;
    }

    int create_regressions(C_sco2_rc_csp_template *model, C_sco2_rc_csp_template::S_des_par &mut_params,
        util::matrix_t<double> &T_htf_ptcs, util::matrix_t<double> &T_amb_ptcs, util::matrix_t<double> &m_dot_htf_ptcs) {

        int out_type = -1;
        std::string out_msg = "";

        // Get or calculate user-defined power cycle parameters
        double T_htf_hot_low = mut_params.m_T_htf_hot_in - 273.15 - 20.0;		//[C]
        if (is_assigned("T_htf_hot_low"))
        {
            T_htf_hot_low = as_double("T_htf_hot_low");		//[C]
        }
        assign("T_htf_hot_low", T_htf_hot_low);

        double T_htf_hot_high = mut_params.m_T_htf_hot_in - 273.15 + 15.0;	//[C]
        if (is_assigned("T_htf_hot_high"))
        {
            T_htf_hot_high = as_double("T_htf_hot_high");	//[C]
        }
        assign("T_htf_hot_high", T_htf_hot_high);

        int n_T_htf_hot_in = 5;
        if (is_assigned("n_T_htf_hot"))
        {
            n_T_htf_hot_in = as_integer("n_T_htf_hot");			//[-]
        }
        assign("n_T_htf_hot", n_T_htf_hot_in);

        double T_amb_low = 0.0;
        if (is_assigned("T_amb_low"))
        {
            T_amb_low = as_double("T_amb_low");				//[C]
        }
        assign("T_amb_low", T_amb_low);

        double T_amb_high = std::max(45.0, mut_params.m_T_amb_des - 273.15 + 5.0);
        if (is_assigned("T_amb_high"))
        {
            T_amb_high = as_double("T_amb_high");			//[C]
        }
        assign("T_amb_high", T_amb_high);

        int n_T_amb_in = 10;
        if (is_assigned("n_T_amb"))
        {
            n_T_amb_in = as_integer("n_T_amb");					//[-]
        }
        assign("n_T_amb", n_T_amb_in);

        double m_dot_htf_ND_high = 1.05;
        if (is_assigned("m_dot_htf_ND_high"))
        {
            m_dot_htf_ND_high = as_double("m_dot_htf_ND_high");
        }
        assign("m_dot_htf_ND_high", m_dot_htf_ND_high);

        int n_m_dot_htf_ND_in = 10;
        if (is_assigned("n_m_dot_htf_ND"))
        {
            n_m_dot_htf_ND_in = as_integer("n_m_dot_htf_ND");
        }
        assign("n_m_dot_htf_ND", n_m_dot_htf_ND_in);

        if (n_T_htf_hot_in < 3 || n_T_amb_in < 3 || n_m_dot_htf_ND_in < 3)
        {
            throw exec_error("sco2_csp_ud_pc_tables", "Need at 3 three points for each independent variable");
        }

        double m_dot_htf_ND_low = as_double("m_dot_htf_ND_low");
        try
        {
            model->generate_ud_pc_tables(T_htf_hot_low, T_htf_hot_high, n_T_htf_hot_in,
                T_amb_low, T_amb_high, n_T_amb_in,
                m_dot_htf_ND_low, m_dot_htf_ND_high, n_m_dot_htf_ND_in,
                T_htf_ptcs, T_amb_ptcs, m_dot_htf_ptcs);
        }
        catch (C_csp_exception &csp_exception)
        {
            // Report warning before exiting with error
            while (model->mc_messages.get_message(&out_type, &out_msg))
            {
                log(out_msg);
            }

            throw exec_error("sco2_csp_system", csp_exception.m_error_message);
        }

        // If all calls were successful, log to SSC any messages from sco2_recomp_csp
        while (model->mc_messages.get_message(&out_type, &out_msg))
        {
            log(out_msg);
        }
    }

    int output_regressions(util::matrix_t<double> &T_htf_ptcs, util::matrix_t<double> &T_amb_ptcs, util::matrix_t<double> &m_dot_htf_ptcs) {
        int n_T_htf_hot = (int)T_htf_ptcs.nrows();
        int n_T_amb = (int)T_amb_ptcs.nrows();
        int n_m_dot_htf_ND = (int)m_dot_htf_ptcs.nrows();

        int ncols = (int)T_htf_ptcs.ncols();

        ssc_number_t *p_T_htf_ind = allocate("T_htf_ind", n_T_htf_hot, ncols);
        for (int i = 0; i < n_T_htf_hot; i++)
        {
            for (int j = 0; j < ncols; j++)
            {
                p_T_htf_ind[i*ncols + j] = (ssc_number_t)T_htf_ptcs(i, j);
            }
        }

        ssc_number_t *p_T_amb_ind = allocate("T_amb_ind", n_T_amb, ncols);
        for (int i = 0; i < n_T_amb; i++)
        {
            for (int j = 0; j < ncols; j++)
            {
                p_T_amb_ind[i*ncols + j] = (ssc_number_t)T_amb_ptcs(i, j);
            }
        }

        ssc_number_t *p_m_dot_htf_ND_ind = allocate("m_dot_htf_ND_ind", n_m_dot_htf_ND, ncols);
        for (int i = 0; i < n_m_dot_htf_ND; i++)
        {
            for (int j = 0; j < ncols; j++)
            {
                p_m_dot_htf_ND_ind[i*ncols + j] = (ssc_number_t)m_dot_htf_ptcs(i, j);
            }
        }

        /*
        // Possible cleaner way
        int ncols = (int)indep_1_w_3.ncols();
        int n_indep_1_w_3 = (int)indep_1_w_3.nrows();
        int n_indep_2_w_1 = (int)indep_2_w_1.nrows();
        int n_indep_3_w_2 = (int)indep_3_w_2.nrows();

        ssc_number_t *p_indep_1_w_3 = allocate("indep_1_w_3", n_indep_1_w_3, ncols);
        std::copy(indep_1_w_3.data(), indep_1_w_3.data() + n_indep_1_w_3 * ncols, p_indep_1_w_3);
        ssc_number_t *p_indep_2_w_1 = allocate("indep_2_w_1", n_indep_2_w_1, ncols);
        std::copy(indep_2_w_1.data(), indep_2_w_1.data() + n_indep_2_w_1 * ncols, p_indep_2_w_1);
        ssc_number_t *p_indep_3_w_2 = allocate("indep_3_w_2", n_indep_3_w_2, ncols);
        std::copy(indep_3_w_2.data(), indep_3_w_2.data() + n_indep_3_w_2 * ncols, p_indep_3_w_2);

        // If all calls were successful, log to SSC any messages from sco2_recomp_csp
        while (mut->mc_messages.get_message(&out_type, &out_msg)) {
        log(out_msg);
        }
        log("\n Tested model tables complete");
        */

        return 0;
    }

    int output_design_vals(C_sco2_rc_csp_template *model) {
        double m_dot_htf_design = model->get_phx_des_par()->m_m_dot_hot_des;	//[kg/s]
        double T_htf_cold_calc = model->get_design_solved()->ms_phx_des_solved.m_T_h_out;		//[K]
        assign("T_htf_cold_des", (ssc_number_t)(T_htf_cold_calc - 273.15));		//[C] convert from K
        assign("m_dot_htf_des", (ssc_number_t)m_dot_htf_design);				//[kg/s]
        assign("eta_thermal_calc", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.m_eta_thermal);	//[-]
        assign("m_dot_co2_full", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.m_m_dot_t);		//[kg/s]
        assign("recomp_frac", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.m_recomp_frac);		//[-]
            // Compressor
        assign("P_comp_in", (ssc_number_t)(model->get_design_solved()->ms_rc_cycle_solved.m_pres[1 - 1] / 1000.0));		//[MPa] convert from kPa
        assign("P_comp_out", (ssc_number_t)(model->get_design_solved()->ms_rc_cycle_solved.m_pres[2 - 1] / 1000.0));		//[MPa] convert from kPa
        assign("mc_phi_des", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.ms_mc_ms_des_solved.m_phi_des);
        assign("mc_tip_ratio_des", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.ms_mc_ms_des_solved.m_w_tip_ratio);		//[-]
        assign("mc_n_stages", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.ms_mc_ms_des_solved.m_n_stages);	//[-]
        assign("mc_N_des", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.ms_mc_ms_des_solved.m_N_design);	//[rpm]
        assign("mc_D", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.ms_mc_ms_des_solved.m_D_rotor);			//[m]
        assign("mc_phi_surge", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.ms_mc_ms_des_solved.m_phi_surge);	//[-]
            // Recompressor
        assign("rc_phi_des", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.ms_rc_ms_des_solved.m_phi_des);	//[-]
        assign("rc_tip_ratio_des", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.ms_rc_ms_des_solved.m_w_tip_ratio);	//[-]
        assign("rc_n_stages", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.ms_rc_ms_des_solved.m_n_stages);	//[-]
        assign("rc_N_des", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.ms_rc_ms_des_solved.m_N_design);	//[rpm]
        assign("rc_D", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.ms_rc_ms_des_solved.m_D_rotor);		//[m] 
        assign("rc_phi_surge", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.ms_rc_ms_des_solved.m_phi_surge);//[-]
            // Turbine
        assign("t_nu_des", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.ms_t_des_solved.m_nu_design);           //[-]
        assign("t_tip_ratio_des", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.ms_t_des_solved.m_w_tip_ratio);  //[-]
        assign("t_N_des", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.ms_t_des_solved.m_N_design);			   //[rpm]
        assign("t_D", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.ms_t_des_solved.m_D_rotor);                  //[m]
            // Recuperator
        double UA_LTR = model->get_design_solved()->ms_rc_cycle_solved.m_UA_LT;		//[kW/K]
        double UA_HTR = model->get_design_solved()->ms_rc_cycle_solved.m_UA_HT;		//[kW/K]
        assign("UA_recup_total", (ssc_number_t)(UA_LTR + UA_HTR));		//[kW/K]
        assign("UA_LTR", (ssc_number_t)UA_LTR);						//[kW/K]
        assign("eff_LTR", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.ms_LT_recup_des_solved.m_eff_design);		//[-]
        assign("NTU_LTR", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.ms_LT_recup_des_solved.m_NTU_design);		//[-]
        assign("UA_HTR", (ssc_number_t)UA_HTR);						//[kW/K]
        assign("eff_HTR", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.ms_HT_recup_des_solved.m_eff_design);		//[-]
        assign("NTU_HTR", (ssc_number_t)model->get_design_solved()->ms_rc_cycle_solved.ms_HT_recup_des_solved.m_NTU_design);		//[-]
            // PHX
        assign("UA_PHX", (ssc_number_t)model->get_design_solved()->ms_phx_des_solved.m_UA_design_total);			//[kW/K]
        assign("eff_PHX", (ssc_number_t)model->get_design_solved()->ms_phx_des_solved.m_eff_design);				//[-]
        assign("NTU_PHX", (ssc_number_t)model->get_design_solved()->ms_phx_des_solved.m_NTU_design);				//[-]
            // Air Cooler

            // State Points
        int n_sp = C_RecompCycle::RC_OUT + 1;
        ssc_number_t *p_T_co2_des = allocate("T_co2_des", n_sp);
        ssc_number_t *p_P_co2_des = allocate("P_co2_des", n_sp);
        for (int i = 0; i < n_sp; i++)
        {
            p_T_co2_des[i] = (ssc_number_t)(model->get_design_solved()->ms_rc_cycle_solved.m_temp[i] - 273.15);		//[C]
            p_P_co2_des[i] = (ssc_number_t)(model->get_design_solved()->ms_rc_cycle_solved.m_pres[i] / 1.E3);		//[MPa]
        }

        double sco2_f_min = 0.5;
        if (!model->get_design_solved()->ms_rc_cycle_solved.m_is_rc)
            sco2_f_min = 0.7;

        double m_dot_htf_ND_low = sco2_f_min;
        if (is_assigned("m_dot_htf_ND_low"))
        {
            if (as_boolean("is_apply_default_htf_mins"))
                m_dot_htf_ND_low = std::max(sco2_f_min, as_double("m_dot_htf_ND_low"));	//[-]
            else
                m_dot_htf_ND_low = as_double("m_dot_htf_ND_low");
        }

        assign("m_dot_htf_ND_low", m_dot_htf_ND_low);

        return 0;
    }

    // The following are helper functions for Rankine model tables
    double P_regr(double T_htf_ND, double T_amb_ND, double m_dot_htf_ND, double P_ref, std::map<std::string, int> &r_map,
        const util::matrix_t<double> &T_htf_ptcs, const util::matrix_t<double> &T_amb_ptcs, const util::matrix_t<double> &m_dot_htf_ptcs) {

        // Notation:
        // A = T_htf_ptcs
        // B = T_amb_ptcs
        // C = m_dot_htf_ptcs

        //double P_A, P_AB, P_B, P_BC, P_C;

        //// extract rows for A, B, C and PA, PB, PC, PAC, PAB, PBC
        //util::matrix_t<double> mat;
        //mat = T_htf_ptcs.row(r_map.find("A")->second);
        //std::vector<double> xvec(mat.data(), mat.data() + mat.ncells());
        //mat = T_htf_ptcs.row(r_map.find("PA")->second);
        //std::vector<double> yvec(mat.data(), mat.data() + mat.ncells());
        //P_A = interpolate(xvec, yvec, T_htf_ND, false);

        //P_B = interpolate(T_amb_ptcs.row(r_map.find("B")->second, ))

        //mat = T_htf_ptcs.row(r_map.find("A")->second);
        //std::vector<double> A(mat.data(), mat.data() + mat.ncells());
        //mat = T_htf_ptcs.row(r_map.find("A")->second);
        //std::vector<double> A(mat.data(), mat.data() + mat.ncells());
        //mat = T_htf_ptcs.row(r_map.find("A")->second);
        //std::vector<double> A(mat.data(), mat.data() + mat.ncells());
        //mat = T_htf_ptcs.row(r_map.find("A")->second);
        //std::vector<double> A(mat.data(), mat.data() + mat.ncells());
        //mat = T_htf_ptcs.row(r_map.find("A")->second);
        //std::vector<double> A(mat.data(), mat.data() + mat.ncells());
        //mat = T_htf_ptcs.row(r_map.find("A")->second);
        //std::vector<double> A(mat.data(), mat.data() + mat.ncells());
        //mat = T_htf_ptcs.row(r_map.find("A")->second);
        //std::vector<double> A(mat.data(), mat.data() + mat.ncells());
        //mat = T_htf_ptcs.row(r_map.find("A")->second);
        //std::vector<double> A(mat.data(), mat.data() + mat.ncells());
        //
        //P = P_ref * ( (1 + ))
    }

    std::map<std::string, int> regr_map() {
        std::map<std::string, int> rmap;

        rmap["A"] = 0;       // independent
        rmap["PA"] = 1;
        rmap["QA"] = 2;
        rmap["B"] = 3;       // independent
        rmap["PB"] = 4;
        rmap["QB"] = 5;
        rmap["C"] = 6;       // independent
        rmap["PC"] = 7;
        rmap["QC"] = 8;
        rmap["AC"] = 9;       // independent
        rmap["PAC"] = 10;
        rmap["QAC"] = 11;
        rmap["AB"] = 12;       // independent
        rmap["PAB"] = 13;
        rmap["QAB"] = 14;
        rmap["BC"] = 15;       // independent
        rmap["PBC"] = 16;
        rmap["QBC"] = 17;

        return rmap;
    }

    double interpolate(const std::vector<double> &xData, const std::vector<double> &yData, const double x, const bool extrapolate) {
        int size = xData.size();
        bool ascending;
        xData[1] > xData[0] ? ascending = true : ascending = false;

        int i = 0;                                                                          // find left end of interval for interpolation
        if ((ascending && x >= xData[size - 2]) || (!ascending && x <= xData[size - 2]))   // special case: beyond right end
        {
            i = size - 2;
        }
        else if (ascending)
        {
            while (x > xData[i + 1]) i++;
        }
        else
        {
            while (x < xData[i + 1]) i++;
        }
        double xL = xData[i], yL = yData[i], xR = xData[i + 1], yR = yData[i + 1];      // points on either side (unless beyond ends)
        if (!extrapolate)                                                         // if beyond ends of array and not extrapolating
        {
            if ((ascending && x < xL) || (!ascending && x > xL)) yR = yL;
            if ((ascending && x > xR) || (!ascending && x < xR)) yL = yR;
        }

        double dydx = (yR - yL) / (xR - xL);                                    // gradient

        return yL + dydx * (x - xL);                                              // linear interpolation
    }

};

DEFINE_MODULE_ENTRY(validate_pc_tables, "Regression model validation framework", 0)
