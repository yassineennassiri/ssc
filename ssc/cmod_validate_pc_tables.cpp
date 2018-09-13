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
#include <fstream>
#include <iostream>
#include <iomanip>
//#include <algorithm>
#include "csp_solver_util.h"
#include "csp_common.h"
#include "sco2_pc_csp_int.h"
#include "interpolation_routines.h"
//#include "stochastic.h"

static var_info _cm_vtab_validate_pc_tables[] = {
//   VARTYPE        DATATYPE        NAME                    LABEL                                                   UNITS          META  GROUP   REQUIRED_IF  CONSTRAINTS  UI_HINTS*/
    // System design
	{ SSC_INPUT,    SSC_STRING,     "model_name",           "Name of model to test (e.g., 'sco2_recomp_csp_scale')",  "",           "",    "",      "*",     "",       "" },
    { SSC_INPUT,    SSC_ARRAY,      "indep_levels",         "Levels of the indep. variables T_htf, m_dot, and T_amb", "",           "",    "",      "*",     "",       "" },
    { SSC_INPUT,    SSC_NUMBER,     "sample_type",          "0 = uniform, 1 = random (rect. distr.)",                 "",           "",    "",      "*",     "",       "" },
    { SSC_INPUT,    SSC_NUMBER,     "load_me_tables",       "Load saved main effect tables?",                         "",           "",    "",      "*",     "",       "" },
    { SSC_INPUT,    SSC_NUMBER,     "load_training_data",   "Load training data set from basis model?",               "",           "",    "",      "*",     "",       "" },
    // Cycle Design
    { SSC_INPUT,    SSC_NUMBER,     "cycle_config",         "Cycle configuration, 1=recompression, 2=partial cooling", "-",         "",    "",      "*",     "",       "" },
    // PHX Design
    // Air Cooler Design
    // Off Design UDPC Options
//    { SSC_INPUT,    SSC_NUMBER,     "is_generate_udpc",     "1 = generate udpc tables, 0 = only calculate design point cyle", "",   "",    "",      "?=1",    "",      "" }, //
    { SSC_INPUT,    SSC_NUMBER,     "is_apply_default_htf_mins", "1 = yes (0.5 rc, 0.7 simple), 0 = no, only use 'm_dot_ND_low'", "", "", "",   "?=1",    "",      "" },     //
    // User Defined Power Cycle Table Inputs
    { SSC_INOUT,    SSC_NUMBER,     "T_htf_hot_low",        "Lower level of HTF hot temperature",					  "C",          "",    "",      "",       "",      "" },
    { SSC_INOUT,    SSC_NUMBER,     "T_htf_hot_high",	    "Upper level of HTF hot temperature",					  "C",		    "",    "",      "",       "",      "" },
    { SSC_INOUT,    SSC_NUMBER,     "n_T_htf_hot",		    "Number of HTF hot temperature parametric runs",		   "",		    "",    "",      "",       "",      "" },
    { SSC_INOUT,    SSC_NUMBER,     "T_amb_low",		    "Lower level of ambient temperature",					  "C",		    "",    "",      "",       "",      "" },
    { SSC_INOUT,    SSC_NUMBER,     "T_amb_high",		    "Upper level of ambient temperature",					  "C",		    "",    "",      "",       "",      "" },
    { SSC_INOUT,    SSC_NUMBER,     "n_T_amb",			    "Number of ambient temperature parametric runs",		   "",		    "",    "",      "",       "",      "" },
    { SSC_INOUT,    SSC_NUMBER,     "m_dot_ND_low",	        "Lower level of normalized HTF mass flow rate",			   "",		    "",    "",      "",       "",      "" },
    { SSC_INOUT,    SSC_NUMBER,     "m_dot_ND_high",        "Upper level of normalized HTF mass flow rate",			   "",		    "",    "",      "",       "",      "" },
    { SSC_INOUT,    SSC_NUMBER,     "n_m_dot_ND",	        "Number of normalized HTF mass flow rate parametric runs", "",		    "",    "",      "",       "",      "" },

    // ** Design OUTPUTS **
    // System Design Solution
    // Compressor																															
    // Recompressor																															
    // Turbine																																
    // Recuperators																				 											
    // PHX Design Solution																													
    // Air Cooler Design
    // ?????
    // State Points

    // Power Cycle Tables
    { SSC_INOUT,   SSC_MATRIX,      "T_htf_me",             "Main FX of HTF temperature w/ ND HTF mass flow rate levels", "",       "",    "",      "?=[[0,1,2,3,4,5,6,7,8,9,10,11,12][0,1,2,3,4,5,6,7,8,9,10,11,12]]",     "",       "" },
    { SSC_INOUT,   SSC_MATRIX,      "T_amb_me",             "Main FX of ambient temp w/ HTF temp levels",         "",     "",       "",             "?=[[0,1,2,3,4,5,6,7,8,9,10,11,12][0,1,2,3,4,5,6,7,8,9,10,11,12]]",     "",       "" },
    { SSC_INOUT,   SSC_MATRIX,      "m_dot_ND_me",          "Main FX of ND HTF mass flow rate w/ ambient temp levels",    "",       "",    "",      "?=[[0,1,2,3,4,5,6,7,8,9,10,11,12][0,1,2,3,4,5,6,7,8,9,10,11,12]]",     "",       "" },

    // Interpolation Training Data
    { SSC_INOUT,   SSC_ARRAY,       "T_htf_hot_ff",         "Training data from basis model, T_htf",                      "",       "",    "",               "",     "",       "" },
    { SSC_INOUT,   SSC_ARRAY,       "m_dot_ND_ff",          "Training data from basis model, m_dot",                      "",       "",    "",               "",     "",       "" },
    { SSC_INOUT,   SSC_ARRAY,       "T_amb_ff",             "Training data from basis model, T_amb",                      "",       "",    "",               "",     "",       "" },
    { SSC_INOUT,   SSC_ARRAY,       "Q_dot_basis_ff",       "Training data from basis model, Q_dot",                      "",       "",    "",               "",     "",       "" },
    { SSC_INOUT,   SSC_ARRAY,       "W_dot_basis_ff",       "Training data from basis model, W_dot",                      "",       "",    "",               "",     "",       "" },

    // Regression vs. basis model comparison metrics
    { SSC_OUTPUT,   SSC_ARRAY,      "T_htf_hot_ff",         "Sample of HTF temp. for full-factorial model runs",         "C",       "",    "",              "",     "",       "" },
    { SSC_OUTPUT,   SSC_ARRAY,      "m_dot_ND_ff",          "Sample of mass flow used for full-factorial model runs",     "",       "",    "",              "",     "",       "" },
    { SSC_OUTPUT,   SSC_ARRAY,      "T_amb_ff",             "Sample of ambient temp. for full-factorial model runs",     "C",       "",    "",              "",     "",       "" },

    { SSC_OUTPUT,   SSC_ARRAY,      "Q_dot_basis_ff",       "Basis model cycle input heat",                            "MWt",       "",    "",              "",     "",       "" },
    { SSC_OUTPUT,   SSC_ARRAY,      "Q_dot_regr_ff",        "Regression model cycle input heat",                       "MWt",       "",    "",              "",     "",       "" },
    { SSC_OUTPUT,   SSC_ARRAY,      "Q_dot_interp_ff",      "Interpolation model cycle input heat",                    "MWt",       "",    "",              "",     "",       "" },
    { SSC_OUTPUT,   SSC_ARRAY,      "W_dot_basis_ff",       "Basis model cycle output power",                          "MWe",       "",    "",              "",     "",       "" },
    { SSC_OUTPUT,   SSC_ARRAY,      "W_dot_regr_ff",        "Regression model cycle output power",                     "MWe",       "",    "",              "",     "",       "" },
    { SSC_OUTPUT,   SSC_ARRAY,      "W_dot_interp_ff",      "Interpolation model cycle output power",                  "MWe",       "",    "",              "",     "",       "" },

    //{ SSC_OUTPUT,   SSC_ARRAY,      "dQ_dot_regr_mns_basis",   "Regression model - basis model cycle input heat",      "MWt",       "",    "",              "",     "",       "" },
    //{ SSC_OUTPUT,   SSC_ARRAY,      "dW_dot_regr_mns_basis",   "Regression model - basis model cycle output power",    "MWe",       "",    "",              "",     "",       "" },
    //{ SSC_OUTPUT,   SSC_ARRAY,      "pcdQ_dot_regr_mns_basis", "Perc. diff. regr. vs. basis model cycle input heat",   "MWt",       "",    "",              "",     "",       "" },
    //{ SSC_OUTPUT,   SSC_ARRAY,      "pcdW_dot_regr_mns_basis", "Perc. diff. regr. vs. basis model cycle output power", "MWe",       "",    "",              "",     "",       "" },

var_info_invalid };


class cm_validate_pc_tables : public compute_module
{
public:
    cm_validate_pc_tables()
    {
		add_var_info(vtab_sco2_design);					// add common variables
		add_var_info(_cm_vtab_validate_pc_tables);
    }

    void exec() throw(general_error)
    {
        // Select and initialize basis model
        string model_name = as_string("model_name");
        C_sco2_recomp_csp mut;
        //C_sco2_recomp_csp sco2_recomp_csp_direct;
        //C_sco2_recomp_csp_10MWe_scale sco2_recomp_csp_scale;
        //if (model_name == "sco2_recomp_csp_direct") {
        //    mut = &sco2_recomp_csp_direct;
        //}
        //else if (model_name == "sco2_recomp_csp_scale") {
        //    mut = &sco2_recomp_csp_scale;
        //}
        //else {
        //    throw exec_error("model_under_test", "model name not found");
        //}
        
        // Compile basis model parameters from SSC inputs
        C_sco2_recomp_csp::S_des_par mut_par;
        if (compile_params(mut_par)) {
            throw exec_error("model_under_test", "error in model parameters");
        };

		int sco2_des_err;
        try
        {
			sco2_des_err = sco2_design_cmod_common(this, mut);
        }
        catch (C_csp_exception &csp_exception)
        {
            throw exec_error("sco2_csp_system", "sco2_design_cmod_common failed");
        }
		if (sco2_des_err != 0) { return; }

        // Set design outputs and state points
        output_design_vals(&mut);

        // Obtain training/validation data from basis model
        std::vector<double> T_htf_hot_ff, m_dot_ND_ff, T_amb_ff;
        std::vector<double> Q_dot_basis_ff, W_dot_basis_ff;
        int n_ff;
        if (as_boolean("load_training_data") && is_assigned("T_htf_hot_ff")
            && is_assigned("m_dot_ND_ff") && is_assigned("T_amb_ff")
            && is_assigned("Q_dot_basis_ff") && is_assigned("W_dot_basis_ff"))
        {    
            // Load training data
            T_htf_hot_ff = as_vector_double("T_htf_hot_ff");
            m_dot_ND_ff = as_vector_double("m_dot_ND_ff");
            T_amb_ff = as_vector_double("T_amb_ff");
            Q_dot_basis_ff = as_vector_double("Q_dot_basis_ff");
            W_dot_basis_ff = as_vector_double("W_dot_basis_ff");
            n_ff = T_htf_hot_ff.size();
        }
        else
        {
            // Generate sample set for running basis model
            // TODO - add option for orthogonal samples (e.g., Latin Hypercubes)
            std::vector<int> indep_levels = as_vector_integer("indep_levels");
            assert(indep_levels.size() == 3);
            int sample_type = as_integer("sample_type");
            generate_ff_samples(indep_levels, sample_type, T_htf_hot_ff, m_dot_ND_ff, T_amb_ff);
            n_ff = T_htf_hot_ff.size();

            // Run basis model with sample set
            C_sco2_recomp_csp::S_od_par mut_od_par;
            int od_strategy = C_sco2_recomp_csp::E_TARGET_POWER_ETA_MAX;
            int off_design_code = -1;
            Q_dot_basis_ff.reserve(n_ff), W_dot_basis_ff.reserve(n_ff);

            // Log to file
            std::ofstream log;
            std::cout << std::setprecision(3) << std::fixed;    // 3 decimal places
            log.open("validate_pc_tables_log.dat");
            log << "T_htf_hot [K]" << "\t" << "m_dot_ND [-]" << "\t" << "T_amb [K]" << "\t" << "Q_dot [MWt]" << "\t" << "W_dot [MWe]" << "\n";

            for (std::vector<int>::size_type i = 0; i != T_htf_hot_ff.size(); i++) {
                mut_od_par.m_T_htf_hot = T_htf_hot_ff.at(i) + 273.15;
                mut_od_par.m_m_dot_htf = m_dot_ND_ff.at(i) * as_number("m_dot_des");  // ND -> kg/s
                mut_od_par.m_T_amb = T_amb_ff.at(i) + 273.15;

                // Log to file
                log << T_htf_hot_ff.at(i) << "\t" << m_dot_ND_ff.at(i) << "\t" << T_amb_ff.at(i) << "\t";
                log.flush();

                try {
                    off_design_code = mut.optimize_off_design(mut_od_par, od_strategy);
                    Q_dot_basis_ff.push_back(mut.get_od_solved()->ms_rc_cycle_od_solved.m_Q_dot / 1000.);          // kWt -> MWt
                    W_dot_basis_ff.push_back(mut.get_od_solved()->ms_rc_cycle_od_solved.m_W_dot_net / 1000.);      // kWe -> MWe

                    // Log to file
                    log << Q_dot_basis_ff.at(i) << "\t" << W_dot_basis_ff.at(i) << "\n";
                    log.flush();
                }
                catch(...) {
                    Q_dot_basis_ff.push_back(-999);
                    W_dot_basis_ff.push_back(-999);
                    log << -999 << "\t" << -999 << "\n";
                    log.flush();
                }
            }
            log.close();

            // Save basis model data set
            ssc_number_t *T_htf_hot_ff_cm = allocate("T_htf_hot_ff", T_htf_hot_ff.size());
            std::copy(T_htf_hot_ff.begin(), T_htf_hot_ff.end(), T_htf_hot_ff_cm);
            ssc_number_t *m_dot_ND_ff_cm = allocate("m_dot_ND_ff", m_dot_ND_ff.size());
            std::copy(m_dot_ND_ff.begin(), m_dot_ND_ff.end(), m_dot_ND_ff_cm);
            ssc_number_t *T_amb_ff_cm = allocate("T_amb_ff", T_amb_ff.size());
            std::copy(T_amb_ff.begin(), T_amb_ff.end(), T_amb_ff_cm);
            ssc_number_t *Q_dot_basis_ff_cm = allocate("Q_dot_basis_ff", n_ff);
            std::copy(Q_dot_basis_ff.begin(), Q_dot_basis_ff.end(), Q_dot_basis_ff_cm);
            ssc_number_t *W_dot_basis_ff_cm = allocate("W_dot_basis_ff", n_ff);
            std::copy(W_dot_basis_ff.begin(), W_dot_basis_ff.end(), W_dot_basis_ff_cm);
        }
        // Remove training data where model did not complete (= -999)
        int i = 0;
        while(i < Q_dot_basis_ff.size()) {
            if (Q_dot_basis_ff.at(i) == -999 || W_dot_basis_ff.at(i) == -999) {
                T_htf_hot_ff.erase(T_htf_hot_ff.begin() + i);
                m_dot_ND_ff.erase(m_dot_ND_ff.begin() + i);
                T_amb_ff.erase(T_amb_ff.begin() + i);
                Q_dot_basis_ff.erase(Q_dot_basis_ff.begin() + i);
                W_dot_basis_ff.erase(W_dot_basis_ff.begin() + i);
            }
            else {
                i++;
            }
        }
        n_ff = T_htf_hot_ff.size();


        // Calculate regression model heat and power from sample set and output to SSC
        // Get main effects tables from basis model
        util::matrix_t<double> T_htf_me, T_amb_me, m_dot_ND_me;
        if (!as_boolean("load_me_tables") || !is_assigned("T_htf_me") || !is_assigned("T_amb_me") || !is_assigned("m_dot_ND_me")) {
            // Generate regression models, get and output main effect tables
            if (create_regressions(&mut, mut_par, T_htf_me, T_amb_me, m_dot_ND_me)) {
                //throw exec_error("model_under_test", "regressions failed");
            }
            output_regressions(T_htf_me, T_amb_me, m_dot_ND_me);  // output main effect tables
        }
        else {
            // Load main effects tables from previous validation run to save compute time
            T_htf_me = as_matrix("T_htf_me");
            T_amb_me = as_matrix("T_amb_me");
            m_dot_ND_me = as_matrix("m_dot_ND_me");
        }

        //// Generate interaction effects tables
        //double
        //    T_htf_hot_des = as_double("T_htf_hot_des"),
        //    T_htf_hot_low = as_double("T_htf_hot_low"),
        //    T_htf_hot_high = as_double("T_htf_hot_high"),
        //    T_amb_des = as_double("T_amb_des"),
        //    T_amb_low = as_double("T_amb_low"),
        //    T_amb_high = as_double("T_amb_high"),
        //    m_dot_ND_des = 1.0,
        //    m_dot_ND_low = as_double("m_dot_ND_low"),
        //    m_dot_ND_high = as_double("m_dot_ND_high");
        //C_ud_power_cycle custom_pc;
        //try {
        //    custom_pc.init(                 // the temperature parameters are in C
        //        T_htf_me, T_htf_hot_des, T_htf_hot_low, T_htf_hot_high,
        //        T_amb_me, T_amb_des, T_amb_low, T_amb_high,
        //        m_dot_ND_me, m_dot_ND_des, m_dot_ND_low, m_dot_ND_high);
        //}
        //catch(C_csp_exception csp_e) {
        //    cout << "Exception " << csp_e.m_error_code << " in " << csp_e.m_code_location << ": " << csp_e.m_error_message;
        //}
        //catch (...) {
        //    cout << "Unknown exception initializing custom power cycle for generating interaction effect tables.";
        //}

        //// Evaluate regression model
        //// TODO - change evaluation to use LHS
        //double Q_dot_des, W_dot_des;
        //W_dot_des = as_double("W_dot_net_des");
        //Q_dot_des = W_dot_des / as_double("eta_thermal_des");
        //std::vector<double> Q_dot_regr_ff, W_dot_regr_ff;
        //Q_dot_regr_ff.reserve(n_ff), W_dot_regr_ff.reserve(n_ff);
        //for (std::vector<int>::size_type i = 0; i != T_htf_hot_ff.size(); i++) {
        //    Q_dot_regr_ff.push_back(custom_pc.get_Q_dot_HTF_ND(T_htf_hot_ff.at(i), T_amb_ff.at(i), m_dot_ND_ff.at(i)) * Q_dot_des);     // MWt
        //    W_dot_regr_ff.push_back(custom_pc.get_W_dot_gross_ND(T_htf_hot_ff.at(i), T_amb_ff.at(i), m_dot_ND_ff.at(i)) * W_dot_des);   // MWe
        //}
        //ssc_number_t *Q_dot_regr_ff_cm = allocate("Q_dot_regr_ff", n_ff);
        //std::copy(Q_dot_regr_ff.begin(), Q_dot_regr_ff.end(), Q_dot_regr_ff_cm);
        //ssc_number_t *W_dot_regr_ff_cm = allocate("W_dot_regr_ff", n_ff);
        //std::copy(W_dot_regr_ff.begin(), W_dot_regr_ff.end(), W_dot_regr_ff_cm);


        //// Calculate interpolation model
        //// TODO - add cooling parasitics and water use
        //MatDoub IndepVars;
        //VectDoub Q_dot_interpT, W_dot_interpT;

        //// Populate interpolation training data
        //if (false) {
        //    // From main effect tables:
        //    interp_inputs_from_maineffects(T_htf_me, m_dot_ND_me, T_amb_me, IndepVars, Q_dot_interpT, W_dot_interpT);
        //}
        //else {
        //    // From basis model training data set
        //    IndepVars.reserve(n_ff);
        //    for (std::vector<int>::size_type i = 0; i != n_ff; i++) {
        //        IndepVars.push_back(vector<double>(3, 0.));
        //        IndepVars.back().at(0) = T_htf_hot_ff.at(i);
        //        IndepVars.back().at(1) = m_dot_ND_ff.at(i);
        //        IndepVars.back().at(2) = T_amb_ff.at(i);
        //    }
        //}

        //// Train interpolation model
        //double interp_beta = 1.5;       // try 1.99 too
        //double interp_nug = 0;
        //Powvargram W_dot_vgram(IndepVars, W_dot_basis_ff, interp_beta, interp_nug);                   // W_dot
        //GaussMarkov *W_dot_interp_table = new GaussMarkov(IndepVars, W_dot_basis_ff, W_dot_vgram);
        //Powvargram Q_dot_vgram(IndepVars, Q_dot_basis_ff, interp_beta, interp_nug);                   // Q_dot
        //GaussMarkov *Q_dot_interp_table = new GaussMarkov(IndepVars, Q_dot_basis_ff, Q_dot_vgram);

        //// Evaluate interpolation model using sample set and output values
        //// TODO - THE SAME DATA SET IS BEING USED TO TRAIN **AND** EVALUATE -> change this
        //std::vector<double> Q_dot_interp_ff, W_dot_interp_ff;
        //Q_dot_interp_ff.reserve(n_ff), W_dot_interp_ff.reserve(n_ff);
        //std::vector<double> indep_vars_test(3, 0.);
        //for (std::vector<int>::size_type i = 0; i != T_htf_hot_ff.size(); i++) {
        //    indep_vars_test.at(0) = T_htf_hot_ff.at(i);
        //    indep_vars_test.at(1) = m_dot_ND_ff.at(i);
        //    indep_vars_test.at(2) = T_amb_ff.at(i);

        //    Q_dot_interp_ff.push_back( Q_dot_interp_table->interp(indep_vars_test) ); // MWt
        //    W_dot_interp_ff.push_back( W_dot_interp_table->interp(indep_vars_test) ); // MWe
        //}
        //ssc_number_t *Q_dot_interp_ff_cm = allocate("Q_dot_interp_ff", n_ff);
        //std::copy(Q_dot_interp_ff.begin(), Q_dot_interp_ff.end(), Q_dot_interp_ff_cm);
        //ssc_number_t *W_dot_interp_ff_cm = allocate("W_dot_interp_ff", n_ff);
        //std::copy(W_dot_interp_ff.begin(), W_dot_interp_ff.end(), W_dot_interp_ff_cm);

    }

    int compile_params(C_sco2_recomp_csp::S_des_par &mut_params) {
        mut_params.m_cycle_config = as_integer("cycle_config");     // cycle configuration, 1=recompression, 2=partial cooling

        return 0;
    }

    int create_regressions(C_sco2_recomp_csp *model, C_sco2_recomp_csp::S_des_par &mut_params,
        util::matrix_t<double> &T_htf_ptcs, util::matrix_t<double> &T_amb_ptcs, util::matrix_t<double> &m_dot_ptcs) {

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

        double m_dot_ND_high = 1.05;
        if (is_assigned("m_dot_ND_high"))
        {
            m_dot_ND_high = as_double("m_dot_ND_high");
        }
        assign("m_dot_ND_high", m_dot_ND_high);

        int n_m_dot_ND_in = 10;
        if (is_assigned("n_m_dot_ND"))
        {
            n_m_dot_ND_in = as_integer("n_m_dot_ND");
        }
        assign("n_m_dot_ND", n_m_dot_ND_in);

        if (n_T_htf_hot_in < 3 || n_T_amb_in < 3 || n_m_dot_ND_in < 3)
        {
            throw exec_error("sco2_csp_ud_pc_tables", "Need at 3 three points for each independent variable");
        }

        double m_dot_ND_low = as_double("m_dot_ND_low");
        try
        {
            model->generate_ud_pc_tables(T_htf_hot_low, T_htf_hot_high, n_T_htf_hot_in,
                T_amb_low, T_amb_high, n_T_amb_in,
                m_dot_ND_low, m_dot_ND_high, n_m_dot_ND_in,
                T_htf_ptcs, T_amb_ptcs, m_dot_ptcs);
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

    int output_regressions(util::matrix_t<double> &T_htf_ptcs, util::matrix_t<double> &T_amb_ptcs, util::matrix_t<double> &m_dot_ptcs) {
        int n_T_htf_hot = (int)T_htf_ptcs.nrows();
        int n_T_amb = (int)T_amb_ptcs.nrows();
        int n_m_dot_ND = (int)m_dot_ptcs.nrows();

        int ncols = (int)T_htf_ptcs.ncols();

        ssc_number_t *p_T_htf_me = allocate("T_htf_me", n_T_htf_hot, ncols);
        for (int i = 0; i < n_T_htf_hot; i++)
        {
            for (int j = 0; j < ncols; j++)
            {
                p_T_htf_me[i*ncols + j] = (ssc_number_t)T_htf_ptcs(i, j);
            }
        }

        ssc_number_t *p_T_amb_me = allocate("T_amb_me", n_T_amb, ncols);
        for (int i = 0; i < n_T_amb; i++)
        {
            for (int j = 0; j < ncols; j++)
            {
                p_T_amb_me[i*ncols + j] = (ssc_number_t)T_amb_ptcs(i, j);
            }
        }

        ssc_number_t *p_m_dot_ND_me = allocate("m_dot_ND_me", n_m_dot_ND, ncols);
        for (int i = 0; i < n_m_dot_ND; i++)
        {
            for (int j = 0; j < ncols; j++)
            {
                p_m_dot_ND_me[i*ncols + j] = (ssc_number_t)m_dot_ptcs(i, j);
            }
        }

        // Possible cleaner way
        //int ncols = (int)indep_1_w_3.ncols();
        //int n_indep_1_w_3 = (int)indep_1_w_3.nrows();
        //int n_indep_2_w_1 = (int)indep_2_w_1.nrows();
        //int n_indep_3_w_2 = (int)indep_3_w_2.nrows();

        //ssc_number_t *p_indep_1_w_3 = allocate("indep_1_w_3", n_indep_1_w_3, ncols);
        //std::copy(indep_1_w_3.data(), indep_1_w_3.data() + n_indep_1_w_3 * ncols, p_indep_1_w_3);
        //ssc_number_t *p_indep_2_w_1 = allocate("indep_2_w_1", n_indep_2_w_1, ncols);
        //std::copy(indep_2_w_1.data(), indep_2_w_1.data() + n_indep_2_w_1 * ncols, p_indep_2_w_1);
        //ssc_number_t *p_indep_3_w_2 = allocate("indep_3_w_2", n_indep_3_w_2, ncols);
        //std::copy(indep_3_w_2.data(), indep_3_w_2.data() + n_indep_3_w_2 * ncols, p_indep_3_w_2);

        //// If all calls were successful, log to SSC any messages from sco2_recomp_csp
        //while (mut->mc_messages.get_message(&out_type, &out_msg)) {
        //log(out_msg);
        //}
        //log("\n Tested model tables complete");

        return 0;
    }

    int output_design_vals(C_sco2_recomp_csp *model) {
        if (!is_assigned("eta_thermal_des") || as_number("eta_thermal_des") <= 0.0) assign("eta_thermal_des", as_number("eta_thermal_calc"));
        return 0;
    }

    int generate_ff_samples(std::vector<int> nSamples, int sample_type, std::vector<double> &T_htf_hot_ff,
        std::vector<double> &m_dot_ND_ff, std::vector<double> &T_amb_ff) {
        // For generating full-factorial samples

        double T_htf_hot_des = as_double("T_htf_hot_des");
        double T_htf_hot_low = as_double("T_htf_hot_low");
        double T_htf_hot_high = as_double("T_htf_hot_high");
        double m_dot_ND_des = 1.0;
        double m_dot_ND_low = as_double("m_dot_ND_low");
        double m_dot_ND_high = as_double("m_dot_ND_high");
        double T_amb_des = as_double("T_amb_des");
        double T_amb_low = as_double("T_amb_low");
        double T_amb_high = as_double("T_amb_high");

        int n_T_htf_hot = nSamples.at(0);
        int n_m_dot_ND = nSamples.at(1);
        int n_T_amb = nSamples.at(2);
        std::vector<double> T_htf_hot(n_T_htf_hot), m_dot_ND(n_m_dot_ND), T_amb(n_T_amb);

        if (sample_type == 0) {             // uniform sample
            double d, inc;
            inc = (T_htf_hot_high - T_htf_hot_low) / max(1, (n_T_htf_hot - 1));        // fill T_htf_hot
            d = T_htf_hot_low - inc;      // subtract inc so first value equals T_htf_hot_low
            generate(T_htf_hot.begin(), T_htf_hot.end(), [&d, inc] { return d += inc; });
            inc = (m_dot_ND_high - m_dot_ND_low) / max(1, (n_m_dot_ND - 1));           // fill m_dot_ND_low
            d = m_dot_ND_low - inc;
            generate(m_dot_ND.begin(), m_dot_ND.end(), [&d, inc] { return d += inc; });
            inc = (T_amb_high - T_amb_low) / max(1, (n_T_amb - 1));                    // fill T_amb
            d = T_amb_low - inc;
            generate(T_amb.begin(), T_amb.end(), [&d, inc] { return d += inc; });
        }
        else if (sample_type == 1) {        // random from rectangular distribution
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine generator(seed);
            std::uniform_real_distribution<double> distr_T_htf_hot(T_htf_hot_low, T_htf_hot_high);
            std::uniform_real_distribution<double> distr_m_dot_ND(m_dot_ND_low, m_dot_ND_high);
            std::uniform_real_distribution<double> distr_T_amb(T_amb_low, T_amb_high);
            for (int i = 0; i < n_T_htf_hot; i++) { T_htf_hot.at(i) = distr_T_htf_hot(generator); }
            for (int i = 0; i < n_m_dot_ND; i++) { m_dot_ND.at(i) = distr_m_dot_ND(generator); }
            for (int i = 0; i < n_T_amb; i++) { T_amb.at(i) = distr_T_amb(generator); }
        }
        else {
            throw(C_csp_exception("No such sample type"));
        }
        // Expand samples into a full factorial
        double n_ff = n_T_htf_hot * n_m_dot_ND * n_T_amb;
        T_htf_hot_ff.reserve(n_ff), m_dot_ND_ff.reserve(n_ff), T_amb_ff.reserve(n_ff);
        for (std::vector<int>::size_type i = 0; i != T_htf_hot.size(); i++) {
            for (std::vector<int>::size_type j = 0; j != m_dot_ND.size(); j++) {
                for (std::vector<int>::size_type k = 0; k != T_amb.size(); k++) {
                    T_htf_hot_ff.push_back(T_htf_hot.at(i));
                    m_dot_ND_ff.push_back(m_dot_ND.at(j));
                    T_amb_ff.push_back(T_amb.at(k));
                }
            }
        }
        // Output T_htf_hot, T_amb and m_dot_ND values
        ssc_number_t *T_htf_hot_ff_cm = allocate("T_htf_hot_ff", n_ff);
        std::copy(T_htf_hot_ff.begin(), T_htf_hot_ff.end(), T_htf_hot_ff_cm);
        ssc_number_t *m_dot_ND_ff_cm = allocate("m_dot_ND_ff", n_ff);
        std::copy(m_dot_ND_ff.begin(), m_dot_ND_ff.end(), m_dot_ND_ff_cm);
        ssc_number_t *T_amb_ff_cm = allocate("T_amb_ff", n_ff);
        std::copy(T_amb_ff.begin(), T_amb_ff.end(), T_amb_ff_cm);
    }

    int interp_inputs_from_maineffects(const util::matrix_t<double> &T_htf_me, const util::matrix_t<double> &m_dot_ND_me,
        const util::matrix_t<double> &T_amb_me, MatDoub &IndepVars, VectDoub &Q_dot_me, VectDoub &W_dot_me) {
        
        double
            T_htf_hot_des = as_double("T_htf_hot_des"),
            T_htf_hot_low = as_double("T_htf_hot_low"),
            T_htf_hot_high = as_double("T_htf_hot_high"),
            T_amb_des = as_double("T_amb_des"),
            T_amb_low = as_double("T_amb_low"),
            T_amb_high = as_double("T_amb_high"),
            m_dot_ND_des = 1.0,
            m_dot_ND_low = as_double("m_dot_ND_low"),
            m_dot_ND_high = as_double("m_dot_ND_high");

        // Main Effect table structure:
        //    Independent |    Gross Power Output   |   HTF Thermal Power	|   Cooling Parasitics  |	 Water Use 
        // 0)  Variable   |  1) -   2) 0     3) +   |  4) -   5) 0    6) +  |  7) -    8) 0    9) + | 10) -  11) 0   12) +

        int nObs = 3 * (T_htf_me.nrows() + m_dot_ND_me.nrows() + T_amb_me.nrows());
        IndepVars.reserve(nObs), Q_dot_me.reserve(nObs), W_dot_me.reserve(nObs);

        // Populate from T_htf table
        for (int i = 0; i < T_htf_me.nrows(); i++) {
            // Low-level
            IndepVars.push_back(vector<double>(3, 0.));     // [T_htf, m_dot, T_amb]
            IndepVars.back().at(0) = T_htf_me.at(i, 0);
            IndepVars.back().at(1) = m_dot_ND_low;
            IndepVars.back().at(2) = T_amb_des;
            W_dot_me.push_back(T_htf_me.at(i, 1));
            Q_dot_me.push_back(T_htf_me.at(i, 4));
            // Design-level
            IndepVars.push_back(vector<double>(3, 0.));
            IndepVars.back().at(0) = T_htf_me.at(i, 0);
            IndepVars.back().at(1) = m_dot_ND_des;
            IndepVars.back().at(2) = T_amb_des;
            W_dot_me.push_back(T_htf_me.at(i, 2));
            Q_dot_me.push_back(T_htf_me.at(i, 5));
            // High-level
            IndepVars.push_back(vector<double>(3, 0.));
            IndepVars.back().at(0) = T_htf_me.at(i, 0);
            IndepVars.back().at(1) = m_dot_ND_high;
            IndepVars.back().at(2) = T_amb_des;
            W_dot_me.push_back(T_htf_me.at(i, 3));
            Q_dot_me.push_back(T_htf_me.at(i, 6));
        }

        // Populate from m_dot_ND table
        for (int i = 0; i < m_dot_ND_me.nrows(); i++) {
            // Low-level
            IndepVars.push_back(vector<double>(3, 0.));     // [T_htf, m_dot, T_amb]
            IndepVars.back().at(0) = T_htf_hot_des;
            IndepVars.back().at(1) = m_dot_ND_me.at(i, 0);
            IndepVars.back().at(2) = T_amb_low;
            W_dot_me.push_back(m_dot_ND_me.at(i, 1));
            Q_dot_me.push_back(m_dot_ND_me.at(i, 4));
            // Design-level
            IndepVars.push_back(vector<double>(3, 0.));
            IndepVars.back().at(0) = T_htf_hot_des;
            IndepVars.back().at(1) = m_dot_ND_me.at(i, 0);
            IndepVars.back().at(2) = T_amb_des;
            W_dot_me.push_back(m_dot_ND_me.at(i, 2));
            Q_dot_me.push_back(m_dot_ND_me.at(i, 5));
            // High-level
            IndepVars.push_back(vector<double>(3, 0.));
            IndepVars.back().at(0) = T_htf_hot_des;
            IndepVars.back().at(1) = m_dot_ND_me.at(i, 0);
            IndepVars.back().at(2) = T_amb_high;
            W_dot_me.push_back(m_dot_ND_me.at(i, 3));
            Q_dot_me.push_back(m_dot_ND_me.at(i, 6));
        }

        // Populate from T_amb table
        for (int i = 0; i < T_amb_me.nrows(); i++) {
            // Low-level
            IndepVars.push_back(vector<double>(3, 0.));     // [T_htf, m_dot, T_amb]
            IndepVars.back().at(0) = T_htf_hot_low;
            IndepVars.back().at(1) = m_dot_ND_des;
            IndepVars.back().at(2) = T_amb_me.at(i, 0);
            W_dot_me.push_back(T_amb_me.at(i, 1));
            Q_dot_me.push_back(T_amb_me.at(i, 4));
            // Design-level
            IndepVars.push_back(vector<double>(3, 0.));
            IndepVars.back().at(0) = T_htf_hot_des;
            IndepVars.back().at(1) = m_dot_ND_des;
            IndepVars.back().at(2) = T_amb_me.at(i, 0);
            W_dot_me.push_back(T_amb_me.at(i, 2));
            Q_dot_me.push_back(T_amb_me.at(i, 5));
            // High-level
            IndepVars.push_back(vector<double>(3, 0.));
            IndepVars.back().at(0) = T_htf_hot_high;
            IndepVars.back().at(1) = m_dot_ND_des;
            IndepVars.back().at(2) = T_amb_me.at(i, 0);
            W_dot_me.push_back(T_amb_me.at(i, 3));
            Q_dot_me.push_back(T_amb_me.at(i, 6));
        }

        return 0;
    }

    // The following are helper functions for Rankine model tables
    double P_regr(double T_htf_ND, double T_amb_ND, double m_dot_ND, double P_ref, std::map<std::string, int> &r_map,
        const util::matrix_t<double> &T_htf_ptcs, const util::matrix_t<double> &T_amb_ptcs, const util::matrix_t<double> &m_dot_ptcs) {

        // Notation:
        // A = T_htf_ptcs
        // B = T_amb_ptcs
        // C = m_dot_ptcs

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
