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
//   VARTYPE        DATATYPE        NAME                        LABEL                                                       UNITS   META    GROUP   REQUIRED_IF     CONSTRAINTS     UI_HINTS*/
    { SSC_INPUT,    SSC_NUMBER,     "mut",                      "name of TCS compute module to be tested",                  "",     "",     "",     "*",            "",             "" },
    { SSC_INPUT,    SSC_NUMBER,     "indep_1_lo",               "lower bound of range for first var. when indep."           "",     "",     "",     "*",            "",             "" },
    { SSC_INPUT,    SSC_NUMBER,     "indep_1_hi",               "upper bound of range for first var. when indep.",          "",     "",     "",     "*",            "",             "" },
    { SSC_INPUT,    SSC_NUMBER,     "indep_1_n",                "num. of range divis. for first var. when indep.",          "",     "",     "",     "*",            "",             "" },
    { SSC_INPUT,    SSC_NUMBER,     "indep_2_lo",               "lower bound of range for secnd var. when indep.",          "",     "",     "",     "*",            "",             "" },
    { SSC_INPUT,    SSC_NUMBER,     "indep_2_hi",               "upper bound of range for secnd var. when indep.",          "",     "",     "",     "*",            "",             "" },
    { SSC_INPUT,    SSC_NUMBER,     "indep_2_n",                "num. of range divis. for secnd var. when indep.",          "",     "",     "",     "*",            "",             "" },
    { SSC_INPUT,    SSC_NUMBER,     "indep_3_lo",               "lower bound of range for third var. when indep.",          "",     "",     "",     "*",            "",             "" },
    { SSC_INPUT,    SSC_NUMBER,     "indep_3_hi",               "upper bound of range for third var. when indep.",          "",     "",     "",     "*",            "",             "" },
    { SSC_INPUT,    SSC_NUMBER,     "indep_3_n",                "num. of range divis. for third var. when indep.",          "",     "",     "",     "*",            "",             "" },

    { SSC_OUTPUT,   SSC_MATRIX,     "indep_1_w_3",              "Parametric of 1st indep. var. w/ 3rd indep. var. levels",  "",     "",     "",     "*",            "",             "" },
    { SSC_OUTPUT,   SSC_MATRIX,     "indep_2_w_1",              "Parametric of 2nd indep. var. w/ 1st indep. var. levels",  "",     "",     "",     "*",            "",             "" },
    { SSC_OUTPUT,   SSC_MATRIX,     "indep_3_w_2",              "Parametric of 3rd indep. var. w/ 2nd indep. var. levels",  "",     "",     "",     "*",            "",             "" },

    var_info_invalid };


class cm_validate_pc_tables : public compute_module
{
public:
    cm_validate_pc_tables()
    {
        add_var_info(_cm_vtab_validate_pc_tables);
    }

    C_sco2_rc_csp_template *mut;            // model under test, TODO - generalize the type of model
    int out_type = -1;
    std::string out_msg = "";

    void exec() throw(general_error)
    {
        util::matrix_t<double> indep_1_w_3, indep_2_w_1, indep_3_w_2;
        
        // Get SSC inputs
        double indep_1_lo = as_double("indep_1_lo");
        double indep_1_hi = as_double("indep_1_hi");
        double indep_1_n = as_double("indep_1_n");
        double indep_2_lo = as_double("indep_2_lo");
        double indep_2_hi = as_double("indep_2_hi");
        double indep_2_n = as_double("indep_2_n");
        double indep_3_lo = as_double("indep_3_lo");
        double indep_3_hi = as_double("indep_3_hi");
        double indep_3_n = as_double("indep_3_n");

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
};

DEFINE_MODULE_ENTRY(validate_pc_tables, "Regression model validation framework", 0)
