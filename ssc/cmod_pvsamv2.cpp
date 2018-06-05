/*******************************************************************************************************
*  Copyright 2017 Alliance for Sustainable Energy, LLC
*
*  NOTICE: This software was developed at least in part by Alliance for Sustainable Energy, LLC
*  (�Alliance�) under Contract No. DE-AC36-08GO28308 with the U.S. Department of Energy and the U.S.
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
*  the underlying software originally provided by Alliance as �System Advisor Model� or �SAM�. Except
*  to comply with the foregoing, the terms �System Advisor Model�, �SAM�, or any confusingly similar
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

#include "cmod_pvsamv2.h"

// comment following define if do not want shading database validation outputs
//#define SHADE_DB_OUTPUTS

static var_info _cm_vtab_pvsamv2[] = {
/*   VARTYPE           DATATYPE         NAME                                            LABEL                                                   UNITS      META                             GROUP                  REQUIRED_IF                 CONSTRAINTS                      UI_HINTS*/
	{ SSC_INPUT,        SSC_STRING,      "solar_resource_file",                         "Weather file in TMY2, TMY3, EPW, or SAM CSV.",         "",         "",                              "pvsamv1",              "?",                        "",                              "" },
	{ SSC_INPUT,        SSC_TABLE,       "solar_resource_data",                         "Weather data",                                         "",         "lat,lon,tz,elev,year,month,hour,minute,gh,dn,df,poa,tdry,twet,tdew,rhum,pres,snow,alb,aod,wspd,wdir",    "pvsamv1",              "?",                        "",                              "" },
	
	// transformer model percent of rated ac output
	{ SSC_INPUT, SSC_NUMBER, "transformer_no_load_loss", "Power transformer no load loss", "%", "", "pvsamv1", "?=0", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "transformer_load_loss", "Power transformer load loss", "%", "", "pvsamv1", "?=0", "", "" },
	
	// optional for lifetime analysis
	{ SSC_INPUT,        SSC_NUMBER,      "system_use_lifetime_output",                  "PV lifetime simulation",                               "0/1",      "",                              "pvsamv1",             "?=0",                        "INTEGER,MIN=0,MAX=1",          "" },
	{ SSC_INPUT,        SSC_NUMBER,      "analysis_period",                             "Lifetime analysis period",                             "years",    "",                              "pvsamv1",             "system_use_lifetime_output=1",   "",                             "" },
	{ SSC_INPUT,        SSC_ARRAY,       "dc_degradation",                              "Annual module degradation",                            "%/year",   "",                              "pvsamv1",             "system_use_lifetime_output=1",   "",                             "" },
//	{ SSC_INPUT,        SSC_ARRAY,       "ac_degradation",                              "Annual AC degradation",                                "%/year",   "",                              "pvsamv1",             "system_use_lifetime_output=1",   "",                             "" },
	{ SSC_OUTPUT,       SSC_ARRAY,       "dc_degrade_factor",                           "Annual module degrade factor",                         "",         "",                              "pvsamv1",             "system_use_lifetime_output=1",   "",                             "" },
//	{ SSC_OUTPUT,       SSC_ARRAY,       "ac_degrade_factor",                           "Annual AC degrade factor",                             "",         "",                              "pvsamv1",             "system_use_lifetime_output=1",   "",                             "" },
	{ SSC_INPUT,        SSC_NUMBER,      "en_dc_lifetime_losses",                       "Enable lifetime daily DC losses",                      "0/1",      "",                              "pvsamv1",             "?=0",                        "INTEGER,MIN=0,MAX=1",          "" },
	{ SSC_INPUT,        SSC_ARRAY,       "dc_lifetime_losses",                          "Lifetime daily DC losses",                             "%",        "",                              "pvsamv1",             "en_dc_lifetime_losses=1",    "",                             "" },
	{ SSC_INPUT,        SSC_NUMBER,      "en_ac_lifetime_losses",                       "Enable lifetime daily AC losses",                      "0/1",      "",                              "pvsamv1",             "?=0",                        "INTEGER,MIN=0,MAX=1",          "" },
	{ SSC_INPUT,        SSC_ARRAY,       "ac_lifetime_losses",                          "Lifetime daily AC losses",                             "%",        "",                              "pvsamv1",             "en_ac_lifetime_losses=1",    "",                             "" },
														                        																		                             
	//SEV: Activating the snow model							                        																		                             
	{ SSC_INPUT,        SSC_NUMBER,      "en_snow_model",                               "Toggle snow loss estimation",                          "0/1",      "",                              "snowmodel",            "?=0",                       "BOOLEAN",                      "" },
	{ SSC_INPUT,        SSC_NUMBER,      "system_capacity",                             "Nameplate capacity",                                   "kW",       "",                              "pvsamv1",              "*",                         "",                             "" },
	
	{ SSC_INPUT,        SSC_NUMBER,      "use_wf_albedo",                               "Use albedo in weather file if provided",               "0/1",      "",                              "pvsamv1",              "?=1",                      "BOOLEAN",                       "" },
	{ SSC_INPUT,        SSC_ARRAY,       "albedo",                                      "User specified ground albedo",                         "0..1",     "",                              "pvsamv1",              "*",						  "LENGTH=12",					  "" },
	{ SSC_INPUT,        SSC_NUMBER,      "irrad_mode",                                  "Irradiance input translation mode",                    "",         "0=beam&diffuse,1=total&beam,2=total&diffuse,3=poa_reference,4=poa_pyranometer",   "pvsamv1",              "?=0",      "INTEGER,MIN=0,MAX=4",           "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sky_model",                                   "Diffuse sky model",                                    "",         "0=isotropic,1=hkdr,2=perez",    "pvsamv1",              "?=2",                      "INTEGER,MIN=0,MAX=2",           "" },
	 
	{ SSC_INPUT,        SSC_NUMBER,      "modules_per_string",                          "Modules per string",                                    "",        "",                              "pvsamv1",              "*",                        "INTEGER,POSITIVE",              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "strings_in_parallel",                         "String in parallel",                                    "",        "",                              "pvsamv1",              "*",                        "INTEGER,POSITIVE",              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inverter_count",                              "Number of inverters",                                   "",        "",                              "pvsamv1",              "*",                        "INTEGER,POSITIVE",              "" },
	
	{ SSC_INPUT,        SSC_NUMBER,      "enable_mismatch_vmax_calc",                   "Enable mismatched subarray Vmax calculation",           "",        "",                              "pvsamv1",              "?=0",                      "BOOLEAN",                       "" },

	{ SSC_INPUT,        SSC_NUMBER,      "subarray1_tilt",                              "Sub-array 1 Tilt",                                      "deg",     "0=horizontal,90=vertical",      "pvsamv1",              "naof:subarray1_tilt_eq_lat", "MIN=0,MAX=90",                "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray1_tilt_eq_lat",                       "Sub-array 1 Tilt=latitude override",                    "0/1",     "",                              "pvsamv1",              "na:subarray1_tilt",          "BOOLEAN",                     "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray1_azimuth",                           "Sub-array 1 Azimuth",                                   "deg",     "0=N,90=E,180=S,270=W",          "pvsamv1",              "*",                        "MIN=0,MAX=359.9",               "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray1_track_mode",                        "Sub-array 1 Tracking mode",                             "",        "0=fixed,1=1axis,2=2axis,3=azi,4=monthly", "pvsamv1",    "*",                        "INTEGER,MIN=0,MAX=4",           "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray1_rotlim",                            "Sub-array 1 Tracker rotation limit",                    "deg",     "",                              "pvsamv1",              "?=45",                     "MIN=5,MAX=85",                  "" },
	{ SSC_INPUT,		SSC_NUMBER,		 "subarray1_shade_mode",				     	"Sub-array 1 shading mode (fixed tilt or 1x tracking)",	 "0/1/2",   "0=none,1=standard(non-linear),2=thin film(linear)",  "pvsamv1",			     "*",                        "INTEGER,MIN=0,MAX=2",		      "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray1_gcr",                               "Sub-array 1 Ground coverage ratio",                     "0..1",    "",                              "pvsamv1",              "?=0.3",                    "MIN=0,MAX=3",               "" },
	{ SSC_INPUT,        SSC_ARRAY,       "subarray1_monthly_tilt",                      "Sub-array 1 monthly tilt input",                        "deg",     "",                              "pvsamv1",              "subarray1_track_mode=4",   "LENGTH=12",                     "" },
//	{ SSC_INPUT, SSC_ARRAY, "subarray1_shading:hourly", "Sub-array 1 Hourly beam shading losses", "%", "", "pvsamv1", "?", "", "" },
//	{ SSC_INPUT, SSC_NUMBER, "subarray1_shading:shading_db_lookup", "Sub-array 1 enable shading database lookup", "", "", "pvsamv1", "?=0", "BOOLEAN", "" },
//	{ SSC_INPUT, SSC_NUMBER, "subarray1_shading:string_option", "Sub-array 1 shading string option", "", "0=shadingdb,1=average,2=maximum,3=minimum", "pvsamv1", "?=-1", "INTEGER,MIN=-1,MAX=3", "" },
	{ SSC_INPUT, SSC_NUMBER, "subarray1_shading:string_option", "Sub-array 1 shading string option", "", "0=shadingdb,1=shadingdb_notc,2=average,3=maximum,4=minimum", "pvsamv1", "?=-1", "INTEGER,MIN=-1,MAX=4", "" },
	{ SSC_INPUT, SSC_MATRIX, "subarray1_shading:timestep", "Sub-array 1 timestep beam shading losses", "%", "", "pvsamv1", "?", "", "" },
	{ SSC_INPUT, SSC_MATRIX, "subarray1_shading:mxh", "Sub-array 1 Month x Hour beam shading losses", "%", "", "pvsamv1", "?", "", "" },
	{ SSC_INPUT,        SSC_MATRIX,      "subarray1_shading:azal",                      "Sub-array 1 Azimuth x altitude beam shading losses",    "%",       "",                              "pvsamv1",              "?",                        "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray1_shading:diff",                      "Sub-array 1 Diffuse shading loss",                      "%",       "",                              "pvsamv1",              "?",                        "",                              "" },
	{ SSC_INPUT,        SSC_ARRAY,       "subarray1_soiling",                           "Sub-array 1 Monthly soiling loss",                      "%",       "",                              "pvsamv1",              "*",                        "LENGTH=12",                      "" },         

	// loss diagram outputs, also used to calculate total dc derate
	{ SSC_INPUT, SSC_NUMBER, "subarray1_mismatch_loss", "Sub-array 1 DC mismatch loss", "%", "", "pvsamv1", "*", "MIN=0,MAX=100", "" },
	{ SSC_INPUT, SSC_NUMBER, "subarray1_diodeconn_loss", "Sub-array 1 DC diodes and connections loss", "%", "", "pvsamv1", "*", "MIN=0,MAX=100", "" },
	{ SSC_INPUT, SSC_NUMBER, "subarray1_dcwiring_loss", "Sub-array 1 DC wiring loss", "%", "", "pvsamv1", "*", "MIN=0,MAX=100", "" },
	{ SSC_INPUT, SSC_NUMBER, "subarray1_tracking_loss", "Sub-array 1 DC tracking error loss", "%", "", "pvsamv1", "*", "MIN=0,MAX=100", "" },
	{ SSC_INPUT, SSC_NUMBER, "subarray1_nameplate_loss", "Sub-array 1 DC nameplate loss", "%", "", "pvsamv1", "*", "MIN=-5,MAX=100", "" },

	{ SSC_INPUT, SSC_NUMBER, "subarray2_mismatch_loss", "Sub-array 2 DC mismatch loss", "%", "", "pvsamv1", "?", "MIN=0,MAX=100", "" },
	{ SSC_INPUT, SSC_NUMBER, "subarray2_diodeconn_loss", "Sub-array 2 DC diodes and connections loss", "%", "", "pvsamv1", "?", "MIN=0,MAX=100", "" },
	{ SSC_INPUT, SSC_NUMBER, "subarray2_dcwiring_loss", "Sub-array 2 DC wiring loss", "%", "", "pvsamv1", "?", "MIN=0,MAX=100", "" },
	{ SSC_INPUT, SSC_NUMBER, "subarray2_tracking_loss", "Sub-array 2 DC tracking error loss", "%", "", "pvsamv1", "?", "MIN=0,MAX=100", "" },
	{ SSC_INPUT, SSC_NUMBER, "subarray2_nameplate_loss", "Sub-array 2 DC nameplate loss", "%", "", "pvsamv1", "?", "MIN=-5,MAX=100", "" },

	{ SSC_INPUT, SSC_NUMBER, "subarray3_mismatch_loss", "Sub-array 3 DC mismatch loss", "%", "", "pvsamv1", "?", "MIN=0,MAX=100", "" },
	{ SSC_INPUT, SSC_NUMBER, "subarray3_diodeconn_loss", "Sub-array 3 DC diodes and connections loss", "%", "", "pvsamv1", "?", "MIN=0,MAX=100", "" },
	{ SSC_INPUT, SSC_NUMBER, "subarray3_dcwiring_loss", "Sub-array 3 DC wiring loss", "%", "", "pvsamv1", "?", "MIN=0,MAX=100", "" },
	{ SSC_INPUT, SSC_NUMBER, "subarray3_tracking_loss", "Sub-array 3 DC tracking error loss", "%", "", "pvsamv1", "?", "MIN=0,MAX=100", "" },
	{ SSC_INPUT, SSC_NUMBER, "subarray3_nameplate_loss", "Sub-array 3 DC nameplate loss", "%", "", "pvsamv1", "?", "MIN=-5,MAX=100", "" },

	{ SSC_INPUT, SSC_NUMBER, "subarray4_mismatch_loss", "Sub-array 4 DC mismatch loss", "%", "", "pvsamv1", "?", "MIN=0,MAX=100", "" },
	{ SSC_INPUT, SSC_NUMBER, "subarray4_diodeconn_loss", "Sub-array 4 DC diodes and connections loss", "%", "?", "pvsamv1", "?", "MIN=0,MAX=100", "" },
	{ SSC_INPUT, SSC_NUMBER, "subarray4_dcwiring_loss", "Sub-array 4 DC wiring loss", "%", "", "pvsamv1", "?", "MIN=0,MAX=100", "" },
	{ SSC_INPUT, SSC_NUMBER, "subarray4_tracking_loss", "Sub-array 4 DC tracking error loss", "%", "", "pvsamv1", "?", "MIN=0,MAX=100", "" },
	{ SSC_INPUT, SSC_NUMBER, "subarray4_nameplate_loss", "Sub-array 4 DC nameplate loss", "%", "", "pvsamv1", "?", "MIN=-5,MAX=100", "" },

	//this is a DC loss that is applied uniformly to all subarrays
	{ SSC_INPUT, SSC_NUMBER, "dcoptimizer_loss", "DC power optimizer loss", "%", "", "pvsamv1", "*", "MIN=0,MAX=100", "" },
	//AC losses are also applied uniformly to all subarrays
	{ SSC_INPUT, SSC_NUMBER, "acwiring_loss", "AC wiring loss", "%", "", "pvsamv1", "*", "MIN=0,MAX=100", "" },
//	{ SSC_INPUT, SSC_NUMBER, "transformer_loss", "AC step-up transformer loss", "%", "", "pvsamv1", "*", "MIN=0,MAX=100", "" },


	//

	{ SSC_INPUT,        SSC_NUMBER,      "subarray1_mod_orient",                        "Sub-array 1 Module orientation for self-shading",         "0/1",    "0=portrait,1=landscape",        "pvsamv1",              "subarray1_shade_mode>0", "INTEGER,MIN=0,MAX=1",           "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray1_nmodx",                             "Sub-array 1 no. of modules along bottom for self-shading","",       "",                              "pvsamv1",              "subarray1_shade_mode>0", "INTEGER,POSITIVE",              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray1_nmody",                             "Sub-array 1 no. of modules along side for self-shading",  "",       "",                              "pvsamv1",              "subarray1_shade_mode>0", "INTEGER,POSITIVE",              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray1_backtrack",                         "Sub-array 1 Backtracking enabled",                        "",       "0=no backtracking,1=backtrack", "pvsamv1",              "subarray1_track_mode=1",   "BOOLEAN",                       "" },

	{ SSC_INPUT,        SSC_NUMBER,      "subarray2_enable",                            "Sub-array 2 Enable",                                      "0/1",    "0=disabled,1=enabled",          "pvsamv1",              "?=0",                      "BOOLEAN",                       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray2_nstrings",                          "Sub-array 2 Number of parallel strings",                  "",       "",                              "pvsamv1",              "subarray2_enable=1",       "INTEGER",                       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray2_tilt",                              "Sub-array 2 Tilt",                                        "deg",    "0=horizontal,90=vertical",      "pvsamv1",              "naof:subarray2_tilt_eq_lat", "MIN=0,MAX=90",                "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray2_tilt_eq_lat",                       "Sub-array 2 Tilt=latitude override",                      "0/1",    "",                              "pvsamv1",              "na:subarray2_tilt",          "BOOLEAN",                     "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray2_azimuth",                           "Sub-array 2 Azimuth",                                     "deg",    "0=N,90=E,180=S,270=W",          "pvsamv1",              "subarray2_enable=1",       "MIN=0,MAX=359.9",               "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray2_track_mode",                        "Sub-array 2 Tracking mode",                               "",       "0=fixed,1=1axis,2=2axis,3=azi,4=monthly", "pvsamv1",    "subarray2_enable=1",       "INTEGER,MIN=0,MAX=4",           "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray2_rotlim",                            "Sub-array 2 Tracker rotation limit",                      "deg",    "",                              "pvsamv1",              "?=45",                     "MIN=5,MAX=85",                  "" },
	{ SSC_INPUT,		SSC_NUMBER,		 "subarray2_shade_mode",				     	"Sub-array 2 shading mode (fixed tilt or 1x tracking)",	   "0/1/2",   "0=none,1=standard(non-linear),2=thin film(linear)",  "pvsamv1",		      "*",                        "INTEGER,MIN=0,MAX=2",		   "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray2_gcr",                               "Sub-array 2 Ground coverage ratio",                       "0..1",   "",                              "pvsamv1",              "?=0.3",                    "MIN=0,MAX=3",               "" },
	{ SSC_INPUT,        SSC_ARRAY,       "subarray2_monthly_tilt",                      "Sub-array 2 monthly tilt input",                          "deg",    "",                              "pvsamv1",              "subarray2_track_mode=4",   "LENGTH=12",                     "" },
//	{ SSC_INPUT,        SSC_ARRAY,       "subarray2_shading:hourly",                    "Sub-array 2 Hourly beam shading losses",                 "%",       "",                              "pvsamv1",              "?",                        "",                              "" },
//	{ SSC_INPUT, SSC_NUMBER, "subarray2_shading:shading_db_lookup", "Sub-array 2 enable shading database lookup", "", "", "pvsamv1", "?=0", "BOOLEAN", "" },
//	{ SSC_INPUT, SSC_NUMBER, "subarray2_shading:string_option", "Sub-array 2 shading string option", "", "0=shadingdb,1=average,2=maximum,3=minimum", "pvsamv1", "?=-1", "INTEGER,MIN=-1,MAX=3", "" },
	{ SSC_INPUT, SSC_NUMBER, "subarray2_shading:string_option", "Sub-array 2 shading string option", "", "0=shadingdb,1=shadingdb_notc,2=average,3=maximum,4=minimum", "pvsamv1", "?=-1", "INTEGER,MIN=-1,MAX=4", "" },
	{ SSC_INPUT, SSC_MATRIX, "subarray2_shading:timestep", "Sub-array 2 timestep beam shading losses", "%", "", "pvsamv1", "?", "", "" },
	{ SSC_INPUT, SSC_MATRIX, "subarray2_shading:mxh", "Sub-array 2 Month x Hour beam shading losses", "%", "", "pvsamv1", "?", "", "" },
	{ SSC_INPUT,        SSC_MATRIX,      "subarray2_shading:azal",                      "Sub-array 2 Azimuth x altitude beam shading losses",     "%",       "",                              "pvsamv1",              "?",                        "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray2_shading:diff",                      "Sub-array 2 Diffuse shading loss",                       "%",       "",                              "pvsamv1",              "?",                        "",                              "" },
	{ SSC_INPUT,        SSC_ARRAY,       "subarray2_soiling",                           "Sub-array 2 Monthly soiling loss",                       "%",   "",                              "pvsamv1",              "subarray2_enable=1",       "LENGTH=12",                     "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray2_mod_orient",                        "Sub-array 2 Module orientation for self-shading",         "0/1",    "0=portrait,1=landscape",        "pvsamv1",              "subarray2_shade_mode>0",  "INTEGER,MIN=0,MAX=1",           "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray2_nmodx",                             "Sub-array 2 no. of modules along bottom for self-shading","",       "",                              "pvsamv1",              "subarray2_shade_mode>0",  "INTEGER,POSITIVE",              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray2_nmody",                             "Sub-array 2 no. of modules along side for self-shading",  "",       "",                              "pvsamv1",              "subarray2_shade_mode>0",  "INTEGER,POSITIVE",              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray2_backtrack",                         "Sub-array 2 Backtracking enabled",                        "",       "0=no backtracking,1=backtrack", "pvsamv1",              "subarray2_track_mode=1",   "BOOLEAN",                       "" },

	{ SSC_INPUT,        SSC_NUMBER,      "subarray3_enable",                            "Sub-array 3 Enable",                                      "0/1",    "0=disabled,1=enabled",          "pvsamv1",              "?=0",                      "BOOLEAN",                       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray3_nstrings",                          "Sub-array 3 Number of parallel strings",                  "",       "",                              "pvsamv1",              "subarray3_enable=1",       "INTEGER",                       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray3_tilt",                              "Sub-array 3 Tilt",                                        "deg",    "0=horizontal,90=vertical",      "pvsamv1",              "naof:subarray3_tilt_eq_lat", "MIN=0,MAX=90",                "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray3_tilt_eq_lat",                       "Sub-array 3 Tilt=latitude override",                      "0/1",    "",                              "pvsamv1",              "na:subarray3_tilt",          "BOOLEAN",                     "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray3_azimuth",                           "Sub-array 3 Azimuth",                                     "deg",    "0=N,90=E,180=S,270=W",          "pvsamv1",              "subarray3_enable=1",       "MIN=0,MAX=359.9",               "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray3_track_mode",                        "Sub-array 3 Tracking mode",                               "",       "0=fixed,1=1axis,2=2axis,3=azi,4=monthly", "pvsamv1",    "subarray3_enable=1",       "INTEGER,MIN=0,MAX=4",           "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray3_rotlim",                            "Sub-array 3 Tracker rotation limit",                      "deg",    "",                              "pvsamv1",              "?=45",                     "MIN=5,MAX=85",                  "" },
	{ SSC_INPUT,		SSC_NUMBER,		 "subarray3_shade_mode",				     	"Sub-array 3 shading mode (fixed tilt or 1x tracking)",	   "0/1/2",   "0=none,1=standard(non-linear),2=thin film(linear)", "pvsamv1",			  "*",                        "INTEGER,MIN=0,MAX=2",		   "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray3_gcr",                               "Sub-array 3 Ground coverage ratio",                       "0..1",   "",                              "pvsamv1",              "?=0.3",                    "MIN=0,MAX=3",               "" },
	{ SSC_INPUT,        SSC_ARRAY,       "subarray3_monthly_tilt",                      "Sub-array 3 monthly tilt input",                          "deg",    "",                              "pvsamv1",              "subarray3_track_mode=4",   "LENGTH=12",                     "" },
//	{ SSC_INPUT,        SSC_ARRAY,       "subarray3_shading:hourly",                    "Sub-array 3 Hourly beam shading losses",                 "%",       "",                              "pvsamv1",              "?",                        "",                              "" },
//	{ SSC_INPUT, SSC_NUMBER, "subarray3_shading:shading_db_lookup", "Sub-array 3 enable shading database lookup", "", "", "pvsamv1", "?=0", "BOOLEAN", "" },
//	{ SSC_INPUT, SSC_NUMBER, "subarray3_shading:string_option", "Sub-array 3 shading string option", "", "0=shadingdb,1=average,2=maximum,3=minimum", "pvsamv1", "?=-1", "INTEGER,MIN=-1,MAX=3", "" },
	{ SSC_INPUT, SSC_NUMBER, "subarray3_shading:string_option", "Sub-array 3 shading string option", "", "0=shadingdb,1=shadingdb_notc,2=average,3=maximum,4=minimum", "pvsamv1", "?=-1", "INTEGER,MIN=-1,MAX=4", "" },
	{ SSC_INPUT, SSC_MATRIX, "subarray3_shading:timestep", "Sub-array 3 timestep beam shading losses", "%", "", "pvsamv1", "?", "", "" },
	{ SSC_INPUT, SSC_MATRIX, "subarray3_shading:mxh", "Sub-array 3 Month x Hour beam shading losses", "%", "", "pvsamv1", "?", "", "" },
	{ SSC_INPUT,        SSC_MATRIX,      "subarray3_shading:azal",                      "Sub-array 3 Azimuth x altitude beam shading losses",     "%",       "",                              "pvsamv1",              "?",                        "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray3_shading:diff",                      "Sub-array 3 Diffuse shading loss",                       "%",       "",                              "pvsamv1",              "?",                        "",                              "" },
	{ SSC_INPUT,        SSC_ARRAY,       "subarray3_soiling",                           "Sub-array 3 Monthly soiling loss",                       "%",   "",                              "pvsamv1",              "subarray3_enable=1",       "LENGTH=12",                     "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray3_mod_orient",                        "Sub-array 3 Module orientation for self-shading",         "0/1",    "0=portrait,1=landscape",        "pvsamv1",              "subarray1_shade_mode>0", "INTEGER,MIN=0,MAX=1",           "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray3_nmodx",                             "Sub-array 3 no. of modules along bottom for self-shading","",       "",                              "pvsamv1",              "subarray3_shade_mode>0", "INTEGER,POSITIVE",              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray3_nmody",                             "Sub-array 3 no. of modules along side for self-shading",  "",       "",                              "pvsamv1",              "subarray3_shade_mode>0", "INTEGER,POSITIVE",              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray3_backtrack",                         "Sub-array 3 Backtracking enabled",                        "",       "0=no backtracking,1=backtrack", "pvsamv1",              "subarray3_track_mode=1",   "BOOLEAN",                       "" },

	{ SSC_INPUT,        SSC_NUMBER,      "subarray4_enable",                            "Sub-array 4 Enable",                                      "0/1",    "0=disabled,1=enabled",          "pvsamv1",              "?=0",                      "BOOLEAN",                       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray4_nstrings",                          "Sub-array 4 Number of parallel strings",                  "",       "",                              "pvsamv1",              "subarray4_enable=1",       "INTEGER",                       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray4_tilt",                              "Sub-array 4 Tilt",                                        "deg",    "0=horizontal,90=vertical",      "pvsamv1",              "naof:subarray4_tilt_eq_lat", "MIN=0,MAX=90",                "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray4_tilt_eq_lat",                       "Sub-array 4 Tilt=latitude override",                      "0/1",    "",                              "pvsamv1",              "na:subarray4_tilt",          "BOOLEAN",                     "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray4_azimuth",                           "Sub-array 4 Azimuth",                                     "deg",    "0=N,90=E,180=S,270=W",          "pvsamv1",              "subarray4_enable=1",       "MIN=0,MAX=359.9",               "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray4_track_mode",                        "Sub-array 4 Tracking mode",                               "",       "0=fixed,1=1axis,2=2axis,3=azi,4=monthly", "pvsamv1",    "subarray4_enable=1",       "INTEGER,MIN=0,MAX=4",           "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray4_rotlim",                            "Sub-array 4 Tracker rotation limit",                      "deg",    "",                              "pvsamv1",              "?=45",                     "MIN=5,MAX=85",                  "" },
	{ SSC_INPUT,		SSC_NUMBER,		 "subarray4_shade_mode",				     	"Sub-array 4 shading mode (fixed tilt or 1x tracking)",	   "0/1/2",  "0=none,1=standard(non-linear),2=thin film(linear)",  "pvsamv1",			  "*",                        "INTEGER,MIN=0,MAX=2",		   "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray4_gcr",                               "Sub-array 4 Ground coverage ratio",                       "0..1",   "",                              "pvsamv1",              "?=0.3",                    "MIN=0,MAX=3",               "" },
	{ SSC_INPUT,        SSC_ARRAY,       "subarray4_monthly_tilt",                      "Sub-array 4 monthly tilt input",                          "deg",    "",                              "pvsamv1",              "subarray4_track_mode=4",   "LENGTH=12",                     "" },
//	{ SSC_INPUT,        SSC_ARRAY,       "subarray4_shading:hourly",                    "Sub-array 4 Hourly beam shading losses",                 "%",       "",                              "pvsamv1",              "?",                        "",                              "" },
//	{ SSC_INPUT, SSC_NUMBER, "subarray4_shading:shading_db_lookup", "Sub-array 4 enable shading database lookup", "", "", "pvsamv1", "?=0", "BOOLEAN", "" },
//	{ SSC_INPUT, SSC_NUMBER, "subarray4_shading:string_option", "Sub-array 4 shading string option", "", "0=shadingdb,1=average,2=maximum,3=minimum", "pvsamv1", "?=-1", "INTEGER,MIN=-1,MAX=3", "" },
	{ SSC_INPUT, SSC_NUMBER, "subarray4_shading:string_option", "Sub-array 4 shading string option", "", "0=shadingdb,1=shadingdb_notc,2=average,3=maximum,4=minimum", "pvsamv1", "?=-1", "INTEGER,MIN=-1,MAX=4", "" },
	{ SSC_INPUT, SSC_MATRIX, "subarray4_shading:timestep", "Sub-array 4 timestep beam shading losses", "%", "", "pvsamv1", "?", "", "" },
	{ SSC_INPUT, SSC_MATRIX, "subarray4_shading:mxh", "Sub-array 4 Month x Hour beam shading losses", "%", "", "pvsamv1", "?", "", "" },
	{ SSC_INPUT,        SSC_MATRIX,      "subarray4_shading:azal",                      "Sub-array 4 Azimuth x altitude beam shading losses",     "%",       "",                              "pvsamv1",              "?",                        "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray4_shading:diff",                      "Sub-array 4 Diffuse shading loss",                       "%",       "",                              "pvsamv1",              "?",                        "",                              "" },
	{ SSC_INPUT,        SSC_ARRAY,       "subarray4_soiling",                           "Sub-array 4 Monthly soiling loss",                       "%",   "",                              "pvsamv1",              "subarray4_enable=1",       "LENGTH=12",                     "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray4_mod_orient",                        "Sub-array 4 Module orientation for self-shading",         "0/1",    "0=portrait,1=landscape",        "pvsamv1",              "subarray4_shade_mode>0", "INTEGER,MIN=0,MAX=1",           "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray4_nmodx",                             "Sub-array 4 no. of modules along bottom for self-shading","",       "",                              "pvsamv1",              "subarray4_shade_mode>0", "INTEGER,POSITIVE",              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray4_nmody",                             "Sub-array 4 no. of modules along side for self-shading",  "",       "",                              "pvsamv1",              "subarray4_shade_mode>0", "INTEGER,POSITIVE",              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "subarray4_backtrack",                         "Sub-array 4 Backtracking enabled",                        "",       "0=no backtracking,1=backtrack", "pvsamv1",              "subarray4_track_mode=1",   "BOOLEAN",                       "" },

	{ SSC_INPUT,        SSC_NUMBER,      "module_model",                                "Photovoltaic module model specifier",                     "",       "0=spe,1=cec,2=6par_user,3=snl,4=sd11-iec61853", "pvsamv1",              "*",                        "INTEGER,MIN=0,MAX=4",           "" },
	{ SSC_INPUT,        SSC_NUMBER,      "module_aspect_ratio",                         "Module aspect ratio",                                     "",       "",                              "pvsamv1",              "?=1.7",                    "",                              "POSITIVE" },
	{ SSC_INPUT,        SSC_NUMBER,      "spe_area",                                    "Module area",                                             "m2",     "",                              "pvsamv1",              "module_model=0",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "spe_rad0",                                    "Irradiance level 0",                                      "W/m2",   "",                              "pvsamv1",              "module_model=0",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "spe_rad1",                                    "Irradiance level 1",                                      "W/m2",   "",                              "pvsamv1",              "module_model=0",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "spe_rad2",                                    "Irradiance level 2",                                      "W/m2",   "",                              "pvsamv1",              "module_model=0",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "spe_rad3",                                    "Irradiance level 3",                                      "W/m2",   "",                              "pvsamv1",              "module_model=0",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "spe_rad4",                                    "Irradiance level 4",                                      "W/m2",   "",                              "pvsamv1",              "module_model=0",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "spe_eff0",                                    "Efficiency at irradiance level 0",                        "%",      "",                              "pvsamv1",              "module_model=0",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "spe_eff1",                                    "Efficiency at irradiance level 1",                        "%",      "",                              "pvsamv1",              "module_model=0",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "spe_eff2",                                    "Efficiency at irradiance level 2",                        "%",      "",                              "pvsamv1",              "module_model=0",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "spe_eff3",                                    "Efficiency at irradiance level 3",                        "%",      "",                              "pvsamv1",              "module_model=0",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "spe_eff4",                                    "Efficiency at irradiance level 4",                        "%",      "",                              "pvsamv1",              "module_model=0",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "spe_reference",                               "Reference irradiance level",                              "",       "",                              "pvsamv1",              "module_model=0",           "INTEGER,MIN=0,MAX=4",           "" },
	{ SSC_INPUT,        SSC_NUMBER,      "spe_module_structure",                        "Mounting and module structure",                           "",       "0=glass/cell/polymer sheet - open rack,1=glass/cell/glass - open rack,2=polymer/thin film/steel - open rack,3=Insulated back, building-integrated PV,4=close roof mount,5=user-defined",                      "pvsamv1",       "module_model=0",                    "INTEGER,MIN=0,MAX=5",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "spe_a",                                       "Cell temp parameter a",                                   "",       "",                              "pvsamv1",              "module_model=0",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "spe_b",                                       "Cell temp parameter b",                                   "",       "",                              "pvsamv1",              "module_model=0",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "spe_dT",                                      "Cell temp parameter dT",                                  "",       "",                              "pvsamv1",              "module_model=0",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "spe_temp_coeff",                              "Temperature coefficient",                                 "%/C",    "",                              "pvsamv1",              "module_model=0",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "spe_fd",                                      "Diffuse fraction",                                        "0..1",   "",                              "pvsamv1",              "module_model=0",           "MIN=0,MAX=1",                   "" },
	{ SSC_INPUT,        SSC_NUMBER,      "spe_vmp",                                     "Nominal max power voltage",                               "V",      "",                              "pvsamv1",              "module_model=0",           "POSITIVE",                      "" },
	{ SSC_INPUT,        SSC_NUMBER,      "spe_voc",                                     "Nominal open circuit voltage",                            "V",      "",                              "pvsamv1",              "module_model=0",           "POSITIVE",                      "" },

	{ SSC_INPUT,        SSC_NUMBER,      "cec_area",                                    "Module area",                                             "m2",     "",                              "pvsamv1",              "module_model=1",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_a_ref",                                   "Nonideality factor a",                                    "",       "",                              "pvsamv1",              "module_model=1",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_adjust",                                  "Temperature coefficient adjustment",                      "%",      "",                              "pvsamv1",              "module_model=1",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_alpha_sc",                                "Short circuit current temperature coefficient",           "A/C",    "",                              "pvsamv1",              "module_model=1",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_beta_oc",                                 "Open circuit voltage temperature coefficient",            "V/C",    "",                              "pvsamv1",              "module_model=1",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_gamma_r",                                 "Maximum power point temperature coefficient",             "%/C",    "",                              "pvsamv1",              "module_model=1",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_i_l_ref",                                 "Light current",                                           "A",      "",                              "pvsamv1",              "module_model=1",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_i_mp_ref",                                "Maximum power point current",                             "A",      "",                              "pvsamv1",              "module_model=1",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_i_o_ref",                                 "Saturation current",                                      "A",      "",                              "pvsamv1",              "module_model=1",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_i_sc_ref",                                "Short circuit current",                                   "A",      "",                              "pvsamv1",              "module_model=1",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_n_s",                                     "Number of cells in series",                               "",       "",                              "pvsamv1",              "module_model=1",           "POSITIVE",                      "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_r_s",                                     "Series resistance",                                       "ohm",    "",                              "pvsamv1",              "module_model=1",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_r_sh_ref",                                "Shunt resistance",                                        "ohm",    "",                              "pvsamv1",              "module_model=1",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_t_noct",                                  "Nominal operating cell temperature",                      "C",      "",                              "pvsamv1",              "module_model=1",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_v_mp_ref",                                "Maximum power point voltage",                             "V",      "",                              "pvsamv1",              "module_model=1",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_v_oc_ref",                                "Open circuit voltage",                                    "V",      "",                              "pvsamv1",              "module_model=1",           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_temp_corr_mode",                          "Cell temperature model selection",                        "",       "0=noct,1=mc",                   "pvsamv1",              "module_model=1",           "INTEGER,MIN=0,MAX=1",           "" },
	
	{ SSC_INPUT,        SSC_NUMBER,      "cec_standoff",                                "Standoff mode",                                           "",       "0=bipv,1=>3.5in,2=2.5-3.5in,3=1.5-2.5in,4=0.5-1.5in,5=<0.5in,6=ground/rack",  "pvsamv1",       "module_model=1",                           "INTEGER,MIN=0,MAX=6",       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_height",                                  "Array mounting height",                                   "",       "0=one story,1=two story",                                           "pvsamv1",       "module_model=1",                           "INTEGER,MIN=0,MAX=1",       "" },

	{ SSC_INPUT,        SSC_NUMBER,      "cec_mounting_config",                         "Mounting configuration",                                  "",       "0=rack,1=flush,2=integrated,3=gap",                                 "pvsamv1",       "module_model=1&cec_temp_corr_mode=1",      "INTEGER,MIN=0,MAX=3",       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_heat_transfer",                           "Heat transfer dimensions",                                "",       "0=module,1=array",                                                  "pvsamv1",       "module_model=1&cec_temp_corr_mode=1",      "INTEGER,MIN=0,MAX=1",       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_mounting_orientation",                    "Mounting structure orientation",                          "",       "0=do not impede flow,1=vertical supports,2=horizontal supports",    "pvsamv1",       "module_model=1&cec_temp_corr_mode=1",      "INTEGER,MIN=0,MAX=2",       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_gap_spacing",                             "Gap spacing",                                             "m",      "",                                                                  "pvsamv1",       "module_model=1&cec_temp_corr_mode=1",      "",                          "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_module_width",                            "Module width",                                            "m",      "",                                                                  "pvsamv1",       "module_model=1&cec_temp_corr_mode=1",      "",                          "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_module_length",                           "Module height",                                           "m",      "",                                                                  "pvsamv1",       "module_model=1&cec_temp_corr_mode=1",      "",                          "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_array_rows",                              "Rows of modules in array",                                "",       "",                                                                  "pvsamv1",       "module_model=1&cec_temp_corr_mode=1",      "",                          "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_array_cols",                              "Columns of modules in array",                             "",       "",                                                                  "pvsamv1",       "module_model=1&cec_temp_corr_mode=1",      "",                          "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cec_backside_temp",                           "Module backside temperature",                             "C",      "",                                                                  "pvsamv1",       "module_model=1&cec_temp_corr_mode=1",      "POSITIVE",                  "" },
		
	{ SSC_INPUT,        SSC_NUMBER,      "6par_celltech",                               "Solar cell technology type",                              "",       "monoSi=0,multiSi=1,CdTe=2,CIS=3,CIGS=4,Amorphous=5",                "pvsamv1",       "module_model=2",                           "INTEGER,MIN=0,MAX=5",       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "6par_vmp",                                    "Maximum power point voltage",                             "V",      "",                                                                  "pvsamv1",       "module_model=2",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "6par_imp",                                    "Imp",                                                     "A",      "",                                                                  "pvsamv1",       "module_model=2",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "6par_voc",                                    "Voc",                                                     "V",      "",                                                                  "pvsamv1",       "module_model=2",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "6par_isc",                                    "Isc",                                                     "A",      "",                                                                  "pvsamv1",       "module_model=2",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "6par_bvoc",                                   "Short circuit current temperature coefficient",           "V/C",    "",                                                                  "pvsamv1",       "module_model=2",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "6par_aisc",                                   "Open circuit voltage temperature coefficient",            "A/C",    "",                                                                  "pvsamv1",       "module_model=2",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "6par_gpmp",                                   "Maximum power point temperature coefficient",             "%/C",    "",                                                                  "pvsamv1",       "module_model=2",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "6par_nser",                                   "Nseries",                                                 "",       "",                                                                  "pvsamv1",       "module_model=2",                           "INTEGER,POSITIVE",              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "6par_area",                                   "Module area",                                             "m2",     "",                                                                  "pvsamv1",       "module_model=2",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "6par_tnoct",                                  "Nominal operating cell temperature",                      "C",      "",                                                                  "pvsamv1",       "module_model=2",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "6par_standoff",                               "Standoff mode",                                           "",       "0=bipv,1=>3.5in,2=2.5-3.5in,3=1.5-2.5in,4=0.5-1.5in,5=<0.5in,6=ground/rack",  "pvsamv1",       "module_model=2",                           "INTEGER,MIN=0,MAX=6",           "" },
	{ SSC_INPUT,        SSC_NUMBER,      "6par_mounting",                               "Array mounting height",                                   "",       "0=one story,1=two story",                                           "pvsamv1",       "module_model=2",                           "INTEGER,MIN=0,MAX=1",           "" },
	
	{ SSC_INPUT,        SSC_NUMBER,      "snl_module_structure",                        "Module and mounting structure configuration",             "",       "0=Use Database Values,1=glass/cell/polymer sheet - open rack,2=glass/cell/glass - open rack,3=polymer/thin film/steel - open rack,4=Insulated back building-integrated PV,5=close roof mount,6=user-defined",                      "pvsamv1",       "module_model=3",                    "INTEGER,MIN=0,MAX=6",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_a",                                       "Temperature coefficient a",                               "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_b",                                       "Temperature coefficient b",                               "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_dtc",                                     "Temperature coefficient dT",                              "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_ref_a",                                   "User-specified a",                                        "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_ref_b",                                   "User-specified b",                                        "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_ref_dT",                                  "User-specified dT",                                       "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_fd",                                      "Diffuse fraction",                                        "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_a0",                                      "Air mass polynomial coeff 0",                             "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_a1",                                      "Air mass polynomial coeff 1",                             "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_a2",                                      "Air mass polynomial coeff 2",                             "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_a3",                                      "Air mass polynomial coeff 3",                             "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_a4",                                      "Air mass polynomial coeff 4",                             "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_aimp",                                    "Max power point current temperature coefficient",         "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_aisc",                                    "Short circuit current temperature coefficient",           "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_area",                                    "Module area",                                             "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_b0",                                      "Incidence angle modifier polynomial coeff 0",             "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_b1",                                      "Incidence angle modifier polynomial coeff 1",             "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_b2",                                      "Incidence angle modifier polynomial coeff 2",             "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_b3",                                      "Incidence angle modifier polynomial coeff 3",             "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_b4",                                      "Incidence angle modifier polynomial coeff 4",             "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_b5",                                      "Incidence angle modifier polynomial coeff 5",             "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_bvmpo",                                   "Max power point voltage temperature coefficient",         "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_bvoco",                                   "Open circuit voltage temperature coefficient",            "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_c0",                                      "C0",                                                      "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_c1",                                      "C1",                                                      "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_c2",                                      "C2",                                                      "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_c3",                                      "C3",                                                      "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_c4",                                      "C4",                                                      "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_c5",                                      "C5",                                                      "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_c6",                                      "C6",                                                      "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_c7",                                      "C7",                                                      "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_impo",                                    "Max power point current",                                 "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_isco",                                    "Short circuit current",                                   "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_ixo",                                     "Ix midpoint current",                                     "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_ixxo",                                    "Ixx midpoint current",                                    "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_mbvmp",                                   "Irradiance dependence of Vmp temperature coefficient",    "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_mbvoc",                                   "Irradiance dependence of Voc temperature coefficient",    "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_n",                                       "Diode factor",                                            "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_series_cells",                            "Number of cells in series",                               "",       "",                      "pvsamv1",       "module_model=3",                    "INTEGER",                       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_vmpo",                                    "Max power point voltage",                                 "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "snl_voco",                                    "Open circuit voltage",                                    "",       "",                      "pvsamv1",       "module_model=3",                    "",                              "" },

	//{ SSC_INPUT,        SSC_NUMBER,      "sd11par_type",                                "Cell technology type",                                    "",       "monoSi=0,multiSi=1,CdTe=2,CIS=3,CIGS=4,Amorphous=5",                "pvsamv1",       "module_model=4",                           "INTEGER,MIN=0,MAX=5",       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_nser",                                "Nseries",                                                 "",       "",                                                                  "pvsamv1",       "module_model=4",                           "INTEGER,POSITIVE",              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_area",                                "Module area",                                             "m2",     "",                                                                  "pvsamv1",       "module_model=4",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_AMa0",                                "Air mass modifier coeff 0",                               "",       "",                                                                  "pvsamv1",       "module_model=4",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_AMa1",                                "Air mass modifier coeff 1",                               "",       "",                                                                  "pvsamv1",       "module_model=4",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_AMa2",                                "Air mass modifier coeff 2",                               "",       "",                                                                  "pvsamv1",       "module_model=4",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_AMa3",                                "Air mass modifier coeff 3",                               "",       "",                                                                  "pvsamv1",       "module_model=4",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_AMa4",                                "Air mass modifier coeff 4",                               "",       "",                                                                  "pvsamv1",       "module_model=4",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_glass",                               "Cover glass type",                                        "",       "0=normal,1=AR glass",                                               "pvsamv1",       "module_model=4",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_tnoct",                               "Nominal operating cell temperature",                      "C",      "",                                                                  "pvsamv1",       "module_model=4",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_standoff",                            "Standoff mode",                                           "",       "0=bipv,1=>3.5in,2=2.5-3.5in,3=1.5-2.5in,4=0.5-1.5in,5=<0.5in,6=ground/rack",  "pvsamv1",       "module_model=4",                 "INTEGER,MIN=0,MAX=6",           "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_mounting",                            "Array mounting height",                                   "",       "0=one story,1=two story",                                           "pvsamv1",       "module_model=4",                           "INTEGER,MIN=0,MAX=1",           "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_Vmp0",                                "Vmp (STC)",                                               "V",      "",                                                                  "pvsamv1",       "module_model=4",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_Imp0",                                "Imp (STC)",                                               "A",      "",                                                                  "pvsamv1",       "module_model=4",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_Voc0",                                "Voc (STC)",                                               "V",      "",                                                                  "pvsamv1",       "module_model=4",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_Isc0",                                "Isc (STC)",                                               "A",      "",                                                                  "pvsamv1",       "module_model=4",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_alphaIsc",                            "Open circuit voltage temperature coefficient",            "A/C",    "",                                                                  "pvsamv1",       "module_model=4",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_n",                                   "Diode nonideality factor",                                "",       "",                                                                  "pvsamv1",       "module_model=4",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_Il",                                  "Light current",                                           "A",      "",                                                                  "pvsamv1",       "module_model=4",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_Io",                                  "Saturation current",                                      "A",      "",                                                                  "pvsamv1",       "module_model=4",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_Egref",                               "Bandgap voltage",                                         "eV",     "",                                                                  "pvsamv1",       "module_model=4",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_d1",                                  "Rs fit parameter 1",                                      "",       "",                                                                  "pvsamv1",       "module_model=4",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_d2",                                  "Rs fit parameter 2",                                      "",       "",                                                                  "pvsamv1",       "module_model=4",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_d3",                                  "Rs fit parameter 3",                                      "",       "",                                                                  "pvsamv1",       "module_model=4",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_c1",                                  "Rsh fit parameter 1",                                     "",       "",                                                                  "pvsamv1",       "module_model=4",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_c2",                                  "Rsh fit parameter 2",                                     "",       "",                                                                  "pvsamv1",       "module_model=4",                           "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "sd11par_c3",                                  "Rsh fit parameter 3",                                     "",       "",                                                                  "pvsamv1",       "module_model=4",                           "",                              "" },
	
// inverter model
	{ SSC_INPUT,        SSC_NUMBER,      "inverter_model",                              "Inverter model specifier",                                "",        "0=cec,1=datasheet,2=partload,3=coefficientgenerator",        "pvsamv1",               "*",                         "INTEGER,MIN=0,MAX=3",           "" },
	{ SSC_INPUT,        SSC_NUMBER,      "mppt_low_inverter",                           "Minimum inverter MPPT voltage window",                    "Vdc",     "",                     "pvsamv1",       "",                    "?=0",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "mppt_hi_inverter",                            "Maximum inverter MPPT voltage window",                    "Vdc",     "",                     "pvsamv1",       "",                    "?=0",                              "" },

	{ SSC_INPUT,        SSC_NUMBER,      "inv_snl_c0",                                  "Curvature between AC power and DC power at ref",          "1/W",     "",                     "pvsamv1",       "inverter_model=0",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_snl_c1",                                  "Coefficient of Pdco variation with DC input voltage",     "1/V",     "",                     "pvsamv1",       "inverter_model=0",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_snl_c2",                                  "Coefficient of Pso variation with DC input voltage",      "1/V",     "",                     "pvsamv1",       "inverter_model=0",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_snl_c3",                                  "Coefficient of Co variation with DC input voltage",       "1/V",     "",                     "pvsamv1",       "inverter_model=0",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_snl_paco",                                "AC maximum power rating",                                 "Wac",     "",                     "pvsamv1",       "inverter_model=0",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_snl_pdco",                                "DC input power at which AC power rating is achieved",     "Wdc",     "",                     "pvsamv1",       "inverter_model=0",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_snl_pnt",                                 "AC power consumed by inverter at night",                  "Wac",     "",                     "pvsamv1",       "inverter_model=0",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_snl_pso",                                 "DC power required to enable the inversion process",       "Wdc",     "",                     "pvsamv1",       "inverter_model=0",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_snl_vdco",                                "DC input voltage for the rated AC power rating",          "Vdc",     "",                     "pvsamv1",       "inverter_model=0",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_snl_vdcmax",                              "Maximum DC input operating voltage",                      "Vdc",     "",                     "pvsamv1",       "inverter_model=0",                    "",                              "" },

	{ SSC_INPUT, SSC_NUMBER, "inv_cec_cg_c0", "Curvature between AC power and DC power at ref", "1/W", "", "pvsamv1", "inverter_model=3", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_cec_cg_c1", "Coefficient of Pdco variation with DC input voltage", "1/V", "", "pvsamv1", "inverter_model=3", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_cec_cg_c2", "Coefficient of Pso variation with DC input voltage", "1/V", "", "pvsamv1", "inverter_model=3", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_cec_cg_c3", "Coefficient of Co variation with DC input voltage", "1/V", "", "pvsamv1", "inverter_model=3", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_cec_cg_paco", "AC maximum power rating", "Wac", "", "pvsamv1", "inverter_model=3", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_cec_cg_pdco", "DC input power at which AC power rating is achieved", "Wdc", "", "pvsamv1", "inverter_model=3", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_cec_cg_pnt", "AC power consumed by inverter at night", "Wac", "", "pvsamv1", "inverter_model=3", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_cec_cg_psco", "DC power required to enable the inversion process", "Wdc", "", "pvsamv1", "inverter_model=3", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_cec_cg_vdco", "DC input voltage for the rated AC power rating", "Vdc", "", "pvsamv1", "inverter_model=3", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_cec_cg_vdcmax", "Maximum DC input operating voltage", "Vdc", "", "pvsamv1", "inverter_model=3", "", "" },

	{ SSC_INPUT,        SSC_NUMBER,      "inv_ds_paco",                                "AC maximum power rating",                                 "Wac",     "",                     "pvsamv1",       "inverter_model=1",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_ds_eff",                                 "Weighted or Peak or Nominal Efficiency",     "Wdc",     "",                     "pvsamv1",       "inverter_model=1",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_ds_pnt",                                 "AC power consumed by inverter at night",                  "Wac",     "",                     "pvsamv1",       "inverter_model=1",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_ds_pso",                                 "DC power required to enable the inversion process",       "Wdc",     "",                     "pvsamv1",       "inverter_model=1",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_ds_vdco",                                "DC input voltage for the rated AC power rating",          "Vdc",     "",                     "pvsamv1",       "inverter_model=1",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_ds_vdcmax",                              "Maximum DC input operating voltage",                      "Vdc",     "",                     "pvsamv1",       "inverter_model=1",                    "",                              "" },

	{ SSC_INPUT,        SSC_NUMBER,      "inv_pd_paco",                                "AC maximum power rating",                                 "Wac",     "",                     "pvsamv1",       "inverter_model=2",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_pd_pdco",                                "DC input power at which AC power rating is achieved",     "Wdc",     "",                     "pvsamv1",       "inverter_model=2",                    "",                              "" },
	{ SSC_INPUT,        SSC_ARRAY,       "inv_pd_partload",                            "Partload curve partload values",                          "%",       "",                     "pvsamv1",       "inverter_model=2",                    "",                              "" },
	{ SSC_INPUT,        SSC_ARRAY,       "inv_pd_efficiency",                          "Partload curve efficiency values",                        "%",       "",                     "pvsamv1",       "inverter_model=2",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_pd_pnt",                                 "AC power consumed by inverter at night",                  "Wac",     "",                     "pvsamv1",       "inverter_model=2",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_pd_vdco",                                "DC input voltage for the rated AC power rating",          "Vdc",     "",                     "pvsamv1",       "inverter_model=2",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_pd_vdcmax",                              "Maximum DC input operating voltage",                      "Vdc",     "",                     "pvsamv1",       "inverter_model=2",                    "",                              "" },
	
	// battery storage and dispatch
	{ SSC_INPUT,        SSC_NUMBER,      "en_batt",                                    "Enable battery storage model",                            "0/1",     "",                     "Battery",       "?=0",                                 "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "batt_replacement_option",                    "Enable battery replacement?",                             "0=none,1=capacity based,2=user schedule", "", "Battery", "?=0", "INTEGER,MIN=0,MAX=2", "" },
	{ SSC_INPUT,        SSC_ARRAY,       "batt_replacement_schedule",                  "Battery bank replacements per year (user specified)",     "number/year", "", "Battery", "batt_replacement_option=2", "", "" },

	{ SSC_INPUT,        SSC_ARRAY,       "load",                                       "Electricity load (year 1)",                         "kW", "", "Battery", "?", "", "" },
	
	// NOTE:  other battery storage model inputs and outputs are defined in batt_common.h/batt_common.cpp
	
	// outputs

/* environmental conditions */
	// irradiance data from weather file
	{ SSC_OUTPUT,        SSC_ARRAY,      "gh",                                         "Irradiance GHI from weather file",                                     "W/m2",   "",                      "Time Series",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "dn",                                         "Irradiance DNI from weather file",                                     "W/m2",   "",                      "Time Series",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "df",                                         "Irradiance DHI from weather file",                                     "W/m2",   "",                      "Time Series",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "wfpoa",                                      "Irradiance POA from weather file",                                     "W/m2",   "",                      "Time Series",       "",                     "",                              "" },
	
	//not all of these three calculated values will be reported, based on irrad_mode selection
	{ SSC_OUTPUT,        SSC_ARRAY,      "gh_calc",                                    "Irradiance GHI calculated",                                       "W/m2",   "",                      "Time Series",       "",                     "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "dn_calc",                                    "Irradiance DNI calculated",                                       "W/m2",   "",                      "Time Series",       "",                     "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "df_calc",                                    "Irradiance DHI calculated",                                       "W/m2",   "",                      "Time Series",       "",                     "",                              "" },

	// non-irradiance data from weather file
	{ SSC_OUTPUT,        SSC_ARRAY,      "wspd",                                       "Weather file wind speed",                                                        "m/s",    "",                      "Time Series",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "tdry",                                       "Weather file ambient temperature",                                               "C",      "",                      "Time Series",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "alb",                                        "Weather file albedo",							                                 "",       "",                     "Time Series",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "snowdepth",                                  "Weather file snow depth",							                            "cm",       "",                    "Time Series",       "",                    "",                              "" },

	// calculated sun position data
	{ SSC_OUTPUT,        SSC_ARRAY,      "sol_zen",                                    "Sun zenith angle",                                                  "deg",    "",                      "Time Series",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "sol_alt",                                    "Sun altitude angle",                                                "deg",    "",                      "Time Series",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "sol_azi",                                    "Sun azimuth angle",                                                 "deg",    "",                      "Time Series",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "sunup",                                      "Sun up over horizon",                                               "0/1/2/3", "",                     "Time Series",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "sunpos_hour",                                "Sun position time",                                     "hour",   "",                      "Time Series",       "",                     "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "airmass",                                    "Absolute air mass",                                                 "",       "",                      "Time Series",       "*",                    "",                              "" },

	/* sub-array level outputs */
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_surf_tilt",                  "Subarray 1 Surface tilt",                                              "deg",    "", "Time Series (Subarray 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_surf_azi",                   "Subarray 1 Surface azimuth",                                           "deg",    "", "Time Series (Subarray 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_aoi",                        "Subarray 1 Angle of incidence",                                        "deg",    "", "Time Series (Subarray 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_axisrot",                    "Subarray 1 Axis rotation for 1 axis trackers",                         "deg",    "", "Time Series (Subarray 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_idealrot",                   "Subarray 1 Axis rotation ideal for 1 axis trackers",                   "deg",    "", "Time Series (Subarray 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_poa_eff_beam",               "Subarray 1 POA beam irradiance after shading and soiling",             "W/m2",   "", "Time Series (Subarray 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_poa_eff_diff",               "Subarray 1 POA diffuse irradiance after shading and soiling",          "W/m2",   "", "Time Series (Subarray 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_poa_nom",                    "Subarray 1 POA total irradiance nominal",                              "W/m2",   "", "Time Series (Subarray 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_poa_shaded",                 "Subarray 1 POA total irradiance after shading only",                   "W/m2",   "", "Time Series (Subarray 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_poa_eff",                    "Subarray 1 POA total irradiance after shading and soiling",            "W/m2",   "", "Time Series (Subarray 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_soiling_derate",             "Subarray 1 Soiling beam irradiance factor",                            "frac",   "", "Time Series (Subarray 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_beam_shading_factor",        "Subarray 1 External shading and soiling beam irradiance factor",       "frac",   "", "Time Series (Subarray 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_linear_derate",              "Subarray 1 Self-shading linear beam irradiance factor",                "frac",   "", "Time Series (Subarray 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_ss_diffuse_derate",          "Subarray 1 Self-shading non-linear sky diffuse irradiance factor",     "frac",   "", "Time Series (Subarray 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_ss_reflected_derate",        "Subarray 1 Self-shading non-linear ground diffuse irradiance factor",  "frac",   "", "Time Series (Subarray 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_ss_derate",                  "Subarray 1 Self-shading non-linear DC factor",                         "frac",   "", "Time Series (Subarray 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "shadedb_subarray1_shade_frac",         "Subarray 1 Partial external shading DC factor",                        "frac",   "", "Time Series (Subarray 1)",       "",                     "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_snow_coverage",              "Subarray 1 Snow cover",                                                "0..1",   "", "Time Series (Subarray 1)",       "",                     "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_snow_loss",                  "Subarray 1 Snow cover DC power loss",                                  "kW",     "", "Time Series (Subarray 1)",       "",                     "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_modeff",                     "Subarray 1 Module efficiency",                                         "%",      "", "Time Series (Subarray 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_celltemp",                   "Subarray 1 Cell temperature",                                          "C",      "", "Time Series (Subarray 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_dc_voltage",                 "Subarray 1 Operating voltage",                                         "V",      "", "Time Series (Subarray 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_voc",                        "Subarray 1 Open circuit voltage",                                      "V",      "", "Time Series (Subarray 1)",       "",                     "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray1_isc",                        "Subarray 1 Short circuit current",                                     "A",      "", "Time Series (Subarray 1)",       "",                     "",                              "" },


	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_surf_tilt",                  "Subarray 2 Surface tilt",                                              "deg",    "", "Time Series (Subarray 2)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_surf_azi",                   "Subarray 2 Surface azimuth",                                           "deg",    "", "Time Series (Subarray 2)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_aoi",                        "Subarray 2 Angle of incidence",                                        "deg",    "", "Time Series (Subarray 2)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_axisrot",                    "Subarray 2 Axis rotation for 1 axis trackers",                         "deg",    "", "Time Series (Subarray 2)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_idealrot",                   "Subarray 2 Axis rotation ideal for 1 axis trackers",                   "deg",    "", "Time Series (Subarray 2)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_poa_eff_beam",               "Subarray 2 POA beam irradiance after shading and soiling",             "W/m2",   "", "Time Series (Subarray 2)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_poa_eff_diff",               "Subarray 2 POA diffuse irradiance after shading and soiling",          "W/m2",   "", "Time Series (Subarray 2)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_poa_nom",                    "Subarray 2 POA total irradiance nominal",                              "W/m2",   "", "Time Series (Subarray 2)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_poa_shaded",                 "Subarray 2 POA total irradiance after shading only",                   "W/m2",   "", "Time Series (Subarray 2)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_poa_eff",                    "Subarray 2 POA total irradiance after shading and soiling",            "W/m2",   "", "Time Series (Subarray 2)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_soiling_derate",             "Subarray 2 Soiling beam irradiance factor",                            "frac",   "", "Time Series (Subarray 2)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_beam_shading_factor",        "Subarray 2 External shading and soiling beam irradiance factor",       "frac",   "", "Time Series (Subarray 2)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_linear_derate",              "Subarray 2 Self-shading linear beam irradiance factor",                "frac",   "", "Time Series (Subarray 2)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_ss_diffuse_derate",          "Subarray 2 Self-shading non-linear sky diffuse irradiance factor",     "frac",   "", "Time Series (Subarray 2)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_ss_reflected_derate",        "Subarray 2 Self-shading non-linear ground diffuse irradiance factor",  "frac",   "", "Time Series (Subarray 2)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_ss_derate",                  "Subarray 2 Self-shading non-linear DC factor",                         "frac",   "", "Time Series (Subarray 2)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "shadedb_subarray2_shade_frac",         "Subarray 2 Partial shading DC factor",                                 "frac",   "", "Time Series (Subarray 2)",       "",                     "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_snow_coverage",				 "Subarray 2 Snow cover",                                                "0..1",   "", "Time Series (Subarray 2)",       "",                     "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_snow_loss",					 "Subarray 2 Snow cover DC power loss",                                  "kW",     "", "Time Series (Subarray 2)",       "",                     "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_modeff",                     "Subarray 2 Module efficiency",                                         "%",      "", "Time Series (Subarray 2)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_celltemp",                   "Subarray 2 Cell temperature",                                          "C",      "", "Time Series (Subarray 2)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_dc_voltage",                 "Subarray 2 Operating voltage",                                         "V",      "", "Time Series (Subarray 2)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_voc",                        "Subarray 2 Open circuit voltage",                                      "V",      "", "Time Series (Subarray 2)",       "",                     "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray2_isc",                        "Subarray 2 Short circuit current",                                     "A",      "", "Time Series (Subarray 2)",       "",                     "",                              "" },

	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_surf_tilt",                  "Subarray 3 Surface tilt",                                              "deg",    "", "Time Series (Subarray 3)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_surf_azi",                   "Subarray 3 Surface azimuth",                                           "deg",    "", "Time Series (Subarray 3)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_aoi",                        "Subarray 3 Angle of incidence",                                        "deg",    "", "Time Series (Subarray 3)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_axisrot",                    "Subarray 3 Axis rotation for 1 axis trackers",                         "deg",    "", "Time Series (Subarray 3)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_idealrot",                   "Subarray 3 Axis rotation ideal for 1 axis trackers",                   "deg",    "", "Time Series (Subarray 3)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_poa_eff_beam",               "Subarray 3 POA beam irradiance after shading and soiling",             "W/m2",   "", "Time Series (Subarray 3)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_poa_eff_diff",               "Subarray 3 POA diffuse irradiance after shading and soiling",          "W/m2",   "", "Time Series (Subarray 3)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_poa_nom",                    "Subarray 3 POA total irradiance nominal",                              "W/m2",   "", "Time Series (Subarray 3)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_poa_shaded",                 "Subarray 3 POA total irradiance after shading only",                   "W/m2",   "", "Time Series (Subarray 3)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_poa_eff",                    "Subarray 3 POA total irradiance after shading and soiling",            "W/m2",   "", "Time Series (Subarray 3)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_soiling_derate",             "Subarray 3 Soiling beam irradiance factor",                            "frac",   "", "Time Series (Subarray 3)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_beam_shading_factor",        "Subarray 3 External shading and soiling beam irradiance factor",       "frac",   "", "Time Series (Subarray 3)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_linear_derate",              "Subarray 3 Self-shading linear beam irradiance factor",                "frac",   "", "Time Series (Subarray 3)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_ss_diffuse_derate",          "Subarray 3 Self-shading non-linear sky diffuse irradiance factor",     "frac",   "", "Time Series (Subarray 3)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_ss_reflected_derate",        "Subarray 3 Self-shading non-linear ground diffuse irradiance factor",  "frac",   "", "Time Series (Subarray 3)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_ss_derate",                  "Subarray 3 Self-shading non-linear DC factor",                         "frac",   "", "Time Series (Subarray 3)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "shadedb_subarray3_shade_frac",         "Subarray 3 Partial external shading DC factor",                        "frac",   "", "Time Series (Subarray 3)",       "",                     "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_snow_coverage",				 "Subarray 3 Snow cover",                                                "0..1",   "", "Time Series (Subarray 3)",       "",                     "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_snow_loss",					 "Subarray 3 Snow cover DC power loss",			                         "kW",     "", "Time Series (Subarray 3)",       "",                     "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_modeff",                     "Subarray 3 Module efficiency",                                         "%",      "", "Time Series (Subarray 3)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_celltemp",                   "Subarray 3 Cell temperature",                                          "C",      "", "Time Series (Subarray 3)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_dc_voltage",                 "Subarray 3 Operating voltage",                                         "V",      "", "Time Series (Subarray 3)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_voc",                        "Subarray 3 Open circuit voltage",                                      "V",      "", "Time Series (Subarray 3)",       "",                     "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray3_isc",                        "Subarray 3 Short circuit current",                                     "A",      "", "Time Series (Subarray 3)",       "",                     "",                              "" },

	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_surf_tilt",                  "Subarray 4 Surface tilt",                                              "deg",    "", "Time Series (Subarray 4)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_surf_azi",                   "Subarray 4 Surface azimuth",                                           "deg",    "", "Time Series (Subarray 4)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_aoi",                        "Subarray 4 Angle of incidence",                                        "deg",    "", "Time Series (Subarray 4)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_axisrot",                    "Subarray 4 Axis rotation for 1 axis trackers",                         "deg",    "", "Time Series (Subarray 4)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_idealrot",                   "Subarray 4 Axis rotation ideal for 1 axis trackers",                   "deg",    "", "Time Series (Subarray 4)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_poa_eff_beam",               "Subarray 4 POA beam irradiance after shading and soiling",             "W/m2",   "", "Time Series (Subarray 4)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_poa_eff_diff",               "Subarray 4 POA diffuse irradiance after shading and soiling",          "W/m2",   "", "Time Series (Subarray 4)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_poa_nom",                    "Subarray 4 POA total irradiance nominal",                              "W/m2",   "", "Time Series (Subarray 4)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_poa_shaded",                 "Subarray 4 POA total irradiance after shading only",                   "W/m2",   "", "Time Series (Subarray 4)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_poa_eff",                    "Subarray 4 POA total irradiance after shading and soiling",            "W/m2",   "", "Time Series (Subarray 4)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_soiling_derate",             "Subarray 4 Soiling beam irradiance factor",                            "frac",   "", "Time Series (Subarray 4)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_beam_shading_factor",        "Subarray 4 External shading and soiling beam irradiance factor",       "frac",   "", "Time Series (Subarray 4)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_linear_derate",              "Subarray 4 Self-shading linear beam irradiance factor",                "frac",   "", "Time Series (Subarray 4)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_ss_diffuse_derate",          "Subarray 4 Self-shading non-linear sky diffuse irradiance factor",     "frac",   "", "Time Series (Subarray 4)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_ss_reflected_derate",        "Subarray 4 Self-shading non-linear ground diffuse irradiance factor",  "frac",   "", "Time Series (Subarray 4)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_ss_derate",                  "Subarray 4 Self-shading non-linear DC factor",                         "frac",   "", "Time Series (Subarray 4)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "shadedb_subarray4_shade_frac",         "Subarray 4 Partial external shading DC factor",                        "frac",   "", "Time Series (Subarray 4)",       "",                     "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_snow_coverage",				 "Subarray 4 Snow cover",                                                "0..1",   "", "Time Series (Subarray 4)",       "",                     "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_snow_loss",					 "Subarray 4 Snow cover DC power loss",                                  "kW",     "", "Time Series (Subarray 4)",       "",                     "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_modeff",                     "Subarray 4 Module efficiency",                                         "%",      "", "Time Series (Subarray 4)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_celltemp",                   "Subarray 4 Cell temperature",                                          "C",      "", "Time Series (Subarray 4)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_dc_voltage",                 "Subarray 4 Operating voltage",                                         "V",      "", "Time Series (Subarray 4)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_voc",                        "Subarray 4 Open circuit voltage",                                      "V",      "", "Time Series (Subarray 4)",       "",                     "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "subarray4_isc",                        "Subarray 4 Short circuit current",                                     "A",      "", "Time Series (Subarray 4)",       "",                     "",                              "" },

/* aggregate array level outputs */
	{ SSC_OUTPUT,        SSC_ARRAY,      "poa_nom",                              "Array POA total radiation nominal",                    "kW",   "",  "Time Series (Array)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "poa_beam_nom",                         "Array POA beam radiation nominal",                     "kW",   "",  "Time Series (Array)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "poa_shaded",                           "Array POA total radiation after shading only",         "kW",   "",  "Time Series (Array)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "poa_eff",                              "Array POA total radiation after shading and soiling",  "kW",   "",  "Time Series (Array)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "poa_beam_eff",                         "Array POA beam radiation after shading and soiling",   "kW",   "",  "Time Series (Array)",       "*",                    "",                              "" },

	//SEV: total dc snow loss time series (not a required output) 
	{ SSC_OUTPUT,        SSC_ARRAY,      "dc_snow_loss",                         "Array DC power loss due to snow",						 "kW",   "",   "Time Series (Array)",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "dc_net",                               "Array DC power",                                       "kW",   "",   "Time Series (Array)",       "*",                    "",                              "" },
	
	//inverter outputs
	{ SSC_OUTPUT,        SSC_ARRAY,      "inverter_dc_voltage",                  "Inverter DC input voltage",                            "V",    "",  "Time Series (Inverter)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "inv_eff",                              "Inverter efficiency",                                  "%",    "",  "Time Series (Inverter)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "dc_invmppt_loss",                      "Inverter clipping loss DC MPPT voltage limits",         "kW",  "",  "Time Series (Inverter)",       "*",                    "",                              "" },
    { SSC_OUTPUT,        SSC_ARRAY,      "inv_cliploss",                         "Inverter clipping loss AC power limit",                "kW",   "",  "Time Series (Inverter)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "inv_psoloss",                          "Inverter power consumption loss",                      "kW",   "",  "Time Series (Inverter)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "inv_pntloss",                          "Inverter night time loss",                             "kW",   "",  "Time Series (Inverter)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "ac_wiring_loss",                       "AC wiring loss",                                       "kW",   "",   "Time Series (Inverter)",              "*",                        "",                   "" },

	// transformer model outputs
	{ SSC_OUTPUT,        SSC_ARRAY,      "xfmr_nll_ts",                          "Transformer no load loss",                              "kW", "",    "Time Series (Transformer)", "", "", "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "xfmr_ll_ts",                           "Transformer load loss",                                 "kW", "",    "Time Series (Transformer)", "", "", "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "xfmr_loss_ts",                         "Transformer total loss",                                "kW", "",    "Time Series (Transformer)", "", "", "" },

	//total losses- not part of loss diagram but now outputs instead of inputs JMF 11/25/15
	{ SSC_OUTPUT,        SSC_NUMBER,      "ac_loss",                             "AC wiring loss",                                       "%",   "",    "Annual (Year 1)",              "*",                        "",                   "" },

	// monthly and annual outputs

	{ SSC_OUTPUT, SSC_NUMBER, "annual_energy", "Annual energy", "kWh", "", "Annual (Year 1)", "*", "", "" },

	{ SSC_OUTPUT,        SSC_NUMBER,     "annual_dc_invmppt_loss",                      "Inverter clipping loss DC MPPT voltage limits",          "kWh/yr",    "",                      "Annual (Year 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_NUMBER,     "annual_inv_cliploss",                         "Inverter clipping loss AC power limit",                  "kWh/yr",    "",                      "Annual (Year 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_NUMBER,     "annual_inv_psoloss",                          "Inverter power consumption loss",                        "kWh/yr",    "",                      "Annual (Year 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_NUMBER,     "annual_inv_pntloss",                          "Inverter night time loss",                               "kWh/yr",    "",                      "Annual (Year 1)",       "*",                    "",                              "" },

	{ SSC_OUTPUT,        SSC_NUMBER,      "subarray1_dcloss",                    "Subarray 1 Total DC power loss",                                       "%",      "", "Annual (Year 1)",              "*",                        "",                   "" },
	{ SSC_OUTPUT,        SSC_NUMBER,      "subarray2_dcloss",                    "Subarray 2 Total DC power loss",                                       "%",      "", "Annual (Year 1)",              "",                        "",                   "" },
	{ SSC_OUTPUT,        SSC_NUMBER,      "subarray3_dcloss",                    "Subarray 3 Total DC power loss",                                       "%",      "", "Annual (Year 1)",              "",                        "",                   "" },
	{ SSC_OUTPUT,        SSC_NUMBER,      "subarray4_dcloss",                    "Subarray 4 Total DC power loss",                                       "%",      "", "Annual (Year 1)",              "",                        "",                   "" },

	{ SSC_OUTPUT,        SSC_NUMBER,     "xfmr_nll_year1",                              "Transformer no load loss",                               "kWh/yr", "", "Annual (Year 1)", "", "", "" },
	{ SSC_OUTPUT,        SSC_NUMBER,     "xfmr_ll_year1",                               "Transformer load loss",                                  "kWh/yr", "", "Annual (Year 1)", "", "", "" },
	{ SSC_OUTPUT,        SSC_NUMBER,     "xfmr_loss_year1",                             "Transformer total loss",                                 "kWh/yr", "", "Annual (Year 1)", "", "", "" },

	{ SSC_OUTPUT,        SSC_ARRAY,      "monthly_poa_nom",                             "POA irradiance total nominal",                          "kWh/mo",    "",                      "Monthly",       "*",                    "LENGTH=12",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "monthly_poa_beam_nom",                        "POA irradiance beam nominal",                           "kWh/mo",    "",                      "Monthly",       "*",                    "LENGTH=12",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "monthly_poa_eff",                             "POA irradiance total after shading and soiling",          "kWh/mo",    "",                      "Monthly",       "*",                    "LENGTH=12",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "monthly_poa_beam_eff",                        "POA irradiance beam after shading and soiling",           "kWh/mo",    "",                      "Monthly",       "*",                    "LENGTH=12",                              "" },

	{ SSC_OUTPUT,        SSC_ARRAY,      "monthly_dc",                                  "PV array DC energy",                                   "kWh/mo",    "",                      "Monthly",       "*",                    "LENGTH=12",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "monthly_energy",                              "System AC energy",                                     "kWh/mo",    "",                      "Monthly",       "*",                    "LENGTH=12",                              "" },

	{ SSC_OUTPUT,        SSC_NUMBER,     "annual_gh",                                   "Annual GHI",                                              "Wh/m2/yr", "",                      "Annual (Year 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_NUMBER,     "annual_poa_nom",                              "POA irradiance total nominal",                          "kWh/yr",    "",                      "Annual (Year 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_NUMBER,     "annual_poa_beam_nom",                         "POA irradiance beam nominal",                           "kWh/yr",    "",                      "Annual (Year 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_NUMBER,     "annual_poa_shaded",                           "POA irradiancetotal after shading only",                 "kWh/yr",    "",                      "Annual (Year 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_NUMBER,     "annual_poa_eff",                              "POA irradiancetotal after shading and soiling",          "kWh/yr",    "",                      "Annual (Year 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_NUMBER,     "annual_poa_beam_eff",                         "POA irradiancebeam after shading and soiling",           "kWh/yr",    "",                      "Annual (Year 1)",       "*",                    "",                              "" },

	{ SSC_OUTPUT,        SSC_NUMBER,     "annual_dc_nominal",                           "Annual DC energy nominal",                           "kWh/yr",    "",                      "Annual (Year 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_NUMBER,     "annual_dc_gross",                             "Annual DC energy gross",                             "kWh/yr",    "",                      "Annual (Year 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_NUMBER,     "annual_dc_net",                               "Annual DC energy",                                   "kWh/yr",    "",                      "Annual (Year 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_NUMBER,     "annual_ac_gross",                             "Annual AC energy gross",                               "kWh/yr",    "",                      "Annual (Year 1)",       "*",                    "",                              "" },


	//SEV: total dc snow loss monthy array and annual value (not a required output) 
	{ SSC_OUTPUT,        SSC_ARRAY,      "monthly_snow_loss",                    "Snow DC energy loss",					       "kWh/mo",    "",                       "Monthly",       "",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_NUMBER,     "annual_snow_loss",                     "Snow DC energy loss",						   "kWh/yr",    "",                       "Annual (Year 1)",       "",                    "",                              "" },

	// loss diagram - order applied
	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray1_dc_gross", "Subarray 1 gross DC energy", "kWh", "", "Annual (Year 1)", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray1_dc_mismatch_loss", "Subarray 1 DC mismatch loss", "kWh", "", "Annual (Year 1)", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray1_dc_diodes_loss", "Subarray 1 DC diodes and connections loss", "kWh", "", "Annual (Year 1)", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray1_dc_wiring_loss", "Subarray 1 DC wiring loss", "kWh", "", "Annual (Year 1)", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray1_dc_tracking_loss", "Subarray 1 DC tracking loss", "kWh", "", "Annual (Year 1)", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray1_dc_nameplate_loss", "Subarray 1 DC nameplate loss", "kWh", "", "Annual (Year 1)", "*", "", "" },

	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray2_dc_gross", "Subarray 2 gross DC energy", "kWh", "", "Annual (Year 1)", "", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray2_dc_mismatch_loss", "Subarray 2 DC mismatch loss", "kWh", "", "Annual (Year 1)", "", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray2_dc_diodes_loss", "Subarray 2 DC diodes and connections loss", "kWh", "", "Annual (Year 1)", "", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray2_dc_wiring_loss", "Subarray 2 DC wiring loss", "kWh", "", "Annual (Year 1)", "", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray2_dc_tracking_loss", "Subarray 2 DC tracking loss", "kWh", "", "Annual (Year 1)", "", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray2_dc_nameplate_loss", "Subarray 2 DC nameplate loss", "kWh", "", "Annual (Year 1)", "", "", "" },

	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray3_dc_gross", "Subarray 3 gross DC energy", "kWh", "", "Annual (Year 1)", "", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray3_dc_mismatch_loss", "Subarray 3 DC mismatch loss", "kWh", "", "Annual (Year 1)", "", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray3_dc_diodes_loss", "Subarray 3 DC diodes and connections loss", "kWh", "", "Annual (Year 1)", "", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray3_dc_wiring_loss", "Subarray 3 DC wiring loss", "kWh", "", "Annual (Year 1)", "", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray3_dc_tracking_loss", "Subarray 3 DC tracking loss", "kWh", "", "Annual (Year 1)", "", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray3_dc_nameplate_loss", "Subarray 3 DC nameplate loss", "kWh", "", "Annual (Year 1)", "", "", "" },

	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray4_dc_gross", "Subarray 4 gross DC energy", "kWh", "", "Annual (Year 1)", "", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray4_dc_mismatch_loss", "Subarray 4 DC mismatch loss", "kWh", "", "Annual (Year 1)", "", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray4_dc_diodes_loss", "Subarray 4 DC diodes and connections loss", "kWh", "", "Annual (Year 1)", "", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray4_dc_wiring_loss", "Subarray 4 DC wiring loss", "kWh", "", "Annual (Year 1)", "", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray4_dc_tracking_loss", "Subarray 4 DC tracking loss", "kWh", "", "Annual (Year 1)", "", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_subarray4_dc_nameplate_loss", "Subarray 4 DC nameplate loss", "kWh", "", "Annual (Year 1)", "", "", "" },

	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_mismatch_loss", "DC mismatch loss", "kWh", "", "Annual (Year 1)", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_diodes_loss", "DC diodes and connections loss", "kWh", "", "Annual (Year 1)", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_wiring_loss", "DC wiring loss", "kWh", "", "Annual (Year 1)", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_tracking_loss", "DC tracking loss", "kWh", "", "Annual (Year 1)", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_nameplate_loss", "DC nameplate loss", "kWh", "", "Annual (Year 1)", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_optimizer_loss", "DC power optimizer loss", "kWh", "", "Annual (Year 1)", "*", "", "" },

	// loss diagram energy outputs: nominal poa, nominal array at STC, net dc, net ac, system output annual_poa_nom, annual_dc_nominal, annual_dc_net, annual_ac_net, annual_energy
	// loss diagram % losses: annual_poa_nom

	{ SSC_OUTPUT, SSC_NUMBER, "annual_poa_shading_loss_percent", "POA shading loss", "%", "", "Loss", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_poa_soiling_loss_percent", "POA soiling loss", "%", "", "Loss", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_module_loss_percent", "DC module modeled loss", "%", "", "Loss", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_snow_loss_percent", "DC snow loss", "%", "", "Loss", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_mppt_clip_loss_percent", "DC inverter MPPT clipping loss", "%", "", "Loss", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_mismatch_loss_percent", "DC mismatch loss", "%", "", "Loss", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_diodes_loss_percent", "DC diodes and connections loss", "%", "", "Loss", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_wiring_loss_percent", "DC wiring loss", "%", "", "Loss", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_tracking_loss_percent", "DC tracking loss", "%", "", "Loss", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_nameplate_loss_percent", "DC nameplate loss", "%", "", "Loss", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_optimizer_loss_percent", "DC power optimizer loss", "%", "", "Loss", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_perf_adj_loss_percent", "DC performance adjustment loss", "%", "", "Loss", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_lifetime_loss_percent", "Lifetime daily DC loss- year 1", "%", "", "Loss", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_battery_loss_percent", "DC connected battery loss- year 1", "%", "", "Loss", "*", "", "" },

	//annual_dc_net
	{ SSC_OUTPUT, SSC_NUMBER, "annual_ac_inv_clip_loss_percent", "AC inverter power clipping loss", "%", "", "Loss", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_ac_inv_pso_loss_percent", "AC inverter power consumption loss", "%", "", "Loss", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_ac_inv_pnt_loss_percent", "AC inverter night tare loss", "%", "", "Loss", "*", "", "" },

	// annual_ac_gross
	{ SSC_OUTPUT, SSC_NUMBER, "annual_ac_inv_eff_loss_percent", "AC inverter efficiency loss", "%", "", "Loss", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_ac_wiring_loss_percent", "AC wiring loss", "%", "", "Loss", "*", "", "" },
//	{ SSC_OUTPUT, SSC_NUMBER, "annual_ac_transformer_loss_percent", "AC step-up transformer loss", "%", "", "Loss", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_ac_lifetime_loss_percent", "Lifetime daily AC loss- year 1", "%", "", "Loss", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_ac_battery_loss_percent", "AC connected battery loss- year 1", "%", "", "Loss", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_xfmr_loss_percent", "Transformer loss percent", "%", "", "Loss", "", "", "" },

	// annual_ac_net
	{ SSC_OUTPUT, SSC_NUMBER, "annual_ac_perf_adj_loss_percent", "AC performance adjustment loss", "%", "", "Loss", "*", "", "" },
	// annual_energy

	/*
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_after_mismatch_loss", "DC output after mismatch loss", "kWh", "", "Annual (Year 1)", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_after_diodes_loss", "DC output after diodes and connections loss", "kWh", "", "Annual (Year 1)", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_after_wiring_loss", "DC output after wiring loss", "kWh", "", "Annual (Year 1)", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_after_tracking_loss", "DC output after tracking loss", "kWh", "", "Annual (Year 1)", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_after_nameplate_loss", "DC output after nameplate loss", "kWh", "", "Annual (Year 1)", "*", "", "" },

	{ SSC_OUTPUT, SSC_NUMBER, "annual_ac_after_inv_cliploss", "AC output after inverter clipping loss", "kWh", "", "Annual (Year 1)", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_ac_after_inv_psoloss", "AC output after inverter power consumption loss", "kWh", "", "Annual (Year 1)", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_ac_after_inv_pntloss", "AC output after inverter night tare loss", "kWh", "", "Annual (Year 1)", "*", "", "" },
	*/
	{ SSC_OUTPUT, SSC_NUMBER, "annual_ac_wiring_loss", "AC wiring loss", "kWh", "", "Annual (Year 1)", "*", "", "" },
//	{ SSC_OUTPUT, SSC_NUMBER, "annual_ac_transformer_loss", "AC step-up transformer loss", "kWh", "", "Annual (Year 1)", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_dc_optimizer_loss", "DC power optimizer loss", "kWh", "", "Annual (Year 1)", "*", "", "" },

	/*
	{ SSC_OUTPUT, SSC_NUMBER, "annual_ac_after_wiring_loss", "AC output after wiring loss", "kWh", "", "Annual (Year 1)", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "annual_ac_after_transformer_loss", "AC output after step-up transformer loss", "kWh", "", "Annual (Year 1)", "*", "", "" },
	*/

	//

	{ SSC_OUTPUT,        SSC_NUMBER,     "6par_a",                                      "CEC 6-parameter: a",        "",       "", "Module CEC 6-parameter model parameters",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_NUMBER,     "6par_Io",                                     "CEC 6-parameter: Io",       "",       "", "Module CEC 6-parameter model parameters",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_NUMBER,     "6par_Il",                                     "CEC 6-parameter: Il",       "",       "", "Module CEC 6-parameter model parameters",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_NUMBER,     "6par_Rs",                                     "CEC 6-parameter: Rs",       "",       "", "Module CEC 6-parameter model parameters",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_NUMBER,     "6par_Rsh",                                    "CEC 6-parameter: Rsh",      "",       "", "Module CEC 6-parameter model parameters",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_NUMBER,     "6par_Adj",                                    "CEC 6-parameter: Adj",      "",       "", "Module CEC 6-parameter model parameters",       "*",                    "",                              "" },
																												     
	{ SSC_OUTPUT,        SSC_NUMBER,     "performance_ratio",                           "Performance ratio",         "",       "",  "Annual (Year 1)",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_NUMBER,     "capacity_factor",                             "Capacity factor",           "%",      "",  "Annual (Year 1)", "*", "", "" },
	{ SSC_OUTPUT,        SSC_NUMBER,     "kwh_per_kw",                                  "First year kWh/kW",         "kWh/kW", "",	"Annual (Year 1)", "*", "", "" },

	//miscellaneous outputs
	{ SSC_OUTPUT,        SSC_NUMBER,      "ts_shift_hours",                            "Sun position time offset",   "hours",  "",  "Miscellaneous", "*",                       "",                          "" },
	{ SSC_OUTPUT,        SSC_NUMBER,     "nameplate_dc_rating",                        "System nameplate DC rating", "kW",     "",  "Miscellaneous",       "*",                    "",                              "" },

// test outputs
#ifdef SHADE_DB_OUTPUTS
	// ShadeDB validation

	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray1_gpoa", "ShadeDB subarray 1 global poa input", "W/m2", "", "Time Series (Subarray 1)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray1_dpoa", "ShadeDB subarray 1 diffuse poa input", "W/m2", "", "Time Series (Subarray 1)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray1_pv_cell_temp", "ShadeDB subarray 1 pv cell temp input", "C", "", "Time Series (Subarray 1)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray1_mods_per_str", "ShadeDB subarray 1 modules per string input", "", "", "Time Series (Subarray 1)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray1_str_vmp_stc", "ShadeDB subarray 1 string Vmp at STC input", "V", "", "Time Series (Subarray 1)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray1_mppt_lo", "ShadeDB subarray 1 MPPT low input", "V", "", "Time Series (Subarray 1)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray1_mppt_hi", "ShadeDB subarray 1 MPPT high input", "V", "", "Time Series (Subarray 1)", "", "", "" },

	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray2_gpoa", "ShadeDB subarray 2 global poa input", "W/m2", "", "Time Series (Subarray 2)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray2_dpoa", "ShadeDB subarray 2 diffuse poa input", "W/m2", "", "Time Series (Subarray 2)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray2_pv_cell_temp", "ShadeDB subarray 2 pv cell temp input", "C", "", "Time Series (Subarray 2)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray2_mods_per_str", "ShadeDB subarray 2 modules per string input", "", "", "Time Series (Subarray 2)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray2_str_vmp_stc", "ShadeDB subarray 2 string Vmp at STC input", "V", "", "Time Series (Subarray 2)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray2_mppt_lo", "ShadeDB subarray 2 MPPT low input", "V", "", "Time Series (Subarray 2)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray2_mppt_hi", "ShadeDB subarray 2 MPPT high input", "V", "", "Time Series (Subarray 2)", "", "", "" },

	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray3_gpoa", "ShadeDB subarray 3 global poa input", "W/m2", "", "Time Series (Subarray 3)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray3_dpoa", "ShadeDB subarray 3 diffuse poa input", "W/m2", "", "Time Series (Subarray 3)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray3_pv_cell_temp", "ShadeDB subarray 3 pv cell temp input", "C", "", "Time Series (Subarray 3)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray3_mods_per_str", "ShadeDB subarray 3 modules per string input", "", "", "Time Series (Subarray 3)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray3_str_vmp_stc", "ShadeDB subarray 3 string Vmp at STC input", "V", "", "Time Series (Subarray 3)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray3_mppt_lo", "ShadeDB subarray 3 MPPT low input", "V", "", "Time Series (Subarray 3)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray3_mppt_hi", "ShadeDB subarray 3 MPPT high input", "V", "", "Time Series (Subarray 3)", "", "", "" },


	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray4_gpoa", "ShadeDB subarray 4 global poa input", "W/m2", "", "Time Series (Subarray 4)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray4_dpoa", "ShadeDB subarray 4 diffuse poa input", "W/m2", "", "Time Series (Subarray 4)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray4_pv_cell_temp", "ShadeDB subarray 4 pv cell temp input", "C", "", "Time Series (Subarray 4)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray4_mods_per_str", "ShadeDB subarray 4 modules per string input", "", "", "Time Series (Subarray 4)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray4_str_vmp_stc", "ShadeDB subarray 4 string Vmp at STC input", "V", "", "Time Series (Subarray 4)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray4_mppt_lo", "ShadeDB subarray 4 MPPT low input", "V", "", "Time Series (Subarray 4)", "", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "shadedb_subarray4_mppt_hi", "ShadeDB subarray 4 MPPT high input", "V", "", "Time Series (Subarray 4)", "", "", "" },

#endif

	// a couple debugging outputs
	/*
	{ SSC_OUTPUT,        SSC_ARRAY,      "p_nonlinear_dc_derate0",                      "SS1x dc derate",                                          "",    "",                      "pvsamv1",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "p_nonlinear_derate_X",                        "SS1x X",                                          "",    "",                      "pvsamv1",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "p_nonlinear_derate_S",                        "SS1x S",                                          "",    "",                      "pvsamv1",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "p_nonlinear_shad1xf",                         "SS1x shade fraction",                                          "",    "",                      "pvsamv1",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "p_nonlinear_Ee_ratio",                        "SS1x Ee ratio",                                          "",    "",                      "pvsamv1",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "p_nonlinear_skyd1xf",                         "SS1x skydiff derate",                                          "",    "",                      "pvsamv1",       "*",                    "",                              "" },
	{ SSC_OUTPUT,        SSC_ARRAY,      "p_nonlinear_gndd1xf",                         "SS1x gnddiff derate",                                          "",    "",                      "pvsamv1",       "*",                    "",                              "" },
	*/

var_info_invalid };

cm_pvsamv2::cm_pvsamv2()
{
	add_var_info( _cm_vtab_pvsamv2 );
	add_var_info(vtab_adjustment_factors);
	add_var_info(vtab_dc_adjustment_factors);
	add_var_info(vtab_technology_outputs);
	add_var_info(vtab_battery_inputs);
	add_var_info(vtab_battery_outputs);
}

void cm_pvsamv2::exec( ) throw (compute_module::general_error)
{	
	// Create the IO (input/output) class for the PV simulation, this handles all inputs and outputs for the entire compute module
	std::unique_ptr<PVIOManager> IOManager(new PVIOManager(this));
	std::unique_ptr<PVSystem> PVSystem(IOManager); //sorry Nick for all the things that are broken! I think we'd only need to pass in the IOManager?

	// Temporal loop governing the simulation
	size_t weatherFileIndex = 0; /// Index in the weather file, which repeats every year of the simulation
	size_t outputIndex = 0; ///current index in the timeseries outputs, which can go for the lifetime of the system in lifetime mode
	
	for (size_t year = 0; year < IOManager->SimulationIO->numberOfYears; year++) {
		for (size_t hour = 0; hour < util::hours_per_year; hour++) {
			for (size_t step = 0; step < IOManager->getSimulationIO->stepsPerHour; step++)
			{
				//read in and error check the weather file data
				PVSystem->WeatherFile->ReadAndCheckLine();
				
				//calculate the solar position
				PVSystem->SolarPosition->CalculateSolarPosition();

				//DC Side of System
				PVSystem->PVDCController->RunSingleStep();

				//AC Side of System
				PVSystem->PVACController->calculate(); //this would do inverter temperature, inverter power, and AC battery?

				//Transformer
				PVSystem->Transformer->calculate();

				//AC Curtailment
				PVSystem->Grid->curtailment();


				//Increment indices-- I think these indices need to be passed in above?
				weatherFileIndex++;
				outputIndex++;
			}

		}

	}

}
	
DEFINE_MODULE_ENTRY( pvsamv2, "Photovoltaic performance model, SAM component models V.2", 1 )
