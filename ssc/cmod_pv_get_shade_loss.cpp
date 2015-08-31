// To duplicate functionality of GetShadeLoss.m from Sara for shade loss database
/* from GetShadeLoss.m:
GetShadeLoss(G,D,Tc,ModsPerString,StrShade,VMaxSTCStrUnshaded,VStrMPPT,ShadeDB )
% This searches the shade database and returns the %loss from partial shading
% G is the global POA irradiance, D is the diffuse irradiance, Tc is PV cell
% temperature, StrShade is a vector with each string's shaded fraction (like 24, 55, 12, etc preferably in terms of byp diode substrs),
%gammaPmp is the temperature coefficient of max power,
% reported in datasheet, VMaxSTCStrUnshaded is the unshaded Vmp of the string at STC,
% VStrMPPT is the lower and upper bounds of the inverter's MPPT range, and
% Shade DB is the database of shading losses (created by the DBX scripts at NREL)
*/
#include <functional>   // std::greater
#include <algorithm>    // std::sort

#include "core.h"
#include "lib_util.h"
#include "lib_pv_shade_loss_db.h"

static var_info _cm_vtab_pv_get_shade_loss[] = {
/*   VARTYPE           DATATYPE         NAME                           LABEL                                UNITS     META                      GROUP                      REQUIRED_IF                 CONSTRAINTS                      UI_HINTS*/
	{ SSC_INPUT, SSC_NUMBER, "global_poa_irrad", "Global POA irradiance", "", "", "PV Shade Loss DB", "*", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "diffuse_irrad", "Diffuse irradiance", "", "", "PV Shade Loss DB", "*", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "pv_cell_temp", "PV cell temperature", "", "", "PV Shade Loss DB", "*", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "mods_per_string", "Modules per string", "", "", "PV Shade Loss DB", "*", "INTEGER,MIN=1", "" },
	{ SSC_INPUT, SSC_ARRAY, "str_shade_fracs", "Shading fractions for each string", "", "", "PV Shade Loss DB", "*", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "str_vmp_stc", "Unshaded Vmp of the string at STC", "", "", "PV Shade Loss DB", "*", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "v_mppt_low", "Lower bound of inverter MPPT range", "", "", "PV Shade Loss DB", "*", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "v_mppt_high", "Upper bound of inverter MPPT range", "", "", "PV Shade Loss DB", "*", "", "" },

	// testing indices from lookup
	{ SSC_OUTPUT, SSC_NUMBER, "N", "N", "", "", "PV Shade Loss DB", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "d", "d", "", "", "PV Shade Loss DB", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "t", "t", "", "", "PV Shade Loss DB", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "S", "S", "", "", "PV Shade Loss DB", "*", "", "" },


	/* can extend to full 8760 shading for each string*/
	{ SSC_OUTPUT, SSC_NUMBER, "shade_loss", "Shade loss fraction", "", "", "PV Shade Loss DB", "*", "", "" },


var_info_invalid };

class cm_pv_get_shade_loss : public compute_module
{
public:

	cm_pv_get_shade_loss()
	{
		add_var_info(_cm_vtab_pv_get_shade_loss);
	}

	void exec() throw(general_error)
	{



		int N = 0;
		int d = 0;
		int t = 0;
		int S = 0;


// start here

		double diffuse = as_double("diffuse_irrad");
		double global = as_double("global_poa_irrad");


		ssc_number_t shade_loss = 0;
		std::vector<double> dbl_str_shade = as_doublevec("str_shade_fracs");
		size_t num_strings = dbl_str_shade.size();

		if ((num_strings > 0) && (global > 0))
		{
			//Sort in descending order of shading
			std::sort(dbl_str_shade.begin(), dbl_str_shade.end(), std::greater<double>());
			//Need to round them to 10s (note should be integer)
			for (size_t i = 0; i < num_strings;i++)
				dbl_str_shade[i] /= 10.0;
			std::vector<int> str_shade;
			for (size_t i = 0; i < num_strings; i++)
				str_shade.push_back((int)round(dbl_str_shade[i]));
//			str_shade.push_back((int)dbl_str_shade[i]);
			int s_max = -1; // = str_shade[0]
			int s_sum = 0; // = str_shade[0] that is if first element zero then sum shoudl be zero
			for (size_t i = 0; i < num_strings; i++)
			{
				if (str_shade[i] > s_max) s_max = str_shade[i];
				s_sum += str_shade[i];
			}
			//Now get the indices for the DB
			if (s_sum > 0)
			{
				int diffuse_frac = (int)round(diffuse / global * 10.0);
				if (diffuse_frac < 1) diffuse_frac = 1;
				int counter = 1;
				bool found = false;
				if (num_strings > 1)
				{
					counter = 0;
					for (int i2 = 0; i2 <= s_max; i2++)
					{
						if (num_strings == 2)
						{
							counter++;
							std::vector<int> cur_case{ s_max, i2 };
							if (str_shade == cur_case)
								found = true;
						}
						else
						{
							for (int i3 = 0; i3 <= i2; i3++)
							{
								if (num_strings == 3)
								{
									counter++;
									std::vector<int> cur_case{ s_max, i2, i3 };
									if (str_shade == cur_case)
										found = true;
								}
								else
								{
									for (int i4 = 0; i4 <= i3; i4++)
									{
										if (num_strings == 4)
										{
											counter++;
											std::vector<int> cur_case{ s_max, i2, i3, i4 };
											if (str_shade == cur_case)
												found = true;
										}
										else
										{
											for (int i5 = 0; i5 <= i4; i5++)
											{
												if (num_strings == 5)
												{
													counter++;
													std::vector<int> cur_case{ s_max, i2, i3, i4, i5 };
													if (str_shade == cur_case)
														found = true;
												}
												else
												{
													for (int i6 = 0; i6 <= i5; i6++)
													{
														if (num_strings == 6)
														{
															counter++;
															std::vector<int> cur_case{ s_max, i2, i3, i4, i5, i6 };
															if (str_shade == cur_case)
																found = true;
														}
														else
														{
															for (int i7 = 0; i7 <= i6; i7++)
															{
																if (num_strings == 7)
																{
																	counter++;
																	std::vector<int> cur_case{ s_max, i2, i3, i4, i5, i6, i7 };
																	if (str_shade == cur_case)
																		found = true;
																}
																else
																{
																	for (int i8 = 0; i8 <= i7; i8++)
																	{
																		if (num_strings == 8)
																		{
																			counter++;
																			std::vector<int> cur_case{ s_max, i2, i3, i4, i5, i6, i7, i8 };
																			if (str_shade == cur_case)
																				found = true;
																		}
																		else
																		{
																			// error message or throw error
																			counter = 0;
																		}
																	} // for i7
																	if (found) break;

																}
															} // for i7
															if (found) break;
														}
													} // for i6
													if (found) break;
												}
											} // for i5
											if (found) break;
										}
									} // for i4
									if (found) break;
								}
							} // for i3
							if (found) break;
						}
						if (found) break;
					} // for i2
				} // (num_strings > 1)


				N = num_strings;
				d = diffuse_frac;
				t = s_max;
				S = counter;

				DB8 db8;
				db8.init();
				std::vector<double>test_vmpp = db8.get_vector(N, d, t, S, DB8::VMPP);
				std::vector<double>test_impp = db8.get_vector(N, d, t, S, DB8::IMPP);
				std::vector<double>test_vs = db8.get_vector(N, d, t, S, DB8::VS);
				std::vector<double>test_is = db8.get_vector(N, d, t, S, DB8::IS);

			} //(sum >0)
		} //  ((num_strings > 0) && (global > 0))

		// testing indices returned
		assign("N", (ssc_number_t)N);
		assign("d", (ssc_number_t)d);
		assign("t", (ssc_number_t)t);
		assign("S", (ssc_number_t)S);


		assign("shade_loss", shade_loss);

		/*

		

		MyStruct = DB{ NumStrings }.d{ DiffuseFrac }.t{ Smax }.shade(counter, :, : );
		Vmps = squeeze(double(MyStruct(1, 1, :)) / 1000);
		Imps = squeeze(double(MyStruct(1, 2, :)) / 1000);
		Vs = squeeze(double(MyStruct(1, 3, :)) / 1000);
		Is = squeeze(double(MyStruct(1, 4, :)) / 1000);
		%Look at the power fractions
			PmpFracs = Vmps.*Imps;
		PFracs = Vs.*Is;
		[PmaxFrac, Pmaxind] = max(PmpFracs);
		%Try scaling the voltages using the Sandia model.Taking numbers from
			%their database for the Yingli YL230.It's a similar module (mc-si,60 cell, etc)to the 
			%Trina 250 PA05 which the database was build from.But user may need more
			%input into this!!!

			n = 1.263;
		BetaVmp = -0.137*ModsPerString; %mult by ModsPerString because it's in V
			Ns = 60 * ModsPerString; %X modules, each with 60 cells
			C2 = -0.05871;
		C3 = 8.35334;
		k = 1.38066E-23; %J / K, Boltzmann's constant
			q = 1.60218E-19;  % Coulomb, elementary charge

			deltaTc = n*k*(Tc + 273.15) / q; %Thermal voltage
			TcVmpMax = Vmps(Pmaxind)*VMaxSTCStrUnshaded + C2*Ns*deltaTc*log(G / 1000) + C3*Ns*(deltaTc*log(G / 1000)) ^ 2 + BetaVmp*(Tc - 25);
		TcVmpScale = TcVmpMax. / Vmps(Pmaxind) / VMaxSTCStrUnshaded;
		TcVmps = Vmps*VMaxSTCStrUnshaded + C2*Ns*deltaTc*log(G / 1000) + C3*Ns*(deltaTc*log(G / 1000)) ^ 2 + BetaVmp*(Tc - 25);
		TcVs = Vs*TcVmpScale*VMaxSTCStrUnshaded;

		%DO NOT USE - NOT PART OF SANDIA MODEL!!
			%Scale voltage fractions by temperature(gamma is in % / deg C)
			% TcVmps = Vmps*((1 + gammaPmp*(Tc - 25)))*VMaxSTCStrUnshaded;
		%TcVs = Vs*((1 + gammaPmp*(Tc - 25)))*VMaxSTCStrUnshaded;

		%Now want to choose the point with a V in range and highest power
			%First, figure out which max power point gives lowest loss

			Veemax = TcVmps(Pmaxind);
		if and(VStrMPPT(1) <= Veemax, VStrMPPT(2) >= Veemax)
			% The global max power point is in range!
			ShadeLoss = 1 - PmaxFrac;
		elseif isempty(PFracs(and(TcVs >= VStrMPPT(1), TcVs <= VStrMPPT(2))))
			ShadeLoss = 1;
		else
			%    %The global max power point is NOT in range
			ShadeLoss = 1 - max(PFracs(and(TcVs >= VStrMPPT(1), TcVs <= VStrMPPT(2))));
		end
		*/

	}
};

DEFINE_MODULE_ENTRY(pv_get_shade_loss, "PV get shade loss fraction for strings from Sara MacAlpine", 1)
