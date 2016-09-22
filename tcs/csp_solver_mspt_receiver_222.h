#ifndef __csp_solver_mspt_receiver_222_
#define __csp_solver_mspt_receiver_222_

#include "csp_solver_util.h"

#include "htf_props.h"
#include "ngcc_powerblock.h"
#include "csp_solver_core.h"

class C_mspt_receiver_222
{
private:
	HTFProperties field_htfProps;		// Instance of HTFProperties class for field HTF
	HTFProperties tube_material;		// Instance of HTFProperties class for receiver tube material
	HTFProperties ambient_air;			// Instance of HTFProperties class for ambient air

	ngcc_power_cycle cycle_calcs;

	double m_id_tube;
	double m_A_tube;
	int m_n_t;
	double m_A_rec_proj;
	double m_A_node;
	
	double m_Q_dot_piping_loss;		//[Wt] = Constant thermal losses from piping to env. = (THT*length_mult + length_add) * piping_loss_coef

	int m_itermode;
	double m_od_control;
	double m_tol_od;
	double m_m_dot_htf_des;
	//double m_q_rec_min;		// twn: moved to public so upstream code can retrieve

	/* declare storage variables here */
	int m_mode;
	int m_mode_prev;
	double m_E_su;
	double m_E_su_prev;
	double m_t_su;
	double m_t_su_prev;

	util::matrix_t<int> m_flow_pattern;
	int m_n_lines;

	util::matrix_t<double> m_flux_in;

	util::matrix_t<double> m_q_dot_inc;

	util::matrix_t<double> m_T_s_guess;
	util::matrix_t<double> m_T_s;
	util::matrix_t<double> m_T_panel_out_guess;
	util::matrix_t<double> m_T_panel_out;
	util::matrix_t<double> m_T_panel_in_guess;
	util::matrix_t<double> m_T_panel_in;
	util::matrix_t<double> m_T_panel_ave;
	util::matrix_t<double> m_T_panel_ave_guess;
	util::matrix_t<double> m_T_film;
	util::matrix_t<double> m_q_dot_conv;
	util::matrix_t<double> m_q_dot_rad;
	util::matrix_t<double> m_q_dot_loss;
	util::matrix_t<double> m_q_dot_abs;

	double m_m_mixed;
	double m_LoverD;
	double m_RelRough;

	// ISCC-specific
	double m_T_amb_low;
	double m_T_amb_high;
	double m_P_amb_low;
	double m_P_amb_high;
	double m_q_iscc_max;

	// member string for exception messages
	std::string error_msg;

	// track number of calls per timestep, reset = -1 in converged() call
	int m_ncall;

	//Transient model parameters
	int m_startup_mode;		
	int m_startup_mode_initial;
	int m_n_call_circ;
	int m_n_call_circ_initial;
	double m_id_riser;				//[m]
	double m_od_riser;				//[m]
	double m_id_downc;				//[m]
	double m_od_downc;				//[m]
	double m_Rinsconst_riser;		//[K*m/W]
	double m_Rinsconst_downc;		//[K*m/W]
	double m_total_startup_time;
	double m_total_startup_time_initial;
	int m_n_elem;
	int m_nz_tot;
	util::matrix_t<double> m_tm;		 //[J/K/m]
	util::matrix_t<double> m_od;		 //[m]
	util::matrix_t<double> m_id;		 //[m]
	util::matrix_t<int> m_flowelem_type;
	util::matrix_t<double> m_tm_solid;	//[J/K/m]

	struct transient_inputs
	{
		int nelem;
		int nztot;
		int npath;
		util::matrix_t<double> lam1, lam2, cval, length, zpts, tinit;
		util::matrix_t<int> nz, startpt;
	} tinputs;

	struct transient_outputs
	{
		double timeavg_tout;					// Time-averaged downcomer outlet T [K]
		double timeavg_conv_loss;				// Time-averaged convection loss from the receiver panels [W]
		double timeavg_rad_loss;				// Time-averaged radiation loss
		double timeavg_piping_loss;				// Time-averaged thermal loss from piping [W]
		double timeavg_qthermal;				// Average thermal power available during the time step [W]
		double timeavg_eta_therm;				// Time-averaged thermal efficiency of the receiver 
		double max_tout;						// Max downcomer outlet T during the time step [K]
		double min_tout;						// Min downcomer outlet T during the time step [K]
		double max_rec_tout;					// Max receiver outlet T during the time step [K]
		util::matrix_t<double> t_profile;		// Axial temperature profile at the end of the time step[K]
		util::matrix_t<double> timeavg_temp;	// Time-average outlet temperature of each flow element [K]
	} toutputs;

	struct parameter_eval_inputs
	{
		double hfor, T_amb, T_sky, nu_amb, c_htf, rho_htf, mu_htf, k_htf, Pr_htf;
		util::matrix_t<double> tm;
		util::matrix_t<double> Tfeval, Tseval;
	} pinputs;


	void calc_header_size(double presdrop, double mdot, double rhof, double muf, double Lh, double &id_calc, double &th_calc, double &od_calc);
	double interpolate(double x, const util::matrix_t<double> &xarray, const util::matrix_t<double> &yarray, int klow, int khigh);
	double integrate(double xlow, double xhigh, const util::matrix_t<double> &xarray, const util::matrix_t<double> &yarray, int klow, int khigh);
	void calc_ss_profile(double inlet_temp, const transient_inputs &tinputs, util::matrix_t<double> &tprofile);
	void calc_timeavg_temp(double inlet_temp, double tstep, const transient_inputs &tinputs, util::matrix_t<double> &timeavg);
	void calc_axial_profile(double inlet_temp, double tpt, const transient_inputs &tinputs, util::matrix_t<double> &tprofile);
	double calc_single_pt(double inlet_temp, double tpt, double zpt, int flowid, int pathid, const transient_inputs &tinputs);
	void calc_extreme_outlet_values(double inlet_temp, double tstep, const transient_inputs &tinputs, double *tmin, double *tmax);
	void update_pde_parameters(double mflow_tot, const parameter_eval_inputs &pinputs, const util::matrix_t<double> &qinc_panel, util::matrix_t<double> &Rtube, transient_inputs &tinputs);
	void solve_transient_model(double mflow_tot, double inlet_temp, double tstep, double allowable_Trise, util::matrix_t<double> &Rtube, const util::matrix_t<double> &qinc_panel, parameter_eval_inputs &pinputs, transient_inputs &tinputs, transient_outputs &toutputs);
	
	enum startup_modes
	{
		HEAT_TRACE = 0,		// No flux on receiver, riser/downcomer heated with heat tracing
		PREHEAT,			// Low flux on receiver, no HTF flow
		CIRCULATE,			// Full available power on receiver, HTF mass flow rate selected to hit target hot at SS
		HOLD				// Models predict that startup has been completed, but minimum startup time has not yet been reached.  Fluid continues to circulate through the receiver.  
	};


public:
	// Class to save messages for up stream classes
	C_csp_messages csp_messages;

	// Data
	int m_n_panels;					//[-]
	double m_d_rec;					//[m]
	double m_h_rec;					//[m]
	double m_h_tower;				//[m]
	double m_od_tube;				//[mm], convert to [m] in init()
	double m_th_tube;				//[mm], convert to [m] in init()
	double m_epsilon;				//[-]
	double m_hl_ffact;				//[-]
	double m_T_htf_hot_des;			//[C], convert to [K] in init()
	double m_T_htf_cold_des;		//[C], convert to [K] in init()
	double m_f_rec_min;				//[-]
	double m_q_rec_des;				//[MW], convert to [W] in init()
	double m_rec_su_delay;			//[-]
	double m_rec_qf_delay;			//[-]
	double m_m_dot_htf_max;			//[kg/hr], convert to [kg/s] in init()
	double m_A_sf;					//[m2]

	// 8.10.2015 twn: add tower piping thermal losses to receiver performance
	double m_pipe_loss_per_m;		//[Wt/m]
	double m_pipe_length_add;		//[m]
	double m_pipe_length_mult;		//[-]

	// Transient model 
	double m_rec_tm_mult;			//[]
	double m_u_riser;				//[m/s]
	double m_th_riser;				//[mm], convert to [m] in init()
	double m_th_downc;				//[mm], convert to [m] in init()
	double m_th_ins_pipe;			//[cm], convert to [m] in init()
	double m_ins_thermal_cond;		//[W/m/K]
	double m_piping_loss_coeff;		//[W/m2/K]
	double m_riser_tm_mult;			//[]
	double m_downc_tm_mult;			//[]
	double m_heat_trace_power;		//[kW/m], convert to [W/m] in init()
	double m_tube_flux_startup;		//[kW/m2]
	double m_preheat_temp_diff;		//[K]
	double m_startup_temp_diff;		//[K]

	int m_n_flux_x;
	int m_n_flux_y;

	// Calculate in init()
	double m_q_rec_min;				//[W]

		// 4.17.15 twn: former TCS inputs, moved to member data because are constant throughout simulation
	double m_T_salt_hot_target;			//[C], convert to K in init() call
	double m_eta_pump;					//[-]
	int m_night_recirc;					//[-]
	double m_hel_stow_deploy;			//[-]

		// Added for csp_solver/tcs wrappers:
	int m_field_fl;
	util::matrix_t<double> m_field_fl_props;	
	int m_mat_tube;
	int m_flow_type;		

		// ISCC specific
	bool m_is_iscc;
	int m_cycle_config;
	
	struct S_inputs
	{
		double m_field_eff;					//[-] 
		int m_input_operation_mode;			//[-]

		const util::matrix_t<double> *m_flux_map_input;		//[-]

		S_inputs()
		{
			m_field_eff = std::numeric_limits<double>::quiet_NaN();

			m_input_operation_mode = -1;
		}
	};

		// Let's put outputs in a structure...
	struct S_outputs
	{
		
		double m_m_dot_salt_tot;		//[kg/hr] 
		double m_eta_therm;				//[-] RECEIVER thermal efficiency
		double m_W_dot_pump;			//[MW] 
		double m_q_conv_sum;			//[MW] 
		double m_q_rad_sum;				//[MW] 
		double m_Q_thermal;				//[MW] Thermal power delivered to TES/PC: subtracts piping losses (q_dot_rec - q_dot_piping_losses)
		double m_T_salt_hot;			//[C]
		double m_field_eff_adj;			//[-]
		double m_q_dot_rec_inc;			//[MWt] Receiver incident thermal power (after reflection losses)
		double m_q_startup;				//[MWt-hr]
		double m_dP_receiver;			//[bar] receiver pressure drop
		double m_dP_total;				//[bar] total pressure drop
		double m_vel_htf;				//[m/s] HTF flow velocity through receiver tubes
		double m_T_salt_cold;			//[C] 
		double m_m_dot_ss;				//[kg/hr] 
		double m_q_dot_ss;				//[MW] 
		double m_f_timestep;			//[-]
		double m_time_required_su;		//[s]
		double m_q_dot_piping_loss;		//[MWt] Thermal power lost from piping to surroundings 

		double m_max_T_salt_hot;		//[C] Maximum salt outlet T during the time step
		double m_min_T_salt_hot;		//[C] Minimum salt outlet T during the time step
		double m_max_rec_tout;			//[C] Maximum salt T (receiver outlet) during the time step
		double m_Q_heattrace;			//[MWt] Power required for heat tracing

		S_outputs()
		{
			m_m_dot_salt_tot = m_eta_therm = m_W_dot_pump = m_q_conv_sum = m_q_rad_sum = m_Q_thermal =
				m_T_salt_hot = m_field_eff_adj = m_q_dot_rec_inc = m_q_startup = m_dP_receiver = m_dP_total =
				m_vel_htf = m_T_salt_cold = m_m_dot_ss = m_q_dot_ss = m_f_timestep = 
				m_time_required_su = m_q_dot_piping_loss = std::numeric_limits<double>::quiet_NaN();

			m_max_T_salt_hot = m_min_T_salt_hot = m_max_rec_tout = std::numeric_limits<double>::quiet_NaN();
			m_Q_heattrace = std::numeric_limits<double>::quiet_NaN();
		}
	};

	S_outputs ms_outputs;

	void clear_outputs();
	
	// Methods
	C_mspt_receiver_222();

	~C_mspt_receiver_222(){};

	void init();

	int get_operating_state();

	void call(const C_csp_weatherreader::S_outputs &weather, 
		const C_csp_solver_htf_1state &htf_state_in, 
		const C_mspt_receiver_222::S_inputs &inputs,
		const C_csp_solver_sim_info &sim_info);

	void off(const C_csp_weatherreader::S_outputs &weather,
		const C_csp_solver_htf_1state &htf_state_in,
		const C_csp_solver_sim_info &sim_info);

	void converged();

    void calc_pump_performance(double rho_f, double mdot, double ffact, double &PresDrop_calc, double &WdotPump_calc);
	
	double est_startup_time();

    HTFProperties *get_htf_property_object();

};
















#endif // __csp_solver_mspt_receiver_222_