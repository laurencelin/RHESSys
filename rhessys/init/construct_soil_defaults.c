/*--------------------------------------------------------------*/
/* 								*/
/*	construct_soil_defaults					*/
/*								*/
/*	construct_soil_defaults.c - makes patch default		*/
/*			objects.				*/
/*								*/
/*	NAME							*/
/*	construct_soil_defaults.c - makes patch default		*/
/*			objects.				*/
/*								*/
/*	SYNOPSIS						*/
/*	struct soil_default *construct_soil_defaults(           */
/*			num_default_files,			*/
/*			default_files,				*/
/*			grow_flag,				*/
/*			default_object_list )			*/
/*								*/
/*	OPTIONS							*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*	Original code, January 15, 1996.			*/
/*	July 28, 1997	C.Tague					*/
/*	removed capillary rise soil variables i.e rooting depth */
/*	pore size index and suction				*/
/*								*/
/*	Sep 15 97 RAF						*/
/*	added cap rise variables back as we now do a cap rise	*/
/*	routine but removed wilting point which was never used.	*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include "rhessys.h"
#include "params.h"

struct soil_default *construct_soil_defaults(
			int	num_default_files,
			char	**default_files,
			struct command_line_object *command_line) 
{
	/*--------------------------------------------------------------*/
	/*	Local function definition.				*/
	/*--------------------------------------------------------------*/
	void	*alloc(	size_t,
		char	*,
		char	*);
	
	double compute_delta_water(int, double, double,	double, double, double);
	int	parse_albedo_flag( char *);
	
	/*--------------------------------------------------------------*/
	/*	Local variable definition.				*/
	/*--------------------------------------------------------------*/
    int strbufLen = 256;
    int filenameLen = 1024;
    int	i;
	double 	ftmp,soil;
	FILE	*default_file;
    char	strbuf[strbufLen];
    char	outFilename[filenameLen];
	char	*newrecord;
	char	record[MAXSTR];
	struct 	soil_default *default_object_list;
	void	*alloc(	size_t, char *, char *);
    param *paramPtr = NULL;
    int paramCnt = 0;

    double p0;
    double p_decay;
    double p_decay_1;
    double vksat0;
    double vksat_decay;
    double vksat_decay_1;
    double hksat0;
    double hksat_decay;
    double hksat_decay_1;
    double PAE;
    double PAE_1;
    double PSI;
    double soildepth;
    double p3;
    double p4;
    double step;
    int len, ii, jj;
    int rt_len;
    int starting_index;
    double zzz;
    double www;
    double por_numintegrated, x_numintegrated, x_numintegrated2, x_numintegrated3, inverse_value1, inverse_value2;
    double x_numintegrated4;
    double correct_ratio;
    double max_sat_def;
    double layer_por[maxSoilDepthIndex];
    double layer_numintegrated[maxSoilDepthIndex];
    double layer_numintegrated2[maxSoilDepthIndex];
    double layer_numintegrated3[maxSoilDepthIndex];
    
	/*--------------------------------------------------------------*/
	/*	Allocate an array of default objects.			*/
	/*--------------------------------------------------------------*/
	default_object_list = (struct soil_default *)
		alloc(num_default_files * sizeof(struct soil_default),"default_object_list",
		"construct_soil_defaults");
	
	/*--------------------------------------------------------------*/
	/*	Loop through the default files list.			*/
	/*--------------------------------------------------------------*/
	for (i=0 ; i<num_default_files; i++){
		/*--------------------------------------------------------------*/
		/*		Try to open the ith default file.		*/
		/*--------------------------------------------------------------*/
		//if ( (default_file = fopen( default_files[i], "r")) == NULL ){
		//	fprintf(stderr,"FATAL ERROR:in construct_soil_defaults",
		//		"unable to open defaults file %d.\n",i);
		//	exit(EXIT_FAILURE);
		//} /*end if*/

		printf("Reading %s\n", default_files[i]);
                paramCnt = 0;
                if (paramPtr != NULL)
                    free(paramPtr);

                paramPtr = readParamFile(&paramCnt, default_files[i]);

		/*--------------------------------------------------------------*/
		/*		read the ith default file into the ith object.	*/
		/*--------------------------------------------------------------*/
		default_object_list[i].ID = 			getIntParam(&paramCnt, &paramPtr, "patch_default_ID", "%d", 3, 1); 
		default_object_list[i].theta_psi_curve = 	getIntParam(&paramCnt, &paramPtr, "theta_psi_curve", "%d", 1.0, 1);
		default_object_list[i].Ksat_0 = 		getDoubleParam(&paramCnt, &paramPtr, "Ksat_0", "%lf", 3.0, 1);
		default_object_list[i].cap_rise_max = 		getDoubleParam(&paramCnt, &paramPtr, "cap_rise_max", "%lf", 0.8, 1);
		default_object_list[i].m = 			getDoubleParam(&paramCnt, &paramPtr, "m", "%lf", 0.12, 1);
		default_object_list[i].porosity_0 = 		getDoubleParam(&paramCnt, &paramPtr, "porosity_0", "%lf", 0.435, 1);
		default_object_list[i].porosity_decay = 	getDoubleParam(&paramCnt, &paramPtr, "porosity_decay", "%lf", 4000.0, 1);
		default_object_list[i].p3 = 			getDoubleParam(&paramCnt, &paramPtr, "P3", "%lf", 0.0, 1); // param name upper case in param file
		default_object_list[i].pore_size_index = 	getDoubleParam(&paramCnt, &paramPtr, "pore_size_index", "%lf", 0.204, 1);
		default_object_list[i].psi_air_entry = 		getDoubleParam(&paramCnt, &paramPtr, "psi_air_entry", "%lf", 0.218, 1);
		default_object_list[i].psi_max = 		getDoubleParam(&paramCnt, &paramPtr, "psi_max", "%lf", 0.01, 1);
		default_object_list[i].soil_depth = 		getDoubleParam(&paramCnt, &paramPtr, "soil_depth", "%lf", 200.0, 1);
		default_object_list[i].m_z = 			getDoubleParam(&paramCnt, &paramPtr, "m_z", "%lf", 0.4, 1);
		default_object_list[i].detention_store_size = 	getDoubleParam(&paramCnt, &paramPtr, "detention_store_size", "%lf", 0.0, 1);
		default_object_list[i].deltaz = 		getDoubleParam(&paramCnt, &paramPtr, "deltaZ", "%lf", 1.0, 1); // param name contains uppercase "Z" in param file
		default_object_list[i].active_zone_z = 		getDoubleParam(&paramCnt, &paramPtr, "active_zone_z", "%lf", 5.0, 1);
		if(default_object_list[i].active_zone_z > default_object_list[i].soil_depth){
		    default_object_list[i].active_zone_z = default_object_list[i].soil_depth;
		}

		if (abs(default_object_list[i].active_zone_z - default_object_list[i].soil_depth) > 0.5) {
			printf("\nNote that soil depth used for biogeochem cycling (active zone z)");
 			printf("\nis more than 0.5 meter different from hydrologic soil depth");
 			printf("\nfor soil default file: %s\n", default_files[i] );
			}

		default_object_list[i].maximum_snow_energy_deficit = 	getDoubleParam(&paramCnt, &paramPtr, "maximum_snow_energy_deficit", "%lf", -10.0, 1);
		default_object_list[i].snow_water_capacity = 		getDoubleParam(&paramCnt, &paramPtr, "snow_water_capacity", "%lf", 0.0, 1);
		default_object_list[i].snow_light_ext_coef = 		getDoubleParam(&paramCnt, &paramPtr, "snow_light_ext_coef", "%lf", 10000.0, 1);
		default_object_list[i].snow_melt_Tcoef = 		getDoubleParam(&paramCnt, &paramPtr, "snow_melt_Tcoef", "%lf", 0.05, 1);
		default_object_list[i].snow_albedo_flag = 	parse_albedo_flag(getStrParam(&paramCnt, &paramPtr, "snow_albedo_flag", "%s", "age", 1));
		default_object_list[i].bats_b = 		getDoubleParam(&paramCnt, &paramPtr, "bats_b", "%lf", 2.0, 1);
		default_object_list[i].bats_r3 = 		getDoubleParam(&paramCnt, &paramPtr, "bats_r3", "%lf", 0.3, 1);
/*
			 default_object_list[i].snow_melt_Tcoef *= command_line[0].tmp_value; 
*/
		default_object_list[i].max_heat_capacity = 	getDoubleParam(&paramCnt, &paramPtr, "max_heat_capacity", "%lf", 0.0, 1);
		default_object_list[i].min_heat_capacity = 	getDoubleParam(&paramCnt, &paramPtr, "min_heat_capacity", "%lf", 0.0, 1);
		default_object_list[i].albedo = 		getDoubleParam(&paramCnt, &paramPtr, "albedo", "%lf", 0.28, 1);
		default_object_list[i].NO3_adsorption_rate =	getDoubleParam(&paramCnt, &paramPtr, "NO3_adsorption_rate", "%lf", 0.0, 1); 
		default_object_list[i].N_decay_rate = 		getDoubleParam(&paramCnt, &paramPtr, "N_decay", "%lf", 0.12, 1); 
		default_object_list[i].fixed_t_mult = 		getDoubleParam(&paramCnt, &paramPtr, "fixed_t_mult", "%lf", -999, 1); 

		if (default_object_list[i].fixed_t_mult > -999) 
				printf("\n Warming using a fixed values for temperature control on soil/liter decomp %lf", 
				default_object_list[i].fixed_t_mult);
		/*
		if (command_line[0].tmp_value > ZERO)
			default_object_list[i].N_decay_rate *= command_line[0].tmp_value;
		*/
		default_object_list[i].soil_type.sand =		getDoubleParam(&paramCnt, &paramPtr, "sand", "%lf", 0.7, 1);
		default_object_list[i].soil_type.silt =		getDoubleParam(&paramCnt, &paramPtr, "silt", "%lf", 0.2, 1);
		default_object_list[i].soil_type.clay = 	getDoubleParam(&paramCnt, &paramPtr, "clay", "%lf", 0.1, 1);
		soil =  default_object_list[i].soil_type.sand
			+ default_object_list[i].soil_type.silt
			+ default_object_list[i].soil_type.clay;
		if  (abs(soil - 1.0) > ZERO) {
			fprintf(stderr,
				"FATAL ERROR:in construct_soil_defaults\n proportion sand, silt, clay = %f\n\n", soil);
			printf("\n %d -  %f %f %f %f \n",
				default_object_list[i].ID,
				default_object_list[i].N_decay_rate,
				default_object_list[i].soil_type.sand,
				default_object_list[i].soil_type.silt,
				default_object_list[i].soil_type.clay);
		} /*end if*/
		if (command_line[0].gw_flag > 0) {
			default_object_list[i].sat_to_gw_coeff = getDoubleParam(&paramCnt, &paramPtr, "sat_to_gw_coeff", "%lf", 1.0, 1);
			default_object_list[i].sat_to_gw_coeff *= command_line[0].sat_to_gw_coeff_mult;
			}

		/*-----------------------------------------------------------------------------
		 *  Fill and Spill parameters
		 *-----------------------------------------------------------------------------*/
		default_object_list[i].fs_spill =	getDoubleParam(&paramCnt, &paramPtr, "fs_spill", "%lf", 1, 1);
		default_object_list[i].fs_percolation =	getDoubleParam(&paramCnt, &paramPtr, "fs_percolation", "%lf", 1, 1);
		default_object_list[i].fs_threshold = 	getDoubleParam(&paramCnt, &paramPtr, "fs_threshold", "%lf", 0.2, 1);
		
		/*--------------------------------------------------------------*/
		/*	vertical soil m and K are initized using soil default	*/
		/*	but sensitivity analysis -s is not applied to them	*/
		/*      use -sv to change these parameters			*/
		/*--------------------------------------------------------------*/
			default_object_list[i].m_v = default_object_list[i].m;
			default_object_list[i].mz_v = default_object_list[i].m_z;
			default_object_list[i].Ksat_0_v = default_object_list[i].Ksat_0;


		/*--------------------------------------------------------------*/
		/* sensitivity adjustment of vertical drainage  soil paramters	*/
		/*--------------------------------------------------------------*/
		if (command_line[0].vsen_flag > 0) {
				default_object_list[i].m_v *= command_line[0].vsen[M];
				default_object_list[i].mz_v *= command_line[0].vsen[M];
				default_object_list[i].Ksat_0_v *= command_line[0].vsen[K];
		}

		/*--------------------------------------------------------------*/
		/* sensitivity adjustment of soil drainage paramters		*/
		/*--------------------------------------------------------------*/
		if (command_line[0].sen_flag > 0) {
				default_object_list[i].m *= command_line[0].sen[M];
				default_object_list[i].m_z *= command_line[0].sen[M];
				default_object_list[i].Ksat_0 *= command_line[0].sen[K];
				default_object_list[i].soil_depth *= command_line[0].sen[SOIL_DEPTH];
		}

		/*--------------------------------------------------------------*/
		/*      calculate water_equivalent depth of soil                */
		/*--------------------------------------------------------------*/
		default_object_list[i].soil_water_cap = compute_delta_water(
			0, default_object_list[i].porosity_0,
			default_object_list[i].porosity_decay,
			default_object_list[i].soil_depth,
			default_object_list[i].soil_depth,
			0.0);

		/*--------------------------------------------------------------*/
		/* initialization of optional default file parms		*/
		/*--------------------------------------------------------------*/
		default_object_list[i].theta_mean_std_p1 = 	getDoubleParam(&paramCnt, &paramPtr, "theta_mean_std_p1", "%lf", 0.0, 1);
		default_object_list[i].theta_mean_std_p2 = 	getDoubleParam(&paramCnt, &paramPtr, "theta_mean_std_p2", "%lf", 0.0, 1);
		default_object_list[i].gl_c = 			getDoubleParam(&paramCnt, &paramPtr, "gl_c", "%lf", 0.0062, 1);
		default_object_list[i].gsurf_slope = 		getDoubleParam(&paramCnt, &paramPtr, "gsurf_slope ", "%lf", 0.01, 1);
		default_object_list[i].gsurf_intercept = 	getDoubleParam(&paramCnt, &paramPtr, "gsurf_intercept", "%lf", 0.001, 1);
		default_object_list[i].p4 = 			getDoubleParam(&paramCnt, &paramPtr, "p4", "%lf", -1.5, 1);
		default_object_list[i].DOM_decay_rate = 	getDoubleParam(&paramCnt, &paramPtr, "DOM_decay_rate", "%lf", 0.05, 1);
		default_object_list[i].NH4_adsorption_rate =    getDoubleParam(&paramCnt, &paramPtr, "NH4_adsorption_rate", "%lf", 0.000005, 1);
		default_object_list[i].DON_production_rate = 	getDoubleParam(&paramCnt, &paramPtr, "DON_production_rate", "%lf", 0.03, 1);
		default_object_list[i].DOC_adsorption_rate = 	getDoubleParam(&paramCnt, &paramPtr, "DOC_adsorption_rate", "%lf", 0.000023, 1);
		default_object_list[i].DON_adsorption_rate = 	getDoubleParam(&paramCnt, &paramPtr, "DON_adsorption_rate", "%lf", 0.000001, 1);
		default_object_list[i].interval_size = 		getDoubleParam(&paramCnt, &paramPtr, "interval_size", "%lf", INTERVAL_SIZE, 1);
		default_object_list[i].recompute_gamma_flag = 	getIntParam(&paramCnt, &paramPtr, "recompute_gamma_flag", "%d", 1, 1);

		/*--------------------------------------------------------------*/
		/* sensitivity adjustment of vertical drainage  soil paramters	*/
		/* an  scale Pore size index and psi air entry or other parameters	*/
		/* that control moisture retention (if curve 3 is used)		*/
		/*--------------------------------------------------------------*/
		if (command_line[0].vsen_alt_flag > 0) {
			if (default_object_list[i].theta_psi_curve != 3)  {
				default_object_list[i].psi_air_entry *= command_line[0].vsen_alt[PA];
				default_object_list[i].pore_size_index *= command_line[0].vsen_alt[PO];
				if (default_object_list[i].pore_size_index >= 1.0) {
					printf("\n Sensitivity analysis giving Pore Size Index > 1.0, not allowed, setting to 1.0\n");
					default_object_list[i].pore_size_index = 0.999;
					}
			}
			else {
				default_object_list[i].p3 *= command_line[0].vsen_alt[PA];
				default_object_list[i].p4 *= command_line[0].vsen_alt[PO];
			}
		}
	
		/*--------------------------------------------------------------*/
		/*      Fire effect parameters                          	*/
		/*--------------------------------------------------------------*/
		default_object_list[i].overstory_height_thresh = getDoubleParam(&paramCnt, &paramPtr, "overstory_height_thresh", "%lf", 6, 1);
		default_object_list[i].understory_height_thresh = getDoubleParam(&paramCnt, &paramPtr, "understory_height_thresh", "%lf", 4, 1);
	
		/*--------------------------------------------------------------*/
		/*		Close the ith default file.								*/
		/*--------------------------------------------------------------*/

                memset(strbuf, '\0', strbufLen);
                strcpy(strbuf, default_files[i]);
                char *s = strbuf;
                char *y = NULL;
                char *token = NULL;
                char filename[256];
    
                // Store filename portion of path in 't'
                while ((token = strtok(s, "/")) != NULL) {
                    // Save the latest component of the filename
                    strcpy(filename, token);
                    s = NULL;
                } 
    
                // Remove the file extension, if one exists
                memset(strbuf, '\0', strbufLen);
                strcpy(strbuf, filename);
                s = strbuf;
                token = strtok(s, ".");
                if (token != NULL) {
                    strcpy(filename, token);
                }
        
                memset(outFilename, '\0', filenameLen);
            

    
            // Concatenate the output prefix with the filename of the input .def file
            // and "_soil.params"
            if (command_line[0].output_prefix != NULL) {
                strcat(outFilename, command_line[0].output_prefix);
                if (filename != NULL) {
                    strcat(outFilename, "_");
                    strcat(outFilename, filename);
                }
                strcat(outFilename, "_soil.params");
            } 
            else {
                if (filename != NULL) {
                    strcat(outFilename, "_");
                    strcat(outFilename, filename);
                }
                strcat(outFilename, "soil.params");
            }
    
        printParams(paramCnt, paramPtr, outFilename);
        
        
        //-------------------------------- lookup table --------------------------------//
        len = (int)(soildepth*1000) + 2;
        if(len>maxSoilDepthIndex){ len = maxSoilDepthIndex; printf("soil(%d) has soil depth %f m. (too deep)",default_object_list[i].ID, soildepth); }
        for( jj=0; jj<maxSoilDepthIndex; jj++){
        layer_por[jj] = 0.0;
        layer_numintegrated[jj] = 0.0;
        layer_numintegrated2[jj] = 0.0;
        layer_numintegrated3[jj] = 0.0;
        }// for loop

        // using the decay and x0 to calculate the following
        default_object_list[i].sat_def_z = (double*)calloc(1002, sizeof(double));
        default_object_list[i].sat_def = (double*)calloc(1002, sizeof(double));
        default_object_list[i].sat_def_0zm = (double*)calloc(1002, sizeof(double));
        //default_object_list[i].por_z = (double*)calloc(1002, sizeof(double)); // change to nabsorbed related

        default_object_list[i].vksat_0zm = (double*)calloc(1002, sizeof(double));
        default_object_list[i].vksat_z = (double*)calloc(1002, sizeof(double));
        default_object_list[i].exfiltration_coef = (double*)calloc(1002, sizeof(double));

        default_object_list[i].fc1_0z = (double*)calloc(1002, sizeof(double));
        default_object_list[i].fc1_030r = (double*)calloc(1002, sizeof(double));
        default_object_list[i].fc1_025r = (double*)calloc(1002, sizeof(double));
        default_object_list[i].fc1_020r = (double*)calloc(1002, sizeof(double));
        default_object_list[i].fc1_015r = (double*)calloc(1002, sizeof(double));
        default_object_list[i].fc1_010r = (double*)calloc(1002, sizeof(double));
        default_object_list[i].fc1_006r = (double*)calloc(1002, sizeof(double));
        default_object_list[i].fc1_003r = (double*)calloc(1002, sizeof(double));
        default_object_list[i].fc1_00012r = (double*)calloc(1002, sizeof(double));
        default_object_list[i].pot_caprise_0z = (double*)calloc(1002, sizeof(double));
        default_object_list[i].pot_caprise_030r = (double*)calloc(1002, sizeof(double));
        default_object_list[i].pot_caprise_025r = (double*)calloc(1002, sizeof(double));
        default_object_list[i].pot_caprise_020r = (double*)calloc(1002, sizeof(double));
        default_object_list[i].pot_caprise_015r = (double*)calloc(1002, sizeof(double));
        default_object_list[i].pot_caprise_010r = (double*)calloc(1002, sizeof(double));
        default_object_list[i].pot_caprise_006r = (double*)calloc(1002, sizeof(double));
        default_object_list[i].pot_caprise_003r = (double*)calloc(1002, sizeof(double));
        default_object_list[i].pot_caprise_00012r = (double*)calloc(1002, sizeof(double));

        default_object_list[i].transmissivity_maxdailyflux = (double*)calloc(1002, sizeof(double));
        default_object_list[i].transmissivity_dailyflux = (double*)calloc(1002, sizeof(double));
        default_object_list[i].transmissivity_dailyfluxSoluteProp = (double*)calloc(1002, sizeof(double));

        p0 = default_object_list[i].porosity_0;
        p_decay = default_object_list[i].porosity_decay;
        p_decay_1 = -1.0/default_object_list[i].porosity_decay;

        vksat0 = default_object_list[i].Ksat_0_v;
        vksat_decay = default_object_list[i].mz_v;
        vksat_decay_1 = -1.0/default_object_list[i].mz_v;

        hksat0 = default_object_list[i].Ksat_0;
        hksat_decay = default_object_list[i].m; // meter
        hksat_decay_1 = -1.0/default_object_list[i].m; // -1/meter

        PAE = default_object_list[i].psi_air_entry;
        PAE_1 = -1.0/PAE;
        PSI = default_object_list[i].pore_size_index;
        soildepth = default_object_list[i].soil_depth; // already s3 scaled at LINE 2xx above.
        p3 = default_object_list[i].p3;
        p4 = default_object_list[i].p4;
        max_sat_def = p0*p_decay*(1-exp(p_decay_1*soildepth));

        // ------------------- solute profiles in sat zone
        // postive with porosity profile and horizontial ksat profile
        double satSolute_decayRate_1 = hksat_decay_1 + p_decay_1; // -1/meter
        double satSolute_decayRate = -1.0 / satSolute_decayRate_1;
        // .. if we assume sat solute profile is different from porosity, then the solute transmissivity is different from water transmissivity


        default_object_list[i].max_sat_def_1 = 1.0/( max_sat_def );
        default_object_list[i].exfiltration_wilting_point = exp(-PSI*log(2.5/PAE));
        default_object_list[i].exfiltration_S_pow = (1/(2*PSI))+2;

        default_object_list[i].soil_water_cap = p0* p_decay*(1.0-exp(p_decay_1*soildepth));
        default_object_list[i].maxrootdepth_index = (int)(round(default_object_list[i].maxrootdepth*1000));
        default_object_list[i].maxrootdepth_sat_0z = p0* p_decay*(1.0-exp(p_decay_1*default_object_list[i].maxrootdepth));
        default_object_list[i].active_zone_index = (int)(round(default_object_list[i].active_zone_z*1000));
        default_object_list[i].active_zone_sat_0z = p0* p_decay*(1.0-exp(p_decay_1*default_object_list[i].active_zone_z));
        default_object_list[i].active_zone_sat_0z_1 = 1.0/default_object_list[i].active_zone_sat_0z;
        default_object_list[i].active_zone_omProp = (1.0-exp(-default_object_list[i].active_zone_z/default_object_list[i].DOMdecayRate))/(1.0-exp(-default_object_list[i].maxrootdepth/default_object_list[i].DOMdecayRate)); // assume 95% soil OM is within 0 - maxrootz

        double check, check2;
        for( ii=0; ii<1001; ii++){
        //if(ii%100==0) printf("soil default file: %s lookup table (%d)\n", default_files[i], ii);
        // POROSITY must be less than 1.0

        default_object_list[i].sat_def_z[ii] = ii>0? min(soildepth,fabs( p_decay * log((1-ii*0.001) + ii*0.001*exp(p_decay_1* soildepth)) )) : 0.0;

        default_object_list[i].sat_def[ii] = 0.001*ii*max_sat_def;
        //p0*p_decay*(1-exp(p_decay_1*default_object_list[i].sat_def_z[ii]));

        default_object_list[i].sat_def_0zm[ii] = default_object_list[i].sat_def_z[ii]>0? default_object_list[i].sat_def[ii]/default_object_list[i].sat_def_z[ii] : 0.0;


        //default_object_list[i].por_z[ii] = p0* exp(default_object_list[i].sat_def_z[ii]*p_decay_1);
        //p0 = default_object_list[i].porosity_0;
        //p_decay = default_object_list[i].porosity_decay;
        //p_decay_1 = -1.0/default_object_list[i].porosity_decay;


        default_object_list[i].vksat_0zm[ii] = default_object_list[i].sat_def_z[ii]>0? vksat_decay* vksat0*(1-exp(vksat_decay_1*default_object_list[i].sat_def_z[ii]))/default_object_list[i].sat_def_z[ii] : vksat0;

        default_object_list[i].vksat_z[ii] = vksat0*exp(vksat_decay_1*default_object_list[i].sat_def_z[ii]);

        default_object_list[i].exfiltration_coef[ii] = default_object_list[i].sat_def_z[ii]>PAE? sqrt( (8*default_object_list[i].sat_def_0zm[ii] * default_object_list[i].vksat_0zm[ii] * PAE) / (3.0 * (1.0 + 3.0*PSI) * (1.0 + 4*PSI))) : 0.0;

        // ii=0=top=surface; ii=1000=bottom
        // loop from z to 0 surface
        len = (int)(1000*default_object_list[i].sat_def_z[ii]) + 1;
        step = default_object_list[i].sat_def_z[ii] / (floor(1000*default_object_list[i].sat_def_z[ii]) + 1);
        if(step>0){
        for( jj = 0; jj<=len; jj++ ){ // going up from z to 0; jj=0=bottom; jj=len=top
        zzz = default_object_list[i].sat_def_z[ii] - jj*step;
        www = default_object_list[i].sat_def_z[ii] - zzz;
        layer_numintegrated[jj] = p0*exp(p_decay_1*zzz) * min(1.0,exp(PSI*log(PAE/www)));
        layer_numintegrated2[jj] = vksat0*exp(vksat_decay_1*zzz) * exp(PAE_1*www);
        layer_por[jj] = p0*exp(p_decay_1*zzz);
        }// for loop jj
        // debug: when z=0, len=1, step=0; zzz=
        // when z=0.003, len=3+1, step=0.003/4, then [jj=0;zzz=(1)0.003;www=0, jj=1;zzz=(3/4)0.003;www=(1/4)0.003,
        // jj=2;zzz=(2/4);www=(2/4), jj=3;zzz=(1/4);www=(3/4), jj=4;zzz=0;www=(4/4)

        //jj=0 -> zzz=sat_def_z & www=0 --> 1/0 problem
        layer_numintegrated[0] = layer_por[0]; // fixing 1/0 problem at the buttom
        }else{
        // fixing // when ii=0; sat_def_z=0 --> len=1 & step=0; zzz_0=0 & www_0=0;
        layer_numintegrated[0] = p0; // fc
        layer_numintegrated2[0] = vksat0;
        layer_por[0] = p0;
        layer_numintegrated[1] = p0; // fc
        layer_numintegrated2[1] = vksat0;
        layer_por[1] = p0;
        }


        // integrate from z to 0
        starting_index = 0;
        por_numintegrated = 0.0;
        x_numintegrated = 0.0;
        x_numintegrated2 = 0.0;
        for( jj = starting_index; jj<=len; jj++ ){ //jj=0=sat_def_z=bottom; jj=len=0=surface
        por_numintegrated += layer_por[jj];
        x_numintegrated += layer_numintegrated[jj];
        x_numintegrated2 += layer_numintegrated2[jj];
        }// for loop jj
        if(len >= starting_index+2) for( jj = starting_index+1; jj<=len-1; jj++ ){
        por_numintegrated += layer_por[jj];
        x_numintegrated += layer_numintegrated[jj];
        x_numintegrated2 += layer_numintegrated2[jj];
        }// for loop jj
        correct_ratio = (step>0)? default_object_list[i].sat_def[ii]/(0.5*por_numintegrated*step) : 1.0; //z=0, step=0;
        default_object_list[i].fc1_0z[ii] = (0.5*x_numintegrated*step) * correct_ratio;
        default_object_list[i].pot_caprise_0z[ii] = (0.5*x_numintegrated2*step) * correct_ratio;
        inverse_value1 = (default_object_list[i].fc1_0z[ii]>0? 1.0/default_object_list[i].fc1_0z[ii] : 1.0);
        inverse_value2 = (default_object_list[i].pot_caprise_0z[ii]>0? 1.0/default_object_list[i].pot_caprise_0z[ii] : 1.0);
        // cannot determine the fc in rtz because its a patch variable


        // integrate from 3.0 to 0
        if(default_object_list[i].sat_def_z[ii] > 3.0){
        starting_index = len-(int)(3.0/step+1);
        por_numintegrated = 0.0;
        x_numintegrated = 0.0;
        x_numintegrated2 = 0.0;
        for( jj = starting_index; jj<=len; jj++ ){
        por_numintegrated += layer_por[jj];
        x_numintegrated += layer_numintegrated[jj];
        x_numintegrated2 += layer_numintegrated2[jj];
        }// for loop jj
        if(len >= starting_index+2) for( jj = starting_index+1; jj<=len-1; jj++ ){
        por_numintegrated += layer_por[jj];
        x_numintegrated += layer_numintegrated[jj];
        x_numintegrated2 += layer_numintegrated2[jj];
        }// for loop jj
        correct_ratio = p0*p_decay*(1.0-exp(p_decay_1*3.0))/(0.5*por_numintegrated*step);
        default_object_list[i].fc1_030r[ii] = (0.5*x_numintegrated*step) * correct_ratio * inverse_value1;
        default_object_list[i].pot_caprise_030r[ii] = min(1.0,(0.5*x_numintegrated2*step) * correct_ratio * inverse_value2);
        check = default_object_list[i].sat_def[ii] - p0*p_decay*(1.0-exp(p_decay_1*3.0));
        check2 = (1.0-default_object_list[i].fc1_030r[ii])*default_object_list[i].fc1_0z[ii];
        //if( check2 > check )
        //printf("soil constr 3m (%d,%d) %f %f\n", default_object_list[i].ID,ii,check, check2);
        }else{
        default_object_list[i].fc1_030r[ii] = 1.0;
        default_object_list[i].pot_caprise_030r[ii] = 1.0;
        }

        // integrate from 2.5 to 0
        if(default_object_list[i].sat_def_z[ii] > 2.5){
        starting_index = len-(int)(2.5/step+1);
        por_numintegrated = 0.0;
        x_numintegrated = 0.0;
        x_numintegrated2 = 0.0;
        for( jj = starting_index; jj<=len; jj++ ){
        por_numintegrated += layer_por[jj];
        x_numintegrated += layer_numintegrated[jj];
        x_numintegrated2 += layer_numintegrated2[jj];
        }// for loop jj
        if(len >= starting_index+2) for( jj = starting_index+1; jj<=len-1; jj++ ){
        por_numintegrated += layer_por[jj];
        x_numintegrated += layer_numintegrated[jj];
        x_numintegrated2 += layer_numintegrated2[jj];
        }// for loop jj
        correct_ratio = p0*p_decay*(1.0-exp(p_decay_1*2.5))/(0.5*por_numintegrated*step);
        default_object_list[i].fc1_025r[ii] = (0.5*x_numintegrated*step) * correct_ratio * inverse_value1;
        default_object_list[i].pot_caprise_025r[ii] = min(1.0,(0.5*x_numintegrated2*step) * correct_ratio * inverse_value2);
        check = default_object_list[i].sat_def[ii] - p0*p_decay*(1.0-exp(p_decay_1*2.5));
        check2 = (1.0-default_object_list[i].fc1_025r[ii])*default_object_list[i].fc1_0z[ii];
        //if( check2 > check )
        //printf("soil constr 2.5m (%d,%d) %f %f\n", default_object_list[i].ID,ii,check, check2);
        }else{
        default_object_list[i].fc1_025r[ii] = 1.0;
        default_object_list[i].pot_caprise_025r[ii] = 1.0;
        }

        // integrate from 2.0 to 0
        if(default_object_list[i].sat_def_z[ii] > 2.0){
        starting_index = len-(int)(2.0/step+1);
        por_numintegrated = 0.0;
        x_numintegrated = 0.0;
        x_numintegrated2 = 0.0;
        for( jj = starting_index; jj<=len; jj++ ){
        por_numintegrated += layer_por[jj];
        x_numintegrated += layer_numintegrated[jj];
        x_numintegrated2 += layer_numintegrated2[jj];
        }// for loop jj
        if(len >= starting_index+2) for( jj = starting_index+1; jj<=len-1; jj++ ){
        por_numintegrated += layer_por[jj];
        x_numintegrated += layer_numintegrated[jj];
        x_numintegrated2 += layer_numintegrated2[jj];
        }// for loop jj
        correct_ratio = p0*p_decay*(1.0-exp(p_decay_1*2.0))/(0.5*por_numintegrated*step);
        default_object_list[i].fc1_020r[ii] = (0.5*x_numintegrated*step) * correct_ratio * inverse_value1;
        default_object_list[i].pot_caprise_020r[ii] = min(1.0,(0.5*x_numintegrated2*step) * correct_ratio * inverse_value2);
        check = default_object_list[i].sat_def[ii] - p0*p_decay*(1.0-exp(p_decay_1*2.0));
        check2 = (1.0-default_object_list[i].fc1_020r[ii])*default_object_list[i].fc1_0z[ii];
        //if( check2 > check )
        //printf("soil constr 2m (%d,%d) %f %f\n", default_object_list[i].ID,ii,check, check2);
        }else{
        default_object_list[i].fc1_020r[ii] = 1.0;
        default_object_list[i].pot_caprise_020r[ii] = 1.0;
        }

        // integrate from 1.5 to 0
        if(default_object_list[i].sat_def_z[ii] > 1.5){
        starting_index = len-(int)(1.5/step+1);
        por_numintegrated = 0.0;
        x_numintegrated = 0.0;
        x_numintegrated2 = 0.0;
        for( jj = starting_index; jj<=len; jj++ ){
        por_numintegrated += layer_por[jj];
        x_numintegrated += layer_numintegrated[jj];
        x_numintegrated2 += layer_numintegrated2[jj];
        }// for loop jj
        if(len >= starting_index+2) for( jj = starting_index+1; jj<=len-1; jj++ ){
        por_numintegrated += layer_por[jj];
        x_numintegrated += layer_numintegrated[jj];
        x_numintegrated2 += layer_numintegrated2[jj];
        }// for loop jj
        correct_ratio = p0*p_decay*(1.0-exp(p_decay_1*1.5))/(0.5*por_numintegrated*step);
        default_object_list[i].fc1_015r[ii] = (0.5*x_numintegrated*step) * correct_ratio * inverse_value1;
        default_object_list[i].pot_caprise_015r[ii] = min(1.0,(0.5*x_numintegrated2*step) * correct_ratio * inverse_value2);
        check = default_object_list[i].sat_def[ii] - p0*p_decay*(1.0-exp(p_decay_1*1.5));
        check2 = (1.0-default_object_list[i].fc1_015r[ii])*default_object_list[i].fc1_0z[ii];
        //if( check2 > check )
        //printf("soil constr 1.5m (%d,%d) %f %f\n", default_object_list[i].ID,ii,check, check2);
        }else{
        default_object_list[i].fc1_015r[ii] = 1.0;
        default_object_list[i].pot_caprise_015r[ii] = 1.0;
        }

        // integrate from 1 to 0
        if(default_object_list[i].sat_def_z[ii] > 1.0){
        starting_index = len-(int)(1.0/step+1);
        por_numintegrated = 0.0;
        x_numintegrated = 0.0;
        x_numintegrated2 = 0.0;
        for( jj = starting_index; jj<=len; jj++ ){
        por_numintegrated += layer_por[jj];
        x_numintegrated += layer_numintegrated[jj];
        x_numintegrated2 += layer_numintegrated2[jj];
        }// for loop jj
        if(len >= starting_index+2) for( jj = starting_index+1; jj<=len-1; jj++ ){
        por_numintegrated += layer_por[jj];
        x_numintegrated += layer_numintegrated[jj];
        x_numintegrated2 += layer_numintegrated2[jj];
        }// for loop jj
        correct_ratio = p0*p_decay*(1.0-exp(p_decay_1))/(0.5*por_numintegrated*step);
        default_object_list[i].fc1_010r[ii] = (0.5*x_numintegrated*step) * correct_ratio * inverse_value1;
        default_object_list[i].pot_caprise_010r[ii] = min(1.0,(0.5*x_numintegrated2*step) * correct_ratio * inverse_value2);
        check = default_object_list[i].sat_def[ii] - p0*p_decay*(1.0-exp(p_decay_1*1.0));
        check2 = (1.0-default_object_list[i].fc1_010r[ii])*default_object_list[i].fc1_0z[ii];
        //if( check2 > check )
        //printf("soil constr 31.0m (%d,%d) %f %f\n", default_object_list[i].ID,ii,check, check2);
        }else{
        default_object_list[i].fc1_010r[ii] = 1.0;
        default_object_list[i].pot_caprise_010r[ii] = 1.0;
        }

        // integrate from 0.6 to 0 --> top 60cm
        if(default_object_list[i].sat_def_z[ii] > 0.6){
        starting_index = len-(int)(0.6/step+1);
        por_numintegrated = 0.0;
        x_numintegrated = 0.0;
        x_numintegrated2 = 0.0;
        for( jj = starting_index; jj<=len; jj++ ){
        por_numintegrated += layer_por[jj];
        x_numintegrated += layer_numintegrated[jj];
        x_numintegrated2 += layer_numintegrated2[jj];
        }// for loop jj
        if(len >= starting_index+2) for( jj = starting_index+1; jj<=len-1; jj++ ){
        por_numintegrated += layer_por[jj];
        x_numintegrated += layer_numintegrated[jj];
        x_numintegrated2 += layer_numintegrated2[jj];
        }// for loop jj
        correct_ratio = p0*p_decay*(1.0-exp(p_decay_1*0.6))/(0.5*por_numintegrated*step);
        default_object_list[i].fc1_006r[ii] = (0.5*x_numintegrated*step) * correct_ratio * inverse_value1;
        default_object_list[i].pot_caprise_006r[ii] = (default_object_list[i].sat_def_z[ii] > PAE) ? min(1.0,(0.5*x_numintegrated2*step) * correct_ratio * inverse_value2) : 1.0;
        check = default_object_list[i].sat_def[ii] - p0*p_decay*(1.0-exp(p_decay_1*0.6));
        check2 = (1.0-default_object_list[i].fc1_006r[ii])*default_object_list[i].fc1_0z[ii];
        if( check2 > check ){
        //printf("soil constr 0.6m (%d,%d) %f %f\n", default_object_list[i].ID,ii,check, check2);
        //fc*(1-y)=fc*x
        //fc*x*(check/check2) = fc*(1-Y)
        //1-(1-y)*(check/check2) = Y
        default_object_list[i].fc1_006r[ii] = 1.0 - (1.0-default_object_list[i].fc1_006r[ii])*(check/check2);
        }
        }else{
        default_object_list[i].fc1_006r[ii] = 1.0;
        default_object_list[i].pot_caprise_006r[ii] = 1.0;
        }

        // integrate from 0.3 to 0 --> top 30cm
        if(default_object_list[i].sat_def_z[ii] > 0.3){
        starting_index = len-(int)(0.3/step+1);
        por_numintegrated = 0.0;
        x_numintegrated = 0.0;
        x_numintegrated2 = 0.0;
        for( jj = starting_index; jj<=len; jj++ ){
        por_numintegrated += layer_por[jj];
        x_numintegrated += layer_numintegrated[jj];
        x_numintegrated2 += layer_numintegrated2[jj];
        }// for loop jj
        if(len >= starting_index+2) for( jj = starting_index+1; jj<=len-1; jj++ ){
        por_numintegrated += layer_por[jj];
        x_numintegrated += layer_numintegrated[jj];
        x_numintegrated2 += layer_numintegrated2[jj];
        }// for loop jj
        correct_ratio = p0*p_decay*(1.0-exp(p_decay_1*0.3))/(0.5*por_numintegrated*step);
        default_object_list[i].fc1_003r[ii] = (0.5*x_numintegrated*step) * correct_ratio * inverse_value1;
        default_object_list[i].pot_caprise_003r[ii] = (default_object_list[i].sat_def_z[ii] > PAE) ? min(1.0,(0.5*x_numintegrated2*step) * correct_ratio * inverse_value2) : 1.0;
        check = default_object_list[i].sat_def[ii] - p0*p_decay*(1.0-exp(p_decay_1*0.3));
        check2 = (1.0-default_object_list[i].fc1_003r[ii])*default_object_list[i].fc1_0z[ii];
        if( check2 > check ){
        //printf("soil constr 0.2m (%d,%d) %f %f\n", default_object_list[i].ID,ii,check, check2);
        //fc*(1-y)=fc*x
        //fc*x*(check/check2) = fc*(1-Y)
        //1-(1-y)*(check/check2) = Y
        default_object_list[i].fc1_003r[ii] = 1.0 - (1.0-default_object_list[i].fc1_003r[ii])*(check/check2);
        }

        }else{
        default_object_list[i].fc1_003r[ii] = 1.0;
        default_object_list[i].pot_caprise_003r[ii] = 1.0;
        }

        // integrate from 0.012 to 0 --> top 12cm
        if(default_object_list[i].sat_def_z[ii] > 0.012){
        starting_index = len-(int)(0.012/step+1);
        por_numintegrated = 0.0;
        x_numintegrated = 0.0;
        x_numintegrated2 = 0.0;
        for( jj = starting_index; jj<=len; jj++ ){
        por_numintegrated += layer_por[jj];
        x_numintegrated += layer_numintegrated[jj];
        x_numintegrated2 += layer_numintegrated2[jj];
        }// for loop jj
        if(len >= starting_index+2) for( jj = starting_index+1; jj<=len-1; jj++ ){
        por_numintegrated += layer_por[jj];
        x_numintegrated += layer_numintegrated[jj];
        x_numintegrated2 += layer_numintegrated2[jj];
        }// for loop jj
        correct_ratio = p0*p_decay*(1.0-exp(p_decay_1*0.012))/(0.5*por_numintegrated*step);
        default_object_list[i].fc1_00012r[ii] = (0.5*x_numintegrated*step) * correct_ratio * inverse_value1;
        default_object_list[i].pot_caprise_00012r[ii] = (default_object_list[i].sat_def_z[ii] > PAE) ? min(1.0,(0.5*x_numintegrated2*step) * correct_ratio * inverse_value2) : 1.0;
        check = default_object_list[i].sat_def[ii] - p0*p_decay*(1.0-exp(p_decay_1*0.012));
        check2 = (1.0-default_object_list[i].fc1_00012r[ii])*default_object_list[i].fc1_0z[ii];
        if( check2 > check ){
        //printf("soil constr 0.01m (%d,%d) %f %f\n", default_object_list[i].ID,ii,check, check2);
        default_object_list[i].fc1_00012r[ii] = 1.0 - (1.0-default_object_list[i].fc1_00012r[ii])*(check/check2);
        }
        }else{
        default_object_list[i].fc1_00012r[ii] = 1.0;
        default_object_list[i].pot_caprise_00012r[ii] = 1.0;
        }





        // loop from z to soildepth
        len = (int)(1000*(soildepth-default_object_list[i].sat_def_z[ii])) + 1;
        step = max(0.0,(soildepth-default_object_list[i].sat_def_z[ii]) / (floor(1000*(soildepth-default_object_list[i].sat_def_z[ii])) + 1));//this could be round up
        if(ii<1000){
        for( jj = 0; jj<=len; jj++ ){
        zzz = soildepth - jj*step;
        www = soildepth - zzz;
        layer_numintegrated[jj] = jj>0? p0*exp(p_decay_1*zzz) * ( 1.0-min(1.0,exp(PSI*log(PAE/www))) ) : 0.0; // movable water in pores
        layer_numintegrated2[jj] = min(hksat0*exp(hksat_decay_1*zzz), layer_numintegrated[jj]); // conductivity vs. movable water in pores
        layer_numintegrated3[jj] = exp(satSolute_decayRate_1*zzz); // mass of the layer solute; n0 = satSolute/(1.0-exp(satSolute_decayRate_1*soildepth))
        layer_por[jj] = p0*exp(p_decay_1*zzz); // volumn of the layer water
        }// for loop jj
        // when jj=0 --> zzz=soil_depth & www=0
        }else{
        //fixing
        // when sat_def_z=0 -> not a problem
        // when sat_def_z=soildepth -> len=1 & step=0
        layer_numintegrated[0] = 0.0;
        layer_numintegrated2[0] = 0.0;
        layer_numintegrated3[0] = 0.0;
        layer_por[0] = p0*exp(p_decay_1*soildepth);
        layer_numintegrated[1] = 0.0;
        layer_numintegrated2[1] = 0.0;
        layer_numintegrated3[1] = 0.0;
        layer_por[1] = p0*exp(p_decay_1*soildepth);
        }

        // jj=0=soildepth=bottom; jj=len=z=top
        // integrate from z to soildepth
        starting_index = 0;
        por_numintegrated = 0.0;
        x_numintegrated = 0.0;
        x_numintegrated2 = 0.0;
        x_numintegrated3 = 0.0;
        x_numintegrated4 = 0.0;
        for( jj = starting_index; jj<=len; jj++ ){
        por_numintegrated += layer_por[jj];
        x_numintegrated += layer_numintegrated[jj];
        x_numintegrated2 += layer_numintegrated2[jj];
        x_numintegrated4 += layer_numintegrated3[jj];
        x_numintegrated3 += (layer_numintegrated[jj]>0? layer_numintegrated3[jj] * layer_numintegrated2[jj]/layer_numintegrated[jj] : 0.0);
        }// for loop jj
        if(len >= starting_index+2) for( jj = starting_index+1; jj<=len-1; jj++ ){
        por_numintegrated += layer_por[jj];
        x_numintegrated += layer_numintegrated[jj];
        x_numintegrated2 += layer_numintegrated2[jj];
        x_numintegrated4 += layer_numintegrated3[jj];
        x_numintegrated3 += (layer_numintegrated[jj]>0? layer_numintegrated3[jj] * layer_numintegrated2[jj]/layer_numintegrated[jj] : 0.0);
        }// for loop jj
        correct_ratio = (step>0) ? p0*p_decay*(exp(p_decay_1*default_object_list[i].sat_def_z[ii])-exp(p_decay_1*soildepth) ) /(0.5*por_numintegrated*step) : 1.0;
        default_object_list[i].transmissivity_maxdailyflux[ii] = (0.5*x_numintegrated*step) * correct_ratio;//volumn constraint; from wtz to soildepth
        default_object_list[i].transmissivity_dailyflux[ii] = (0.5*x_numintegrated2*step) * correct_ratio;//volumn and ksat constraint; from wtz to soildepth
        //transmissivity_dailyflux[ii] last one is -inf / nan sometimes
        default_object_list[i].transmissivity_dailyfluxSoluteProp[ii] = (x_numintegrated4>0? x_numintegrated3 / x_numintegrated4 : 0.0);

        }// for loop ii 0.1% loop



        rt_len = (int)(soildepth*1000) + 1; //(int)(default_object_list[i].maxrootdepth*1000)+1;
        default_object_list[i].soildepthLen = rt_len;
        default_object_list[i].rtz2sat_def_0z = (double*)calloc(rt_len, sizeof(double));
        default_object_list[i].rtz2sat_def_pct_index = (int*)calloc(rt_len, sizeof(int));
        default_object_list[i].rtz2NO3prop = (double*)calloc(rt_len, sizeof(double));
        default_object_list[i].rtz2NH4prop = (double*)calloc(rt_len, sizeof(double));
        default_object_list[i].rtz2DOMprop = (double*)calloc(rt_len, sizeof(double));


        // ------------------- solute profiles in soil column
        // assume most (95%) solutes and OM are bounded in between 0 - maxrootdepth (horizon C)
        // maxrootdepth (horizon C) >= active_zone_z (horizon B)

        // ... assume the profile is inverse of vKsat profile below the maxrootdepth zone
        double rtz_soil_vdeay_1 = 1.0/default_object_list[i].mz_v; // inverse vKsat profile (vksat_decay_1)
        double rtz_soil_coef = exp(rtz_soil_vdeay_1*default_object_list[i].maxrootdepth);

        // ... estimate the om% deeper than maxrootdepth
        double rtz_soil_om;
        if(fabs(soilcDecay)>0) rtz_soil_om = 0.1*(exp(-default_object_list[i].active_zone_z/soilcDecay) -
        exp(-default_object_list[i].maxrootdepth/soilcDecay)) / (1.0 - exp(-default_object_list[i].maxrootdepth/soilcDecay));
        // 10% of OM in between active_zone_z and maxrootdepth
        else rtz_soil_om = 0.05;
        double solute_0_rtz_soil = rtz_soil_om/(rtz_soil_coef-exp(rtz_soil_vdeay_1*soildepth));// negative when rtz_soil_vdeay_1>0

        // ... profile decay rates
        if(fabs(soilcDecay)>0){
        // use soil om profile to characterize solute profile
        default_object_list[i].NO3decayRate_1 = -1.0/soilcDecay + 1.0/default_object_list[i].mz_v;
        default_object_list[i].NH4decayRate_1 = -1.0/soilcDecay + 1.0/default_object_list[i].mz_v;
        default_object_list[i].DOMdecayRate_1 = -1.0/soilcDecay + 1.0/default_object_list[i].mz_v;

        default_object_list[i].NO3decayRate = -1.0/default_object_list[i].NO3decayRate_1;
        default_object_list[i].NH4decayRate = -1.0/default_object_list[i].NH4decayRate_1;
        default_object_list[i].DOMdecayRate = -1.0/default_object_list[i].DOMdecayRate_1;
        }else{
        // specify values
        if(fabs(default_object_list[i].NO3decayRate)>0){
        // NO3decayRate is non-zero
        default_object_list[i].NO3decayRate_1 = -1.0/default_object_list[i].NO3decayRate;
        }else{
        // NO3decayRate is zero
        default_object_list[i].NO3decayRate = -default_object_list[i].maxrootdepth/log(rtz_soil_om); //m
        default_object_list[i].NO3decayRate_1 = -1.0/default_object_list[i].NO3decayRate;
        }

        if(fabs(default_object_list[i].NH4decayRate)>0){
        // NH4decayRate is non-zero
        default_object_list[i].NH4decayRate_1 = -1.0/default_object_list[i].NH4decayRate;
        }else{
        // NH4decayRate is zero
        default_object_list[i].NH4decayRate = -default_object_list[i].maxrootdepth/log(rtz_soil_om); //m
        default_object_list[i].NH4decayRate_1 = -1.0/default_object_list[i].NH4decayRate;
        }

        if(fabs(default_object_list[i].DOMdecayRate)>0){
        // DOMdecayRate is non-zero
        default_object_list[i].DOMdecayRate_1 = -1.0/default_object_list[i].DOMdecayRate;
        }else{
        // DOMdecayRate is zero
        default_object_list[i].DOMdecayRate = -default_object_list[i].maxrootdepth/log(rtz_soil_om); //m
        default_object_list[i].DOMdecayRate_1 = -1.0/default_object_list[i].DOMdecayRate;
        }
        }// end of if
        if(fabs(default_object_list[i].NO3decayRate_1)>0){}else{ default_object_list[i].NO3decayRate_1 = -0.00025; }
        if(fabs(default_object_list[i].NH4decayRate_1)>0){}else{ default_object_list[i].NH4decayRate_1 = -0.00025; }
        if(fabs(default_object_list[i].DOMdecayRate_1)>0){}else{ default_object_list[i].DOMdecayRate_1 = -0.00025; }
        double NO3_0 = 1.0/(1.0-exp(default_object_list[i].NO3decayRate_1*default_object_list[i].maxrootdepth));
        double NH4_0 = 1.0/(1.0-exp(default_object_list[i].NH4decayRate_1*default_object_list[i].maxrootdepth));
        double DOM_0 = 1.0/(1.0-exp(default_object_list[i].DOMdecayRate_1*default_object_list[i].maxrootdepth));
        // N0 * (1-exp(z)) = soilN;

        for( ii=0; ii<rt_len; ii++){
        zzz = ii*0.001;
        default_object_list[i].rtz2sat_def_0z[ii] = p0*p_decay*(1.0-exp(p_decay_1*zzz));
        default_object_list[i].rtz2sat_def_pct_index[ii] = (int)(default_object_list[i].rtz2sat_def_0z[ii]*default_object_list[i].max_sat_def_1*1000);
        //default_object_list[i].rtz2sat_def_pct_indexM[ii] = 1000*(default_object_list[i].rtz2sat_def_0z[ii] - default_object_list[i].rtz2sat_def_pct_index[ii]*0.001);

        if( zzz > default_object_list[i].maxrootdepth){
        // (1-rtz_soil_om) is the om% in between 0 and maxrootdepth
        default_object_list[i].rtz2NO3prop[ii] = (1-rtz_soil_om)+(rtz_soil_coef-exp(rtz_soil_vdeay_1*zzz)) * solute_0_rtz_soil;
        default_object_list[i].rtz2NH4prop[ii] = (1-rtz_soil_om)+(rtz_soil_coef-exp(rtz_soil_vdeay_1*zzz)) * solute_0_rtz_soil;
        default_object_list[i].rtz2DOMprop[ii] = (1-rtz_soil_om)+(rtz_soil_coef-exp(rtz_soil_vdeay_1*zzz)) * solute_0_rtz_soil;
        }else{
        default_object_list[i].rtz2NO3prop[ii] = (1.0-exp(default_object_list[i].NO3decayRate_1*zzz)) * NO3_0 * (1-rtz_soil_om);
        default_object_list[i].rtz2NH4prop[ii] = (1.0-exp(default_object_list[i].NH4decayRate_1*zzz)) * NH4_0 * (1-rtz_soil_om);
        default_object_list[i].rtz2DOMprop[ii] = (1.0-exp(default_object_list[i].DOMdecayRate_1*zzz)) * DOM_0 * (1-rtz_soil_om);
        }
        }// for loop ii
        default_object_list[i].rtz2sat_def_0z[0] = 0.0;


        default_object_list[i].sat_def_z[1001] = default_object_list[i].sat_def_z[1000];
        default_object_list[i].sat_def[1001] = default_object_list[i].sat_def[1000];
        default_object_list[i].sat_def_0zm[1001] = default_object_list[i].sat_def_0zm[1000];
        //default_object_list[i].sat_zZ[1001] = default_object_list[i].sat_zZ[1000];

        default_object_list[i].vksat_0zm[1001] = default_object_list[i].vksat_0zm[1000];
        default_object_list[i].vksat_z[1001] = default_object_list[i].vksat_z[1000];
        default_object_list[i].exfiltration_coef[1001] = default_object_list[i].exfiltration_coef[1000];

        default_object_list[i].fc1_0z[1001] = default_object_list[i].fc1_0z[1000];
        default_object_list[i].fc1_030r[1001] = default_object_list[i].fc1_030r[1000];
        default_object_list[i].fc1_025r[1001] = default_object_list[i].fc1_025r[1000];
        default_object_list[i].fc1_020r[1001] = default_object_list[i].fc1_020r[1000];
        default_object_list[i].fc1_015r[1001] = default_object_list[i].fc1_015r[1000];
        default_object_list[i].fc1_010r[1001] = default_object_list[i].fc1_010r[1000];
        default_object_list[i].fc1_006r[1001] = default_object_list[i].fc1_006r[1000];
        default_object_list[i].fc1_003r[1001] = default_object_list[i].fc1_003r[1000];
        default_object_list[i].fc1_00012r[1001] = default_object_list[i].fc1_00012r[1000];
        default_object_list[i].pot_caprise_0z[1001] = default_object_list[i].pot_caprise_0z[1000];
        default_object_list[i].pot_caprise_030r[1001] = default_object_list[i].pot_caprise_030r[1000];
        default_object_list[i].pot_caprise_025r[1001] = default_object_list[i].pot_caprise_025r[1000];
        default_object_list[i].pot_caprise_020r[1001] = default_object_list[i].pot_caprise_020r[1000];
        default_object_list[i].pot_caprise_015r[1001] = default_object_list[i].pot_caprise_015r[1000];
        default_object_list[i].pot_caprise_010r[1001] = default_object_list[i].pot_caprise_010r[1000];
        default_object_list[i].pot_caprise_006r[1001] = default_object_list[i].pot_caprise_006r[1000];
        default_object_list[i].pot_caprise_003r[1001] = default_object_list[i].pot_caprise_003r[1000];
        default_object_list[i].pot_caprise_00012r[1001] = default_object_list[i].pot_caprise_00012r[1000];

        default_object_list[i].transmissivity_maxdailyflux[1001] = default_object_list[i].transmissivity_maxdailyflux[1000];
        default_object_list[i].transmissivity_dailyflux[1001] = default_object_list[i].transmissivity_dailyflux[1000];
        default_object_list[i].transmissivity_dailyfluxSoluteProp[1001] = default_object_list[i].transmissivity_dailyfluxSoluteProp[1000];
        printf("soil,%d,%f,%f,%f\n",default_object_list[i].ID, soildepth, default_object_list[i].vksat_z[1000], p0*exp(p_decay_1* soildepth) );
        
        
	} /*end for*/

	if (paramPtr != NULL)
            free(paramPtr);
  return(default_object_list);
} /*end construct_soil_defaults*/
