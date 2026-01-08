/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"
#include "../BioFVM/BioFVM.h"  
using namespace BioFVM;


#include "rrc_api.h"
#include "rrc_types.h"

#include <iostream>
#include <algorithm>
#include <math.h>
#include <random>
#include <chrono>
#ifndef MPI
#define MPI 3.14159265358979323846
#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// #include "rrc_utilities.h"
extern "C" rrc::RRHandle createRRInstance();

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = NULL;
	//cell_defaults.functions.update_velocity = standard_update_cell_velocity;
	//cell_defaults.functions.update_velocity = drag_update_cell_vel;
	cell_defaults.functions.update_velocity = drag_update_cell_velocity_Ines;



	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
    setup_signal_behavior_dictionaries(); 
    
    
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{	

    static int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" );
    static int glucose_substrate_index = microenvironment.find_density_index( "glucose" ); 
    static int lactate_substrate_index = microenvironment.find_density_index( "lactate");
	static int pyr_substrate_index = microenvironment.find_density_index( "Pyr" );
    static int nad_substrate_index = microenvironment.find_density_index( "NAD" ); 
    static int nadh_substrate_index = microenvironment.find_density_index( "NADH");
    
    
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create cells 
    
	Cell* pCell;
	
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.8 * 2.0 * cell_radius; 
	//double initial_tumor_radius = 1.2*cell_radius;//Tumor Radius, proportional to number of cells
	//double initial_tumor_radius = 12;//Tumor Radius, proportional to number of cells
	double initial_tumor_radius = 21;//Tumor Radius, proportional to number of cells
    double retval;
    //2D
	//std::vector<std::vector<double>> positions = create_cell_circle_positions(cell_radius,initial_tumor_radius);
    
    //3D
    std::vector<std::vector<double>> positions = create_cell_circle_positions_3D(cell_radius,initial_tumor_radius);

    std::cout << "NUMBER OF CELLS : " << positions.size() << " __________" << std::endl;
    for( int i=0; i < positions.size(); i++ )
    {
        pCell = create_cell(get_cell_definition("default")); 
        pCell->assign_position( positions[i] );

        int i_Oxy_i = pCell->custom_data.find_variable_index( "intra_oxy" );
        int i_Glu_i = pCell->custom_data.find_variable_index( "intra_glu" );
        int i_Lac_i = pCell->custom_data.find_variable_index( "intra_lac" );
        int energy_vi = pCell->custom_data.find_variable_index( "intra_energy" );
        int i_NAD_i = pCell->custom_data.find_variable_index( "intra_NAD" );
        int i_NADH_i = pCell->custom_data.find_variable_index( "intra_NADH" );
        int i_Pyr_i = pCell->custom_data.find_variable_index( "intra_pyr" );

		pCell->custom_data[i_Oxy_i] = parameters.doubles("initial_internal_oxygen");
        pCell->custom_data[i_Glu_i] = parameters.doubles("initial_internal_glucose");
        pCell->custom_data[i_Lac_i] = parameters.doubles("initial_internal_lactate");
        pCell->custom_data[energy_vi] = parameters.doubles("initial_energy");
        pCell->custom_data[i_NAD_i] = parameters.doubles("initial_internal_NAD");
        pCell->custom_data[i_NADH_i] = parameters.doubles("initial_internal_NADH");
        pCell->custom_data[i_Pyr_i] = parameters.doubles("initial_internal_pyr");
		
		
		printf("Custom Data Oxy:%lf Glu:%lf Lac:%lf Ener:%lf Pyr:%lf Nad:%lf Nadh:%lf\n",pCell->custom_data[i_Oxy_i],
				pCell->custom_data[i_Glu_i],
				pCell->custom_data[i_Lac_i],
				pCell->custom_data[energy_vi],
				pCell->custom_data[i_NAD_i], 
				pCell->custom_data[i_NADH_i],
				pCell->custom_data[i_Pyr_i] );
		
/*         pCell->custom_data[i_Oxy_i] = parameters.doubles("initial_internal_oxygen");
        pCell->custom_data[i_Glu_i] = parameters.doubles("initial_internal_glucose");
        pCell->custom_data[i_Lac_i] = parameters.doubles("initial_internal_lactate");
        pCell->custom_data[energy_vi] = parameters.doubles("initial_energy"); */
        
        double cell_volume = pCell->phenotype.volume.total;
        
        //std::cout << "oxygen custom data : " << pCell->custom_data[i_Oxy_i] << std::endl;
        //std::cout << "oxygen custom data : SIGNAL" << get_single_signal( pCell, "custom:intra_oxy") << std::endl;
       
        
        pCell->phenotype.molecular.internalized_total_substrates[oxygen_substrate_index]= get_single_signal( pCell, "custom:intra_oxy") * cell_volume;
        pCell->phenotype.molecular.internalized_total_substrates[glucose_substrate_index]= get_single_signal( pCell, "custom:intra_glu") * cell_volume;
        pCell->phenotype.molecular.internalized_total_substrates[lactate_substrate_index]= get_single_signal( pCell, "custom:intra_lac") * cell_volume;
		pCell->phenotype.molecular.internalized_total_substrates[pyr_substrate_index]= get_single_signal( pCell, "custom:intra_pyr") * cell_volume;
        pCell->phenotype.molecular.internalized_total_substrates[nad_substrate_index]= get_single_signal( pCell, "custom:intra_NAD") * cell_volume;
        pCell->phenotype.molecular.internalized_total_substrates[nadh_substrate_index]= get_single_signal( pCell, "custom:intra_NADH") * cell_volume;
        pCell->phenotype.intracellular->start();
        (*all_cells)[i]->phenotype.intracellular->set_parameter_value("Energy",get_single_signal( pCell, "custom:intra_energy"));

		
		printf("Intracell Oxy:%lf Glu:%lf Lac:%lf Ener:%lf Pyr:%lf Nad:%lf Nadh:%lf\n",pCell->phenotype.molecular.internalized_total_substrates[oxygen_substrate_index],
			pCell->phenotype.molecular.internalized_total_substrates[glucose_substrate_index],
			pCell->phenotype.molecular.internalized_total_substrates[lactate_substrate_index],
			pCell->phenotype.intracellular->get_parameter_value("Energy"),
			pCell->phenotype.molecular.internalized_total_substrates[pyr_substrate_index], 
			pCell->phenotype.molecular.internalized_total_substrates[nad_substrate_index],
			pCell->phenotype.molecular.internalized_total_substrates[nadh_substrate_index]);
		
		/*
		// Get Intracellular Concentrations
        double oxy_val_ini = parameters.doubles("initial_internal_oxygen"); 
        double glu_val_ini = parameters.doubles("initial_internal_glucose"); 
        double lac_val_ini = parameters.doubles("initial_internal_lactate"); 
		double ene_val_ini = parameters.doubles("initial_energy"); 
		double pyr_val_ini = parameters.doubles("initial_internal_pyr"); 
        double nad_val_ini = parameters.doubles("initial_internal_NAD"); 
        double nadh_val_ini = parameters.doubles("initial_internal_NADH"); 
		// Update SBML 
        (*all_cells)[i]->phenotype.intracellular->set_parameter_value("Oxygen",oxy_val_ini);
        (*all_cells)[i]->phenotype.intracellular->set_parameter_value("Glucose",glu_val_ini);
        (*all_cells)[i]->phenotype.intracellular->set_parameter_value("Lactate",lac_val_ini);
		(*all_cells)[i]->phenotype.intracellular->set_parameter_value("Energy",ene_val_ini);
		(*all_cells)[i]->phenotype.intracellular->set_parameter_value("Pyr",pyr_val_ini);
        (*all_cells)[i]->phenotype.intracellular->set_parameter_value("NAD",nad_val_ini);
        (*all_cells)[i]->phenotype.intracellular->set_parameter_value("NADH",nadh_val_ini);
		
		printf("Celula %d SBML INICIAL:\n Oxy:%lf Glu:%lf Lac:%lf Ener:%lf Pyr:%lf Nad:%lf Nadh:%lf\n", i, 
			(*all_cells)[i]->phenotype.intracellular->get_parameter_value("Oxygen"),
			(*all_cells)[i]->phenotype.intracellular->get_parameter_value("Glucose"), 
			(*all_cells)[i]->phenotype.intracellular->get_parameter_value("Lactate"), 
			(*all_cells)[i]->phenotype.intracellular->get_parameter_value("Energy"), 
			(*all_cells)[i]->phenotype.intracellular->get_parameter_value("Pyr"), 
			(*all_cells)[i]->phenotype.intracellular->get_parameter_value("NAD"),  
			(*all_cells)[i]->phenotype.intracellular->get_parameter_value("NADH"));
		*/
		
	}

	return; 
}

void update_intracellular()
{	/*
	for( int i=0; i < (*all_cells).size(); i++ )
    {
		if((*all_cells)[i]->phenotype.death.dead==true)
			printf("Celula %d muerta despues sbml\n",i);
	}
	*/
	
    // BioFVM Indices
    static int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" );
    static int glucose_substrate_index = microenvironment.find_density_index( "glucose" ); 
    static int lactate_substrate_index = microenvironment.find_density_index( "lactate");
	static int pyr_substrate_index = microenvironment.find_density_index( "Pyr" );
    static int nad_substrate_index = microenvironment.find_density_index( "NAD" ); 
    static int nadh_substrate_index = microenvironment.find_density_index( "NADH");

    #pragma omp parallel for 
    for( int i=0; i < (*all_cells).size(); i++ )
    {
        if( (*all_cells)[i]->is_out_of_domain == false)
        {
            // Cell Volume
            double cell_volume = (*all_cells)[i]->phenotype.volume.total;
            
            // Get Intracellular Concentrations
            double oxy_val_int = get_single_signal((*all_cells)[i], "intracellular oxygen"); 
            double glu_val_int = get_single_signal((*all_cells)[i], "intracellular glucose"); 
            double lac_val_int = get_single_signal((*all_cells)[i], "intracellular lactate"); 
			double pyr_val_int = get_single_signal((*all_cells)[i], "intracellular Pyr"); 
            double nad_val_int = get_single_signal((*all_cells)[i], "intracellular NAD"); 
            double nadh_val_int = get_single_signal((*all_cells)[i], "intracellular NADH"); 
			//printf("Celula %d Intracellular concentrations:\n Oxy:%lf Glu:%lf Lac:%lf Pyr:%lf Nad:%lf Nadh:%lf\n", i, oxy_val_int, glu_val_int, lac_val_int, pyr_val_int, nad_val_int,  nadh_val_int); 
            
            // Update SBML 
            (*all_cells)[i]->phenotype.intracellular->set_parameter_value("Oxygen",oxy_val_int);
            (*all_cells)[i]->phenotype.intracellular->set_parameter_value("Glucose",glu_val_int);
            (*all_cells)[i]->phenotype.intracellular->set_parameter_value("Lactate",lac_val_int);
			(*all_cells)[i]->phenotype.intracellular->set_parameter_value("Pyr",pyr_val_int);
            (*all_cells)[i]->phenotype.intracellular->set_parameter_value("NAD",nad_val_int);
            (*all_cells)[i]->phenotype.intracellular->set_parameter_value("NADH",nadh_val_int);
			
			
			/*
            printf("Celula %d SBML Entry:\n Oxy:%lf Glu:%lf Lac:%lf Ener:%lf  Pyr:%lf Nad:%lf Nadh:%lf\n", i, 
				(*all_cells)[i]->phenotype.intracellular->get_parameter_value("Oxygen"),
				(*all_cells)[i]->phenotype.intracellular->get_parameter_value("Glucose"), 
				(*all_cells)[i]->phenotype.intracellular->get_parameter_value("Lactate"), 
				(*all_cells)[i]->phenotype.intracellular->get_parameter_value("Energy"), 
				(*all_cells)[i]->phenotype.intracellular->get_parameter_value("Pyr"), 
				(*all_cells)[i]->phenotype.intracellular->get_parameter_value("NAD"),  
				(*all_cells)[i]->phenotype.intracellular->get_parameter_value("NADH"));
			*/		
			/*
			//Are the cells dead?
			if((*all_cells)[i]->phenotype.death.dead==true)
				printf("Celula %d muerta antes sbml\n",i);
			*/
            // SBML Simulation
            (*all_cells)[i]->phenotype.intracellular->update();
            
            // Phenotype Simulation
            (*all_cells)[i]->phenotype.intracellular->update_phenotype_parameters((*all_cells)[i]->phenotype);
			/*
			//Are the cells dead?
			if((*all_cells)[i]->phenotype.death.dead==true)
				printf("Celula %d muerta despues sbml\n",i);
			*/
			/*
			printf("Celula %d SBML Exit:\n Oxy:%lf Glu:%lf Lac:%lf Ener:%lf Pyr:%lf Nad:%lf Nadh:%lf\n", i, (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Oxygen"),
				(*all_cells)[i]->phenotype.intracellular->get_parameter_value("Glucose"), 
				(*all_cells)[i]->phenotype.intracellular->get_parameter_value("Lactate"), 
				(*all_cells)[i]->phenotype.intracellular->get_parameter_value("Energy"), 
				(*all_cells)[i]->phenotype.intracellular->get_parameter_value("Pyr"), 
				(*all_cells)[i]->phenotype.intracellular->get_parameter_value("NAD"),  
				(*all_cells)[i]->phenotype.intracellular->get_parameter_value("NADH")); 
			*/
            // Internalized Chemical Update After SBML Simulation
            (*all_cells)[i]->phenotype.molecular.internalized_total_substrates[oxygen_substrate_index] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Oxygen") * cell_volume;
            (*all_cells)[i]->phenotype.molecular.internalized_total_substrates[glucose_substrate_index] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Glucose") * cell_volume;
            (*all_cells)[i]->phenotype.molecular.internalized_total_substrates[lactate_substrate_index] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Lactate") * cell_volume;
			(*all_cells)[i]->phenotype.molecular.internalized_total_substrates[pyr_substrate_index] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Pyr") * cell_volume;
            (*all_cells)[i]->phenotype.molecular.internalized_total_substrates[nad_substrate_index] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("NAD") * cell_volume;
            (*all_cells)[i]->phenotype.molecular.internalized_total_substrates[nadh_substrate_index] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("NADH") * cell_volume;
			//printf("Cellular volume: %lf\n", cell_volume);
			
			/*
			printf("Celula %d Internalized sustances:\n Oxy:%lf Glu:%lf Lac:%lf Pyr:%lf Nad:%lf Nadh:%lf\n", i, 
				(*all_cells)[i]->phenotype.molecular.internalized_total_substrates[oxygen_substrate_index],
				(*all_cells)[i]->phenotype.molecular.internalized_total_substrates[glucose_substrate_index], 
				(*all_cells)[i]->phenotype.molecular.internalized_total_substrates[lactate_substrate_index], 
				(*all_cells)[i]->phenotype.molecular.internalized_total_substrates[pyr_substrate_index], 
				(*all_cells)[i]->phenotype.molecular.internalized_total_substrates[nad_substrate_index],  
				(*all_cells)[i]->phenotype.molecular.internalized_total_substrates[nadh_substrate_index]); 
			*/
            //Save custom data
            set_single_behavior( (*all_cells)[i] , "custom:intra_oxy" , (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Oxygen") ); 
            set_single_behavior( (*all_cells)[i] , "custom:intra_glu" , (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Glucose") ); 
            set_single_behavior( (*all_cells)[i] , "custom:intra_lac" , (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Lactate") ); 
            set_single_behavior( (*all_cells)[i] , "custom:intra_energy" , (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Energy") ); 
			set_single_behavior( (*all_cells)[i] , "custom:intra_pyr" , (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Pyr") ); 
            set_single_behavior( (*all_cells)[i] , "custom:intra_NAD" , (*all_cells)[i]->phenotype.intracellular->get_parameter_value("NAD") ); 
            set_single_behavior( (*all_cells)[i] , "custom:intra_NADH" , (*all_cells)[i]->phenotype.intracellular->get_parameter_value("NADH") );
			
			/*
			int i_Oxy_i = (*all_cells)[i]->custom_data.find_variable_index( "intra_oxy" );
       		int i_Glu_i = (*all_cells)[i]->custom_data.find_variable_index( "intra_glu" );
        	int i_Lac_i = (*all_cells)[i]->custom_data.find_variable_index( "intra_lac" );
        	int energy_vi = (*all_cells)[i]->custom_data.find_variable_index( "intra_energy" );
        	int i_NAD_i = (*all_cells)[i]->custom_data.find_variable_index( "intra_NAD" );
        	int i_NADH_i = (*all_cells)[i]->custom_data.find_variable_index( "intra_NADH" );
        	int i_Pyr_i = (*all_cells)[i]->custom_data.find_variable_index( "intra_pyr" );
            
			printf("Celula %d Custom data exit:\n Oxy:%lf Glu:%lf Lac:%lf Ener:%lf Pyr:%lf Nad:%lf Nadh:%lf\n",i , 
				(*all_cells)[i]->custom_data[i_Oxy_i],
				(*all_cells)[i]->custom_data[i_Glu_i],
				(*all_cells)[i]->custom_data[i_Lac_i],
				(*all_cells)[i]->custom_data[energy_vi],
				(*all_cells)[i]->custom_data[i_Pyr_i],
				(*all_cells)[i]->custom_data[i_NAD_i],
				(*all_cells)[i]->custom_data[i_NADH_i]); 
			
			printf("\n");
			*/
        }
    }
    
}


void feed_the_cells(){
	
	double tolerance = 1e-6;  // Use double for floating-point precision
	//printf("Feed the cells:%lf \n", PhysiCell_globals.current_time);
	if (fabs(PhysiCell_globals.current_time - 2880.0) < tolerance || fabs(PhysiCell_globals.current_time - 5760.0) < tolerance || fabs(PhysiCell_globals.current_time - 8640.0) < tolerance || fabs(PhysiCell_globals.current_time - 11520.0) < tolerance) {
    	setup_microenvironment();
	}

}

		
std::vector<std::string> my_coloring_function( Cell* pCell )
{
	
    // Get Energy index
	static int energy_vi = pCell->custom_data.find_variable_index( "intra_energy" );
    
    // start with flow cytometry coloring 
	std::vector<std::string> output = false_cell_coloring_cytometry(pCell); 
	
	// color
    // proliferative cell
	if( pCell->phenotype.death.dead == false && pCell->type == 0 && pCell->custom_data[energy_vi] > 1000)
	{
		output[0] = "rgb(255,255,0)";
		output[2] = "rgb(125,125,0)";
	}

    // arrested cell
	if( pCell->phenotype.death.dead == false && pCell->type == 0 && pCell->custom_data[energy_vi] <= 1000)
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
	}     
    
    // dead cell
	if( pCell->phenotype.death.dead == true && pCell->type == 0)
	{
		output[0] = "rgb(20,20,20)";
		output[2] = "rgb(10,10,10)";
	}  
	
	return output; 
}


std::vector<std::vector<double>> create_cell_circle_positions(double cell_radius, double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*sqrt(3);

	std::vector<double> tempPoint(3,0.0);
	
	for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
	{
		for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
		{
			tempPoint[1]=y + (xc%2) * cell_radius;
			tempPoint[0]=x;
			tempPoint[2]=0;
			if(sqrt(norm_squared(tempPoint))< sphere_radius)
			{ cells.push_back(tempPoint); }
		}
	}
	return cells;
}
std::vector<std::vector<double>> create_cell_circle_positions_3D(double cell_radius, double sphere_radius)
{
    printf("3D");
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*sqrt(3);
    double z_spacing= cell_radius*sqrt(3);

	std::vector<double> tempPoint(3,0.0);
	
	for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
	{
		for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
		{
            for(double z=-sphere_radius;z<sphere_radius;z+=z_spacing, zc++)
		    {
            tempPoint[1]=y + (xc%2) * cell_radius;
			tempPoint[0]=x;
			tempPoint[2]=z + (xc%2)* cell_radius  ;
			if(sqrt(norm_squared(tempPoint))< sphere_radius)
			{ cells.push_back(tempPoint); }
            printf("X:%f\t, Y:%f\t, Z:%f\n", tempPoint[0], tempPoint[1], tempPoint[2]);
		    }
		}
	}
	return cells;
}
std::vector<std::vector<double>> create_cell_square_positions_3D(double cell_radius, double size)
{
    printf("3D");
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*sqrt(3);
    double z_spacing= cell_radius*sqrt(3);

	std::vector<double> tempPoint(3,0.0);
	
	for(double x=-size;x<size;x+=x_spacing, xc++)
	{
		for(double y=-size;y<size;y+=y_spacing, yc++)
		{
            for(double z=-size;z<size;z+=z_spacing, zc++)
		    {
            tempPoint[1]=y ;
			tempPoint[0]=x ;
			tempPoint[2]=z ;
            printf("X:%f\t, Y:%f\t, Z:%f\n", tempPoint[0], tempPoint[1], tempPoint[2]);
			//if(sqrt(norm_squared(tempPoint))< sphere_radius)
			//{ cells.push_back(tempPoint); }
            cells.push_back(tempPoint);
		    }
		}
	}
	return cells;
}

double generateRandomNumberNormal( ){
    #define mean 15
    #define stdDev 0.2

    double u1, u2;
    double z0;


    u1 = (rand() + 1.0) / (RAND_MAX + 2.0);
    u2 = (rand() + 1.0) / (RAND_MAX + 2.0);
    z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * 3.1415926535897932 * u2);

    return mean + (stdDev * z0);
}

double UniformRandom_personal(){
    /*srand(time(NULL));

    int randomInteger = rand() % 100; // Generate a random integer between 0 and 99
    //printf("Random Integer: %d\n", randomInteger);

    double randomFloat = (double)rand() / RAND_MAX; // Generate a random float between 0 and 1
    //printf("Random Float: %f\n", randomFloat);

    return randomFloat;
	*/
    static std::random_device rd; // Create a random device to seed the random number generator
    static std::mt19937 gen(rd()); // Create a Mersenne Twister random number engine
    std::uniform_real_distribution<double> dis(0.0, 1.0); // Create a uniform distribution between 0 and 1

    return dis(gen); // Generate a random number using the engine and distribution
}
double generateNormalRandom(){
    // Set seed based on current time
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);

    // Define mean and standard deviation
    double mean_e = 0.15;
    double stddev = 0.2;

    // Create a normal distribution
    std::normal_distribution<double> distribution(mean_e, stddev);

    // Generate random number following the normal distribution
    return distribution(generator);
}
double random_between_0_and_1() {
    return (double) rand() / RAND_MAX;
}
double locomotive_forces_generator_Ines()
{
    // random number generator to define cell locomotive forces
    // based on anempirically obtained velocity distribution
    double random_value, force_value;
    random_value = UniformRandom();
    force_value =
    0.1157 * pow(random_value, 3) + 0.2422 * pow(random_value, 2) +
    0.0053 * random_value + 0.0048;
    force_value = 13.5*force_value;
    return force_value;
}
double locomotive_forces_generator_Ines_mod()
{
    // random number generator to define cell locomotive forces
    // based on anempirically obtained velocity distribution
    double random_value, force_value;
    random_value = UniformRandom_personal();
	//printf("Random value:%lf\n", random_value);
    force_value =
    0.1157 * pow(random_value, 3) + 0.2422 * pow(random_value, 2) +
    0.0053 * random_value + 0.0048;
    force_value = 13.5*force_value;
    return force_value;
}
double locomotive_forces_generator()
{
    // random number generator to define cell locomotive forces
    // based on anempirically obtained velocity distribution
    double random_value, force_value;
	
    do{
    random_value = generateNormalRandom();
	//random_value = UniformRandom();
    }while(random_value<0);    
    //printf("Random_value:%lf\n", random_value);
    
	//Ines force
	/*
	force_value =
    0.1157 * pow(random_value, 3) + 0.2422 * pow(random_value, 2) +
    0.0053 * random_value + 0.0048;
    force_value = 13.5*force_value;
	*/
    //Force created to fit
    
	force_value =0.1157 * pow(random_value, 3) + 0.2422 * pow(random_value, 2) + 0.0053 * random_value + 0.0048;
    force_value = (920*force_value);

	//Gnuplot fit
	//force_value =0.004987 * pow(random_value, 2) + 0.00049943* random_value + 0.000369;

	//Python fit exponential
	/*
	int a;
	int b;
	a=0.3048498793346672;
	b=-0.2642112540565802;
	force_value =a*exp(b*random_value);
	*/

	//Python fit recta
	
	//double a;
	//double b;
	//double c;
	//a=662.276643627848;
	//b=-263.4036172565524;
	//c=-26.009033062411003;
	//force_value = (a*pow(random_value,2))*(b*random_value)+c;
	

    //printf("Random value:%lf\t m:%lf\t n:%lf\t orce value:%lf\n", random_value, m, n, force_value);
    
	return force_value;
}

/*
double locomotive_forces_generator( )
{
    // random number generator to define cell locomotive forces
    // based on anempirically obtained velocity distribution
    double random_value, force_value;
    random_value = UniformRandom();
	printf("random:%lf\n",random_value);
    force_value =
    0.1157 * pow(random_value, 3) + 0.2422 * pow(random_value, 2) +
    0.0053 * random_value + 0.0048;
    force_value = 13.5*force_value;
    return force_value;
}
*/



void drag_update_cell_velocity_Ines( Cell* pCell, Phenotype& phenotype, double dt ) {
    // sample ECM
    int ECM_density_index = microenvironment.find_density_index( "ECM" );
    double ECM_density = pCell->nearest_density_vector()[ECM_density_index];
    double dyn_viscosity;
    // get viscosity based on concentration
    if(ECM_density == 2.5) {
        //dyn_viscosity = 7.96;
		dyn_viscosity = 7.96;
    }
    else if(ECM_density == 4.0){
        //dyn_viscosity = 18.42;
		dyn_viscosity = 18.42;
    }
    else if(ECM_density == 6.0) {
        //dyn_viscosity = 39.15;
		dyn_viscosity = 39.15;
    }
    // update the speed value
    pCell->phenotype.motility.migration_speed = 
        locomotive_forces_generator_Ines_mod();
    // update velocity
	//printf("Migration speed:%lf\n", pCell->phenotype.motility.migration_speed);
    standard_update_cell_velocity(pCell, phenotype, dt);
    // include the 1/vu (1/ECM density) term to consider friction
	//printf("Vel Antes:%lf\n", pCell->velocity);
    pCell->velocity /= dyn_viscosity;
	//printf("Vel Despues:%lf\n", pCell->velocity);
    return;
}


void drag_update_cell_vel( Cell* pCell, Phenotype& phenotype, double dt ) {
    // sample ECM
	
    int ECM_density_index = microenvironment.find_density_index( "ECM" );
    double ECM_density = pCell->nearest_density_vector()[ECM_density_index];
    double dyn_viscosity;
    // get viscosity based on concentration of the collagen of the extracellular matrix
    if(ECM_density == 2.5) {
        //dyn_viscosity = 7.96;
		dyn_viscosity = 7.96/60.0;
    }
    else if(ECM_density == 4.0){
        //dyn_viscosity = 18.42;
		dyn_viscosity = 18.42/60.0;
    }
    else if(ECM_density == 6.0) {
        //dyn_viscosity = 39.15;
		dyn_viscosity = 39.15/60.0;
    }    

    // calculate the velocity generated by the random locomotive force that drives the cell. Store it in the migration speed of the cell
    pCell->phenotype.motility.migration_speed = 
        locomotive_forces_generator();
    //printf("Vel Antes:%lf\n", pCell->phenotype.motility.migration_speed);

	// update velocity with the standard model. The drag velocity still needs to be added.
    standard_update_cell_velocity(pCell, phenotype, dt);	

	//Add the drag velocity
	double viscosity_factor;//This will be the factor needed to add the drag  
	viscosity_factor= 6*MPI*phenotype.geometry.radius*dyn_viscosity;


	
	//printf("Vel Antes:%lf\n", pCell->velocity);
    pCell->velocity/=(viscosity_factor);
	//pCell->velocity /= dyn_viscosity;
	//printf("Vel Despues:%lf\n", pCell->velocity);
	
    return;
}
