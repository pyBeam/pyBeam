/*
 * pyBeam, a Beam Solver
 *
 * Copyright (C) 2018 Tim Albring, Ruben Sanchez, Rauno Cavallaro, Rocco Bombardieri
 * 
 * Developers: Tim Albring, Ruben Sanchez (SciComp, TU Kaiserslautern)
 *             Rauno Cavallaro, Rocco Bombardieri (Carlos III University Madrid)
 *
 * This file is part of pyBeam.
 *
 * pyBeam is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License
 * as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * pyBeam is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU Affero General Public License for more details.
 * You should have received a copy of the GNU Affero
 * General Public License along with pyBeam.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <iostream>
#include <fstream>
#include <chrono>


#include "../include/input.h"

using namespace std;

CInput::CInput(void) {
  
}

void CInput::SetParameters(case_filename){
	
    /*--- Parsing the config file  ---*/

    SetConfig_Parsing(case_filename);	
	
	//##################     Numerical Inputs     ###########################

	//nNodes = 101; 			
	// number of overall nodes along the wing (no collapsed)
	addUnsignedLongOption("N_NODES", nNodes, 0);  
	// TO ADD WARNING FOR DEFAULT VALUE	

	//nDOF = 6;                
	// number of rigid modes to be calculated
	addUnsignedShortOption("N_DOF", nDOF, 0);
	// TO ADD WARNING FOR DEFAULT VALUE
		
	//load = 5000; 			
	// [N];
	addDoubleOption("LOAD", load, 0);
		
	//follower_flag = 0;		
	// (0) Nonfollower (1) follower (2) approx follower
	addDoubleOption("FOLLOWER_FLAG", follower_flag, "NON_FOLLOWER");
	
	
	//loadSteps = 1;			
	// Number of load steps
	addUnsignedLongOption("LOAD_STEPS", loadSteps, 1);
	
	
	//nIter = 30;			
	// Number of iterations 
    addUnsignedLongOption("N_STRUCT_ITER", nIter, 1);
    
    //t = thickness;		
    // web & flange thickness [m]
    addDoubleOption("W_THICKNESS", t, 0);
    // TO ADD WARNING FOR DEFAULT VALUE
    
	//h = 40*1e-2;			
	// web height [m]    
    addDoubleOption("W_HEIGHT", h, 0);
    // TO ADD WARNING FOR DEFAULT VALUE
    
	//b = 20*1e-2;			
	// flange width [m]    
    addDoubleOption("F_WIDTH", b, 0);
    // TO ADD WARNING FOR DEFAULT VALUE
    
	//E = 70*1e9; 			
	// Elastic modulus [GPa]    
    addDoubleOption("Y_MODULUS", E, 0);
    // TO ADD WARNING FOR DEFAULT VALUE
    
	//Poiss = 0.3; 			
	// Poisson Ratio    
    addDoubleOption("POISSON", Poiss, 0);
    // TO ADD WARNING FOR DEFAULT VALUE
    
	//ro = 2.7e3;				
	// Beam Density [kg/m^3]
    addDoubleOption("RHO", ro, 0);
    // TO ADD WARNING FOR DEFAULT VALUE    
    
	//l = 30; 			
	//// Wing Length [m] 
	addDoubleOption("B_LENGTH", l, 0);	   
    // TO ADD WARNING FOR DEFAULT VALUE   
    	
	//################     Convergence Parameters    ###########################

	//convCriteria = 1e-4;	
	addDoubleOption("CONV_CRITERIUM", convCriteria, 1e-4);	   
	
    
    nFEM = nNodes - 1;
	//##############    Wing Inputs  ###########################
	// Units Sys: SI

	G = E/(2*(1+Poiss) );	// Shear modulus

	A = t*h+2*t*b;			// cross section area
	As_z = t*h; 			// z Effective shear area
	As_y = 5/6.0*2*t*b;		// y Effective shear area

	Iyy = (pow(t,3)*h)/12.0 + 2.0* (  b*pow(t,3)/12.0 + b*t*pow( (h+t)/2.0 , 2)    );  //cross section Izz [m^4]
	Izz = ( h*pow(t,3) + 2*t*pow(b,3))/12.0;
	Jx = Iyy+Izz; 			//Polar Moment of Inertia


	Mwing = A*l*ro;		//Wing's mass [kg]
	EIy = E*Iyy;
	EIz = E*Izz;
	GJ = G*Jx;
	AE = A*E;

	Clalpha = 2*   atan(1)*4;  //  recall pi = atan(1)*4;
	Cldelta = 1;


	//#################    Elements properties    ############################

	le = l/(nNodes-1);      	//element length
	m_e = Mwing/(nNodes-1); 	//Element mass
	m = m_e; 				//Element's mass

	m_w = ro*le*t*h; 		//web mass
	m_f = ro*le*t*b; 		//flange mass
	
	Iz = ( m_w*(pow(t,2)+pow(le,2)) +
		 2*m_f*(pow(b,2)+pow(le,2)) )/12.0;     //[kg*m^2]
		 
	Ix = m_w/12.0*(pow(h,2)+pow(t,2))+
		 2*m_f*( (pow(t,2)+pow(b,2))/12+
		 pow((t/2+h/2),2)) ; 		//[kg*m^2]

    
}

void CInput::SetConfig_Parsing(char case_filename) {
  string text_line, option_name;
  ifstream case_file;
  vector<string> option_value;
  
  /*--- Read the configuration file ---*/
  
  case_file.open(case_filename, ios::in);

  if (case_file.fail()) {
	  
    printf("The beam configuration file (.cfg) is missing!!");
    exit (EXIT_FAILURE);
  }

  string errorString;

  int  err_count = 0;  // How many errors have we found in the config file
  int max_err_count = 30; // Maximum number of errors to print before stopping

  map<string, bool> included_options;

  /*--- Parse the configuration file and set the options ---*/
  
  while (getline (case_file, text_line)) {
    
    if (err_count >= max_err_count) {
      errorString.append("too many errors. Stopping parse");

      cout << errorString << endl;
      throw(1);
    }
    
    if (TokenizeString(text_line, option_name, option_value)) {
      
      /*--- See if it's a python option ---*/

      /*if (option_map.find(option_name) == option_map.end()) {
          string newString;
          newString.append(option_name);
          newString.append(": invalid option name");
          newString.append(". Check current SU2 options in config_template.cfg.");
          newString.append("\n");
          //if (!option_name.compare("AD_COEFF_FLOW")) newString.append("AD_COEFF_FLOW= (1st, 2nd, 4th) is now JST_SENSOR_COEFF= (2nd, 4th).\n");
          nue;
      }*/

      /*--- Option exists, check if the option has already been in the config file ---*/
      
      if (included_options.find(option_name) != included_options.end()) {
        string newString;
        newString.append(option_name);
        newString.append(": option appears twice");
        newString.append("\n");
        errorString.append(newString);
        err_count++;
        continue;
      }


      /*--- New found option. Add it to the map, and delete from all options ---*/
      
      included_options.insert(pair<string, bool>(option_name, true));
      all_options.erase(option_name);

      /*--- Set the value and check error ---*/
      
      string out = option_map[option_name]->SetValue(option_value);
      if (out.compare("") != 0) {
        errorString.append(out);
        errorString.append("\n");
        err_count++;
      }
    }
  }

  /*--- See if there were any errors parsing the config file ---*/
      
  if (errorString.size() != 0) {
    printf(errorString);
    exit (EXIT_FAILURE);
  }

  /*--- Set the default values for all of the options that weren't set ---*/
      
  for (map<string, bool>::iterator iter = all_options.begin(); iter != all_options.end(); ++iter) {
    option_map[iter->first]->SetDefault();
  }

  case_file.close();
  
}



CInput::~CInput(void) {
	
}


