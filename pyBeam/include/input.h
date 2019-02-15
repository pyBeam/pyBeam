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

#pragma once
#include <math.h>
#include <map>
#include <stdio.h>      /* printf, fopen */
/*We can think of using #include "./mpi_structure.hpp" and the command SU2_MPI::Error
 * once fully integrated with SU2 core
 */ 
#include <stdlib.h>     /* exit, EXIT_FAILURE */ 

#include "../include/types.h"





class CInput
{

private:
	
protected:

	//##################     Numerical Inputs     ###########################

	unsigned long nNodes;	// Number of overall nodes along the wing (no collapsed)
	unsigned long nFEM;		// Number of finite elements

	unsigned short nDOF; 	// Number of degrees of freedom
	
	addouble load; 			// [N];
	int follower_flag;		// (0) Nonfollower (1) follower (2) approx follower
	unsigned long loadSteps;			// Number of load steps
	unsigned long nIter;				// Number of iterations

	addouble end_time;		// [sec] for SS calculation
	addouble dt;   			// [sec] time increment for SS calculation

	//##############    Wing Inputs  ###########################
	// Units Sys: SI

	addouble t; 				// web & flange thickness [m]
	addouble h;				// web height [m]
	addouble b;				// flange width [m]
	addouble E; 				// Elastic modulus [GPa]
	addouble Poiss; 			// Poisson Ratio
	addouble ro;				// Beam Density [kg/m^3]
	addouble G;				// Shear modulus
	addouble l; 				// Wing Length [m]
	addouble A;				// cross section area
	addouble As_z; 			// z Effective shear area
	addouble As_y;			// y Effective shear area

	addouble Iyy, Izz;
	addouble Jx; 				//Polar Moment of Inertia


	addouble Mwing;			//Wing's mass [kg]
	addouble EIy, EIz, GJ, AE;

	addouble Clalpha;  		//  recall pi = atan(1)*4;
	addouble Cldelta;

	//#################    Elements properties    ############################

	addouble le;      	//element length
	addouble m, m_e; 		//Element's mass
	addouble m_w, m_f; 	//web and flange mass
	addouble Ix, Iz; 		//[kg*m^2]

	//################     Convergence Parameters    ###########################

	addouble convCriteria;

	
public:

  CInput(void);
  
  virtual ~CInput(void);
  
  void SetWebThickness(passivedouble thickness) {t = thickness; }

  void SetWebHeight(passivedouble height) {h = height; }  
  
  void SetFlangeWidth(passivedouble FlangeWidth) {b = FlangeWidth; }   
 
  void SetYoungModulus(passivedouble YoungModulus) {E = YoungModulus; } 

  void SetPoisson(passivedouble Poisson) {Poiss = Poisson; }  
  
  void SetDensity(passivedouble Density) {ro = Density; }
    
  void SetBeamLength(passivedouble BeamLength) {l = BeamLength; } 
  
  void SetLoad(passivedouble Load) {load = Load; }  
  
  void SetFollowerFlag(int FollowerFlag) {follower_flag = FollowerFlag; }  
  
  void SetLoadSteps(unsigned long LoadSteps) { loadSteps = LoadSteps; }   
  
  void SetNStructIter(unsigned long NStructIter) {nIter = NStructIter; }   
  
  void SetConvCriterium(passivedouble ConvCriterium) {convCriteria = ConvCriterium; }   


  void SetParameters();
  
  unsigned long Get_nNodes(void) { return nNodes; }  
    
  unsigned long Get_nFEM(void) { return nFEM; }
  
  unsigned short Get_nDOF(void) { return nDOF; }
  
  unsigned short Get_FollowerFlag(void) { return follower_flag; }  
  
  unsigned long Get_LoadSteps(void) { return loadSteps;}
  
  unsigned long Get_nIter(void) { return nIter;}  
  
  addouble Get_Thickness(void) { return t;}
  
  addouble Get_Load(void) { return load;}  
  
  addouble Get_l(void) { return l; }  
  
  addouble Get_le(void) { return le; }
  
  addouble Get_Jx(void) { return Jx; }
  
  addouble Get_m_e(void) { return m_e; } 
  
  addouble Get_A(void) { return A; } 
      
  addouble Get_EIz(void) { return EIz; } 
        
  addouble Get_EIy(void) { return EIy; } 
  
  addouble Get_GJ(void) { return GJ; } 
  
  addouble Get_AE(void) { return AE; } 
  
  addouble Get_m(void) { return m; }
  
  addouble Get_Iyy(void) { return Iyy; }     
          
  addouble Get_Izz(void) { return Izz; }  
    

  addouble Get_ConvCriteria(void) { return convCriteria; }  

    

  
    
  /* Please see config_structure.hpp in common/include of SU2 suite for more info about the fuction
   * */

  /*--- all_options is a map containing all of the options. This is used during config file parsing
   to track the options which have not been set (so the default values can be used). Without this map
   there would be no list of all the config file options. ---*/
  
  map<string, bool> all_options;
  
  /*--- brief param is a map from the option name (config file string) to its decoder (the specific child
   class of COptionBase that turns the string into a value) ---*/
  
  map<string, COptionBase*> option_map;   
     
  // All of the addXxxOptions take in the name of the option, and a refernce to the field of that option
  // in the option structure. Depending on the specific type, it may take in a default value, and may
  // take in extra options. The addXxxOptions mostly follow the same pattern, so please see addDoubleOption
  // for detailed comments.     
     
  /*!<\brief addDoubleOption creates a config file parser for an option with the given name whose
   value can be represented by a addouble.*/    
    
  void addDoubleOption(const string name, addouble & option_field, addouble default_value) {
    // Check if the key is already in the map. If this fails, it is coder error
    // and not user error, so throw.
    assert(option_map.find(name) == option_map.end());
    
    // Add this option to the list of all the options
    all_options.insert(pair<string, bool>(name, true));
    
    // Create the parser for a su2double option with a reference to the option_field and the desired
    // default value. This will take the string in the config file, convert it to a su2double, and
    // place that su2double in the memory location specified by the reference.
    COptionBase* val = new COptionDouble(name, option_field, default_value);
    
    // Create an association between the option name ("CFL") and the parser generated above.
    // During configuration, the parsing script will get the option name, and use this map
    // to find how to parse that option.
    option_map.insert(pair<string, COptionBase *>(name, val));
  }
  
  void addStringOption(const string name, string & option_field, string default_value) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionString(name, option_field, default_value);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }  
  
  void addUnsignedLongOption(const string name, unsigned long & option_field, unsigned long default_value) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionULong(name, option_field, default_value);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }  
  
  void addUnsignedShortOption(const string name, unsigned short & option_field, unsigned short default_value) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionUShort(name, option_field, default_value);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }  
  
  void addIntegerOption(const string name, int & option_field, int default_value) {
    assert(option_map.find(name) == option_map.end());
    all_options.insert(pair<string, bool>(name, true));
    COptionBase* val = new COptionInt(name, option_field, default_value);
    option_map.insert(pair<string, COptionBase *>(name, val));
  }  
    
};


