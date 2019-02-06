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

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include <iostream>
#include <vector>
#include <string>

#include "../CoDiPack/include/codi.hpp"

using namespace std;

typedef codi::RealReverse addouble;
typedef double passivedouble;

typedef Eigen::Matrix<addouble, Eigen::Dynamic, Eigen::Dynamic> MatrixXdDiff; // MatrixXd
typedef Eigen::Matrix<addouble, Eigen::Dynamic, 1> VectorXdDiff;       // VectorXd
typedef Eigen::Matrix<addouble, 3, 3> Matrix3dDiff;             // Matrix3d
typedef Eigen::Matrix<addouble, 3, 1> Vector3dDiff;             // Vector3d

//const unsigned int MAX_STRING_SIZE = 200;    /*!< \brief Maximum number of domains. */

class COptionBase {
private:
public:
  COptionBase() {};
  virtual  ~COptionBase() = 0;
  //  virtual string SetValue(string) {SU2MPI::PrintAndFinalize("shouldn't be here"); return "";};
  virtual string SetValue(vector<string>) = 0;
  virtual void SetDefault() = 0;

  string optionCheckMultipleValues(vector<string> & option_value, string type_id, string option_name) {
    if (option_value.size() != 1) {
      string newString;
      newString.append(option_name);
      newString.append(": multiple values for type ");
      newString.append(type_id);
      return newString;
    }
    return "";
  }
  
  string badValue(vector<string> & option_value, string type_id, string option_name) {
    string newString;
    newString.append(option_name);
    newString.append(": improper option value for type ");
    newString.append(type_id);
    return newString;
  }
};  

class COptionDouble : public COptionBase {
  addouble & field; // Reference to the fieldname
  addouble def; // Default value
  string name; // identifier for the option

public:
  COptionDouble(string option_field_name, addouble & option_field, addouble default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionDouble() {};
  string SetValue(vector<string> option_value) {
    // check if there is more than one value
    string out = optionCheckMultipleValues(option_value, "addouble", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    istringstream is(option_value[0]);
    addouble val;
    if (is >> val) {
      this->field = val;
      return "";
    }
    return badValue(option_value, "addouble", this->name);
  }
  void SetDefault() {
    this->field = this->def;
  }
};

class COptionString : public COptionBase {
  string & field; // Reference to the fieldname
  string def; // Default value
  string name; // identifier for the option

public:
  COptionString(string option_field_name, string & option_field, string default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionString() {};
  string SetValue(vector<string> option_value) {
    // check if there is more than one value
    string out = optionCheckMultipleValues(option_value, "addouble", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    this->field.assign(option_value[0]);
    return "";
  }
  void SetDefault() {
    this->field = this->def;
  }
};

class COptionULong : public COptionBase {
  unsigned long & field; // Reference to the feildname
  unsigned long def; // Default value
  string name; // identifier for the option

public:
  COptionULong(string option_field_name, unsigned long & option_field, unsigned long default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionULong() {};
  string SetValue(vector<string> option_value) {
    string out = optionCheckMultipleValues(option_value, "unsigned long", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    istringstream is(option_value[0]);
    unsigned long val;
    if (is >> val) {
      this->field = val;
      return "";
    }
    return badValue(option_value, "unsigned long", this->name);
  }
  void SetDefault() {
    this->field = this->def;
  }
};

class COptionUShort : public COptionBase {
  unsigned short & field; // Reference to the feildname
  unsigned short def; // Default value
  string name; // identifier for the option

public:
  COptionUShort(string option_field_name, unsigned short & option_field, unsigned short default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionUShort() {};
  string SetValue(vector<string> option_value) {
    string out = optionCheckMultipleValues(option_value, "unsigned short", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    istringstream is(option_value[0]);
    unsigned short val;
    if (is >> val) {
      this->field = val;
      return "";
    }
    return badValue(option_value, "unsigned short", this->name);
  }
  void SetDefault() {
    this->field = this->def;
  }
};

class COptionInt : public COptionBase {
  int & field; // Reference to the feildname
  int def; // Default value
  string name; // identifier for the option

public:
  COptionInt(string option_field_name, int & option_field, int default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionInt() {};
  string SetValue(vector<string> option_value) {
    string out = optionCheckMultipleValues(option_value, "int", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    istringstream is(option_value[0]);
    int val;
    if (is >> val) {
      this->field = val;
      return "";
    }
    return badValue(option_value, "int", this->name);
  }
  void SetDefault() {
    this->field = this->def;
  }
};

/*!
 * \brief utility function for converting strings to uppercase
 * \param[in, out] str - string we want to convert
 */
inline void StringToUpperCase(string & str) {
  std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}