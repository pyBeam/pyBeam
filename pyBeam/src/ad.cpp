/*
 * pyBeam, a Beam Solver
 *
 * Copyright (C) 2018 Tim Albring, Ruben Sanchez, Rauno Cavallaro
 *
 * Developers: Tim Albring, Ruben Sanchez (SciComp, TU Kaiserslautern)
 *             Rauno Cavallaro (Carlos III University Madrid)
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



#ifdef CODI_REVERSE_TYPE
#include "../include/datatypes/ad_reverse.hpp"
#elif CODI_FORWARD_TYPE
#include "../include/datatypes/ad_forward.hpp"
#else
#include "../include/datatypes/ad_passive.hpp"
#endif

namespace AD {
#ifdef CODI_REVERSE_TYPE

  /*--- Initialization of the tape ---*/
  addouble::TapeType& globalTape = addouble::getGlobalTape();

#endif
}
