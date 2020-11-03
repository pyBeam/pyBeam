#!/usr/bin/env python3
#
# pyBeam, a Beam Solver
#
# Copyright (C) 2018 Ruben Sanchez, Rocco Bombardieri, Rauno Cavallaro
#
# Developers: Ruben Sanchez (SciComp, TU Kaiserslautern)
#             Rocco Bombardieri, Rauno Cavallaro (Carlos III University Madrid)
#
# This file is part of pyBeam.
#
# pyBeam is free software: you can redistribute it and/or
# modify it under the terms of the GNU Affero General Public License
# as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# pyBeam is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# See the GNU Affero General Public License for more details.
# You should have received a copy of the GNU Affero
# General Public License along with pyBeam.
# If not, see <http://www.gnu.org/licenses/>.
#


import numpy as np
import os,os.path
import openmdao.api as om
import sys
from pyBeamLib import pyBeamSolver
from pyBeamOpenMDAO import pyBeamOpt
from pprint import pprint
import sqlitedict
import matplotlib.pyplot as plt



class Weight(om.ExplicitComponent):

    def setup(self):
        self.add_input('C_wb',units='mm')
        self.add_input('h',units='mm')
        self.add_input('t_sk',units='mm')
        self.add_input('t_sp',units='mm')
        self.add_input('A_fl',units='mm**2')
        self.add_input('n_stiff',units='mm')
        self.add_input('A_stiff',units='mm**2')


        self.add_output('Obj_f')

        # Finite difference all partials.
        self.declare_partials(of='Obj_f', wrt='C_wb', method='fd',form='central')
        self.declare_partials(of='Obj_f', wrt='h', method='fd',form='central')
        self.declare_partials(of='Obj_f', wrt='t_sk', method='fd',form='central')
        self.declare_partials(of='Obj_f', wrt='t_sp', method='fd',form='central')
        self.declare_partials(of='Obj_f', wrt='A_fl', method='fd',form='central')
        self.declare_partials(of='Obj_f', wrt='A_stiff', method='fd',form='central')

    def compute(self, inputs, outputs):


        C_wb = inputs['C_wb']
        h    = inputs['h']
        t_sk = inputs['t_sk']
        t_sp = inputs['t_sp']
        A_fl = inputs['A_fl']
        n_stiff = inputs['n_stiff']
        A_stiff = inputs['A_stiff']


        outputs['Obj_f'] = beam_opt.ComputeWeight_opt( (float(C_wb), float(h), float(t_sk), float(t_sp), float(A_fl),  int(n_stiff), float(A_stiff)))





class KS_constraint(om.ExplicitComponent):

    def setup(self):
        self.add_input('C_wb',units='mm')
        self.add_input('h',units='mm')
        self.add_input('t_sk',units='mm')
        self.add_input('t_sp',units='mm')
        self.add_input('A_fl',units='mm**2')
        self.add_input('n_stiff',units='mm')
        self.add_input('A_stiff',units='mm**2')


        self.add_output('Const_KS')

        # Finite difference all partials.
        self.declare_partials(of='Const_KS', wrt='C_wb', method='fd',form='central')
        self.declare_partials(of='Const_KS', wrt='h', method='fd',form='central')
        self.declare_partials(of='Const_KS', wrt='t_sk', method='fd',form='central')
        self.declare_partials(of='Const_KS', wrt='t_sp', method='fd',form='central')
        self.declare_partials(of='Const_KS', wrt='A_fl', method='fd',form='central')
        self.declare_partials(of='Const_KS', wrt='A_stiff', method='fd',form='central')


    def compute(self, inputs, outputs):


        C_wb = inputs['C_wb']
        h    = inputs['h']
        t_sk = inputs['t_sk']
        t_sp = inputs['t_sp']
        A_fl = inputs['A_fl']
        n_stiff = inputs['n_stiff']
        A_stiff = inputs['A_stiff']



        outputs['Const_KS'] =   beam_opt.ComputeResponseKSStress_opt_lin( (float(C_wb), float(h), float(t_sk), float(t_sp), float(A_fl), int(n_stiff), float(A_stiff)))

        """NonLinear Analysis"""
        #outputs['Const_KS'] = beam_opt.ComputeResponseKSStress_opt_non_lin( (float(C_wb), float(h), float(t_sk), float(t_sp), float(A_fl), int(n_stiff), float(A_stiff)))








if __name__ == "__main__":

    prob = om.Problem()

    Loads = (60, 150, 5000, 5000)
    file_dir = os.path.dirname(os.path.realpath(__file__))
    #beam = pyBeamSolver(file_dir, 'config_lin.cfg')
    beam_opt = pyBeamOpt(file_dir, 'config_lin.cfg', Loads)

    """NonLinear Analysis"""
    #beam = pyBeamSolver(file_dir, 'config_non_lin.cfg')
    #beam_opt = pyBeamOpt(file_dir, 'config_non_lin.cfg', Loads)


    #DVs = beam.GetDesignVariables()

    C_wb = 3000
    h = 500
    t_sk = 3
    t_sp = 5
    A_fl = 200
    n_stiff = 6
    A_stiff = 50

    prob.model.add_subsystem('weight_comp', Weight(),
                             promotes_inputs=['C_wb', 'h', 't_sk', 't_sp', 'A_fl','n_stiff','A_stiff'])

    # define the component whose output will be constrained
    prob.model.add_subsystem('KS_comp', KS_constraint(),
                             promotes_inputs=['C_wb', 'h', 't_sk', 't_sp', 'A_fl','n_stiff','A_stiff'])



    # setup the optimization

    #prob.driver = om.pyOptSparseDriver()
    prob.driver = om.ScipyOptimizeDriver()
    prob.driver.options['optimizer'] = 'SLSQP'



    #prob.driver.options['maxiter'] = 100
    prob.driver.options['debug_print']=['desvars','objs','nl_cons','totals']



    prob.model.add_design_var('C_wb',lower=1)
    prob.model.add_design_var('h',lower=1)
    prob.model.add_design_var('t_sk',lower=1)
    prob.model.add_design_var('t_sp',lower=1)
    prob.model.add_design_var('A_fl',lower=1)
    prob.model.add_design_var('A_stiff',lower = 1)


    filename = "opt_results.sql"
    recorder = om.SqliteRecorder(filename)
    prob.driver.add_recorder(recorder)
    prob.driver.recording_options['record_desvars'] = True
    prob.driver.recording_options['record_objectives'] = True
    prob.driver.recording_options['record_constraints'] = True
    prob.driver.recording_options['record_derivatives'] = True
    prob.driver.recording_options['includes'] = []
    prob.driver.recording_options['excludes'] = []

    prob.model.add_objective('weight_comp.Obj_f')
    prob.model.add_constraint('KS_comp.Const_KS',upper=0)

    prob.setup()

    prob.set_val('C_wb', C_wb,units='mm')
    prob.set_val('h', h,units='mm')
    prob.set_val('t_sk', t_sk,units='mm')
    prob.set_val('t_sp', t_sp,units='mm')
    prob.set_val('A_fl', A_fl,units='mm**2')
    prob.set_val('n_stiff', n_stiff,units='mm')
    prob.set_val('A_stiff', A_stiff,units='mm**2')

    prob.run_driver()

    print("Minimum Weight =", prob.get_val('weight_comp.Obj_f'))
    print("Constraint =", prob.get_val('KS_comp.Const_KS'))

    prob.cleanup()


    """  Record Data : """


    cr = om.CaseReader("opt_results.sql")

    case_names = cr.list_cases(out_stream=None)

    path = os.path.abspath(os.path.join(file_dir, "../Optimization"))
    final_directory = os.path.join(path, r'opt_iter')
    if not os.path.exists(final_directory):
        os.makedirs(final_directory)



    path_iter = os.path.abspath(os.path.join(file_dir, "../Optimization/opt_iter"))






    for it in range(0,len(case_names),1):

            case = cr.get_case(it)
            derivs=cr.get_case(it).derivatives

            final_directory_iter = os.path.join(path_iter, r'000'+ str(it))
            if not os.path.exists(final_directory_iter):
                os.makedirs(final_directory_iter)

            if derivs is not None :
             Weight_deriv = os.path.join(final_directory_iter, "Weight_sens" + ".txt")

             Deriv_obj = open(Weight_deriv, "w")
             Deriv_obj.write('Weight wrt C_wb = ' + str(derivs['weight_comp.Obj_f','C_wb']) + '\n'
                             'Weight wrt h = ' + str(derivs['weight_comp.Obj_f','h']) + '\n'
                             'Weight wrt t_sk = ' + str(derivs['weight_comp.Obj_f','t_sk']) + '\n'
                             'Weight wrt t_sp = ' + str(derivs['weight_comp.Obj_f','t_sp']) + '\n'
                             'Weight wrt A_fl = ' + str(derivs['weight_comp.Obj_f','A_fl']) + '\n'
                             'Weight wrt A_stiff = ' + str(derivs['weight_comp.Obj_f','A_stiff']))
             Deriv_obj.close()

             KS_deriv = os.path.join(final_directory_iter, "KS_sens" + ".txt")

             Deriv_KS = open(KS_deriv, "w")
             Deriv_KS.write('KS wrt C_wb = ' + str(derivs['KS_comp.Const_KS', 'C_wb']) + '\n'
                                                                                               'KS wrt h = ' + str(
                 derivs['KS_comp.Const_KS', 'h']) + '\n'
                        'KS wrt t_sk = ' + str(
                 derivs['KS_comp.Const_KS', 't_sk']) + '\n'
                        'KS wrt t_sp = ' + str(
                 derivs['KS_comp.Const_KS', 't_sp']) + '\n'
                        'KS wrt A_fl = ' + str(derivs['KS_comp.Const_KS', 'A_fl']) + '\n'
             'KS wrt A_stiff = ' + str(derivs['KS_comp.Const_KS', 'A_stiff']))
             Deriv_KS.close()

            if derivs is None:
                pass

            Obj_func = os.path.join(final_directory_iter,  "Weight" + ".txt")

            outputs_obj = open(Obj_func , "w")
            outputs_obj.write(str(case['weight_comp.Obj_f']))
            outputs_obj.close()

            KS = os.path.join(final_directory_iter, "KS_constraint" + ".txt")

            outputs_KS = open(KS, "w")
            outputs_KS.write(str(case['KS_comp.Const_KS']))
            outputs_KS.close()

            Box_width = os.path.join(final_directory_iter, "C_wb" + ".txt")

            outputs_C_wb = open(Box_width, "w")
            outputs_C_wb.write(str(case['C_wb']))
            outputs_C_wb.close()

            height= os.path.join(final_directory_iter, "h" + ".txt")

            outputs_h = open(height, "w")
            outputs_h.write(str(case['h']))
            outputs_h.close()

            t_skin = os.path.join(final_directory_iter, "t_sk" + ".txt")

            outputs_t_sk = open(t_skin, "w")
            outputs_t_sk.write(str(case['t_sk']))
            outputs_t_sk.close()

            t_spar= os.path.join(final_directory_iter, "t_sp" + ".txt")

            outputs_t_sp = open(t_spar, "w")
            outputs_t_sp.write(str(case['t_sp']))
            outputs_t_sp.close()


            A_flanges = os.path.join(final_directory_iter, "A_fl" + ".txt")

            outputs_A_fl = open(A_flanges, "w")
            outputs_A_fl.write(str(case['A_fl']))
            outputs_A_fl.close()

            A_stiff = os.path.join(final_directory_iter, "A_stiff" + ".txt")

            outputs_A_stiff = open(A_stiff, "w")
            outputs_A_stiff.write(str(case['A_stiff']))
            outputs_A_stiff.close()







    plt.figure()
    for i in range(0, len(case_names), 1):
        case_plot = cr.get_case(i)
        plt.plot(i, float(case_plot['weight_comp.Obj_f']), 'ro', label='Weight')
    plt.legend(['Weight'])
    plt.xlabel('Driver Iterations', fontsize=18)
    plt.ylabel('Weight', fontsize=16)


    plt.figure()
    for i in range(0, len(case_names), 1):
        case_plot = cr.get_case(i)
        plt.plot(i, float(case_plot['KS_comp.Const_KS']), 'ro')
    plt.legend(['KS_constraint'])
    plt.xlabel('Driver Iterations', fontsize=18)
    plt.ylabel('KS Constraint', fontsize=16)

    plt.figure()
    for i in range(0, len(case_names), 1):
        case_plot = cr.get_case(i)

        plt.plot(i, float(case_plot['C_wb']), 'ro', label='C_wb')
        plt.plot(i, float(case_plot['h']), 'bs', label='h')
        plt.plot(i, float(case_plot['t_sk']), 'g^', label='t_sk')
        plt.plot(i, float(case_plot['t_sp']), 'yo', label='t_sp')
        plt.plot(i, float(case_plot['A_fl']), 'k+', label='A_fl')
        plt.plot(i, float(case_plot['A_stiff']), 'co', label='A_stiff')
    plt.legend(['C_wb', 'h', 't_sk', 't_sp', 'A_fl','A_stiff'], bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xlabel('Driver Iterations', fontsize=18)
    plt.ylabel('Design Variables', fontsize=16)

    plt.show()