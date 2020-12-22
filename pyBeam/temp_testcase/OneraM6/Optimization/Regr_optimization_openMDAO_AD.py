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
import shutil
import openmdao.api as om
import sys
from PIL import Image
from pyBeamLibAD import pyBeamSolverAD
from pyBeamOpenMDAO import pyBeamOpt
from pprint import pprint
import sqlitedict
import matplotlib.pyplot as plt


class Weight(om.ExplicitComponent):

    def setup(self):

        self.add_input('C_wb',shape = len(C_wb),units='mm')
        self.add_input('h', shape = len(h),units='mm')
        self.add_input('t_sk',shape = len(t_sk),units='mm')
        self.add_input('t_sp',shape = len(t_sp),units='mm')
        self.add_input('A_fl',shape = len(A_fl),units='mm**2')
        self.add_input('n_stiff',shape = len(n_stiff),units='mm')
        self.add_input('A_stiff',shape = len(A_stiff),units='mm**2')



        self.add_output('Obj_f')

        # Finite difference all partials.
        self.declare_partials('Obj_f', 'C_wb')
        self.declare_partials('Obj_f', 'h')
        self.declare_partials('Obj_f', 't_sk')
        self.declare_partials('Obj_f', 't_sp')
        self.declare_partials('Obj_f', 'A_fl')
        self.declare_partials('Obj_f', 'A_stiff')

    def compute(self, inputs, outputs):


        C_wb = inputs['C_wb']
        h    = inputs['h']
        t_sk = inputs['t_sk']
        t_sp = inputs['t_sp']
        A_fl = inputs['A_fl']
        n_stiff = inputs['n_stiff']
        A_stiff = inputs['A_stiff']

        DesignVArs = np.array([C_wb, h, t_sk, t_sp, A_fl, n_stiff, A_stiff]).transpose()
        DesignVArs = np.array(DesignVArs.tolist())
        DesignVArs = np.reshape(DesignVArs, DesignVArs.size)




        outputs['Obj_f'] = beam_opt.ComputeWeight_opt(DesignVArs)



    def compute_partials(self,inputs,partials):

        C_wb = inputs['C_wb']
        h = inputs['h']
        t_sk = inputs['t_sk']
        t_sp = inputs['t_sp']
        A_fl = inputs['A_fl']
        n_stiff = inputs['n_stiff']
        A_stiff = inputs['A_stiff']

        DesignVArs = np.array([C_wb, h, t_sk, t_sp, A_fl, n_stiff, A_stiff]).transpose()
        DesignVArs = np.array(DesignVArs.tolist())
        DesignVArs = np.reshape(DesignVArs, DesignVArs.size)

        Weight_adj=beam_opt.ComputeAdjointWeight_opt(DesignVArs)

        partials['Obj_f', 'C_wb']=Weight_adj[::6]
        partials['Obj_f', 'h'] = Weight_adj[1::6]
        partials['Obj_f', 't_sk'] = Weight_adj[2::6]
        partials['Obj_f', 't_sp'] = Weight_adj[3::6]
        partials['Obj_f', 'A_fl'] = Weight_adj[4::6]
        partials['Obj_f', 'A_stiff'] = Weight_adj[5::6]




class KSStress_constraint(om.ExplicitComponent):

    def setup(self):
        self.add_input('C_wb',shape=len(C_wb),units='mm')
        self.add_input('h',shape=len(h),units='mm')
        self.add_input('t_sk',shape=len(t_sk),units='mm')
        self.add_input('t_sp',shape=len(t_sp),units='mm')
        self.add_input('A_fl',shape=len(A_fl),units='mm**2')
        self.add_input('n_stiff',shape=len(n_stiff),units='mm')
        self.add_input('A_stiff',shape=len(A_stiff),units='mm**2')


        self.add_output('Const_KS')

        # Finite difference all partials.

        self.declare_partials('Const_KS', 'C_wb')
        self.declare_partials('Const_KS', 'h')
        self.declare_partials('Const_KS', 't_sk')
        self.declare_partials('Const_KS', 't_sp')
        self.declare_partials('Const_KS', 'A_fl')
        self.declare_partials('Const_KS', 'A_stiff')

    def compute(self, inputs, outputs):


        C_wb = inputs['C_wb']
        h    = inputs['h']
        t_sk = inputs['t_sk']
        t_sp = inputs['t_sp']
        A_fl = inputs['A_fl']
        n_stiff = inputs['n_stiff']
        A_stiff = inputs['A_stiff']

        DesignVArs = np.array([C_wb, h, t_sk, t_sp, A_fl, n_stiff, A_stiff]).transpose()
        DesignVArs = np.array(DesignVArs.tolist())
        DesignVArs = np.reshape(DesignVArs, DesignVArs.size)

        if flag_lin == 'lin':
         outputs['Const_KS'] = beam_opt.ComputeResponseKSStress_opt_lin(DesignVArs)
        elif flag_lin == 'nonlin':
            outputs['Const_KS'] = beam_opt.ComputeResponseKSStress_opt_non_lin(DesignVArs)


    def compute_partials(self,inputs,partials):

        C_wb = inputs['C_wb']
        h = inputs['h']
        t_sk = inputs['t_sk']
        t_sp = inputs['t_sp']
        A_fl = inputs['A_fl']
        n_stiff = inputs['n_stiff']
        A_stiff = inputs['A_stiff']

        DesignVArs = np.array([C_wb, h, t_sk, t_sp, A_fl, n_stiff, A_stiff]).transpose()
        DesignVArs = np.array(DesignVArs.tolist())
        DesignVArs = np.reshape(DesignVArs, DesignVArs.size)

        if flag_lin == 'lin':
            KS_adj = beam_opt.ComputeAdjointKSstresses_opt_lin(DesignVArs)
        elif flag_lin == 'nonlin':
            KS_adj=beam_opt.ComputeAdjointKSstresses_opt_non_lin(DesignVArs)


        partials['Const_KS', 'C_wb'] = KS_adj[0::6]
        partials['Const_KS', 'h']    = KS_adj[1::6]
        partials['Const_KS', 't_sk'] = KS_adj[2::6]
        partials['Const_KS', 't_sp'] = KS_adj[3::6]
        partials['Const_KS', 'A_fl'] = KS_adj[4::6]
        partials['Const_KS', 'A_stiff'] = KS_adj[5::6]



class KSBuckling_constraint(om.ExplicitComponent):

            def setup(self):
                self.add_input('C_wb', shape=len(C_wb), units='mm')
                self.add_input('h', shape=len(h), units='mm')
                self.add_input('t_sk', shape=len(t_sk), units='mm')
                self.add_input('t_sp', shape=len(t_sp), units='mm')
                self.add_input('A_fl', shape=len(A_fl), units='mm**2')
                self.add_input('n_stiff', shape=len(n_stiff), units='mm')
                self.add_input('A_stiff', shape=len(A_stiff), units='mm**2')

                self.add_output('Const_KS_buckl')

                # Finite difference all partials.

                self.declare_partials('Const_KS_buckl', 'C_wb')
                self.declare_partials('Const_KS_buckl', 'h')
                self.declare_partials('Const_KS_buckl', 't_sk')
                self.declare_partials('Const_KS_buckl', 't_sp')
                self.declare_partials('Const_KS_buckl', 'A_fl')
                self.declare_partials('Const_KS_buckl', 'A_stiff')

            def compute(self, inputs, outputs):

                C_wb = inputs['C_wb']
                h = inputs['h']
                t_sk = inputs['t_sk']
                t_sp = inputs['t_sp']
                A_fl = inputs['A_fl']
                n_stiff = inputs['n_stiff']
                A_stiff = inputs['A_stiff']

                DesignVArs = np.array([C_wb, h, t_sk, t_sp, A_fl, n_stiff, A_stiff]).transpose()
                DesignVArs = np.array(DesignVArs.tolist())
                DesignVArs = np.reshape(DesignVArs, DesignVArs.size)

                if flag_lin == 'lin':
                    outputs['Const_KS_buckl'] = beam_opt.ComputeResponseKSBuckling_opt_lin(DesignVArs)
                elif flag_lin == 'nonlin':
                    outputs['Const_KS_buckl'] = beam_opt.ComputeResponseKSBuckling_opt_non_lin(DesignVArs)

            def compute_partials(self, inputs, partials):

                C_wb = inputs['C_wb']
                h = inputs['h']
                t_sk = inputs['t_sk']
                t_sp = inputs['t_sp']
                A_fl = inputs['A_fl']
                n_stiff = inputs['n_stiff']
                A_stiff = inputs['A_stiff']

                DesignVArs = np.array([C_wb, h, t_sk, t_sp, A_fl, n_stiff, A_stiff]).transpose()
                DesignVArs = np.array(DesignVArs.tolist())
                DesignVArs = np.reshape(DesignVArs, DesignVArs.size)

                if flag_lin == 'lin':
                    KS_adj = beam_opt.ComputeAdjointKSBuckling_opt_lin(DesignVArs)
                elif flag_lin == 'nonlin':
                    KS_adj = beam_opt.ComputeAdjointKSBuckling_opt_non_lin(DesignVArs)

                partials['Const_KS_buckl', 'C_wb'] = KS_adj[0::6]
                partials['Const_KS_buckl', 'h'] = KS_adj[1::6]
                partials['Const_KS_buckl', 't_sk'] = KS_adj[2::6]
                partials['Const_KS_buckl', 't_sp'] = KS_adj[3::6]
                partials['Const_KS_buckl', 'A_fl'] = KS_adj[4::6]
                partials['Const_KS_buckl', 'A_stiff'] = KS_adj[5::6]







if __name__ == "__main__":




    prob = om.Problem()

    flag_lin = input('Which optimization do you want to run ? (lin or nonlin) \n')

    Loads = (19,100,5000,50000)
    file_dir = os.path.dirname(os.path.realpath(__file__))

    path_opt = os.path.abspath(os.path.join(file_dir, ".."))
    final_directory_opt = os.path.join(path_opt, r'Optimization')
    if os.path.exists(final_directory_opt):
        shutil.rmtree(final_directory_opt)
    os.makedirs(final_directory_opt)




    if  flag_lin == 'lin':

        beam_opt = pyBeamOpt(file_dir, 'config_lin.cfg', Loads)
    elif flag_lin =='nonlin':

        beam_opt = pyBeamOpt(file_dir, 'config_non_lin.cfg', Loads)
    else:
        raise ValueError("the type of analysis to be performed is missing . Execution aborted")




    DVs   = (beam_opt.SetInitialParameters()[0])
    nDVs  = beam_opt.SetInitialParameters()[1]
    nProp = beam_opt.SetInitialParameters()[2]
    L     =beam_opt.SetInitialParameters()[3]
    print("Length", L)

    C_wb = DVs[::nDVs]
    h = DVs[1::nDVs]
    t_sk = DVs[2::nDVs]
    t_sp = DVs[3::nDVs]
    A_fl = DVs[4::nDVs]
    n_stiff = DVs[5::nDVs]
    A_stiff = DVs[6::nDVs]

    prob.model.add_subsystem('weight_comp', Weight(),
                             promotes_inputs=['C_wb', 'h', 't_sk', 't_sp', 'A_fl', 'n_stiff', 'A_stiff'])

    prob.model.add_subsystem('KS_comp_stress', KSStress_constraint(),
                             promotes_inputs=['C_wb', 'h', 't_sk', 't_sp', 'A_fl','n_stiff','A_stiff'])

    prob.model.add_subsystem('KS_comp_buckl', KSBuckling_constraint(),
                             promotes_inputs=['C_wb', 'h', 't_sk', 't_sp', 'A_fl', 'n_stiff', 'A_stiff'])
    if not n_stiff.all() == 0:
        prob.model.add_subsystem('Astiff_Afl_comp', om.ExecComp('g1 = A_stiff-A_fl', g1=np.ones(nProp),
                                                               A_stiff=np.ones(nProp), A_fl=np.ones(nProp)),promotes=['*'])

    # If properties > 1

    if nProp > 1:
        prob.model.add_subsystem('t_sk_con', om.ExecComp('g2 = t_sk[0:-1] - t_sk[1:]', g2=np.ones(nProp-1), t_sk=np.ones(nProp)),
                                                           promotes=['*'])
        prob.model.add_subsystem('t_sp_con', om.ExecComp('g3 = t_sp[0:-1] - t_sp[1:]', g3=np.ones(nProp-1), t_sp=np.ones(nProp)),
                                                          promotes=['*'])
        prob.model.add_subsystem('Afl_con',om.ExecComp('g4 = A_fl[0:-1] - A_fl[1:]', g4=np.ones(nProp- 1),A_fl=np.ones(nProp)),
                                 promotes=['*'])
        if not n_stiff.all() == 0:
            prob.model.add_subsystem('Astiff_con',om.ExecComp('g5 = A_stiff[0:-1] - A_stiff[1:]', g5=np.ones((nProp)-1),
                                                 A_stiff=np.ones(nProp)),promotes=['*'])


    # setup the optimization

    prob.driver = om.pyOptSparseDriver()
    #prob.driver = om.ScipyOptimizeDriver()
    prob.driver.options['optimizer'] = 'SLSQP'


    #prob.driver.options['maxiter'] = 10
    prob.driver.options['debug_print']=['desvars','objs','nl_cons','totals']

    prob.model.set_input_defaults('C_wb', C_wb, units='mm')
    prob.model.set_input_defaults('h', h, units='mm')
    prob.model.set_input_defaults('t_sk', t_sk, units='mm')
    prob.model.set_input_defaults('t_sp', t_sp, units='mm')
    prob.model.set_input_defaults('A_fl', A_fl, units='mm**2')
    prob.model.set_input_defaults('n_stiff', n_stiff, units='mm')
    prob.model.set_input_defaults('A_stiff', A_stiff, units='mm**2')




    filename = "opt_results.sql"
    recorder = om.SqliteRecorder(filename)
    prob.driver.add_recorder(recorder)
    prob.driver.recording_options['record_desvars'] = True
    prob.driver.recording_options['record_objectives'] = True
    prob.driver.recording_options['record_constraints'] = True
    prob.driver.recording_options['record_derivatives'] = True
    prob.driver.recording_options['includes'] = []
    prob.driver.recording_options['excludes'] = []


    """"""""""""""""""""""""""""""" ADD constraints , Objective function  and design variables """""""""""""""""""""""""""""""""

    #prob.model.add_design_var('C_wb')
    #prob.model.add_design_var('h')
    prob.model.add_design_var('t_sk', lower=0.5)
    prob.model.add_design_var('t_sp', lower=0.5)
    prob.model.add_design_var('A_fl')
    if not n_stiff.all() == 0:
        prob.model.add_design_var('A_stiff')

    prob.model.add_objective('weight_comp.Obj_f')

    prob.model.add_constraint('KS_comp_stress.Const_KS', upper=0.0)

    prob.model.add_constraint('KS_comp_buckl.Const_KS_buckl', upper=0.0)

    if not n_stiff.all() == 0:
        prob.model.add_constraint('g1', upper=0)
    if nProp > 1:
        prob.model.add_constraint('g2', equals=0)
        prob.model.add_constraint('g3', equals=0)
        prob.model.add_constraint('g4', equals=0)
        if not n_stiff.all() == 0:
            prob.model.add_constraint('g5', equals=0)

    prob.set_solver_print(level=0)

    prob.setup(check=False, mode='rev')



    prob.run_driver()

    print("Minimum Weight =", prob.get_val('weight_comp.Obj_f'))
    print("Constraint_stress =", prob.get_val('KS_comp_stress.Const_KS'))
    print("Constraint_buckling =", prob.get_val('KS_comp_buckl.Const_KS_buckl'))

    prob.cleanup()


    """  Record Data : """


    cr = om.CaseReader("opt_results.sql")

    case_names = cr.list_cases(out_stream=None)

    path = os.path.abspath(os.path.join(file_dir, "../Optimization"))
    final_directory = os.path.join(path, r'opt_iter')
    final_directory_figures = os.path.join(path, r'Figures')
    if not os.path.exists(final_directory):
        os.makedirs(final_directory)
    if not os.path.exists(final_directory_figures):
        os.makedirs(final_directory_figures)

    path_iter = os.path.abspath(os.path.join(file_dir, "../Optimization/opt_iter"))
    path_figures = os.path.abspath(os.path.join(file_dir, "../Optimization/Figures"))





    for it in range(0,len(case_names),1):

            case = cr.get_case(it)
            derivs=cr.get_case(it).derivatives

            final_directory_iter = os.path.join(path_iter, r'000'+ str(it))
            if not os.path.exists(final_directory_iter):
                os.makedirs(final_directory_iter)

            if derivs is not None :
             Weight_deriv = os.path.join(final_directory_iter, "Weight_sens" + ".txt")

             Deriv_obj = open(Weight_deriv, "w")
             Deriv_obj.write(
                             'Weight wrt t_sk = ' + str(derivs['weight_comp.Obj_f', 't_sk']) + '\n'+
                             'Weight wrt t_sp = ' + str(derivs['weight_comp.Obj_f', 't_sp']) + '\n'+
                             'Weight wrt A_fl = ' + str(derivs['weight_comp.Obj_f', 'A_fl']) + '\n')
             if not n_stiff.all() == 0:
                Deriv_obj.write('Weight wrt A_stiff = ' + str(derivs['weight_comp.Obj_f', 'A_stiff']))
             Deriv_obj.close()


             KS_deriv_buckl = os.path.join(final_directory_iter, "KSBuckling_sens" + ".txt")
             Deriv_KS_buckl = open(KS_deriv_buckl, "w")
             Deriv_KS_buckl.write(
                            'KS_buckl wrt t_sk = ' + str(derivs['KS_comp_buckl.Const_KS_buckl', 't_sk']) + '\n'
                            'KS_buckl wrt t_sp = ' + str(derivs['KS_comp_buckl.Const_KS_buckl', 't_sp'])+ '\n')
             if not n_stiff.all() == 0:
                 Deriv_KS_buckl.write('KS_buckl wrt A_stiff = ' + str(derivs['KS_comp_buckl.Const_KS_buckl', 'A_stiff']))
             Deriv_KS_buckl.close()



             KS_deriv = os.path.join(final_directory_iter, "KSStress_sens" + ".txt")
             Deriv_KS = open(KS_deriv, "w")
             Deriv_KS.write(
                 'KS wrt t_sk = ' + str(derivs['KS_comp_stress.Const_KS', 't_sk']) + '\n'
                                                                              'KS wrt t_sp = ' + str(
                     derivs['KS_comp_stress.Const_KS', 't_sp']) + '\n')
             if not n_stiff.all() == 0:
                 Deriv_KS.write('KS wrt A_stiff = ' + str(derivs['KS_comp_stress.Const_KS', 'A_stiff']))
             Deriv_KS.close()


            if derivs is None:
                pass

            Obj_func = os.path.join(final_directory_iter,  "Weight" + ".txt")

            outputs_obj = open(Obj_func , "w")
            outputs_obj.write(str(case['weight_comp.Obj_f']))
            outputs_obj.close()

            KS = os.path.join(final_directory_iter, "KSstress_constraint" + ".txt")

            outputs_KS = open(KS, "w")
            outputs_KS.write(str(case['KS_comp_stress.Const_KS']))
            outputs_KS.close()

            KS_buckl = os.path.join(final_directory_iter, "KSBuckling_constraint" + ".txt")

            outputs_KS_buckl = open(KS_buckl, "w")
            outputs_KS_buckl.write(str(case['KS_comp_buckl.Const_KS_buckl']))
            outputs_KS_buckl.close()



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

            if not n_stiff.all() == 0 :
                A_stiffeners = os.path.join(final_directory_iter, "A_stiff" + ".txt")

                outputs_A_stiff = open(A_stiffeners, "w")
                outputs_A_stiff.write(str(case['A_stiff']))
                outputs_A_stiff.close()



    """                          PLOT                        """


    """Weight"""
    fig1= plt.figure()
    for i in range(0, len(case_names), 1):
        case_plot = cr.get_case(i)
        plt.plot(i, float(case_plot['weight_comp.Obj_f']), 'ro', label='Weight')
    plt.legend(['Weight'])
    plt.xlabel('Driver Iterations', fontsize=18)
    plt.ylabel('Weight', fontsize=16)

    plt.savefig('../Optimization/Figures/ONERA_Weight_' + str(int(
            n_stiff[0])) + '_stiff_' + str(flag_lin)+'.pdf')  # save the figure to file
    plt.close(fig1)

    """KS constraint"""
    fig2=plt.figure()
    for i in range(0, len(case_names), 1):
        case_plot = cr.get_case(i)
        plt.plot(i, float(case_plot['KS_comp_stress.Const_KS']), 'ro')
    plt.legend(['KS_Stresses_constraint'])
    plt.xlabel('Driver Iterations', fontsize=18)
    plt.ylabel('KS_stresses Constraint', fontsize=16)
    plt.savefig('../Optimization/Figures/ONERA_KSStressconstraint_' + str(int(
            n_stiff[0]))  + '_stiff_' + str(flag_lin)+'.pdf')  # save the figure to file
    plt.close(fig2)



    """KS_Buckling constraint"""
    fig3 = plt.figure()
    for i in range(0, len(case_names), 1):
        case_plot = cr.get_case(i)
        plt.plot(i, float(case_plot['KS_comp_buckl.Const_KS_buckl']), 'ro')
    plt.legend(['KS_Buckling_constraint'])
    plt.xlabel('Driver Iterations', fontsize=18)
    plt.ylabel('KS_Buckling Constraint', fontsize=16)
    plt.savefig('../Optimization/Figures/ONERA_KSBucklingconstraint_' + str(int(
        n_stiff[0])) + '_stiff_' + str(flag_lin) + '.pdf')  # save the figure to file
    plt.close(fig3)


    """Design Variables """

    if nProp > 1:

        """Design Variables property > 1"""

        fig_root = plt.figure()
        """Root"""

        for i in range(0, len(case_names), 1):
            case_plot = cr.get_case(i)
            plt.plot(i, float(case_plot['t_sk'][0]), 'g^', label='t_sk')
            plt.plot(i, float(case_plot['t_sp'][0]), 'yo', label='t_sp')
            plt.plot(i, float(case_plot['A_fl'][0]), 'k+', label='A_fl')
            if not n_stiff.all() == 0:
             plt.plot(i, float(case_plot['A_stiff'][0]), 'co', label='A_stiff')



        if not n_stiff.all() == 0:
            plt.legend(['t_sk', 't_sp', 'A_fl','A_stiff'], bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.legend([ 't_sk', 't_sp', 'A_fl'], bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.xlabel('Driver Iterations', fontsize=18)
        plt.ylabel('Design Variables', fontsize=16)
        plt.tight_layout()
        plt.savefig('../Optimization/Figures/ONERA_DesignVAriables_at_root_' + str(int(
            n_stiff[0])) + '_stiff_' + str(flag_lin)+'.pdf')  # save the figure to file
        plt.close(fig_root)

        fig_tip = plt.figure()
        """Tip"""
        for i in range(0, len(case_names), 1):
            case_plot = cr.get_case(i)
            plt.plot(i, float(case_plot['t_sk'][-1]), 'g^', label='t_sk')
            plt.plot(i, float(case_plot['t_sp'][-1]), 'yo', label='t_sp')
            plt.plot(i, float(case_plot['A_fl'][-1]), 'k+', label='A_fl')
            if not n_stiff.all() == 0:
                plt.plot(i, float(case_plot['A_stiff'][-1]), 'co', label='A_stiff')

        if not n_stiff.all() == 0:
            plt.legend(['t_sk', 't_sp', 'A_fl', 'A_stiff'], bbox_to_anchor=(1.05, 1), loc='upper left')
        else:
            plt.legend(['t_sk', 't_sp', 'A_fl'], bbox_to_anchor=(1.05, 1), loc='upper left')

        plt.xlabel('Driver Iterations', fontsize=18)
        plt.ylabel('Design Variables', fontsize=16)
        plt.tight_layout()
        plt.savefig('../Optimization/Figures/ONERA_DesignVAriables_at_tip_' + str(int(
            n_stiff[0])) + '_stiff_' + str(flag_lin)+'.pdf')  # save the figure to file
        plt.close(fig_tip)


        """Design variables Vs Span otimization  """
        fig_DVs_span_opt = plt.figure()

        span = np.linspace(0,L,nProp)
        case_plot = cr.get_case(-1)

        plt.plot(span, (C_wb), 'ro', label='C_wb')
        plt.plot(span, (h), 'b^', label='h')
        plt.plot(span, (case_plot['t_sk'][:]), 'g^', label='t_sk')
        plt.plot(span, (case_plot['t_sp'][:]), 'yo', label='t_sp')
        plt.plot(span, (case_plot['A_fl'][:]), 'k+', label='A_fl')
        if not n_stiff.all() == 0:
            plt.plot(span, (case_plot['A_stiff'][:]), 'co', label='A_stiff')
            plt.legend(['C_wb', 'h','t_sk', 't_sp', 'A_fl', 'A_stiff'], bbox_to_anchor=(1.05, 1), loc='upper left')
        else:
         plt.legend([ 'C_wb', 'h','t_sk', 't_sp', 'A_fl'], bbox_to_anchor=(1.05, 1), loc='upper left')

        plt.xlabel('Span', fontsize=18)
        plt.ylabel('Design Variables', fontsize=16)
        plt.tight_layout()
        plt.savefig('../Optimization/Figures/ONERA_DesignVAriables_Vs_Span_Optimization' + str(int(
            n_stiff[0])) + '_stiff_' + str(flag_lin) + '.pdf')  # save the figure to file
        plt.close(fig_DVs_span_opt)


        """Design variables Vs Span Initial parameters """
        fig_DVs_span_init_design = plt.figure()
        plt.plot(span, (C_wb), 'ro', label='C_wb')
        plt.plot(span, (h), 'b^', label='h')
        plt.plot(span, (t_sk), 'g^', label='t_sk')
        plt.plot(span, (t_sp), 'yo', label='t_sp')
        plt.plot(span, (A_fl), 'k+', label='A_fl')
        if not n_stiff.all() == 0:
            plt.plot(span, (A_stiff), 'co', label='A_stiff')
            plt.legend(['C_wb', 'h','t_sk', 't_sp', 'A_fl', 'A_stiff'], bbox_to_anchor=(1.05, 1), loc='upper left')
        else:
            plt.legend(['C_wb', 'h','t_sk', 't_sp', 'A_fl'], bbox_to_anchor=(1.05, 1), loc='upper left')

        plt.xlabel('Span', fontsize=18)
        plt.ylabel('Design Variables', fontsize=16)
        plt.tight_layout()
        plt.savefig('../Optimization/Figures/ONERA_DesignVAriables_Vs_Span_initial_design' + str(int(
            n_stiff[0])) + '_stiff_' + str(flag_lin) + '.pdf')  # save the figure to file
        plt.close(fig_DVs_span_init_design)

    else:

        """ Design variables property =1 """
        fig4 = plt.figure()
        for i in range(0, len(case_names), 1):
            case_plot = cr.get_case(i)


            plt.plot(i, float(case_plot['t_sk']), 'g^', label='t_sk')
            plt.plot(i, float(case_plot['t_sp']), 'yo', label='t_sp')
            plt.plot(i, float(case_plot['A_fl']), 'k+', label='A_fl')
            if not n_stiff.all() == 0:
             plt.plot(i, float(case_plot['A_stiff']), 'co', label='A_stiff')


        if not n_stiff.all() == 0:
            plt.legend([ 't_sk', 't_sp', 'A_fl','A_stiff'], bbox_to_anchor=(1.05, 1), loc='upper left')
        else:
            plt.legend(['t_sk', 't_sp', 'A_fl'], bbox_to_anchor=(1.05, 1), loc='upper left')

        plt.xlabel('Driver Iterations', fontsize=18)
        plt.ylabel('Design Variables', fontsize=16)
        plt.tight_layout()
        plt.savefig('../Optimization/Figures/ONERA_DesignVAriables' + str(int(n_stiff)) + '_stiff_' + str(flag_lin)+'.pdf')  # save the figure to file
        plt.close(fig4)





