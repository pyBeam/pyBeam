from pyBeam import CInput

def parseInput(BEAM_config, inputs,Constr, nConstr):
    
    inputs.SetBeamLength(BEAM_config['B_LENGTH'])
    inputs.SetWebThickness(BEAM_config['W_THICKNESS'])
    inputs.SetWebHeight(BEAM_config['W_HEIGHT'])
    inputs.SetFlangeWidth(BEAM_config['F_WIDTH'])
    inputs.SetYoungModulus(BEAM_config['Y_MODULUS'])
    inputs.SetPoisson(BEAM_config['POISSON'])
    inputs.SetDensity(BEAM_config['RHO'])
    inputs.SetLoad(BEAM_config['LOAD'])
    inputs.SetFollowerFlag(BEAM_config['FOLLOWER_FLAG'])
    inputs.SetLoadSteps(BEAM_config['LOAD_STEPS'])
    inputs.SetNStructIter(BEAM_config['N_STRUCT_ITER'])
    inputs.SetConvCriterium(BEAM_config['CONV_CRITERIUM'])
    
    #Now setting the constraints
    inputs.SetnConstr(nConstr)
    for iConstr in range(nConstr): 
        inputs.SetSingleConstr( iConstr, int(Constr[iConstr,0]), int(Constr[iConstr,1]) )
