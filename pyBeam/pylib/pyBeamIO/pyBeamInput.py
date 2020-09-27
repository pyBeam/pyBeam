#!/usr/bin/env python
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
import pdb
import os, sys, shutil, copy
import numpy as np
import scipy as sp
import scipy.linalg as linalg
from math import *

# ----------------------------------------------------------------------
#  Input file loading
# ----------------------------------------------------------------------
class Point:
  """ Description. """

  def __init__(self):
    self.ID = 0
    self.Coord0 = np.zeros((3,1)) # marks the initial position of the nodes: to this we superimpose the mode shape      
    self.Coord = np.zeros((3,1))
    self.Vel = np.zeros((3,1))
    self.Force = np.zeros((3,1))
    
  def GetCoord0(self):
    return self.Coord0
 

  def GetCoord(self):
    return self.Coord

  def GetVel(self):
    return self.Vel

  def GetForce(self):
    return self.Force

  def GetID(self):
    return self.ID

  def SetCoord0(self, val_Coord):
    x, y, z = val_Coord
    self.Coord0[0] = x
    self.Coord0[1] = y
    self.Coord0[2] = z


  def SetCoord(self, val_Coord):
    x, y, z = val_Coord
    self.Coord[0] = x
    self.Coord[1] = y
    self.Coord[2] = z

  def SetVel(self, val_Vel):
    vx, vy, vz = val_Vel
    self.Vel[0] = vx
    self.Vel[1] = vy
    self.Vel[2] = vz

  def SetForce(self, val_Force):
    fx, fy, fz = val_Force
    self.Force[0] = fx
    self.Force[1] = fy
    self.Force[2] = fz
   
  def SetID(self,ID):
    self.ID = int(ID)

class Element:  # for boundary elements
    
  def __init__(self):    
    self.Conn = np.zeros((2,1), dtype=int)
    self.ID = 0
    self.Property = 0  
    self.AuxVect = np.zeros((3,1))
       
  def SetID(self,ID):    
    self.ID = ID     
      
  def SetProperty(self,Property):    
    self.Property = Property         
    
  def SetConnectivity(self,val_Conn): # line element
    node1, node2, = val_Conn
    self.Conn[0] = node1
    self.Conn[1] = node2 
    
  def SetAuxVector(self,Auxval_Coord):
    x, y, z = Auxval_Coord
    self.AuxVect[0] = x
    self.AuxVect[1] = y
    self.AuxVect[2] = z
      
  def GetNodes(self):
    return self.Conn

  def GetProperty(self):
    return self.Property

  def GetID(self):
    return self.ID

  def GetAuxVector(self):
    return self.AuxVect

class RBE2_elem:  # for boundary elements
    
  def __init__(self):
    self.Conn = np.zeros((2,1), dtype=int)
    self.ID = 0

  def SetConnectivity(self,val_Conn): # line element
    node1, node2, = val_Conn
    self.Conn[0,0] = node1
    self.Conn[1,0] = node2
    
  def GetNodes(self):
    return self.Conn

  def SetID(self,ID):
    self.ID = ID




class Property:
  """ Description. """
    
  def __init__(self):
    self.A = 0
    self.Iyy = 0
    self.Izz = 0
    self.Jt = 0
    self.t_sk = 0;       # skin  thickness 
    self.t_sp = 0;           # spar thickness
    self.A_stiff = 0;        # stiffener Area
    self.A_fl = 0;           # flanges Area 
    self.h = 0;              # box tot height 
    self.C_wb = 0;           # box tot length 
    self.n_stiff = 0;        # number of stiffeners
    self.format = "N";         # N tipical one in which inertias are defined
                             # S on in which wing box sizes are defined

  def SetA(self,A):
    self.A = A 
    
  def SetIyy(self,Iyy):
    self.Iyy = Iyy
    
  def SetIzz(self,Izz):
    self.Izz = Izz
    
  def SetJt(self,Jt):
    self.Jt = Jt

  def Sett_sk(self,t_sk):
    self.t_sk = t_sk
    
  def Sett_sp(self,t_sp):
    self.t_sp = t_sp   
    
  def SetA_stiff(self,A_stiff):
    self.A_stiff = A_stiff    

  def SetA_fl(self,A_fl):
    self.A_fl = A_fl  
    
  def Seth(self,h):
    self.h = h  

  def SetC_wb(self,C_wb):
    self.C_wb = C_wb  
    
  def Setn_stiff(self,n_stiff):
    self.n_stiff = n_stiff 
    
  def SetFormat(self,Pformat):
    self.Format = Pformat  
    
  def GetFormat(self):
    return self.Format  
    
  def GetA(self):
    return self.A

  def GetIyy(self):
    return self.Iyy

  def GetIzz(self):
    return self.Izz

  def GetJt(self):
    return self.Jt

  def GetC_wb(self):
    return self.C_wb 

  def Geth(self):
    return self.h 

  def Gett_sk(self):
    return self.t_sk
    
  def Gett_sp(self):
    return self.t_sp

  def GetA_fl(self):
    return self.A_fl 

  def Getn_stiff(self):
    return self.n_stiff 
    
  def GetA_stiff(self):
    return self.A_stiff   

    

def readDimension(Mesh_file):

    nDim = 0

    with open(Mesh_file, 'r') as meshfile:
      print('--> Reading mesh file: ' + Mesh_file + '.')
      while 1:
        line = meshfile.readline()
        if not line:
          break
        if line.strip():
          if (line[0] == '%'):
            continue	  
        pos = line.find('NDIM')
        if pos != -1:
          line = line.strip('\r\n')
          line = line.split("=",1)
          nDim = (int(line[1]))
          continue

    return nDim      

def readMesh(Mesh_file,nDim):
	
    nPoint = 0
    node = []
     
    with open(Mesh_file, 'r') as meshfile:
      #print('Opened mesh file ' + Mesh_file + '.')
      while 1:
        line = meshfile.readline()
        if not line:
          break	
        if line.strip():
          if (line[0] == '%'):
            continue
        pos = line.find('NPOIN')
        if pos != -1:
          line = line.strip('\r\n')
          line = line.split("=",1)
          nPoint = int(line[1])
          for iPoint in range(nPoint):
            node.append(Point())
            line = meshfile.readline()
            line = line.strip('\r\n')
            line = line.split() ## important modification in case the formatting includes tabs
            x = float(line[0])
            y = float(line[1])
            z = 0.0
            if nDim == 3:
              z = float(line[2])
            node[iPoint].SetCoord((x,y,z))
            node[iPoint].SetCoord0((x,y,z))
            node[iPoint].SetID(int(iPoint+1)) # sequential ID
          continue	

    return node, nPoint	

def readConstr(Mesh_file):	  

    nConstr = 0         
     
    with open(Mesh_file, 'r') as meshfile:
      #print('Opened mesh file ' + Mesh_file + '.')
      while 1:
        line = meshfile.readline()
        if not line:
          break	
        if line.strip():
          if (line[0] == '%'):
            continue
        pos = line.find('NCONSTR')
        if pos != -1:
          line = line.strip('\r\n')
          line = line.split("=",1)
          nConstr = int(line[1])
          Constr = np.zeros((nConstr,2), dtype=int)
          for iConstr in range(nConstr):
            line = meshfile.readline()
            line = line.strip('\r\n')
            line = line.split() ## important modification in case the formatting includes tabs
            nid = int(line[0])
            dofid = int(line[1])
            Constr[iConstr,0] = nid;Constr[iConstr,1] = dofid; 
          continue	

    return Constr, nConstr	

def readConnectivity(Mesh_file):	

    Elem = []
          
    with open(Mesh_file, 'r') as meshfile:
      #print('Opened mesh file ' + Mesh_file + '.')
      while 1:
        line = meshfile.readline()
        if not line:
          break	
        if line.strip():
          if (line[0] == '%'):
            continue	
        pos = line.find('NELEM')
        if pos != -1:
          line = line.strip('\r\n')
          line = line.replace(" ", "")
          line = line.split("=",1)
          nElem = int(line[1])
          for iElem in range(nElem):
              Elem.append(Element())
              line = meshfile.readline()
              line = line.strip('\r\n')
              line = line.split() ## important modification in case the formatting includes tabs
              nodes = line[0:3]#.split()  ## important modification in case the formatting includes tabs
              AuxVector = line[3:6]
              Elem[iElem].SetConnectivity([ int(nodes[0]), int(nodes[1]) ])   
              Elem[iElem].SetID(iElem)   
              Elem[iElem].SetProperty(int(nodes[2]))   
              Elem[iElem].SetAuxVector([ float(AuxVector[0]) , float(AuxVector[1]), float(AuxVector[2]) ])
          continue
        else:
          continue	
        
    return Elem, nElem	

def Check_RBE2(RBE2, nRBE2):
    
    if nRBE2 !=0:
       masters = np.zeros((nRBE2,1))
       slaves = np.zeros((nRBE2,1))
    
       for i in range(0,nRBE2):
          masters[i,0] = RBE2[i].GetNodes()[0,0];
          slaves[i,0] = RBE2[i].GetNodes()[1,0];
       
       check = np.isin(slaves,masters)
    
       if np.any(check) == True:
          raise ValueError('RBE2 issue: a slave cannot be master as well. Execution aborted') 


def readRBE2(Mesh_file):	

    RBE2 = []
          
    with open(Mesh_file, 'r') as meshfile:
      #print('Opened mesh file ' + Mesh_file + '.')
      while 1:
        line = meshfile.readline()
        if not line:
          break	
        if line.strip():
          if (line[0] == '%'):
            continue	
        pos = line.find('NRBE2')
        if pos != -1:
          line = line.strip('\r\n')
          line = line.replace(" ", "")
          line = line.split("=",1)
          nRBE2 = int(line[1])
          for iRBE2 in range(nRBE2):
              RBE2.append(RBE2_elem())
              line = meshfile.readline()
              line = line.strip('\r\n')
              line = line.split() ## important modification in case the formatting includes tabs
              nodes = line[0:3]#.split()  ## important modification in case the formatting includes tabs
              AuxVector = line[3:6]
              RBE2[iRBE2].SetConnectivity([ int(nodes[0]), int(nodes[1]) ])   
              RBE2[iRBE2].SetID(iRBE2)   

          continue
        else:
          continue	    
    Check_RBE2(RBE2, nRBE2)    
    
    return RBE2, nRBE2	

def readProp(Prop_file):	

    Prop = []    
              
    with open(Prop_file, 'r') as propfile:
      print('--> Reading property file: ' + Prop_file + '.')
      
      #for line in propfile:
      #  # For each line, check if line contains the string
      #  if string_to_search in line:
      #          return True          
          
      while 1:
        line = propfile.readline()
        if not line:
          break	
        if line.strip():
          if (line[0] == '%'):
            continue	
                
        pos = line.find('NPROP')
        #print ("NPROPS", pos)
        pos2 = line.find('NPROPS')
        #print ("NPROPS2", pos2)
        
        if pos != -1 and pos2 ==-1:
          print('--> OLD FORMAT FOR SPECIFYING PROPERTIES. WILL BE DISCONTINUED')
          line = line.strip('\r\n')
          line = line.replace(" ", "")
          line = line.split("=",1)
          nProp = int(line[1])
          for iProp in range(nProp):   
              print("Storing Property", iProp)
              Prop.append(Property())
              line = propfile.readline()
              line = line.strip('\r\n')
              line = line.split() ## important modification in case the formatting includes tabs    
              A = float(line[0]); Iyy = float(line[1]); Izz = float(line[2]); Jt = float(line[3]); 
              Prop[iProp].SetA(A); Prop[iProp].SetIyy(Iyy); Prop[iProp].SetIzz(Izz); Prop[iProp].SetJt(Jt); 
              print("Success!",iProp)
          
#        pos = line.find('NPROPS')
        elif pos != -1 and pos2 !=-1:
          print('--> NEW FORMAT FOR SPECIFYING PROPERTIES.')
          #if iOLD == 1:
          #  raise ValueError('Cannot specify properties in new and old format. Execution aborted') 
          line = line.strip('\r\n')
          line = line.replace(" ", "")
          line = line.split("=",1)
          nProp = int(line[1])
          for iProp in range(nProp):    
              Prop.append(Property())
              line = propfile.readline()
              line = line.strip('\r\n')
              line = line.split() ## important modification in case the formatting includes tabs    
              Pformat = line[0]
              Prop[iProp].SetFormat(Pformat)
              line = propfile.readline()
              line = line.strip('\r\n')
              line = line.split() ## important modification in case the formatting includes tabs   
                
              if (Pformat == 'N'):        
                A = float(line[0]); Iyy = float(line[1]); Izz = float(line[2]); Jt = float(line[3]); 
                Prop[iProp].SetA(A); Prop[iProp].SetIyy(Iyy); Prop[iProp].SetIzz(Izz); Prop[iProp].SetJt(Jt)              
              elif (Pformat == 'S'):
                C_wb = float(line[0]);  h = float(line[1]);  t_sk = float(line[2]);  t_sp = float(line[3]); 
                A_fl = float(line[4]);  n_stiff =  int(line[5]);   A_stiff = float(line[6]);
                Prop[iProp].SetC_wb(C_wb); 
                Prop[iProp].Seth(h); 
                Prop[iProp].Sett_sk(t_sk); 
                Prop[iProp].Sett_sp(t_sp);
                Prop[iProp].SetA_fl(A_fl);                   
                Prop[iProp].Setn_stiff(n_stiff);   
                Prop[iProp].SetA_stiff(A_stiff);  
              else: 
                raise ValueError("Unknown paramter for Property CARD input. Execution aborted") 
                
    return Prop, nProp


