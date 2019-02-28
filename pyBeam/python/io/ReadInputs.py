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
    
  #def SetConnectivity_tria(self,val_Conn): # triangular element
  #  node1, node2, node3 = val_Conn
  #  self.Conn[0] = node1
  #  self.Conn[1] = node2          
  #  self.Conn[2] = node3    
    
  def GetNodes(self):
    return self.Conn

  def GetProperty(self):
    return self.Property

  def GetID(self):
    return self.ID

  def GetAuxVector(self):
    return self.AuxVect

class Property:
  """ Description. """
    
  def __init__(self):
    self.A = 0
    self.Iyy = 0
    self.Izz = 0
    self.Jt = 0

  def SetA(self,A):
    self.A = A 
    
  def SetIyy(self,Iyy):
    self.Iyy = Iyy    
    
  def SetIzz(self,Izz):
    self.Izz = Izz    
    
  def SetJt(self,Jt):
    self.Jt = Jt  
    
  def GetA(self):
    return self.A

  def GetIyy(self):
    return self.Iyy

  def GetIzz(self):
    return self.Izz

  def GetJt(self):
    return self.Jt


def readDimension(Mesh_file):

    nDim = 0

    with open(Mesh_file, 'r') as meshfile:
      print('Opened Structural mesh file ' + Mesh_file + '.')
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
              Elem[iElem].SetConnectivity([ int(nodes[0]), int(nodes[1]) ])
              Elem[iElem].SetID(iElem)   
              Elem[iElem].SetProperty(int(nodes[2]))   
              Elem[iElem].SetAuxVector([ float(AuxVector[0]) , float(AuxVector[1]), float(AuxVector[2]) ])
          continue
        else:
          continue	
        
    return Elem, nElem	

def readProp(Prop_file):	

    Prop = []    
    
    with open(Prop_file, 'r') as propfile:
      print('Opened property file ' + Prop_file + '.')
      while 1:
        line = propfile.readline()
        if not line:
          break	
        if line.strip():
          if (line[0] == '%'):
            continue	
        pos = line.find('NPROP')
        if pos != -1:
          line = line.strip('\r\n')
          line = line.replace(" ", "")
          line = line.split("=",1)
          nProp = int(line[1])
          for iProp in range(nProp):    
              Prop.append(Property())
              line = propfile.readline()
              line = line.strip('\r\n')
              line = line.split() ## important modification in case the formatting includes tabs    
              A = float(line[0]); Iyy = float(line[1]); Izz = float(line[2]); Jt = float(line[3]); 
              Prop[iProp].SetA(A); Prop[iProp].SetIyy(Iyy); Prop[iProp].SetIzz(Izz); Prop[iProp].SetJt(Jt); 
              
    return Prop, nProp

