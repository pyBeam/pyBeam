/*
 * pyBeam, an open-source Beam Solver
 *
 * Copyright (C) 2019 by the authors
 *
 * File developers: Rocco Bombardieri (Carlos III University Madrid)
 *                  Rauno Cavallaro (Carlos III University Madrid)
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



#include "../include/property.h"
#include<cstdlib>
#include<iostream>

using namespace std;
   


void CProperty::SetSectionProperties2()
{
    
    int n_stiff=6;
    addouble C_wb=2;
    addouble t_sk=2*pow(10,-3);
    addouble t_sp=5*pow(10,-3);
    addouble h=0.7;
    addouble A_stiff=5*pow(10,-4);
    addouble A_fl=12*pow(10,-4);
    addouble b=(C_wb)/((n_stiff/2)-1);  //distance within stiffeners
    addouble a=0;  
    addouble summ_ys=0;
    
    
    
  
    if (n_stiff=0 ){
        
        summ_ys=0;                    //stiffeners moment of inertia respect vertical axis = 0
                                        
    } else
    {
    int r= ((n_stiff)/2)%2;
    if (r==0)
    {
       for (int j=0;j<=(n_stiff/4)-1;j+=1)
       {
           
           a=a + pow(0.5+j,2);
  
       }
          
    } else 
    {
      for (int i=1;i<=(((n_stiff/2)-1)/2);i+=1)
      {
        a=a+pow(i,2);  
        cout<<"a="<<a<<endl;
      }
    }
     addouble summ_ys=4*pow(b,2)*a;
    }   
  
   cout<<"summ"<<summ_ys;

   addouble Izz_fl=4*A_fl*pow(C_wb/2,2);
   
    addouble A_skin=C_wb*t_sk;
    addouble A_sp=t_sp*h;
 
    addouble  A=A_skin+A_sp+n_stiff*A_stiff+4*A_fl;
  
    addouble  Iyy=2*(C_wb*(pow(t_sk,3)/12)+(C_wb*t_sk)*(pow(h/2,2)))+2*(t_sp*(pow(h,3)/12))+n_stiff*((A_stiff+A_fl)*pow(h/2,2));
    addouble  Izz=2*(t_sk*(pow(C_wb,3)/12))+2*(h*(pow(t_sp,3)/12)+A_sp*(pow(C_wb/2,2)))+(A_stiff*summ_ys+Izz_fl);
    addouble  Jt=(2*t_sp*t_sk*pow(C_wb,2)*pow(h,2))/(C_wb*t_sp+h*t_sk);
    
    
    
    
    cout<<"Iyy"<<Iyy<<endl;
    cout<<"Izz"<<Izz<<endl;
    cout<<"Jt"<<Jt<<endl;
    
   
      
}
