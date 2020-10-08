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

/**  This function overloads the other function with the same name but here as inputs are passed
 wingbox thicknesses instead of inertias. Hence, inertias are evaluated.*/

void CProperty::SetSectionProperties(passivedouble C_wb_, passivedouble h_, passivedouble t_sk_,  
                              passivedouble t_sp_, passivedouble A_fl_, 
                              int n_stiff_, passivedouble A_stiff_)
{
    C_wb = C_wb_;
    h = h_;
    t_sk = t_sk_;
    t_sp = t_sp_;
    A_fl = A_fl_;
    n_stiff = n_stiff_;
    A_stiff = A_stiff_;
    
    addouble b=(C_wb)/(((n_stiff+4)/2)-1);      //distance within stiffeners
    
    int r= ((n_stiff)/2)%2;
    
    addouble a=0;  
    addouble summ_ys=0;
    
     if (n_stiff == 0 ){  //std::cout << "--> no stiffeners "<< std::endl;
        
        summ_ys=0;}     //stiffeners moment of inertia respect vertical axis = 0                                
   
     else if (r ==0)  {  //std::cout << "--> Even Number of stiffeners"<< std::endl;
        for (int j=0;j<=(n_stiff/4)-1;j+=1){
           a=a + pow(0.5+j,2);}   }
     else  { //std::cout << "--> Odd Number of stiffeners"<< std::endl;

        for (int i=1;i<=(((n_stiff/2)-1)/2);i+=1) {
            a=a+pow(i,2);}   }
    
    summ_ys=4*pow(b,2)*a;   
   
    addouble A_skin = (C_wb+t_sp)*t_sk;                              //  skin area
    addouble A_sp   = (t_sp*(h-t_sk));                               //  spar area
       
    addouble Iyy_skin= (C_wb+t_sp)*(pow(t_sk,3)/12)+A_skin*(pow(((h/2)),2));  // Iyy of the skin
    addouble Iyy_spar= t_sp*(pow((h-t_sk),3)/12);                            //Iyy of the spar
    addouble Izz_skin= t_sk*(pow((C_wb-t_sp),3)/12);                         // Izz of the skin 
    addouble Izz_spar= (h-t_sk)*(pow(t_sp,3)/12)+A_sp*(pow(((C_wb/2)),2));   // Izz of the spar
    addouble Izz_fl = 4*A_fl*pow((C_wb/2),2); // Izz of the flanges 
    
    A_b= n_stiff*A_stiff + 4*A_fl;  /// Number of stiffeners
       
    Iyy_b=(n_stiff)*((A_stiff)*pow((h/2),2))  + 4*A_fl*pow((h/2),2);  // Iyy of the booms system for ideal shell theory
    
    Izz_b=(A_stiff*summ_ys+Izz_fl);                                  // Izz of the booms system for ideal shell theory
    
    A = 2*A_skin + 2*A_sp  +  n_stiff*A_stiff  +  4*A_fl;       //  total area
           
    Iyy=2*Iyy_skin  + 2*Iyy_spar + (n_stiff)*((A_stiff)*pow((h/2),2))  + 4*A_fl*pow((h/2),2);   
         
    Izz=2*Izz_skin  +   2*Izz_spar + (A_stiff*summ_ys+Izz_fl);
        
    Jt=(2*t_sp*t_sk*pow(C_wb,2)*pow(h,2))/(C_wb*t_sp+h*t_sk);
   
    J0=Iyy+Izz;
   
   //addouble Sy= A*(h/2);   // Y static moment of the section
   
   //addouble Sz=A*(C_wb/2);  // Z static moment of the section   
   
   isWBDV = 1;
   
//   std::cout << "Inertias of the boom never used. Static Moments neither" << std::endl;
    
}


void CProperty::RegisterInput_WB(void) {
    if (isWBDV == 1){
        AD::RegisterInput(C_wb);
        AD::RegisterInput(h);
        AD::RegisterInput(t_sk);
        AD::RegisterInput(t_sp);
        AD::RegisterInput(A_fl);    
        AD::RegisterInput(A_stiff);      
    }
    else if  (isWBDV == 0){
        AD::RegisterInput(A);
        AD::RegisterInput(Iyy);
        AD::RegisterInput(Izz);
        AD::RegisterInput(Jt); 
    }
}


void CProperty::GetGradient_WB(void) {
    if (isWBDV == 1){
        AD::RegisterInput(C_wb);
        AD::RegisterInput(h);
        AD::RegisterInput(t_sk);
        AD::RegisterInput(t_sp);
        AD::RegisterInput(A_fl);    
        AD::RegisterInput(A_stiff);      
    }
    else if  (isWBDV == 0){
        AD::RegisterInput(A);
        AD::RegisterInput(Iyy);
        AD::RegisterInput(Izz);
        AD::RegisterInput(Jt); 
    }
}




/*void CProperty::SetSectionProperties2()
{
    // input parameters 
    n_stiff =0;
    C_wb    =3000;
    t_sk    =3;
    t_sp    =5;
     h       =500;
    addouble A_stiff=50 ;
    addouble A_fl    =200;
    
    addouble b=(C_wb)/(((n_stiff+4)/2)-1);      //distance within stiffeners
    //cout<<"b"<<b<<endl;
    
    
   
    int r= ((n_stiff)/2)%2;
    addouble a=0;  
    addouble summ_ys=0;
    
     if (n_stiff == 0 ){
         
        //std::cout << "--> no stiffeners "<< std::endl;
        
        summ_ys=0;                           //stiffeners moment of inertia respect vertical axis = 0
                                        
     } else if (r==0) // Even number 
     {
        //std::cout << "--> Even Number of stiffeners"<< std::endl;
       for (int j=0;j<=(n_stiff/4)-1;j+=1)
        {
           
           a=a + pow(0.5+j,2);
  
        }
          
     } else  //odd nnumber 
    {
     //std::cout << "--> Odd Number of stiffeners"<< std::endl;
      for (int i=1;i<=(((n_stiff/2)-1)/2);i+=1)
      {
        a=a+pow(i,2);  
       
      }
    }
    summ_ys=4*pow(b,2)*a;   
    
   
    addouble A_skin = (C_wb+t_sp)*t_sk;                                        //  skin area
    addouble A_sp = (t_sp*(h-t_sk));                                          //  spar area
 
    addouble Iyy_skin=(C_wb+t_sp)*(pow(t_sk,3)/12)+A_skin*(pow(((h/2)),2));                // Iyy of the skin
    addouble Iyy_spar= t_sp*(pow((h-t_sk),3)/12);                                        //Iyy of the spar
    addouble Izz_skin= t_sk*(pow((C_wb-t_sp),3)/12);                                               // Izz of the skin 
    addouble Izz_spar= (h-t_sk)*(pow(t_sp,3)/12)+A_sp*(pow(((C_wb/2)),2));      // Izz of the spar
    addouble Izz_fl = 4*A_fl*pow((C_wb/2),2);                                  // Izz of the flanges 
    
    
     addouble Iyy_b=(n_stiff)*((A_stiff)*pow((h/2),2))  + 4*A_fl*pow((h/2),2);  // Iyy of the booms system for ideal shell theory
     addouble Izz_b=(A_stiff*summ_ys+Izz_fl);                                  // Izz of the booms system for ideal shell theory
    
    A=2*A_skin + 2*A_sp  +  n_stiff*A_stiff  +  4*A_fl;       //  total area
             
                  
             
    Iyy=2*Iyy_skin  + 2*Iyy_spar + (n_stiff)*((A_stiff)*pow((h/2),2))  + 4*A_fl*pow((h/2),2);   
   
      
    Izz=2*Izz_skin  +   2*Izz_spar + (A_stiff*summ_ys+Izz_fl);
    
    
    Jt=(2*t_sp*t_sk*pow(C_wb,2)*pow(h,2))/(C_wb*t_sp+h*t_sk);
   

    J0=Iyy+Izz;
   
   addouble Sy= A*(h/2);   // Y static moment of the section
   
   addouble Sz=A*(C_wb/2);  // Z static moment of the section
    
   
    
    cout<<"Atot="<< setprecision(20) <<A<<endl;
    cout<<"Iyy="<< setprecision(20)<<Iyy<<endl;
    cout<<"Izz="<< setprecision(20)<<Izz<<endl;
    cout<<"Jt="<< setprecision(20) <<Jt<<endl;
    //cout<<"Sy="<<Sy<<endl;
    //cout<<"Sz"<<Sz<<endl;
    cout<< "I_yy_b"<< setprecision(20) << Iyy_b<<endl;
    cout<< "I_zz_b"<< setprecision(20) << Izz_b<<endl;
    
   
      
}
 * */
