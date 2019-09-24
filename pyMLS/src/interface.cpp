/*
 * pyMLS, an open-source Moving Least Squares library
 *
 * Copyright (C) 2019 by the authors
 *
 * File developers: Rocco Bombardieri (Carlos III University Madrid)
 *                  Rauno Cavallaro (Carlos III University Madrid)
 *                  Ruben Sanchez (SciComp, TU Kaiserslautern)
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


#include "../include/interface.hpp"


void annInterface(Eigen::MatrixXd &dist, Eigen::MatrixXi &nnpos, int nDim, long int pointsUsed,
                  long int aeroNodes, long int strNodes, Eigen::MatrixXd strData, Eigen::MatrixXd aeroData){

    /*--- Initialize the ANN objects ---*/
    ANNkd_tree* annTree;
    ANNpointArray annData;

    /*--- Add some screen output to track what we are doing at each time ---*/
    std::cout << std::endl;
    std::cout << "Solving a " << nDim << "D problem." << std::endl;
    std::cout << "Total number of aerodynamic nodes: " << aeroNodes <<"." << std::endl;
    std::cout << "Total number of structural nodes (base mesh): " << strNodes <<"." << std::endl;
    std::cout << "Using " << pointsUsed <<" structural nodes for the search." << std::endl;

    double *aeroNode, *testDist;
    aeroNode = new double[nDim];
    testDist = new double[pointsUsed];

    int* testIndex;
    testIndex = new int[pointsUsed];

    int iNode, jNode;
    short iDim;

    /*--- Allocate the structural points to the ANN data structure ---*/
    annData = annAllocPts(strNodes, nDim);
    for (iNode=0; iNode < strNodes; iNode++) {
      for (iDim=0; iDim < nDim; iDim++) {
        annData[iNode][iDim] = strData(iNode,iDim);
      }
    }

    /*--- Create the tree with the structural nodes ---*/
    annTree = new ANNkd_tree(annData, strNodes, nDim);

    /*--- Loop over the fluid nodes and find positions and distances for the number of points used ---*/
    for (iNode = 0; iNode < aeroNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        aeroNode[iDim] = aeroData(iNode,iDim);
      }
      annTree->annkSearch(aeroNode, pointsUsed, testIndex, testDist, 1.0E-16);
      for (jNode=0; jNode < pointsUsed; jNode++) {
        nnpos(iNode,jNode) = testIndex[jNode];
        dist(iNode,jNode) = testDist[jNode];
      }
    }

    /*--- Delete the tree and remaining objects ---*/
    delete annTree;
    delete annData;
    delete [] aeroNode;
    delete [] testDist;
    delete [] testIndex;

}

long int ipow(long int base, long int exp)
{
    int result = 1;
    while (exp != 0)
    {
        if ((exp & 1) == 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }

    return result;
};

void mls_poly(Eigen::MatrixXd &P, Eigen::VectorXd &pp, Eigen::MatrixXd str_data, int type, Eigen::VectorXd aer_data){
    
    long int N = str_data.rows();
    long int dim = str_data.cols();
    long int i, j, k,l;
    
    // Eigen::MatrixXf P(N, 1+(dim)^type);
    // Eigen::VectorXf pp( 1+(dim)^type);
    // aer_data here is the single line of the aero_data matrix
    
    //re-scale nodes coordinates to avoid matrix ill-conditioning
    Eigen::VectorXd cg = str_data.colwise().sum(); 
    cg = cg/N;
    
    //std::cout << "cg = " << cg << "\n";
    Eigen::MatrixXd str = Eigen::MatrixXd::Zero(N,dim);
    for (j=0;j<N;j++)
    {
        for (i=0;i<str_data.cols();i++)
        {
        str(j,i) = str_data(j,i) - cg(i);
        }
    }
   
    Eigen::VectorXd aer = aer_data - cg;
    
    switch (type)
    {
        case 1:
            P.topLeftCorner(N, 1) = Eigen::MatrixXd::Ones(N,1); //setIdentity();
            P.block(0,1,N, dim) = str;
            pp(0) = 1.;
            for(i=1;i<dim+1;i++)
            {
            pp(i) = aer(i-1);
            }
            break;
        case 2:
            P.topLeftCorner(N, 1) = Eigen::MatrixXd::Ones(N,1); //setIdentity();
            P.block(0,1,N, dim) = str;
            pp(0) = 1.;
            for(i=1;i<dim+1;i++)
            {
            pp(i) = aer(i-1);
            }
            k = dim ;
            for (i=0;i<dim;i++)
            {
                for (j=i;j<dim;j++)
                {
                    k = k+1;
                    for (l=0;l<str.rows();l++)
                    {
                        P(l,k) = str(l,i) * str(l,j);
                    }
                    //P.col(k) = str.col(i)*str.col(j);
                    pp(k) = aer(i)* aer(j);
                     
                }
            }
            break;
            
    }
    
    //std::cout << "P = " << P <<"\n";
    //std::cout << "pp =" << pp <<"\n";
};

void   mls_wgt(Eigen::MatrixXd &W, Eigen::VectorXd dist, int type, double delta)
{
    double mx = delta * dist.maxCoeff();
    Eigen::VectorXd aux;
    Eigen::VectorXd aux2;
    Eigen::VectorXd one = Eigen::VectorXd::Ones(dist.size());
    
    switch (type)
    {
        case 1:
            aux = 1.*one - dist/mx;
            aux = aux.array().pow(2).matrix();
            W = aux.asDiagonal();
            break;
        case 2:
            aux = 1.*one - dist/mx;
            aux = aux.array().pow(4).matrix();
            aux = (aux.array()*(4.*dist/mx+1.*one).array()).matrix();
            W = aux.asDiagonal();
            break;
        case 3:
            aux = (1.*one - dist/mx).array().pow(6).matrix();
            aux2 = 35./3.*(dist/mx).array().pow(2).matrix() + 18./3.*dist/mx + 1.*one;
            //aux3 = 18./3.*dist/mx + 1.*one;
            aux = (aux.array()* aux2.array()).matrix();
            W = aux.asDiagonal();
            break;
        case 4:
            aux = (1.*one - dist/mx).array().pow(8).matrix();
            aux2 = 32.*(dist/mx).array().pow(3).matrix() + 25.*(dist/mx).array().pow(2).matrix() + 8.*dist/mx + 1.*one;
            //aux3 = 25.*(dist/mx).array().pow(2).matrix();
            //aux4 = 8.*dist/mx + 1.*one;
            aux = (aux.array()* aux2.array()).matrix();
            W = aux.asDiagonal();
            break;
    }
    
 //std::cout << "W = " << W <<"\n";

    
};

void check_MLS( Eigen::Map< Eigen::VectorXd > &norm_err, Eigen::MatrixXd str_data , Eigen::MatrixXd interpolation_matrix , Eigen::MatrixXd aero_data )
{
    /*% "check_MLS" performs a check on the consistency of the interface matrices
     % of MSL (in this case of H_VR_str_to_aero, but it can be considered
     % representative also for all the other matrices used to transfer FUNCTIONS
     % (and not DERIVATIVE of functions)
     
     % The check consists in a comparison between the coordinates of the rings
     % evaluated with the MLS approximation (RVRNODES_ext) and the real ones (RVRNODES).
     
     % norm_err is a vector with the error for each ring */
    //norm_err = Eigen::VectorXd(aero_data.rows())
    
    int i;
    Eigen::MatrixXd err = Eigen::MatrixXd::Zero(aero_data.rows(),3);
    norm_err = Eigen::VectorXd::Zero(aero_data.rows());
    Eigen::MatrixXd aero_ext =  interpolation_matrix * str_data;
    
    
    for (i=0;i<aero_data.rows();i++){
        err.row(i)    = aero_ext.row(i) - aero_data.row(i);
        
        norm_err(i) = err.row(i).norm();
    }
}

void pseudo_inv_eig(Eigen::MatrixXd &pinvA, Eigen::MatrixXd A, double toll)
{
    
    int m = A.rows();
    int n = A.cols();
    Eigen::MatrixXd U(m,m);
    Eigen::MatrixXd V(n,n);
    Eigen::MatrixXd S;

    //Eigen::VectorXd s_lap = Eigen::VectorXd::Zero(MD);
    //Eigen::MatrixXd S_lap(MD,MD);
   
    // Very careful!!! JacobiSVD is very efficient for small matrices and very slow for large
    //moreover here we ask for the full U and V matrices: if the code is slow think about switching to other Eigen algorithms
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
    
    /* It is chosen to relay on the Lapack's method DGESVD which is *supposedly* the same as used by matlab which I trust
     the method requires a 'linking' to the routine*/
    
    // Get the optimum work size.
    //dgesvd_ ( 'A' , 'A', m , n, A.data(), LDA , s_lap.data(), U.data(), m, VT.data(), n, WORKDUMMY, LWORK, INFO )
    
    /*LWORK = int(WORKDUMMY) + 32;
    Eigen::VectorXd WORK(LWORK);
    dgesvd_ ( 'A' , 'A', m , n, A.data(), LDA , s.data(), U.data(), m, VT.data(), n, WORK, LWORK, INFO )
    */
          
    U = svd.matrixU();
    V = svd.matrixV();
    S = svd.singularValues().asDiagonal();
    //check (to be commented)
    //std::cout << "Check: difference between A and U*S*V' = " << A - U*S*(V.transpose()) << "\n";    
    int r =0;   // Bug 04/2019

    Eigen::VectorXd s = S.diagonal();
    Eigen::VectorXd s_norm = Eigen::VectorXd::Zero(s.size());
    double tol = s.minCoeff();
    double maxs = s(0);
      
    for (int i=0;i<s.size();i++) 
    {           
            s_norm(i) = s(i)/maxs;       
            if (s_norm(i) > toll){
                r = r + 1;
            }
    }

     
    
    
    
    if (r == 0)  // also valid for the else case as for up
    {
        pinvA = Eigen::MatrixXd::Zero(m,n);
    }     
    else 
    {
             Eigen::MatrixXd S_bis = Eigen::MatrixXd::Zero(r,r);
    
             for (int i=0;i< r ;i++){
                  S_bis(i,i) = 1./s(i); 
              } 
             pinvA = V.block(0,0,V.rows(),r) * S_bis * U.block(0,0,U.rows(),r).transpose();
    }    
       /* DUBIGGING /////////////////////////////////////////////////////////
        std::cout << "S_norm = " << s_norm<< "\n"; 
        std::cout << "U = " << U << "\n";
        std::cout << "V = " << V << "\n";
        std::cout << "S = " << S << "\n";
        
        std::cout << "r = " << r << "\n";
        std::cout << "Press a key to continue...";
        std::cin.get();
    /////////////////////////////////////////////////////////    */
         
}

void mls_interface (std::vector<double> &interpolation_matrix_std, std::vector<double> &norm_err_std,
                    int str_nodenumb, int aero_nodenumb,
                    std::vector<double> str_data_std, std::vector<double> aero_data_std ,
                    int poly, int weight, long int points, double rmax,double delta, double toll) {
    
    /*
     % the first 3 outputs enable to evaluate the gradient of a function from
     % the values of the function
     
     %     -   h is the vector which have, for each quaery node, the coefficients of
     %           the shape functions such that f_(quaery_node)=H*f_(center)
     %            this is a vector with dimension (PARAM.MLS.poly1 * n_quaery_points)
     
     %    -      nnpos is a matrix   [n_quaery_points x PARAM.MLS.points] :
     %           for each quaery_node it gives the ID of the "points" closest
     %           centers: this is necessary to build the final interface matrix because
     %           it provides, for each row, which column has to be filled)
     
     */
    // INTIIALIZATION ----------------------------------
    
    int i, j, coeff;
    // conversion of input arrays into matrices  (same memory position of the preceding array)
    double *str_data_pointer = &str_data_std[0];
    double *aero_data_pointer = &aero_data_std[0];
    
    
    Eigen::Map< Eigen::MatrixXd > str_data(str_data_pointer,str_nodenumb,3)  ;
    Eigen::Map< Eigen::MatrixXd > aero_data(aero_data_pointer,aero_nodenumb,3)  ;

    //std::cout << " str_data =" << str_data << "\n";
           
    double *interpolation_matrix_pointer = &interpolation_matrix_std[0];
    double *norm_err_pointer = &norm_err_std[0];
    Eigen::Map< Eigen::MatrixXd > interpolation_matrix(interpolation_matrix_pointer,aero_nodenumb,str_nodenumb)  ;
    Eigen::Map< Eigen::VectorXd > norm_err(norm_err_pointer,aero_nodenumb)  ;

    int min_p;
    int DEFAULT_POLY = 2;
    int DEFAULT_WEIGHT = 2;
    double DEFAULT_RADIUS = pow(10.,+16);
    long int strN = str_data.rows();
    long int aerN = aero_data.rows();
    int dim = str_data.cols();
    long int dim2 = aero_data.cols();
    int nloc;
    long int aux_vector_size;
    
    //Eigen::MatrixXi nnpos(aerN= length(aero_data),points);
    Eigen::MatrixXi nnpos;
    //Eigen::VectorXd h = Eigen::VectorXd::Zero(aerN*points);
    Eigen::VectorXd h;
    
    Eigen::MatrixXd dist(aerN,points);
    
    // P = [ones(N,1), str]; //initialize auxiliary MatrixX
    Eigen::MatrixXd P(strN, 1 + ipow(dim,poly));
    Eigen::VectorXd pp(1 + ipow(dim,poly));
    Eigen::MatrixXd aux;
    Eigen::MatrixXd b;
    Eigen::MatrixXd A;
    //Eigen::MatrixXd pinvA;
    //----------------------------------------------------------  
    if (dim !=dim2)
    {
        printf("\nInconsistent data space dimensions.");      
    }
    
    if ((poly <= 0) || (poly > 2))
    {
        printf("\nWrong polynomial order. Correct values are 1 (linear) or 2 (parabolic). Value set to %d by default.",DEFAULT_POLY);
        poly = DEFAULT_POLY;
        std::cout << "OVERRIDEN. TO MAKE THE PROCESS MORE USER-AWARE, THE RUN IS TERMINATED. \n ERROR." << std::endl;
        throw -3;        
    }
    
    if ((weight <= 0) || (weight > 4))
    {
        printf("\nWrong Radial Basis weight order. Correct values are between 1 and 4. Value set to %d by default.",DEFAULT_WEIGHT);
        weight = DEFAULT_WEIGHT;
        std::cout << "OVERRIDEN. TO MAKE THE PROCESS MORE USER-AWARE, THE RUN IS TERMINATED. \n ERROR." << std::endl;
        throw -4;        
    }
    
    if (rmax == 0)
    {
        printf("\n No rmax value required. Value set to Inf.");
        rmax = DEFAULT_RADIUS;
        std::cout << "OVERRIDEN. TO MAKE THE PROCESS MORE USER-AWARE, THE RUN IS TERMINATED. \n ERROR." << std::endl;
        throw -5;        
    }
    

    switch (poly)
    {
        case 1:
            min_p = 1 + dim; break;
        case 2:
            min_p = 1 + dim*dim;  break;
    }
    
    if (min_p > strN) {  //(NEW CHECK!! 23/05/18)
        printf("\nRequired search points larger than available one in the structural mesh. ERROR. CODE IS STOPPED.");
        throw 113;
    }
    
    if (points < min_p)
    {
        //printf("\n Required search points too small. Value set to %d.", min_p);
        printf("\nDefault neighbor points smaller than minimum required (for the order of polynomial). Value set to %d.", min_p);    
        points = min_p;
    }
    
//    if (points > strN)    (OLD CHECK!! 23/05/18)
//    {
//        printf("\n Required too many search points. Value set to %ld.", strN);
//        points = strN;
//    }
    
    // further allocation needed
    Eigen::VectorXd h_internal(points*aerN);
    
    /* At this point the original matlab routine was writing some input files
     to be read by the C++ part of the code. Now the whole opeartion is done
     inside the C++ environment so we only need to pass the correct variables defined upside
     */
    
    // i need to give a dimension to matrix before entering the interface routine
    nnpos = Eigen::MatrixXi::Zero(aerN,points); // note that points won't change anymore from now on
    dist = Eigen::MatrixXd::Zero(aerN,points); // note that points won't change anymore from now on
    
    // Now the ANN interface routine is called

    annInterface( dist, nnpos, dim, points, aerN, strN, str_data, aero_data);
    
    for (i=0;i<dist.rows();i++)
    {
        for(j=0;j<dist.cols();j++)
        {
            dist(i,j) = sqrt(dist(i,j));
        }
    }
    
    //std::cout << "dist\n" << dist << "\n";
    //std::cout << "nnpos\n" << nnpos << "\n";
    
    coeff = 0;
    
    // for the methodology used to build the functions, see "Phenomenology of Nonlinear Aeroelastic Responses"
    for ( i=0;i<aerN;i++ )
    {
        nloc = 0;
        for ( j=0;j<dist.cols();j++ )
        {
            if (dist(i,j) < rmax)
            {
                nloc += 1;
            }
        }
        
        if (nloc < min_p)
        {

            nloc = min_p;
            printf("\nInserted radius too small (rmax parameter). Compact support enlarged to minimum size.\n");
        }
        
        /* what I have to do is (MATLAB) "str_data(nnpos(i,1:nloc),:)"
        I need to allocate something in loco*/
    
        Eigen::VectorXi aux_vector(nloc);
        Eigen::VectorXd aux_vector_dist(nloc);
        for ( j=0;j<nloc;j++ )
        {
            aux_vector(j) = nnpos(i,j);
            aux_vector_dist(j) = dist(i,j);
        }
        
        aux_vector_size = aux_vector_dist.size();
        Eigen::MatrixXd aux_matrix_str = Eigen::MatrixXd::Zero(aux_vector_size, str_data.cols());
        Eigen::MatrixXd W = Eigen::MatrixXd::Zero(aux_vector_size,aux_vector_size);
        igl::slice(str_data,aux_vector,1,aux_matrix_str); //igl::
        Eigen::VectorXd aero_data_row(aero_data.cols()) ;
        for ( j=0;j<aero_data.cols();j++ )
        {
          aero_data_row(j)  = aero_data(i,j);
        }
        
        mls_poly(P, pp, aux_matrix_str, poly, aero_data_row);

        mls_wgt( W, aux_vector_dist, weight, delta);
        //std::cout << "P = \n" << P << "\n";
        //std::cout << "W = \n" << W << "\n";
        b = P.block(0,0,aux_vector_size, 1 + ipow(dim,poly)).transpose() * W;
        A = P.block(0,0,aux_vector_size, 1 + ipow(dim,poly)).transpose() * W * P.block(0,0,aux_vector_size, 1 + ipow(dim,poly));
        //std::cout << "b = \n" << b << "\n";
        //std::cout << "A = \n" << A << "\n";

        //
        //pinvA = A.completeOrthogonalDecomposition().pseudoInverse();    // This is done trhough the QR method which is different from the SVD
        Eigen::MatrixXd pinvA;
        pseudo_inv_eig( pinvA, A, toll);
        
        
        Eigen::MatrixXd temp = pp.transpose() * (pinvA * b); //std::cout << "pp.transpose() * (pinvA * b) = \n" << pp.transpose() * (pinvA * b) << "\n";
            
        //std::cout << "coeff = \n" << coeff << "\n";
        //std::cout << "nloc = \n" << nloc << "\n";
        for ( j=coeff;j<coeff+nloc;j++ )
        {
        h_internal(j) = temp(0,j-coeff);
        }
        /* debugging //////////////////////////////////////
        if (i == 1-1){
        std::cout << "i = " << i <<"\n";
        std::cout << "pinvA = " << pinvA << "\n";
        std::ofstream file("test.txt");
        file <<  pinvA << '\n';
        std::cout << "Press a key to continue...";
        std::cin.get();
        }
        //////////////////////////////////////////////// */
        coeff += nloc;
        
        P.setZero();
        pp.setZero();
        temp.resize(0,0);
        pinvA.resize(0,0);
        aux_vector.resize(0); // this may be a problem
        aux_vector_dist.resize(0);
        aux_matrix_str.resize(0,0);
        aero_data_row.resize(0);
    }

    h = h_internal.segment(0,coeff);
    
    //std::cout << "h = \n" << h << "\n";
    //std::cout << "nnpos = \n" << nnpos << "\n";

    coeff = 0;

    //std::cout << interpolation_matrix << std::endl;

    for (i=0;i<nnpos.rows();i++){
        for (j=0;j<nnpos.cols();j++){
            interpolation_matrix(i,nnpos(i,j)) = h(coeff);
            coeff +=1;
        }
    }
    //std::cout << "interpolation_matrix = \n" << interpolation_matrix << "\n";

    check_MLS( norm_err,  str_data , interpolation_matrix , aero_data );
     
    
    for (int i=0;i<aero_data.cols()+1;i++){
    printf("*(norm_err + %d) : %f\n",  i, *(norm_err_pointer + i) );
    }
    

};
