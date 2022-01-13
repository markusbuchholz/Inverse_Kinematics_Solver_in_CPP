// (column, row)

//g++ 2R_2D_Kinematics.cpp -o t -I/usr/include/python3.8 -lpython3.8

#include "matplotlibcpp.h"

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include<vector>


using Eigen::MatrixXd;
using Eigen::MatrixXf;

namespace plt = matplotlibcpp;



MatrixXf FWD_Kinemeatic (float theta1, float theta2){

    float c1 = std::cos(theta1);
    float c2 = std::cos(theta2);
    float s1 = std::sin(theta1);
    float s2 = std::sin(theta2);

    float c12 = std::cos(theta1 + theta2);
    float s12 = std::sin(theta1 + theta2);
    
    float l1 = 1;
    float l2 = 1;


    MatrixXf FWD(2,1);

    FWD(0,0) = l1 * c1 + l2 * c12;
    FWD(1,0) = l1 * s1 + l2 * s12;




    return FWD;
}




MatrixXf computePseudoInverse(float theta1, float theta2)
{

    MatrixXf J(2,2);
    
    float c1 = std::cos(theta1);
    float c2 = std::cos(theta2);
    float c12 = std::cos(theta1 + theta2);
    
    float s1 = std::sin(theta1);
    float s2 = std::sin(theta2);
    float s12 = std::sin(theta1 + theta2);
    
    float l1 = 1;
    float l2 = 1;

    J(0,0) = -l1 * s1 - l2 * s12; J(0,1) = -l2 * s12;
    J(1,0) =  l1 * c1 + l2 * c12; J(1,1) =  l2 * c12;

    Eigen::MatrixXf pinvJ = J.completeOrthogonalDecomposition().pseudoInverse();


    return pinvJ;
}

MatrixXf invJacobian( float theta1, float theta2){

    MatrixXf invJ (2,2);

    float c1 = std::cos(theta1);
    float c2 = std::cos(theta2);
    float c12 = std::cos(theta1 + theta2);
    
    float s1 = std::sin(theta1);
    float s2 = std::sin(theta2);
    float s12 = std::sin(theta1 + theta2);
    
    float l1 = 1;
    float l2 = 1;

    float factor = 1/(l1 * l2 * s2);

    invJ(0,0) = l2 * c12; 
    invJ(0,1) = l2 * s12;
    invJ(1,0) = -l1 * c1 - l2 * c12;
    invJ(1,1) = -l1 * s1 - l2 * s12;



    return factor*invJ;


}

MatrixXf applyAlgorithm(){



   
    MatrixXf e(2,1);
    MatrixXf Xd (2,1);
    MatrixXf xiFWD (2,1);
    MatrixXf FWD (2,2);

    MatrixXf i_1_theta (2,1);
    MatrixXf i_theta (2,1);

    i_theta(0,0) =  M_PI/3;
    i_theta(1,0) = M_PI/4; 

 
   
    FWD = FWD_Kinemeatic(i_theta(0,0), i_theta(1,0));


    Xd(0,0) = 0.366; Xd(1,0) = 1.366; 
    e = Xd - FWD;





    while((std::abs(e(0,0)) > 0.000001) || (std::abs(e(1,0)) > 0.000001))
    {

        MatrixXf invJ = computePseudoInverse(i_theta(0,0), i_theta(1,0));

 
        
        i_1_theta = i_theta + invJ * e;

        FWD = FWD_Kinemeatic(i_1_theta(0,0), i_1_theta(1,0));
       

        e = Xd - FWD;
 

        i_theta = i_1_theta;
        
        std::cout << i_theta* (180/M_PI) << std::endl;

        

    }

    return i_theta;
}



void plotRobot(MatrixXf theta){

std::vector<float> xX;
std::vector<float> yY;

float theta1 = theta(0,0);
float theta2 = theta(1,0);

float l1 = 1;
float l2 = 1;

float x0 = 0;
float y0 = 0;
float x1 = l1 * std::cos(theta1);
float y1 = l1 * std::sin(theta1);


float x2 = l1 * std::cos(theta1) + l2 * std::cos(theta1 + theta2);
float y2 = l1 * std::sin(theta1) + l2 * std::sin(theta1 + theta2);

xX.push_back(x0);
xX.push_back(x1);
xX.push_back(x2);

yY.push_back(y0);
yY.push_back(y1);
yY.push_back(y2);

plt::figure_size(600, 400);
plt::title("2R XY Robot");

plt::plot(xX,yY);

plt::show();



}

int main()
{



    MatrixXf thetaRobot = applyAlgorithm();

  

  
    MatrixXf check = FWD_Kinemeatic(M_PI/6, M_PI/2);
    std::cout<<check<<"\n";
    plotRobot(thetaRobot);

}