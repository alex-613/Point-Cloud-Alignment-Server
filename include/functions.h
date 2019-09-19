/**
* @file functions.h
* @author Alex Ranne/Hisham Iqbal
* @date 18th of September, 2019
* @brief Functions header file of the project, declears all functions used.
*
*
*
*
*/

#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <vector>
#include <WS2tcpip.h>

#pragma comment(lib,"WS2_32")

using namespace Eigen;

class functions 
{
	public:
		//Least squares optimisation
		static MatrixXd matrixGen(int rows, int cols);
		static MatrixXd quaternionConv2R(double ANGLE, int x, int y, int z);
		static MatrixXd predefinedTransformation(Eigen::MatrixXd points, Eigen::MatrixXd R, int translation[3]);
		static MatrixXd centroidFinding(Eigen::MatrixXd points, Eigen::MatrixXd pointsNew);
		static MatrixXd findH(Eigen::MatrixXd points, Eigen::MatrixXd pointsNew, Eigen::MatrixXd centroids);
		static MatrixXd singleValueDecomposition(Eigen::MatrixXd H);
		static MatrixXd findTranslation(Eigen::MatrixXd centroids, Eigen::MatrixXd X);

		////WINSOCK
		static void initialization();
		static void binding(SOCKET &in);
		static std::vector<Eigen::MatrixXd> receivingloop(SOCKET in, int pcln);
		static void termination(SOCKET &in);
};