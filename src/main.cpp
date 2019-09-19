/**
* @file main.cpp
* @author Alex Ranne/Hisham Iqbal
* @date 18th of September, 2019
* @brief This is the main file of the project
*
* \mainpage Description
* In the field of Augmented Reality and mixed reality, where lots of decives are involved in creating a believable interavtive environment for the user, it is important that engineers and programmers create a program that
* enables interaction between those devices. In this program, a connectionless UDP (user datagram protocol) is use to transmit small amounts of critical data, in order to optimize the hologram projections by the hololens.
* The optimization is done through the use of the method of least squares. This is a common method used for coordinate mapping and optimizing positions of various coordinate systems within a 3D space. Given two or more point clouds, expressed in their coordinate system
* and the geometric relation between them is unknown, then a mamized rotation matrix and translation vector must be found using quaternion mathematics, SVD and matrix algebra, eigen value method, or iteration(trial and error).
* This program is written to weave the sending/receiving and calculation process together. Note that this is the server (receiving) side of the connection.
*
* \section library_sec Libraries used in this project
* \subsection library1 The Winsock WS2tcpip.h library
* This library enables the program to utilize the powerful windows sockets functions, such as socket binding, conversion from big to little endian, sending and receiving char arrays from specific ip addresses and more.
*
* \subsection library2 Eigen
* The Eigen library is a very powerful, header only library that has lots of built in matrix manipulation tools. Such as transposing, inversing, conducting SVD, enabling functions to return the Matrix type, etc etc.
*
* \section code_sec Code in main.cpp breakdown
* @code
* functions::initialization();
* @endcode
* Call the init function which starts winsock, allowing us to create, bind sockets. See functions.cpp
*
* @code
*SOCKET in = socket(AF_INET, SOCK_DGRAM, 0);
*sockaddr_in serverHint;
*functions::binding(in);
*@endcode
* Create a UDP style socket called in, and initialize the serverHint sock address structure to be filled in. The binding function binds the ip address with the socket/port number, and creates a tuple for later use.
* See functions.cpp
*
* @code
*sockaddr_in client;
*int clientLength = sizeof(client);
*int pcln = 2; //Number of point clouds
*
*std::vector<Eigen::MatrixXd> matrixArray = functions::receivingloop(in, pcln);
*@endcode
* Initialize another sockaddr stucture to be filled in later.
* The use should modify the number of point clouds accordingly, depending on how many point clouds we are receiving.
* The function returns a vector of matrices, stored for further analysis by the rest of the program.
* @code
* 	MatrixXd centroids = functions::centroidFinding(points, pointsNew);
*
*	std::cout << "\n" << centroids << std::endl;
* @endcode
*
*  This step uses one of the generated functions to find the centroid of two point clouds, which is sourced from the
*  microsoft hololens.
*
* @code
*  	MatrixXd H = functions::findH(points, pointsNew, centroids);
*
*	std::cout << "\n" << H << std::endl;
*
*	MatrixXd X = functions::singleValueDecomposition(H);
*
*	std::cout << "\n" << X << std::endl;
* @endcode
*
* The H matrix is a very common matrix that has appeared lots in the documentations regarding the mathematics behind the
* least squares optimization method. In short, the H matrix is formed by firstly finding the difference matrix
* (the difference of each coordinate to the cnetroid in each point cloud), then taking one row of the difference matrix for one point , transpose it,
* and multiple by the same row in a the difference matrix of the second point cloud. Then summing of these matrices together.
*
* Following that step, we then apply SVD (single value decomposiiton, see outside reference) to this matrix. And recover the rotational matrix.
*
* @code
* 	MatrixXd t = functions::findTranslation(centroids, X);
*
*	std::cout << "\n" << t << std::endl;
* @endcode
*
* Finally, we recover the translation matrix using the rotation matrix.
*
*/
#include <iostream>
#include <WS2tcpip.h>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <cmath>
#include <string>
#include <sstream>
#include <algorithm>

#include <vector>
#include <stdio.h>

#include "functions.h"

int main() 
{
	//////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////START OF UDP TRANSMISSION////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////

	/**
	Initialize winsock
	*/

	functions::initialization();


	/**
	Socket creation
	*/

	SOCKET in = socket(AF_INET, SOCK_DGRAM, 0);
	sockaddr_in serverHint;

	functions::binding(in);

	/**
	Main setup, the server receives
	*/
	sockaddr_in client;
	
	int pcln = 2; //Number of point clouds

	std::vector<Eigen::MatrixXd> matrixArray = functions::receivingloop(in, pcln);

	/**
	Close winsock and clean up
	*/
	functions::termination(in);

	//////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////END OF UDP TRANSMISSION//////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////

	/**
	Extract the data from the matrix, not including the commas
	*/
	std::cout << "\n" << matrixArray[0] << std::endl;
	std::cout << "\n" << matrixArray[1] << std::endl;

	////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////Least squares optimization//////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////

	// Find the centroids

	Eigen::MatrixXd centroids = functions::centroidFinding(matrixArray[0], matrixArray[1]);

	std::cout << "\n" << centroids << std::endl;

	// Find H matrix

	Eigen::MatrixXd H = functions::findH(matrixArray[0], matrixArray[1], centroids);

	std::cout << "\n" << H << std::endl;

	Eigen::MatrixXd X = functions::singleValueDecomposition(H);

	std::cout << "\n" << X << std::endl;

	//Find the t vector

	Eigen::MatrixXd t = functions::findTranslation(centroids, X);

	std::cout << "\n" << t << std::endl;
	std::cout << "Any key to exit" << std::endl;
	std::cin.get();
	return 0;
}