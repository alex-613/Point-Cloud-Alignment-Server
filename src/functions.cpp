/**
* @file functions.cpp
* @author Alex Ranne/Hisham Iqbal
* @data 18th of September, 2019
* @brief This file consists of all the function definitions. Most definitions uses the Eigen library to manipulate
* martices and performs matrix algebra. And the Winsock library for socket programming.
*
*
*
*/

#include "functions.h"
#pragma comment(lib,"WS2_32")

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////Least Squares Optimization//////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

//Matrix generator function

Eigen::MatrixXd functions::matrixGen(int rows, int cols) {
	Eigen::MatrixXd points = Eigen::MatrixXd::Random(rows, cols);

	return points;
}

//Quaternion converter

Eigen::MatrixXd functions::quaternionConv2R(double ANGLE, int x, int y, int z) {

	double q[4] = { cos(ANGLE), sin(ANGLE) * x, sin(ANGLE) * y, sin(ANGLE) * z };

	double xx = q[1] * q[1];
	double xy = q[1] * q[2];
	double xz = q[1] * q[3];
	double xw = q[1] * q[0];

	double yy = q[2] * q[2];
	double yz = q[2] * q[3];
	double yw = q[2] * q[0];

	double zz = q[3] * q[3];
	double zw = q[3] * q[0];

	double m00 = 1 - 2 * (yy + zz);
	double m01 = 2 * (xy - zw);
	double m02 = 2 * (xz + yw);

	double m10 = 2 * (xy + zw);
	double m11 = 1 - 2 * (xx + zz);
	double m12 = 2 * (yz - xw);

	double m20 = 2 * (xz - yw);
	double m21 = 2 * (yz + xw);
	double m22 = 1 - 2 * (xx + yy);

	Eigen::MatrixXd R(3, 3);

	R(0, 0) = m00;
	R(0, 1) = m01;
	R(0, 2) = m02;

	R(1, 0) = m10;
	R(1, 1) = m11;
	R(1, 2) = m12;

	R(2, 0) = m20;
	R(2, 1) = m21;
	R(2, 2) = m22;

	return R;
}

//Predefined transformation
Eigen::MatrixXd functions::predefinedTransformation(Eigen::MatrixXd points, Eigen::MatrixXd R, int translation[3]) {
	Eigen::MatrixXd pointsNew = points * R;

	for (int i = 0; i < pointsNew.rows(); i++) {
		pointsNew(i, 0) += translation[0];
		pointsNew(i, 1) -= translation[1];
		pointsNew(i, 2) += translation[2];
	}
	return pointsNew;

}
//Finding centroids of two data clouds
Eigen::MatrixXd functions::centroidFinding(Eigen::MatrixXd points, Eigen::MatrixXd pointsNew) {
	double sum_1x = 0;
	double sum_1y = 0;
	double sum_1z = 0;

	double sum_2x = 0;
	double sum_2y = 0;
	double sum_2z = 0;

	for (int i = 0; i < points.rows(); i++) {
		sum_1x += points(i, 0);
		sum_1y += points(i, 1);
		sum_1z += points(i, 2);

		sum_2x += pointsNew(i, 0);
		sum_2y += pointsNew(i, 1);
		sum_2z += pointsNew(i, 2);
	}
	double avg_cloud_1x = sum_1x / points.rows();
	double avg_cloud_1y = sum_1y / points.rows();
	double avg_cloud_1z = sum_1z / points.rows();

	double avg_cloud_2x = sum_2x / points.rows();
	double avg_cloud_2y = sum_2y / points.rows();
	double avg_cloud_2z = sum_2z / points.rows();

	Eigen::MatrixXd centroids(2, 3);

	centroids(0, 0) = avg_cloud_1x;
	centroids(0, 1) = avg_cloud_1y;
	centroids(0, 2) = avg_cloud_1z;

	centroids(1, 0) = avg_cloud_2x;
	centroids(1, 1) = avg_cloud_2y;
	centroids(1, 2) = avg_cloud_2z;


	return centroids;

}

Eigen::MatrixXd functions::findH(Eigen::MatrixXd points, Eigen::MatrixXd pointsNew, Eigen::MatrixXd centroids) {
	Eigen::MatrixXd diff1 = Eigen::MatrixXd::Constant(points.rows(), 3, 0);
	Eigen::MatrixXd diff2 = Eigen::MatrixXd::Constant(points.rows(), 3, 0);

	Eigen::MatrixXd H = Eigen::MatrixXd::Constant(3, 3, 0);
	Eigen::MatrixXd a = Eigen::MatrixXd::Constant(3, 1, 0);
	Eigen::MatrixXd b = Eigen::MatrixXd::Constant(1, 3, 0);

	for (int i = 0; i < points.rows(); i++) {
		diff1(i, 0) = points(i, 0) - centroids(0, 0); //x
		diff1(i, 1) = points(i, 1) - centroids(0, 1); // y
		diff1(i, 2) = points(i, 2) - centroids(0, 2); // z

		diff2(i, 0) = pointsNew(i, 0) - centroids(1, 0);// x
		diff2(i, 1) = pointsNew(i, 1) - centroids(1, 1);// y
		diff2(i, 2) = pointsNew(i, 2) - centroids(1, 2);// z

		a << diff1(i, 0),
			diff1(i, 1),
			diff1(i, 2);

		b << diff2(i, 0), diff2(i, 1), diff2(i, 2);
		H = H + a * b;

	}
	return H;

}

//Jacobian SVD decomposition of matrix to recover R

Eigen::MatrixXd functions::singleValueDecomposition(Eigen::MatrixXd H) {
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(H, Eigen::ComputeThinU | Eigen::ComputeThinV);

	Eigen::MatrixXd V = svd.matrixV();
	Eigen::MatrixXd U = svd.matrixU();
	Eigen::MatrixXd S = svd.singularValues();

	Eigen::MatrixXd X = U * V.transpose();

	Eigen::MatrixXd correction(3, 3);

	correction << 1, 0, 0,
		0, 1, 0,
		0, 0, X.determinant();


	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {

			if (abs(X(i, j)) < 1e-5) {
				X(i, j) = 0.00;
			}

		}
	}

	if (X.determinant() == 1) {
		std::cout << "Congratz" << std::endl;


	}
	else if (X.determinant() == -1) {
		std::cout << "Humm, check with Admin" << std::endl;
		X = U * correction * V.transpose();
	}
	return X;
}

Eigen::MatrixXd functions::findTranslation(Eigen::MatrixXd centroids, Eigen::MatrixXd X) {
	Eigen::MatrixXd t = Eigen::MatrixXd::Constant(3, 1, 0);
	Eigen::MatrixXd centroid_2_Vector = Eigen::MatrixXd::Constant(3, 1, 0);
	Eigen::MatrixXd centroid_1_Vector = Eigen::MatrixXd::Constant(3, 1, 0);

	centroid_1_Vector << round(centroids(0, 0)),
		round(centroids(0, 1)),
		round(centroids(0, 2));

	centroid_2_Vector << round(centroids(1, 0)),
		round(centroids(1, 1)),
		round(centroids(1, 2));

	t = centroid_2_Vector - X * centroid_1_Vector;
	return t;
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////UDP socket programming//////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

void functions::initialization() {
	WSADATA data;
	WORD version = MAKEWORD(2, 2);

	int wsOK = WSAStartup(version, &data);

	if (wsOK != 0) {
		//Error check
		std::cout << "Can't start Winsock!" << wsOK << std::endl;
		return;
	}

}

void functions::binding(SOCKET &in) {
	sockaddr_in serverHint;
	//Fill the serverHint structure
	serverHint.sin_addr.S_un.S_addr = ADDR_ANY;
	serverHint.sin_family = AF_INET;
	serverHint.sin_port = htons(54000);

	if (bind(in, (sockaddr*)& serverHint, sizeof(serverHint)) == SOCKET_ERROR) {
		std::cout << "Can't bind socket!" << WSAGetLastError() << std::endl;
		return;
	}
}

std::vector<Eigen::MatrixXd> functions::receivingloop(SOCKET in, int pcln) {
	char buf[1024];
	std::vector<Eigen::MatrixXd> matrixArray;

	sockaddr_in client;
	int clientLength = sizeof(client);

	for (int i = 0; i < pcln; i++) {
		//Clear anything in these variables
		ZeroMemory(&client, clientLength);
		ZeroMemory(buf, 1024);

		//Wait for message
		std::cout << "Waiting to receive a message (index no: " << i+1 << ")" << std::endl;
		int bytesIn = recvfrom(in, buf, 1024, 0, (sockaddr*)& client, &clientLength);
		if (bytesIn == SOCKET_ERROR) {
			std::cout << "Error receiving from client " << WSAGetLastError() << std::endl;
			continue;
		}

		// Display message and client info
		char clientIp[256]; // Create enough space to convert the address byte array
		ZeroMemory(clientIp, 256); // to string of characters

		// Convert from byte array to chars
		inet_ntop(AF_INET, &client.sin_addr, clientIp, 256);
		//std::cout << "Message recv from " << clientIp << " : " << buf << std::endl;

		// The system receives the buf
		std::cout << "Message recv from " << clientIp << std::endl;
		std::cout << buf << std::endl;

		//=============================================================//

		//Step 1:convert buf to string
		std::string bufStr(buf);

		//Step 2:remove all null characters, convert to double array
		bufStr.erase(std::find(bufStr.begin(), bufStr.end(), '\0'), bufStr.end());
		std::vector<double> bufDoublevect;

		std::stringstream ss(bufStr);

		double j;

		while (ss >> j) {
			bufDoublevect.push_back(j);

			if (ss.peek() == ',')
				ss.ignore();
		}

		//Step3: Convert the bufDoublevect into an Eigen matrix
		int pointsets = bufDoublevect.size() / 3;
		Eigen::MatrixXd bufMatrix(pointsets, 3);
		int counter = 0;

		for (int k = 0; k < bufMatrix.rows(); k++) {
			for (int h = 0; h < bufMatrix.cols(); h++) {
				bufMatrix(k, h) = bufDoublevect[counter];
				counter++;
			}
		}

		//Step4: Append the matrix into a matrix array

		matrixArray.push_back(bufMatrix);

	}

	return matrixArray;

}

void functions::termination(SOCKET &in) {

	closesocket(in);
	WSACleanup();
}

