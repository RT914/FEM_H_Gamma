#include <Eigen/Dense>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <vector>

#include "NewtonRaphsonMethod.h"
#include "Square.h"
#include "FEM.h"
#include "Export.h"
#include "CalKernel.h"

// Interpolation Function Sophia
// [\alpha][\beta][\gamma][i-xi][j-xj][k-xk]
double SophiaX[3][3][3][3][3][3] = {};
// [\alpha][\beta][\gamma][i-yi][j-yj][k-yk]
double SophiaY[3][3][3][3][3][3] = {};
// [\alpha][\beta][\gamma][i-zi][j-zj][k-zk]
double SophiaZ[3][3][3][3][3][3] = {};

// Interpolation Function Chloe
// [i-xi][i-yi][i-zi]
double Chloe[3][3][3] = {};

// Interpolation Function Aria
// [i-xi][i-yi][i-zi][j-xj][j-yj][j-zj]
double Aria[3][3][3][3][3][3] = {};

// Interpolation Function Mia
// [\alpha][\beta][\gamma][i-xi][j-xj][k-xk][l-xl]
double MiaX[3][3][3][3][3][3][3] = {};
// [\alpha][\beta][\gamma][i-yi][j-yj][k-yk][l-yl]
double MiaY[3][3][3][3][3][3][3] = {};
// [\alpha][\beta][\gamma][i-zi][j-zj][k-zk][l-zl]
double MiaZ[3][3][3][3][3][3][3] = {};

double kappa = 1.0;

// �s��P�̑傫����[3�~NumberOfParticles][NumberOfParticles]�@�Ƃ���D
Eigen::MatrixXd MatrixP(3 * NumberOfParticles, NumberOfParticles);

// �s��Q�̑傫����[NumberOfParticles][NumberOfParticles]�@�Ƃ���D
Eigen::MatrixXd MatrixQ(NumberOfParticles, NumberOfParticles);

// �s��R�̑傫����[NumberOfParticles][NumberOfParticles]�@�Ƃ���D
Eigen::MatrixXd MatrixR(NumberOfParticles, NumberOfParticles);

// �x�N�g��b�̑傫���� [NumberOfParticles]�@�Ƃ���D
Eigen::VectorXd Vectorb(NumberOfParticles);

// �x�N�g��c�̑傫���� [NumberOfParticles]�@�Ƃ���D
Eigen::VectorXd Vectorc(NumberOfParticles);

// �s��T�̑傫����[4�~NumberOfParticles][2�~NumberOfParticles]�@�Ƃ���D
Eigen::MatrixXd MatrixS(4 * NumberOfParticles, 2 * NumberOfParticles);

// �s��T�̑傫����[4�~NumberOfParticles][4�~NumberOfParticles]�@�Ƃ���D
Eigen::MatrixXd MatrixT(4 * NumberOfParticles, 4 * NumberOfParticles);

// �x�N�g��d�̑傫���� [2�~NumberOfParticles]�@�Ƃ���D
Eigen::VectorXd Vectord(2 * NumberOfParticles);

// �x�N�g��e�̑傫���� [4�~NumberOfParticles]�@�Ƃ���D
Eigen::VectorXd Vectore(4 * NumberOfParticles);

Eigen::VectorXd VectorDeltaPhi(3 * NumberOfParticles);

Eigen::VectorXd VectorDeltaTheta(NumberOfParticles);

// For H
// �s��P�̑傫����[NumberOfParticles][NumberOfParticles]�@�Ƃ���D
Eigen::MatrixXd MatrixP_H(3 * NumberOfParticles, NumberOfParticles);

// �x�N�g��b�̑傫���� [NumberOfParticles]�@�Ƃ���D
Eigen::VectorXd Vectorb_H(NumberOfParticles);

// �s��R�̑傫����[NumberOfParticles][NumberOfParticles]�@�Ƃ���D
Eigen::MatrixXd MatrixR_H(3 * NumberOfParticles, 3 * NumberOfParticles);

// �x�N�g���̑傫���� [3�~NumberOfParticles]�@�Ƃ���D
Eigen::VectorXd Vectorc_H(3 * NumberOfParticles);

Eigen::VectorXd Newton(Square square) {
	double NormVectorDeltaPhi = 1.0;
	double NormVectorDeltaTheta = 1.0;
	double NormVectorDelta = 1.0;
	vector <double>norm;
	Eigen::MatrixXd MatrixS(4 * NumberOfParticles, 2 * NumberOfParticles); // �s��S���`
	Eigen::VectorXd Vectord(2 * NumberOfParticles); // �x�N�g��d���`
	int SquarePointsNumber = square.points.size();

	// �x�N�g���̑傫���� [3�~NumberOfParticles]�@�Ƃ���D
	Eigen::VectorXd barphi(3 * NumberOfParticles);
	// �x�N�g���̑傫���� [3�~NumberOfParticles]�@�Ƃ���D
	Eigen::VectorXd bartheta(NumberOfParticles);
	// �x�N�g���̑傫����[3�~NumberOfParticles]�@�Ƃ���D�Q�Ɨp
	Eigen::VectorXd re_barphi(3 * NumberOfParticles);
	// �x�N�g���̑傫���� [3�~NumberOfParticles]�@�Ƃ���D�����T���@�ŗp����D
	Eigen::VectorXd barphi_prime(3 * NumberOfParticles);
	// �x�N�g���̑傫���� [NumberOfParticles]�@�Ƃ���D�����T���@�ŗp����D
	Eigen::VectorXd bartheta_prime(NumberOfParticles);


	// ���W�̎擾
	for (int i = 0; i < SquarePointsNumber; i++) {
		barphi(3 * i) = square.points[i].position[0];
		barphi(3 * i + 1) = square.points[i].position[1];
		barphi(3 * i + 2) = square.points[i].position[2];
		bartheta(i) = square.points[i].theta;
	}

	// �Q�Ɨp�̍��W�̎擾
	for (int i = 0; i < SquarePointsNumber; i++) {
		re_barphi(3 * i) = square.points[i].reference_position[0];
		re_barphi(3 * i + 1) = square.points[i].reference_position[1];
		re_barphi(3 * i + 2) = square.points[i].reference_position[2];
	}

	// std::cout << barphi << "\n";

	int looptimes = 0;

	while (NormVectorDelta > 1.0e-6) {
		MatrixP = calMatrixP(square, barphi, bartheta); // �s��P���v�Z
		exportMatrixP(MatrixP);
		MatrixQ = calMatrixQ(square, barphi, bartheta); // �s��Q���v�Z
		if (looptimes == 0) {
			exportMatrixQ(MatrixQ);
		}
		MatrixR = calMatrixR(square, barphi, bartheta); // �s��R���v�Z
		if (looptimes == 0){ exportMatrixR(MatrixR); }
		Vectorb = calVectorb(square, barphi, bartheta); // �x�N�g��b���v�Z
		exportVectorb(Vectorb);
		Vectorc = calVectorc(square, barphi, bartheta); // �x�N�g��c���v�Z
		exportVectorc(Vectorc);

		
		for (int i = 0; i < NumberOfParticles;i++) {
			for (int j = 0; j < NumberOfParticles; j ++ ) {
				// ��s��S�̍쐬 ----------------------------------------
				// �s��S�ɍs��P���Z�b�g
				MatrixS(3 * i, j) = MatrixP(3 * i, j);
				MatrixS(3 * i + 1, j) = MatrixP(3 * i + 1, j);
				MatrixS(3 * i + 2, j) = MatrixP(3 * i + 2, j);

				// �s��S�ɍs��Q���Z�b�g
				MatrixS(i + 3 * NumberOfParticles, j) = MatrixQ(i, j);

				// �s��S�ɍs��R���Z�b�g
				MatrixS(i + 3 * NumberOfParticles, j + NumberOfParticles) = MatrixR(i, j);
			}

			// �s��S�̏o��
			exportMatrixS(MatrixS);
			// std::cout << MatrixS << std::endl;

			// ��x�N�g��d�̍쐬 ------------------------------------
			// �x�N�g��d�Ƀx�N�g��b���Z�b�g
			Vectord(i) = Vectorb(i);

			// �x�N�g��d�Ƀx�N�g��c���Z�b�g
			Vectord(i + NumberOfParticles) = Vectorc(i);

			// �x�N�g��d�̏o��
			exportVectord(Vectord);
		}

		MatrixT = calMatrixT(MatrixS); // �s��T���v�Z
		exportMatrixT(MatrixT);
		Vectore = calVectore(Vectord, MatrixS); // �x�N�g��e���v�Z
		exportVectore(Vectore);

		Eigen::FullPivLU<Eigen::MatrixXd> LU(MatrixT);
		Eigen::VectorXd VectorDelta = LU.solve(Vectore);

		// phi��theta�̕���
		for (int i = 0; i < NumberOfParticles; i++) {
			// Phi
			VectorDeltaPhi(3 * i) = VectorDelta(3 * i);
			VectorDeltaPhi(3 * i + 1) = VectorDelta(3 * i + 1);
			VectorDeltaPhi(3 * i + 2) = VectorDelta(3 * i + 2);

			// Theta
			VectorDeltaTheta(i) = VectorDelta(i + 3 * NumberOfParticles);
		}

		exportVectorDeltaPhi(VectorDeltaPhi);
		exportVectorDeltaTheta(VectorDeltaTheta);

		NormVectorDelta = VectorDelta.norm();
		norm.emplace_back(NormVectorDelta);
		NormVectorDeltaPhi = VectorDeltaPhi.norm();
		NormVectorDeltaTheta = VectorDeltaTheta.norm();

		std::cout << " NormVectorDelta : " << NormVectorDelta << std::endl;

		//Line Search with Wolfe
		double sigma = 0.5;
		double eps1 = 0.01;
		double eps2 = 0.9;
		double lambda = 1.0;
		int line_search_times = 1;

		barphi_prime = VectorDeltaPhi * lambda + barphi;
		bartheta_prime = VectorDeltaTheta * lambda + bartheta;

		Eigen::VectorXd Vectorb_prime = -1 * calVectorb(square, barphi_prime, bartheta_prime);
		Eigen::VectorXd Vectorc_prime = -1 * calVectorc(square, barphi_prime, bartheta_prime);
		Eigen::VectorXd Vectord_prime(2 * NumberOfParticles);
		Eigen::MatrixXd MatrixP_prime = calMatrixP(square, barphi_prime, bartheta_prime);
		Eigen::MatrixXd MatrixQ_prime = calMatrixQ(square, barphi_prime, bartheta_prime);
		Eigen::MatrixXd MatrixR_prime = calMatrixR(square, barphi_prime, bartheta_prime);
		Eigen::MatrixXd MatrixS_prime(4 * NumberOfParticles, 2 * NumberOfParticles);


		for (int i = 0; i < NumberOfParticles; i++)
		{
			for (int j = 0; j < NumberOfParticles; j++) {
				// ��s��S_prime�̍쐬 ----------------------------------------
				// �s��S_prime�ɍs��P_prime���Z�b�g
				MatrixS_prime(3 * i, j) = MatrixP_prime(3 * i, j);
				MatrixS_prime(3 * i + 1, j) = MatrixP_prime(3 * i + 1, j);
				MatrixS_prime(3 * i + 2, j) = MatrixP_prime(3 * i + 2, j);

				// �s��S_prime�ɍs��Q_prime���Z�b�g
				MatrixS_prime(i + 3 * NumberOfParticles, j) = MatrixQ_prime(i, j);

				// �s��S_prime�ɍs��R_prime���Z�b�g
				MatrixS_prime(i + 3 * NumberOfParticles, j + NumberOfParticles) = MatrixR_prime(i, j);
			}

			// ��x�N�g��d_prime�̍쐬 ------------------------------------
			// �x�N�g��d_prime�Ƀx�N�g��b_prime���Z�b�g
			Vectord_prime(i) = Vectorb_prime(i);
			// �x�N�g��d_prime�Ƀx�N�g��c_prime���Z�b�g
			Vectord_prime(i + NumberOfParticles) = Vectorc_prime(i);
		}


		//Armijo
		Eigen::VectorXd f_prime = Vectord_prime;
		Eigen::VectorXd f = Vectord;
		Eigen::VectorXd nabla_f1 = eps1 * lambda * MatrixS.transpose() * VectorDelta;
		Eigen::VectorXd right_f = f + nabla_f1;

		//curvature
		Eigen::VectorXd nabla_f_prime = MatrixS_prime.transpose() * VectorDelta;
		Eigen::VectorXd nabla_f2 = eps2 * MatrixS.transpose() * VectorDelta;

		while (f_prime.norm() > right_f.norm() || nabla_f_prime.norm() < nabla_f2.norm()) {
			// �ɁCbarphi_prime�Cbartheta_prime�̍X�V
			// lambda = 2.0 / (2 + line_search_times);
			lambda = sigma * lambda;
			// std::cout << line_search_times << std::endl;
			barphi_prime = VectorDeltaPhi * lambda + barphi;
			bartheta_prime = VectorDeltaTheta * lambda + bartheta;

			//f_new = calVectorb(square, barphi_prime);
			// Vectord_prime�CMatrixS_prime�̍X�V
			Vectorb_prime = calVectorb(square, barphi_prime, bartheta_prime);
			Vectorc_prime = calVectorc(square, barphi_prime, bartheta_prime);
			MatrixP_prime = calMatrixP(square, barphi_prime, bartheta_prime);
			MatrixQ_prime = calMatrixQ(square, barphi_prime, bartheta_prime);
			MatrixR_prime = calMatrixR(square, barphi_prime, bartheta_prime);
			for (int i = 0; i < NumberOfParticles; i++)
			{
				for (int j = 0; j < NumberOfParticles; j++) {
					// ��s��S_prime�̍쐬 ----------------------------------------
					// �s��S_prime�ɍs��P_prime���Z�b�g
					MatrixS_prime(3 * i, j) = MatrixP_prime(3 * i, j);
					MatrixS_prime(3 * i + 1, j) = MatrixP_prime(3 * i + 1, j);
					MatrixS_prime(3 * i + 2, j) = MatrixP_prime(3 * i + 2, j);

					// �s��S_prime�ɍs��Q_prime���Z�b�g
					MatrixS_prime(i + 3 * NumberOfParticles, j) = MatrixQ_prime(i, j);

					// �s��S_prime�ɍs��R_prime���Z�b�g
					MatrixS_prime(i + 3 * NumberOfParticles, j + NumberOfParticles) = MatrixR_prime(i, j);
				}

				// ��x�N�g��d_prime�̍쐬 ------------------------------------
				// �x�N�g��d_prime�Ƀx�N�g��b_prime���Z�b�g
				Vectord_prime(i) = Vectorb_prime(i);
				// �x�N�g��d_prime�Ƀx�N�g��c_prime���Z�b�g
				Vectord_prime(i + NumberOfParticles) = Vectorc_prime(i);
			}

			//Armijo
			f_prime = Vectord_prime;
			f = Vectord;
			nabla_f1 = eps1 * lambda * MatrixS.transpose() * VectorDelta;
			right_f = f + nabla_f1;

			//curvature
			nabla_f_prime = MatrixS_prime.transpose() * VectorDelta;
			nabla_f2 = eps2 * MatrixS.transpose() * VectorDelta;

			line_search_times++;
			std::cout << "lambda : " << lambda << "\n";
		}

		std::cout << "Result Lambda : " << lambda << "\n";
		std::cout << "\n";

		// ���W�̍X�V
		barphi += lambda * VectorDeltaPhi;

		// �̐ϕω����̍X�V
		bartheta += lambda * VectorDeltaTheta;

		looptimes++;
	}

	exportNorm(norm);

	// std::cout << barphi << "\n";
	return barphi;
}

Eigen::VectorXd Newton_H_onetime(Square square) {
	double NormVectorDeltaPhi = 1.0;

	vector <double>norm;

	int SquarePointsNumber = square.points.size();

	// �x�N�g���̑傫���� [3�~NumberOfParticles]�@�Ƃ���D
	Eigen::VectorXd barphi(3 * NumberOfParticles);
	// �x�N�g���̑傫���� [NumberOfParticles]�@�Ƃ���D
	Eigen::VectorXd bartheta(NumberOfParticles);
	// �x�N�g���̑傫����[3�~NumberOfParticles]�@�Ƃ���D�Q�Ɨp
	Eigen::VectorXd re_barphi(3 * NumberOfParticles);
	// �x�N�g���̑傫���� [3�~NumberOfParticles]�@�Ƃ���D�����T���@�ŗp����D
	Eigen::VectorXd barphi_prime(3 * NumberOfParticles);


	// ���W�̎擾
	for (int i = 0; i < SquarePointsNumber; i++) {
		barphi(3 * i) = square.points[i].position[0];
		barphi(3 * i + 1) = square.points[i].position[1];
		barphi(3 * i + 2) = square.points[i].position[2];
		// bartheta(i) = square.points[i].theta;
	}

	// �Q�Ɨp�̍��W�̎擾
	for (int i = 0; i < SquarePointsNumber; i++) {
		re_barphi(3 * i) = square.points[i].reference_position[0];
		re_barphi(3 * i + 1) = square.points[i].reference_position[1];
		re_barphi(3 * i + 2) = square.points[i].reference_position[2];
	}

	bartheta = 0.5 * calTheta(square, barphi);

	// std::cout << barphi << "\n";

	MatrixP_H = calMatrixP(square, barphi, bartheta); // �s��P���v�Z
	// exportMatrixP(MatrixP);
	Vectorb_H = calVectorb(square, barphi, bartheta); // �x�N�g��b���v�Z�i��ԗʂ𗝑z���j
	// std::cout << Vectorb_H << "\n";
	// exportVectorb_Convenient(Vectorb);

	MatrixR_H = MatrixP_H * MatrixP_H.transpose();
	Vectorc_H = Vectorb_H.transpose() * MatrixP_H.transpose();

	Eigen::FullPivLU<Eigen::MatrixXd> LU(MatrixR_H);
	Eigen::VectorXd VectorDeltaPhi = LU.solve(Vectorc_H);


	NormVectorDeltaPhi = VectorDeltaPhi.norm();
	norm.emplace_back(NormVectorDeltaPhi);
	std::cout << "NormVectorDeltaPhi = " << NormVectorDeltaPhi << "\n";
	// exportVectorDeltaPhi(VectorDeltaPhi);

	//Line Search with Wolfe
	double sigma = 0.5;
	double eps1 = 0.01;
	double eps2 = 0.1;
	double lambda = 1.0;
	int line_search_times = 1;

	barphi_prime = VectorDeltaPhi * lambda + barphi;

	//Armijo
	Eigen::VectorXd f_prime = calVectorb_Convenient(square, barphi_prime, re_barphi, bartheta);;
	Eigen::VectorXd f = Vectorb_H;
	Eigen::VectorXd nabla_f1 = eps1 * lambda * MatrixP_H.transpose() * VectorDeltaPhi;
	Eigen::VectorXd right_f = f + nabla_f1;

	//curvature
	Eigen::VectorXd nabla_f_prime = calMatrixP(square, barphi_prime, bartheta).transpose() * VectorDeltaPhi;
	Eigen::VectorXd nabla_f2 = eps2 * MatrixP_H.transpose() * VectorDeltaPhi;

	
	while (f_prime.norm() > right_f.norm() || abs(nabla_f_prime.norm()) < abs(nabla_f2.norm())) {
		// �ɁCbarphi_prime�̍X�V
		lambda = sigma * lambda;
		barphi_prime = VectorDeltaPhi * lambda + barphi;

		//Armijo
		f_prime = calVectorb(square, barphi_prime, bartheta);
		nabla_f1 = eps1 * lambda * MatrixP_H.transpose() * VectorDeltaPhi;
		right_f = f + nabla_f1;

		//curvature
		nabla_f_prime = calMatrixP(square, barphi_prime, bartheta).transpose() * VectorDeltaPhi;
		nabla_f2 = eps2 * MatrixP_H.transpose() * VectorDeltaPhi;

		// line_search_times++;
	}


	std::cout << "lambda : " << lambda << "\n";
	


	// ���W�̍X�V
	barphi += lambda * VectorDeltaPhi;

	return barphi;
}

Eigen::VectorXd Newton_Convenient(Square square) {
	double NormVectorDeltaPhi = 1.0;
	double NormVectorDeltaTheta = 1.0;
	double NormVectorDelta = 1.0;

	vector <double>norm;
	
	int SquarePointsNumber = square.points.size();

	// �x�N�g���̑傫���� [3�~NumberOfParticles]�@�Ƃ���D
	Eigen::VectorXd barphi(3 * NumberOfParticles);
	// �x�N�g���̑傫���� [NumberOfParticles]�@�Ƃ���D
	Eigen::VectorXd bartheta(NumberOfParticles);
	// �x�N�g���̑傫����[3�~NumberOfParticles]�@�Ƃ���D�Q�Ɨp
	Eigen::VectorXd re_barphi(3 * NumberOfParticles);
	// �x�N�g���̑傫���� [3�~NumberOfParticles]�@�Ƃ���D�����T���@�ŗp����D
	Eigen::VectorXd barphi_prime(3 * NumberOfParticles);
	// �x�N�g���̑傫���� [3�~NumberOfParticles]�@�Ƃ���D�����T���@�ŗp����D
	Eigen::VectorXd bartheta_prime(NumberOfParticles);

	// ���W�̎擾
	for (int i = 0; i < SquarePointsNumber; i++) {
		barphi(3 * i) = square.points[i].position[0];
		barphi(3 * i + 1) = square.points[i].position[1];
		barphi(3 * i + 2) = square.points[i].position[2];
		bartheta(i) = square.points[i].theta;
	}

	// �Q�Ɨp�̍��W�̎擾
	for (int i = 0; i < SquarePointsNumber; i++) {
		re_barphi(3 * i) = square.points[i].reference_position[0];
		re_barphi(3 * i + 1) = square.points[i].reference_position[1];
		re_barphi(3 * i + 2) = square.points[i].reference_position[2];
	}

	// std::cout << barphi << "\n";

	while(NormVectorDelta > 1.0e-6) {
		MatrixP = calMatrixP(square, barphi, bartheta); // �s��P���v�Z
		exportMatrixP(MatrixP);
		MatrixQ = calMatrixQ(square, barphi, bartheta); // �s��Q���v�Z
		exportMatrixQ(MatrixQ);
		MatrixR = calMatrixR(square, barphi, bartheta); // �s��R���v�Z
		exportMatrixR(MatrixR);
		Vectorb = calVectorb_Convenient(square, barphi, re_barphi, bartheta); // �x�N�g��b���v�Z�i��ԗʂ𗝑z���j
		exportVectorb_Convenient(Vectorb);
		Vectorc = calVectorc(square, barphi, bartheta); // �x�N�g��c���v�Z
		exportVectorc(Vectorc);
		// std::cout << "Vectorc :" << "\n" << Vectorc << "\n";


		for (int i = 0; i < NumberOfParticles; i++) {
			for (int j = 0; j < NumberOfParticles; j++) {
				// ��s��S�̍쐬 ----------------------------------------
				// �s��S�ɍs��P���Z�b�g
				MatrixS(3 * i, j) = MatrixP(3 * i, j);
				MatrixS(3 * i + 1, j) = MatrixP(3 * i + 1, j);
				MatrixS(3 * i + 2, j) = MatrixP(3 * i + 2, j);

				// �s��S�ɍs��Q���Z�b�g
				MatrixS(i + 3 * NumberOfParticles, j) = MatrixQ(i, j);

				// �s��S�ɍs��R���Z�b�g
				MatrixS(i + 3 * NumberOfParticles, j + NumberOfParticles) = MatrixR(i, j);
			}
			// std::cout << MatrixS << std::endl;

			// ��x�N�g��d�̍쐬 ------------------------------------
			// �x�N�g��d�Ƀx�N�g��b���Z�b�g
			Vectord(i) = Vectorb(i);

			// �x�N�g��d�Ƀx�N�g��c���Z�b�g
			Vectord(i + NumberOfParticles) = Vectorc(i);

			
		}

		// �s��S�̏o��
		exportMatrixS(MatrixS);
		// �x�N�g��d�̏o��
		exportVectord(Vectord);
		//std::cout << "Vectord :" << "\n" << Vectord << "\n";

		MatrixT = calMatrixT(MatrixS); // �s��T���v�Z
		exportMatrixT(MatrixT);
		Vectore = calVectore(Vectord, MatrixS); // �x�N�g��e���v�Z
		exportVectore(Vectore);

		Eigen::FullPivLU<Eigen::MatrixXd> LU(MatrixT);
		Eigen::VectorXd VectorDelta = LU.solve(Vectore);
		// std::cout << "DeltaVector :" << "\n" << VectorDelta << "\n";
		NormVectorDelta = VectorDelta.norm();
		norm.emplace_back(NormVectorDelta);
		std::cout << "NormDeltaVector = " << NormVectorDelta << "\n";

		// phi��theta�̕���
		for (int i = 0; i < NumberOfParticles; i++) {
			// Phi
			VectorDeltaPhi(3 * i) = VectorDelta(3 * i);
			VectorDeltaPhi(3 * i + 1) = VectorDelta(3 * i + 1);
			VectorDeltaPhi(3 * i + 2) = VectorDelta(3 * i + 2);

			// Theta
			VectorDeltaTheta(i) = VectorDelta(i + 3 * NumberOfParticles);
		}

		// DeltaPhi��DeltaTheta�̏o��
		// std::cout << "DeltaVectorPhi :" << "\n" << VectorDeltaPhi << "\n";
		// std::cout << "DeltaVectorTheta :" << "\n" << VectorDeltaTheta << "\n";

		exportVectorDeltaPhi(VectorDeltaPhi);
		exportVectorDeltaTheta(VectorDeltaTheta);

		NormVectorDeltaPhi = VectorDeltaPhi.norm();
		// std::cout << "NormDeltaVectorPhi = " << NormVectorDeltaPhi << "\n";
		NormVectorDeltaTheta = VectorDeltaTheta.norm();
		// std::cout << "NormDeltaVectorTheta = " << NormVectorDeltaTheta << "\n";

		// std::cout << "NormDeltaVector : " << VectorDelta.norm() << "\n";
		
		
		//Line Search with Wolfe
		double sigma = 0.5;
		double eps1 = 0.01;
		double eps2 = 0.1;
		double lambda = 1.0;
		int line_search_times = 1;

		barphi_prime = VectorDeltaPhi * lambda + barphi;
		bartheta_prime = VectorDeltaTheta * lambda + bartheta;

		Eigen::VectorXd Vectorb_prime = -1 * calVectorb_Convenient(square, barphi_prime, re_barphi, bartheta_prime);
		Eigen::VectorXd Vectorc_prime = -1 * calVectorc(square, barphi_prime, bartheta_prime);
		Eigen::VectorXd Vectord_prime(2 * NumberOfParticles);
		Eigen::MatrixXd MatrixP_prime = calMatrixP(square, barphi_prime, bartheta_prime);
		Eigen::MatrixXd MatrixQ_prime = calMatrixQ(square, barphi_prime, bartheta_prime);
		Eigen::MatrixXd MatrixR_prime = calMatrixR(square, barphi_prime, bartheta_prime);
		Eigen::MatrixXd MatrixS_prime(4 * NumberOfParticles, 2 * NumberOfParticles);


		for (int i = 0; i < NumberOfParticles; i++) 
		{
			for (int j = 0; j < NumberOfParticles; j++) {
				// ��s��S_prime�̍쐬 ----------------------------------------
				// �s��S_prime�ɍs��P_prime���Z�b�g
				MatrixS_prime(3 * i, j) = MatrixP_prime(3 * i, j);
				MatrixS_prime(3 * i + 1, j) = MatrixP_prime(3 * i + 1, j);
				MatrixS_prime(3 * i + 2, j) = MatrixP_prime(3 * i + 2, j);

				// �s��S_prime�ɍs��Q_prime���Z�b�g
				MatrixS_prime(i + 3 * NumberOfParticles, j) = MatrixQ_prime(i, j);

				// �s��S_prime�ɍs��R_prime���Z�b�g
				MatrixS_prime(i + 3 * NumberOfParticles, j + NumberOfParticles) = MatrixR_prime(i, j);
			}

			// ��x�N�g��d_prime�̍쐬 ------------------------------------
			// �x�N�g��d_prime�Ƀx�N�g��b_prime���Z�b�g
			Vectord_prime(i) = Vectorb_prime(i);
			// �x�N�g��d_prime�Ƀx�N�g��c_prime���Z�b�g
			Vectord_prime(i + NumberOfParticles) = Vectorc_prime(i);
		}

		
		//Armijo
		Eigen::VectorXd f_prime = Vectord_prime;
		Eigen::VectorXd f = Vectord;
		Eigen::VectorXd nabla_f1 = eps1 * lambda * MatrixS.transpose() * VectorDelta;
		Eigen::VectorXd right_f = f + nabla_f1;

		//curvature
		Eigen::VectorXd nabla_f_prime = MatrixS_prime.transpose() * VectorDelta;
		Eigen::VectorXd nabla_f2 = eps2 * MatrixS.transpose() * VectorDelta;
		
		while (f_prime.norm() > right_f.norm() || nabla_f_prime.norm() < nabla_f2.norm()) {
			// �ɁCbarphi_prime�Cbartheta_prime�̍X�V
			// lambda = 2.0 / (2 + line_search_times);
			lambda = sigma * lambda;
			// std::cout << line_search_times << std::endl;
			barphi_prime = VectorDeltaPhi * lambda + barphi;
			bartheta_prime = VectorDeltaTheta * lambda + bartheta;

			//f_new = calVectorb(square, barphi_prime);
			// Vectord_prime�CMatrixS_prime�̍X�V
			Vectorb_prime = calVectorb_Convenient(square, barphi_prime, re_barphi, bartheta_prime);
			Vectorc_prime = calVectorc(square, barphi_prime, bartheta_prime);
			MatrixP_prime = calMatrixP(square, barphi_prime, bartheta_prime);
			MatrixQ_prime = calMatrixQ(square, barphi_prime, bartheta_prime);
			MatrixR_prime = calMatrixR(square, barphi_prime, bartheta_prime);
			for (int i = 0; i < NumberOfParticles; i++)
			{
				for (int j = 0; j < NumberOfParticles; j++) {
					// ��s��S_prime�̍쐬 ----------------------------------------
					// �s��S_prime�ɍs��P_prime���Z�b�g
					MatrixS_prime(3 * i, j) = MatrixP_prime(3 * i, j);
					MatrixS_prime(3 * i + 1, j) = MatrixP_prime(3 * i + 1, j);
					MatrixS_prime(3 * i + 2, j) = MatrixP_prime(3 * i + 2, j);

					// �s��S_prime�ɍs��Q_prime���Z�b�g
					MatrixS_prime(i + 3 * NumberOfParticles, j) = MatrixQ_prime(i, j);

					// �s��S_prime�ɍs��R_prime���Z�b�g
					MatrixS_prime(i + 3 * NumberOfParticles, j + NumberOfParticles) = MatrixR_prime(i, j);
				}

				// ��x�N�g��d_prime�̍쐬 ------------------------------------
				// �x�N�g��d_prime�Ƀx�N�g��b_prime���Z�b�g
				Vectord_prime(i) = Vectorb_prime(i);
				// �x�N�g��d_prime�Ƀx�N�g��c_prime���Z�b�g
				Vectord_prime(i + NumberOfParticles) = Vectorc_prime(i);
			}

			//Armijo
			f_prime = Vectord_prime;
			f = Vectord;
			nabla_f1 = eps1 * lambda * MatrixS.transpose() * VectorDelta;
			right_f = f + nabla_f1;

			//curvature
			nabla_f_prime = MatrixS_prime.transpose() * VectorDelta;
			nabla_f2 = eps2 * MatrixS.transpose() * VectorDelta;

			line_search_times++;
		}
		

		std::cout << "lambda : " << lambda << "\n";

		// ���W�̍X�V
		barphi += lambda * VectorDeltaPhi;

		// �̐ϕω����̍X�V
		bartheta += lambda * VectorDeltaTheta;
	}

	// std::cout << barphi.size() << "\n";

	/*
	for (int i = 0; i < norm.size(); i++) {
		std::cout << norm[i] << std::endl;
	}
	*/

	exportNorm(norm);

	return barphi;
}

// ���̒l���Œ�(1�Ȃ�)
Eigen::VectorXd Newton_H(Square square) {

	double NormVectorDeltaPhi = 1.0;

	vector <double>norm;

	int SquarePointsNumber = square.points.size();

	// �x�N�g���̑傫���� [3�~NumberOfParticles]�@�Ƃ���D
	Eigen::VectorXd barphi(3 * NumberOfParticles);
	// �x�N�g���̑傫���� [NumberOfParticles]�@�Ƃ���D
	Eigen::VectorXd bartheta(NumberOfParticles);
	// �x�N�g���̑傫����[3�~NumberOfParticles]�@�Ƃ���D�Q�Ɨp
	Eigen::VectorXd re_barphi(3 * NumberOfParticles);
	// �x�N�g���̑傫���� [3�~NumberOfParticles]�@�Ƃ���D�����T���@�ŗp����D
	Eigen::VectorXd barphi_prime(3 * NumberOfParticles);


	// ���W�̎擾
	for (int i = 0; i < SquarePointsNumber; i++) {
		barphi(3 * i) = square.points[i].position[0];
		barphi(3 * i + 1) = square.points[i].position[1];
		barphi(3 * i + 2) = square.points[i].position[2];
		// bartheta(i) = square.points[i].theta;
	}

	// �Q�Ɨp�̍��W�̎擾
	for (int i = 0; i < SquarePointsNumber; i++) {
		re_barphi(3 * i) = square.points[i].reference_position[0];
		re_barphi(3 * i + 1) = square.points[i].reference_position[1];
		re_barphi(3 * i + 2) = square.points[i].reference_position[2];
	}

	bartheta = 0.5 * calTheta(square, barphi);

	// std::cout << barphi << "\n";

	int looptimes = 0;

	while (NormVectorDeltaPhi > 1.0e-6) {
		MatrixP_H = calMatrixP(square, barphi, bartheta); // �s��P���v�Z
		// exportMatrixP(MatrixP);
		Vectorb_H = calVectorb(square, barphi, bartheta); // �x�N�g��b���v�Z�i��ԗʂ𗝑z���j
		// std::cout << Vectorb_H << "\n";
		// exportVectorb_Convenient(Vectorb);

		if (looptimes == 0) {
			exportMatrixP(MatrixP_H);
		}

		MatrixR_H = MatrixP_H * MatrixP_H.transpose();
		Vectorc_H = Vectorb_H.transpose() * MatrixP_H.transpose();

		Eigen::FullPivLU<Eigen::MatrixXd> LU(MatrixR_H);
		Eigen::VectorXd VectorDeltaPhi = LU.solve(Vectorc_H);
		

		NormVectorDeltaPhi = VectorDeltaPhi.norm();
		norm.emplace_back(NormVectorDeltaPhi);
		std::cout << "NormVectorDeltaPhi = " << NormVectorDeltaPhi << "\n";
		exportVectorDeltaPhi(VectorDeltaPhi);

		//Line Search with Wolfe
		double sigma = 0.5;
		double eps1 = 0.01;
		double eps2 = 0.1;
		double lambda = 1.0;
		int line_search_times = 1;

		barphi_prime = VectorDeltaPhi * lambda + barphi;

		//Armijo
		Eigen::VectorXd f_prime = calVectorb_Convenient(square, barphi_prime, re_barphi, bartheta);;
		Eigen::VectorXd f = Vectorb_H;
		Eigen::VectorXd nabla_f1 = eps1 * lambda * MatrixP_H.transpose() * VectorDeltaPhi;
		Eigen::VectorXd right_f = f + nabla_f1;

		//curvature
		Eigen::VectorXd nabla_f_prime = calMatrixP(square, barphi_prime, bartheta).transpose() * VectorDeltaPhi;
		Eigen::VectorXd nabla_f2 = eps2 * MatrixP_H.transpose() * VectorDeltaPhi;

		while (f_prime.norm() > right_f.norm() || abs(nabla_f_prime.norm()) < abs(nabla_f2.norm()) ) {
			// �ɁCbarphi_prime�̍X�V
			lambda = sigma * lambda;
			barphi_prime = VectorDeltaPhi * lambda + barphi;

			//Armijo
			f_prime = calVectorb(square, barphi_prime, bartheta);
			nabla_f1 = eps1 * lambda * MatrixP_H.transpose() * VectorDeltaPhi;
			right_f = f + nabla_f1;

			//curvature
			nabla_f_prime = calMatrixP(square, barphi_prime, bartheta).transpose() * VectorDeltaPhi;
			nabla_f2 = eps2 * MatrixP_H.transpose() * VectorDeltaPhi;

			// line_search_times++;
			std::cout << "lambda : " << lambda << "\n";
		}

		std::cout << "Result Lambda : " << lambda << "\n";
		std::cout << "\n";
		

		// ���W�̍X�V
		barphi += lambda * VectorDeltaPhi;

		looptimes++;
	}

	// std::cout << barphi.size() << "\n";

	exportNorm(norm);

	return barphi;
}

void calInterpolationSophia() 
{
	Kernel_Array ka = createKernelArray();

	for (int l = 0; l < dimensions; l++) {
		for (int m = 0; m < dimensions; m++) {
			for (int n = 0; n < dimensions; n++) {

				/*
				* �e�O���b�h�m�[�h�̔������ꂽ��Ԋ֐���x,y,z�̂ǂ̎����Ŕ�������Ă��邩
				* x = 1, y = 2, z = 3
				* ki = kernel index
				*/
				int ki_i = l + 1;  // grid i �ɂ����Ĕ�������kernel
				int	ki_j = m + 1;  // grid j �ɂ����Ĕ�������kernel
				int	ki_k = n + 1;  // grid k �ɂ����Ĕ�������kernel

				for (int i = 0; i < dimensions; i++) {
					for (int j = 0; j < dimensions; j++) {
						for (int k = 0; k < dimensions; k++) {

							// grid node's difference
							// -1 �` 1
							int gn_i = i - 1;
							int gn_j = j - 1;
							int gn_k = k - 1;

							// case x
							ka = setCoefficientSophiaX(ka, ki_i, ki_j, ki_k, gn_i, gn_j, gn_k);
							SophiaX[l][m][n][i][j][k] = productSophiaX(ka);


							// case y
							ka = setCoefficientSophiaY(ka, ki_i, ki_j, ki_k, gn_i, gn_j, gn_k);
							SophiaY[l][m][n][i][j][k] = productSophiaY(ka);


							// case z
							ka = setCoefficientSophiaZ(ka, ki_i, ki_j, ki_k, gn_i, gn_j, gn_k);
							SophiaZ[l][m][n][i][j][k] = productSophiaZ(ka);

						}
					}
				}

			}
		}
	}

	// ���ʂ��o��
	exportInterpolationSophia(SophiaX, SophiaY, SophiaZ);

}

void calInterpolationChloe() {
	Kernel_Array ka = createKernelArray();

	// grid node
	for (int ix = 0; ix < dimensions; ix++) {
		for (int iy = 0; iy < dimensions; iy++) {
			for (int iz = 0; iz < dimensions; iz++) {

				// -1 �` 1
				int gn_ix = ix - 1;
				int gn_iy = iy - 1;
				int gn_iz = iz - 1;

				// std::cout << "grid node index : ";
				// std::cout << "x " << gn_ix << ",y " << gn_iy << ",z " << gn_iz << "\n";

				ka = setCoefficientChloeX(ka, gn_ix);
				ka = setCoefficientChloeY(ka, gn_iy);
				ka = setCoefficientChloeZ(ka, gn_iz);

				double el = productChloeX(ka) * productChloeY(ka) * productChloeZ(ka);

				// std::cout << "result : " << el << "\n" << "\n";

				Chloe[ix][iy][iz] = { el };
			}
		}
	}

	// ���ʂ��o��
	exportInterpolationChloe(Chloe);
}

void calInterpolationAria() 
{
	Kernel_Array ka = createKernelArray();

	// grid node
	for (int ix = 0; ix < dimensions; ix++) {
		for (int iy = 0; iy < dimensions; iy++) {
			for (int iz = 0; iz < dimensions; iz++) {

				// -1 �` 1
				int gn_ix = ix - 1;
				int gn_iy = iy - 1;
				int gn_iz = iz - 1;

				// std::cout << "grid node index : ";
				// std::cout << "x " << gn_ix << ",y " << gn_iy << ",z " << gn_iz << "\n";

				for (int jx = 0; jx < dimensions; jx++) {
					for (int jy = 0; jy < dimensions; jy++) {
						for (int jz = 0; jz < dimensions; jz++) {

							// -1 �` 1
							int gn_jx = jx - 1;
							int gn_jy = jy - 1;
							int gn_jz = jz - 1;

							// std::cout << "grid node index : ";
							// std::cout << "x " << gn_jx << ",y " << gn_jy << ",z " << gn_jz << "\n";

							ka = setCoefficientAriaX(ka, gn_ix, gn_jx);
							ka = setCoefficientAriaY(ka, gn_iy, gn_jy);
							ka = setCoefficientAriaZ(ka, gn_iz, gn_jz);

							double el = productAriaX(ka) * productAriaY(ka) * productAriaZ(ka);

							/*
							if (ix == 0 && iy == 0 && iz == 0 && jx == 0 && jy == 0 && jz == 0) {
								std::cout << "result : " << productAriaX(ka) << "\n";
								std::cout << "result : " << productAriaY(ka) << "\n";
								std::cout << "result : " << productAriaZ(ka) << "\n";
							}
							*/

							// std::cout << "result : " << el << "\n" << "\n";

							Aria[ix][iy][iz][jx][jy][jz] = { el };

							// std::cout << "result : " << Aria[ix][iy][iz][jx][jy][jz] << "\n" << "\n";
						}
					}
				}

			}
		}
	}

	// ���ʂ��o��
	exportInterpolationAria(Aria);
}

void calInterpolationMia()
{
	Kernel_Array ka = createKernelArray();

	for (int p = 0; p < dimensions; p++) {
		for (int q = 0; q < dimensions; q++) {
			for (int r = 0; r < dimensions; r++) {

				/*
				* �e�O���b�h�m�[�h�̔������ꂽ��Ԋ֐���x,y,z�̂ǂ̎����Ŕ�������Ă��邩
				* x = 1, y = 2, z = 3
				* ki = kernel index
				*/
				int ki_i = p + 1;  // grid i �ɂ����Ĕ�������kernel
				int	ki_j = q + 1;  // grid j �ɂ����Ĕ�������kernel
				int	ki_k = r + 1;  // grid k �ɂ����Ĕ�������kernel

				// i, j, k, l��1,-1,0��3�����l��
				for (int i = 0; i < dimensions; i++) {
					for (int j = 0; j < dimensions; j++) {
						for (int k = 0; k < dimensions; k++) {
							for (int l = 0; l < dimensions; l++) {

								// grid node's difference
								// -1 �` 1
								int gn_i = i - 1;
								int gn_j = j - 1;
								int gn_k = k - 1;
								int gn_l = l - 1;

								// case x
								ka = setCoefficientMiaX(ka, ki_i, ki_j, ki_k, gn_i, gn_j, gn_k, gn_l);
								MiaX[p][q][r][i][j][k][l] = productMiaX(ka);

								
								// case y
								ka = setCoefficientMiaY(ka, ki_i, ki_j, ki_k, gn_i, gn_j, gn_k, gn_l);
								MiaY[p][q][r][i][j][k][l] = productMiaY(ka);


								// case z
								ka = setCoefficientMiaZ(ka, ki_i, ki_j, ki_k, gn_i, gn_j, gn_k, gn_l);
								MiaZ[p][q][r][i][j][k][l] = productMiaZ(ka);
								
							}
						}
					}
				}

			}
		}
	}

	
	// ���ʂ��o��
	exportInterpolationMia(MiaX, MiaY, MiaZ);
	
}

// For Equation of Theta
Eigen::VectorXd calTheta(Square square, Eigen::VectorXd phi) {
	Eigen::MatrixXd W(NumberOfParticles, NumberOfParticles);
	bool FlagCalGrid;

	// calculation Matrx W
	for (int xi = 0; xi < NumberOfParticles; xi++) {
		for (int i = 0; i < NumberOfParticles; i++) {

			FlagCalGrid = false;

			Eigen::Vector3i grid_xi = FlatToGrid(xi);
			Eigen::Vector3i grid_i = FlatToGrid(i);
			Eigen::Vector3i i_minus_xi = grid_i - grid_xi;

			W(i, xi) = 0.0;

			// i - xi�� -1�`1�͈̔͂ł��邩�̔���
			// i - xi
			if ((abs(i_minus_xi[0]) <= 1) && (abs(i_minus_xi[1]) <= 1) && (abs(i_minus_xi[2]) <= 1)) {
				FlagCalGrid = true;
				//std::cout << "true\n";
			}

			if (FlagCalGrid) {
				W(i, xi) = Chloe[i_minus_xi[0] + 1][i_minus_xi[1] + 1][i_minus_xi[2] + 1];
			}

		}
	}

	/*
	* �s��̏o��
	for (int xi = 0; xi < NumberOfParticles; xi++) {
		for (int i = 0; i < NumberOfParticles; i++) {
			std::cout << W(i, xi);
		}
		std::cout << std::endl;
	}
	*/

	
	// Vector d
	Eigen::VectorXd D = Eigen::VectorXd::Zero(NumberOfParticles);

	for (int xi = 0; xi < NumberOfParticles; xi++) {
		Eigen::Vector3i grid_xi = FlatToGrid(xi);
		
		for (int i = 0; i < NumberOfParticles; i++) {
			Eigen::Vector3i grid_i = FlatToGrid(i);
			Eigen::Vector3i i_minus_xi = grid_i - grid_xi;
			double bar_phi_i_x = phi(3 * i);
			double bar_phi_i_y = phi(3 * i + 1);
			double bar_phi_i_z = phi(3 * i + 2);
			
			for (int j = 0; j < NumberOfParticles; j++) {
				Eigen::Vector3i grid_j = FlatToGrid(j);
				Eigen::Vector3i j_minus_xi = grid_j - grid_xi;
				double bar_phi_j_x = phi(3 * j);
				double bar_phi_j_y = phi(3 * j + 1);
				double bar_phi_j_z = phi(3 * j + 2);
				// std::cout << D(j) << std::endl;

				for (int k = 0; k < NumberOfParticles; k++) {
					Eigen::Vector3i grid_k = FlatToGrid(k);
					Eigen::Vector3i k_minus_xi = grid_k - grid_xi;
					double bar_phi_k_x = phi(3 * k);
					double bar_phi_k_y = phi(3 * k + 1);
					double bar_phi_k_z = phi(3 * k + 2);
					// std::cout << D(k) << std::endl;
					FlagCalGrid = false;

					// i - xi, j - xi, k - xi �� -1�`1�͈̔͂ł��邩�̔���
					//i - xi
					if ((abs(i_minus_xi[0]) <= 1) && (abs(i_minus_xi[1]) <= 1) && (abs(i_minus_xi[2]) <= 1)) {
						//j - xi
						if ((abs(j_minus_xi[0]) <= 1) && (abs(j_minus_xi[1]) <= 1) && (abs(j_minus_xi[2]) <= 1)) {
							//k - xi
							if ((abs(k_minus_xi[0]) <= 1) && (abs(k_minus_xi[1]) <= 1) && (abs(k_minus_xi[2]) <= 1)) {
								FlagCalGrid = true;
								//std::cout << "true\n";
							}
						}
					}

					// ����Z�Ȃ�v�Z
					if (FlagCalGrid) {
						double w1w2w3 = SophiaX[0][1][2][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* SophiaY[0][1][2][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* SophiaZ[0][1][2][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];


						double Phi1 =
							bar_phi_i_x * (bar_phi_j_y * bar_phi_k_z - bar_phi_j_z * bar_phi_k_y)
							+ bar_phi_i_y * (bar_phi_j_z * bar_phi_k_x - bar_phi_j_x * bar_phi_k_z)
							+ bar_phi_i_z * (bar_phi_j_x * bar_phi_k_y - bar_phi_j_y * bar_phi_k_x);

						// i, j, k�Řa���Ƃ�
						D(xi) += Phi1 * w1w2w3;
					}

				}
			}
		}
	}

	// std::cout << D.size() << std::endl;
	// std::cout << W.rows() << std::endl;
	// std::cout << W.cols() << std::endl;

	// std::cout << W.determinant() << std::endl;
	
	Eigen::PartialPivLU<Eigen::MatrixXd> LU(W * pow(square.dx, 3));
	Eigen::VectorXd Theta = LU.solve(D);

	/*
	std::cout << "Theta value:" << std::endl;

	for (int i = 0; i < Theta.size(); i++) {
		std::cout << Theta(i) << std::endl;
	}
	
	*/
	
	return Theta;
}


Eigen::MatrixXd calMatrixP(Square square, Eigen::VectorXd phi, Eigen::VectorXd theta)
{
	int las = 0;
	bool FlagCalGrid;
	Eigen::MatrixXd OutputM(3 * NumberOfParticles, NumberOfParticles);

	//calculate P matrix
	for (int xi = 0; xi < NumberOfParticles; xi++) {
		for (int k = 0; k < NumberOfParticles; k++) {

			Eigen::Vector3i grid_k = FlatToGrid(k);
			Eigen::Vector3i grid_xi = FlatToGrid(xi);
			Eigen::Vector3i k_minus_xi = grid_k - grid_xi;

			OutputM(3 * k, xi) = 0.0;
			OutputM(3 * k + 1, xi) = 0.0;
			OutputM(3 * k + 2, xi) = 0.0;

			double m_x[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
			double m_y[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
			double m_z[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

			double Q1 = 0.0;
			double Q2 = 0.0;
			double Q3 = 0.0;

			int n = 0;

			for (int i = 0; i < NumberOfParticles; i++) {
				for (int j = 0; j < NumberOfParticles; j++) {

					FlagCalGrid = false;

					Eigen::Vector3i grid_i = FlatToGrid(i);
					Eigen::Vector3i grid_j = FlatToGrid(j);

					Eigen::Vector3i i_minus_xi = grid_i - grid_xi;
					Eigen::Vector3i j_minus_xi = grid_j - grid_xi;


					// \Phi = \bar{\varphi2}^i * \bar{\varphi3}^j - \bar{\varphi3}^i * \bar{\varphi2}^j
					double Phi_2_3 = (phi(3 * i + 1) * phi(3 * j + 2));
					double Phi_3_2 = (phi(3 * i + 2) * phi(3 * j + 1));


					// \Phi = \bar{\varphi3}^i * \bar{\varphi1}^j - \bar{\varphi1}^i * \bar{\varphi3}^j
					double Phi_3_1 = (phi(3 * i + 2) * phi(3 * j));
					double Phi_1_3 = (phi(3 * i) * phi(3 * j + 2));


					// \Phi = \bar{\varphi1}^i * \bar{\varphi2}^j - \bar{\varphi2}^i * \bar{\varphi1}^j
					double Phi_1_2 = (phi(3 * i) * phi(3 * j + 1));
					double Phi_2_1 = (phi(3 * i + 1) * phi(3 * j));


					// i - xi, j - xi, k - xi �� -1�`1�͈̔͂ł��邩�̔���
					//i - xi
					if ((abs(i_minus_xi[0]) <= 1) && (abs(i_minus_xi[1]) <= 1) && (abs(i_minus_xi[2]) <= 1)) {
						//j - xi
						if ((abs(j_minus_xi[0]) <= 1) && (abs(j_minus_xi[1]) <= 1) && (abs(j_minus_xi[2]) <= 1)) {
							//k - xi
							if ((abs(k_minus_xi[0]) <= 1) && (abs(k_minus_xi[1]) <= 1) && (abs(k_minus_xi[2]) <= 1)) {
								FlagCalGrid = true;
								//std::cout << "true\n";
							}
						}
					}

					// ����Z�Ȃ�v�Z
					if (FlagCalGrid) {
						
						/*
						double w2w3w1 = SophiaX[1][2][0][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* SophiaY[1][2][0][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* SophiaZ[1][2][0][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						double w1w3w2 = SophiaX[0][2][1][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* SophiaY[0][2][1][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* SophiaZ[0][2][1][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						double w1w2w3 = SophiaX[0][1][2][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* SophiaY[0][1][2][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* SophiaZ[0][1][2][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						m_x[0] += Phi_2_3 * w2w3w1;
						m_x[1] += Phi_3_2 * w2w3w1;
						m_x[2] += Phi_3_2 * w1w3w2;
						m_x[3] += Phi_2_3 * w1w3w2;
						m_x[4] += Phi_2_3 * w1w2w3;
						m_x[5] += Phi_3_2 * w1w2w3;

						m_y[0] += Phi_3_1 * w2w3w1;
						m_y[1] += Phi_1_3 * w2w3w1;
						m_y[2] += Phi_1_3 * w1w3w2;
						m_y[3] += Phi_3_1 * w1w3w2;
						m_y[4] += Phi_3_1 * w1w2w3;
						m_y[5] += Phi_1_3 * w1w2w3;

						m_z[0] += Phi_1_2 * w2w3w1;
						m_z[1] += Phi_2_1 * w2w3w1;
						m_z[2] += Phi_2_1 * w1w3w2;
						m_z[3] += Phi_1_2 * w1w3w2;
						m_z[4] += Phi_1_2 * w1w2w3;
						m_z[5] += Phi_2_1 * w1w2w3;
						*/

						
						vector<Eigen::Vector3i> v = { i_minus_xi, j_minus_xi, k_minus_xi };

						double w2w3w1 = Riemann_Sum_for_Sophia(v, { 1, 2, 0 }, square.dx);

						double w1w3w2 = Riemann_Sum_for_Sophia(v, { 0, 2, 1 }, square.dx);

						double w1w2w3 = Riemann_Sum_for_Sophia(v, { 0, 1, 2 }, square.dx);

						m_x[0] += Phi_2_3 * w2w3w1;
						m_x[1] += Phi_3_2 * w2w3w1;
						m_x[2] += Phi_3_2 * w1w3w2;
						m_x[3] += Phi_2_3 * w1w3w2;
						m_x[4] += Phi_2_3 * w1w2w3;
						m_x[5] += Phi_3_2 * w1w2w3;

						m_y[0] += Phi_3_1 * w2w3w1;
						m_y[1] += Phi_1_3 * w2w3w1;
						m_y[2] += Phi_1_3 * w1w3w2;
						m_y[3] += Phi_3_1 * w1w3w2;
						m_y[4] += Phi_3_1 * w1w2w3;
						m_y[5] += Phi_1_3 * w1w2w3;

						m_z[0] += Phi_1_2 * w2w3w1;
						m_z[1] += Phi_2_1 * w2w3w1;
						m_z[2] += Phi_2_1 * w1w3w2;
						m_z[3] += Phi_1_2 * w1w3w2;
						m_z[4] += Phi_1_2 * w1w2w3;
						m_z[5] += Phi_2_1 * w1w2w3;
						
					}


					
				}
			}



			/*
			if (abs(m_x[1]) > 1.0e-6) {
				std::cout << "xi : " << xi << ", k : " << k << ", m_x[1] : " << m_x[1] << "\n";
			}
			*/

			/*
			OutputM(3 * k, xi) = Q1;
			OutputM(3 * k + 1, xi) = Q2;
			OutputM(3 * k + 2, xi) = Q3;
			*/


			OutputM(3 * k, xi) = m_x[0] - m_x[1] + m_x[2] - m_x[3] + m_x[4] - m_x[5];
			OutputM(3 * k + 1, xi) = m_y[0] - m_y[1] + m_y[2] - m_y[3] + m_y[4] - m_y[5];
			OutputM(3 * k + 2, xi) = m_z[0] - m_z[1] + m_z[2] - m_z[3] + m_z[4] - m_z[5];

		}
	}

	return OutputM;
}

Eigen::MatrixXd calMatrixQ(Square square, Eigen::VectorXd phi, Eigen::VectorXd theta)
{
	int las = 0;
	bool FlagCalGrid;
	Eigen::MatrixXd OutputQ(NumberOfParticles, NumberOfParticles);

	//calculate Q matrix
	for (int xi = 0; xi < NumberOfParticles; xi++) {
		for (int i = 0; i < NumberOfParticles; i++) {
			FlagCalGrid = false;
			OutputQ(i, xi) = 0.0;

			Eigen::Vector3i grid_i = FlatToGrid(i);
			Eigen::Vector3i grid_xi = FlatToGrid(xi);
			Eigen::Vector3i i_minus_xi = grid_i - grid_xi;

			// i - xi, j - xi, k - xi �� -1�`1�͈̔͂ł��邩�̔���
			//i - xi
			if ((abs(i_minus_xi[0]) <= 1) && (abs(i_minus_xi[1]) <= 1) && (abs(i_minus_xi[2]) <= 1)) {
				FlagCalGrid = true;
				//std::cout << "true\n";
			}

			if (FlagCalGrid) {
				// ��͓I
				// OutputQ(i, xi) = - pow(square.dx, 3) * Chloe[i_minus_xi[0] + 1][i_minus_xi[1] + 1][i_minus_xi[2] + 1];

				// ���l�I
				OutputQ(i, xi)
					= - Riemann_Sum_for_Chloe(i_minus_xi[0], square.dx)
					* Riemann_Sum_for_Chloe(i_minus_xi[1], square.dx)
					* Riemann_Sum_for_Chloe(i_minus_xi[2], square.dx);
				
				
			}

		}
	}



	return OutputQ;
}

Eigen::MatrixXd calMatrixR(Square square, Eigen::VectorXd phi, Eigen::VectorXd theta)
{
	bool FlagCalGrid1;
	bool FlagCalGrid2;
	Eigen::MatrixXd OutputR(NumberOfParticles, NumberOfParticles);
	Eigen::MatrixXd OutputR1(NumberOfParticles, NumberOfParticles);
	Eigen::MatrixXd OutputR2(NumberOfParticles, NumberOfParticles);

	//calculate R matrix
	for (int k = 0; k < NumberOfParticles; k++) {
		for (int i = 0; i < NumberOfParticles; i++) {

			FlagCalGrid1 = false;
			Eigen::Vector3i grid_i = FlatToGrid(i);
			Eigen::Vector3i grid_k = FlatToGrid(k);
			Eigen::Vector3i i_minus_k = grid_i - grid_k;
			OutputR1(i, k) = 0.0;
			OutputR2(i, k) = 0.0;
			
			// i - k, j - k �� -1�`1�͈̔͂ł��邩�̔���
				//i - k
			if ((abs(i_minus_k[0]) <= 1) && (abs(i_minus_k[1]) <= 1) && (abs(i_minus_k[2]) <= 1)) {
				FlagCalGrid1 = true;
				//std::cout << "true\n";
			}

			// ����Z�Ȃ�v�Z
			if (FlagCalGrid1) {
				OutputR1(i, k) = pow(square.dx, 3) * (kappa / 2) * Chloe[i_minus_k[0] + 1][i_minus_k[1] + 1][i_minus_k[2] + 1];
			}

			for (int j = 0; j < NumberOfParticles; j++) {

				FlagCalGrid2 = false;
				Eigen::Vector3i grid_j = FlatToGrid(j);
				Eigen::Vector3i j_minus_k = grid_j - grid_k;

				// i - k�� -1�`1�͈̔͂ł��邩�̔���@���@j - k��-1�`1�͈̔͂ł��邩�̔��� 
				//i - k
				if ((abs(i_minus_k[0]) <= 1) && (abs(i_minus_k[1]) <= 1) && (abs(i_minus_k[2]) <= 1)) {
					//j - k
					if ((abs(j_minus_k[0]) <= 1) && (abs(j_minus_k[1]) <= 1) && (abs(j_minus_k[2]) <= 1)) {
						FlagCalGrid2 = true;
						//std::cout << "true\n";
					}
				}
				

				// ����Z�Ȃ�v�Z
				if (FlagCalGrid2) {
					
					/*
					// ��͓I
					OutputR2(i, k) += pow(square.dx, 3) * (1 / pow(theta[i], 2) * (kappa / 2))
						* Aria[i_minus_k[0] + 1][i_minus_k[1] + 1][i_minus_k[2] + 1][j_minus_k[0] + 1][j_minus_k[1] + 1][j_minus_k[2] + 1];
					*/
				}
				
			}

			if ((abs(i_minus_k[0]) <= 1) && (abs(i_minus_k[1]) <= 1) && (abs(i_minus_k[2]) <= 1)) {
				// ���l�I
				OutputR2(i, k) += (kappa / 2)
					* Riemann_Sum_for_Luna(i_minus_k[0], k, theta, square.dx)
					* Riemann_Sum_for_Luna(i_minus_k[1], k, theta, square.dx)
					* Riemann_Sum_for_Luna(i_minus_k[2], k, theta, square.dx);
			}
			

			/*
			for (int a = 0; a < NumberOfParticles; a++) {
				for (int b = 0; b < NumberOfParticles; b++) {
					std::cout << OutputR2(a, b);
				}
				std::cout << std::endl;
			}
			*/
			

			OutputR(i, k) = OutputR1(i, k) + OutputR2(i, k);

		}
	}



	return OutputR;
}

Eigen::VectorXd calVectorb(Square square, Eigen::VectorXd phi, Eigen::VectorXd theta)
{
	bool FlagCalGrid;
	Eigen::VectorXd OutputV(NumberOfParticles);
	Eigen::VectorXd Vectorb1(NumberOfParticles);
	Eigen::VectorXd Vectorb2(NumberOfParticles);

	double bar_phi_i_x;
	double bar_phi_i_y;
	double bar_phi_i_z;

	double bar_phi_j_x;
	double bar_phi_j_y;
	double bar_phi_j_z;

	double bar_phi_k_x;
	double bar_phi_k_y;
	double bar_phi_k_z;

	//calculate b vector
	for (int xi = 0; xi < NumberOfParticles; xi++) {

		Eigen::Vector3i grid_xi = FlatToGrid(xi);

		Vectorb1(xi) = 0.0;
		Vectorb2(xi) = 0.0;

		for (int i = 0; i < NumberOfParticles; i++) {

			// �x�N�g��b1�̌v�Z���� ------------------------------------------------------------------------------
			FlagCalGrid = false;
			Eigen::Vector3i grid_i = FlatToGrid(i);

			Eigen::Vector3i i_minus_xi = grid_i - grid_xi;


			// i - xi�� -1�`1�͈̔͂ł��邩�̔���
			//i - xi
			if ((abs(i_minus_xi[0]) <= 1) && (abs(i_minus_xi[1]) <= 1) && (abs(i_minus_xi[2]) <= 1)) {
				FlagCalGrid = true;
			}

			if (FlagCalGrid == true) {
				double ww = Chloe[i_minus_xi[0] + 1][i_minus_xi[1] + 1][i_minus_xi[2] + 1];
				Vectorb1(xi) += theta(i) * ww;
			}

			// �x�N�g��b2�̌v�Z���� ------------------------------------------------------------------------------
			
			bar_phi_i_x = phi(3 * i);
			bar_phi_i_y = phi(3 * i + 1);
			bar_phi_i_z = phi(3 * i + 2);

			for (int j = 0; j < NumberOfParticles; j++) {

				Eigen::Vector3i grid_j = FlatToGrid(j);
				Eigen::Vector3i j_minus_xi = grid_j - grid_xi;

				bar_phi_j_x = phi(3 * j);
				bar_phi_j_y = phi(3 * j + 1);
				bar_phi_j_z = phi(3 * j + 2);

				for (int k = 0; k < NumberOfParticles; k++) {

					Eigen::Vector3i grid_k = FlatToGrid(k);
					Eigen::Vector3i k_minus_xi = grid_k - grid_xi;

					bar_phi_k_x = phi(3 * k);
					bar_phi_k_y = phi(3 * k + 1);
					bar_phi_k_z = phi(3 * k + 2);

					FlagCalGrid = false;

					// i - xi, j - xi, k - xi �� -1�`1�͈̔͂ł��邩�̔���
					//i - xi
					if ((abs(i_minus_xi[0]) <= 1) && (abs(i_minus_xi[1]) <= 1) && (abs(i_minus_xi[2]) <= 1)) {
						//j - xi
						if ((abs(j_minus_xi[0]) <= 1) && (abs(j_minus_xi[1]) <= 1) && (abs(j_minus_xi[2]) <= 1)) {
							//k - xi
							if ((abs(k_minus_xi[0]) <= 1) && (abs(k_minus_xi[1]) <= 1) && (abs(k_minus_xi[2]) <= 1)) {
								FlagCalGrid = true;
							}
						}
					}

					if (FlagCalGrid == true) {

						double Phi1 =
							bar_phi_i_x * (bar_phi_j_y * bar_phi_k_z - bar_phi_j_z * bar_phi_k_y)
							+ bar_phi_i_y * (bar_phi_j_z * bar_phi_k_x - bar_phi_j_x * bar_phi_k_z)
							+ bar_phi_i_z * (bar_phi_j_x * bar_phi_k_y - bar_phi_j_y * bar_phi_k_x);

						/*
						// ��͓I
						double sophia = SophiaX[0][1][2][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* SophiaY[0][1][2][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* SophiaZ[0][1][2][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						Vectorb2(xi) += Phi1 * sophia;
						*/

						// ���l�I
						vector<Eigen::Vector3i> v = { i_minus_xi, j_minus_xi, k_minus_xi };
						double w1w2w3 = Riemann_Sum_for_Sophia(v, { 0, 1, 2 }, square.dx);

						Vectorb2(xi) += Phi1 * w1w2w3;

						//std::cout << "No.  " << xi << "  b2 : " << "Phi = " << Phi << ", w1w2w3 = " << w1w2w3 << ", Phi1 * w1w2w3 = " << Phi1 * w1w2w3 << "\n";
						
					}

				}
			}


		}

		

		OutputV(xi) = (pow(square.dx, 3) * Vectorb1(xi)) - Vectorb2(xi);

		// std::cout << "No.  " << xi << "  b1 : " << pow(square.dx, 3) * Vectorb1(xi) << "\n";
		// std::cout << "No.  " << xi << "  b2 : " << Vectorb2(xi) << "\n";
		// std::cout << "No.  " << xi << "  V : " << OutputV(xi) << "\n";

	}

	// std::cout << Vectorb1 << "\n";

	return OutputV;
}

Eigen::VectorXd calVectorb_Convenient(Square square, Eigen::VectorXd phi, Eigen::VectorXd re_phi, Eigen::VectorXd theta)
{
	bool FlagCalGrid;
	Eigen::VectorXd OutputV(NumberOfParticles);
	Eigen::VectorXd Vectorb1(NumberOfParticles);
	Eigen::VectorXd Vectorb2(NumberOfParticles);

	double bar_phi_i_x = 0.0;
	double bar_phi_i_y = 0.0;
	double bar_phi_i_z = 0.0;

	double bar_phi_j_x = 0.0;
	double bar_phi_j_y = 0.0;
	double bar_phi_j_z = 0.0;

	double bar_phi_k_x = 0.0;
	double bar_phi_k_y = 0.0;
	double bar_phi_k_z = 0.0;

	//calculate b vector
	for (int xi = 0; xi < NumberOfParticles; xi++) {

		Eigen::Vector3i grid_xi = FlatToGrid(xi);

		Vectorb1(xi) = 0.0;
		Vectorb2(xi) = 0.0;

		// �x�N�g��b1�̌v�Z
		for (int i = 0; i < NumberOfParticles; i++) {

			Eigen::Vector3i grid_i = FlatToGrid(i);
			Eigen::Vector3i i_minus_xi = grid_i - grid_xi;

			double re_bar_phi_i_x = re_phi(3 * i);
			double re_bar_phi_i_y = re_phi(3 * i + 1);
			double re_bar_phi_i_z = re_phi(3 * i + 2);

			for (int j = 0; j < NumberOfParticles; j++) {

				Eigen::Vector3i grid_j = FlatToGrid(j);
				Eigen::Vector3i j_minus_xi = grid_j - grid_xi;

				double re_bar_phi_j_x = re_phi(3 * j);
				double re_bar_phi_j_y = re_phi(3 * j + 1);
				double re_bar_phi_j_z = re_phi(3 * j + 2);

				for (int k = 0; k < NumberOfParticles; k++) {

					Eigen::Vector3i grid_k = FlatToGrid(k);
					Eigen::Vector3i k_minus_xi = grid_k - grid_xi;
					
					double re_bar_phi_k_x = re_phi(3 * k);
					double re_bar_phi_k_y = re_phi(3 * k + 1);
					double re_bar_phi_k_z = re_phi(3 * k + 2);

					FlagCalGrid = false;

					// i - xi, j - xi, k - xi �� -1�`1�͈̔͂ł��邩�̔���
						//i - xi
					if ((abs(i_minus_xi[0]) <= 1) && (abs(i_minus_xi[1]) <= 1) && (abs(i_minus_xi[2]) <= 1)) {
						//j - xi
						if ((abs(j_minus_xi[0]) <= 1) && (abs(j_minus_xi[1]) <= 1) && (abs(j_minus_xi[2]) <= 1)) {
							//k - xi
							if ((abs(k_minus_xi[0]) <= 1) && (abs(k_minus_xi[1]) <= 1) && (abs(k_minus_xi[2]) <= 1)) {

								FlagCalGrid = true;
								//std::cout << "true\n";

							}
						}
					}

					if (FlagCalGrid == true) {

						double Phi1_b1 =
							re_bar_phi_i_x * (re_bar_phi_j_y * re_bar_phi_k_z - re_bar_phi_j_z * re_bar_phi_k_y)
							+ re_bar_phi_i_y * (re_bar_phi_j_z * re_bar_phi_k_x - re_bar_phi_j_x * re_bar_phi_k_z)
							+ re_bar_phi_i_z * (re_bar_phi_j_x * re_bar_phi_k_y - re_bar_phi_j_y * re_bar_phi_k_x);


						double Sophia_Convenient = SophiaX[0][1][2][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* SophiaY[0][1][2][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* SophiaZ[0][1][2][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];


						Vectorb1(xi) += Phi1_b1 * Sophia_Convenient;

					}

				}
			}
		}

		for (int i = 0; i < NumberOfParticles; i++) {

			Eigen::Vector3i grid_i = FlatToGrid(i);
			Eigen::Vector3i i_minus_xi = grid_i - grid_xi;

			double bar_phi_i_x = phi(3 * i);
			double bar_phi_i_y = phi(3 * i + 1);
			double bar_phi_i_z = phi(3 * i + 2);

			for (int j = 0; j < NumberOfParticles; j++) {

				Eigen::Vector3i grid_j = FlatToGrid(j);
				Eigen::Vector3i j_minus_xi = grid_j - grid_xi;

				double bar_phi_j_x = phi(3 * j);
				double bar_phi_j_y = phi(3 * j + 1);
				double bar_phi_j_z = phi(3 * j + 2);

				for (int k = 0; k < NumberOfParticles; k++) {

					Eigen::Vector3i grid_k = FlatToGrid(k);
					Eigen::Vector3i k_minus_xi = grid_k - grid_xi;

					double bar_phi_k_x = phi(3 * k);
					double bar_phi_k_y = phi(3 * k + 1);
					double bar_phi_k_z = phi(3 * k + 2);

					FlagCalGrid = false;

					// i - xi, j - xi, k - xi �� -1�`1�͈̔͂ł��邩�̔���
					//i - xi
					if ((abs(i_minus_xi[0]) <= 1) && (abs(i_minus_xi[1]) <= 1) && (abs(i_minus_xi[2]) <= 1)) {
						//j - xi
						if ((abs(j_minus_xi[0]) <= 1) && (abs(j_minus_xi[1]) <= 1) && (abs(j_minus_xi[2]) <= 1)) {
							//k - xi
							if ((abs(k_minus_xi[0]) <= 1) && (abs(k_minus_xi[1]) <= 1) && (abs(k_minus_xi[2]) <= 1)) {

								FlagCalGrid = true;
								//std::cout << "true\n";

							}
						}
					}

					if (FlagCalGrid == true) {

						double Phi1 =
							bar_phi_i_x * (bar_phi_j_y * bar_phi_k_z - bar_phi_j_z * bar_phi_k_y)
							+ bar_phi_i_y * (bar_phi_j_z * bar_phi_k_x - bar_phi_j_x * bar_phi_k_z)
							+ bar_phi_i_z * (bar_phi_j_x * bar_phi_k_y - bar_phi_j_y * bar_phi_k_x);

						double Sophia = SophiaX[0][1][2][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* SophiaY[0][1][2][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* SophiaZ[0][1][2][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						Vectorb2(xi) += Phi1 * Sophia;

					}

				}
			}
		}

		OutputV(xi) = Vectorb1(xi) - Vectorb2(xi);
	}

	// std::cout << Vectorb1 << std::endl;

	return OutputV;
}

Eigen::VectorXd calVectorc(Square square, Eigen::VectorXd phi, Eigen::VectorXd theta)
{
	bool FlagCalGrid;
	Eigen::VectorXd Outputc(NumberOfParticles);// square.size()
	Eigen::VectorXd Outputc1(NumberOfParticles);// square.size()
	Eigen::VectorXd Outputc2(NumberOfParticles);// square.size()

	
	for (int k = 0; k < NumberOfParticles; k++) {
		Eigen::Vector3i grid_k = FlatToGrid(k);

		// calculate c1 vector
		Outputc1(k) = 0.0;

		for (int i = 0; i < NumberOfParticles; i++) {
			FlagCalGrid = false;
			Eigen::Vector3i grid_i = FlatToGrid(i);
			Eigen::Vector3i i_minus_k = grid_i - grid_k;
			
			// i - k�� -1�`1�͈̔͂ł��邩�̔���
			//i - k
			if ((abs(i_minus_k[0]) <= 1) && (abs(i_minus_k[1]) <= 1) && (abs(i_minus_k[2]) <= 1)) {
				FlagCalGrid = true;
				//std::cout << "true\n";
			}

			// ����Z�Ȃ�v�Z
			if (FlagCalGrid) {
				Outputc1(k) += pow(square.dx, 3) * ( square.points[i].power - kappa / 2 * theta[i] ) * Chloe[i_minus_k[0] + 1][i_minus_k[1] + 1][i_minus_k[2] + 1];
			}

		}

		// calculate c2 vector
		Outputc2(k) = 0.0;

		for (int i = 0; i < NumberOfParticles; i++) {
			FlagCalGrid = false;
			Eigen::Vector3i grid_i = FlatToGrid(i);
			Eigen::Vector3i i_minus_k = grid_i - grid_k;

			// i - k�� -1�`1�͈̔͂ł��邩�̔���
			//i - k
			if ((abs(i_minus_k[0]) <= 1) && (abs(i_minus_k[1]) <= 1) && (abs(i_minus_k[2]) <= 1)) {
				FlagCalGrid = true;
				//std::cout << "true\n";
			}

			// ����Z�Ȃ�v�Z
			if (FlagCalGrid) {

				/*
				Outputc2(k) += kappa / 2 * (1 / theta[i])
					* Riemann_Sum_for_Mary(i_minus_k[0], square.dx)
					* Riemann_Sum_for_Mary(i_minus_k[1], square.dx) 
					* Riemann_Sum_for_Mary(i_minus_k[2], square.dx);
				*/

				Outputc2(k) += pow(square.dx, 3) * kappa / 2 * (1 / theta[i]) * Chloe[i_minus_k[0] + 1][i_minus_k[1] + 1][i_minus_k[2] + 1];

			}

		}

		Outputc(k) = Outputc1(k) + Outputc2(k);
		// std::cout << "No.  " << k << "  V1 : " << Outputc1(k) << "\n";
		// std::cout << "No.  " << k << "  V2 : " << Outputc2(k) << "\n";
	}

	return Outputc;
}

Eigen::MatrixXd calMatrixT(Eigen::MatrixXd M) {
	return (M) * (M.transpose());
}

Eigen::VectorXd calVectore(Eigen::VectorXd V, Eigen::MatrixXd M) {
	// �{���Cvector b �͓]�u���Ȃ����C��x�N�g���ł��邽�߁C�s�x�N�g���ɂ����D
	return V.transpose() * (M.transpose());
}

int GridToFlat(Eigen::Vector3i grid_index)
{
	int flat_index;
	int grid_x = grid_index[0];
	int grid_y = grid_index[1];
	int grid_z = grid_index[2];

	flat_index = grid_x * int(pow(NumberOfOneDemensionParticles, 2)) + grid_y * NumberOfOneDemensionParticles + grid_z;

	return flat_index;
};

Eigen::Vector3i FlatToGrid(int flat_index)
{
	Eigen::Vector3i grid_index = {};

	grid_index[0] = flat_index / int(pow(NumberOfOneDemensionParticles, 2));
	grid_index[1] = (flat_index % int(pow(NumberOfOneDemensionParticles, 2))) / NumberOfOneDemensionParticles;
	grid_index[2] = ((flat_index % int(pow(NumberOfOneDemensionParticles, 2))) % NumberOfOneDemensionParticles);

	return grid_index;
};

double HatFunction(double x)
{
	double x_num = (abs(x) < 1.0) ? (1.0 - abs(x)) : 0.0;

	return x_num;
}

double DifferentialHatFunction(double x)
{
	if (x >= -1.0 && x < 0.0) {
		return 1.0;
	}
	else if (x <= 1.0 && x > 0.0) {
		return -1.0;
	}
	else {
		return 0.0;
	}
}

double Riemann_Sum_for_Mary(int i, double h) {

	int d = 10; // �ЂƔ͈͂̕�����
	double w = 1.0 / d; // ������
	// �e���̑�������

	vector<double> cal_point;

	for (int a = 0; a < d; a++) {
		cal_point.emplace_back(1.0 / (2.0 * d) + a * w);
	}

	double totalSum = 0.0;

	for (int offset = -2; offset <= 1; offset++) {
		double sum = 0.0;
		for (int a = 0; a < d; a++) {
			double cal = static_cast<double>(offset) + cal_point[a];
			double denominator = HatFunction(cal - i);
			

			// �[�����Z������邽�߂̃`�F�b�N
			if (std::abs(denominator) < 1e-10) { // 1e-10 �͏�����臒l�ł�
				continue; // �P�ɉ��Z���X�L�b�v���܂�
			}
			else {
				
				sum += (1.0 / d) * h * (1 / denominator) * HatFunction(cal);
			}
		}
		// �eoffset�̍��v��S�̂̍��v�ɒǉ�
		totalSum += sum;
	}

	cal_point.clear();

	// printf("%lf\n", totalSum);

	return totalSum;
}

double Riemann_Sum_for_Luna(int i_minus_k, int k, Eigen::VectorXd theta, double h) {

	int d = 10; // �ЂƔ͈͂̕�����
	double w = 1.0 / d; // ������
	// �e���̑�������

	vector<double> cal_point;

	for (int a = 0; a < d; a++) {
		cal_point.emplace_back(1.0 / (2.0 * d) + a * w);
	}

	double totalSum = 0.0;

	for (int offset = -2; offset <= 1; offset++) {
		double sum = 0.0;
		for (int a = 0; a < d; a++) {
			double cal = static_cast<double>(offset) + cal_point[a];
			double sum_j = 0.0;
			for (int j = 0; j < NumberOfParticles; j++) {
				int j_minus_k = j - k;

				if(abs(j_minus_k) <= 1) {
					sum_j += theta[j] * HatFunction(cal - j_minus_k); // theta�̔z��̐ς��l��
				}
				
			}

			// �[�����Z������邽�߂̃`�F�b�N
			if (abs(sum_j) < 1e-10) { // 1e-10 �͏�����臒l�ł�
				continue; // �P�ɉ��Z���X�L�b�v���܂�
			}
			else {
				sum += w * h * (1 / pow(sum_j, 2)) * HatFunction(cal - i_minus_k) * HatFunction(cal);
			}
		}
	}

	cal_point.clear();

	return totalSum;
}

double Riemann_Sum_for_Chloe(int i, double h) {

	int d = 10; // �ЂƔ͈͂̕�����
	double w = 1.0 / d; // ������
	// �e���̑�������

	vector<double> cal_point;

	for (int a = 0; a < d; a++) {
		cal_point.emplace_back(1.0 / (2.0 * d) + a * w);
	}

	double totalSum = 0.0;

	for (int offset = -2; offset <= 1; offset++) {
		double sum = 0.0;
		for (int a = 0; a < d; a++) {
			double cal = static_cast<double>(offset) + cal_point[a];
			double denominator = HatFunction(cal - i);

			// �[�����Z������邽�߂̃`�F�b�N
			if (std::abs(denominator) < 1e-10) { // 1e-10 �͏�����臒l�ł�
				continue; // �P�ɉ��Z���X�L�b�v���܂�
			}
			else {
				sum += w * h * denominator * HatFunction(cal);
			}
		}

		// �eoffset�̍��v��S�̂̍��v�ɒǉ�
		totalSum += sum;

	}

	cal_point.clear();

	return totalSum;
}

double Riemann_Sum_Example(double h) {

	int d = 10; // �ЂƔ͈͂̕�����
	double w = 1.0 / d; // ������
	// �e���̑�������

	vector<double> cal_point;

	for (int a = 0; a < d; a++) {
		cal_point.emplace_back(1.0 / (2.0 * d) + a * w);
	}

	double totalSum = 0.0;

	for (int offset = -2; offset <= 1; offset++) {
		double sum = 0.0;
		for (int a = 0; a < d; a++) {
			double cal = static_cast<double>(offset) + cal_point[a];
			double denominator = HatFunction(cal);
			// cout << denominator << endl;
			// �[�����Z������邽�߂̃`�F�b�N
			if (std::abs(denominator) < 1e-10) { // 1e-10 �͏�����臒l�ł�
				continue; // �P�ɉ��Z���X�L�b�v���܂�
			}
			else {
				sum += w * h * denominator;
			}
		}

		// �eoffset�̍��v��S�̂̍��v�ɒǉ�
		totalSum += sum;

	}

	cal_point.clear();

	return totalSum;
}

double Riemann_Sum_for_Sophia(vector<Eigen::Vector3i> v, vector<int> axis, double h) {
	// axis = {1,2,0}�Ȃǂ�����

	int d = 10; // �ЂƔ͈͂̕�����
	double w = 1.0 / d; // ������
	// �e���̑�������

	vector<double> cal_point;

	for (int a = 0; a < d; a++) {
		cal_point.emplace_back(1.0 / (2.0 * d) + a * w);
	}

	double AllDimSum = 1.0;

	// �e�������v�Z
	for (int dim = 0; dim < 3; dim++) {

		// �e��Ԃ��ƂɌv�Z
		double totalSum = 0.0;
		for (int offset = -2; offset <= 1; offset++) {

			double sum = 0.0;
			for (int a = 0; a < d; a++) {
				double denominator = 1.0;
				double cal = static_cast<double>(offset) + cal_point[a];
				if (dim == axis[0]) {
					denominator *= (1 / h) * DifferentialHatFunction(cal - v[0][dim]) * HatFunction(cal - v[1][dim]) * HatFunction(cal - v[2][dim]);
				}
				else if (dim == axis[1]) {
					denominator *= (1 / h) * HatFunction(cal - v[0][dim]) * DifferentialHatFunction(cal - v[1][dim]) * HatFunction(cal - v[2][dim]);
				}
				else if (dim == axis[2]) {
					denominator *= (1 / h) * HatFunction(cal - v[0][dim]) * HatFunction(cal - v[1][dim]) * DifferentialHatFunction(cal - v[2][dim]);
				}
				else {
					denominator = 0.0;
				}

				// �[�����Z������邽�߂̃`�F�b�N
				if (std::abs(denominator) < 1e-10) { // 1e-10 �͏�����臒l�ł�
					continue; // �P�ɉ��Z���X�L�b�v���܂�
				}
				else {
					sum += w * h * denominator * HatFunction(cal);
				}
			}

			// �eoffset�̍��v��S�̂̍��v�ɒǉ�
			totalSum += sum;

		}
		AllDimSum *= totalSum;
	}

	
	cal_point.clear();

	return AllDimSum;
}