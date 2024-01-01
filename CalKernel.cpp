#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include "CalKernel.h"

// カーネル計算に必要な各区間の方程式の係数を格納する行列配列を作成する構造体を生成する関数
Kernel_Array createKernelArray() {

	// grid node i
	Eigen::MatrixXd c_ix = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd c_iy = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd c_iz = Eigen::MatrixXd::Zero(4, 2);

	// grid node j
	Eigen::MatrixXd c_jx = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd c_jy = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd c_jz = Eigen::MatrixXd::Zero(4, 2);

	// grid node k
	Eigen::MatrixXd c_kx = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd c_ky = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd c_kz = Eigen::MatrixXd::Zero(4, 2);

	// grid node l
	Eigen::MatrixXd c_lx = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd c_ly = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd c_lz = Eigen::MatrixXd::Zero(4, 2);

	// original point grid
	Eigen::MatrixXd c_o = Eigen::MatrixXd::Zero(4, 2);

	// 原点固定のkernel
	// [1+x][1-x]
	c_o(1, 0) = 1;
	c_o(1, 1) = 1;
	c_o(2, 0) = 1;
	c_o(2, 1) = -1;

	/*
	* 計算したいカーネルは「Psi*Psi*Psi*w(Sophia)」である
	*/
	Kernel_Array ka = Kernel_Array(
		c_ix, c_iy, c_iz,
		c_jx, c_jy, c_jz,
		c_kx, c_ky, c_kz,
		c_lx, c_ly, c_lz, c_o);

	return ka;
};


/*
* xyzのインデックスやグリッドノードによって，行列に適切な係数をセットする関数 ----------------
*/

Kernel_Array setCoefficientSophiaX(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_ix, int gn_jx, int gn_kx) {

	// 値の初期化
	ka.coefficient_i_x = Eigen::MatrixXd::Zero(4, 2);
	ka.coefficient_j_x = Eigen::MatrixXd::Zero(4, 2);
	ka.coefficient_k_x = Eigen::MatrixXd::Zero(4, 2);

	// grid node i について
	if (gn_ix == 1) {
		if (frag_i == 1)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.coefficient_i_x(2, 0) = 1;
			ka.coefficient_i_x(3, 0) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.coefficient_i_x(2, 0) = 1 - gn_ix;
			ka.coefficient_i_x(2, 1) = 1;
			ka.coefficient_i_x(3, 0) = 1 + gn_ix;
			ka.coefficient_i_x(3, 1) = -1;

		}
		else if (frag_i == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.coefficient_i_x(2, 0) = 1 - gn_ix;
			ka.coefficient_i_x(2, 1) = 1;
			ka.coefficient_i_x(3, 0) = 1 + gn_ix;
			ka.coefficient_i_x(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_ix == 0) {
		if (frag_i == 1)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.coefficient_i_x(1, 0) = 1;
			ka.coefficient_i_x(2, 0) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.coefficient_i_x(1, 0) = 1 - gn_ix;
			ka.coefficient_i_x(1, 1) = 1;
			ka.coefficient_i_x(2, 0) = 1 + gn_ix;
			ka.coefficient_i_x(2, 1) = -1;
		}
		else if (frag_i == 3)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.coefficient_i_x(1, 0) = 1 - gn_ix;
			ka.coefficient_i_x(1, 1) = 1;
			ka.coefficient_i_x(2, 0) = 1 + gn_ix;
			ka.coefficient_i_x(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_ix == -1)
	{
		if (frag_i == 1)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.coefficient_i_x(0, 0) = 1;
			ka.coefficient_i_x(1, 0) = -1;
		}
		else if (frag_i == 2)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.coefficient_i_x(0, 0) = 1 - gn_ix;
			ka.coefficient_i_x(0, 1) = 1;
			ka.coefficient_i_x(1, 0) = 1 + gn_ix;
			ka.coefficient_i_x(1, 1) = -1;
		}
		else if (frag_i == 3)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.coefficient_i_x(0, 0) = 1 - gn_ix;
			ka.coefficient_i_x(0, 1) = 1;
			ka.coefficient_i_x(1, 0) = 1 + gn_ix;
			ka.coefficient_i_x(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node ix value reference.\n";
		std::exit(1);
	}


	// grid node j について
	if (gn_jx == 1) {
		if (frag_j == 1)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.coefficient_j_x(2, 0) = 1;
			ka.coefficient_j_x(3, 0) = -1;

		}
		else if (frag_j == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.coefficient_j_x(2, 0) = 1 - gn_jx;
			ka.coefficient_j_x(2, 1) = 1;
			ka.coefficient_j_x(3, 0) = 1 + gn_jx;
			ka.coefficient_j_x(3, 1) = -1;

		}
		else if (frag_j == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.coefficient_j_x(2, 0) = 1 - gn_jx;
			ka.coefficient_j_x(2, 1) = 1;
			ka.coefficient_j_x(3, 0) = 1 + gn_jx;
			ka.coefficient_j_x(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jx == 0) {
		if (frag_j == 1)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.coefficient_j_x(1, 0) = 1;
			ka.coefficient_j_x(2, 0) = -1;

		}
		else if (frag_j == 2)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.coefficient_j_x(1, 0) = 1 - gn_jx;
			ka.coefficient_j_x(1, 1) = 1;
			ka.coefficient_j_x(2, 0) = 1 + gn_jx;
			ka.coefficient_j_x(2, 1) = -1;
		}
		else if (frag_j == 3)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.coefficient_j_x(1, 0) = 1 - gn_jx;
			ka.coefficient_j_x(1, 1) = 1;
			ka.coefficient_j_x(2, 0) = 1 + gn_jx;
			ka.coefficient_j_x(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jx == -1)
	{
		if (frag_j == 1)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.coefficient_j_x(0, 0) = 1;
			ka.coefficient_j_x(1, 0) = -1;
		}
		else if (frag_j == 2)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.coefficient_j_x(0, 0) = 1 - gn_jx;
			ka.coefficient_j_x(0, 1) = 1;
			ka.coefficient_j_x(1, 0) = 1 + gn_jx;
			ka.coefficient_j_x(1, 1) = -1;
		}
		else if (frag_j == 3)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.coefficient_j_x(0, 0) = 1 - gn_jx;
			ka.coefficient_j_x(0, 1) = 1;
			ka.coefficient_j_x(1, 0) = 1 + gn_jx;
			ka.coefficient_j_x(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node jx value reference.\n";
		std::exit(1);
	}

	// grid node k について
	if (gn_kx == 1) {
		if (frag_k == 1)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.coefficient_k_x(2, 0) = 1;
			ka.coefficient_k_x(3, 0) = -1;

		}
		else if (frag_k == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.coefficient_k_x(2, 0) = 1 - gn_kx;
			ka.coefficient_k_x(2, 1) = 1;
			ka.coefficient_k_x(3, 0) = 1 + gn_kx;
			ka.coefficient_k_x(3, 1) = -1;

		}
		else if (frag_k == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.coefficient_k_x(2, 0) = 1 - gn_kx;
			ka.coefficient_k_x(2, 1) = 1;
			ka.coefficient_k_x(3, 0) = 1 + gn_kx;
			ka.coefficient_k_x(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_kx == 0) {
		if (frag_k == 1)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.coefficient_k_x(1, 0) = 1;
			ka.coefficient_k_x(2, 0) = -1;

		}
		else if (frag_k == 2)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.coefficient_k_x(1, 0) = 1 - gn_kx;
			ka.coefficient_k_x(1, 1) = 1;
			ka.coefficient_k_x(2, 0) = 1 + gn_kx;
			ka.coefficient_k_x(2, 1) = -1;
		}
		else if (frag_k == 3)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.coefficient_k_x(1, 0) = 1 - gn_kx;
			ka.coefficient_k_x(1, 1) = 1;
			ka.coefficient_k_x(2, 0) = 1 + gn_kx;
			ka.coefficient_k_x(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_kx == -1)
	{
		if (frag_k == 1)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.coefficient_k_x(0, 0) = 1;
			ka.coefficient_k_x(1, 0) = -1;
		}
		else if (frag_k == 2)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.coefficient_k_x(0, 0) = 1 - gn_kx;
			ka.coefficient_k_x(0, 1) = 1;
			ka.coefficient_k_x(1, 0) = 1 + gn_kx;
			ka.coefficient_k_x(1, 1) = -1;
		}
		else if (frag_k == 3)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.coefficient_k_x(0, 0) = 1 - gn_kx;
			ka.coefficient_k_x(0, 1) = 1;
			ka.coefficient_k_x(1, 0) = 1 + gn_kx;
			ka.coefficient_k_x(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node kx value reference.\n";
		std::exit(1);
	}

	return ka;
}

Kernel_Array setCoefficientSophiaY(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_iy, int gn_jy, int gn_ky) {

	// 値の初期化
	ka.coefficient_i_y = Eigen::MatrixXd::Zero(4, 2);
	ka.coefficient_j_y = Eigen::MatrixXd::Zero(4, 2);
	ka.coefficient_k_y = Eigen::MatrixXd::Zero(4, 2);

	// grid node i について
	if (gn_iy == 1) {
		if (frag_i == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.coefficient_i_y(2, 0) = 1 - gn_iy;
			ka.coefficient_i_y(2, 1) = 1;
			ka.coefficient_i_y(3, 0) = 1 + gn_iy;
			ka.coefficient_i_y(3, 1) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.coefficient_i_y(2, 0) = 1;
			ka.coefficient_i_y(3, 0) = -1;

		}
		else if (frag_i == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.coefficient_i_y(2, 0) = 1 - gn_iy;
			ka.coefficient_i_y(2, 1) = 1;
			ka.coefficient_i_y(3, 0) = 1 + gn_iy;
			ka.coefficient_i_y(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_iy == 0) {
		if (frag_i == 1)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.coefficient_i_y(1, 0) = 1 - gn_iy;
			ka.coefficient_i_y(1, 1) = 1;
			ka.coefficient_i_y(2, 0) = 1 + gn_iy;
			ka.coefficient_i_y(2, 1) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.coefficient_i_y(1, 0) = 1;
			ka.coefficient_i_y(2, 0) = -1;
		}
		else if (frag_i == 3)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.coefficient_i_y(1, 0) = 1 - gn_iy;
			ka.coefficient_i_y(1, 1) = 1;
			ka.coefficient_i_y(2, 0) = 1 + gn_iy;
			ka.coefficient_i_y(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_iy == -1)
	{
		if (frag_i == 1)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.coefficient_i_y(0, 0) = 1 - gn_iy;
			ka.coefficient_i_y(0, 1) = 1;
			ka.coefficient_i_y(1, 0) = 1 + gn_iy;
			ka.coefficient_i_y(1, 1) = -1;
		}
		else if (frag_i == 2)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.coefficient_i_y(0, 0) = 1;
			ka.coefficient_i_y(1, 0) = -1;
		}
		else if (frag_i == 3)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.coefficient_i_y(0, 0) = 1 - gn_iy;
			ka.coefficient_i_y(0, 1) = 1;
			ka.coefficient_i_y(1, 0) = 1 + gn_iy;
			ka.coefficient_i_y(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node iy value reference.\n";
		std::exit(1);
	}


	// grid node j について
	if (gn_jy == 1) {
		if (frag_j == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.coefficient_j_y(2, 0) = 1 - gn_jy;
			ka.coefficient_j_y(2, 1) = 1;
			ka.coefficient_j_y(3, 0) = 1 + gn_jy;
			ka.coefficient_j_y(3, 1) = -1;

		}
		else if (frag_j == 2)
		{

			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.coefficient_j_y(2, 0) = 1;
			ka.coefficient_j_y(3, 0) = -1;

		}
		else if (frag_j == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.coefficient_j_y(2, 0) = 1 - gn_jy;
			ka.coefficient_j_y(2, 1) = 1;
			ka.coefficient_j_y(3, 0) = 1 + gn_jy;
			ka.coefficient_j_y(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jy == 0) {
		if (frag_j == 1)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.coefficient_j_y(1, 0) = 1 - gn_jy;
			ka.coefficient_j_y(1, 1) = 1;
			ka.coefficient_j_y(2, 0) = 1 + gn_jy;
			ka.coefficient_j_y(2, 1) = -1;
		}
		else if (frag_j == 2)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.coefficient_j_y(1, 0) = 1;
			ka.coefficient_j_y(2, 0) = -1;
		}
		else if (frag_j == 3)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.coefficient_j_y(1, 0) = 1 - gn_jy;
			ka.coefficient_j_y(1, 1) = 1;
			ka.coefficient_j_y(2, 0) = 1 + gn_jy;
			ka.coefficient_j_y(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jy == -1)
	{
		if (frag_j == 1)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.coefficient_j_y(0, 0) = 1 - gn_jy;
			ka.coefficient_j_y(0, 1) = 1;
			ka.coefficient_j_y(1, 0) = 1 + gn_jy;
			ka.coefficient_j_y(1, 1) = -1;
		}
		else if (frag_j == 2)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.coefficient_j_y(0, 0) = 1;
			ka.coefficient_j_y(1, 0) = -1;
		}
		else if (frag_j == 3)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.coefficient_j_y(0, 0) = 1 - gn_jy;
			ka.coefficient_j_y(0, 1) = 1;
			ka.coefficient_j_y(1, 0) = 1 + gn_jy;
			ka.coefficient_j_y(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node jy value reference.\n";
		std::exit(1);
	}

	// grid node k について
	if (gn_ky == 1) {
		if (frag_k == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.coefficient_k_y(2, 0) = 1 - gn_ky;
			ka.coefficient_k_y(2, 1) = 1;
			ka.coefficient_k_y(3, 0) = 1 + gn_ky;
			ka.coefficient_k_y(3, 1) = -1;
		}
		else if (frag_k == 2)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.coefficient_k_y(2, 0) = 1;
			ka.coefficient_k_y(3, 0) = -1;
		}
		else if (frag_k == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.coefficient_k_y(2, 0) = 1 - gn_ky;
			ka.coefficient_k_y(2, 1) = 1;
			ka.coefficient_k_y(3, 0) = 1 + gn_ky;
			ka.coefficient_k_y(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_ky == 0) {
		if (frag_k == 1)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.coefficient_k_y(1, 0) = 1 - gn_ky;
			ka.coefficient_k_y(1, 1) = 1;
			ka.coefficient_k_y(2, 0) = 1 + gn_ky;
			ka.coefficient_k_y(2, 1) = -1;

		}
		else if (frag_k == 2)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.coefficient_k_y(1, 0) = 1;
			ka.coefficient_k_y(2, 0) = -1;
		}
		else if (frag_k == 3)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.coefficient_k_y(1, 0) = 1 - gn_ky;
			ka.coefficient_k_y(1, 1) = 1;
			ka.coefficient_k_y(2, 0) = 1 + gn_ky;
			ka.coefficient_k_y(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_ky == -1)
	{
		if (frag_k == 1)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.coefficient_k_y(0, 0) = 1 - gn_ky;
			ka.coefficient_k_y(0, 1) = 1;
			ka.coefficient_k_y(1, 0) = 1 + gn_ky;
			ka.coefficient_k_y(1, 1) = -1;
		}
		else if (frag_k == 2)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.coefficient_k_y(0, 0) = 1;
			ka.coefficient_k_y(1, 0) = -1;
		}
		else if (frag_k == 3)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.coefficient_k_y(0, 0) = 1 - gn_ky;
			ka.coefficient_k_y(0, 1) = 1;
			ka.coefficient_k_y(1, 0) = 1 + gn_ky;
			ka.coefficient_k_y(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node ky value reference.\n";
		std::exit(1);
	}

	return ka;
}

Kernel_Array setCoefficientSophiaZ(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_iz, int gn_jz, int gn_kz) {

	// 値の初期化
	ka.coefficient_i_z = Eigen::MatrixXd::Zero(4, 2);
	ka.coefficient_j_z = Eigen::MatrixXd::Zero(4, 2);
	ka.coefficient_k_z = Eigen::MatrixXd::Zero(4, 2);

	// grid node i について
	if (gn_iz == 1) {
		if (frag_i == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.coefficient_i_z(2, 0) = 1 - gn_iz;
			ka.coefficient_i_z(2, 1) = 1;
			ka.coefficient_i_z(3, 0) = 1 + gn_iz;
			ka.coefficient_i_z(3, 1) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.coefficient_i_z(2, 0) = 1 - gn_iz;
			ka.coefficient_i_z(2, 1) = 1;
			ka.coefficient_i_z(3, 0) = 1 + gn_iz;
			ka.coefficient_i_z(3, 1) = -1;
		}
		else if (frag_i == 3)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.coefficient_i_z(2, 0) = 1;
			ka.coefficient_i_z(3, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_iz == 0) {
		if (frag_i == 1)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.coefficient_i_z(1, 0) = 1 - gn_iz;
			ka.coefficient_i_z(1, 1) = 1;
			ka.coefficient_i_z(2, 0) = 1 + gn_iz;
			ka.coefficient_i_z(2, 1) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.coefficient_i_z(1, 0) = 1 - gn_iz;
			ka.coefficient_i_z(1, 1) = 1;
			ka.coefficient_i_z(2, 0) = 1 + gn_iz;
			ka.coefficient_i_z(2, 1) = -1;
		}
		else if (frag_i == 3)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.coefficient_i_z(1, 0) = 1;
			ka.coefficient_i_z(2, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_iz == -1)
	{
		if (frag_i == 1)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.coefficient_i_z(0, 0) = 1 - gn_iz;
			ka.coefficient_i_z(0, 1) = 1;
			ka.coefficient_i_z(1, 0) = 1 + gn_iz;
			ka.coefficient_i_z(1, 1) = -1;
		}
		else if (frag_i == 2)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.coefficient_i_z(0, 0) = 1 - gn_iz;
			ka.coefficient_i_z(0, 1) = 1;
			ka.coefficient_i_z(1, 0) = 1 + gn_iz;
			ka.coefficient_i_z(1, 1) = -1;
		}
		else if (frag_i == 3)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.coefficient_i_z(0, 0) = 1;
			ka.coefficient_i_z(1, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node i value reference.\n";
		std::exit(1);
	}


	// grid node j について
	if (gn_jz == 1) {
		if (frag_j == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.coefficient_j_z(2, 0) = 1 - gn_jz;
			ka.coefficient_j_z(2, 1) = 1;
			ka.coefficient_j_z(3, 0) = 1 + gn_jz;
			ka.coefficient_j_z(3, 1) = -1;

		}
		else if (frag_j == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.coefficient_j_z(2, 0) = 1 - gn_jz;
			ka.coefficient_j_z(2, 1) = 1;
			ka.coefficient_j_z(3, 0) = 1 + gn_jz;
			ka.coefficient_j_z(3, 1) = -1;
		}
		else if (frag_j == 3)
		{

			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.coefficient_j_z(2, 0) = 1;
			ka.coefficient_j_z(3, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jz == 0) {
		if (frag_j == 1)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.coefficient_j_z(1, 0) = 1 - gn_jz;
			ka.coefficient_j_z(1, 1) = 1;
			ka.coefficient_j_z(2, 0) = 1 + gn_jz;
			ka.coefficient_j_z(2, 1) = -1;
		}
		else if (frag_j == 2)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.coefficient_j_z(1, 0) = 1 - gn_jz;
			ka.coefficient_j_z(1, 1) = 1;
			ka.coefficient_j_z(2, 0) = 1 + gn_jz;
			ka.coefficient_j_z(2, 1) = -1;
		}
		else if (frag_j == 3)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.coefficient_j_z(1, 0) = 1;
			ka.coefficient_j_z(2, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jz == -1)
	{
		if (frag_j == 1)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.coefficient_j_z(0, 0) = 1 - gn_jz;
			ka.coefficient_j_z(0, 1) = 1;
			ka.coefficient_j_z(1, 0) = 1 + gn_jz;
			ka.coefficient_j_z(1, 1) = -1;
		}
		else if (frag_j == 2)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.coefficient_j_z(0, 0) = 1 - gn_jz;
			ka.coefficient_j_z(0, 1) = 1;
			ka.coefficient_j_z(1, 0) = 1 + gn_jz;
			ka.coefficient_j_z(1, 1) = -1;
		}
		else if (frag_j == 3)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.coefficient_j_z(0, 0) = 1;
			ka.coefficient_j_z(1, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node j value reference.\n";
		std::exit(1);
	}

	// grid node k について
	if (gn_kz == 1) {
		if (frag_k == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.coefficient_k_z(2, 0) = 1 - gn_kz;
			ka.coefficient_k_z(2, 1) = 1;
			ka.coefficient_k_z(3, 0) = 1 + gn_kz;
			ka.coefficient_k_z(3, 1) = -1;
		}
		else if (frag_k == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.coefficient_k_z(2, 0) = 1 - gn_kz;
			ka.coefficient_k_z(2, 1) = 1;
			ka.coefficient_k_z(3, 0) = 1 + gn_kz;
			ka.coefficient_k_z(3, 1) = -1;
		}
		else if (frag_k == 3)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.coefficient_k_z(2, 0) = 1;
			ka.coefficient_k_z(3, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_kz == 0) {
		if (frag_k == 1)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.coefficient_k_z(1, 0) = 1 - gn_kz;
			ka.coefficient_k_z(1, 1) = 1;
			ka.coefficient_k_z(2, 0) = 1 + gn_kz;
			ka.coefficient_k_z(2, 1) = -1;

		}
		else if (frag_k == 2)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.coefficient_k_z(1, 0) = 1 - gn_kz;
			ka.coefficient_k_z(1, 1) = 1;
			ka.coefficient_k_z(2, 0) = 1 + gn_kz;
			ka.coefficient_k_z(2, 1) = -1;
		}
		else if (frag_k == 3)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.coefficient_k_z(1, 0) = 1;
			ka.coefficient_k_z(2, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_kz == -1)
	{
		if (frag_k == 1)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.coefficient_k_z(0, 0) = 1 - gn_kz;
			ka.coefficient_k_z(0, 1) = 1;
			ka.coefficient_k_z(1, 0) = 1 + gn_kz;
			ka.coefficient_k_z(1, 1) = -1;
		}
		else if (frag_k == 2)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.coefficient_k_z(0, 0) = 1 - gn_kz;
			ka.coefficient_k_z(0, 1) = 1;
			ka.coefficient_k_z(1, 0) = 1 + gn_kz;
			ka.coefficient_k_z(1, 1) = -1;
		}
		else if (frag_k == 3)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.coefficient_k_z(0, 0) = 1;
			ka.coefficient_k_z(1, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node k value reference.\n";
		std::exit(1);
	}

	return ka;
}


Kernel_Array setCoefficientChloeX(Kernel_Array ka, int gn_ix) {

	// 値の初期化
	ka.coefficient_i_x = Eigen::MatrixXd::Zero(4, 2);

	// grid node ix について
	if (gn_ix == 1)
	{
		ka.coefficient_i_x(2, 0) = 1 - gn_ix;
		ka.coefficient_i_x(2, 1) = 1;
		ka.coefficient_i_x(3, 0) = 1 + gn_ix;
		ka.coefficient_i_x(3, 1) = -1;
	}
	else if (gn_ix == 0)
	{
		ka.coefficient_i_x(1, 0) = 1 - gn_ix;
		ka.coefficient_i_x(1, 1) = 1;
		ka.coefficient_i_x(2, 0) = 1 + gn_ix;
		ka.coefficient_i_x(2, 1) = -1;
	}
	else if (gn_ix == -1)
	{
		ka.coefficient_i_x(0, 0) = 1 - gn_ix;
		ka.coefficient_i_x(0, 1) = 1;
		ka.coefficient_i_x(1, 0) = 1 + gn_ix;
		ka.coefficient_i_x(1, 1) = -1;
	}
	else
	{
		std::cout << "Impossible grid node x value reference.\n";
		std::exit(1);
	}

	return ka;
}

Kernel_Array setCoefficientChloeY(Kernel_Array ka, int gn_iy) {

	// 値の初期化
	ka.coefficient_i_y = Eigen::MatrixXd::Zero(4, 2);

	// grid node ix について
	if (gn_iy == 1)
	{
		ka.coefficient_i_y(2, 0) = 1 - gn_iy;
		ka.coefficient_i_y(2, 1) = 1;
		ka.coefficient_i_y(3, 0) = 1 + gn_iy;
		ka.coefficient_i_y(3, 1) = -1;
	}
	else if (gn_iy == 0)
	{
		ka.coefficient_i_y(1, 0) = 1 - gn_iy;
		ka.coefficient_i_y(1, 1) = 1;
		ka.coefficient_i_y(2, 0) = 1 + gn_iy;
		ka.coefficient_i_y(2, 1) = -1;
	}
	else if (gn_iy == -1)
	{
		ka.coefficient_i_y(0, 0) = 1 - gn_iy;
		ka.coefficient_i_y(0, 1) = 1;
		ka.coefficient_i_y(1, 0) = 1 + gn_iy;
		ka.coefficient_i_y(1, 1) = -1;
	}
	else
	{
		std::cout << "Impossible grid node y value reference.\n";
		std::exit(1);
	}

	return ka;
}

Kernel_Array setCoefficientChloeZ(Kernel_Array ka, int gn_iz) {

	// 値の初期化
	ka.coefficient_i_z = Eigen::MatrixXd::Zero(4, 2);

	// grid node ix について
	if (gn_iz == 1)
	{
		ka.coefficient_i_z(2, 0) = 1 - gn_iz;
		ka.coefficient_i_z(2, 1) = 1;
		ka.coefficient_i_z(3, 0) = 1 + gn_iz;
		ka.coefficient_i_z(3, 1) = -1;
	}
	else if (gn_iz == 0)
	{
		ka.coefficient_i_z(1, 0) = 1 - gn_iz;
		ka.coefficient_i_z(1, 1) = 1;
		ka.coefficient_i_z(2, 0) = 1 + gn_iz;
		ka.coefficient_i_z(2, 1) = -1;
	}
	else if (gn_iz == -1)
	{
		ka.coefficient_i_z(0, 0) = 1 - gn_iz;
		ka.coefficient_i_z(0, 1) = 1;
		ka.coefficient_i_z(1, 0) = 1 + gn_iz;
		ka.coefficient_i_z(1, 1) = -1;
	}
	else
	{
		std::cout << "Impossible grid node z value reference.\n";
		std::exit(1);
	}

	return ka;
}


Kernel_Array setCoefficientAriaX(Kernel_Array ka, int gn_ix, int gn_jx) {

	// 値の初期化
	ka.coefficient_i_x = Eigen::MatrixXd::Zero(4, 2);
	ka.coefficient_j_x = Eigen::MatrixXd::Zero(4, 2);

	// grid node ix について
	if (gn_ix == 1)
	{
		ka.coefficient_i_x(2, 0) = 1 - gn_ix;
		ka.coefficient_i_x(2, 1) = 1;
		ka.coefficient_i_x(3, 0) = 1 + gn_ix;
		ka.coefficient_i_x(3, 1) = -1;
	}
	else if (gn_ix == 0)
	{
		ka.coefficient_i_x(1, 0) = 1 - gn_ix;
		ka.coefficient_i_x(1, 1) = 1;
		ka.coefficient_i_x(2, 0) = 1 + gn_ix;
		ka.coefficient_i_x(2, 1) = -1;
	}
	else if (gn_ix == -1)
	{
		ka.coefficient_i_x(0, 0) = 1 - gn_ix;
		ka.coefficient_i_x(0, 1) = 1;
		ka.coefficient_i_x(1, 0) = 1 + gn_ix;
		ka.coefficient_i_x(1, 1) = -1;
	}
	else
	{
		std::cout << "Impossible grid node x value reference.\n";
		std::exit(1);
	}

	// grid node jx について
	if (gn_jx == 1)
	{
		ka.coefficient_j_x(2, 0) = 1 - gn_jx;
		ka.coefficient_j_x(2, 1) = 1;
		ka.coefficient_j_x(3, 0) = 1 + gn_jx;
		ka.coefficient_j_x(3, 1) = -1;
	}
	else if (gn_jx == 0)
	{
		ka.coefficient_j_x(1, 0) = 1 - gn_jx;
		ka.coefficient_j_x(1, 1) = 1;
		ka.coefficient_j_x(2, 0) = 1 + gn_jx;
		ka.coefficient_j_x(2, 1) = -1;
	}
	else if (gn_jx == -1)
	{
		ka.coefficient_j_x(0, 0) = 1 - gn_jx;
		ka.coefficient_j_x(0, 1) = 1;
		ka.coefficient_j_x(1, 0) = 1 + gn_jx;
		ka.coefficient_j_x(1, 1) = -1;
	}
	else
	{
		std::cout << "Impossible grid node x value reference.\n";
		std::exit(1);
	}

	return ka;
}

Kernel_Array setCoefficientAriaY(Kernel_Array ka, int gn_iy, int gn_jy) {

	// 値の初期化
	ka.coefficient_i_y = Eigen::MatrixXd::Zero(4, 2);
	ka.coefficient_j_y = Eigen::MatrixXd::Zero(4, 2);

	// grid node iy について
	if (gn_iy == 1)
	{
		ka.coefficient_i_y(2, 0) = 1 - gn_iy;
		ka.coefficient_i_y(2, 1) = 1;
		ka.coefficient_i_y(3, 0) = 1 + gn_iy;
		ka.coefficient_i_y(3, 1) = -1;
	}
	else if (gn_iy == 0)
	{
		ka.coefficient_i_y(1, 0) = 1 - gn_iy;
		ka.coefficient_i_y(1, 1) = 1;
		ka.coefficient_i_y(2, 0) = 1 + gn_iy;
		ka.coefficient_i_y(2, 1) = -1;
	}
	else if (gn_iy == -1)
	{
		ka.coefficient_i_y(0, 0) = 1 - gn_iy;
		ka.coefficient_i_y(0, 1) = 1;
		ka.coefficient_i_y(1, 0) = 1 + gn_iy;
		ka.coefficient_i_y(1, 1) = -1;
	}
	else
	{
		std::cout << "Impossible grid node y value reference.\n";
		std::exit(1);
	}

	// grid node jy について
	if (gn_jy == 1)
	{
		ka.coefficient_j_y(2, 0) = 1 - gn_jy;
		ka.coefficient_j_y(2, 1) = 1;
		ka.coefficient_j_y(3, 0) = 1 + gn_jy;
		ka.coefficient_j_y(3, 1) = -1;
	}
	else if (gn_jy == 0)
	{
		ka.coefficient_j_y(1, 0) = 1 - gn_jy;
		ka.coefficient_j_y(1, 1) = 1;
		ka.coefficient_j_y(2, 0) = 1 + gn_jy;
		ka.coefficient_j_y(2, 1) = -1;
	}
	else if (gn_jy == -1)
	{
		ka.coefficient_j_y(0, 0) = 1 - gn_jy;
		ka.coefficient_j_y(0, 1) = 1;
		ka.coefficient_j_y(1, 0) = 1 + gn_jy;
		ka.coefficient_j_y(1, 1) = -1;
	}
	else
	{
		std::cout << "Impossible grid node y value reference.\n";
		std::exit(1);
	}

	return ka;
}

Kernel_Array setCoefficientAriaZ(Kernel_Array ka, int gn_iz, int gn_jz) {

	// 値の初期化
	ka.coefficient_i_z = Eigen::MatrixXd::Zero(4, 2);
	ka.coefficient_j_z = Eigen::MatrixXd::Zero(4, 2);

	// grid node iz について
	if (gn_iz == 1)
	{
		ka.coefficient_i_z(2, 0) = 1 - gn_iz;
		ka.coefficient_i_z(2, 1) = 1;
		ka.coefficient_i_z(3, 0) = 1 + gn_iz;
		ka.coefficient_i_z(3, 1) = -1;
	}
	else if (gn_iz == 0)
	{
		ka.coefficient_i_z(1, 0) = 1 - gn_iz;
		ka.coefficient_i_z(1, 1) = 1;
		ka.coefficient_i_z(2, 0) = 1 + gn_iz;
		ka.coefficient_i_z(2, 1) = -1;
	}
	else if (gn_iz == -1)
	{
		ka.coefficient_i_z(0, 0) = 1 - gn_iz;
		ka.coefficient_i_z(0, 1) = 1;
		ka.coefficient_i_z(1, 0) = 1 + gn_iz;
		ka.coefficient_i_z(1, 1) = -1;
	}
	else
	{
		std::cout << "Impossible grid node z value reference.\n";
		std::exit(1);
	}

	// grid node jz について
	if (gn_jz == 1)
	{
		ka.coefficient_j_z(2, 0) = 1 - gn_jz;
		ka.coefficient_j_z(2, 1) = 1;
		ka.coefficient_j_z(3, 0) = 1 + gn_jz;
		ka.coefficient_j_z(3, 1) = -1;
	}
	else if (gn_jz == 0)
	{
		ka.coefficient_j_z(1, 0) = 1 - gn_jz;
		ka.coefficient_j_z(1, 1) = 1;
		ka.coefficient_j_z(2, 0) = 1 + gn_jz;
		ka.coefficient_j_z(2, 1) = -1;
	}
	else if (gn_jz == -1)
	{
		ka.coefficient_j_z(0, 0) = 1 - gn_jz;
		ka.coefficient_j_z(0, 1) = 1;
		ka.coefficient_j_z(1, 0) = 1 + gn_jz;
		ka.coefficient_j_z(1, 1) = -1;
	}
	else
	{
		std::cout << "Impossible grid node z value reference.\n";
		std::exit(1);
	}

	return ka;
}


Kernel_Array setCoefficientMiaX(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_ix, int gn_jx, int gn_kx, int gn_lx) {

	// 値の初期化
	ka.coefficient_i_x = Eigen::MatrixXd::Zero(4, 2);
	ka.coefficient_j_x = Eigen::MatrixXd::Zero(4, 2);
	ka.coefficient_k_x = Eigen::MatrixXd::Zero(4, 2);
	ka.coefficient_l_x = Eigen::MatrixXd::Zero(4, 2);

	// grid node i について
	if (gn_ix == 1) {
		if (frag_i == 1)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.coefficient_i_x(2, 0) = 1;
			ka.coefficient_i_x(3, 0) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.coefficient_i_x(2, 0) = 1 - gn_ix;
			ka.coefficient_i_x(2, 1) = 1;
			ka.coefficient_i_x(3, 0) = 1 + gn_ix;
			ka.coefficient_i_x(3, 1) = -1;

		}
		else if (frag_i == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.coefficient_i_x(2, 0) = 1 - gn_ix;
			ka.coefficient_i_x(2, 1) = 1;
			ka.coefficient_i_x(3, 0) = 1 + gn_ix;
			ka.coefficient_i_x(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_ix == 0) {
		if (frag_i == 1)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.coefficient_i_x(1, 0) = 1;
			ka.coefficient_i_x(2, 0) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.coefficient_i_x(1, 0) = 1 - gn_ix;
			ka.coefficient_i_x(1, 1) = 1;
			ka.coefficient_i_x(2, 0) = 1 + gn_ix;
			ka.coefficient_i_x(2, 1) = -1;
		}
		else if (frag_i == 3)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.coefficient_i_x(1, 0) = 1 - gn_ix;
			ka.coefficient_i_x(1, 1) = 1;
			ka.coefficient_i_x(2, 0) = 1 + gn_ix;
			ka.coefficient_i_x(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_ix == -1)
	{
		if (frag_i == 1)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.coefficient_i_x(0, 0) = 1;
			ka.coefficient_i_x(1, 0) = -1;
		}
		else if (frag_i == 2)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.coefficient_i_x(0, 0) = 1 - gn_ix;
			ka.coefficient_i_x(0, 1) = 1;
			ka.coefficient_i_x(1, 0) = 1 + gn_ix;
			ka.coefficient_i_x(1, 1) = -1;
		}
		else if (frag_i == 3)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.coefficient_i_x(0, 0) = 1 - gn_ix;
			ka.coefficient_i_x(0, 1) = 1;
			ka.coefficient_i_x(1, 0) = 1 + gn_ix;
			ka.coefficient_i_x(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node ix value reference.\n";
		std::exit(1);
	}


	// grid node j について
	if (gn_jx == 1) {
		if (frag_j == 1)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.coefficient_j_x(2, 0) = 1;
			ka.coefficient_j_x(3, 0) = -1;

		}
		else if (frag_j == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.coefficient_j_x(2, 0) = 1 - gn_jx;
			ka.coefficient_j_x(2, 1) = 1;
			ka.coefficient_j_x(3, 0) = 1 + gn_jx;
			ka.coefficient_j_x(3, 1) = -1;

		}
		else if (frag_j == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.coefficient_j_x(2, 0) = 1 - gn_jx;
			ka.coefficient_j_x(2, 1) = 1;
			ka.coefficient_j_x(3, 0) = 1 + gn_jx;
			ka.coefficient_j_x(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jx == 0) {
		if (frag_j == 1)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.coefficient_j_x(1, 0) = 1;
			ka.coefficient_j_x(2, 0) = -1;

		}
		else if (frag_j == 2)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.coefficient_j_x(1, 0) = 1 - gn_jx;
			ka.coefficient_j_x(1, 1) = 1;
			ka.coefficient_j_x(2, 0) = 1 + gn_jx;
			ka.coefficient_j_x(2, 1) = -1;
		}
		else if (frag_j == 3)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.coefficient_j_x(1, 0) = 1 - gn_jx;
			ka.coefficient_j_x(1, 1) = 1;
			ka.coefficient_j_x(2, 0) = 1 + gn_jx;
			ka.coefficient_j_x(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jx == -1)
	{
		if (frag_j == 1)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.coefficient_j_x(0, 0) = 1;
			ka.coefficient_j_x(1, 0) = -1;
		}
		else if (frag_j == 2)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.coefficient_j_x(0, 0) = 1 - gn_jx;
			ka.coefficient_j_x(0, 1) = 1;
			ka.coefficient_j_x(1, 0) = 1 + gn_jx;
			ka.coefficient_j_x(1, 1) = -1;
		}
		else if (frag_j == 3)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.coefficient_j_x(0, 0) = 1 - gn_jx;
			ka.coefficient_j_x(0, 1) = 1;
			ka.coefficient_j_x(1, 0) = 1 + gn_jx;
			ka.coefficient_j_x(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node jx value reference.\n";
		std::exit(1);
	}

	// grid node k について
	if (gn_kx == 1) {
		if (frag_k == 1)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.coefficient_k_x(2, 0) = 1;
			ka.coefficient_k_x(3, 0) = -1;

		}
		else if (frag_k == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.coefficient_k_x(2, 0) = 1 - gn_kx;
			ka.coefficient_k_x(2, 1) = 1;
			ka.coefficient_k_x(3, 0) = 1 + gn_kx;
			ka.coefficient_k_x(3, 1) = -1;

		}
		else if (frag_k == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.coefficient_k_x(2, 0) = 1 - gn_kx;
			ka.coefficient_k_x(2, 1) = 1;
			ka.coefficient_k_x(3, 0) = 1 + gn_kx;
			ka.coefficient_k_x(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_kx == 0) {
		if (frag_k == 1)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.coefficient_k_x(1, 0) = 1;
			ka.coefficient_k_x(2, 0) = -1;

		}
		else if (frag_k == 2)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.coefficient_k_x(1, 0) = 1 - gn_kx;
			ka.coefficient_k_x(1, 1) = 1;
			ka.coefficient_k_x(2, 0) = 1 + gn_kx;
			ka.coefficient_k_x(2, 1) = -1;
		}
		else if (frag_k == 3)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.coefficient_k_x(1, 0) = 1 - gn_kx;
			ka.coefficient_k_x(1, 1) = 1;
			ka.coefficient_k_x(2, 0) = 1 + gn_kx;
			ka.coefficient_k_x(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_kx == -1)
	{
		if (frag_k == 1)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.coefficient_k_x(0, 0) = 1;
			ka.coefficient_k_x(1, 0) = -1;
		}
		else if (frag_k == 2)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.coefficient_k_x(0, 0) = 1 - gn_kx;
			ka.coefficient_k_x(0, 1) = 1;
			ka.coefficient_k_x(1, 0) = 1 + gn_kx;
			ka.coefficient_k_x(1, 1) = -1;
		}
		else if (frag_k == 3)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.coefficient_k_x(0, 0) = 1 - gn_kx;
			ka.coefficient_k_x(0, 1) = 1;
			ka.coefficient_k_x(1, 0) = 1 + gn_kx;
			ka.coefficient_k_x(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node kx value reference.\n";
		std::exit(1);
	}

	// grid node l について
	if (gn_lx == 1) {
		ka.coefficient_l_x(2, 0) = 1 - gn_lx;
		ka.coefficient_l_x(2, 1) = 1;
		ka.coefficient_l_x(3, 0) = 1 + gn_lx;
		ka.coefficient_l_x(3, 1) = -1;
	}
	else if (gn_lx == 0) {
		ka.coefficient_l_x(1, 0) = 1 - gn_lx;
		ka.coefficient_l_x(1, 1) = 1;
		ka.coefficient_l_x(2, 0) = 1 + gn_lx;
		ka.coefficient_l_x(2, 1) = -1;
	}
	else if (gn_lx == -1)
	{
		ka.coefficient_l_x(0, 0) = 1 - gn_lx;
		ka.coefficient_l_x(0, 1) = 1;
		ka.coefficient_l_x(1, 0) = 1 + gn_lx;
		ka.coefficient_l_x(1, 1) = -1;
	}
	else
	{
		std::cout << "Impossible grid node lx value reference.\n";
		std::exit(1);
	}

	/*
	if (frag_i == 1 && frag_j == 2 && frag_k == 3) {

		std::cout << "formula " << std::endl;
		for (int i = 0; i < 4; i++) {
			std::cout << "range : " << i << std::endl;
			std::cout << ka.coefficient_i_x(i, 0) << std::endl;
			std::cout << ka.coefficient_i_x(i, 1) << std::endl;
		}

	}
	*/

	return ka;
}

Kernel_Array setCoefficientMiaY(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_iy, int gn_jy, int gn_ky, int gn_ly) {

	// 値の初期化
	ka.coefficient_i_y = Eigen::MatrixXd::Zero(4, 2);
	ka.coefficient_j_y = Eigen::MatrixXd::Zero(4, 2);
	ka.coefficient_k_y = Eigen::MatrixXd::Zero(4, 2);
	ka.coefficient_l_y = Eigen::MatrixXd::Zero(4, 2);

	// grid node i について
	if (gn_iy == 1) {
		if (frag_i == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.coefficient_i_y(2, 0) = 1 - gn_iy;
			ka.coefficient_i_y(2, 1) = 1;
			ka.coefficient_i_y(3, 0) = 1 + gn_iy;
			ka.coefficient_i_y(3, 1) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.coefficient_i_y(2, 0) = 1;
			ka.coefficient_i_y(3, 0) = -1;

		}
		else if (frag_i == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.coefficient_i_y(2, 0) = 1 - gn_iy;
			ka.coefficient_i_y(2, 1) = 1;
			ka.coefficient_i_y(3, 0) = 1 + gn_iy;
			ka.coefficient_i_y(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_iy == 0) {
		if (frag_i == 1)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.coefficient_i_y(1, 0) = 1 - gn_iy;
			ka.coefficient_i_y(1, 1) = 1;
			ka.coefficient_i_y(2, 0) = 1 + gn_iy;
			ka.coefficient_i_y(2, 1) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.coefficient_i_y(1, 0) = 1;
			ka.coefficient_i_y(2, 0) = -1;
		}
		else if (frag_i == 3)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.coefficient_i_y(1, 0) = 1 - gn_iy;
			ka.coefficient_i_y(1, 1) = 1;
			ka.coefficient_i_y(2, 0) = 1 + gn_iy;
			ka.coefficient_i_y(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_iy == -1)
	{
		if (frag_i == 1)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.coefficient_i_y(0, 0) = 1 - gn_iy;
			ka.coefficient_i_y(0, 1) = 1;
			ka.coefficient_i_y(1, 0) = 1 + gn_iy;
			ka.coefficient_i_y(1, 1) = -1;
		}
		else if (frag_i == 2)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.coefficient_i_y(0, 0) = 1;
			ka.coefficient_i_y(1, 0) = -1;
		}
		else if (frag_i == 3)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.coefficient_i_y(0, 0) = 1 - gn_iy;
			ka.coefficient_i_y(0, 1) = 1;
			ka.coefficient_i_y(1, 0) = 1 + gn_iy;
			ka.coefficient_i_y(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node iy value reference.\n";
		std::exit(1);
	}


	// grid node j について
	if (gn_jy == 1) {
		if (frag_j == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.coefficient_j_y(2, 0) = 1 - gn_jy;
			ka.coefficient_j_y(2, 1) = 1;
			ka.coefficient_j_y(3, 0) = 1 + gn_jy;
			ka.coefficient_j_y(3, 1) = -1;

		}
		else if (frag_j == 2)
		{

			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.coefficient_j_y(2, 0) = 1;
			ka.coefficient_j_y(3, 0) = -1;

		}
		else if (frag_j == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.coefficient_j_y(2, 0) = 1 - gn_jy;
			ka.coefficient_j_y(2, 1) = 1;
			ka.coefficient_j_y(3, 0) = 1 + gn_jy;
			ka.coefficient_j_y(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jy == 0) {
		if (frag_j == 1)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.coefficient_j_y(1, 0) = 1 - gn_jy;
			ka.coefficient_j_y(1, 1) = 1;
			ka.coefficient_j_y(2, 0) = 1 + gn_jy;
			ka.coefficient_j_y(2, 1) = -1;
		}
		else if (frag_j == 2)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.coefficient_j_y(1, 0) = 1;
			ka.coefficient_j_y(2, 0) = -1;
		}
		else if (frag_j == 3)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.coefficient_j_y(1, 0) = 1 - gn_jy;
			ka.coefficient_j_y(1, 1) = 1;
			ka.coefficient_j_y(2, 0) = 1 + gn_jy;
			ka.coefficient_j_y(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jy == -1)
	{
		if (frag_j == 1)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.coefficient_j_y(0, 0) = 1 - gn_jy;
			ka.coefficient_j_y(0, 1) = 1;
			ka.coefficient_j_y(1, 0) = 1 + gn_jy;
			ka.coefficient_j_y(1, 1) = -1;
		}
		else if (frag_j == 2)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.coefficient_j_y(0, 0) = 1;
			ka.coefficient_j_y(1, 0) = -1;
		}
		else if (frag_j == 3)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.coefficient_j_y(0, 0) = 1 - gn_jy;
			ka.coefficient_j_y(0, 1) = 1;
			ka.coefficient_j_y(1, 0) = 1 + gn_jy;
			ka.coefficient_j_y(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node jy value reference.\n";
		std::exit(1);
	}

	// grid node k について
	if (gn_ky == 1) {
		if (frag_k == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.coefficient_k_y(2, 0) = 1 - gn_ky;
			ka.coefficient_k_y(2, 1) = 1;
			ka.coefficient_k_y(3, 0) = 1 + gn_ky;
			ka.coefficient_k_y(3, 1) = -1;
		}
		else if (frag_k == 2)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.coefficient_k_y(2, 0) = 1;
			ka.coefficient_k_y(3, 0) = -1;
		}
		else if (frag_k == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.coefficient_k_y(2, 0) = 1 - gn_ky;
			ka.coefficient_k_y(2, 1) = 1;
			ka.coefficient_k_y(3, 0) = 1 + gn_ky;
			ka.coefficient_k_y(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_ky == 0) {
		if (frag_k == 1)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.coefficient_k_y(1, 0) = 1 - gn_ky;
			ka.coefficient_k_y(1, 1) = 1;
			ka.coefficient_k_y(2, 0) = 1 + gn_ky;
			ka.coefficient_k_y(2, 1) = -1;

		}
		else if (frag_k == 2)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.coefficient_k_y(1, 0) = 1;
			ka.coefficient_k_y(2, 0) = -1;
		}
		else if (frag_k == 3)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.coefficient_k_y(1, 0) = 1 - gn_ky;
			ka.coefficient_k_y(1, 1) = 1;
			ka.coefficient_k_y(2, 0) = 1 + gn_ky;
			ka.coefficient_k_y(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_ky == -1)
	{
		if (frag_k == 1)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.coefficient_k_y(0, 0) = 1 - gn_ky;
			ka.coefficient_k_y(0, 1) = 1;
			ka.coefficient_k_y(1, 0) = 1 + gn_ky;
			ka.coefficient_k_y(1, 1) = -1;
		}
		else if (frag_k == 2)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.coefficient_k_y(0, 0) = 1;
			ka.coefficient_k_y(1, 0) = -1;
		}
		else if (frag_k == 3)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.coefficient_k_y(0, 0) = 1 - gn_ky;
			ka.coefficient_k_y(0, 1) = 1;
			ka.coefficient_k_y(1, 0) = 1 + gn_ky;
			ka.coefficient_k_y(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node ky value reference.\n";
		std::exit(1);
	}

	// grid node l について
	if (gn_ly == 1) {
		ka.coefficient_l_y(2, 0) = 1 - gn_ly;
		ka.coefficient_l_y(2, 1) = 1;
		ka.coefficient_l_y(3, 0) = 1 + gn_ly;
		ka.coefficient_l_y(3, 1) = -1;
	}
	else if (gn_ly == 0) {
		ka.coefficient_l_y(1, 0) = 1 - gn_ly;
		ka.coefficient_l_y(1, 1) = 1;
		ka.coefficient_l_y(2, 0) = 1 + gn_ly;
		ka.coefficient_l_y(2, 1) = -1;
	}
	else if (gn_ly == -1)
	{
		ka.coefficient_l_y(0, 0) = 1 - gn_ly;
		ka.coefficient_l_y(0, 1) = 1;
		ka.coefficient_l_y(1, 0) = 1 + gn_ly;
		ka.coefficient_l_y(1, 1) = -1;
	}
	else
	{
		std::cout << "Impossible grid node ly value reference.\n";
		std::exit(1);
	}

	return ka;
}

Kernel_Array setCoefficientMiaZ(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_iz, int gn_jz, int gn_kz, int gn_lz) {

	// 値の初期化
	ka.coefficient_i_z = Eigen::MatrixXd::Zero(4, 2);
	ka.coefficient_j_z = Eigen::MatrixXd::Zero(4, 2);
	ka.coefficient_k_z = Eigen::MatrixXd::Zero(4, 2);
	ka.coefficient_l_z = Eigen::MatrixXd::Zero(4, 2);

	// grid node i について
	if (gn_iz == 1) {
		if (frag_i == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.coefficient_i_z(2, 0) = 1 - gn_iz;
			ka.coefficient_i_z(2, 1) = 1;
			ka.coefficient_i_z(3, 0) = 1 + gn_iz;
			ka.coefficient_i_z(3, 1) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.coefficient_i_z(2, 0) = 1 - gn_iz;
			ka.coefficient_i_z(2, 1) = 1;
			ka.coefficient_i_z(3, 0) = 1 + gn_iz;
			ka.coefficient_i_z(3, 1) = -1;
		}
		else if (frag_i == 3)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.coefficient_i_z(2, 0) = 1;
			ka.coefficient_i_z(3, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_iz == 0) {
		if (frag_i == 1)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.coefficient_i_z(1, 0) = 1 - gn_iz;
			ka.coefficient_i_z(1, 1) = 1;
			ka.coefficient_i_z(2, 0) = 1 + gn_iz;
			ka.coefficient_i_z(2, 1) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.coefficient_i_z(1, 0) = 1 - gn_iz;
			ka.coefficient_i_z(1, 1) = 1;
			ka.coefficient_i_z(2, 0) = 1 + gn_iz;
			ka.coefficient_i_z(2, 1) = -1;
		}
		else if (frag_i == 3)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.coefficient_i_z(1, 0) = 1;
			ka.coefficient_i_z(2, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_iz == -1)
	{
		if (frag_i == 1)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.coefficient_i_z(0, 0) = 1 - gn_iz;
			ka.coefficient_i_z(0, 1) = 1;
			ka.coefficient_i_z(1, 0) = 1 + gn_iz;
			ka.coefficient_i_z(1, 1) = -1;
		}
		else if (frag_i == 2)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.coefficient_i_z(0, 0) = 1 - gn_iz;
			ka.coefficient_i_z(0, 1) = 1;
			ka.coefficient_i_z(1, 0) = 1 + gn_iz;
			ka.coefficient_i_z(1, 1) = -1;
		}
		else if (frag_i == 3)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.coefficient_i_z(0, 0) = 1;
			ka.coefficient_i_z(1, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node i value reference.\n";
		std::exit(1);
	}


	// grid node j について
	if (gn_jz == 1) {
		if (frag_j == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.coefficient_j_z(2, 0) = 1 - gn_jz;
			ka.coefficient_j_z(2, 1) = 1;
			ka.coefficient_j_z(3, 0) = 1 + gn_jz;
			ka.coefficient_j_z(3, 1) = -1;

		}
		else if (frag_j == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.coefficient_j_z(2, 0) = 1 - gn_jz;
			ka.coefficient_j_z(2, 1) = 1;
			ka.coefficient_j_z(3, 0) = 1 + gn_jz;
			ka.coefficient_j_z(3, 1) = -1;
		}
		else if (frag_j == 3)
		{

			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.coefficient_j_z(2, 0) = 1;
			ka.coefficient_j_z(3, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jz == 0) {
		if (frag_j == 1)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.coefficient_j_z(1, 0) = 1 - gn_jz;
			ka.coefficient_j_z(1, 1) = 1;
			ka.coefficient_j_z(2, 0) = 1 + gn_jz;
			ka.coefficient_j_z(2, 1) = -1;
		}
		else if (frag_j == 2)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.coefficient_j_z(1, 0) = 1 - gn_jz;
			ka.coefficient_j_z(1, 1) = 1;
			ka.coefficient_j_z(2, 0) = 1 + gn_jz;
			ka.coefficient_j_z(2, 1) = -1;
		}
		else if (frag_j == 3)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.coefficient_j_z(1, 0) = 1;
			ka.coefficient_j_z(2, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jz == -1)
	{
		if (frag_j == 1)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.coefficient_j_z(0, 0) = 1 - gn_jz;
			ka.coefficient_j_z(0, 1) = 1;
			ka.coefficient_j_z(1, 0) = 1 + gn_jz;
			ka.coefficient_j_z(1, 1) = -1;
		}
		else if (frag_j == 2)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.coefficient_j_z(0, 0) = 1 - gn_jz;
			ka.coefficient_j_z(0, 1) = 1;
			ka.coefficient_j_z(1, 0) = 1 + gn_jz;
			ka.coefficient_j_z(1, 1) = -1;
		}
		else if (frag_j == 3)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.coefficient_j_z(0, 0) = 1;
			ka.coefficient_j_z(1, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node j value reference.\n";
		std::exit(1);
	}

	// grid node k について
	if (gn_kz == 1) {
		if (frag_k == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.coefficient_k_z(2, 0) = 1 - gn_kz;
			ka.coefficient_k_z(2, 1) = 1;
			ka.coefficient_k_z(3, 0) = 1 + gn_kz;
			ka.coefficient_k_z(3, 1) = -1;
		}
		else if (frag_k == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.coefficient_k_z(2, 0) = 1 - gn_kz;
			ka.coefficient_k_z(2, 1) = 1;
			ka.coefficient_k_z(3, 0) = 1 + gn_kz;
			ka.coefficient_k_z(3, 1) = -1;
		}
		else if (frag_k == 3)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.coefficient_k_z(2, 0) = 1;
			ka.coefficient_k_z(3, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_kz == 0) {
		if (frag_k == 1)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.coefficient_k_z(1, 0) = 1 - gn_kz;
			ka.coefficient_k_z(1, 1) = 1;
			ka.coefficient_k_z(2, 0) = 1 + gn_kz;
			ka.coefficient_k_z(2, 1) = -1;

		}
		else if (frag_k == 2)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.coefficient_k_z(1, 0) = 1 - gn_kz;
			ka.coefficient_k_z(1, 1) = 1;
			ka.coefficient_k_z(2, 0) = 1 + gn_kz;
			ka.coefficient_k_z(2, 1) = -1;
		}
		else if (frag_k == 3)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.coefficient_k_z(1, 0) = 1;
			ka.coefficient_k_z(2, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_kz == -1)
	{
		if (frag_k == 1)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.coefficient_k_z(0, 0) = 1 - gn_kz;
			ka.coefficient_k_z(0, 1) = 1;
			ka.coefficient_k_z(1, 0) = 1 + gn_kz;
			ka.coefficient_k_z(1, 1) = -1;
		}
		else if (frag_k == 2)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.coefficient_k_z(0, 0) = 1 - gn_kz;
			ka.coefficient_k_z(0, 1) = 1;
			ka.coefficient_k_z(1, 0) = 1 + gn_kz;
			ka.coefficient_k_z(1, 1) = -1;
		}
		else if (frag_k == 3)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.coefficient_k_z(0, 0) = 1;
			ka.coefficient_k_z(1, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node k value reference.\n";
		std::exit(1);
	}

	// grid node l について
	if (gn_lz == 1) {
		ka.coefficient_l_z(2, 0) = 1 - gn_lz;
		ka.coefficient_l_z(2, 1) = 1;
		ka.coefficient_l_z(3, 0) = 1 + gn_lz;
		ka.coefficient_l_z(3, 1) = -1;
	}
	else if (gn_lz == 0) {
		ka.coefficient_l_z(1, 0) = 1 - gn_lz;
		ka.coefficient_l_z(1, 1) = 1;
		ka.coefficient_l_z(2, 0) = 1 + gn_lz;
		ka.coefficient_l_z(2, 1) = -1;
	}
	else if (gn_lz == -1)
	{
		ka.coefficient_l_z(0, 0) = 1 - gn_lz;
		ka.coefficient_l_z(0, 1) = 1;
		ka.coefficient_l_z(1, 0) = 1 + gn_lz;
		ka.coefficient_l_z(1, 1) = -1;
	}
	else
	{
		std::cout << "Impossible grid node lz value reference.\n";
		std::exit(1);
	}

	return ka;
}

// -------------------------------------------------------------------------------------------


// 定積分計算
double integral(int start_p, int stop_p, std::vector<double> array)
{
	double ans = 0.0;
	for (int i = 0; i < array.size(); i++) {
		ans += array[i] * pow(stop_p, i + 1) / (i + 1) - array[i] * pow(start_p, i + 1) / (i + 1);
	}

	return ans;
}

// 変数式の積
std::vector<double> product(std::vector<double> a1, std::vector<double> a2)
{
	int length1 = a1.size();
	int length2 = a2.size();

	int lem = length1 + length2 - 1;

	std::vector<double> ans_array;

	for (int i = 0; i < lem; i++)
	{
		ans_array.emplace_back(0.0);
	}

	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {

			ans_array[i + j] += a1[i] * a2[j];
		}
	}

	return ans_array;
}


/*
* 各内装関数の積 -------------------------------------------
*/

double productSophiaX(Kernel_Array ka) {
	double ans_x = 0.0;
	std::vector<double> formula;
	std::vector<double> a_i_vector;
	std::vector<double> a_j_vector;
	std::vector<double> a_k_vector;
	std::vector<double> a_o_vector;

	// std::cout << "X ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_i_vector.clear();
		a_j_vector.clear();
		a_k_vector.clear();
		a_o_vector.clear();

		a_i_vector.emplace_back(ka.coefficient_i_x(a, 0));
		a_i_vector.emplace_back(ka.coefficient_i_x(a, 1));

		a_j_vector.emplace_back(ka.coefficient_j_x(a, 0));
		a_j_vector.emplace_back(ka.coefficient_j_x(a, 1));

		a_k_vector.emplace_back(ka.coefficient_k_x(a, 0));
		a_k_vector.emplace_back(ka.coefficient_k_x(a, 1));

		a_o_vector.emplace_back(ka.coefficient_o(a, 0));
		a_o_vector.emplace_back(ka.coefficient_o(a, 1));

		formula = product(product(a_i_vector, a_j_vector), product(a_k_vector, a_o_vector));
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "x dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/* std::cout << "formula : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/

		ans_x += ans;
		///std::abs(ans);
	}

	return ans_x;
}

double productSophiaY(Kernel_Array ka) {
	double ans_y = 0.0;
	std::vector<double> formula;
	std::vector<double> a_i_vector;
	std::vector<double> a_j_vector;
	std::vector<double> a_k_vector;
	std::vector<double> a_o_vector;

	// std::cout << "Y ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_i_vector.clear();
		a_j_vector.clear();
		a_k_vector.clear();
		a_o_vector.clear();

		a_i_vector.emplace_back(ka.coefficient_i_y(a, 0));
		a_i_vector.emplace_back(ka.coefficient_i_y(a, 1));

		a_j_vector.emplace_back(ka.coefficient_j_y(a, 0));
		a_j_vector.emplace_back(ka.coefficient_j_y(a, 1));

		a_k_vector.emplace_back(ka.coefficient_k_y(a, 0));
		a_k_vector.emplace_back(ka.coefficient_k_y(a, 1));

		a_o_vector.emplace_back(ka.coefficient_o(a, 0));
		a_o_vector.emplace_back(ka.coefficient_o(a, 1));

		formula = product(product(a_i_vector, a_j_vector), product(a_k_vector, a_o_vector));
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "y dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/*
		std::cout << "formula : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/

		ans_y += ans;
		//std::abs(ans);
	}

	return ans_y;
}

double productSophiaZ(Kernel_Array ka) {
	double ans_z = 0.0;
	std::vector<double> formula;
	std::vector<double> a_i_vector;
	std::vector<double> a_j_vector;
	std::vector<double> a_k_vector;
	std::vector<double> a_o_vector;

	// std::cout << "Z ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_i_vector.clear();
		a_j_vector.clear();
		a_k_vector.clear();
		a_o_vector.clear();

		a_i_vector.emplace_back(ka.coefficient_i_z(a, 0));
		a_i_vector.emplace_back(ka.coefficient_i_z(a, 1));

		a_j_vector.emplace_back(ka.coefficient_j_z(a, 0));
		a_j_vector.emplace_back(ka.coefficient_j_z(a, 1));

		a_k_vector.emplace_back(ka.coefficient_k_z(a, 0));
		a_k_vector.emplace_back(ka.coefficient_k_z(a, 1));

		a_o_vector.emplace_back(ka.coefficient_o(a, 0));
		a_o_vector.emplace_back(ka.coefficient_o(a, 1));

		formula = product(product(a_i_vector, a_j_vector), product(a_k_vector, a_o_vector));
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "z dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/*
		std::cout << "formula : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/

		ans_z += ans;
		//std::abs(ans);
	}

	return ans_z;
}


double productChloeX(Kernel_Array ka) {
	double ans_x = 0.0;
	std::vector<double> formula;
	std::vector<double> a_i_vector;
	std::vector<double> a_o_vector;

	// std::cout << "X ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_i_vector.clear();
		a_o_vector.clear();

		a_i_vector.emplace_back(ka.coefficient_i_x(a, 0));
		a_i_vector.emplace_back(ka.coefficient_i_x(a, 1));
		a_o_vector.emplace_back(ka.coefficient_o(a, 0));
		a_o_vector.emplace_back(ka.coefficient_o(a, 1));

		formula = product(a_i_vector, a_o_vector);
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "x dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/*
		std::cout << "formula : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/

		ans_x += ans;
		//std::abs(ans);
	}

	return ans_x;
}

double productChloeY(Kernel_Array ka) {
	double ans_y = 0.0;
	std::vector<double> formula;
	std::vector<double> a_i_vector;
	std::vector<double> a_o_vector;

	// std::cout << "Y ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_i_vector.clear();
		a_o_vector.clear();

		a_i_vector.emplace_back(ka.coefficient_i_y(a, 0));
		a_i_vector.emplace_back(ka.coefficient_i_y(a, 1));
		a_o_vector.emplace_back(ka.coefficient_o(a, 0));
		a_o_vector.emplace_back(ka.coefficient_o(a, 1));

		formula = product(a_i_vector, a_o_vector);
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "y dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/*
		std::cout << "formula : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/

		ans_y += ans;
		//std::abs(ans);
	}

	return ans_y;
}

double productChloeZ(Kernel_Array ka) {
	double ans_z = 0.0;
	std::vector<double> formula;
	std::vector<double> a_i_vector;
	std::vector<double> a_o_vector;

	// std::cout << "Z ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_i_vector.clear();
		a_o_vector.clear();

		a_i_vector.emplace_back(ka.coefficient_i_z(a, 0));
		a_i_vector.emplace_back(ka.coefficient_i_z(a, 1));
		a_o_vector.emplace_back(ka.coefficient_o(a, 0));
		a_o_vector.emplace_back(ka.coefficient_o(a, 1));

		formula = product(a_i_vector, a_o_vector);
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "z dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/* std::cout << "formula : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/

		ans_z += ans;
		//std::abs(ans);
	}

	return ans_z;
}


double productAriaX(Kernel_Array ka) {
	double ans_x = 0.0;
	std::vector<double> formula;
	std::vector<double> a_i_vector;
	std::vector<double> a_j_vector;
	std::vector<double> a_o_vector;

	// std::cout << "X ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_i_vector.clear();
		a_j_vector.clear();
		a_o_vector.clear();

		a_i_vector.emplace_back(ka.coefficient_i_x(a, 0));
		a_i_vector.emplace_back(ka.coefficient_i_x(a, 1));
		a_j_vector.emplace_back(ka.coefficient_j_x(a, 0));
		a_j_vector.emplace_back(ka.coefficient_j_x(a, 1));
		a_o_vector.emplace_back(ka.coefficient_o(a, 0));
		a_o_vector.emplace_back(ka.coefficient_o(a, 1));

		formula = product(product(a_i_vector, a_j_vector), a_o_vector);
		
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "x dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/*
		std::cout << "formula X : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/

		ans_x += ans;
		//std::abs(ans);
	}

	// std::cout << "\n";

	return ans_x;
}

double productAriaY(Kernel_Array ka) {
	double ans_y = 0.0;
	std::vector<double> formula;
	std::vector<double> a_i_vector;
	std::vector<double> a_j_vector;
	std::vector<double> a_o_vector;

	// std::cout << "Y ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_i_vector.clear();
		a_j_vector.clear();
		a_o_vector.clear();

		a_i_vector.emplace_back(ka.coefficient_i_y(a, 0));
		a_i_vector.emplace_back(ka.coefficient_i_y(a, 1));
		a_j_vector.emplace_back(ka.coefficient_j_y(a, 0));
		a_j_vector.emplace_back(ka.coefficient_j_y(a, 1));
		a_o_vector.emplace_back(ka.coefficient_o(a, 0));
		a_o_vector.emplace_back(ka.coefficient_o(a, 1));

		formula = product(product(a_i_vector, a_j_vector), a_o_vector);
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "y dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/*
		std::cout << "formula Y : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/

		ans_y += ans;
		//std::abs(ans);
	}

	return ans_y;
}

double productAriaZ(Kernel_Array ka) {
	double ans_z = 0.0;
	std::vector<double> formula;
	std::vector<double> a_i_vector;
	std::vector<double> a_j_vector;
	std::vector<double> a_o_vector;

	// std::cout << "Z ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_i_vector.clear();
		a_j_vector.clear();
		a_o_vector.clear();

		a_i_vector.emplace_back(ka.coefficient_i_z(a, 0));
		a_i_vector.emplace_back(ka.coefficient_i_z(a, 1));
		a_j_vector.emplace_back(ka.coefficient_j_z(a, 0));
		a_j_vector.emplace_back(ka.coefficient_j_z(a, 1));
		a_o_vector.emplace_back(ka.coefficient_o(a, 0));
		a_o_vector.emplace_back(ka.coefficient_o(a, 1));

		formula = product(product(a_i_vector, a_j_vector), a_o_vector);
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "z dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/* std::cout << "formula Z : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/

		ans_z += ans;
		//std::abs(ans);
	}

	return ans_z;
}


double productMiaX(Kernel_Array ka) {
	double ans_x = 0.0;
	std::vector<double> formula;
	std::vector<double> a_i_vector;
	std::vector<double> a_j_vector;
	std::vector<double> a_k_vector;
	std::vector<double> a_l_vector;
	std::vector<double> a_o_vector;

	// std::cout << "X ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_i_vector.clear();
		a_j_vector.clear();
		a_k_vector.clear();
		a_l_vector.clear();
		a_o_vector.clear();

		a_i_vector.emplace_back(ka.coefficient_i_x(a, 0));
		a_i_vector.emplace_back(ka.coefficient_i_x(a, 1));

		a_j_vector.emplace_back(ka.coefficient_j_x(a, 0));
		a_j_vector.emplace_back(ka.coefficient_j_x(a, 1));

		a_k_vector.emplace_back(ka.coefficient_k_x(a, 0));
		a_k_vector.emplace_back(ka.coefficient_k_x(a, 1));

		a_l_vector.emplace_back(ka.coefficient_l_x(a, 0));
		a_l_vector.emplace_back(ka.coefficient_l_x(a, 1));

		a_o_vector.emplace_back(ka.coefficient_o(a, 0));
		a_o_vector.emplace_back(ka.coefficient_o(a, 1));


		formula = product(product(product(a_i_vector, a_j_vector), product(a_k_vector, a_l_vector)), a_o_vector);
		
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "x dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		/*
		// output formula
		std::cout << "formula ";
		std::cout << a << " : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/
		
		

		ans_x += ans;
		///std::abs(ans);
	}

	return ans_x;
}

double productMiaY(Kernel_Array ka) {
	double ans_y = 0.0;
	std::vector<double> formula;
	std::vector<double> a_i_vector;
	std::vector<double> a_j_vector;
	std::vector<double> a_k_vector;
	std::vector<double> a_l_vector;
	std::vector<double> a_o_vector;

	// std::cout << "Y ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_i_vector.clear();
		a_j_vector.clear();
		a_k_vector.clear();
		a_l_vector.clear();
		a_o_vector.clear();

		a_i_vector.emplace_back(ka.coefficient_i_y(a, 0));
		a_i_vector.emplace_back(ka.coefficient_i_y(a, 1));

		a_j_vector.emplace_back(ka.coefficient_j_y(a, 0));
		a_j_vector.emplace_back(ka.coefficient_j_y(a, 1));

		a_k_vector.emplace_back(ka.coefficient_k_y(a, 0));
		a_k_vector.emplace_back(ka.coefficient_k_y(a, 1));

		a_l_vector.emplace_back(ka.coefficient_l_y(a, 0));
		a_l_vector.emplace_back(ka.coefficient_l_y(a, 1));

		a_o_vector.emplace_back(ka.coefficient_o(a, 0));
		a_o_vector.emplace_back(ka.coefficient_o(a, 1));

		formula = product(product(product(a_i_vector, a_j_vector), product(a_k_vector, a_l_vector)), a_o_vector);
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "y dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/*
		std::cout << "formula : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/

		ans_y += ans;
		//std::abs(ans);
	}

	return ans_y;
}

double productMiaZ(Kernel_Array ka) {
	double ans_z = 0.0;
	std::vector<double> formula;
	std::vector<double> a_i_vector;
	std::vector<double> a_j_vector;
	std::vector<double> a_k_vector;
	std::vector<double> a_l_vector;
	std::vector<double> a_o_vector;

	// std::cout << "Z ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_i_vector.clear();
		a_j_vector.clear();
		a_k_vector.clear();
		a_l_vector.clear();
		a_o_vector.clear();

		a_i_vector.emplace_back(ka.coefficient_i_z(a, 0));
		a_i_vector.emplace_back(ka.coefficient_i_z(a, 1));

		a_j_vector.emplace_back(ka.coefficient_j_z(a, 0));
		a_j_vector.emplace_back(ka.coefficient_j_z(a, 1));

		a_k_vector.emplace_back(ka.coefficient_k_z(a, 0));
		a_k_vector.emplace_back(ka.coefficient_k_z(a, 1));

		a_l_vector.emplace_back(ka.coefficient_l_z(a, 0));
		a_l_vector.emplace_back(ka.coefficient_l_z(a, 1));

		a_o_vector.emplace_back(ka.coefficient_o(a, 0));
		a_o_vector.emplace_back(ka.coefficient_o(a, 1));

		formula = product(product(product(a_i_vector, a_j_vector), product(a_k_vector, a_l_vector)), a_o_vector);
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "z dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/*
		std::cout << "formula : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/

		ans_z += ans;
		//std::abs(ans);
	}

	return ans_z;
}

// ---------------------------------------------------------