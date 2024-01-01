#ifndef __CALKERNEL_H__
#define __CALKERNEL_H__

#include <Eigen/Dense>
#include <vector>

struct Kernel_Array
{
	/*
	* åvéZÇµÇΩÇ¢ÉJÅ[ÉlÉãÇÕÅuPsi*Psi*Psi*w(Sophia)ÅvÇ≈Ç†ÇÈ
	*/
	Kernel_Array(Eigen::MatrixXd& c_ix, Eigen::MatrixXd& c_iy, Eigen::MatrixXd& c_iz,
		Eigen::MatrixXd& c_jx, Eigen::MatrixXd& c_jy, Eigen::MatrixXd& c_jz,
		Eigen::MatrixXd& c_kx, Eigen::MatrixXd& c_ky, Eigen::MatrixXd& c_kz,
		Eigen::MatrixXd& c_lx, Eigen::MatrixXd& c_ly, Eigen::MatrixXd& c_lz, Eigen::MatrixXd& c_o) :
		coefficient_i_x(c_ix), coefficient_i_y(c_iy), coefficient_i_z(c_iz),
		coefficient_j_x(c_jx), coefficient_j_y(c_jy), coefficient_j_z(c_jz),
		coefficient_k_x(c_kx), coefficient_k_y(c_ky), coefficient_k_z(c_kz),
		coefficient_l_x(c_lx), coefficient_l_y(c_ly), coefficient_l_z(c_lz), coefficient_o(c_o) {};

	/*
	* -2Ç©ÇÁ2Ç‹Ç≈ÇÃ4ãÊä‘ÇÃï˚íˆéÆÇÃåXÇ´Ç∆êÿï–(åWêî)Çäiî[Ç∑ÇÈçsóÒ
	*/
	Eigen::MatrixXd coefficient_i_x; //Coefficient   grid node i   division x
	Eigen::MatrixXd coefficient_i_y; //Coefficient   grid node i   division y
	Eigen::MatrixXd coefficient_i_z; //Coefficient   grid node i   division z
	Eigen::MatrixXd coefficient_j_x; //Coefficient   grid node j   division x
	Eigen::MatrixXd coefficient_j_y; //Coefficient   grid node j   division y
	Eigen::MatrixXd coefficient_j_z; //Coefficient   grid node j   division z
	Eigen::MatrixXd coefficient_k_x; //Coefficient   grid node k   division x
	Eigen::MatrixXd coefficient_k_y; //Coefficient   grid node k   division y
	Eigen::MatrixXd coefficient_k_z; //Coefficient   grid node k   division z
	Eigen::MatrixXd coefficient_l_x; //Coefficient   grid node l   division x
	Eigen::MatrixXd coefficient_l_y; //Coefficient   grid node l   division y
	Eigen::MatrixXd coefficient_l_z; //Coefficient   grid node l   division z
	Eigen::MatrixXd coefficient_o;   //Coefficient   origin point

};

Kernel_Array createKernelArray();

Kernel_Array setCoefficientSophiaX(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_ix, int gn_jx, int gn_kx);
Kernel_Array setCoefficientSophiaY(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_iy, int gn_jy, int gn_ky);
Kernel_Array setCoefficientSophiaZ(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_iz, int gn_jz, int gn_kz);

Kernel_Array setCoefficientChloeX(Kernel_Array ka, int gn_ix);
Kernel_Array setCoefficientChloeY(Kernel_Array ka, int gn_iy);
Kernel_Array setCoefficientChloeZ(Kernel_Array ka, int gn_iz);

Kernel_Array setCoefficientAriaX(Kernel_Array ka, int gn_ix, int gn_jx);
Kernel_Array setCoefficientAriaY(Kernel_Array ka, int gn_iy, int gn_jy);
Kernel_Array setCoefficientAriaZ(Kernel_Array ka, int gn_iz, int gn_jz);

Kernel_Array setCoefficientMiaX(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_ix, int gn_jx, int gn_kx, int gn_lx);
Kernel_Array setCoefficientMiaY(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_iy, int gn_jy, int gn_ky, int gn_ly);
Kernel_Array setCoefficientMiaZ(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_iz, int gn_jz, int gn_kz, int gn_lz);

std::vector<double> product(std::vector<double> a1, std::vector<double> a2);

double integral(int start_p, int stop_p, std::vector<double> array);

// For Sophia
double productSophiaX(Kernel_Array ka);
double productSophiaY(Kernel_Array ka);
double productSophiaZ(Kernel_Array ka);

// For Chloe
double productChloeX(Kernel_Array ka);
double productChloeY(Kernel_Array ka);
double productChloeZ(Kernel_Array ka);

// For Aria
double productAriaX(Kernel_Array ka);
double productAriaY(Kernel_Array ka);
double productAriaZ(Kernel_Array ka);

// For Mia
double productMiaX(Kernel_Array ka);
double productMiaY(Kernel_Array ka);
double productMiaZ(Kernel_Array ka);

#endif
