#ifndef __NEWTONRAPHSONMETHOD_H__
#define __NEWTONRAPHSONMETHOD_H__

#include "Square.h"
#include <Eigen/Dense>

Eigen::VectorXd Newton(Square square);
Eigen::VectorXd Newton_H_onetime(Square square);
Eigen::VectorXd Newton_Convenient(Square square);
Eigen::VectorXd Newton_H(Square square);
void calInterpolationSophia();
void calInterpolationChloe();
void calInterpolationAria();
void calInterpolationMia();

Eigen::VectorXd calTheta(Square square, Eigen::VectorXd phi);

Eigen::MatrixXd calMatrixP(Square square, Eigen::VectorXd phi, Eigen::VectorXd theta); // �s��P�̌v�Z
Eigen::MatrixXd calMatrixQ(Square square, Eigen::VectorXd phi, Eigen::VectorXd theta); // �s��Q�̌v�Z
Eigen::MatrixXd calMatrixR(Square square, Eigen::VectorXd phi, Eigen::VectorXd theta); // �s��R�̌v�Z
Eigen::VectorXd calVectorb(Square square, Eigen::VectorXd phi, Eigen::VectorXd theta); // �x�N�g��b�̌v�Z
Eigen::VectorXd calVectorb_Convenient(Square square, Eigen::VectorXd phi, Eigen::VectorXd re_phi, Eigen::VectorXd theta); // newton_Convenient���\�b�h�p�̃x�N�g��b�̌v�Z
Eigen::VectorXd calVectorc(Square square, Eigen::VectorXd phi, Eigen::VectorXd theta); // �x�N�g��c�̌v�Z
Eigen::MatrixXd calMatrixT(Eigen::MatrixXd M); // �s��T�̌v�Z           
Eigen::VectorXd calVectore(Eigen::VectorXd V, Eigen::MatrixXd M); // �x�N�g��e�̌v�Z

int GridToFlat(Eigen::Vector3i grid_index); // �O���b�h�C���f�b�N�X����t���b�g�C���f�b�N�X�ւ̕ϊ�
Eigen::Vector3i FlatToGrid(int flat_index); // �t���b�g�C���f�b�N�X����O���b�h�C���f�b�N�X�ւ̕ϊ�

double HatFunction(double x); // ���l�v�Z�̍ۂɎg�p
double DifferentialHatFunction(double x); // ���l�v�Z�̍ۂɎg�p

// �敪�����@�̐��l�v�Z
double Riemann_Sum_for_Chloe(int i, double h);
double Riemann_Sum_for_Mary(Eigen::Vector3i grid_k, Eigen::VectorXd theta, double h);
double Riemann_Sum_for_Luna(Eigen::Vector3i i_minus_k, Eigen::Vector3i grid_k, Eigen::VectorXd theta, double h);
double Riemann_Sum_for_Sophia(vector<Eigen::Vector3i> v, vector<int> axis, double h);
double Riemann_Sum_Example(double h);

#endif

