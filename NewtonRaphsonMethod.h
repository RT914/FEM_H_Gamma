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

Eigen::MatrixXd calMatrixP(Square square, Eigen::VectorXd phi, Eigen::VectorXd theta); // 行列Pの計算
Eigen::MatrixXd calMatrixQ(Square square, Eigen::VectorXd phi, Eigen::VectorXd theta); // 行列Qの計算
Eigen::MatrixXd calMatrixR(Square square, Eigen::VectorXd phi, Eigen::VectorXd theta); // 行列Rの計算
Eigen::VectorXd calVectorb(Square square, Eigen::VectorXd phi, Eigen::VectorXd theta); // ベクトルbの計算
Eigen::VectorXd calVectorb_Convenient(Square square, Eigen::VectorXd phi, Eigen::VectorXd re_phi, Eigen::VectorXd theta); // newton_Convenientメソッド用のベクトルbの計算
Eigen::VectorXd calVectorc(Square square, Eigen::VectorXd phi, Eigen::VectorXd theta); // ベクトルcの計算
Eigen::MatrixXd calMatrixT(Eigen::MatrixXd M); // 行列Tの計算           
Eigen::VectorXd calVectore(Eigen::VectorXd V, Eigen::MatrixXd M); // ベクトルeの計算

int GridToFlat(Eigen::Vector3i grid_index); // グリッドインデックスからフラットインデックスへの変換
Eigen::Vector3i FlatToGrid(int flat_index); // フラットインデックスからグリッドインデックスへの変換

double HatFunction(double x); // 数値計算の際に使用
double DifferentialHatFunction(double x); // 数値計算の際に使用

// 区分級数法の数値計算
double Riemann_Sum_for_Chloe(int i, double h);
double Riemann_Sum_for_Mary(int i, double h);
double Riemann_Sum_for_Luna(int i, int j, double h);
double Riemann_Sum_for_Sophia(vector<Eigen::Vector3i> v, vector<int> axis, double h);

#endif

