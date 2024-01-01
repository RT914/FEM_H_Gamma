#ifndef __EXPORT_H__
#define __EXPORT_H__

#include <Eigen/Dense>
#include <vector>

void exportInterpolationSophia(double X[3][3][3][3][3][3], double Y[3][3][3][3][3][3], double Z[3][3][3][3][3][3]);

void exportInterpolationChloe(double W[3][3][3]);

void exportInterpolationAria(double W[3][3][3][3][3][3]);

void exportInterpolationMia(double X[3][3][3][3][3][3][3], double Y[3][3][3][3][3][3][3], double Z[3][3][3][3][3][3][3]);

void exportMatrixP(Eigen::MatrixXd M);

void exportMatrixQ(Eigen::MatrixXd M);

void exportMatrixR(Eigen::MatrixXd M);

void exportVectorb(Eigen::VectorXd V);

void exportVectorb_Convenient(Eigen::VectorXd V);

void exportVectorc(Eigen::VectorXd V);

void exportMatrixS(Eigen::MatrixXd M);

void exportVectord(Eigen::VectorXd V);

void exportMatrixT(Eigen::MatrixXd M);

void exportVectore(Eigen::VectorXd V);

void exportVectorDelta(Eigen::VectorXd V);

void exportVectorDeltaPhi(Eigen::VectorXd V);

void exportVectorDeltaTheta(Eigen::VectorXd V);

void exportNorm(vector<double> V);

#endif
