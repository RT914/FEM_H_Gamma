#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>

#include "fem.h"

void exportInterpolationSophia(double X[3][3][3][3][3][3], double Y[3][3][3][3][3][3], double Z[3][3][3][3][3][3]) {

	std::string output_csv_file_name = "csv/Sophia.csv";
	std::ofstream data_file(output_csv_file_name);

	// Header
	data_file << "ix-ƒÌx" << "," << "iy-ƒÌy" << "," << "iz-ƒÌz" << ",";
	data_file << "jx-ƒÌx" << "," << "jy-ƒÌy" << "," << "jz-ƒÌz" << ",";
	data_file << "kx-ƒÌx" << "," << "ky-ƒÌy" << "," << "kz-ƒÌz" << ",";

	// kernel index
	for (int wi = 0; wi < dimensions; wi++) {
		for (int wj = 0; wj < dimensions; wj++) {
			for (int wk = 0; wk < dimensions; wk++) {

				data_file << "w" << wi + 1 << "w" << wj + 1 << "w" << wk + 1 << ",";

			}
		}
	}

	data_file << std::endl;


	// Grid index & Value
	// grid node					
	for (int ix = 0; ix < dimensions; ix++) {
		for (int iy = 0; iy < dimensions; iy++) {
			for (int iz = 0; iz < dimensions; iz++) {

				// case i
				// -1 ` 1
				int gn_i_x = ix - 1;
				int gn_i_y = iy - 1;
				int gn_i_z = iz - 1;

				for (int jx = 0; jx < dimensions; jx++) {
					for (int jy = 0; jy < dimensions; jy++) {
						for (int jz = 0; jz < dimensions; jz++) {

							// case j
							// -1 ` 1
							int gn_j_x = jx - 1;
							int gn_j_y = jy - 1;
							int gn_j_z = jz - 1;

							for (int kx = 0; kx < dimensions; kx++) {
								for (int ky = 0; ky < dimensions; ky++) {
									for (int kz = 0; kz < dimensions; kz++) {

										// case k
										// -1 ` 1

										int gn_k_x = kx - 1;
										int gn_k_y = ky - 1;
										int gn_k_z = kz - 1;


										data_file << gn_i_x << "," << gn_i_y << "," << gn_i_z << ","
											<< gn_j_x << "," << gn_j_y << "," << gn_j_z << ","
											<< gn_k_x << "," << gn_k_y << "," << gn_k_z << ",";


										for (int wi = 0; wi < dimensions; wi++) {
											for (int wj = 0; wj < dimensions; wj++) {
												for (int wk = 0; wk < dimensions; wk++) {

													data_file << X[wi][wj][wk][ix][jx][kx] * Y[wi][wj][wk][iy][jy][ky] * Z[wi][wj][wk][iz][jz][kz] << ",";

												}
											}
										}



										data_file << std::endl;


									}
								}
							}


						}
					}
				}



			}
		}
	}

	data_file.close();

};

void exportInterpolationChloe(double W[3][3][3]) {

	std::string output_csv_file_name = "csv/Chloe.csv";

	std::ofstream data_file(output_csv_file_name);

	// Header
	data_file << "ix-ƒÌx" << "," << "iy-ƒÌy" << "," << "iz-ƒÌz" << ",";

	data_file << std::endl;


	// Grid index & Value
	// grid node
	for (int ix = 0; ix < dimensions; ix++) {
		for (int iy = 0; iy < dimensions; iy++) {
			for (int iz = 0; iz < dimensions; iz++) {

				// -1 ` 1
				int gn_ix = ix - 1;
				int gn_iy = iy - 1;
				int gn_iz = iz - 1;


				data_file << gn_ix << "," << gn_iy << "," << gn_iz << ",";
				data_file << W[ix][iy][iz] << ",";

				data_file << std::endl;
			}
		}
	}

	data_file.close();

}

void exportInterpolationAria(double W[3][3][3][3][3][3]) {

	std::string output_csv_file_name = "csv/Aria.csv";

	std::ofstream data_file(output_csv_file_name);

	// Header
	data_file << "ix-ƒÌx" << "," << "iy-ƒÌy" << "," << "iz-ƒÌz" << ",";
	data_file << "jx-ƒÌx" << "," << "jy-ƒÌy" << "," << "jz-ƒÌz" << ",";

	data_file << std::endl;


	// Grid index & Value
	// grid node
	for (int ix = 0; ix < dimensions; ix++) {
		for (int iy = 0; iy < dimensions; iy++) {
			for (int iz = 0; iz < dimensions; iz++) {

				// -1 ` 1
				int gn_ix = ix - 1;
				int gn_iy = iy - 1;
				int gn_iz = iz - 1;


				for (int jx = 0; jx < dimensions; jx++) {
					for (int jy = 0; jy < dimensions; jy++) {
						for (int jz = 0; jz < dimensions; jz++) {

							// -1 ` 1
							int gn_jx = jx - 1;
							int gn_jy = jy - 1;
							int gn_jz = jz - 1;


							data_file << gn_ix << "," << gn_iy << "," << gn_iz << ",";
							data_file << gn_jx << "," << gn_jy << "," << gn_jz << ",";
							data_file << W[ix][iy][iz][jx][jy][jz] << ",";

							data_file << std::endl;
						}
					}
				}


			}
		}
	}

	data_file.close();

}

void exportInterpolationMia(double X[3][3][3][3][3][3][3], double Y[3][3][3][3][3][3][3], double Z[3][3][3][3][3][3][3]) {

	std::string output_csv_file_name = "csv/Mia.csv";
	std::ofstream data_file(output_csv_file_name);

	// Header
	data_file << "ix-ƒÌx" << "," << "iy-ƒÌy" << "," << "iz-ƒÌz" << ",";
	data_file << "jx-ƒÌx" << "," << "jy-ƒÌy" << "," << "jz-ƒÌz" << ",";
	data_file << "kx-ƒÌx" << "," << "ky-ƒÌy" << "," << "kz-ƒÌz" << ",";
	data_file << "lx-ƒÌx" << "," << "ly-ƒÌy" << "," << "lz-ƒÌz" << ",";

	// kernel index
	for (int wi = 0; wi < dimensions; wi++) {
		for (int wj = 0; wj < dimensions; wj++) {
			for (int wk = 0; wk < dimensions; wk++) {

				data_file << "w" << wi + 1 << "w" << wj + 1 << "w" << wk + 1 << ",";

			}
		}
	}

	data_file << std::endl;


	// Grid index & Value
	// grid node					
	for (int ix = 0; ix < dimensions; ix++) {
		for (int iy = 0; iy < dimensions; iy++) {
			for (int iz = 0; iz < dimensions; iz++) {

				// case i
				// -1 ` 1
				int gn_i_x = ix - 1;
				int gn_i_y = iy - 1;
				int gn_i_z = iz - 1;

				for (int jx = 0; jx < dimensions; jx++) {
					for (int jy = 0; jy < dimensions; jy++) {
						for (int jz = 0; jz < dimensions; jz++) {

							// case j
							// -1 ` 1
							int gn_j_x = jx - 1;
							int gn_j_y = jy - 1;
							int gn_j_z = jz - 1;

							for (int kx = 0; kx < dimensions; kx++) {
								for (int ky = 0; ky < dimensions; ky++) {
									for (int kz = 0; kz < dimensions; kz++) {

										// case k
										// -1 ` 1

										int gn_k_x = kx - 1;
										int gn_k_y = ky - 1;
										int gn_k_z = kz - 1;

										for (int lx = 0; lx < dimensions; lx++) {
											for (int ly = 0; ly < dimensions; ly++) {
												for (int lz = 0; lz < dimensions; lz++) {

													// case k
													// -1 ` 1

													int gn_l_x = lx - 1;
													int gn_l_y = ly - 1;
													int gn_l_z = lz - 1;

													data_file << gn_i_x << "," << gn_i_y << "," << gn_i_z << ","
														<< gn_j_x << "," << gn_j_y << "," << gn_j_z << ","
														<< gn_k_x << "," << gn_k_y << "," << gn_k_z << ","
														<< gn_l_x << "," << gn_l_y << "," << gn_l_z << ",";

													for (int wi = 0; wi < dimensions; wi++) {
														for (int wj = 0; wj < dimensions; wj++) {
															for (int wk = 0; wk < dimensions; wk++) {

																data_file << X[wi][wj][wk][ix][jx][kx][lx] * Y[wi][wj][wk][iy][jy][ky][ly] * Z[wi][wj][wk][iz][jz][kz][lz] << ",";

															}
														}
													}

													data_file << std::endl;

												}
											}
										}


									}
								}
							}


						}
					}
				}



			}
		}
	}

	data_file.close();

};

void exportMatrixP(Eigen::MatrixXd M) {
	std::string output_csv_file_name = "csv/MatrixP.csv";
	std::ofstream data_file(output_csv_file_name);
	// std::cout << "Matrix P ‚Ìs” : " << M.rows() << std::endl;
	// std::cout << "Matrix P ‚Ì—ñ” : " << M.cols() << std::endl;

	// Header
	data_file << "" << ",";
	for (int i = 0; i < M.rows(); i++) {
		data_file << i + 1 << ",";
	}
	data_file << std::endl;

	for (int i = 0; i < M.cols(); i++) {
		data_file << i + 1 << ",";
		for (int j = 0; j < M.rows(); j++) {
			data_file << M(j, i) << ",";
		}
		data_file << std::endl;
	}

	data_file.close();
}

void exportMatrixQ(Eigen::MatrixXd M) {
	std::string output_csv_file_name = "csv/MatrixQ.csv";
	std::ofstream data_file(output_csv_file_name);
	// std::cout << "Matrix Q ‚Ìs” : " << M.rows() << std::endl;
	// std::cout << "Matrix Q ‚Ì—ñ” : " << M.cols() << std::endl;

	// Header
	data_file << "" << ",";
	for (int i = 0; i < M.rows(); i++) {
		data_file << i + 1 << ",";
	}
	data_file << std::endl;

	for (int i = 0; i < M.cols(); i++) {
		data_file << i + 1 << ",";
		for (int j = 0; j < M.rows(); j++) {
			data_file << M(j, i) << ",";
		}
		data_file << std::endl;
	}

	data_file.close();
}

void exportMatrixR(Eigen::MatrixXd M) {
	std::string output_csv_file_name = "csv/MatrixR.csv";
	std::ofstream data_file(output_csv_file_name);
	// std::cout << "Matrix R ‚Ìs” : " << M.rows() << std::endl;
	// std::cout << "Matrix R ‚Ì—ñ” : " << M.cols() << std::endl;

	// Header
	data_file << "" << ",";
	for (int i = 0; i < M.rows(); i++) {
		data_file << i + 1 << ",";
	}
	data_file << std::endl;

	for (int i = 0; i < M.cols(); i++) {
		data_file << i + 1 << ",";
		for (int j = 0; j < M.rows(); j++) {
			data_file << M(j, i) << ",";
		}
		data_file << std::endl;
	}

	data_file.close();
}

void exportVectorb(Eigen::VectorXd V) {
	std::string output_csv_file_name = "csv/Vectorb.csv";
	std::ofstream data_file(output_csv_file_name);
	// std::cout << "Vector b ‚Ì‘å‚«‚³ : " << V.size() << std::endl;

	// Header
	for (int i = 0; i < V.size(); i++) {
		data_file << i + 1 << ",";
		data_file << V(i) << ",";
		data_file << std::endl;
	}

	data_file.close();
}

void exportVectorb_Convenient(Eigen::VectorXd V) {
	std::string output_csv_file_name = "csv/Vectorb_Convenient.csv";
	std::ofstream data_file(output_csv_file_name);
	// std::cout << "Vector b_Convenient ‚Ì‘å‚«‚³ : " << V.size() << std::endl;

	// Header
	for (int i = 0; i < V.size(); i++) {
		data_file << i + 1 << ",";
		data_file << V(i) << ",";
		data_file << std::endl;
	}

	data_file.close();
}

void exportVectorc(Eigen::VectorXd V) {
	std::string output_csv_file_name = "csv/Vectorc.csv";
	std::ofstream data_file(output_csv_file_name);
	// std::cout << "Vector c ‚Ì‘å‚«‚³ : " << V.size() << std::endl;

	// Header
	for (int i = 0; i < V.size(); i++) {
		data_file << i + 1 << ",";
		data_file << V(i) << ",";
		data_file << std::endl;
	}

	data_file.close();
}

void exportMatrixS(Eigen::MatrixXd M) {
	std::string output_csv_file_name = "csv/MatrixS.csv";
	std::ofstream data_file(output_csv_file_name);
	// std::cout << "Matrix S ‚Ìs” : " << M.rows() << std::endl;
	// std::cout << "Matrix S ‚Ì—ñ” : " << M.cols() << std::endl;

	// Header
	data_file << "" << ",";
	for (int i = 0; i < M.rows(); i++) {
		data_file << i + 1 << ",";
	}
	data_file << std::endl;

	for (int i = 0; i < M.cols(); i++) {
		data_file << i + 1 << ",";
		for (int j = 0; j < M.rows(); j++) {
			data_file << M(j, i) << ",";
		}
		data_file << std::endl;
	}

	data_file.close();
}

void exportVectord(Eigen::VectorXd V) {
	std::string output_csv_file_name = "csv/Vectord.csv";
	std::ofstream data_file(output_csv_file_name);
	// std::cout << "Vector d ‚Ì‘å‚«‚³ : " << V.size() << std::endl;

	// Header
	for (int i = 0; i < V.size(); i++) {
		data_file << i + 1 << ",";
		data_file << V(i) << ",";
		data_file << std::endl;
	}

	data_file.close();
}

void exportMatrixT(Eigen::MatrixXd M) {
	std::string output_csv_file_name = "csv/MatrixT.csv";
	std::ofstream data_file(output_csv_file_name);
	// std::cout << "Matrix S ‚Ìs” : " << M.rows() << std::endl;
	// std::cout << "Matrix S ‚Ì—ñ” : " << M.cols() << std::endl;

	// Header
	data_file << "" << ",";
	for (int i = 0; i < M.rows(); i++) {
		data_file << i + 1 << ",";
	}
	data_file << std::endl;

	for (int i = 0; i < M.cols(); i++) {
		data_file << i + 1 << ",";
		for (int j = 0; j < M.rows(); j++) {
			data_file << M(j, i) << ",";
		}
		data_file << std::endl;
	}

	data_file.close();
}

void exportVectore(Eigen::VectorXd V) {
	std::string output_csv_file_name = "csv/Vectore.csv";
	std::ofstream data_file(output_csv_file_name);
	// std::cout << "Vector e ‚Ì‘å‚«‚³ : " << V.size() << std::endl;

	// Header
	for (int i = 0; i < V.size(); i++) {
		data_file << i + 1 << ",";
		data_file << V(i) << ",";
		data_file << std::endl;
	}

	data_file.close();
}

void exportVectorDelta(Eigen::VectorXd V) {
	std::string output_csv_file_name = "csv/VectorDelta.csv";
	std::ofstream data_file(output_csv_file_name);
	// std::cout << "Vector Delta ‚Ì‘å‚«‚³ : " << V.size() << std::endl;

	// Header
	for (int i = 0; i < V.size(); i++) {
		data_file << i + 1 << ",";
		data_file << V(i) << ",";
		data_file << std::endl;
	}

	data_file.close();
}

void exportVectorDeltaPhi(Eigen::VectorXd V) {
	std::string output_csv_file_name = "csv/VectorDeltaPhi.csv";
	std::ofstream data_file(output_csv_file_name);
	// std::cout << "Vector DeltaPhi ‚Ì‘å‚«‚³ : " << V.size() << std::endl;

	// Header
	for (int i = 0; i < V.size(); i++) {
		data_file << i + 1 << ",";
		data_file << V(i) << ",";
		data_file << std::endl;
	}

	data_file.close();
}

void exportVectorDeltaTheta(Eigen::VectorXd V) {
	std::string output_csv_file_name = "csv/VectorDeltaTheta.csv";
	std::ofstream data_file(output_csv_file_name);
	// std::cout << "Vector DeltaTheta ‚Ì‘å‚«‚³ : " << V.size() << std::endl;

	// Header
	for (int i = 0; i < V.size(); i++) {
		data_file << i + 1 << ",";
		data_file << V(i) << ",";
		data_file << std::endl;
	}

	data_file.close();
}

void exportNorm(std::vector <double> V) {
	std::string output_csv_file_name = "csv/Norm.csv";
	std::ofstream data_file(output_csv_file_name);
	// std::cout << "Vector DeltaTheta ‚Ì‘å‚«‚³ : " << V.size() << std::endl;

	// Header
	for (int i = 0; i < V.size(); i++) {
		data_file << i + 1 << ",";
		data_file << V[i] << ",";
		data_file << std::endl;
	}

	data_file.close();
}