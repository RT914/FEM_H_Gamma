#include <Eigen/Dense>
#include <stdio.h>
#include "square.h"
#include <iostream>
#include <math.h>

double generation_random(double min, double max);

Square createSquare(int N)
{
	// —§•û‘Ì‚Ì’†SÀ•W
	Eigen::Vector3d pos;
	pos << 0.0, 0.0, 0.0;
	int one_d_point_num = N;
	double range = (N - 1) / 2.0;
	double dx = 2 * range / (N - 1);
	Square square(pos, dx, one_d_point_num);
	Eigen::Vector3d velocity;
	Eigen::Vector3d position;
	Eigen::Vector3d re_position;
	Eigen::Vector3i grid_node;
	double theta = 1.0;
	velocity << 0.0, 0.0, 0.0;

	double square_x = square.position(0);
	double square_y = square.position(1);
	double square_z = square.position(2);
	Eigen::Vector3d base_point;
	base_point << pos.x() - range, pos.y() - range, pos.z() - range;
	srand(5);

	double c = 1.0 * 0.01;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				double x = i * dx + base_point.x();// +c * ((rand() % 2) - 1);
				double y = j * dx + base_point.y();// +c * ((rand() % 2) - 1);
				double z = k * dx + base_point.z();// +c * ((rand() % 2) - 1);

				/*
				if (i == 2 && j == 1 && k == 1) {
					x = i * dx + i * i + base_point.x() + 4;
					y = j * dx + j * j + base_point.y() + 4;
					z = k * dx + k * k + base_point.z() + 4;

				}

				if (i == 1 && j == 1 && k == 2) {
					x = i * dx + i * i + base_point.x() + 1;
					y = j * dx + j * j + base_point.y() + 1;
					z = k * dx + k * k + base_point.z() + 1;

				}


				if (i == 1 && j == 2 && k == 1) {
					x = i * dx + i * i + base_point.x() + 6;
					y = j * dx + j * j + base_point.y() + 6;
					z = k * dx + k * k + base_point.z() + 6;

				}
				*/

				/*
				if (i == 1) {
					x = i * dx + i * i + base_point.x();
					y = j * dx + j * j + base_point.y();
					z = k * dx + k * k + base_point.z() - 1.5;

				}
				*/

				double re_x = i * dx + base_point.x();
				double re_y = j * dx + base_point.y();
				double re_z = k * dx + base_point.z();
				// printf("%d‰ñ–Ú\n", i);
				position << x, y, z;
				// printf("x:%f, y:%f, z:%f\n", x, y, z);
				re_position << re_x, re_y, re_z;
				// printf("x:%f, y:%f, z:%f\n", re_x, re_y, re_z);
				grid_node << i, j, k;
				Point p = Point(position, re_position, velocity, theta, grid_node);
				square.points.emplace_back(p);
			}
		}
	}

	return square;
};

//----------------------------------------------------
// ”ÍˆÍ“à‚Ì—”‚Ì¶¬
//----------------------------------------------------
double generation_random(double min, double max)
{
	return((max - min) * ((float)rand() / RAND_MAX)) + min;
};