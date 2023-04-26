#ifndef __SQUARE_H__
#define __SQUARE_H__

#include <Eigen/Dense>
#include <vector>

using namespace std;

struct Point
{
	Point(Eigen::Vector3d& pos, Eigen::Vector3d& re_pos, Eigen::Vector3d& vel, double& t, Eigen::Vector3i& gri)
		: position(pos), reference_position(re_pos), velocity(vel), theta(t), grid_node(gri) {};

	Eigen::Vector3d position;
	Eigen::Vector3d reference_position;
	Eigen::Vector3d velocity;
	double theta;
	Eigen::Vector3i grid_node;
};

struct Square
{
	Square(Eigen::Vector3d& pos, double delta_x, int one_d_point_num) : position(pos), dx(delta_x), one_dimension_point_number(one_d_point_num) {};

	Eigen::Vector3d position;
	vector<Point> points;
	double dx;
	int one_dimension_point_number;
};

Square createSquare(int particles);

#endif

