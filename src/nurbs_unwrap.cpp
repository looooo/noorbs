#include "unwrap.h"


NurbsFlat::NurbsFlat(
    RowMat<double, 3> poles, Eigen::VectorXi u_knots, 
    Eigen::VectorXi v_knots, Eigen::VectorXd weights,
    int degree_u=3, int degree_v=3)
{
	this->poles = poles;
	this->base = NurbsBase(u_knots, v_knots, weights, u_degree, v_degree);
	// we have to make a map from u_i, v_i to matrix entry and the other direction
	// find the straight direction!
}

void NurbsFlat::lscm()
{

}