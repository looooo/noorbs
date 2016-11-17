// LeastSquareConformalMapping + fem relaxing
// ------------------------------------------
// 
#ifndef UNWRAP_H
#define UNWRAP_H

// 1: local coordinates 2d representation  q_l_0
// 2: least square conformal map -> flat_vertices_0
// 3: local coordinates of mapped mesh q_l_1
// 4: diff in local coordinates -> forces R.B^T.B.(x1-x0)
// 5: stiffnes mat K
// 6: K.u=forces ->u
// 7: x1, y1 += w * u

#include <vector>
#include <memory>
#include <tuple>

#include <Eigen/Geometry>
#include <Eigen/IterativeLinearSolvers>

#include "nurbs.h"

typedef Eigen::SparseMatrix<double> spMat;

namespace nurbs{

typedef Eigen::Vector3d Vector3;
typedef Eigen::Vector2d Vector2;

template <typename type, unsigned int size>
using ColMat = Eigen::Matrix<type, Eigen::Dynamic, size>;

template <typename type, unsigned int size>
using RowMat = Eigen::Matrix<type, size, Eigen::Dynamic>;

class LscmRelax{
private:
    ColMat<double, 3> q_l_g;  // the position of the 3d triangles at there locale coord sys
    ColMat<double, 3> q_l_m;  // the mapped position in local coord sys

    void set_q_l_g();
    void set_q_l_m();
    void set_fixed_pins();
    void set_position(Eigen::VectorXd);
    void set_shift(Eigen::VectorXd);

    std::vector<long> new_order;
    std::vector<long> old_order;

    Eigen::Matrix<double, 3, 3> C;
    Eigen::VectorXd sol;

    std::vector<long> get_fem_fixed_pins();

    void init(
        RowMat<double, 3> vertices, 
        RowMat<long, 3> triangles,
        std::vector<long> fixed_pins);
public:
    LscmRelax(
        RowMat<double, 3> vertices, 
        RowMat<long, 3> triangles,
        std::vector<long> fixed_pins);

    LscmRelax(
        std::vector<std::array<double, 3>> vertices, 
        std::vector<std::array<long, 3>> triangles,
        std::vector<long> fixed_pins={});

    std::vector<long> fixed_pins;
    RowMat<double, 3> vertices;
    RowMat<long, 3> triangles;
    RowMat<double, 2> flat_vertices;
    ColMat<double, 1> rhs;
    Eigen::MatrixXd MATRIX;

    double nue=0.0;
    double elasticity=1.;

    void lscm();
    void relax(double);

    ColMat<double, 3> get_flat_vertices_3D();

    void rotate_by_min_bound_area();
    void transform(bool scale=false);
    
    double get_area();
    double get_flat_area();

};

struct NurbsFlat{

    NurbsFlat(
        RowMat<double, 3> poles,
        Eigen::VectorXi u_knots,
        Eigen::VectorXi v_knots,
        Eigen::VectorXd weights,
        int degree_u=3, int degree_v=3);

    NurbsBase base;
    std::vector<long> fixed_poles;
    RowMat<double, 3> poles;
    RowMat<double, 2> flat_poles;

    double nue=0.0;
    double elasticity=1.;

    void lscm();
    void relax(double);

};


}


#endif