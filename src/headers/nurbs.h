
#ifndef NURBS_H
#define NURBS_H

#include <Eigen/Geometry>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>

namespace nurbs{

typedef Eigen::Triplet<double> trip;
typedef Eigen::SparseMatrix<double> spMat;

struct NurbsBase3D
{
    //
    NurbsBase3D(Eigen::VectorXd u_knots, Eigen::VectorXd v_knots, Eigen::VectorXd w_knots, 
              Eigen::VectorXd weights,
              int degree_u=3, int degree_v=3, int degree_w=3);
    int degree_u;
    int degree_v;
    int degree_w;
    Eigen::VectorXd u_knots;
    Eigen::VectorXd v_knots;
    Eigen::VectorXd w_knots;

    Eigen::VectorXd weights;

    std::vector<std::function<double(double)>> u_functions;
    std::vector<std::function<double(double)>> v_functions;
    std::vector<std::function<double(double)>> w_functions;

    std::vector<std::function<double(double)>> Du_functions;
    std::vector<std::function<double(double)>> Dv_functions;
    std::vector<std::function<double(double)>> Dw_functions;

    std::vector<std::function<double(double)>> DDu_functions;
    std::vector<std::function<double(double)>> DDv_functions;
    std::vector<std::function<double(double)>> DDw_functions;

    void computeFirstDerivatives();
    void computeSecondDerivatives();

    Eigen::VectorXd getInfluenceVector(Eigen::Vector3d u);
    spMat getInfluenceMatrix(Eigen::Matrix<double, Eigen::Dynamic, 3> U);

    Eigen::VectorXd getDuVector(Eigen::Vector3d u);
    spMat getDuMatrix(Eigen::Matrix<double, Eigen::Dynamic, 3> U);

    Eigen::VectorXd getDvVector(Eigen::Vector3d u);
    spMat getDvMatrix(Eigen::Matrix<double, Eigen::Dynamic, 3> U);

    Eigen::VectorXd getDwVector(Eigen::Vector3d u);
    spMat getDwMatrix(Eigen::Matrix<double, Eigen::Dynamic, 3> U);
};


struct NurbsBase2D
{
    //
    NurbsBase2D(Eigen::VectorXd u_knots, Eigen::VectorXd v_knots,
              Eigen::VectorXd weights,
              int degree_u=3, int degree_v=3);
    int degree_u;
    int degree_v;
    Eigen::VectorXd u_knots;
    Eigen::VectorXd v_knots;
    Eigen::VectorXd weights;

    std::vector<std::function<double(double)>> u_functions;
    std::vector<std::function<double(double)>> v_functions;

    std::vector<std::function<double(double)>> Du_functions;
    std::vector<std::function<double(double)>> Dv_functions;

    std::vector<std::function<double(double)>> DDu_functions;
    std::vector<std::function<double(double)>> DDv_functions;

    void computeFirstDerivatives();
    void computeSecondDerivatives();

    Eigen::VectorXd getInfluenceVector(Eigen::Vector2d u);
    spMat getInfluenceMatrix(Eigen::Matrix<double, Eigen::Dynamic, 2> U);

    Eigen::VectorXd getDuVector(Eigen::Vector2d u);
    spMat getDuMatrix(Eigen::Matrix<double, Eigen::Dynamic, 2> U);

    Eigen::VectorXd getDvVector(Eigen::Vector2d u);
    spMat getDvMatrix(Eigen::Matrix<double, Eigen::Dynamic, 2> U);
};

struct NurbsBase1D
{
    NurbsBase1D(Eigen::VectorXd u_knots, Eigen::VectorXd weights, int degree_u=3);
    int degree_u;
    Eigen::VectorXd u_knots;
    Eigen::VectorXd weights;
    std::vector<std::function<double(double)>> u_functions;
    std::vector<std::function<double(double)>> Du_functions;
    std::vector<std::function<double(double)>> DDu_functions;

    void computeFirstDerivatives();
    void computeSecondDerivatives();

    Eigen::VectorXd getInfluenceVector(double u);
    spMat getInfluenceMatrix(Eigen::VectorXd u);

    Eigen::VectorXd getDuVector(double u);
    spMat getDuMatrix(Eigen::VectorXd u);

};

}

#endif