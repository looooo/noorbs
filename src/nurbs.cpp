#include "nurbs.h"
#include <iostream>
#include "math.h"

namespace nurbs{

double divide(double a, double b)
{
    if (fabs(b) < 10e-14)
        return 0;
    else
        return a / b;
}

// DE BOOR ALGORITHM FROM OPENGLIDER
std::function<double(double)> get_basis(int degree, int i, Eigen::VectorXd knots)
    // Return a basis_function for the given degree """
{
    if (degree == 0)
    {
        return [degree, i, knots](double t)
        {
            // The basis function for degree = 0 as per eq. 7
            double t_this = knots[i];
            double t_next = knots[i+1];
            if (t == knots[0])
                return (t_next >= t and t >= t_this);
            return (t_next >= t and t > t_this);
        };
    }
    else
    {
        return [degree, i, knots](double t)
        {
            // The basis function for degree > 0 as per eq. 8
            double out = 0.;
            double t_this = knots[i];
            double t_next = knots[i + 1];
            double t_precog  = knots[i + degree];
            double t_horizon = knots[i + degree + 1];
            if (t_this == t_horizon)
                return 0.;

            double bottom = (t_precog - t_this);
            out = divide(t - t_this, bottom) * get_basis(degree - 1, i, knots)(t);

            bottom = (t_horizon - t_next);
            out += divide(t_horizon - t, bottom) * get_basis(degree - 1, i + 1, knots)(t);
            return out;
        };
    }
};


std::function<double(double)> get_basis_derivative(int order, int degree, int i, Eigen::VectorXd knots)
    // Return the derivation of the basis function
    // order of basis function
    // degree of basis function
    // 
    // knots sequence
{
    if (order == 1)
    {
        return [degree, i, knots, order](double t)
        {
            double out = 0;
            if (not (knots[i + degree] - knots[i] == 0))
            {
                out +=  get_basis(degree - 1, i, knots)(t) *
                        degree / (knots[i + degree] - knots[i]);
            }
            if (not (knots[i + degree + 1] - knots[i + 1] == 0))
            {
                out -=  get_basis(degree - 1, i + 1, knots)(t) *
                        degree / (knots[i + degree + 1] - knots[i + 1]);
            }
            return out; 
        };
    }
    else
    {   
        return [degree, i, knots, order](double t)
        {
            double out = 0;
            if (not (knots[i + degree] - knots[i] == 0))
            {
                out +=  get_basis_derivative(order - 1, degree - 1, i, knots)(t) *
                        degree / (knots[i + degree] - knots[i]);
            }
            if (not (knots[i + degree + 1] - knots[i + 1] == 0))
            {
                out -=  get_basis_derivative(order - 1, degree - 1, i + 1, knots)(t) *
                        degree / (knots[i + degree + 1] - knots[i + 1]);
            }
            return out; 
        };
    }
}

void add_triplets(Eigen::VectorXd values, double row, std::vector<trip> &triplets)
{
    for (int i=0; i < values.size(); i++)
    {
        if (values(i) != 0.)
            triplets.push_back(trip(row, i, values(i)));
    }
}


NurbsBase3D::NurbsBase3D(Eigen::VectorXd u_knots, Eigen::VectorXd v_knots, Eigen::VectorXd w_knots,
                     Eigen::VectorXd weights, int degree_u, int degree_v, int degree_w)
{
    // assert(weights.size() == u_knots.size() * v_knots.size());
    this->u_knots = u_knots;
    this->v_knots = v_knots;
    this->w_knots = w_knots;
    this->weights = weights;
    this->degree_u = degree_u;
    this->degree_v = degree_v;
    this->degree_w = degree_w;
    for (int u_i = 0; u_i < u_knots.size() - degree_u - 1; u_i ++)
        this->u_functions.push_back(get_basis(degree_u, u_i, u_knots));
    for (int v_i = 0; v_i < v_knots.size() - degree_v - 1; v_i ++)
        this->v_functions.push_back(get_basis(degree_v, v_i, v_knots));
    for (int w_i = 0; w_i < w_knots.size() - degree_w - 1; w_i ++)
        this->w_functions.push_back(get_basis(degree_w, w_i, w_knots));
}

Eigen::VectorXd NurbsBase3D::getInfluenceVector(Eigen::Vector3d u)
{
    Eigen::VectorXd n_u, n_v, n_w;
    double sum_weights = 0;
    Eigen::VectorXd infl(this->u_functions.size() * this->v_functions.size() * this->w_functions.size());
    int i = 0;
    int u_i = 0;
    int v_i = 0;
    int w_i = 0;

    n_u.resize(this->u_functions.size());
    n_v.resize(this->v_functions.size());
    n_w.resize(this->w_functions.size());

    for (int i = 0; i < this->u_functions.size(); i ++)
        n_u[i] = this->u_functions[i](u.x());
    for (int i = 0; i < this->v_functions.size(); i ++)
        n_v[i] = this->v_functions[i](u.y());
    for (int i = 0; i < this->w_functions.size(); i ++)
        n_w[i] = this->w_functions[i](u.z());

    for (int u_i = 0; u_i < this->u_functions.size(); u_i++)
    {
        for (int v_i = 0; v_i < this->v_functions.size(); v_i++)
        {
            for (int w_i = 0; w_i < this->w_functions.size(); w_i++)
            {
                infl[i] = weights[i] * n_u[u_i] * n_v[v_i] * n_w[w_i];
                sum_weights += infl[i];
                i ++;
            }
        }
    }
    return infl / sum_weights;
}


spMat NurbsBase3D::getInfluenceMatrix(Eigen::Matrix<double, Eigen::Dynamic, 3> U)
{
    std::vector<trip> triplets;
    for (int row_index=0; row_index < U.rows(); row_index++)
        add_triplets(this->getInfluenceVector(U.row(row_index)), row_index, triplets);
    spMat mat(U.rows(), this->u_functions.size() * this->v_functions.size() * this->w_functions.size());
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
}

void NurbsBase3D::computeFirstDerivatives()
{
    for (int u_i = 0; u_i < u_functions.size(); u_i ++)
        this->Du_functions.push_back(get_basis_derivative(1, this->degree_u, u_i, this->u_knots));
    for (int v_i = 0; v_i < v_functions.size(); v_i ++)
        this->Dv_functions.push_back(get_basis_derivative(1, this->degree_v, v_i, this->v_knots));
    for (int w_i = 0; w_i < w_functions.size(); w_i ++)
        this->Dw_functions.push_back(get_basis_derivative(1, this->degree_w, w_i, this->w_knots));
}

void NurbsBase3D::computeSecondDerivatives()
{
    for (int u_i = 0; u_i < u_functions.size(); u_i ++)
        this->DDu_functions.push_back(get_basis_derivative(2, this->degree_u, u_i, this->u_knots));
    for (int v_i = 0; v_i < v_functions.size(); v_i ++)
        this->DDv_functions.push_back(get_basis_derivative(2, this->degree_v, v_i, this->v_knots));
    for (int w_i = 0; w_i < w_functions.size(); w_i ++)
        this->DDw_functions.push_back(get_basis_derivative(2, this->degree_w, w_i, this->w_knots));
}

Eigen::VectorXd NurbsBase3D::getDuVector(Eigen::Vector3d u)
{
    Eigen::VectorXd A1, A2;
    double C1, C2;
    A1.resize(this->u_functions.size() * this->v_functions.size() * this->w_functions.size());
    A2.resize(this->u_functions.size() * this->v_functions.size() * this->w_functions.size());
    A1.setZero();
    A2.setZero();
    double A3 = 0;
    double A5 = 0;
    int i = 0;
    int u_i = 0;
    int v_i = 0;
    Eigen::VectorXd n_u, n_v, n_w, Dn_u;
    n_u.resize(this->u_functions.size());
    n_v.resize(this->v_functions.size());
    n_w.resize(this->w_functions.size());
    Dn_u.resize(this->u_functions.size());
    for (int u_i=0; u_i < this->u_functions.size(); u_i++)
    {
        n_u[u_i] = this->u_functions[u_i](u.x());
        Dn_u[u_i] = this->Du_functions[u_i](u.x());
    }
    for (int v_i=0; v_i < this->v_functions.size(); v_i++)
    {
        n_v[v_i] = this->v_functions[v_i](u.y());
    }
    for (int w_i=0; w_i < this->w_functions.size(); w_i++)
    {
        n_w[w_i] = this->w_functions[w_i](u.z());
    }

    for (int u_i=0; u_i < this->u_functions.size(); u_i++)
    {
        for (int v_i=0; v_i < this->v_functions.size(); v_i++)
        {
            for (int w_i=0; w_i < this->w_functions.size(); w_i++)
            {
                C1 = weights[i] * Dn_u[u_i] * n_v[v_i] * n_w[w_i];
                C2 = weights[i] * n_u[u_i] * n_v[v_i] * n_w[w_i];
                A1[i] = C1;
                A2[i] = C2;
                A3 += C2;
                A5 += C1;
                i ++;
            }
        }
    }
    return (A1 * A3 - A2 * A5) / A3 / A3;
}

Eigen::VectorXd NurbsBase3D::getDvVector(Eigen::Vector3d u)
{
    Eigen::VectorXd A1, A2;
    double C1, C2;
    A1.resize(this->u_functions.size() * this->v_functions.size() * this->w_functions.size());
    A2.resize(this->u_functions.size() * this->v_functions.size() * this->w_functions.size());
    A1.setZero();
    A2.setZero();
    double A3 = 0;
    double A5 = 0;
    int i = 0;
    Eigen::VectorXd n_u, n_v, n_w, Dn_v;
    n_u.resize(this->u_functions.size());
    n_v.resize(this->v_functions.size());
    n_w.resize(this->w_functions.size());
    Dn_v.resize(this->u_functions.size());
    for (int u_i=0; u_i < this->u_functions.size(); u_i++)
    {
        n_u[u_i] = this->u_functions[u_i](u.x());
    }
    for (int v_i=0; v_i < this->v_functions.size(); v_i++)
    {
        n_v[v_i] = this->v_functions[v_i](u.y());
        Dn_v[v_i] = this->Dv_functions[v_i](u.y());
    }
    for (int w_i=0; w_i < this->w_functions.size(); w_i++)
    {
        n_w[w_i] = this->w_functions[w_i](u.z());
    }

    for (int u_i=0; u_i < this->u_functions.size(); u_i++)
    {
        for (int v_i=0; v_i < this->v_functions.size(); v_i++)
        {
            for (int w_i=0; w_i < this->w_functions.size(); w_i++)
            {
                C1 = weights[i] * n_u[u_i] * Dn_v[v_i] * n_w[w_i];
                C2 = weights[i] * n_u[u_i] * n_v[v_i] * n_w[w_i];
                A1[i] = C1;
                A2[i] = C2;
                A3 += C2;
                A5 += C1;
                i ++;
            }
        }
    }
    return (A1 * A3 - A2 * A5) / A3 / A3;
}

Eigen::VectorXd NurbsBase3D::getDwVector(Eigen::Vector3d u)
{
    Eigen::VectorXd A1, A2;
    double C1, C2;
    A1.resize(this->u_functions.size() * this->v_functions.size() * this->w_functions.size());
    A2.resize(this->u_functions.size() * this->v_functions.size() * this->w_functions.size());
    A1.setZero();
    A2.setZero();
    double A3 = 0;
    double A5 = 0;
    int i = 0;
    int u_i = 0;
    int v_i = 0;
    Eigen::VectorXd n_u, n_v, n_w, Dn_w;
    n_u.resize(this->u_functions.size());
    n_v.resize(this->v_functions.size());
    n_w.resize(this->w_functions.size());
    Dn_w.resize(this->w_functions.size());
    for (int u_i=0; u_i < this->u_functions.size(); u_i++)
    {
        n_u[u_i] = this->u_functions[u_i](u.x());
    }
    for (int v_i=0; v_i < this->v_functions.size(); v_i++)
    {
        n_v[v_i] = this->v_functions[v_i](u.y());
    }
    for (int w_i=0; w_i < this->w_functions.size(); w_i++)
    {
        n_w[w_i] = this->w_functions[w_i](u.z());
        Dn_w[w_i] = this->Dw_functions[w_i](u.z());
    }

    for (int u_i=0; u_i < this->u_functions.size(); u_i++)
    {
        for (int v_i=0; v_i < this->v_functions.size(); v_i++)
        {
            for (int w_i=0; w_i < this->w_functions.size(); w_i++)
            {
                C1 = weights[i] * n_u[u_i] * n_v[v_i] * Dn_w[w_i];
                C2 = weights[i] * n_u[u_i] * n_v[v_i] * n_w[w_i];
                A1[i] = C1;
                A2[i] = C2;
                A3 += C2;
                A5 += C1;
                i ++;
            }
        }
    }
    return (A1 * A3 - A2 * A5) / A3 / A3;
}


spMat NurbsBase3D::getDuMatrix(Eigen::Matrix<double, Eigen::Dynamic, 3> U)
{
    std::vector<trip> triplets;
    for (int row_index=0; row_index < U.rows(); row_index++)
        add_triplets(this->getDuVector(U.row(row_index)), row_index, triplets);
    spMat mat(U.rows(), this->u_functions.size() * this->v_functions.size() * this->w_functions.size());
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
}

spMat NurbsBase3D::getDvMatrix(Eigen::Matrix<double, Eigen::Dynamic, 3> U)
{
    std::vector<trip> triplets;
    for (int row_index=0; row_index < U.rows(); row_index++)
        add_triplets(this->getDvVector(U.row(row_index)), row_index, triplets);
    spMat mat(U.rows(), this->u_functions.size() * this->v_functions.size() * this->w_functions.size());
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
}

spMat NurbsBase3D::getDwMatrix(Eigen::Matrix<double, Eigen::Dynamic, 3> U)
{
    std::vector<trip> triplets;
    for (int row_index=0; row_index < U.rows(); row_index++)
        add_triplets(this->getDwVector(U.row(row_index)), row_index, triplets);
    spMat mat(U.rows(), this->u_functions.size() * this->v_functions.size() * this->w_functions.size());
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
}


NurbsBase2D::NurbsBase2D(Eigen::VectorXd u_knots, Eigen::VectorXd v_knots,
                     Eigen::VectorXd weights,
                     int degree_u, int degree_v)
{
    // assert(weights.size() == u_knots.size() * v_knots.size());
    this->u_knots = u_knots;
    this->v_knots = v_knots;
    this->weights = weights;
    this->degree_u = degree_u;
    this->degree_v = degree_v;
    for (int u_i = 0; u_i < u_knots.size() - degree_u - 1; u_i ++)
        this->u_functions.push_back(get_basis(degree_u, u_i, u_knots));
    for (int v_i = 0; v_i < v_knots.size() - degree_v - 1; v_i ++)
        this->v_functions.push_back(get_basis(degree_v, v_i, v_knots));
}

Eigen::VectorXd NurbsBase2D::getInfluenceVector(Eigen::Vector2d u)
{
    Eigen::VectorXd n_u, n_v;
    double sum_weights = 0;
    Eigen::VectorXd infl(this->u_functions.size() * this->v_functions.size());
    int i = 0;
    int u_i = 0;
    int v_i = 0;

    n_u.resize(this->u_functions.size());
    n_v.resize(this->v_functions.size());

    for (int i = 0; i < this->u_functions.size(); i ++)
        n_u[i] = this->u_functions[i](u.x());
    for (int i = 0; i < this->v_functions.size(); i ++)
        n_v[i] = this->v_functions[i](u.y());

    for (int u_i = 0; u_i < this->u_functions.size(); u_i++)
    {
        for (int v_i = 0; v_i < this->v_functions.size(); v_i++)
        {
            sum_weights += weights[i] * n_u[u_i] * n_v[v_i];
            infl[i] = weights[i] * n_u[u_i] * n_v[v_i];
            i ++;
        }
    }
    return infl / sum_weights;
}


spMat NurbsBase2D::getInfluenceMatrix(Eigen::Matrix<double, Eigen::Dynamic, 2> U)
{
    std::vector<trip> triplets;
    for (int row_index=0; row_index < U.rows(); row_index++)
        add_triplets(this->getInfluenceVector(U.row(row_index)), row_index, triplets);
    spMat mat(U.rows(), this->u_functions.size() * this->v_functions.size());
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
}

void NurbsBase2D::computeFirstDerivatives()
{
    for (int u_i = 0; u_i < u_functions.size(); u_i ++)
        this->Du_functions.push_back(get_basis_derivative(1, this->degree_u, u_i, this->u_knots));
    for (int v_i = 0; v_i < v_functions.size(); v_i ++)
        this->Dv_functions.push_back(get_basis_derivative(1, this->degree_v, v_i, this->v_knots));
}

void NurbsBase2D::computeSecondDerivatives()
{
    for (int u_i = 0; u_i < u_functions.size(); u_i ++)
        this->DDu_functions.push_back(get_basis_derivative(2, this->degree_u, u_i, this->u_knots));
    for (int v_i = 0; v_i < v_functions.size(); v_i ++)
        this->DDv_functions.push_back(get_basis_derivative(2, this->degree_v, v_i, this->v_knots));
}

Eigen::VectorXd NurbsBase2D::getDuVector(Eigen::Vector2d u)
{
    Eigen::VectorXd A1, A2;
    double C1, C2;
    A1.resize(this->u_functions.size() * this->v_functions.size());
    A2.resize(this->u_functions.size() * this->v_functions.size());
    double A3 = 0;
    double A5 = 0;
    int i = 0;
    int u_i = 0;
    int v_i = 0;
    Eigen::VectorXd n_u, n_v, Dn_u;
    n_u.resize(this->u_functions.size());
    Dn_u.resize(this->u_functions.size());
    n_v.resize(this->v_functions.size());
    for (int u_i=0; u_i < this->u_functions.size(); u_i++)
    {
        // std::cout << "u_i: " << u_i << " , n_u: " << n_u.size() 
        //           << " , Dn_u: " << Dn_u.size() << std::endl;
        n_u[u_i] = this->u_functions[u_i](u.x());
        Dn_u[u_i] = this->Du_functions[u_i](u.x());
    }
    for (int v_i=0; v_i < this->v_functions.size(); v_i++)
    {
        n_v[v_i] = this->v_functions[v_i](u.y());
        // std::cout << v_i << std::endl;
    }

    for (int u_i=0; u_i < this->u_functions.size(); u_i++)
    {
        for (int v_i=0; v_i < this->v_functions.size(); v_i++)
        {
            C1 = weights[i] * Dn_u[u_i] * n_v[v_i];
            C2 = weights[i] * n_u[u_i] * n_v[v_i];
            A1[i] = C1;
            A2[i] = C2;
            A3 += C2;
            A5 += C1;
            i ++;
        }
    }
    return (A1 * A3 - A2 * A5) / A3 / A3;
}

Eigen::VectorXd NurbsBase2D::getDvVector(Eigen::Vector2d u)
{
    Eigen::VectorXd A1, A2;
    double C1, C2;
    A1.resize(this->u_functions.size() * this->v_functions.size());
    A2.resize(this->u_functions.size() * this->v_functions.size());
    double A3 = 0;
    double A5 = 0;
    int i = 0;
    int u_i = 0;
    int v_i = 0;
    Eigen::VectorXd n_u, n_v, Dn_v;
    n_u.resize(this->u_functions.size());
    Dn_v.resize(this->v_functions.size());
    n_v.resize(this->v_functions.size());
    for (int u_i=0; u_i < this->u_functions.size(); u_i++)
    {
        n_u[u_i] = this->u_functions[u_i](u.x());
    }
    for (int v_i=0; v_i < this->v_functions.size(); v_i++)
    {
        n_v[v_i] = this->v_functions[v_i](u.y());
        Dn_v[v_i] = this->Dv_functions[v_i](u.y());
    }

    for (int u_i=0; u_i < this->u_functions.size(); u_i++)
    {
        for (int v_i=0; v_i < this->v_functions.size(); v_i++)
        {
            C1 = weights[i] * Dn_v[v_i] * n_u[u_i];
            C2 = weights[i] * n_v[v_i] * n_u[u_i];
            A1[i] = C1; 
            A2[i] = C2;
            A3 += C2; 
            A5 += C1;
            i ++;
        }
    }
    return (A1 * A3 - A2 * A5) / A3 / A3;
}


spMat NurbsBase2D::getDuMatrix(Eigen::Matrix<double, Eigen::Dynamic, 2> U)
{
    std::vector<trip> triplets;
    for (int row_index=0; row_index < U.rows(); row_index++)
        add_triplets(this->getDuVector(U.row(row_index)), row_index, triplets);
    spMat mat(U.rows(), this->u_functions.size() * this->v_functions.size());
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
}

spMat NurbsBase2D::getDvMatrix(Eigen::Matrix<double, Eigen::Dynamic, 2> U)
{
    std::vector<trip> triplets;
    for (int row_index=0; row_index < U.rows(); row_index++)
        add_triplets(this->getDvVector(U.row(row_index)), row_index, triplets);
    spMat mat(U.rows(), this->u_functions.size() * this->v_functions.size());
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
}

NurbsBase1D::NurbsBase1D(Eigen::VectorXd u_knots, Eigen::VectorXd weights, int degree_u)
{
    this->u_knots = u_knots;
    this->weights = weights;
    this->degree_u = degree_u;
    for (int u_i = 0; u_i < u_knots.size() - degree_u - 1; u_i ++)
        this->u_functions.push_back(get_basis(degree_u, u_i, u_knots));
}

Eigen::VectorXd NurbsBase1D::getInfluenceVector(double u)
{
    Eigen::VectorXd n_u;
    double sum_weights = 0;
    Eigen::VectorXd infl(this->u_functions.size());
    int u_i = 0;

    n_u.resize(this->u_functions.size());

    for (int i = 0; i < this->u_functions.size(); i ++)
        n_u[i] = this->u_functions[i](u);

    for (int u_i = 0; u_i < this->u_functions.size(); u_i++)
    {
        sum_weights += weights[u_i] * n_u[u_i];
        infl[u_i] = weights[u_i] * n_u[u_i];
    }
    return infl / sum_weights;
}

spMat NurbsBase1D::getInfluenceMatrix(Eigen::VectorXd u)
{
    std::vector<trip> triplets;
    for (int row_index=0; row_index < u.size(); row_index++)
        add_triplets(this->getInfluenceVector(u[row_index]), row_index, triplets);
    spMat mat(u.size(), this->u_functions.size());
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
}

void NurbsBase1D::computeFirstDerivatives()
{
    for (int u_i = 0; u_i < u_functions.size(); u_i ++)
        this->Du_functions.push_back(get_basis_derivative(1, this->degree_u, u_i, this->u_knots));
}

void NurbsBase1D::computeSecondDerivatives()
{
    for (int u_i = 0; u_i < u_functions.size(); u_i ++)
        this->DDu_functions.push_back(get_basis_derivative(2, this->degree_u, u_i, this->u_knots));
}


Eigen::VectorXd NurbsBase1D::getDuVector(double u)
{
    Eigen::VectorXd A1, A2;
    double C1, C2;
    double C3 = 0;
    double C4 = 0;
    int i = 0;
    int u_i = 0;
    A1.resize(this->u_functions.size());
    A2.resize(this->u_functions.size());
    Eigen::VectorXd n_u, Dn_u;
    n_u.resize(this->u_functions.size());
    Dn_u.resize(this->u_functions.size());

    for (int u_i=0; u_i < this->u_functions.size(); u_i++)
    {
        n_u[u_i] = this->u_functions[u_i](u);
        Dn_u[u_i] = this->Du_functions[u_i](u);
    }

    for (int u_i=0; u_i < this->Du_functions.size(); u_i++)
    {
        C1 = weights[i] * Dn_u[u_i];
        C2 = weights[i] * n_u[u_i];
        C3 += C1; 
        C4 += C2;

        A1[i] = C1; 
        A2[i] = C2;
        i ++;
    }
    return (A1 * C4 - A2 * C3) / C4 / C4 ;
}


spMat NurbsBase1D::getDuMatrix(Eigen::VectorXd U)
{
    std::vector<trip> triplets;
    for (int row_index=0; row_index < U.size(); row_index++)
        add_triplets(this->getDuVector(U[row_index]), row_index, triplets);
    spMat mat(U.size(), this->u_functions.size());
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
}

}