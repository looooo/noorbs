#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/eigen.h>

#include "unwrap.h"
#include "nurbs.h"


PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);


namespace py = pybind11;

void init_nurbs(py::module &m){
    py::class_<nurbs::LscmRelax>(m, "LscmRelax")
        .def(py::init<std::vector<std::array<double, 3>>, std::vector<std::array<long, 3>>, std::vector<long>>())
        .def(py::init<nurbs::ColMat<double, 3>, nurbs::ColMat<long, 3>, std::vector<long>>())
        .def("lscm", &nurbs::LscmRelax::lscm)
        .def("relax", &nurbs::LscmRelax::relax)
        .def("rotate_by_min_bound_area", &nurbs::LscmRelax::rotate_by_min_bound_area)
        .def("transform", &nurbs::LscmRelax::transform)
        .def_readonly("rhs", &nurbs::LscmRelax::rhs)
        .def_readonly("MATRIX", &nurbs::LscmRelax::MATRIX)
        .def_property_readonly("area", &nurbs::LscmRelax::get_area)
        .def_property_readonly("flat_area", &nurbs::LscmRelax::get_flat_area)
        .def_property_readonly("flat_vertices", [](nurbs::LscmRelax& L){return L.flat_vertices.transpose();}, py::return_value_policy::copy)
        .def_property_readonly("flat_vertices_3D", &nurbs::LscmRelax::get_flat_vertices_3D);

    py::class_<nurbs::NurbsBase>(m, "NurbsBase")
        .def(py::init<Eigen::VectorXi, Eigen::VectorXi, Eigen::VectorXd, int, int>())
        .def("computeFirstDerivatives", &nurbs::NurbsBase::computeFirstDerivatives)
        .def("getInfluenceVector", &nurbs::NurbsBase::getInfluenceVector)
        .def("getInfluenceMatrix", &nurbs::NurbsBase::getInfluenceMatrix)
        .def("getDuVector", &nurbs::NurbsBase::getDuVector)
        .def("getDuMatrix", &nurbs::NurbsBase::getDuMatrix)
        .def("getDvVector", &nurbs::NurbsBase::getDvVector)
        .def("getDvMatrix", &nurbs::NurbsBase::getDvMatrix);
}


PYBIND11_PLUGIN(_nurbs){
    py::module m("_nurbs");
    init_nurbs(m);
    return m.ptr();
};