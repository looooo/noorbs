#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/eigen.h>

#include "nurbs.h"


PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);


namespace py = pybind11;

void init_nurbs(py::module &m){
    py::class_<nurbs::NurbsBase3D>(m, "NurbsBase3D")
        .def(py::init<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, int, int, int>())
        .def("computeFirstDerivatives", &nurbs::NurbsBase3D::computeFirstDerivatives)
        .def("getInfluenceVector", &nurbs::NurbsBase3D::getInfluenceVector)
        .def("getInfluenceMatrix", &nurbs::NurbsBase3D::getInfluenceMatrix)
        .def("getDuVector", &nurbs::NurbsBase3D::getDuVector)
        .def("getDuMatrix", &nurbs::NurbsBase3D::getDuMatrix)
        .def("getDvVector", &nurbs::NurbsBase3D::getDvVector)
        .def("getDvMatrix", &nurbs::NurbsBase3D::getDvMatrix)
        .def("getDwVector", &nurbs::NurbsBase3D::getDwVector)
        .def("getDwMatrix", &nurbs::NurbsBase3D::getDwMatrix);

    py::class_<nurbs::NurbsBase2D>(m, "NurbsBase2D")
        .def(py::init<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, int, int>())
        .def("computeFirstDerivatives", &nurbs::NurbsBase2D::computeFirstDerivatives)
        .def("getInfluenceVector", &nurbs::NurbsBase2D::getInfluenceVector)
        .def("getInfluenceMatrix", &nurbs::NurbsBase2D::getInfluenceMatrix)
        .def("getDuVector", &nurbs::NurbsBase2D::getDuVector)
        .def("getDuMatrix", &nurbs::NurbsBase2D::getDuMatrix)
        .def("getDvVector", &nurbs::NurbsBase2D::getDvVector)
        .def("getDvMatrix", &nurbs::NurbsBase2D::getDvMatrix);

    py::class_<nurbs::NurbsBase1D>(m, "NurbsBase1D")
        .def(py::init<Eigen::VectorXd, Eigen::VectorXd, int>())
        .def("computeFirstDerivatives", &nurbs::NurbsBase1D::computeFirstDerivatives)
        .def("getInfluenceVector", &nurbs::NurbsBase1D::getInfluenceVector)
        .def("getInfluenceMatrix", &nurbs::NurbsBase1D::getInfluenceMatrix)
        .def("getDuVector", &nurbs::NurbsBase1D::getDuVector)
        .def("getDuMatrix", &nurbs::NurbsBase1D::getDuMatrix);
}


PYBIND11_MODULE(_nurbs, m){
    init_nurbs(m);
};