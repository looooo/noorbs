set(Python_ADDITIONAL_VERSIONS ${py})
find_package(PythonInterp ${py})
find_package(PythonLibs ${py})

EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c 
  "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
  OUTPUT_VARIABLE python_libs OUTPUT_STRIP_TRAILING_WHITESPACE )

SET(PYTHON_INSTALL_DIR ${python_libs}/nurbs)

INCLUDE_DIRECTORIES( ${PYTHON_INCLUDE_DIRS})    # python.h, pybind11
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/src/headers)

add_library(_nurbs SHARED _nurbs.cpp)
target_link_libraries(_nurbs nurbs)
set_target_properties(_nurbs PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/nurbs)
# !!!!!!!!!remember
set_target_properties(_nurbs PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
set_target_properties(_nurbs PROPERTIES PREFIX "")


configure_file(__init__.py ${CMAKE_BINARY_DIR}/nurbs/__init__.py COPYONLY)


install(TARGETS _nurbs DESTINATION ${PYTHON_INSTALL_DIR})
install(DIRECTORY "" DESTINATION ${PYTHON_INSTALL_DIR} FILES_MATCHING PATTERN "*.py")