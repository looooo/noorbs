execute_process(COMMAND ${PYTHON_EXECUTABLE} -c 
  "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
  OUTPUT_VARIABLE python_libs OUTPUT_STRIP_TRAILING_WHITESPACE )

message(STATUS "python side dir: ${Python_SITELIB}")
message(STATUS "pybind11_INCLUDE_DIRS: ${pybind11_INCLUDE_DIRS}")

set(PYTHON_INSTALL_DIR ${python_libs}/nurbs)
include_directories(${pybind11_INCLUDE_DIRS})
include_directories(${pybind11_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/src/headers)

add_library(_nurbs SHARED _nurbs.cpp)
target_link_libraries(_nurbs nurbs)
set_target_properties(_nurbs PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/nurbs)

set_target_properties(_nurbs PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
set_target_properties(_nurbs PROPERTIES PREFIX "")


configure_file(__init__.py ${CMAKE_BINARY_DIR}/nurbs/__init__.py COPYONLY)


install(TARGETS _nurbs DESTINATION ${PYTHON_INSTALL_DIR})
install(DIRECTORY "" DESTINATION ${PYTHON_INSTALL_DIR} FILES_MATCHING PATTERN "*.py")