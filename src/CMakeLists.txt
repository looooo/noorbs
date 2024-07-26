include_directories(${CMAKE_SOURCE_DIR}/src/headers)

set(HEADER_FILES ${CMAKE_SOURCE_DIR}/src/headers/nurbs.h )
set(CPP_FILES ${CMAKE_SOURCE_DIR}/src/nurbs.cpp )
add_library(nurbs SHARED ${CPP_FILES} ${HEADER_FILES})

install(TARGETS nurbs DESTINATION lib)

if(BUILD_PYTHON_BINDINGS)
    add_subdirectory("python/")
endif(BUILD_PYTHON_BINDINGS)
