include_directories(${CMAKE_SOURCE_DIR}/src/headers)

set(HEADER_FILES ${CMAKE_SOURCE_DIR}/src/headers/unwrap.h
                 ${CMAKE_SOURCE_DIR}/src/headers/nurbs.h )

set(CPP_FILES ${CMAKE_SOURCE_DIR}/src/unwrap.cpp
              ${CMAKE_SOURCE_DIR}/src/nurbs.cpp )

add_library(nurbs SHARED ${CPP_FILES} ${HEADER_FILES})

install(TARGETS nurbs DESTINATION lib)

IF (BUILD_PYTHON_BINDINGS)
    ADD_SUBDIRECTORY("python/")
ENDIF(BUILD_PYTHON_BINDINGS)
