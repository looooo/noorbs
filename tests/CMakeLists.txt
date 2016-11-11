INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/src/headers)

# unwrap
add_executable(unwrap_test unwrap_test.cpp)
target_link_libraries(unwrap_test nurbs)
