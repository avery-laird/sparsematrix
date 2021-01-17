cmake_minimum_required(VERSION 3.17)

project(sparse_matrix_opt)

add_compile_definitions(
        A_PATH="${A_PATH}"
        B_PATH="${B_PATH}"
        WRITE_PATH="${WRITE_PATH}"
)

find_package(OpenMP REQUIRED)

add_executable(starter starter.cpp utils.cpp matrix.cpp)
add_executable(bench bench.cpp  utils.cpp matrix.cpp)

add_executable(eigentest eigentest.cpp)
target_include_directories(eigentest PUBLIC eigen)

target_link_libraries(bench PRIVATE OpenMP::OpenMP_CXX)

add_custom_target(
        test_c ALL
        COMMAND source ${CMAKE_CURRENT_SOURCE_DIR}/venv/bin/activate && python ${CMAKE_CURRENT_SOURCE_DIR}/gen_triang.py
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        BYPRODUCTS optimized_triang_test0.cpp optimized_triang_test1.cpp optimized_triang_test2.cpp
)

add_dependencies(bench test_c)
