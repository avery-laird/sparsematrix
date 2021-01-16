cmake_minimum_required(VERSION 3.17)

project(sparse_matrix_opt)

add_compile_definitions(
        A_PATH="${A_PATH}"
        B_PATH="${B_PATH}"
        WRITE_PATH="${WRITE_PATH}"
)

add_executable(starter starter.cpp)
add_executable(eigentest eigentest.cpp)
target_include_directories(eigentest PUBLIC eigen)