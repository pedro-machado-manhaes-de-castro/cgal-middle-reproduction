# Table: exhaustive six-point ordering checks in README.md.
cmake_minimum_required(VERSION 3.22)
project(Middle_order_witness LANGUAGES CXX)
find_package(CGAL REQUIRED)
foreach(variant IN ITEMS original corrected)
    foreach(test IN ITEMS example contract_witness)
        add_executable(${test}_${variant} ${test}.cpp)
        target_compile_features(${test}_${variant} PRIVATE cxx_std_17)
        target_link_libraries(${test}_${variant} PRIVATE CGAL::CGAL)
        target_include_directories(${test}_${variant} BEFORE PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/../headers/${variant})
    endforeach()
endforeach()
