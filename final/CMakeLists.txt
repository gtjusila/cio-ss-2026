cmake_minimum_required(VERSION 4.00)

project(cio-day-4)

find_package(SCIP REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)

add_executable(read_mps src/read_mps.cpp src/heur_zirounding.cpp
                        src/branch_driebeekpenalties.cpp
                        src/sepa_liftedknapsack.cpp)
target_link_libraries(read_mps PRIVATE libscip Boost::program_options)
set_property(TARGET read_mps PROPERTY INTERPROCEDURAL_OPTIMIZATION ON)
target_compile_features(read_mps PRIVATE cxx_std_23)
