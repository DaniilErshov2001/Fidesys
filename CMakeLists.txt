cmake_minimum_required(VERSION 3.10)
project(Fidesys)

set(CMAKE_CXX_STANDARD 17)



find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)



foreach(CONFIG Debug Release RelWithDebInfo MinSizeRel)
    string(TOUPPER ${CONFIG} CONFIG_UPPER)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${CONFIG_UPPER} ${CMAKE_SOURCE_DIR}/bin/${CONFIG})
endforeach()



include_directories(include)

file(GLOB SOURCE_FILES "src/*.cpp")

add_executable(fidesys ${SOURCE_FILES})

set_target_properties(fidesys PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

target_link_libraries(fidesys Eigen3::Eigen)
