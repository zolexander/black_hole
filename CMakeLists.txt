cmake_minimum_required(VERSION 3.21)
project(BlackHoleSim)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Dependencies
# Find GLFW
find_package(PkgConfig REQUIRED)
pkg_search_module(GLFW REQUIRED glfw3)
find_package(glm REQUIRED)
find_package(OpenGL REQUIRED)

add_library(glad STATIC
    external/glad/src/glad.c
)
# Common dependencies
target_include_directories(glad
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/external/glad/include
)

# 2D Lensing executable
add_executable(BlackHole2D 2D_lensing.cpp)
target_link_libraries(BlackHole2D PRIVATE
    glad
    glm::glm
    ${GLFW_LIBRARIES}
    dl
    pthread
)
target_include_directories(BlackHole2D PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# 3D Black Hole executable
add_executable(BlackHole3D black_hole.cpp)
target_link_libraries(BlackHole3D PRIVATE
    glad
    glm::glm
    ${GLFW_LIBRARIES}
    dl
    pthread
)
target_include_directories(BlackHole3D PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Shader files (copy to output dir)
file(GLOB SHADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.vert"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.frag"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.comp"
)

foreach(shader ${SHADERS})
    add_custom_command(TARGET BlackHole3D POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${shader}
        $<TARGET_FILE_DIR:BlackHole3D
    )
endforeach()

