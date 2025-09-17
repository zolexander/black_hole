cmake_minimum_required(VERSION 3.21)
project(BlackHoleSim)

# --- Compiler-Einstellungen --- #
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Abhängigkeiten finden --- #
find_package(PkgConfig REQUIRED)
pkg_search_module(GLFW REQUIRED glfw3)
find_package(glm REQUIRED)
find_package(OpenGL REQUIRED)

# --- GLAD-Bibliothek einbinden --- #
add_library(glad STATIC external/glad/src/glad.c)
target_include_directories(glad PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/external/glad/include)

# --- ImGui-Bibliothek einbinden --- #
file(GLOB IMGUI_SRC
    external/imgui/*.cpp
    external/imgui/backends/imgui_impl_glfw.cpp
    external/imgui/backends/imgui_impl_opengl3.cpp
)
add_library(imgui STATIC ${IMGUI_SRC})
target_include_directories(imgui PUBLIC external/imgui external/imgui/backends)

# --- 2D Lensing Executable --- #
add_executable(BlackHole2D 2D_lensing.cpp)
target_link_libraries(BlackHole2D PRIVATE
    glad
    imgui
    glm::glm
    ${GLFW_LIBRARIES}
    dl
    pthread
)
target_include_directories(BlackHole2D PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/external/imgui
    ${CMAKE_CURRENT_SOURCE_DIR}/external/imgui/backends
)

# --- 3D Black Hole Executable --- #
add_executable(BlackHole3D black_hole.cpp)
target_link_libraries(BlackHole3D PRIVATE
    glad
    glm::glm
    ${GLFW_LIBRARIES}
    dl
    pthread
    imgui
)
target_include_directories(BlackHole3D PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# --- Shader-Dateien ins Ausgabeverzeichnis kopieren --- #
file(GLOB SHADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.vert"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.frag"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.comp"
)
foreach(shader ${SHADERS})
    add_custom_command(TARGET BlackHole3D POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${shader}
        $<TARGET_FILE_DIR:BlackHole3D>
    )
endforeach()

