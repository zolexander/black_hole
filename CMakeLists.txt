cmake_minimum_required(VERSION 3.21)
project(BlackHoleSim)
if(NOT CMAKE_BUILD_TYPE) 
set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()
# --- Compiler-Einstellungen --- #
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O0")

# Avoid GLFW including any OpenGL headers; we use GLAD for that
add_compile_definitions(GLFW_INCLUDE_NONE)

# --- Abhängigkeiten finden --- #
find_package(PkgConfig REQUIRED)
pkg_search_module(GLFW REQUIRED glfw3)
find_package(glm REQUIRED)
find_package(OpenGL REQUIRED) 
find_package(OpenMP REQUIRED)

# --- GLAD-Bibliothek einbinden --- #
# Prefer a packaged GLAD (e.g. via vcpkg). If not found, auto-fetch via FetchContent
find_package(glad CONFIG QUIET)
if (glad_FOUND)
    message(STATUS "Using packaged GLAD: glad::glad")
else()
    include(FetchContent)
    FetchContent_Declare(
        glad
        GIT_REPOSITORY https://github.com/Dav1dde/glad.git
        GIT_TAG v0.1.36
    )
    FetchContent_MakeAvailable(glad)
    # The fetched project defines target 'glad'. Create an alias to match packaged usage.
    if (NOT TARGET glad::glad)
        add_library(glad::glad ALIAS glad)
    endif()
endif()

# --- ImGui-Bibliothek einbinden --- #
file(GLOB IMGUI_SRC
    external/imgui/*.cpp
    external/imgui/backends/imgui_impl_glfw.cpp
    external/imgui/backends/imgui_impl_opengl3.cpp
)
add_library(imgui STATIC ${IMGUI_SRC})
target_include_directories(imgui PUBLIC external/imgui external/imgui/backends)
# Ensure ImGui uses GLAD as the OpenGL loader to avoid double-including GL headers
target_compile_definitions(imgui PUBLIC IMGUI_IMPL_OPENGL_LOADER_GLAD)
# Ensure ImGui can find GLAD headers by linking to glad
target_link_libraries(imgui PUBLIC glad::glad)

# --- 2D Lensing Executable --- #
add_executable(BlackHole2D 2D_lensing.cpp)
target_link_libraries(BlackHole2D PRIVATE
    glad::glad
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
    glad::glad
    glm::glm
    ${GLFW_LIBRARIES}
    dl
    pthread
    imgui
)
target_include_directories(BlackHole3D PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(KerrBlackHole kerr_blackhole.cpp)
target_link_libraries(KerrBlackHole PRIVATE
    glad::glad
    ${OpenMP_LIBRARIES}
    glm::glm
    ${GLFW_LIBRARIES}
    dl
    pthread
    imgui
)
target_include_directories(KerrBlackHole PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(mainKerr main.cpp kerrintegrate.cpp blackholesim.cpp engine.cpp blackhole_struct.cpp)

target_link_libraries(mainKerr PRIVATE
    glad::glad
    ${OpenMP_LIBRARIES}
    glm::glm
    ${GLFW_LIBRARIES}
    dl
    pthread
    imgui
)
target_include_directories(mainKerr PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/external/imgui
    ${CMAKE_CURRENT_SOURCE_DIR}/external/imgui/backends
)
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