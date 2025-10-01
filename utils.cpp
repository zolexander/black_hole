#include "utils.hpp"
#include <filesystem>
#include <optional>
#include <fstream>
#include <iostream>
namespace BlackholeSim {
namespace Utils {
// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================
/**
 * @brief Validates that a relative path doesn't traverse outside the base
 * directory
 * @param basepath The base directory path
 * @param relpath The relative path to validate
 * @return Optional absolute path if valid, nullopt if path traversal detected
 */
std::optional<std::filesystem::path>
Shader::abspath_no_traversal(const std::filesystem::path &basepath,
                             const std::filesystem::path &relpath) {
  const auto abspath = std::filesystem::weakly_canonical(basepath / relpath);
      // Check if the resolved path starts with the base path to prevent
      // directory traversal
      const auto index = abspath.string().rfind(basepath.string(), 0);
  if (index != 0) {
    return std::nullopt;
  }
  return abspath;
}
std::string Shader::loadShaderSource(const char* filePath) {
    if (!filePath) {
        std::cerr << "ERROR: Null file path provided\n";
        return "";
    }
    
    if (!std::filesystem::exists(filePath)) {
        std::cerr << "ERROR: File does not exist: " << filePath << "\n";
        return "";
    }
    
    if (abspath_no_traversal(std::filesystem::current_path(), filePath) == std::nullopt) {
        std::cerr << "ERROR: File path traversal detected: " << filePath << "\n";
        return "";
    }
    
    try {
        std::ifstream file(filePath, std::ios::in | std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "ERROR: Cannot open file: " << filePath << "\n";
            return "";
        }
        
        // Get file size for efficient memory allocation
        file.seekg(0, std::ios::end);
        const auto fileSize = file.tellg();
        file.seekg(0, std::ios::beg);
        
        std::string content;
        content.reserve(static_cast<size_t>(fileSize));
        content.assign(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>());
        
        return content;
    }
    catch (const std::exception& e) {
        std::cerr << "ERROR: Failed to read file " << filePath << ": " << e.what() << "\n";
        return "";
    }
}
GLuint Shader::compileShader(const char* source, GLenum type) {
    GLuint shader = glCreateShader(type);
    glShaderSource(shader, 1, &source, nullptr);
    glCompileShader(shader);

    GLint success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        char infoLog[512];
        glGetShaderInfoLog(shader, 512, nullptr, infoLog);
        std::cerr << "ERROR::SHADER::COMPILATION_FAILED (" 
                  << (type == GL_VERTEX_SHADER ? "VERTEX" :
                      type == GL_FRAGMENT_SHADER ? "FRAGMENT" : "GEOMETRY")
                  << ")\n" << infoLog << std::endl;
    }
    return shader;
}
GLint Shader::getUniformLocation(const std::string &name) const {
    if (uniformCache.find(name) != uniformCache.end()) {
        return uniformCache[name];
    }
    GLint location = glGetUniformLocation(ID, name.c_str());
    if (location == -1) {
        std::cerr << "WARNING::SHADER::UNIFORM_NOT_FOUND: " << name << std::endl;
    }

    uniformCache[name] = location;
    return location;
}
Shader::Shader() {}
// --- Graphics shader constructor (vertex+fragment+optional geometry) ---
Shader::Shader(const char* vertexPath, const char* fragmentPath, const char* geometryPath) {
    std::string vCode = loadShaderSource(vertexPath);
    std::string fCode = loadShaderSource(fragmentPath);

    GLuint vertex   = compileShader(vCode.c_str(), GL_VERTEX_SHADER);
    GLuint fragment = compileShader(fCode.c_str(), GL_FRAGMENT_SHADER);

    GLuint geometry = 0;
    if (geometryPath) {
        std::string gCode = loadShaderSource(geometryPath);
        geometry = compileShader(gCode.c_str(), GL_GEOMETRY_SHADER);
    }

    ID = glCreateProgram();
    glAttachShader(ID, vertex);
    glAttachShader(ID, fragment);
    if (geometryPath) glAttachShader(ID, geometry);
    glLinkProgram(ID);

    GLint success;
    glGetProgramiv(ID, GL_LINK_STATUS, &success);
    if (!success) {
        char infoLog[512];
        glGetProgramInfoLog(ID, 512, nullptr, infoLog);
        std::cerr << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
    }

    glDeleteShader(vertex);
    glDeleteShader(fragment);
    if (geometryPath) glDeleteShader(geometry);
}
// --- Compute shader constructor ---
Shader::Shader(const char* computePath) {
    std::string cCode = loadShaderSource(computePath);
    GLuint compute = compileShader(cCode.c_str(), GL_COMPUTE_SHADER);

    ID = glCreateProgram();
    glAttachShader(ID, compute);
    glLinkProgram(ID);

    GLint success;
    glGetProgramiv(ID, GL_LINK_STATUS, &success);
    if (!success) {
        char infoLog[512];
        glGetProgramInfoLog(ID, 512, nullptr, infoLog);
        std::cerr << "ERROR::SHADER::COMPUTE_PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
    }

    glDeleteShader(compute);
}
// --- Use ---
void Shader::use() const {
    glUseProgram(ID);
}

// --- Uniform helpers ---
void Shader::setBool(const std::string &name, bool value) const {
    glUniform1i(getUniformLocation(name), (int)value);
}
void Shader::setInt(const std::string &name, int value) const {
    glUniform1i(getUniformLocation(name), value);
}
void Shader::setFloat(const std::string &name, float value) const {
    glUniform1f(getUniformLocation(name), value);
}
void Shader::setVec3(const std::string &name, const glm::vec3 &value) const {
    glUniform3fv(getUniformLocation(name), 1, glm::value_ptr(value));
}
void Shader::setVec4(const std::string &name, const glm::vec4 &value) const {
    glUniform4fv(getUniformLocation(name), 1, glm::value_ptr(value));
}
void Shader::setMat4(const std::string &name, const glm::mat4 &mat) const {
    glUniformMatrix4fv(getUniformLocation(name), 1, GL_FALSE, glm::value_ptr(mat));
}


} // namespace Utils
} // namespace BlackholeSim