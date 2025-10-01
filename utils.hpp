#ifndef BLACKHOLESIM_UTILS_HPP
#define BLACKHOLESIM_UTILS_HPP
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
namespace BlackholeSim {
    namespace Utils {
       
        std::string readFromFile(const char *filePath);
        class Shader {
            public:
                GLuint ID;
            
                // Constructor (supports optional geometry shader)
                Shader();
                Shader(const char* vertexPath, const char* fragmentPath, const char* geometryPath = nullptr);
                // Compute shader constructor
                Shader(const char* computePath);
                GLint getUniformLocation(const std::string &name) const;

                void use() const;
            
                // Uniform utilities
                void setBool (const std::string &name, bool value)  const;
                void setInt  (const std::string &name, int value)   const;
                void setFloat(const std::string &name, float value) const;
                void setVec3 (const std::string &name, const glm::vec3 &value) const;
                void setVec4 (const std::string &name, const glm::vec4 &value) const;
                void setMat4 (const std::string &name, const glm::mat4 &mat)   const;
            
            private:
                std::string loadShaderSource(const char* path);
                GLuint compileShader(const char* source, GLenum type);
            
                std::optional<std::filesystem::path> abspath_no_traversal(
                    const std::filesystem::path &basepath,
                    const std::filesystem::path &relpath);      
                // Cache for uniform locations
                mutable std::unordered_map<std::string, GLint> uniformCache;
            };
    }
}
#endif


