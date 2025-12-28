#pragma once

#include <string>
#include <string_view>

namespace stellar::render {

class ShaderProgram {
public:
  ShaderProgram() = default;
  ~ShaderProgram();

  ShaderProgram(const ShaderProgram&) = delete;
  ShaderProgram& operator=(const ShaderProgram&) = delete;

  ShaderProgram(ShaderProgram&&) noexcept;
  ShaderProgram& operator=(ShaderProgram&&) noexcept;

  bool build(std::string_view vertexSrc, std::string_view fragmentSrc, std::string* outError = nullptr);

  void bind() const;

  int uniformLocation(const char* name) const;
  void setUniformMat4(const char* name, const float* mat4) const;
  void setUniform1i(const char* name, int v) const;
  void setUniform1f(const char* name, float v) const;
  void setUniform2f(const char* name, float x, float y) const;
  void setUniform3f(const char* name, float x, float y, float z) const;
  void setUniform4f(const char* name, float x, float y, float z, float w) const;

  unsigned int handle() const { return program_; }

private:
  unsigned int program_{0};
};

} // namespace stellar::render
