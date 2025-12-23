#pragma once

#include "stellar/render/Gl.h"

#include <string>

namespace stellar::render {

class ShaderProgram {
public:
  ShaderProgram() = default;
  ~ShaderProgram();

  ShaderProgram(const ShaderProgram&) = delete;
  ShaderProgram& operator=(const ShaderProgram&) = delete;

  bool create(const char* vertexSrc, const char* fragmentSrc, std::string* outError = nullptr);
  void destroy();

  void bind() const;

  gl::GLuint id() const { return m_program; }

  gl::GLint uniformLocation(const char* name) const;

private:
  gl::GLuint m_program = 0;
};

} // namespace stellar::render
