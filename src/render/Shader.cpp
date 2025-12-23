#include "stellar/render/Shader.h"

#include <vector>

namespace stellar::render {
namespace {

bool compileShader(gl::GLuint shader, const char* src, std::string* outError) {
  gl::ShaderSource(shader, 1, &src, nullptr);
  gl::CompileShader(shader);

  gl::GLint ok = 0;
  gl::GetShaderiv(shader, gl::GL_COMPILE_STATUS, &ok);
  if (ok) return true;

  gl::GLint logLen = 0;
  gl::GetShaderiv(shader, gl::GL_INFO_LOG_LENGTH, &logLen);
  std::vector<char> log(static_cast<std::size_t>(logLen) + 1u);
  gl::GetShaderInfoLog(shader, logLen, nullptr, log.data());
  if (outError) *outError = std::string(log.data());
  return false;
}

bool linkProgram(gl::GLuint program, std::string* outError) {
  gl::LinkProgram(program);

  gl::GLint ok = 0;
  gl::GetProgramiv(program, gl::GL_LINK_STATUS, &ok);
  if (ok) return true;

  gl::GLint logLen = 0;
  gl::GetProgramiv(program, gl::GL_INFO_LOG_LENGTH, &logLen);
  std::vector<char> log(static_cast<std::size_t>(logLen) + 1u);
  gl::GetProgramInfoLog(program, logLen, nullptr, log.data());
  if (outError) *outError = std::string(log.data());
  return false;
}

}

ShaderProgram::~ShaderProgram() {
  destroy();
}

bool ShaderProgram::create(const char* vertexSrc, const char* fragmentSrc, std::string* outError) {
  destroy();

  const gl::GLuint vs = gl::CreateShader(gl::GL_VERTEX_SHADER);
  const gl::GLuint fs = gl::CreateShader(gl::GL_FRAGMENT_SHADER);

  std::string err;
  if (!compileShader(vs, vertexSrc, &err)) {
    gl::DeleteShader(vs);
    gl::DeleteShader(fs);
    if (outError) *outError = "Vertex shader compile failed: " + err;
    return false;
  }
  if (!compileShader(fs, fragmentSrc, &err)) {
    gl::DeleteShader(vs);
    gl::DeleteShader(fs);
    if (outError) *outError = "Fragment shader compile failed: " + err;
    return false;
  }

  m_program = gl::CreateProgram();
  gl::AttachShader(m_program, vs);
  gl::AttachShader(m_program, fs);

  if (!linkProgram(m_program, &err)) {
    gl::DeleteShader(vs);
    gl::DeleteShader(fs);
    gl::DeleteProgram(m_program);
    m_program = 0;
    if (outError) *outError = "Program link failed: " + err;
    return false;
  }

  // Shaders can be deleted after a successful link.
  gl::DeleteShader(vs);
  gl::DeleteShader(fs);
  return true;
}

void ShaderProgram::destroy() {
  if (m_program) {
    gl::DeleteProgram(m_program);
    m_program = 0;
  }
}

void ShaderProgram::bind() const {
  gl::UseProgram(m_program);
}

gl::GLint ShaderProgram::uniformLocation(const char* name) const {
  return gl::GetUniformLocation(m_program, name);
}

} // namespace stellar::render
