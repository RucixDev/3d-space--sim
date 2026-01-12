#include "stellar/render/Shader.h"

#include "stellar/render/Gl.h"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <vector>

namespace stellar::render {

static const char* shaderStageName(unsigned int type) {
  switch (type) {
    case GL_VERTEX_SHADER:
      return "vertex";
    case GL_FRAGMENT_SHADER:
      return "fragment";
    default:
      return "shader";
  }
}

static int tryParseGlslErrorLine(std::string_view log) {
  // Common formats:
  //  - "ERROR: 0:63: ..."
  //  - "0:63(12): error: ..."
  //  - "0(63) : error ..."
  auto parseDigits = [](std::string_view s, std::size_t i) -> int {
    int v = 0;
    bool any = false;
    while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i]))) {
      any = true;
      v = v * 10 + (s[i] - '0');
      ++i;
    }
    return any ? v : -1;
  };

  if (auto pos = log.find("0:"); pos != std::string_view::npos) {
    int line = parseDigits(log, pos + 2);
    if (line > 0) return line;
  }
  if (auto pos = log.find("0("); pos != std::string_view::npos) {
    int line = parseDigits(log, pos + 2);
    if (line > 0) return line;
  }
  return -1;
}

static void appendSourceExcerpt(std::string& out, std::string_view src, int lineNo, int radius = 4) {
  if (lineNo <= 0) return;

  const int start = std::max(1, lineNo - radius);
  const int end = lineNo + radius;

  out += "\n--- GLSL source excerpt ---\n";

  std::istringstream iss(std::string(src));
  std::string line;
  int cur = 1;
  while (std::getline(iss, line)) {
    if (cur >= start && cur <= end) {
      out += (cur == lineNo) ? "> " : "  ";
      out += std::to_string(cur);
      out += ": ";
      out += line;
      out += "\n";
    }
    if (cur > end) break;
    ++cur;
  }
}

static unsigned int compile(unsigned int type, std::string_view src, std::string* outError) {
  unsigned int sh = gl::CreateShader(type);
  const char* csrc = src.data();
  int len = static_cast<int>(src.size());
  gl::ShaderSource(sh, 1, &csrc, &len);
  gl::CompileShader(sh);

  int status = 0;
  gl::GetShaderiv(sh, GL_COMPILE_STATUS, &status);
  if (!status) {
    int logLen = 0;
    gl::GetShaderiv(sh, GL_INFO_LOG_LENGTH, &logLen);
    std::vector<char> log(static_cast<std::size_t>(std::max(1, logLen)));
    gl::GetShaderInfoLog(sh, logLen, nullptr, log.data());

    if (outError) {
      std::string msg = std::string("Shader compile failed (") + shaderStageName(type) + "): " + log.data();
      const int lineNo = tryParseGlslErrorLine(log.data());
      appendSourceExcerpt(msg, src, lineNo);
      *outError = std::move(msg);
    }
    gl::DeleteShader(sh);
    return 0;
  }
  return sh;
}

ShaderProgram::~ShaderProgram() {
  if (program_) {
    gl::DeleteProgram(program_);
    program_ = 0;
  }
}

ShaderProgram::ShaderProgram(ShaderProgram&& o) noexcept {
  program_ = o.program_;
  o.program_ = 0;
}

ShaderProgram& ShaderProgram::operator=(ShaderProgram&& o) noexcept {
  if (this == &o) return *this;
  if (program_) gl::DeleteProgram(program_);
  program_ = o.program_;
  o.program_ = 0;
  return *this;
}

bool ShaderProgram::build(std::string_view vertexSrc, std::string_view fragmentSrc, std::string* outError) {
  if (program_) {
    gl::DeleteProgram(program_);
    program_ = 0;
  }

  std::string err;
  unsigned int vs = compile(GL_VERTEX_SHADER, vertexSrc, &err);
  if (!vs) {
    if (outError) *outError = err;
    return false;
  }
  unsigned int fs = compile(GL_FRAGMENT_SHADER, fragmentSrc, &err);
  if (!fs) {
    gl::DeleteShader(vs);
    if (outError) *outError = err;
    return false;
  }

  program_ = gl::CreateProgram();
  gl::AttachShader(program_, vs);
  gl::AttachShader(program_, fs);
  gl::LinkProgram(program_);

  gl::DeleteShader(vs);
  gl::DeleteShader(fs);

  int status = 0;
  gl::GetProgramiv(program_, GL_LINK_STATUS, &status);
  if (!status) {
    int logLen = 0;
    gl::GetProgramiv(program_, GL_INFO_LOG_LENGTH, &logLen);
    std::vector<char> log(static_cast<std::size_t>(std::max(1, logLen)));
    gl::GetProgramInfoLog(program_, logLen, nullptr, log.data());

    if (outError) *outError = std::string("Program link failed: ") + log.data();
    gl::DeleteProgram(program_);
    program_ = 0;
    return false;
  }

  return true;
}

void ShaderProgram::bind() const {
  gl::UseProgram(program_);
}

int ShaderProgram::uniformLocation(const char* name) const {
  return gl::GetUniformLocation(program_, name);
}

void ShaderProgram::setUniformMat4(const char* name, const float* mat4) const {
  const int loc = uniformLocation(name);
  if (loc >= 0) gl::UniformMatrix4fv(loc, 1, GL_FALSE, mat4);
}

void ShaderProgram::setUniform1i(const char* name, int v) const {
  const int loc = uniformLocation(name);
  if (loc >= 0) gl::Uniform1i(loc, v);
}

void ShaderProgram::setUniform1f(const char* name, float v) const {
  const int loc = uniformLocation(name);
  if (loc >= 0) gl::Uniform1f(loc, v);
}

void ShaderProgram::setUniform2f(const char* name, float x, float y) const {
  const int loc = uniformLocation(name);
  if (loc >= 0) gl::Uniform2f(loc, x, y);
}

void ShaderProgram::setUniform3f(const char* name, float x, float y, float z) const {
  const int loc = uniformLocation(name);
  if (loc >= 0) gl::Uniform3f(loc, x, y, z);
}

void ShaderProgram::setUniform4f(const char* name, float x, float y, float z, float w) const {
  const int loc = uniformLocation(name);
  if (loc >= 0) gl::Uniform4f(loc, x, y, z, w);
}

} // namespace stellar::render
