#include "stellar/render/Shader.h"

#include "stellar/render/Gl.h"

#include <algorithm>
#include <cctype>
#include <string>
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

  // Manual line scanning avoids locale/istream quirks and prevents MSVC's
  // occasional std::getline overload resolution issues (C2672) in some build
  // configurations.
  int cur = 1;
  std::size_t i = 0;

  while (i <= src.size()) {
    // Find end-of-line: \n, \r, or end.
    std::size_t e = i;
    while (e < src.size() && src[e] != '\n' && src[e] != '\r') ++e;
    const std::string_view line = src.substr(i, e - i);

    if (cur >= start && cur <= end) {
      out += (cur == lineNo) ? "> " : "  ";
      out += std::to_string(cur);
      out += ": ";
      out.append(line.data(), line.size());
      out += "\n";
    }

    if (cur > end) break;

    if (e >= src.size()) break;

    // Skip newline sequences: \n, \r, or \r\n.
    if (src[e] == '\r' && (e + 1) < src.size() && src[e + 1] == '\n') {
      i = e + 2;
    } else {
      i = e + 1;
    }
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
  clearUniformCache();
}

ShaderProgram::ShaderProgram(ShaderProgram&& o) noexcept {
  program_ = o.program_;
  uniformCache_ = std::move(o.uniformCache_);
  o.program_ = 0;
  o.uniformCache_.clear();
}

ShaderProgram& ShaderProgram::operator=(ShaderProgram&& o) noexcept {
  if (this == &o) return *this;
  if (program_) gl::DeleteProgram(program_);
  program_ = o.program_;
  uniformCache_ = std::move(o.uniformCache_);
  o.program_ = 0;
  o.uniformCache_.clear();
  return *this;
}

bool ShaderProgram::build(std::string_view vertexSrc, std::string_view fragmentSrc, std::string* outError) {
  if (program_) {
    gl::DeleteProgram(program_);
    program_ = 0;
  }
  clearUniformCache();

  if (vertexSrc.empty() || fragmentSrc.empty()) {
    if (outError) *outError = "ShaderProgram::build called with an empty shader source.";
    return false;
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
    clearUniformCache();
    return false;
  }

  return true;
}

void ShaderProgram::bind() const {
  gl::UseProgram(program_);
}

int ShaderProgram::uniformLocation(const char* name) const {
  if (!program_ || !name || !*name) return -1;

  // Tiny linear cache: avoids repeated glGetUniformLocation calls (which can be
  // surprisingly expensive) while keeping compile-time dependencies minimal.
  for (const auto& kv : uniformCache_) {
    if (kv.first == name) return kv.second;
  }

  const int loc = gl::GetUniformLocation(program_, name);
  uniformCache_.emplace_back(name, loc);
  return loc;
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
