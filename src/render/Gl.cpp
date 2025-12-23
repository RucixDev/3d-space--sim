#include "stellar/render/Gl.h"

#include <sstream>

namespace stellar::render::gl {

PFNGLVIEWPORTPROC Viewport = nullptr;
PFNGLCLEARCOLORPROC ClearColor = nullptr;
PFNGLCLEARPROC Clear = nullptr;
PFNGLENABLEPROC Enable = nullptr;
PFNGLDISABLEPROC Disable = nullptr;
PFNGLBLENDFUNCPROC BlendFunc = nullptr;

PFNGLGENVERTEXARRAYSPROC GenVertexArrays = nullptr;
PFNGLBINDVERTEXARRAYPROC BindVertexArray = nullptr;
PFNGLDELETEVERTEXARRAYSPROC DeleteVertexArrays = nullptr;

PFNGLGENBUFFERSPROC GenBuffers = nullptr;
PFNGLBINDBUFFERPROC BindBuffer = nullptr;
PFNGLBUFFERDATAPROC BufferData = nullptr;
PFNGLDELETEBUFFERSPROC DeleteBuffers = nullptr;

PFNGLENABLEVERTEXATTRIBARRAYPROC EnableVertexAttribArray = nullptr;
PFNGLVERTEXATTRIBPOINTERPROC VertexAttribPointer = nullptr;

PFNGLCREATESHADERPROC CreateShader = nullptr;
PFNGLSHADERSOURCEPROC ShaderSource = nullptr;
PFNGLCOMPILESHADERPROC CompileShader = nullptr;
PFNGLGETSHADERIVPROC GetShaderiv = nullptr;
PFNGLGETSHADERINFOLOGPROC GetShaderInfoLog = nullptr;
PFNGLDELETESHADERPROC DeleteShader = nullptr;

PFNGLCREATEPROGRAMPROC CreateProgram = nullptr;
PFNGLATTACHSHADERPROC AttachShader = nullptr;
PFNGLLINKPROGRAMPROC LinkProgram = nullptr;
PFNGLGETPROGRAMIVPROC GetProgramiv = nullptr;
PFNGLGETPROGRAMINFOLOGPROC GetProgramInfoLog = nullptr;
PFNGLUSEPROGRAMPROC UseProgram = nullptr;
PFNGLDELETEPROGRAMPROC DeleteProgram = nullptr;

PFNGLDRAWARRAYSPROC DrawArrays = nullptr;

PFNGLGETUNIFORMLOCATIONPROC GetUniformLocation = nullptr;
PFNGLUNIFORMMATRIX4FVPROC UniformMatrix4fv = nullptr;

namespace {
template <typename T>
bool loadOne(T& outFn, GLGetProcAddressFn getProc, const char* name, std::ostringstream& missing) {
  outFn = reinterpret_cast<T>(getProc(name));
  if (!outFn) {
    missing << "  - " << name << "\n";
    return false;
  }
  return true;
}
}

bool load(GLGetProcAddressFn getProc, std::string* outError) {
  std::ostringstream missing;
  bool ok = true;

  ok &= loadOne(Viewport, getProc, "glViewport", missing);
  ok &= loadOne(ClearColor, getProc, "glClearColor", missing);
  ok &= loadOne(Clear, getProc, "glClear", missing);
  ok &= loadOne(Enable, getProc, "glEnable", missing);
  ok &= loadOne(Disable, getProc, "glDisable", missing);
  ok &= loadOne(BlendFunc, getProc, "glBlendFunc", missing);

  ok &= loadOne(GenVertexArrays, getProc, "glGenVertexArrays", missing);
  ok &= loadOne(BindVertexArray, getProc, "glBindVertexArray", missing);
  ok &= loadOne(DeleteVertexArrays, getProc, "glDeleteVertexArrays", missing);

  ok &= loadOne(GenBuffers, getProc, "glGenBuffers", missing);
  ok &= loadOne(BindBuffer, getProc, "glBindBuffer", missing);
  ok &= loadOne(BufferData, getProc, "glBufferData", missing);
  ok &= loadOne(DeleteBuffers, getProc, "glDeleteBuffers", missing);

  ok &= loadOne(EnableVertexAttribArray, getProc, "glEnableVertexAttribArray", missing);
  ok &= loadOne(VertexAttribPointer, getProc, "glVertexAttribPointer", missing);

  ok &= loadOne(CreateShader, getProc, "glCreateShader", missing);
  ok &= loadOne(ShaderSource, getProc, "glShaderSource", missing);
  ok &= loadOne(CompileShader, getProc, "glCompileShader", missing);
  ok &= loadOne(GetShaderiv, getProc, "glGetShaderiv", missing);
  ok &= loadOne(GetShaderInfoLog, getProc, "glGetShaderInfoLog", missing);
  ok &= loadOne(DeleteShader, getProc, "glDeleteShader", missing);

  ok &= loadOne(CreateProgram, getProc, "glCreateProgram", missing);
  ok &= loadOne(AttachShader, getProc, "glAttachShader", missing);
  ok &= loadOne(LinkProgram, getProc, "glLinkProgram", missing);
  ok &= loadOne(GetProgramiv, getProc, "glGetProgramiv", missing);
  ok &= loadOne(GetProgramInfoLog, getProc, "glGetProgramInfoLog", missing);
  ok &= loadOne(UseProgram, getProc, "glUseProgram", missing);
  ok &= loadOne(DeleteProgram, getProc, "glDeleteProgram", missing);

  ok &= loadOne(DrawArrays, getProc, "glDrawArrays", missing);

  ok &= loadOne(GetUniformLocation, getProc, "glGetUniformLocation", missing);
  ok &= loadOne(UniformMatrix4fv, getProc, "glUniformMatrix4fv", missing);

  if (!ok && outError) {
    std::ostringstream oss;
    oss << "Failed to load required OpenGL symbols:\n" << missing.str();
    *outError = oss.str();
  }

  return ok;
}

} // namespace stellar::render::gl
