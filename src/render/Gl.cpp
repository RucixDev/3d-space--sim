#include "stellar/render/Gl.h"

#include "stellar/core/Log.h"

#include <string>

namespace stellar::render::gl {

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

PFNGLGETUNIFORMLOCATIONPROC GetUniformLocation = nullptr;
PFNGLUNIFORM1IPROC Uniform1i = nullptr;
PFNGLUNIFORM1FPROC Uniform1f = nullptr;
PFNGLUNIFORM2FPROC Uniform2f = nullptr;
PFNGLUNIFORM3FPROC Uniform3f = nullptr;
PFNGLUNIFORM4FPROC Uniform4f = nullptr;
PFNGLUNIFORMMATRIX4FVPROC UniformMatrix4fv = nullptr;

PFNGLGENVERTEXARRAYSPROC GenVertexArrays = nullptr;
PFNGLBINDVERTEXARRAYPROC BindVertexArray = nullptr;
PFNGLDELETEVERTEXARRAYSPROC DeleteVertexArrays = nullptr;

PFNGLGENBUFFERSPROC GenBuffers = nullptr;
PFNGLBINDBUFFERPROC BindBuffer = nullptr;
PFNGLBUFFERDATAPROC BufferData = nullptr;
PFNGLBUFFERSUBDATAPROC BufferSubData = nullptr;
PFNGLDELETEBUFFERSPROC DeleteBuffers = nullptr;

PFNGLMAPBUFFERRANGEPROC MapBufferRange = nullptr;
PFNGLUNMAPBUFFERPROC UnmapBuffer = nullptr;

PFNGLENABLEVERTEXATTRIBARRAYPROC EnableVertexAttribArray = nullptr;
PFNGLVERTEXATTRIBPOINTERPROC VertexAttribPointer = nullptr;
PFNGLVERTEXATTRIBDIVISORPROC VertexAttribDivisor = nullptr;

PFNGLDRAWARRAYSINSTANCEDPROC DrawArraysInstanced = nullptr;
PFNGLDRAWELEMENTSINSTANCEDPROC DrawElementsInstanced = nullptr;

PFNGLACTIVETEXTUREPROC ActiveTexture = nullptr;
PFNGLGENTEXTURESPROC GenTextures = nullptr;
PFNGLBINDTEXTUREPROC BindTexture = nullptr;
PFNGLTEXIMAGE2DPROC TexImage2D = nullptr;
PFNGLTEXSUBIMAGE2DPROC TexSubImage2D = nullptr;
PFNGLTEXPARAMETERIPROC TexParameteri = nullptr;
PFNGLGENERATEMIPMAPPROC GenerateMipmap = nullptr;
PFNGLDELETETEXTURESPROC DeleteTextures = nullptr;

PFNGLGENFRAMEBUFFERSPROC GenFramebuffers = nullptr;
PFNGLBINDFRAMEBUFFERPROC BindFramebuffer = nullptr;
PFNGLFRAMEBUFFERTEXTURE2DPROC FramebufferTexture2D = nullptr;
PFNGLCHECKFRAMEBUFFERSTATUSPROC CheckFramebufferStatus = nullptr;
PFNGLDELETEFRAMEBUFFERSPROC DeleteFramebuffers = nullptr;

PFNGLGENRENDERBUFFERSPROC GenRenderbuffers = nullptr;
PFNGLBINDRENDERBUFFERPROC BindRenderbuffer = nullptr;
PFNGLRENDERBUFFERSTORAGEPROC RenderbufferStorage = nullptr;
PFNGLFRAMEBUFFERRENDERBUFFERPROC FramebufferRenderbuffer = nullptr;
PFNGLDELETERENDERBUFFERSPROC DeleteRenderbuffers = nullptr;

PFNGLGENQUERIESPROC GenQueries = nullptr;
PFNGLDELETEQUERIESPROC DeleteQueries = nullptr;
PFNGLBEGINQUERYPROC BeginQuery = nullptr;
PFNGLENDQUERYPROC EndQuery = nullptr;
PFNGLGETQUERYOBJECTUI64VPROC GetQueryObjectui64v = nullptr;

// --- KHR_debug (optional) ---
PFNGLDEBUGMESSAGECALLBACKPROC DebugMessageCallback = nullptr;
PFNGLDEBUGMESSAGECONTROLPROC DebugMessageControl = nullptr;
PFNGLPUSHDEBUGGROUPPROC PushDebugGroup = nullptr;
PFNGLPOPDEBUGGROUPPROC PopDebugGroup = nullptr;
PFNGLOBJECTLABELPROC ObjectLabel = nullptr;

template <class T>
static T loadProc(const char* name) {
  return reinterpret_cast<T>(SDL_GL_GetProcAddress(name));
}

bool load() {
  CreateShader = loadProc<PFNGLCREATESHADERPROC>("glCreateShader");
  ShaderSource = loadProc<PFNGLSHADERSOURCEPROC>("glShaderSource");
  CompileShader = loadProc<PFNGLCOMPILESHADERPROC>("glCompileShader");
  GetShaderiv = loadProc<PFNGLGETSHADERIVPROC>("glGetShaderiv");
  GetShaderInfoLog = loadProc<PFNGLGETSHADERINFOLOGPROC>("glGetShaderInfoLog");
  DeleteShader = loadProc<PFNGLDELETESHADERPROC>("glDeleteShader");

  CreateProgram = loadProc<PFNGLCREATEPROGRAMPROC>("glCreateProgram");
  AttachShader = loadProc<PFNGLATTACHSHADERPROC>("glAttachShader");
  LinkProgram = loadProc<PFNGLLINKPROGRAMPROC>("glLinkProgram");
  GetProgramiv = loadProc<PFNGLGETPROGRAMIVPROC>("glGetProgramiv");
  GetProgramInfoLog = loadProc<PFNGLGETPROGRAMINFOLOGPROC>("glGetProgramInfoLog");
  UseProgram = loadProc<PFNGLUSEPROGRAMPROC>("glUseProgram");
  DeleteProgram = loadProc<PFNGLDELETEPROGRAMPROC>("glDeleteProgram");

  GetUniformLocation = loadProc<PFNGLGETUNIFORMLOCATIONPROC>("glGetUniformLocation");
  Uniform1i = loadProc<PFNGLUNIFORM1IPROC>("glUniform1i");
  Uniform1f = loadProc<PFNGLUNIFORM1FPROC>("glUniform1f");
  Uniform2f = loadProc<PFNGLUNIFORM2FPROC>("glUniform2f");
  Uniform3f = loadProc<PFNGLUNIFORM3FPROC>("glUniform3f");
  Uniform4f = loadProc<PFNGLUNIFORM4FPROC>("glUniform4f");
  UniformMatrix4fv = loadProc<PFNGLUNIFORMMATRIX4FVPROC>("glUniformMatrix4fv");

  GenVertexArrays = loadProc<PFNGLGENVERTEXARRAYSPROC>("glGenVertexArrays");
  BindVertexArray = loadProc<PFNGLBINDVERTEXARRAYPROC>("glBindVertexArray");
  DeleteVertexArrays = loadProc<PFNGLDELETEVERTEXARRAYSPROC>("glDeleteVertexArrays");

  GenBuffers = loadProc<PFNGLGENBUFFERSPROC>("glGenBuffers");
  BindBuffer = loadProc<PFNGLBINDBUFFERPROC>("glBindBuffer");
  BufferData = loadProc<PFNGLBUFFERDATAPROC>("glBufferData");
  BufferSubData = loadProc<PFNGLBUFFERSUBDATAPROC>("glBufferSubData");
  DeleteBuffers = loadProc<PFNGLDELETEBUFFERSPROC>("glDeleteBuffers");

  // Optional buffer mapping (debug tooling).
  MapBufferRange = loadProc<PFNGLMAPBUFFERRANGEPROC>("glMapBufferRange");
  if (!MapBufferRange) MapBufferRange = loadProc<PFNGLMAPBUFFERRANGEPROC>("glMapBufferRangeARB");
  UnmapBuffer = loadProc<PFNGLUNMAPBUFFERPROC>("glUnmapBuffer");
  if (!UnmapBuffer) UnmapBuffer = loadProc<PFNGLUNMAPBUFFERPROC>("glUnmapBufferARB");

  EnableVertexAttribArray = loadProc<PFNGLENABLEVERTEXATTRIBARRAYPROC>("glEnableVertexAttribArray");
  VertexAttribPointer = loadProc<PFNGLVERTEXATTRIBPOINTERPROC>("glVertexAttribPointer");
  VertexAttribDivisor = loadProc<PFNGLVERTEXATTRIBDIVISORPROC>("glVertexAttribDivisor");

  DrawArraysInstanced = loadProc<PFNGLDRAWARRAYSINSTANCEDPROC>("glDrawArraysInstanced");
  DrawElementsInstanced = loadProc<PFNGLDRAWELEMENTSINSTANCEDPROC>("glDrawElementsInstanced");

  ActiveTexture = loadProc<PFNGLACTIVETEXTUREPROC>("glActiveTexture");
  if (!ActiveTexture) {
    // Some drivers expose this via ARB extension entry points.
    ActiveTexture = loadProc<PFNGLACTIVETEXTUREPROC>("glActiveTextureARB");
  }
  GenTextures = loadProc<PFNGLGENTEXTURESPROC>("glGenTextures");
  BindTexture = loadProc<PFNGLBINDTEXTUREPROC>("glBindTexture");
  TexImage2D = loadProc<PFNGLTEXIMAGE2DPROC>("glTexImage2D");
  TexSubImage2D = loadProc<PFNGLTEXSUBIMAGE2DPROC>("glTexSubImage2D");
  TexParameteri = loadProc<PFNGLTEXPARAMETERIPROC>("glTexParameteri");
  GenerateMipmap = loadProc<PFNGLGENERATEMIPMAPPROC>("glGenerateMipmap");
  if (!GenerateMipmap) {
    // Some older drivers expose mipmap generation as an extension entry point.
    GenerateMipmap = loadProc<PFNGLGENERATEMIPMAPPROC>("glGenerateMipmapEXT");
  }
  DeleteTextures = loadProc<PFNGLDELETETEXTURESPROC>("glDeleteTextures");

  GenFramebuffers = loadProc<PFNGLGENFRAMEBUFFERSPROC>("glGenFramebuffers");
  BindFramebuffer = loadProc<PFNGLBINDFRAMEBUFFERPROC>("glBindFramebuffer");
  FramebufferTexture2D = loadProc<PFNGLFRAMEBUFFERTEXTURE2DPROC>("glFramebufferTexture2D");
  CheckFramebufferStatus = loadProc<PFNGLCHECKFRAMEBUFFERSTATUSPROC>("glCheckFramebufferStatus");
  DeleteFramebuffers = loadProc<PFNGLDELETEFRAMEBUFFERSPROC>("glDeleteFramebuffers");

  GenRenderbuffers = loadProc<PFNGLGENRENDERBUFFERSPROC>("glGenRenderbuffers");
  BindRenderbuffer = loadProc<PFNGLBINDRENDERBUFFERPROC>("glBindRenderbuffer");
  RenderbufferStorage = loadProc<PFNGLRENDERBUFFERSTORAGEPROC>("glRenderbufferStorage");
  FramebufferRenderbuffer = loadProc<PFNGLFRAMEBUFFERRENDERBUFFERPROC>("glFramebufferRenderbuffer");
  DeleteRenderbuffers = loadProc<PFNGLDELETERENDERBUFFERSPROC>("glDeleteRenderbuffers");

  // Optional timer query functions (for GPU profiling).
  GenQueries = loadProc<PFNGLGENQUERIESPROC>("glGenQueries");
  DeleteQueries = loadProc<PFNGLDELETEQUERIESPROC>("glDeleteQueries");
  BeginQuery = loadProc<PFNGLBEGINQUERYPROC>("glBeginQuery");
  EndQuery = loadProc<PFNGLENDQUERYPROC>("glEndQuery");
  GetQueryObjectui64v = loadProc<PFNGLGETQUERYOBJECTUI64VPROC>("glGetQueryObjectui64v");

  if (!GenQueries) GenQueries = loadProc<PFNGLGENQUERIESPROC>("glGenQueriesARB");
  if (!DeleteQueries) DeleteQueries = loadProc<PFNGLDELETEQUERIESPROC>("glDeleteQueriesARB");
  if (!BeginQuery) BeginQuery = loadProc<PFNGLBEGINQUERYPROC>("glBeginQueryARB");
  if (!EndQuery) EndQuery = loadProc<PFNGLENDQUERYPROC>("glEndQueryARB");
  if (!GetQueryObjectui64v) GetQueryObjectui64v = loadProc<PFNGLGETQUERYOBJECTUI64VPROC>("glGetQueryObjectui64vARB");

  // Optional KHR_debug functions (debug output, debug groups, object labels).
  // Per spec, OpenGL entry points are non-suffixed; we also try KHR-suffixed
  // names for extra robustness across loader quirks.
  DebugMessageCallback = loadProc<PFNGLDEBUGMESSAGECALLBACKPROC>("glDebugMessageCallback");
  if (!DebugMessageCallback) DebugMessageCallback = loadProc<PFNGLDEBUGMESSAGECALLBACKPROC>("glDebugMessageCallbackKHR");
  DebugMessageControl = loadProc<PFNGLDEBUGMESSAGECONTROLPROC>("glDebugMessageControl");
  if (!DebugMessageControl) DebugMessageControl = loadProc<PFNGLDEBUGMESSAGECONTROLPROC>("glDebugMessageControlKHR");
  PushDebugGroup = loadProc<PFNGLPUSHDEBUGGROUPPROC>("glPushDebugGroup");
  if (!PushDebugGroup) PushDebugGroup = loadProc<PFNGLPUSHDEBUGGROUPPROC>("glPushDebugGroupKHR");
  PopDebugGroup = loadProc<PFNGLPOPDEBUGGROUPPROC>("glPopDebugGroup");
  if (!PopDebugGroup) PopDebugGroup = loadProc<PFNGLPOPDEBUGGROUPPROC>("glPopDebugGroupKHR");
  ObjectLabel = loadProc<PFNGLOBJECTLABELPROC>("glObjectLabel");
  if (!ObjectLabel) ObjectLabel = loadProc<PFNGLOBJECTLABELPROC>("glObjectLabelKHR");


  // Fallback for entry points that may be exported directly by the OpenGL library
  // rather than returned by SDL_GL_GetProcAddress (common on Windows for GL 1.1/1.3 funcs).
  if (!GenTextures)   GenTextures   = reinterpret_cast<PFNGLGENTEXTURESPROC>(&::glGenTextures);
  if (!BindTexture)   BindTexture   = reinterpret_cast<PFNGLBINDTEXTUREPROC>(&::glBindTexture);
  if (!TexImage2D)    TexImage2D    = reinterpret_cast<PFNGLTEXIMAGE2DPROC>(&::glTexImage2D);
  if (!TexSubImage2D) TexSubImage2D = reinterpret_cast<PFNGLTEXSUBIMAGE2DPROC>(&::glTexSubImage2D);
  if (!TexParameteri) TexParameteri = reinterpret_cast<PFNGLTEXPARAMETERIPROC>(&::glTexParameteri);
  if (!DeleteTextures)DeleteTextures= reinterpret_cast<PFNGLDELETETEXTURESPROC>(&::glDeleteTextures);

  // Note: Do NOT fall back to &::glActiveTexture or &::glGenerateMipmap on Windows.
  // The system OpenGL import library (opengl32.dll) only guarantees prototypes/exports for
  // OpenGL 1.1, so referencing these symbols can fail to compile even if the runtime supports
  // them. These must be loaded via SDL_GL_GetProcAddress.

  const bool ok =
      CreateShader && ShaderSource && CompileShader && GetShaderiv && GetShaderInfoLog && DeleteShader &&
      CreateProgram && AttachShader && LinkProgram && GetProgramiv && GetProgramInfoLog && UseProgram && DeleteProgram &&
      GetUniformLocation && Uniform1i && Uniform1f && Uniform2f && Uniform3f && Uniform4f && UniformMatrix4fv &&
      GenVertexArrays && BindVertexArray && DeleteVertexArrays &&
      GenBuffers && BindBuffer && BufferData && BufferSubData && DeleteBuffers &&
      EnableVertexAttribArray && VertexAttribPointer && VertexAttribDivisor &&
      DrawArraysInstanced && DrawElementsInstanced &&
      ActiveTexture && GenTextures && BindTexture && TexImage2D && TexSubImage2D && TexParameteri && GenerateMipmap && DeleteTextures &&
      GenFramebuffers && BindFramebuffer && FramebufferTexture2D && CheckFramebufferStatus && DeleteFramebuffers &&
      GenRenderbuffers && BindRenderbuffer && RenderbufferStorage && FramebufferRenderbuffer && DeleteRenderbuffers;

  if (!ok) {
    stellar::core::log(stellar::core::LogLevel::Error, "OpenGL loader: missing required functions (context too old?)");
  }

  return ok;
}

const char* glVersionString() {
  const auto* s = glGetString(GL_VERSION);
  return s ? reinterpret_cast<const char*>(s) : "(null)";
}

// --- KHR_debug helpers ---

namespace {

const char* dbgSourceStr(GLenum src) {
  switch (src) {
    case GL_DEBUG_SOURCE_API: return "API";
    case GL_DEBUG_SOURCE_WINDOW_SYSTEM: return "WindowSystem";
    case GL_DEBUG_SOURCE_SHADER_COMPILER: return "ShaderCompiler";
    case GL_DEBUG_SOURCE_THIRD_PARTY: return "ThirdParty";
    case GL_DEBUG_SOURCE_APPLICATION: return "Application";
    case GL_DEBUG_SOURCE_OTHER: return "Other";
    default: return "Unknown";
  }
}

const char* dbgTypeStr(GLenum ty) {
  switch (ty) {
    case GL_DEBUG_TYPE_ERROR: return "Error";
    case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR: return "Deprecated";
    case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR: return "Undefined";
    case GL_DEBUG_TYPE_PORTABILITY: return "Portability";
    case GL_DEBUG_TYPE_PERFORMANCE: return "Performance";
    case GL_DEBUG_TYPE_MARKER: return "Marker";
    case GL_DEBUG_TYPE_PUSH_GROUP: return "PushGroup";
    case GL_DEBUG_TYPE_POP_GROUP: return "PopGroup";
    case GL_DEBUG_TYPE_OTHER: return "Other";
    default: return "Unknown";
  }
}

core::LogLevel dbgSeverityToLogLevel(GLenum sev) {
  switch (sev) {
    case GL_DEBUG_SEVERITY_HIGH: return core::LogLevel::Error;
    case GL_DEBUG_SEVERITY_MEDIUM: return core::LogLevel::Warn;
    case GL_DEBUG_SEVERITY_LOW: return core::LogLevel::Info;
    case GL_DEBUG_SEVERITY_NOTIFICATION: return core::LogLevel::Debug;
    default: return core::LogLevel::Info;
  }
}

void APIENTRY glDebugCallback(GLenum source, GLenum type, GLuint id, GLenum severity,
                             GLsizei length, const GLchar* message, const void* userParam) {
  (void)length;
  (void)userParam;

  // Many drivers already include a prefix; keep ours compact but structured.
  const core::LogLevel lvl = dbgSeverityToLogLevel(severity);

  std::string s;
  s.reserve(128);
  s += "GL[";
  s += dbgSourceStr(source);
  s += "/";
  s += dbgTypeStr(type);
  s += "] id=";
  s += std::to_string((unsigned)id);
  s += ": ";
  if (message) s += message;

  core::log(lvl, s);
}

} // namespace

bool debugOutputSupported() {
  return DebugMessageCallback != nullptr;
}

bool debugGroupsSupported() {
  return PushDebugGroup != nullptr && PopDebugGroup != nullptr;
}

bool debugLabelsSupported() {
  return ObjectLabel != nullptr;
}

void setDebugOutputEnabled(bool enabled, bool synchronous, bool includeNotifications) {
  if (!DebugMessageCallback) return;

  if (!enabled) {
    // Disable callback + state.
    DebugMessageCallback(nullptr, nullptr);
    glDisable(GL_DEBUG_OUTPUT);
    glDisable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
    return;
  }

  glEnable(GL_DEBUG_OUTPUT);
  if (synchronous) glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
  else glDisable(GL_DEBUG_OUTPUT_SYNCHRONOUS);

  // Route messages through our logger.
  DebugMessageCallback(&glDebugCallback, nullptr);

  // Filter out notification spam by default (can be toggled by callers).
  if (DebugMessageControl) {
    // Enable everything first, then optionally disable notifications.
    DebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, nullptr, GL_TRUE);
    if (!includeNotifications) {
      DebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DEBUG_SEVERITY_NOTIFICATION, 0, nullptr, GL_FALSE);
    }
  }
}

} // namespace stellar::render::gl
