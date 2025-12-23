#pragma once

#include <cstdint>
#include <cstddef>
#include <string>

namespace stellar::render::gl {

// Minimal OpenGL declarations.
// We intentionally avoid external loader dependencies to keep the starter project small.

using GLenum = std::uint32_t;
using GLuint = std::uint32_t;
using GLint = std::int32_t;
using GLsizei = std::int32_t;
using GLfloat = float;
using GLchar = char;
using GLboolean = std::uint8_t;
using GLsizeiptr = std::intptr_t;
using GLintptr = std::intptr_t;

#if defined(_WIN32)
  #define STELLAR_GL_APIENTRY __stdcall
#else
  #define STELLAR_GL_APIENTRY
#endif

// Common constants (subset)
constexpr GLenum GL_FALSE = 0;
constexpr GLenum GL_TRUE = 1;

constexpr GLenum GL_COLOR_BUFFER_BIT = 0x00004000;
constexpr GLenum GL_DEPTH_BUFFER_BIT = 0x00000100;

constexpr GLenum GL_BLEND = 0x0BE2;
constexpr GLenum GL_SRC_ALPHA = 0x0302;
constexpr GLenum GL_ONE_MINUS_SRC_ALPHA = 0x0303;

constexpr GLenum GL_DEPTH_TEST = 0x0B71;
constexpr GLenum GL_PROGRAM_POINT_SIZE = 0x8642;

constexpr GLenum GL_ARRAY_BUFFER = 0x8892;
constexpr GLenum GL_DYNAMIC_DRAW = 0x88E8;

constexpr GLenum GL_FLOAT = 0x1406;
constexpr GLenum GL_POINTS = 0x0000;

constexpr GLenum GL_VERTEX_SHADER = 0x8B31;
constexpr GLenum GL_FRAGMENT_SHADER = 0x8B30;
constexpr GLenum GL_COMPILE_STATUS = 0x8B81;
constexpr GLenum GL_LINK_STATUS = 0x8B82;
constexpr GLenum GL_INFO_LOG_LENGTH = 0x8B84;

// Function pointer types
using GLGetProcAddressFn = void* (*)(const char* name);

using PFNGLVIEWPORTPROC = void (STELLAR_GL_APIENTRY *)(GLint x, GLint y, GLsizei width, GLsizei height);
using PFNGLCLEARCOLORPROC = void (STELLAR_GL_APIENTRY *)(GLfloat r, GLfloat g, GLfloat b, GLfloat a);
using PFNGLCLEARPROC = void (STELLAR_GL_APIENTRY *)(GLenum mask);
using PFNGLENABLEPROC = void (STELLAR_GL_APIENTRY *)(GLenum cap);
using PFNGLDISABLEPROC = void (STELLAR_GL_APIENTRY *)(GLenum cap);
using PFNGLBLENDFUNCPROC = void (STELLAR_GL_APIENTRY *)(GLenum sfactor, GLenum dfactor);

using PFNGLGENVERTEXARRAYSPROC = void (STELLAR_GL_APIENTRY *)(GLsizei n, GLuint* arrays);
using PFNGLBINDVERTEXARRAYPROC = void (STELLAR_GL_APIENTRY *)(GLuint array);
using PFNGLDELETEVERTEXARRAYSPROC = void (STELLAR_GL_APIENTRY *)(GLsizei n, const GLuint* arrays);

using PFNGLGENBUFFERSPROC = void (STELLAR_GL_APIENTRY *)(GLsizei n, GLuint* buffers);
using PFNGLBINDBUFFERPROC = void (STELLAR_GL_APIENTRY *)(GLenum target, GLuint buffer);
using PFNGLBUFFERDATAPROC = void (STELLAR_GL_APIENTRY *)(GLenum target, GLsizeiptr size, const void* data, GLenum usage);
using PFNGLDELETEBUFFERSPROC = void (STELLAR_GL_APIENTRY *)(GLsizei n, const GLuint* buffers);

using PFNGLENABLEVERTEXATTRIBARRAYPROC = void (STELLAR_GL_APIENTRY *)(GLuint index);
using PFNGLVERTEXATTRIBPOINTERPROC = void (STELLAR_GL_APIENTRY *)(
  GLuint index,
  GLint size,
  GLenum type,
  GLboolean normalized,
  GLsizei stride,
  const void* pointer);

using PFNGLCREATESHADERPROC = GLuint (STELLAR_GL_APIENTRY *)(GLenum type);
using PFNGLSHADERSOURCEPROC = void (STELLAR_GL_APIENTRY *)(GLuint shader, GLsizei count, const GLchar* const* string, const GLint* length);
using PFNGLCOMPILESHADERPROC = void (STELLAR_GL_APIENTRY *)(GLuint shader);
using PFNGLGETSHADERIVPROC = void (STELLAR_GL_APIENTRY *)(GLuint shader, GLenum pname, GLint* params);
using PFNGLGETSHADERINFOLOGPROC = void (STELLAR_GL_APIENTRY *)(GLuint shader, GLsizei maxLength, GLsizei* length, GLchar* infoLog);
using PFNGLDELETESHADERPROC = void (STELLAR_GL_APIENTRY *)(GLuint shader);

using PFNGLCREATEPROGRAMPROC = GLuint (STELLAR_GL_APIENTRY *)(void);
using PFNGLATTACHSHADERPROC = void (STELLAR_GL_APIENTRY *)(GLuint program, GLuint shader);
using PFNGLLINKPROGRAMPROC = void (STELLAR_GL_APIENTRY *)(GLuint program);
using PFNGLGETPROGRAMIVPROC = void (STELLAR_GL_APIENTRY *)(GLuint program, GLenum pname, GLint* params);
using PFNGLGETPROGRAMINFOLOGPROC = void (STELLAR_GL_APIENTRY *)(GLuint program, GLsizei maxLength, GLsizei* length, GLchar* infoLog);
using PFNGLUSEPROGRAMPROC = void (STELLAR_GL_APIENTRY *)(GLuint program);
using PFNGLDELETEPROGRAMPROC = void (STELLAR_GL_APIENTRY *)(GLuint program);

using PFNGLDRAWARRAYSPROC = void (STELLAR_GL_APIENTRY *)(GLenum mode, GLint first, GLsizei count);

using PFNGLGETUNIFORMLOCATIONPROC = GLint (STELLAR_GL_APIENTRY *)(GLuint program, const GLchar* name);
using PFNGLUNIFORMMATRIX4FVPROC = void (STELLAR_GL_APIENTRY *)(GLint location, GLsizei count, GLboolean transpose, const GLfloat* value);

// Loaded function pointers
extern PFNGLVIEWPORTPROC Viewport;
extern PFNGLCLEARCOLORPROC ClearColor;
extern PFNGLCLEARPROC Clear;
extern PFNGLENABLEPROC Enable;
extern PFNGLDISABLEPROC Disable;
extern PFNGLBLENDFUNCPROC BlendFunc;

extern PFNGLGENVERTEXARRAYSPROC GenVertexArrays;
extern PFNGLBINDVERTEXARRAYPROC BindVertexArray;
extern PFNGLDELETEVERTEXARRAYSPROC DeleteVertexArrays;

extern PFNGLGENBUFFERSPROC GenBuffers;
extern PFNGLBINDBUFFERPROC BindBuffer;
extern PFNGLBUFFERDATAPROC BufferData;
extern PFNGLDELETEBUFFERSPROC DeleteBuffers;

extern PFNGLENABLEVERTEXATTRIBARRAYPROC EnableVertexAttribArray;
extern PFNGLVERTEXATTRIBPOINTERPROC VertexAttribPointer;

extern PFNGLCREATESHADERPROC CreateShader;
extern PFNGLSHADERSOURCEPROC ShaderSource;
extern PFNGLCOMPILESHADERPROC CompileShader;
extern PFNGLGETSHADERIVPROC GetShaderiv;
extern PFNGLGETSHADERINFOLOGPROC GetShaderInfoLog;
extern PFNGLDELETESHADERPROC DeleteShader;

extern PFNGLCREATEPROGRAMPROC CreateProgram;
extern PFNGLATTACHSHADERPROC AttachShader;
extern PFNGLLINKPROGRAMPROC LinkProgram;
extern PFNGLGETPROGRAMIVPROC GetProgramiv;
extern PFNGLGETPROGRAMINFOLOGPROC GetProgramInfoLog;
extern PFNGLUSEPROGRAMPROC UseProgram;
extern PFNGLDELETEPROGRAMPROC DeleteProgram;

extern PFNGLDRAWARRAYSPROC DrawArrays;

extern PFNGLGETUNIFORMLOCATIONPROC GetUniformLocation;
extern PFNGLUNIFORMMATRIX4FVPROC UniformMatrix4fv;

// Loads required GL symbols. Returns true on success.
bool load(GLGetProcAddressFn getProc, std::string* outError = nullptr);

} // namespace stellar::render::gl
