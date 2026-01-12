#pragma once

// Small OpenGL function loader using SDL_GL_GetProcAddress.
// This is intentionally minimal (not a full GLAD replacement).

#include <SDL.h>
#include <SDL_opengl.h>
// SDL provides a portable OpenGL extension header that defines the PFNGL*PROC
// function pointer typedefs on platforms where system OpenGL headers may be
// incomplete (notably Windows' legacy gl.h).
#include <SDL_opengl_glext.h>

// -----------------------------------------------------------------------------
// Fallback PFN typedefs
//
// Some platform/header combinations (notably Windows + legacy OpenGL headers)
// may not provide PFNGL*PROC typedefs for core 1.1 texture entry points.
// We only need a tiny subset for this project, so define minimal fallbacks.
// -----------------------------------------------------------------------------
#ifndef PFNGLGENTEXTURESPROC
typedef void (APIENTRYP PFNGLGENTEXTURESPROC)(GLsizei n, GLuint* textures);
#endif
#ifndef PFNGLBINDTEXTUREPROC
typedef void (APIENTRYP PFNGLBINDTEXTUREPROC)(GLenum target, GLuint texture);
#endif
#ifndef PFNGLTEXIMAGE2DPROC
typedef void (APIENTRYP PFNGLTEXIMAGE2DPROC)(GLenum target, GLint level, GLint internalformat,
                                            GLsizei width, GLsizei height, GLint border,
                                            GLenum format, GLenum type, const void* pixels);
#endif
#ifndef PFNGLTEXSUBIMAGE2DPROC
typedef void (APIENTRYP PFNGLTEXSUBIMAGE2DPROC)(GLenum target, GLint level,
                                               GLint xoffset, GLint yoffset,
                                               GLsizei width, GLsizei height,
                                               GLenum format, GLenum type, const void* pixels);
#endif
#ifndef PFNGLTEXPARAMETERIPROC
typedef void (APIENTRYP PFNGLTEXPARAMETERIPROC)(GLenum target, GLenum pname, GLint param);
#endif
#ifndef PFNGLDELETETEXTURESPROC
typedef void (APIENTRYP PFNGLDELETETEXTURESPROC)(GLsizei n, const GLuint* textures);
#endif
#ifndef PFNGLGENERATEMIPMAPPROC
typedef void (APIENTRYP PFNGLGENERATEMIPMAPPROC)(GLenum target);
#endif

// Timer queries (for GPU profiling). These are optional in our loader; the
// engine will gracefully disable profiling if they are unavailable.
#ifndef PFNGLGENQUERIESPROC
typedef void (APIENTRYP PFNGLGENQUERIESPROC)(GLsizei n, GLuint* ids);
#endif
#ifndef PFNGLDELETEQUERIESPROC
typedef void (APIENTRYP PFNGLDELETEQUERIESPROC)(GLsizei n, const GLuint* ids);
#endif
#ifndef PFNGLBEGINQUERYPROC
typedef void (APIENTRYP PFNGLBEGINQUERYPROC)(GLenum target, GLuint id);
#endif
#ifndef PFNGLENDQUERYPROC
typedef void (APIENTRYP PFNGLENDQUERYPROC)(GLenum target);
#endif
#ifndef PFNGLGETQUERYOBJECTUI64VPROC
typedef void (APIENTRYP PFNGLGETQUERYOBJECTUI64VPROC)(GLuint id, GLenum pname, GLuint64* params);
#endif

// Buffer mapping (used by async texture readback / debug tooling).
#ifndef PFNGLMAPBUFFERRANGEPROC
typedef void* (APIENTRYP PFNGLMAPBUFFERRANGEPROC)(GLenum target, GLintptr offset, GLsizeiptr length, GLbitfield access);
#endif
#ifndef PFNGLUNMAPBUFFERPROC
typedef GLboolean (APIENTRYP PFNGLUNMAPBUFFERPROC)(GLenum target);
#endif

namespace stellar::render::gl {

bool load();

// Debug helpers
const char* glVersionString();

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

extern PFNGLGETUNIFORMLOCATIONPROC GetUniformLocation;
extern PFNGLUNIFORM1IPROC Uniform1i;
extern PFNGLUNIFORM1FPROC Uniform1f;
extern PFNGLUNIFORM2FPROC Uniform2f;
extern PFNGLUNIFORM3FPROC Uniform3f;
extern PFNGLUNIFORM4FPROC Uniform4f;
extern PFNGLUNIFORMMATRIX4FVPROC UniformMatrix4fv;

extern PFNGLGENVERTEXARRAYSPROC GenVertexArrays;
extern PFNGLBINDVERTEXARRAYPROC BindVertexArray;
extern PFNGLDELETEVERTEXARRAYSPROC DeleteVertexArrays;

extern PFNGLGENBUFFERSPROC GenBuffers;
extern PFNGLBINDBUFFERPROC BindBuffer;
extern PFNGLBUFFERDATAPROC BufferData;
extern PFNGLBUFFERSUBDATAPROC BufferSubData;
extern PFNGLDELETEBUFFERSPROC DeleteBuffers;

// Optional buffer mapping functions (nullptr when unsupported).
extern PFNGLMAPBUFFERRANGEPROC MapBufferRange;
extern PFNGLUNMAPBUFFERPROC UnmapBuffer;

extern PFNGLENABLEVERTEXATTRIBARRAYPROC EnableVertexAttribArray;
extern PFNGLVERTEXATTRIBPOINTERPROC VertexAttribPointer;
extern PFNGLVERTEXATTRIBDIVISORPROC VertexAttribDivisor;

extern PFNGLDRAWARRAYSINSTANCEDPROC DrawArraysInstanced;
extern PFNGLDRAWELEMENTSINSTANCEDPROC DrawElementsInstanced;

extern PFNGLACTIVETEXTUREPROC ActiveTexture;
extern PFNGLGENTEXTURESPROC GenTextures;
extern PFNGLBINDTEXTUREPROC BindTexture;
extern PFNGLTEXIMAGE2DPROC TexImage2D;
extern PFNGLTEXSUBIMAGE2DPROC TexSubImage2D;
extern PFNGLTEXPARAMETERIPROC TexParameteri;
extern PFNGLGENERATEMIPMAPPROC GenerateMipmap;
extern PFNGLDELETETEXTURESPROC DeleteTextures;

extern PFNGLGENFRAMEBUFFERSPROC GenFramebuffers;
extern PFNGLBINDFRAMEBUFFERPROC BindFramebuffer;
extern PFNGLFRAMEBUFFERTEXTURE2DPROC FramebufferTexture2D;
extern PFNGLCHECKFRAMEBUFFERSTATUSPROC CheckFramebufferStatus;
extern PFNGLDELETEFRAMEBUFFERSPROC DeleteFramebuffers;

extern PFNGLGENRENDERBUFFERSPROC GenRenderbuffers;
extern PFNGLBINDRENDERBUFFERPROC BindRenderbuffer;
extern PFNGLRENDERBUFFERSTORAGEPROC RenderbufferStorage;
extern PFNGLFRAMEBUFFERRENDERBUFFERPROC FramebufferRenderbuffer;
extern PFNGLDELETERENDERBUFFERSPROC DeleteRenderbuffers;

// Optional timer query functions (nullptr when unsupported).
extern PFNGLGENQUERIESPROC GenQueries;
extern PFNGLDELETEQUERIESPROC DeleteQueries;
extern PFNGLBEGINQUERYPROC BeginQuery;
extern PFNGLENDQUERYPROC EndQuery;
extern PFNGLGETQUERYOBJECTUI64VPROC GetQueryObjectui64v;

} // namespace stellar::render::gl
