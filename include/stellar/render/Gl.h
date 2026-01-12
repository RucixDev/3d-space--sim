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

// KHR_debug (debug message callbacks, groups, labels). Optional but extremely
// helpful for GPU debugging/profiling tools and for printing driver messages.
//
// SDL_opengl_glext.h usually provides these typedefs, but some platform/header
// combinations omit them, so provide minimal fallbacks.
#ifndef GLDEBUGPROC
typedef void (APIENTRY *GLDEBUGPROC)(GLenum source, GLenum type, GLuint id, GLenum severity,
                                    GLsizei length, const GLchar* message, const void* userParam);
#endif
#ifndef PFNGLDEBUGMESSAGECALLBACKPROC
typedef void (APIENTRYP PFNGLDEBUGMESSAGECALLBACKPROC)(GLDEBUGPROC callback, const void* userParam);
#endif
#ifndef PFNGLDEBUGMESSAGECONTROLPROC
typedef void (APIENTRYP PFNGLDEBUGMESSAGECONTROLPROC)(GLenum source, GLenum type, GLenum severity,
                                                     GLsizei count, const GLuint* ids, GLboolean enabled);
#endif
#ifndef PFNGLPUSHDEBUGGROUPPROC
typedef void (APIENTRYP PFNGLPUSHDEBUGGROUPPROC)(GLenum source, GLuint id, GLsizei length, const GLchar* message);
#endif
#ifndef PFNGLPOPDEBUGGROUPPROC
typedef void (APIENTRYP PFNGLPOPDEBUGGROUPPROC)(void);
#endif
#ifndef PFNGLOBJECTLABELPROC
typedef void (APIENTRYP PFNGLOBJECTLABELPROC)(GLenum identifier, GLuint name, GLsizei length, const GLchar* label);
#endif

// Token fallbacks for KHR_debug (some header combos omit them).
#ifndef GL_DEBUG_OUTPUT
#define GL_DEBUG_OUTPUT 0x92E0
#endif
#ifndef GL_DEBUG_OUTPUT_SYNCHRONOUS
#define GL_DEBUG_OUTPUT_SYNCHRONOUS 0x8242
#endif
#ifndef GL_DONT_CARE
#define GL_DONT_CARE 0x1100
#endif

#ifndef GL_DEBUG_SOURCE_API
#define GL_DEBUG_SOURCE_API 0x8246
#endif
#ifndef GL_DEBUG_SOURCE_WINDOW_SYSTEM
#define GL_DEBUG_SOURCE_WINDOW_SYSTEM 0x8247
#endif
#ifndef GL_DEBUG_SOURCE_SHADER_COMPILER
#define GL_DEBUG_SOURCE_SHADER_COMPILER 0x8248
#endif
#ifndef GL_DEBUG_SOURCE_THIRD_PARTY
#define GL_DEBUG_SOURCE_THIRD_PARTY 0x8249
#endif
#ifndef GL_DEBUG_SOURCE_APPLICATION
#define GL_DEBUG_SOURCE_APPLICATION 0x824A
#endif
#ifndef GL_DEBUG_SOURCE_OTHER
#define GL_DEBUG_SOURCE_OTHER 0x824B
#endif

#ifndef GL_DEBUG_TYPE_ERROR
#define GL_DEBUG_TYPE_ERROR 0x824C
#endif
#ifndef GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR
#define GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR 0x824D
#endif
#ifndef GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR
#define GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR 0x824E
#endif
#ifndef GL_DEBUG_TYPE_PORTABILITY
#define GL_DEBUG_TYPE_PORTABILITY 0x824F
#endif
#ifndef GL_DEBUG_TYPE_PERFORMANCE
#define GL_DEBUG_TYPE_PERFORMANCE 0x8250
#endif
#ifndef GL_DEBUG_TYPE_OTHER
#define GL_DEBUG_TYPE_OTHER 0x8251
#endif
#ifndef GL_DEBUG_TYPE_MARKER
#define GL_DEBUG_TYPE_MARKER 0x8268
#endif
#ifndef GL_DEBUG_TYPE_PUSH_GROUP
#define GL_DEBUG_TYPE_PUSH_GROUP 0x8269
#endif
#ifndef GL_DEBUG_TYPE_POP_GROUP
#define GL_DEBUG_TYPE_POP_GROUP 0x826A
#endif

#ifndef GL_DEBUG_SEVERITY_HIGH
#define GL_DEBUG_SEVERITY_HIGH 0x9146
#endif
#ifndef GL_DEBUG_SEVERITY_MEDIUM
#define GL_DEBUG_SEVERITY_MEDIUM 0x9147
#endif
#ifndef GL_DEBUG_SEVERITY_LOW
#define GL_DEBUG_SEVERITY_LOW 0x9148
#endif
#ifndef GL_DEBUG_SEVERITY_NOTIFICATION
#define GL_DEBUG_SEVERITY_NOTIFICATION 0x826B
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

// --- KHR_debug (optional) ---
// These will be nullptr when unsupported by the driver/context.
extern PFNGLDEBUGMESSAGECALLBACKPROC DebugMessageCallback;
extern PFNGLDEBUGMESSAGECONTROLPROC DebugMessageControl;
extern PFNGLPUSHDEBUGGROUPPROC PushDebugGroup;
extern PFNGLPOPDEBUGGROUPPROC PopDebugGroup;
extern PFNGLOBJECTLABELPROC ObjectLabel;

// Query support for KHR_debug features.
bool debugOutputSupported();
bool debugGroupsSupported();
bool debugLabelsSupported();

// Enable/disable driver debug output (glDebugMessageCallback).
//
// Notes:
//  - Many drivers only generate messages for debug contexts; for non-debug
//    contexts, it's valid for the log to stay empty.
//  - Debug groups/labels can still be useful for external GPU tools even when
//    debug output is silent.
void setDebugOutputEnabled(bool enabled, bool synchronous = false, bool includeNotifications = false);

} // namespace stellar::render::gl
