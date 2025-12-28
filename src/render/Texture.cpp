#include "stellar/render/Texture.h"

#include "stellar/render/Gl.h"

#include <cstdint>
#include <fstream>
#include <vector>

namespace stellar::render {

Texture2D::~Texture2D() { destroy(); }

Texture2D::Texture2D(Texture2D&& o) noexcept {
  tex_ = o.tex_;
  w_ = o.w_;
  h_ = o.h_;
  o.tex_ = 0;
  o.w_ = o.h_ = 0;
}

Texture2D& Texture2D::operator=(Texture2D&& o) noexcept {
  if (this == &o) return *this;
  destroy();
  tex_ = o.tex_;
  w_ = o.w_;
  h_ = o.h_;
  o.tex_ = 0;
  o.w_ = o.h_ = 0;
  return *this;
}

void Texture2D::destroy() {
  if (tex_) {
    gl::DeleteTextures(1, &tex_);
    tex_ = 0;
  }
  w_ = h_ = 0;
}

static std::uint16_t readU16(std::istream& in) {
  std::uint8_t b0=0,b1=0;
  in.read(reinterpret_cast<char*>(&b0), 1);
  in.read(reinterpret_cast<char*>(&b1), 1);
  return static_cast<std::uint16_t>(b0 | (b1 << 8));
}

static std::uint32_t readU32(std::istream& in) {
  std::uint8_t b[4]{};
  in.read(reinterpret_cast<char*>(b), 4);
  return (std::uint32_t)b[0] | ((std::uint32_t)b[1] << 8) | ((std::uint32_t)b[2] << 16) | ((std::uint32_t)b[3] << 24);
}

static std::int32_t readI32(std::istream& in) {
  return static_cast<std::int32_t>(readU32(in));
}

bool Texture2D::loadBMP(const std::string& path, std::string* outError) {
  destroy();

  std::ifstream in(path, std::ios::binary);
  if (!in) {
    if (outError) *outError = "Failed to open BMP: " + path;
    return false;
  }

  // BITMAPFILEHEADER
  char sig[2]{};
  in.read(sig, 2);
  if (sig[0] != 'B' || sig[1] != 'M') {
    if (outError) *outError = "Not a BMP file: " + path;
    return false;
  }

  (void)readU32(in); // file size
  (void)readU16(in); // reserved1
  (void)readU16(in); // reserved2
  const std::uint32_t dataOffset = readU32(in);

  // BITMAPINFOHEADER (assume 40 bytes)
  const std::uint32_t infoSize = readU32(in);
  if (infoSize < 40) {
    if (outError) *outError = "Unsupported BMP header size";
    return false;
  }

  const std::int32_t width = readI32(in);
  const std::int32_t heightSigned = readI32(in);
  const bool topDown = (heightSigned < 0);
  const std::int32_t height = topDown ? -heightSigned : heightSigned;

  const std::uint16_t planes = readU16(in);
  const std::uint16_t bpp = readU16(in);
  const std::uint32_t compression = readU32(in);

  if (planes != 1 || (bpp != 24 && bpp != 32) || compression != 0) {
    if (outError) *outError = "Unsupported BMP format (need uncompressed 24/32bpp)";
    return false;
  }

  // Skip rest of header
  in.seekg(static_cast<std::streamoff>(14 + infoSize), std::ios::beg);
  in.seekg(static_cast<std::streamoff>(dataOffset), std::ios::beg);

  w_ = width;
  h_ = height;

  std::vector<std::uint8_t> rgba(static_cast<std::size_t>(w_ * h_ * 4), 255);

  const std::size_t bytesPerPixel = bpp / 8;
  const std::size_t rowStride = ((w_ * bytesPerPixel + 3) / 4) * 4; // 4-byte aligned
  std::vector<std::uint8_t> row(rowStride);

  for (int y = 0; y < h_; ++y) {
    const int dstY = topDown ? y : (h_ - 1 - y);
    in.read(reinterpret_cast<char*>(row.data()), static_cast<std::streamsize>(rowStride));
    if (!in) {
      if (outError) *outError = "Unexpected EOF in BMP pixel data";
      return false;
    }

    for (int x = 0; x < w_; ++x) {
      const std::size_t src = static_cast<std::size_t>(x) * bytesPerPixel;
      const std::uint8_t b = row[src + 0];
      const std::uint8_t g = row[src + 1];
      const std::uint8_t r = row[src + 2];
      const std::uint8_t a = (bytesPerPixel == 4) ? row[src + 3] : 255;

      const std::size_t dst = (static_cast<std::size_t>(dstY) * w_ + static_cast<std::size_t>(x)) * 4;
      rgba[dst + 0] = r;
      rgba[dst + 1] = g;
      rgba[dst + 2] = b;
      rgba[dst + 3] = a;
    }
  }

  gl::GenTextures(1, &tex_);
  gl::BindTexture(GL_TEXTURE_2D, tex_);
  gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

  gl::TexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, w_, h_, 0, GL_RGBA, GL_UNSIGNED_BYTE, rgba.data());
  gl::GenerateMipmap(GL_TEXTURE_2D);

  return true;
}

void Texture2D::createChecker(int w, int h, int checkSize) {
  destroy();
  w_ = w;
  h_ = h;
  std::vector<std::uint8_t> rgba(static_cast<std::size_t>(w_ * h_ * 4), 255);

  for (int y = 0; y < h_; ++y) {
    for (int x = 0; x < w_; ++x) {
      const bool c = ((x / checkSize) + (y / checkSize)) % 2 == 0;
      const std::uint8_t v = c ? 220 : 40;
      const std::size_t i = (static_cast<std::size_t>(y) * w_ + static_cast<std::size_t>(x)) * 4;
      rgba[i+0] = v;
      rgba[i+1] = v;
      rgba[i+2] = v;
      rgba[i+3] = 255;
    }
  }

  gl::GenTextures(1, &tex_);
  gl::BindTexture(GL_TEXTURE_2D, tex_);
  gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

  gl::TexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, w_, h_, 0, GL_RGBA, GL_UNSIGNED_BYTE, rgba.data());
  gl::GenerateMipmap(GL_TEXTURE_2D);
}

void Texture2D::createRGBA(int w, int h, const void* rgbaPixels,
                           bool generateMips, bool nearestFilter, bool clampToEdge) {
  destroy();
  if (w <= 0 || h <= 0 || !rgbaPixels) return;
  w_ = w;
  h_ = h;

  gl::GenTextures(1, &tex_);
  gl::BindTexture(GL_TEXTURE_2D, tex_);

  const GLint minFilter = generateMips
    ? (nearestFilter ? GL_NEAREST_MIPMAP_NEAREST : GL_LINEAR_MIPMAP_LINEAR)
    : (nearestFilter ? GL_NEAREST : GL_LINEAR);
  const GLint magFilter = nearestFilter ? GL_NEAREST : GL_LINEAR;

  gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minFilter);
  gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magFilter);

  const GLint wrap = clampToEdge ? GL_CLAMP_TO_EDGE : GL_REPEAT;
  gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrap);
  gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrap);

  gl::TexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, w_, h_, 0, GL_RGBA, GL_UNSIGNED_BYTE, rgbaPixels);
  if (generateMips) {
    gl::GenerateMipmap(GL_TEXTURE_2D);
  }
}

void Texture2D::bind(int unit) const {
  gl::ActiveTexture(GL_TEXTURE0 + unit);
  gl::BindTexture(GL_TEXTURE_2D, tex_);
}

} // namespace stellar::render
