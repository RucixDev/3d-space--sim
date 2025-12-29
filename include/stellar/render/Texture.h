#pragma once

#include <string>

namespace stellar::render {

class Texture2D {
public:
  Texture2D() = default;
  ~Texture2D();

  Texture2D(const Texture2D&) = delete;
  Texture2D& operator=(const Texture2D&) = delete;

  Texture2D(Texture2D&&) noexcept;
  Texture2D& operator=(Texture2D&&) noexcept;

  bool loadBMP(const std::string& path, std::string* outError = nullptr);
  void createChecker(int w, int h, int checkSize);

  // Create a texture from raw RGBA8 pixels (row-major, 4 bytes per pixel).
  // Useful for procedural textures and UI sprites.
  void createRGBA(int w, int h, const void* rgbaPixels,
                  bool generateMips = false,
                  bool nearestFilter = true,
                  bool clampToEdge = true);

  // Allocate an empty RGBA8 texture (initial contents undefined). Use updateRGBA()
  // to upload regions.
  void allocateRGBA(int w, int h,
                    bool generateMips = false,
                    bool nearestFilter = true,
                    bool clampToEdge = true);

  // Update a sub-region of an RGBA8 texture previously created via createRGBA()
  // or allocateRGBA(). (x,y) is the lower-left corner in texel coordinates.
  void updateRGBA(int x, int y, int w, int h, const void* rgbaPixels);

  void bind(int unit = 0) const;

  unsigned int handle() const { return tex_; }
  int width() const { return w_; }
  int height() const { return h_; }

private:
  void destroy();

  unsigned int tex_{0};
  int w_{0};
  int h_{0};
};

} // namespace stellar::render
