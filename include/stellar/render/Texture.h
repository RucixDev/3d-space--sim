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
