#pragma once

#include <string>
#include <string_view>

namespace stellar::core {

// RFC 4648 compatible Base64 helpers.
//
// Intended use in StellarForge:
//  - Encode arbitrary byte strings (including UTF-8) into a compact, whitespace-free token
//    that can be embedded safely into the line-oriented SaveGame format.
//
// NOTE:
//  - This is standard ("+" and "/" alphabet with "=" padding).
//  - The decoder is tolerant of ASCII whitespace.

std::string base64Encode(std::string_view bytes);

// Decodes Base64 into `outBytes`.
// Returns false if the input contains invalid characters or malformed padding.
bool base64Decode(std::string_view b64, std::string* outBytes);

} // namespace stellar::core
