#include "texture.h"
#include "CGL/color.h"

#include <cmath>
#include <algorithm>

namespace CGL {

  Color Texture::sample(const SampleParams& sp) {
    // TODO: Task 6: Fill this in.

    // problem 6
    // Level sampling: choose from L_ZERO, L_NEAREST, or L_LINEAR.
  if (sp.lsm == L_ZERO) {
    // Always sample from the base level.
    if (sp.psm == P_NEAREST)
      return sample_nearest(sp.p_uv, 0);
    else
      return sample_bilinear(sp.p_uv, 0);
  }
  else if (sp.lsm == L_NEAREST) {
    float level = get_level(sp);
    int level_nearest = static_cast<int>(level + 0.5f);
    if (sp.psm == P_NEAREST)
      return sample_nearest(sp.p_uv, level_nearest);
    else
      return sample_bilinear(sp.p_uv, level_nearest);
  }
  else if (sp.lsm == L_LINEAR) {
    float level = get_level(sp);
    int level_low = static_cast<int>(floor(level));
    int level_high = static_cast<int>(ceil(level));
    float t = level - level_low;
    Color c_low, c_high;
    if (sp.psm == P_NEAREST) {
      c_low  = sample_nearest(sp.p_uv, level_low);
      c_high = sample_nearest(sp.p_uv, level_high);
    } else {
      c_low  = sample_bilinear(sp.p_uv, level_low);
      c_high = sample_bilinear(sp.p_uv, level_high);
    }
    return c_low * (1 - t) + c_high * t;
  }
  // Fallback: return magenta to indicate an error.

// return magenta for invalid level
    return Color(1, 0, 1);
  }

  float Texture::get_level(const SampleParams& sp) {
    // TODO: Task 6: Fill this in.

    // probelm 6
    // Use the full-resolution texture dimensions from mipmap[0].
  float w = static_cast<float>(mipmap[0].width);
  float h = static_cast<float>(mipmap[0].height);

  // Compute derivative vectors in uv space.
  Vector2D duv_dx = sp.p_dx_uv - sp.p_uv;
  Vector2D duv_dy = sp.p_dy_uv - sp.p_uv;

  // Scale the derivatives by the texture resolution.
  duv_dx[0] *= w;  duv_dx[1] *= h;
  duv_dy[0] *= w;  duv_dy[1] *= h;

  // Compute the lengths of the derivative vectors.
  float Lx = sqrt(duv_dx[0] * duv_dx[0] + duv_dx[1] * duv_dx[1]);
  float Ly = sqrt(duv_dy[0] * duv_dy[0] + duv_dy[1] * duv_dy[1]);

  float L = std::max(Lx, Ly);

  // Compute the mipmap level as log2(L).
  float level = log2(L);

  // Clamp level to valid range: [0, num_mip_levels-1].
  level = std::min(level, float(mipmap.size() - 1));
  level = std::max(level, 0.0f);

  return level;
  
    // return 0;
  }

  Color MipLevel::get_texel(int tx, int ty) {
    return Color(&texels[tx * 3 + ty * width * 3]);
  }

  Color Texture::sample_nearest(Vector2D uv, int level) {
    // TODO: Task 5: Fill this in.
    auto& mip = mipmap[level];

    // probelm 5
    // Clamp uv coordinates to [0, 1]
    uv[0] = std::min<float>(1.0f, std::max<float>(0.0f, uv[0]));
    uv[1] = std::min<float>(1.0f, std::max<float>(0.0f, uv[1]));


    // Map uv to texture space (centered on texel centers).
    float x = uv[0] * mip.width - 0.5f;
    float y = uv[1] * mip.height - 0.5f;
    
    int ix = static_cast<int>(std::round(x));
    int iy = static_cast<int>(std::round(y));

    // Clamp to valid indices.
    ix = std::max(0, std::min(ix, static_cast<int>(mip.width) - 1));
    iy = std::max(0, std::min(iy, static_cast<int>(mip.height) - 1));
    return mip.get_texel(ix, iy);

    // return magenta for invalid level
    // return Color(1, 0, 1);
  }

  Color Texture::sample_bilinear(Vector2D uv, int level) {
    // TODO: Task 5: Fill this in.
    auto& mip = mipmap[level];

    // problem 5

    // Clamp uv coordinates to [0, 1]
    float u = std::min<float>(1.0f, std::max<float>(0.0f, uv[0]));
    float v = std::min<float>(1.0f, std::max<float>(0.0f, uv[1]));
    // Map uv to texture space.
    float x = u * mip.width - 0.5f;
    float y = v * mip.height - 0.5f;
    
    int x0 = static_cast<int>(std::floor(x));
    int y0 = static_cast<int>(std::floor(y));
    int x1 = x0 + 1;
    int y1 = y0 + 1;
    
    // Compute the fractional parts.
    float s = x - x0;
    float t = y - y0;
    
    // Clamp texture indices.
    x0 = std::max(0, std::min(x0, static_cast<int>(mip.width) - 1));
    y0 = std::max(0, std::min(y0, static_cast<int>(mip.height) - 1));
    x1 = std::max(0, std::min(x1, static_cast<int>(mip.width) - 1));
    y1 = std::max(0, std::min(y1, static_cast<int>(mip.height) - 1));


    // Fetch the four surrounding texels.
    Color c00 = mip.get_texel(x0, y0);
    Color c10 = mip.get_texel(x1, y0);
    Color c01 = mip.get_texel(x0, y1);
    Color c11 = mip.get_texel(x1, y1);
    
    // Interpolate along the x direction.
    Color c0 = c00 * (1 - s) + c10 * s;
    Color c1 = c01 * (1 - s) + c11 * s;
    
    // Interpolate along the y direction.
    Color result = c0 * (1 - t) + c1 * t;
    
    return result;


    // return magenta for invalid level
    // return Color(1, 0, 1);
  }



  /****************************************************************************/

  // Helpers

  inline void uint8_to_float(float dst[3], unsigned char* src) {
    uint8_t* src_uint8 = (uint8_t*)src;
    dst[0] = src_uint8[0] / 255.f;
    dst[1] = src_uint8[1] / 255.f;
    dst[2] = src_uint8[2] / 255.f;
  }

  inline void float_to_uint8(unsigned char* dst, float src[3]) {
    uint8_t* dst_uint8 = (uint8_t*)dst;
    dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
    dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
    dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
  }

  void Texture::generate_mips(int startLevel) {

    // make sure there's a valid texture
    if (startLevel >= mipmap.size()) {
      std::cerr << "Invalid start level";
    }

    // allocate sublevels
    int baseWidth = mipmap[startLevel].width;
    int baseHeight = mipmap[startLevel].height;
    int numSubLevels = (int)(log2f((float)max(baseWidth, baseHeight)));

    numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
    mipmap.resize(startLevel + numSubLevels + 1);

    int width = baseWidth;
    int height = baseHeight;
    for (int i = 1; i <= numSubLevels; i++) {

      MipLevel& level = mipmap[startLevel + i];

      // handle odd size texture by rounding down
      width = max(1, width / 2);
      //assert (width > 0);
      height = max(1, height / 2);
      //assert (height > 0);

      level.width = width;
      level.height = height;
      level.texels = vector<unsigned char>(3 * width * height);
    }

    // create mips
    int subLevels = numSubLevels - (startLevel + 1);
    for (int mipLevel = startLevel + 1; mipLevel < startLevel + subLevels + 1;
      mipLevel++) {

      MipLevel& prevLevel = mipmap[mipLevel - 1];
      MipLevel& currLevel = mipmap[mipLevel];

      int prevLevelPitch = prevLevel.width * 3; // 32 bit RGB
      int currLevelPitch = currLevel.width * 3; // 32 bit RGB

      unsigned char* prevLevelMem;
      unsigned char* currLevelMem;

      currLevelMem = (unsigned char*)&currLevel.texels[0];
      prevLevelMem = (unsigned char*)&prevLevel.texels[0];

      float wDecimal, wNorm, wWeight[3];
      int wSupport;
      float hDecimal, hNorm, hWeight[3];
      int hSupport;

      float result[3];
      float input[3];

      // conditional differentiates no rounding case from round down case
      if (prevLevel.width & 1) {
        wSupport = 3;
        wDecimal = 1.0f / (float)currLevel.width;
      }
      else {
        wSupport = 2;
        wDecimal = 0.0f;
      }

      // conditional differentiates no rounding case from round down case
      if (prevLevel.height & 1) {
        hSupport = 3;
        hDecimal = 1.0f / (float)currLevel.height;
      }
      else {
        hSupport = 2;
        hDecimal = 0.0f;
      }

      wNorm = 1.0f / (2.0f + wDecimal);
      hNorm = 1.0f / (2.0f + hDecimal);

      // case 1: reduction only in horizontal size (vertical size is 1)
      if (currLevel.height == prevLevel.height) {
        //assert (currLevel.height == 1);

        for (int i = 0; i < currLevel.width; i++) {
          wWeight[0] = wNorm * (1.0f - wDecimal * i);
          wWeight[1] = wNorm * 1.0f;
          wWeight[2] = wNorm * wDecimal * (i + 1);

          result[0] = result[1] = result[2] = 0.0f;

          for (int ii = 0; ii < wSupport; ii++) {
            uint8_to_float(input, prevLevelMem + 3 * (2 * i + ii));
            result[0] += wWeight[ii] * input[0];
            result[1] += wWeight[ii] * input[1];
            result[2] += wWeight[ii] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (3 * i), result);
        }

        // case 2: reduction only in vertical size (horizontal size is 1)
      }
      else if (currLevel.width == prevLevel.width) {
        //assert (currLevel.width == 1);

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          result[0] = result[1] = result[2] = 0.0f;
          for (int jj = 0; jj < hSupport; jj++) {
            uint8_to_float(input, prevLevelMem + prevLevelPitch * (2 * j + jj));
            result[0] += hWeight[jj] * input[0];
            result[1] += hWeight[jj] * input[1];
            result[2] += hWeight[jj] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (currLevelPitch * j), result);
        }

        // case 3: reduction in both horizontal and vertical size
      }
      else {

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          for (int i = 0; i < currLevel.width; i++) {
            wWeight[0] = wNorm * (1.0f - wDecimal * i);
            wWeight[1] = wNorm * 1.0f;
            wWeight[2] = wNorm * wDecimal * (i + 1);

            result[0] = result[1] = result[2] = 0.0f;

            // convolve source image with a trapezoidal filter.
            // in the case of no rounding this is just a box filter of width 2.
            // in the general case, the support region is 3x3.
            for (int jj = 0; jj < hSupport; jj++)
              for (int ii = 0; ii < wSupport; ii++) {
                float weight = hWeight[jj] * wWeight[ii];
                uint8_to_float(input, prevLevelMem +
                  prevLevelPitch * (2 * j + jj) +
                  3 * (2 * i + ii));
                result[0] += weight * input[0];
                result[1] += weight * input[1];
                result[2] += weight * input[2];
              }

            // convert back to format of the texture
            float_to_uint8(currLevelMem + currLevelPitch * j + 3 * i, result);
          }
        }
      }
    }
  }

}
