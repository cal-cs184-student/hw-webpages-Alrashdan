#include "rasterizer.h"

using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)


    // sample_buffer[y * width + x] = c;





    //problem 2
     // For points/lines, fill all the samples for the pixel.
     for (unsigned int s = 0; s < sample_rate; s++) {
      sample_buffer[(y * width + x) * sample_rate + s] = c;
    }


  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    fill_pixel(sx, sy, color);
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
    // TODO: Task 2: Update to implement super-sampled rasterization

  //   // Compute the bounding box of the triangle
  // float min_x = std::min({ x0, x1, x2 });
  // float max_x = std::max({ x0, x1, x2 });
  // float min_y = std::min({ y0, y1, y2 });
  // float max_y = std::max({ y0, y1, y2 });

  // // Clamp the bounding box to the framebuffer dimensions
  // int x_start = std::max(0, static_cast<int>(floor(min_x)));
  // int x_end   = std::min(static_cast<int>(width) - 1, static_cast<int>(ceil(max_x)));
  // int y_start = std::max(0, static_cast<int>(floor(min_y)));
  // int y_end   = std::min(static_cast<int>(height) - 1, static_cast<int>(ceil(max_y)));

  // // Define an edge function lambda
  // // Given edge from (ax, ay) to (bx, by) and a point (cx, cy)
  // // the function returns a signed area. A positive value means (cx,cy) is on one side,
  // // a negative value means the other.
  // auto edge = [](float ax, float ay, float bx, float by, float cx, float cy) -> float {
  //   return (cx - ax) * (by - ay) - (cy - ay) * (bx - ax);
  // };

  // // Iterate over each pixel in the bounding box
  // for (int y = y_start; y <= y_end; y++) {
  //   for (int x = x_start; x <= x_end; x++) {
  //     // Compute the sample point at the center of the pixel
  //     float sample_x = x + 0.5f;
  //     float sample_y = y + 0.5f;

  //     // Evaluate the edge functions for the sample point
  //     float e0 = edge(x1, y1, x2, y2, sample_x, sample_y);
  //     float e1 = edge(x2, y2, x0, y0, sample_x, sample_y);
  //     float e2 = edge(x0, y0, x1, y1, sample_x, sample_y);

  //     // If the sample is on the correct side of all three edges (or exactly on an edge),
  //     // the point is inside the triangle.
  //     if ((e0 >= 0 && e1 >= 0 && e2 >= 0) || (e0 <= 0 && e1 <= 0 && e2 <= 0)) {
  //       fill_pixel(x, y, color);
  //     }
  //   }
  // }
      










  // problem 2
  // Compute the triangle bounding box in pixel coordinates
  float min_x = std::min({ x0, x1, x2 });
  float max_x = std::max({ x0, x1, x2 });
  float min_y = std::min({ y0, y1, y2 });
  float max_y = std::max({ y0, y1, y2 });

  int x_start = std::max(0, static_cast<int>(floor(min_x)));
  int x_end   = std::min(static_cast<int>(width) - 1, static_cast<int>(ceil(max_x)));
  int y_start = std::max(0, static_cast<int>(floor(min_y)));
  int y_end   = std::min(static_cast<int>(height) - 1, static_cast<int>(ceil(max_y)));

  // Determine the number of subpixel samples per dimension.
  int n = static_cast<int>(sqrt(sample_rate));  // e.g., if sample_rate==4, then n==2

  // Lambda for the edge function.
  // This computes the signed area (or twice the area) of the triangle (a,b,c)
  auto edge = [](float ax, float ay, float bx, float by, float cx, float cy) -> float {
    return (cx - ax) * (by - ay) - (cy - ay) * (bx - ax);
  };

  // Loop over each pixel in the bounding box.
  for (int y = y_start; y <= y_end; y++) {
    for (int x = x_start; x <= x_end; x++) {
      // Loop over the n x n subpixel grid
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          // Compute the sample location within the pixel.
          // For example, for n=2, the offsets are 0.25 and 0.75.
          float sample_x = x + (i + 0.5f) / n;
          float sample_y = y + (j + 0.5f) / n;

          // Evaluate the edge functions at the sample location.
          float e0 = edge(x1, y1, x2, y2, sample_x, sample_y);
          float e1 = edge(x2, y2, x0, y0, sample_x, sample_y);
          float e2 = edge(x0, y0, x1, y1, sample_x, sample_y);

          // A point is inside the triangle if it is on the same side (or exactly on) all three edges.
          if ((e0 >= 0 && e1 >= 0 && e2 >= 0) ||
              (e0 <= 0 && e1 <= 0 && e2 <= 0)) {
            // Compute the index into the sample buffer.
            // The samples for pixel (x,y) are stored contiguously.
            int sample_index = (y * width + x) * sample_rate + (j * n + i);
            sample_buffer[sample_index] = color;
          }
        }
      }
    }
  }

}


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle

    // problem 4

    // Compute the bounding box of the triangle in pixel coordinates.
  float min_x = std::min({ x0, x1, x2 });
  float max_x = std::max({ x0, x1, x2 });
  float min_y = std::min({ y0, y1, y2 });
  float max_y = std::max({ y0, y1, y2 });

  int x_start = std::max(0, static_cast<int>(floor(min_x)));
  int x_end   = std::min(static_cast<int>(width) - 1, static_cast<int>(ceil(max_x)));
  int y_start = std::max(0, static_cast<int>(floor(min_y)));
  int y_end   = std::min(static_cast<int>(height) - 1, static_cast<int>(ceil(max_y)));

  // Determine the number of subpixel samples per dimension.
  int n = static_cast<int>(sqrt(sample_rate));  // e.g., sample_rate==4 implies n==2

  // Lambda for the edge function.
  // Returns a signed value proportional to the area of the triangle (ax,ay)-(bx,by)-(cx,cy)
  auto edge = [](float ax, float ay, float bx, float by, float cx, float cy) -> float {
    return (cx - ax) * (by - ay) - (cy - ay) * (bx - ax);
  };

  // Compute the total area of the triangle (twice the area actually).
  float area = edge(x0, y0, x1, y1, x2, y2);
  if (area == 0) return;  // Degenerate triangle

  // Loop over each pixel in the bounding box.
  for (int y = y_start; y <= y_end; y++) {
    for (int x = x_start; x <= x_end; x++) {
      // Loop over the subpixel grid.
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          // Compute the center of this subpixel.
          float sample_x = x + (i + 0.5f) / n;
          float sample_y = y + (j + 0.5f) / n;

          // Evaluate edge functions for this sample point.
          float e0 = edge(x1, y1, x2, y2, sample_x, sample_y);
          float e1 = edge(x2, y2, x0, y0, sample_x, sample_y);
          float e2 = edge(x0, y0, x1, y1, sample_x, sample_y);

          // Check if the sample is inside the triangle.
          // Accept the sample if it lies on the correct side of all edges (or exactly on an edge).
          if ((e0 >= 0 && e1 >= 0 && e2 >= 0) ||
              (e0 <= 0 && e1 <= 0 && e2 <= 0)) {

            // Compute barycentric coordinates.
            // Note: These are proportional to the areas of the subtriangles.
            float alpha = edge(x1, y1, x2, y2, sample_x, sample_y) / area;
            float beta  = edge(x2, y2, x0, y0, sample_x, sample_y) / area;
            float gamma = edge(x0, y0, x1, y1, sample_x, sample_y) / area;

            // Interpolate the color using the barycentric coordinates.
            Color sampleColor = c0 * alpha + c1 * beta + c2 * gamma;

            // Compute the sample index in the sample buffer and store the computed color.
            int sample_index = (y * width + x) * sample_rate + (j * n + i);
            sample_buffer[sample_index] = sampleColor;
          }
        }
      }
    }
  }
}




  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle

//     //probelm 5
 // Compute the bounding box in screen space.
 float min_x = std::min({ x0, x1, x2 });
 float max_x = std::max({ x0, x1, x2 });
 float min_y = std::min({ y0, y1, y2 });
 float max_y = std::max({ y0, y1, y2 });

 int x_start = std::max(0, static_cast<int>(floor(min_x)));
 int x_end   = std::min(static_cast<int>(width) - 1, static_cast<int>(ceil(max_x)));
 int y_start = std::max(0, static_cast<int>(floor(min_y)));
 int y_end   = std::min(static_cast<int>(height) - 1, static_cast<int>(ceil(max_y)));

 // Determine the number of subpixel samples per dimension.
 int n = static_cast<int>(sqrt(sample_rate));

 // Lambda for the edge function.
 auto edge = [](float ax, float ay, float bx, float by, float cx, float cy) -> float {
   return (cx - ax) * (by - ay) - (cy - ay) * (bx - ax);
 };

 // Compute the total area of the triangle (used for barycentrics).
 float area = edge(x0, y0, x1, y1, x2, y2);
 if (area == 0) return; // degenerate triangle

 // Loop over each pixel in the bounding box.
 for (int y = y_start; y <= y_end; y++) {
   for (int x = x_start; x <= x_end; x++) {
     // Loop over each subpixel.
     for (int i = 0; i < n; i++) {
       for (int j = 0; j < n; j++) {
         // Compute the center of the subpixel.
         float sample_x = x + (i + 0.5f) / n;
         float sample_y = y + (j + 0.5f) / n;

         // Evaluate edge functions for point-in-triangle test.
         float e0 = edge(x1, y1, x2, y2, sample_x, sample_y);
         float e1 = edge(x2, y2, x0, y0, sample_x, sample_y);
         float e2 = edge(x0, y0, x1, y1, sample_x, sample_y);

         // Accept the sample if it lies inside the triangle (or exactly on an edge).
         if ((e0 >= 0 && e1 >= 0 && e2 >= 0) ||
             (e0 <= 0 && e1 <= 0 && e2 <= 0)) {

           // Compute barycentric coordinates.
           float alpha = edge(x1, y1, x2, y2, sample_x, sample_y) / area;
           float beta  = edge(x2, y2, x0, y0, sample_x, sample_y) / area;
           float gamma = edge(x0, y0, x1, y1, sample_x, sample_y) / area;

           // Interpolate UV coordinates.
           float u = alpha * u0 + beta * u1 + gamma * u2;
           float v = alpha * v0 + beta * v1 + gamma * v2;
           Vector2D sample_uv(u, v);

           // Use the texture to sample a color.
           Color sampleColor;
           if (psm == P_NEAREST) {
             sampleColor = tex.sample_nearest(sample_uv, 0);
           } else if ( psm == P_LINEAR) {
             sampleColor = tex.sample_bilinear(sample_uv, 0);
           }

           // Store the sample color into the supersample buffer.
           int sample_index = (y * width + x) * sample_rate + (j * n + i);
           sample_buffer[sample_index] = sampleColor;
         }
       }
     }
   }
 }


 //////////// uncomment below and comment above to test task 6

// // probelm 6
// // Compute the bounding box in screen space.
// float min_x = std::min({ x0, x1, x2 });
// float max_x = std::max({ x0, x1, x2 });
// float min_y = std::min({ y0, y1, y2 });
// float max_y = std::max({ y0, y1, y2 });

// int x_start = std::max(0, static_cast<int>(floor(min_x)));
// int x_end   = std::min(static_cast<int>(width) - 1, static_cast<int>(ceil(max_x)));
// int y_start = std::max(0, static_cast<int>(floor(min_y)));
// int y_end   = std::min(static_cast<int>(height) - 1, static_cast<int>(ceil(max_y)));

// // Determine subpixel grid dimensions (assumes sample_rate is a perfect square).
// int n = static_cast<int>(sqrt(sample_rate));

// // Lambda for computing the edge function (twice the signed area).
// auto edge = [](float ax, float ay, float bx, float by, float cx, float cy) -> float {
//   return (cx - ax) * (by - ay) - (cy - ay) * (bx - ax);
// };

// // Compute the triangle’s total area (used for barycentrics).
// float area = edge(x0, y0, x1, y1, x2, y2);
// if (area == 0) return; // Degenerate triangle

// // Loop over each pixel in the bounding box.
// for (int y = y_start; y <= y_end; y++) {
//   for (int x = x_start; x <= x_end; x++) {
//     // Loop over the n x n subpixel grid.
//     for (int i = 0; i < n; i++) {
//       for (int j = 0; j < n; j++) {
//         // Compute the center of this subpixel.
//         float sample_x = x + (i + 0.5f) / n;
//         float sample_y = y + (j + 0.5f) / n;

//         // Evaluate edge functions at (sample_x, sample_y).
//         float e0 = edge(x1, y1, x2, y2, sample_x, sample_y);
//         float e1 = edge(x2, y2, x0, y0, sample_x, sample_y);
//         float e2 = edge(x0, y0, x1, y1, sample_x, sample_y);

//         // Check if the sample is inside the triangle (or exactly on an edge).
//         if ((e0 >= 0 && e1 >= 0 && e2 >= 0) ||
//             (e0 <= 0 && e1 <= 0 && e2 <= 0)) {

//           // Compute barycentrics for (sample_x, sample_y).
//           float alpha = e0 / area;
//           float beta  = e1 / area;
//           float gamma = e2 / area;

//           // Interpolate uv coordinates at (sample_x, sample_y).
//           Vector2D p_uv = alpha * Vector2D(u0, v0) +
//                           beta  * Vector2D(u1, v1) +
//                           gamma * Vector2D(u2, v2);

//           // Compute barycentrics for (sample_x+1, sample_y).
//           float e0_dx = edge(x1, y1, x2, y2, sample_x + 1, sample_y);
//           float e1_dx = edge(x2, y2, x0, y0, sample_x + 1, sample_y);
//           float e2_dx = edge(x0, y0, x1, y1, sample_x + 1, sample_y);
//           float alpha_dx = e0_dx / area;
//           float beta_dx  = e1_dx / area;
//           float gamma_dx = e2_dx / area;
//           Vector2D p_dx_uv = alpha_dx * Vector2D(u0, v0) +
//                              beta_dx  * Vector2D(u1, v1) +
//                              gamma_dx * Vector2D(u2, v2);

//           // Compute barycentrics for (sample_x, sample_y+1).
//           float e0_dy = edge(x1, y1, x2, y2, sample_x, sample_y + 1);
//           float e1_dy = edge(x2, y2, x0, y0, sample_x, sample_y + 1);
//           float e2_dy = edge(x0, y0, x1, y1, sample_x, sample_y + 1);
//           float alpha_dy = e0_dy / area;
//           float beta_dy  = e1_dy / area;
//           float gamma_dy = e2_dy / area;
//           Vector2D p_dy_uv = alpha_dy * Vector2D(u0, v0) +
//                              beta_dy  * Vector2D(u1, v1) +
//                              gamma_dy * Vector2D(u2, v2);

//           // Set up sample parameters.
//           SampleParams sp;
//           sp.p_uv    = p_uv;
//           sp.p_dx_uv = p_dx_uv;
//           sp.p_dy_uv = p_dy_uv;
//           sp.lsm     = this->lsm;  // Level-sampling method from RasterizerImp.
//           sp.psm     = this->psm;  // Pixel-sampling method from RasterizerImp.

//           // Use texture sampling with mipmap level selection.
//           Color sampleColor = tex.sample(sp);

//           // Write the color into the supersample buffer.
//           int sample_index = (y * width + x) * sample_rate + (j * n + i);
//           sample_buffer[sample_index] = sampleColor;
//         }
//       }
//     }
//   }
// }
}

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    // this->sample_rate = rate;

    // this->sample_buffer.resize(width * height, Color::White);



    // problem 2

    this->sample_rate = rate;
    sample_buffer.resize(width * height * sample_rate, Color::White);
  }
  


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // // TODO: Task 2: You may want to update this function for supersampling support

    // this->width = width;
    // this->height = height;
    // this->rgb_framebuffer_target = rgb_framebuffer;


    // this->sample_buffer.resize(width * height, Color::White);




    // problem 2
    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;
    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Clear both the framebuffer and the supersample buffer.
  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }



  // void RasterizerImp::clear_buffers() {
  //   std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
  //   std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  // }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support


    // for (int x = 0; x < width; ++x) {
    //   for (int y = 0; y < height; ++y) {
    //     Color col = sample_buffer[y * width + x];

    //     for (int k = 0; k < 3; ++k) {
    //       this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
    //     }
    //   }
    // }




    // problem 2
    // // Loop over each pixel.
    // for (int y = 0; y < height; ++y) {
    //   for (int x = 0; x < width; ++x) {
    //     Color avg(0, 0, 0);
    //     int base_index = (y * width + x) * sample_rate;
    //     // Sum over all supersamples for this pixel.
    //     for (unsigned int s = 0; s < sample_rate; s++) {
    //       avg = avg + sample_buffer[base_index + s];
    //     }
    //     // Average the color.
    //     avg = avg * (1.0f / sample_rate);

    //     // Write the averaged color to the framebuffer.
    //     int fb_index = 3 * (y * width + x);
    //     rgb_framebuffer_target[fb_index    ] = static_cast<unsigned char>(std::min(1.0f, avg.r) * 255);
    //     rgb_framebuffer_target[fb_index + 1] = static_cast<unsigned char>(std::min(1.0f, avg.g) * 255);
    //     rgb_framebuffer_target[fb_index + 2] = static_cast<unsigned char>(std::min(1.0f, avg.b) * 255);
    //   }
    // }
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        Color avg(0, 0, 0);
        int base_index = (y * width + x) * sample_rate;
        // Sum over all supersamples for this pixel.
        for (unsigned int s = 0; s < sample_rate; s++) {
          avg = avg + sample_buffer[base_index + s];
        }
        // Average the color by multiplying with the reciprocal instead of using division.
        avg = avg * (1.0f / sample_rate);
  
        // Write the averaged color to the framebuffer.
        int fb_index = 3 * (y * width + x);
        rgb_framebuffer_target[fb_index    ] = static_cast<unsigned char>(std::min(1.0f, avg.r) * 255);
        rgb_framebuffer_target[fb_index + 1] = static_cast<unsigned char>(std::min(1.0f, avg.g) * 255);
        rgb_framebuffer_target[fb_index + 2] = static_cast<unsigned char>(std::min(1.0f, avg.b) * 255);
      }
    }


  }

  Rasterizer::~Rasterizer() { }


}// CGL
