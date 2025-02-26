#include "transforms.h"

#include "CGL/matrix3x3.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"

namespace CGL {

Vector2D operator*(const Matrix3x3 &m, const Vector2D &v) {
	Vector3D mv = m * Vector3D(v.x, v.y, 1);
	return Vector2D(mv.x / mv.z, mv.y / mv.z);
}

Matrix3x3 translate(float dx, float dy) {
	// Part 3: Fill this in.


	// return Matrix3x3();

	// probelm 3
	// return Matrix3x3(1, 0, dx,
	// 	0, 1, dy,
	// 	0, 0, 1);

	return Matrix3x3(1, 1, dx,
		0, 1, dy,
		0, 0, 1);
}

Matrix3x3 scale(float sx, float sy) {
	// Part 3: Fill this in.


	// return Matrix3x3();

	// problem 3
	return Matrix3x3(sx, 0, 0,
		0, sy, 0,
		0, 0, 1);

}

// The input argument is in degrees counterclockwise
Matrix3x3 rotate(float deg) {
	// Part 3: Fill this in.


	// return Matrix3x3();

	// problem 3

	// Convert degrees to radians
	// const float PI = 3.14159265358979323846f;
    float rad = deg * 3.1414159265358979323f / 180.0f;
    float c = cos(rad);
    float s = sin(rad);
    
    // Create a rotation matrix:
    // [ cos  -sin  0 ]
    // [ sin   cos  0 ]
    // [  0     0   1 ]
    return Matrix3x3(c, -s, 0,
                     s,  c, 0,
                     0,  0, 1);


}

}
