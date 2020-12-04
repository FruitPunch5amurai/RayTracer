#pragma once
#include <glm.hpp>
#include <utility> 

inline float clamp(const float& lo, const float& hi, const float& v)
{
	return std::max(lo, std::min(hi, v));
}
inline float deg2rad(const float& deg)
{
	return deg * M_PI / 180;
}
inline glm::vec3 multDirMatrix(const glm::mat4& mat, const glm::vec3& src)
{
	float a, b, c;

	a = src[0] * mat[0][0] + src[1] * mat[1][0] + src[2] * mat[2][0];
	b = src[0] * mat[0][1] + src[1] * mat[1][1] + src[2] * mat[2][1];
	c = src[0] * mat[0][2] + src[1] * mat[1][2] + src[2] * mat[2][2];

	return glm::vec3(a, b, c);
}
inline glm::vec3 multVecMatrix(const glm::mat4& mat, const glm::vec3& src)
{
	float a, b, c, w;

	a = src[0] * mat[0][0] + src[1] * mat[1][0] + src[2] * mat[2][0] + mat[3][0];
	b = src[0] * mat[0][1] + src[1] * mat[1][1] + src[2] * mat[2][1] + mat[3][1];
	c = src[0] * mat[0][2] + src[1] * mat[1][2] + src[2] * mat[2][2] + mat[3][2];
	w = src[0] * mat[0][3] + src[1] * mat[1][3] + src[2] * mat[2][3] + mat[3][3];

	return glm::vec3(a, b, c);
}
inline bool solveQuadratic(const float& a, const float& b, const float& c, float& x0, float& x1)
{
	float discr = b * b - 4 * a * c;
	if (discr < 0) return false;
	else if (discr == 0) {
		x0 = x1 = -0.5 * b / a;
	}
	else {
		float q = (b > 0) ?
			-0.5 * (b + sqrt(discr)) :
			-0.5 * (b - sqrt(discr));
		x0 = q / a;
		x1 = c / q;
	}

	return true;
}
inline glm::vec3 mix(const glm::vec3& a, const glm::vec3& b, const float& mixValue)
{
	return a * (1 - mixValue) + b * mixValue;
}

