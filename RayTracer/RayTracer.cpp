// RayTracer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <glm.hpp>
#include <gtx/transform.hpp>

#define M_PI 3.14159265359
const float kInfinity = std::numeric_limits<float>::max();


float deg2rad(const float& deg)
{
	return deg * M_PI / 180;
}
float clamp(const float& lo, const float& hi, const float& v)
{
	return std::max(lo, std::min(hi, v));
}
glm::vec3 multDirMatrix(const glm::mat4& mat, const glm::vec3& src)
{
	float a, b, c;

	a = src[0] * mat[0][0] + src[1] * mat[1][0] + src[2] * mat[2][0];
	b = src[0] * mat[0][1] + src[1] * mat[1][1] + src[2] * mat[2][1];
	c = src[0] * mat[0][2] + src[1] * mat[1][2] + src[2] * mat[2][2];

	return glm::vec3(a, b, c);
}
glm::vec3 multVecMatrix(const glm::mat4& mat, const glm::vec3& src)
{
	float a, b, c, w;

	a = src[0] * mat[0][0] + src[1] * mat[1][0] + src[2] * mat[2][0] + mat[3][0];
	b = src[0] * mat[0][1] + src[1] * mat[1][1] + src[2] * mat[2][1] + mat[3][1];
	c = src[0] * mat[0][2] + src[1] * mat[1][2] + src[2] * mat[2][2] + mat[3][2];
	w = src[0] * mat[0][3] + src[1] * mat[1][3] + src[2] * mat[2][3] + mat[3][3];
	
	return glm::vec3(a, b, c);
}
bool solveQuadratic(const float& a, const float& b, const float& c, float& x0, float& x1)
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
struct Ray {
	glm::vec3 origin;
	glm::vec3 direction;

	Ray(const glm::vec3& o, const glm::vec3& dir) :origin(o), direction(dir) {}
};
class Object {
public:
	Object() {}
	virtual ~Object() {}
};
class Light
{
public:
	Light() {}
};
struct Options {
	uint32_t width;
	uint32_t height;
	float fov;
};
struct Sphere {
	float radius;
	glm::vec3 center;
	Sphere(const glm::vec3& c, const float& r) : center(c), radius(r) {}
	bool Intersect(const Ray& ray, double& t) {

		float t0=0, t1=0;

		glm::vec3 L = ray.origin - center;
		float a = glm::dot(ray.direction, ray.direction);
		float b = 2 * glm::dot(ray.direction,L);
		float c = glm::dot(L,L) - (radius*radius);
		if (!solveQuadratic(a, b, c, t0, t1)) return false;
		if (t0 > t1) std::swap(t0, t1);

		if (t0 < 0) {
			t0 = t1; // if t0 is negative, let's use t1 instead 
			if (t0 < 0) return false; // both t0 and t1 are negative 
		}

		t = t0;

		return true;
	}
};

glm::vec3 CastRay(const glm::vec3& origin,
	const glm::vec3& direction,
	const std::vector<std::unique_ptr<Object>>& objects,
	const std::vector<std::unique_ptr<Light>>& lights,
	const Options& options,
	uint32_t depth)
{
	glm::vec3 hitColor = (direction + glm::vec3(1.0f)) * 0.5f;
	return hitColor;
}

void Render(const Options& options,
	const std::vector<std::unique_ptr<Object>>& objects,
	const std::vector<std::unique_ptr<Light>>& lights)
{
	glm::mat4 cameraToWorld(1.0f);
	glm::vec3* framebuffer = new glm::vec3[options.width * options.height];
	glm::vec3* pix = framebuffer;
	float scale = tan(deg2rad(options.fov * 0.5));
	float imageAspectRatio = options.width / (float)options.height;

	glm::vec3 origin;
	origin = multVecMatrix(cameraToWorld, glm::vec3(0.0f));
	for (int j = 0; j < options.height; ++j) {
		for (int i = 0; i < options.width; ++i) {
			float x = (2 * (i + 0.5) / (float)options.width - 1) * imageAspectRatio * scale;
			float y = (1 - 2 * (j + 0.5) / (float)options.height) * scale;

			glm::vec3 dir = multDirMatrix(cameraToWorld, glm::vec3(x, y, -1));
			dir = glm::normalize(dir);

			*(pix++) = CastRay(origin, dir, objects, lights, options, 0);
		}
	}
	// Save result to a PPM image (keep these flags if you compile under Windows)
	std::ofstream ofs("./out.ppm", std::ios::out | std::ios::binary);
	ofs << "P6\n" << options.width << " " << options.height << "\n255\n";
	for (uint32_t i = 0; i < options.height * options.width; ++i) {
		char r = (char)(255 * clamp(0, 1, framebuffer[i].x));
		char g = (char)(255 * clamp(0, 1, framebuffer[i].y));
		char b = (char)(255 * clamp(0, 1, framebuffer[i].z));
		ofs << r << g << b;
	}

	ofs.close();

	delete[] framebuffer;
}

int main()
{
	// creating the scene (adding objects and lights)
	std::vector<std::unique_ptr<Object>> objects;
	std::vector<std::unique_ptr<Light>> lights;

	// setting up options
	Options options;
	options.width = 640;
	options.height = 480;
	options.fov = 90;

	// finally, render
	Render(options, objects, lights);

	return 0;


	//const int WIDTH = 500;
	//const int HEIGHT = 500;
	//Sphere sphere(glm::vec3(WIDTH * 0.5, HEIGHT * 0.5, 50), 100);
	//glm::vec3 origin(0.0f, 0.0f, 0.0f);
	//std::ofstream out("out.ppm");
	//out << "P3\n" << WIDTH << ' ' << HEIGHT << ' ' << "255\n";

	//double t;

	//for (int y = 0; y < HEIGHT; ++y) {
	//	for (int x = 0; x < WIDTH; ++x) {
	//		
	//		glm::vec3 dir = glm::vec3(0.0f,0.0f, 1.0f);
	//		origin = glm::vec3(x, y, 0.0f);
	//		const Ray ray(origin, dir);

	//		t = 200000;
	//		glm::vec3 pixel_col = glm::vec3(0.0f,0.0f,0.0f);
	//		if (sphere.Intersect(ray, t))
	//		{
	//			pixel_col.r = 255;
	//			pixel_col.g = 255;
	//			pixel_col.b = 255;
	//		}
	//		out << (int)pixel_col.x << ' '
	//			<< (int)pixel_col.y << ' '
	//			<< (int)pixel_col.z << '\n';
	//	}
	//}
	//return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
