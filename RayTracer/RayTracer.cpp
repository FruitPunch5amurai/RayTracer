// RayTracer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <cstdio> 
#include <cstdlib> 
#include <memory> 
#include <vector> 
#include <utility> 
#include <cstdint> 
#include <iostream> 
#include <fstream> 
#include <cmath> 
#include <limits> 
#include <random> 

#include <glm.hpp>
#include <gtx/transform.hpp>

#define M_PI 3.14159265359
const float kInfinity = std::numeric_limits<float>::max();
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0, 1);

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
inline glm::vec3 mix(const glm::vec3& a, const glm::vec3& b, const float& mixValue)
{
	return a * (1 - mixValue) + b * mixValue;
}


struct Ray {
	glm::vec3 origin;
	glm::vec3 direction;

	Ray(const glm::vec3& o, const glm::vec3& dir) :origin(o), direction(dir) {}
};
class Object {
public:
	Object() : color(dis(gen), dis(gen), dis(gen)) {}
	virtual ~Object() {}

	virtual bool Intersect(const glm::vec3& orig, const glm::vec3& dir, float& t) const = 0;
	virtual void GetSurfaceData(const glm::vec3&, glm::vec3&, glm::vec2&) const = 0;
	glm::vec3 color;
};
class Sphere : public Object {
public:
	glm::vec3 center;
	float radius;
	float radius2;
	Sphere(const glm::vec3& origin, const float& r) :center(origin), radius(r), radius2(r* r) {}
	bool Intersect(const glm::vec3& orig, const glm::vec3& dir, float& t) const {
		float t0, t1; // solutions for t if the ray intersects 

					  // analytic solution
		glm::vec3 L = orig - center;
		float a = glm::dot(dir,dir);
		float b = 2 * glm::dot(dir,L);
		float c = glm::dot(L,L) - radius2;
		if (!solveQuadratic(a, b, c, t0, t1)) return false;

		if (t0 > t1) std::swap(t0, t1);

		if (t0 < 0) {
			t0 = t1; // if t0 is negative, let's use t1 instead 
			if (t0 < 0) return false; // both t0 and t1 are negative 
		}

		t = t0;

		return true;
	}
	void GetSurfaceData(const glm::vec3& Phit, glm::vec3& Nhit, glm::vec2& tex) const
	{
		Nhit = Phit - center;
		Nhit = glm::normalize(Nhit);
		// In this particular case, the normal is simular to a point on a unit sphere
		// centred around the origin. We can thus use the normal coordinates to compute
		// the spherical coordinates of Phit.
		// atan2 returns a value in the range [-pi, pi] and we need to remap it to range [0, 1]
		// acosf returns a value in the range [0, pi] and we also need to remap it to the range [0, 1]
		tex.x = (1 + atan2(Nhit.z, Nhit.x) / M_PI) * 0.5;
		tex.y = acosf(Nhit.y) / M_PI;
	}
};
class Light
{
public:
	Light() {}
};
struct Options {
	int width;
	int height;
	float fov;
	glm::mat4 cameraToWorld;
};
bool Trace(const glm::vec3& orig, const glm::vec3& dir, const std::vector<std::unique_ptr<Object>>& objects, float& tNear, const Object*& hitObject)
{
	tNear = kInfinity;
	std::vector<std::unique_ptr<Object>>::const_iterator iter = objects.begin();
	for (; iter != objects.end(); ++iter) {
		float t = kInfinity;
		if ((*iter)->Intersect(orig, dir, t) && t < tNear) {
			hitObject = iter->get();
			tNear = t;
		}
	}

	return (hitObject != nullptr);
}
glm::vec3 CastRay(const glm::vec3& origin,
	const glm::vec3& direction,
	const std::vector<std::unique_ptr<Object>>& objects,
	const std::vector<std::unique_ptr<Light>>& lights,
	const Options& options,
	uint32_t depth)
{
	glm::vec3 hitColor(0.0f);
	const Object* hitObject = nullptr; // this is a pointer to the hit object 
	float t; // this is the intersection distance from the ray origin to the hit point 
	if (Trace(origin, direction, objects, t, hitObject)) {
		glm::vec3 Phit = origin + direction * t;
		glm::vec3 Nhit;
		glm::vec2 tex;
		hitObject->GetSurfaceData(Phit, Nhit, tex);
		// Use the normal and texture coordinates to shade the hit point.
		// The normal is used to compute a simple facing ratio and the texture coordinate
		// to compute a basic checker board pattern
		float scale = 4;
		float pattern = (fmodf(tex.x * scale, 1) > 0.5) ^ (fmodf(tex.y * scale, 1) > 0.5);
		hitColor = std::max(0.f, glm::dot(Nhit,-direction)) 
			* mix(hitObject->color, hitObject->color * 0.8f, pattern);
	}

	return hitColor;
}

void Render(const Options& options,const std::vector<std::unique_ptr<Object>>& objects,const std::vector<std::unique_ptr<Light>>& lights)
{	
	glm::vec3* framebuffer = new glm::vec3[options.width * options.height];
	glm::vec3* pix = framebuffer;
	float scale = tan(deg2rad(options.fov * 0.5));
	float imageAspectRatio = options.width / (float)options.height;

	glm::vec3 origin;
	origin = multVecMatrix(options.cameraToWorld, glm::vec3(0.0f));
	for (int j = 0; j < options.height; ++j) {
		for (int i = 0; i < options.width; ++i) {
			//Normalize pixel Position using the frame's dimensions.. Basically we are converting
			//from Raster space to NDC space to Screen space to World space to Camera space.
			float x = (2 * (i + 0.5) / (float)options.width - 1) * imageAspectRatio * scale;
			float y = (1 - 2 * (j + 0.5) / (float)options.height) * scale;

			//float x = i - options.width / 2.0f;
			//float y = options.height / 2 - j;
			//float z = -(options.height / 2) / scale;

			//Use a Matrix to represent the camera transform...
			//Translate the new direction vector by the camera matrix
			glm::vec3 dir = multDirMatrix(options.cameraToWorld, glm::vec3(x, y, -1));
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

	// generate a scene made of random spheres
	uint32_t numSpheres = 32;
	gen.seed(0);
	for (uint32_t i = 0; i < numSpheres; ++i) {
		glm::vec3 randPos((0.5 - dis(gen)) * 10, (0.5 - dis(gen)) * 10, (0.5 + dis(gen) * 10));
		float randRadius = (0.5 + dis(gen) * 0.5);
		objects.push_back(std::unique_ptr<Object>(new Sphere(randPos, randRadius)));
	}

	// setting up options
	Options options;
	options.width = 640;
	options.height = 480;
	options.fov = 90;
	options.cameraToWorld = glm::mat4(0.945519, 0, -0.325569, 0, -0.179534, 0.834209, -0.521403, 0, 0.271593, 0.551447, 0.78876, 0, 4.208271, 8.374532, 17.932925, 1);

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
