// RayTracer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#define SDL_MAIN_HANDLED
#include "MainWindow.h"
#include "common.h"

#include <cstdio> 
#include <cstdlib> 
#include <memory> 
#include <vector> 
#include <utility> 
#include <cstdint> 
#include <iostream> 
#include <fstream> 
#include <cmath> 
#include <random> 

#include <glm.hpp>
#include <gtx/transform.hpp>


#define M_PI 3.14159265359
const float kInfinity = std::numeric_limits<float>::max();
constexpr float kEpsilon = 1e-8;
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0, 1);
glm::vec3* framebuffer = nullptr;


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
class Plane : public Object {
public:

	glm::vec3 origin;
	glm::vec3 normal;
	Plane(const glm::vec3& o, const glm::vec3& n): origin(o), normal(n){}
	bool Intersect(const glm::vec3& orig, const glm::vec3& dir, float& t) const {

		float denominator = glm::dot(normal, dir);
		if (glm::abs(denominator) > 0.0001f)
		{
			glm::vec3 diff = origin - orig;
			t = glm::dot(diff , normal) / denominator;
			if (t > 0.0001)
				return true;
		}
		return false;
	}
	void GetSurfaceData(const glm::vec3& Phit, glm::vec3& Nhit, glm::vec2& tex) const
	{

	}
};
class Disk : public Object {
public:
	glm::vec3 origin;
	float radius;
	Disk(const glm::vec3& o, const float& r) : origin(o), radius(r) {}
	bool Intersect(const glm::vec3& orig, const glm::vec3& dir, float& t) const {
		return false;
	}
	void GetSurfaceData(const glm::vec3& Phit, glm::vec3& Nhit, glm::vec2& tex) const
	{

	}

};
class Triangle :public Object {
public:
	glm::vec3 v0;
	glm::vec3 v1;
	glm::vec3 v2;
	bool isSingleSided = false;
	Triangle(const glm::vec3& p0, const glm::vec3& p1, const glm::vec3& p2) :
		v0(p0), v1(p1), v2(p2) {}
	bool Intersect(const glm::vec3& orig, const glm::vec3& dir, float& t) const {


		/*
		*	===Fundamentals of Computer Graphic's Way===		
		*/
		float j = v0.x - orig.x;
		float k = v0.y - orig.y;
		float l = v0.z - orig.z;
		glm::vec3 abc = glm::vec3(v0.x - v1.x, v0.y - v1.y, v0.z - v1.z);
		glm::vec3 def = glm::vec3(v0.x - v2.x, v0.y - v2.y, v0.z - v2.z);
		glm::vec3 ghi = glm::vec3(dir.x,dir.y,dir.z);

		float eiMhf = (def.y * ghi.z) - (ghi.y * def.z);
		float gfMdi = (ghi.x * def.z) - (def.x * ghi.z);
		float dhMeg = (def.x * ghi.y) - (def.y * ghi.x);
		float jcMal = (j * abc.z) - (abc.x * l);
		float akMjb = (abc.x * k) - (j * abc.y);
		float blMkc = (abc.y * l) - (k * abc.z);

			float M =
			abc.x * (eiMhf) +
			abc.y * (gfMdi) +
			abc.z * (dhMeg);

		t = (def.z*akMjb + def.y*jcMal + def.x*blMkc)/M;
		float r = (ghi.z * (akMjb)+ghi.y * (jcMal)+ghi.x * (blMkc))/M;
		if (r < 0 || r > 1)
			return false;
		float beta = (j * (eiMhf)+k * (gfMdi)+l * (dhMeg)) / M;
		if (beta < 0 || beta > 1 - r)
			return false;
		return true;
		
		/*OLD WAY*/
//		glm::vec3 newOrig = glm::vec3(orig.x, orig.y, -orig.z);

//		//Compute Plane's Normal
//		glm::vec3 v0v1 = v1 - v0;
//		glm::vec3 v0v2 = v2 - v0;
//
//		glm::vec3 N = glm::cross(v0v1, v0v2);
//		N = glm::normalize(N);
//		float area = glm::length(N); // area of the triangle 
//		/*
//		* 		Double sided test
//		 If this dot product is lower than 0, it means that the two vectors are pointing in opposite directions. 
//		 Thus, the surface is front-facing. If the dot product is greater than 0 the vectors are pointing in the same direction. 
//		 We are looking at the inside of the surface (or the back of it, it's back-facing).
//		 If the primitive was declared single-sided, 
//		 it shouldn't then be visible. In this case, the function returns false.
//		*/
//		if (glm::dot(dir, N) > 0 && isSingleSided)
//			return false; // Back-facing surface
//
//		//Step 1: Find P
//		//Check if Ray and Plane are parallel
//		float NdotRayDirection = glm::dot(N, dir);
//		if (glm::abs(NdotRayDirection) < kEpsilon)//Almost 0
//			return false;//They are paralel so they dont intersect
//
//		//compute d parameter using equation 2
//		float d = glm::dot(N, v0);
//
//		//compute t (equation 3)
//		t = (glm::dot(N, newOrig) + d) / NdotRayDirection;
//
//		//Check if triangle is in behind the ray
//		if (t < 0) 
//			return false; // the triangle is behind
//
//		//Compute the intersection point using Equation 1
//		glm::vec3 P = orig + t * dir;
//
//		//Step2: inside-outside test
//		glm::vec3 C; // vector perpendicular to triangle's plane
//
//		//edge 0
//		glm::vec3 edge0 = v1 - v0;
//		glm::vec3 vp0 = P - v0;
//		C = glm::cross(edge0, vp0);
//		if (glm::dot(N, C) < 0.0f) 
//			return false; // P is on the right side
//
//		//edge 1
//		glm::vec3 edge1 = v2 - v1;
//		glm::vec3 vp1 = P - v1;
//		C = glm::cross(edge1, vp1);
//		float u = C.length() / area;	//Baycentric Coordinates
//		if (glm::dot(N, C) < 0.0f)
//			return false; // P is on the right side
//
//		// edge 2
//		glm::vec3 edge2 = v0 - v2;
//		glm::vec3 vp2 = P - v2;
//		C = glm::cross(edge2,vp2);
//		float v = C.length() / area;	
//		if (glm::dot(N,C) < 0) 
//			return false; // P is on the right side;
//
//
//		return true;
	}
	void GetSurfaceData(const glm::vec3& Phit, glm::vec3& Nhit, glm::vec2& tex) const
	{
		glm::vec3 center;
		//Find center of triangle
		center = (v0 + v1 + v2) / 3.0f;

		Nhit = Phit - center;
		Nhit = glm::normalize(Nhit);

		/*
		*	=================Baycentric Coordinates===============
			P = uA +vB + wC			A,B,C are verts of a triangle.. u,v,w are Barycentric Coordinates
			A point is within the triangle(A,B,C) if 
			0 <= u,v,w <= 1
		*/
		//Compute Plane's Normal
		glm::vec3 v0v1 = v1 - v0;
		glm::vec3 v0v2 = v2 - v0;

		glm::vec3 N = glm::cross(v0v1, v0v2);

		float areaABC = glm::dot(N, glm::cross((v1 - v0), (v2 - v0)));
		float areaPBC = glm::dot(N, glm::cross((v1 - Phit), (v2 - Phit)));
		float areaPCA = glm::dot(N, glm::cross((v2 - Phit), (v0 - Phit)));

		float x = areaPBC / areaABC; // alpha
		float y = areaPCA / areaABC; // beta
		float z = 1.0f - x - y; // gamma

		tex.x = x * 0 + y * 1 + z * .5;
		tex.y = x * 1 + y * 1 + z * 0;

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

		//return glm::vec3(255.0f, 0.0f, 0.0f);

		hitObject->GetSurfaceData(Phit, Nhit, tex);
		// Use the normal and texture coordinates to shade the hit point.
		// The normal is used to compute a simple facing ratio and the texture coordinate
		// to compute a basic checker board pattern
		float scale = 4;
		float pattern = (fmodf(tex.x * scale, 1) > 0.5) ^ (fmodf(tex.y * scale, 1) > 0.5);
		hitColor = mix(hitObject->color, hitObject->color * 0.8f, pattern);
	}

	return hitColor;
}

void Render(const Options& options, const std::vector<std::unique_ptr<Object>>& objects, const std::vector<std::unique_ptr<Light>>& lights)
{
	framebuffer = new glm::vec3[options.width * options.height];
	glm::vec3* pix = framebuffer;
	float scale = tan(deg2rad(options.fov * 0.5));
	float imageAspectRatio = options.width / (float)options.height;

	glm::vec3 origin(0.0f);
	origin = multVecMatrix(options.cameraToWorld, glm::vec3(0.0f));
	for (int j = 0; j < options.height; ++j) {
		for (int i = 0; i < options.width; ++i) {
			//Normalize pixel Position using the frame's dimensions.. Basically we are converting
			//from Raster space to NDC space to Screen space to World space to Camera space.
			float x = (2 * (i + 0.5) / (float)options.width - 1) * imageAspectRatio * scale;
			float y = (1 - 2 * (j + 0.5) / (float)options.height) * scale;

			//Use a Matrix to represent the camera transform...
			//Translate the new direction vector by the camera matrix
			
			//glm::vec3 dir = multDirMatrix(options.cameraToWorld, glm::vec3(x, y, -1));
			glm::vec3 dir(x, y, -1);
			dir = glm::normalize(dir);

			*(pix++) = CastRay(origin, dir, objects, lights, options, 0);
		}
	}
	// Save result to a PPM image (keep these flags if you compile under Windows)
	//std::ofstream ofs("./out.ppm", std::ios::out | std::ios::binary);
	//ofs << "P6\n" << options.width << " " << options.height << "\n255\n";
	//for (uint32_t i = 0; i < options.height * options.width; ++i) {
	//	char r = (char)(255 * clamp(0, 1, framebuffer[i].x));
	//	char g = (char)(255 * clamp(0, 1, framebuffer[i].y));
	//	char b = (char)(255 * clamp(0, 1, framebuffer[i].z));
	//	ofs << r << g << b;
	//}

	//ofs.close();


}

int main()
{

	// creating the scene (adding objects and lights)
	std::vector<std::unique_ptr<Object>> objects;
	std::vector<std::unique_ptr<Light>> lights;

	// generate a scene made of random spheres
	uint32_t numSpheres = 10;
	gen.seed(0);
	//for (uint32_t i = 0; i < numSpheres; ++i) {
	//	glm::vec3 randPos((0.5 - dis(gen)) * 10, (0.5 - dis(gen)) * 10.0f, (0.5 + dis(gen) * 10));
	//	float randRadius = (0.5 + dis(gen) * 0.5);
	//	objects.push_back(std::unique_ptr<Object>(new Sphere(randPos, randRadius)));
	//}
	//objects.push_back(std::unique_ptr<Object>(new Plane(glm::vec3(0.0f, 0.0f, 0.0f),
	//	glm::vec3(0.0f, 0.0f, 1.0f))));
	//objects.push_back(std::unique_ptr<Object>(new Disk(glm::vec3(0.0f, 0.0f, 0.0f),1.0f)));
	objects.push_back(std::unique_ptr<Object>(new Triangle(
		glm::vec3(-1.0f, -1.0f, 9.0f),
		glm::vec3(1.0f, -1.0f,	9.0f),
		glm::vec3(0.0f, 1.0f,	9.0f))));

	// setting up options
	Options options;
	options.width = 640;
	options.height = 480;
	options.fov = 90;
	options.cameraToWorld = glm::translate(glm::mat4(1.0f),
		glm::vec3(0.0f,0.0f,13.0f)); //glm::mat4(0.945519, 0, -0.325569, 0, -0.179534, 0.834209, -0.521403, 0, 0.271593, 0.551447, 0.78876, 0, 4.208271, 8.374532, 17.932925, 1);




	//Render
	Render(options, objects, lights);

	//Create Window
	MainWindow* mainWindow = new MainWindow(options.width, options.height);
	mainWindow->Initialize();

	//Draw Scene
	mainWindow->Draw(framebuffer, options.width, options.height);

	bool quitApplication = false;
	SDL_Event e;
	while (!quitApplication)
	{
		//Handle events on queue
		while (SDL_PollEvent(&e) != 0)
		{
			//User requests quit
			if (e.type == SDL_QUIT)
			{
				quitApplication = true;
			}
		}

	}
	mainWindow->Destroy();

	delete[] framebuffer;
	return 0;

}
