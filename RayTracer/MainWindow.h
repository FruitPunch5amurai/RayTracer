#pragma once
#define SDL_MAIN_HANDLED
#include <SDL2/SDL.h>
#include <glm.hpp>
class MainWindow {
public:
	MainWindow(const int& screenWidth, const int& screenHeight) :
		m_screenWidth(screenWidth), m_screenHeight(screenHeight) {}
	~MainWindow() {}

	int Initialize();
	void Draw(glm::vec3* framebuffer,int screenwidth,int screenheight);
	void Destroy();

private:
	int m_screenWidth;
	int m_screenHeight;

	SDL_Window* m_window = nullptr;
	SDL_Renderer* m_renderer = nullptr;
};
