#include "MainWindow.h"
#include "common.h"
#include <iostream>
	int MainWindow::Initialize()
	{
		if (SDL_Init(SDL_INIT_VIDEO) < 0)
		{
			printf("SDL could not initialize! SDL_Error: %s", SDL_GetError());
			return -1;
		}
		m_window = SDL_CreateWindow("Ray Tracer", 100, 100, m_screenWidth, m_screenHeight, SDL_WINDOW_SHOWN);
		if (m_window == NULL)
		{
			printf("SDL could not initialize! SDL_Error: %s", SDL_GetError());
			return -1;
		}
		//Get window surface
		m_renderer = SDL_CreateRenderer(m_window,-1,0);

		//Fill the surface white
		//SDL_FillRect(m_surface, NULL, 
		//	SDL_MapRGB(m_surface->format, 0xFF, 0xFF, 0xFF));

		//Update the surface
		//SDL_UpdateWindowSurface(m_window);

		//Wait two seconds
		//SDL_Delay(2000);

		return 1;
	}

	void MainWindow::Draw(glm::vec3* framebuffer, int screenwidth, int screenheight)
	{
		SDL_SetRenderDrawColor(m_renderer, 0, 0, 0, 255);
		SDL_RenderClear(m_renderer);

		for (int y = 0; y < screenheight; ++y)
		{
			for (int x = 0; x < screenwidth; ++x)
			{
				int i = x + y * m_screenWidth;
				int r = (255 * clamp(0, 1, framebuffer[i].x));
				int g = (255 * clamp(0, 1, framebuffer[i].y));
				int b = (255 * clamp(0, 1, framebuffer[i].z));
				SDL_SetRenderDrawColor(m_renderer, r, g, b, 255);
				SDL_RenderDrawPoint(m_renderer, x, y);
			}
		}

		SDL_RenderPresent(m_renderer);
	}

	void MainWindow::Destroy()
	{
		SDL_DestroyRenderer(m_renderer);

		//Destroy window
		SDL_DestroyWindow(m_window);

		//Quit SDL subsystems
		SDL_Quit();
	}
