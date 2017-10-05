#ifndef H_MAIN
#define H_MAIN

#pragma once

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <functional>
#include <algorithm>
#include <map>

#include "GL/glut.h"

#include "mapp.hpp"

#define WINDOW_SIZE_ONE 800
#define WINDOW_SIZE_TWO 600
#define WINDOW_POSITION_ONE 100
#define WINDOW_POSITION_TWO 100
#define WINDOW_NAME "Window"

#define RED {255, 0, 0}
#define BLUE {0, 0, 255}
#define GREEN {0, 255, 0}
#define WHITE {255, 255, 255}
#define BLACK {0, 0, 0}
#define YELLOW {255, 255, 0}
#define CYAN {0, 255, 255}
#define MAGENTA {255, 0, 255}

#define DEBUG0(X) { std::cout << "DEBUG MACRO START\n";				\
		std::cout << "PLACE IS " << X << "\n";				\
		std::cout << "x1 y1: " << x1 << "\t" << y1 << "\n"; 		\
		std::cout << "x2 y2: " << x2 << "\t" << y2 << "\n"; 		\
		std::cout << "x y  : " << x << "\t" << y << "\n"; 		\
		std::cout << "d    : " << d << "\n"; 				\
		std::cout << "1/d  : " << 1.0/d << "\n"; 			\
		std::cout << "dx dy: " << deltax << "\t" << deltay << "\n"; 	\
		std::cout << "DEBUG MACRO END\n"; }

#define DEBUG(X, Y) DEBUG0(X)

#define PINFO0(X, Y) { std::cout << "Point info: \n";					\
			std::cout << "(x, y):\t( " << X << ", " << Y << ")\n"; 		\
			std::cout << "x = " << X << "\t;y = " << Y << "\n"; }		\

#define PINFO(X, Y) PINFO0(X, Y)

#undef DINFO

typedef struct
{
	int x;
	int y;
	std::vector<int> color;
	void Coords(const std::vector<int> value = {255, 255, 255})
	{
		x = y = 0;
		color = value;
	}
	void Coords(int a, int b, int c)
	{
		x = y = 0;
		color.push_back(c);
		color.push_back(b);
		color.push_back(a);
	}
	void ReCalc(void)
	{
		x = x - WINDOW_SIZE_ONE / 2;
		y = WINDOW_SIZE_TWO / 2 - y;
	}
	void Calc(void)
	{
		x = -(WINDOW_SIZE_ONE / 2 - x);
		y = -(-WINDOW_SIZE_TWO / 2 + y);
	}
	int operator[](int counter)
	{
		if(counter == 0) return x;
		else return y;
	}
} Coords;

void ActiveMouseFunction(int x, int y);
void PassiveMouseFunction(int x, int y);
void DisableActiveMouse();
void DisablePassiveMouse();
void EnableActiveMouse();
void EnablePassiveMouse();
void DisableMouseButton();
void EnableMouseButton();
void MouseButton(int button, int state, int x, int y);

typedef std::function<void(int, int, int, int, std::vector<int>)> DF0;
typedef std::function<void(std::vector<int>, std::vector<int>)> DF1;
typedef std::function<double(double)> FF;
typedef std::function<int(int)> ff;

inline void SetPoint(int x, int y, const std::vector<int> color);
inline void SetPixel(int x, int y, const std::vector<int> color);
void DrawLine(int x1, int y1, int x2, int y2, const std::vector<int> color);
bool CheckFlag(unsigned int Flags, int num);
unsigned int PrepareFlags(std::vector<unsigned int> Options);
void ClearScreen(const std::vector<int> color);
void SetUserSettings();

void DrawText(std::string text, int x, int y, bool Flag, std::vector<int> color, void *font);
void NamedPointRN(int x, int y, std::vector<int> color, std::vector<int> color2);
void NamedPointLN(int x, int y, std::vector<int> color, std::vector<int> color2);
void NamedPoint2N(int x, int y, int a, int b, std::vector<int> color, std::vector<int> color2);
void NamedLineN(int x, int y, int a, int b, std::vector<int> color, std::vector<int> color2);


enum PointPosition { NONE = 0, LEFT, RIGHT, BEYOND, BEHIND, BETWEEN, ORIGIN, DESTINATION };
enum IntersectType { COLLINEAR = 0, PARALLEL, SKEW, SKEW_CROSS, SKEW_NO_CROSS };
enum TypeMouseSettings { ACTIVE_MOUSE_OFF = 1, PASSIVE_MOUSE_OFF, MOUSE_BUTTON_OFF };


void PrepareMouse(unsigned int Flags);


PointPosition PointPos(double x1, double y1, double x2, double y2, double x, double y);
PointPosition PointPos(const std::vector<int> Points);
PointPosition PointPos(const std::vector<int> Points, double x, double y);
IntersectType Intersect(double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy, double *t);
IntersectType Intersect(const std::vector<double> Points, double *t);
IntersectType CROSS(double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy, double *tab, double *tcd);
IntersectType CROSS(const std::vector<double> Points, double *tab, double *tcd);

void DrawLineN(int x1, int y1, int x2, int y2, const std::vector<int> color);
void DrawPoly(const std::vector<int> Point, const std::vector<int> color, bool Flag);
void DrawLineNV(const std::vector<int> Point, const std::vector<int> color);
void DrawLineV(const std::vector<int> Point, const std::vector<int> color);
void DrawCircle(int xc, int yc, int R, std::vector<int> color);

bool IsConvex(const std::vector<int> Poly);
bool IsConvex2(const std::vector<int> Poly);
bool IsSimple(std::vector<double> Poly);

void DrawBarcodeLine(int x1, int y1, int x2, int y2, double shift, const std::vector<int> color);
void DrawBarcodeLineN(int x1, int y1, int x2, int y2, double shift,const std::vector<int> color);

void DrawLFunI(ff Fnc, const std::vector<int> p,  const std::vector<int> color);
void DrawLFunD(FF Fnc, const std::vector<double> p, const std::vector<int> color);
std::vector<int> DrawSquareN(int x, int y, int d, const std::vector<int> col, double angle);
void Draw_0(std::vector<int> C, int d, const std::vector<int> a, const std::vector<int> b, unsigned int style, bool coords);
void Draw_1(std::vector<int> C, int d, const std::vector<int> a, const std::vector<int> b, unsigned int style, bool coords);
void Draw_2(std::vector<int> c, int d, const std::vector<int> col, const std::vector<int> coll, unsigned int style, bool coords);
double CalcAngle(std::vector<int> C);

void FlatCaps(std::vector<int> c, int d, const std::vector<int> col, unsigned int style, bool coords);
void QuadCaps(std::vector<int> c, int d, const std::vector<int> col, const std::vector<int> coll, unsigned int style, bool coords);

void EvenOdd(std::vector<int> poly, std::vector<int> color, bool Flag);
void NZW(std::vector<int> poly, std::vector<int> color, bool Flag);
void ShowPolyEdges(std::vector<std::vector<int>> edges);

void Test_DrawBezier2(std::vector<int> points, std::vector<int> color, bool Flag);
void DrawBezier2(std::vector<int> points, std::vector<int> color, bool Flag);

#endif
