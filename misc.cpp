#include "main.hpp"

bool CheckFlag(unsigned int Flags, int num)
{
	unsigned int pattern = 0x1;
	pattern = pattern << (num - 1);
	if(Flags & pattern != 0) return true;
	else return false;
}

unsigned int PrepareFlags(std::vector<unsigned int> Options)
{
	unsigned int result = 0;
	for(int i = 0; i < Options.size(); i++)
	{
		unsigned int one = 1;
		one = one << Options[i];
		result = result | one;
	}
	return result;
}

inline void SetPoint(int x, int y, const std::vector<int> color = { 255, 255, 255 })
{
	
	glBegin(GL_POINTS);
		glColor3f((float)color[0]/255.0, (float)color[1]/255.0, (float)color[2]/255.0);
		glVertex3f((float)x/(0.5 * WINDOW_SIZE_ONE), (float)y/(0.5 * WINDOW_SIZE_TWO), 0);
	glEnd();
}

inline void SetPixel(int x, int y, const std::vector<int> color = { 255, 255, 255 })
{
	
	x = x - WINDOW_SIZE_ONE / 2;
	y = WINDOW_SIZE_TWO / 2 - y;

	glBegin(GL_POINTS);
		glColor3f((float)color[0]/255.0, (float)color[1]/255.0, (float)color[2]/255.0);
		glVertex3f((float)x/(0.5 * WINDOW_SIZE_ONE), (float)y/(0.5 * WINDOW_SIZE_TWO), 0);
	glEnd();
}


void DrawLine(int x1, int y1, int x2, int y2, const std::vector<int> color = { 255, 255, 255 })
{
	double Error = -0.5;
	int deltax = x2 - x1;
	int deltay = y2 - y1;
	int x, y;
	int DY;
	
	if(deltax == 0 && deltay == 0)
	{
		SetPixel(x1, y1, color);
		return;
	}

	if(deltax == 0)
	{
		x = x1;
		y = std::min(y1, y2);
		
		for(int i = 0; i <= std::abs(deltay); i++, y++)
			SetPixel(x, y, color);
		
		return;
	}

	if(deltay == 0)
	{
		y = y1;
		x = std::min(x1, x2);
		
		for(int i = 0; i <= std::abs(deltax); i++, x++)
			SetPixel(x, y, color);
		
		return;
	}
	
	x = std::min(x1, x2);
	if(x == x1) { y = y1; DY = -y2 + y1; }
	else	{ y = y2; DY = -y1 + y2; }

	double d = (double)DY/(double)(std::max(x1, x2) - std::min(x1, x2));
	
	
	if(d > 0.0f && d <= 1.0f)
	{

		for(int i = 0; i <= std::abs(deltax); i++)
		{
			Error = Error + d;
			if(Error < 0) y = y;
			if(Error >= 0)
			{
				y = y - 1;
				Error = Error - 1;
			}
			SetPixel(x, y, color);
			x = x + 1;
		}
		return;
	}


	if(d > 1.0f)
	{
		y = std::max(y1, y2);
		if(y == y1) x = x1;
		else x = x2;

		for(int i = 0; i <= std::abs(deltay); i++)
		{
			Error = Error + 1.0/d;
			if(Error < 0) x = x;
			if(Error >= 0)
			{
				x = x + 1;
				Error = Error - 1;
			}
			SetPixel(x, y, color);
			y = y - 1;
		}
		return;
	}
	
	x = std::min(x1, x2);
	if(x == x1) { y = y1, DY = y2 - y1; }
	else { y = y2, DY = y1 - y2; }

	d = (double)DY/(std::min(x1, x2) - std::max(x1, x2));


	if(d < -1.0f)
	{
		int xx;
		
		y = std::max(y1, y2);
		if(y == y1) { x = x1; xx = x2; }
		else { x = x2; xx = x1; }
		
		int sign = 1;
		if(x <= xx) sign = 1;
		else sign = -1;

		for(int i = 0; i <= std::abs(deltay); i++)
		{
			Error = Error - 1.0/d;
			if(Error < 0) x = x;
			if(Error >= 0)
			{
				x = x + sign * 1;
				Error = Error - 1;
			}
			SetPixel(x, y, color);
			y = y - 1;
		}
		return;
	}

	if(d >= -1.0f && d < 0.0f)
	{
		
		int sign;
		
		y = std::max(y1, y2);
		if(y == y1) { x = x1; (x1 > x2) ? sign = -1:sign = 1; }
		else { x = x2; (x2 > x1) ? sign = -1:sign = 1; }		

		for(int i = 0; i <= std::abs(deltax); i++)
		{
			Error = Error - d;
			if(Error < 0) y = y;
			if(Error >= 0)
			{
				y = y - 1;
				Error = Error - 1;
			}
			SetPixel(x, y, color);
			x = x + sign * 1;
		}
		return;
	}
	
}

void DrawLineN(int x1, int y1, int x2, int y2, const std::vector<int> color = { 255, 255, 255 })
{
	DrawLine(x1 + WINDOW_SIZE_ONE / 2, WINDOW_SIZE_TWO / 2 - y1, x2 + WINDOW_SIZE_ONE / 2, WINDOW_SIZE_TWO / 2 - y2, color);
}

void DrawBarcodeLine(int x1, int y1, int x2, int y2, double shift = 2,const std::vector<int> color = {255, 255, 255})
{
	double length = std::sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
	int steps = (int)std::round(length/(shift));
	
	
	std::function<int(double)> X = [=](double t)->int { return (int)std::round(x1 + t*(x2 - x1)); };
	std::function<int(double)> Y = [=](double t)->int { return (int)std::round(y1 + t*(y2 - y1)); };
	int x = X(0);
	int y = Y(0);
	int a, b;
	
	for(int i = 1; i <= steps; i++)
	{
		a = X(1.0/(2 * steps) * (double)i);
		b = Y(1.0/(2 * steps) * (double)i);
		if((i % 2) == 1) DrawLine(x, y, a, b, color);
		x = a;
		y = b;
	}
	
}

void DrawBarcodeLineN(int x1, int y1, int x2, int y2, double shift = 2,const std::vector<int> color = {255, 255, 255})
{
	DrawBarcodeLine(x1 + WINDOW_SIZE_ONE / 2, WINDOW_SIZE_TWO / 2 - y1, x2 + WINDOW_SIZE_ONE / 2, WINDOW_SIZE_TWO / 2 - y2, shift, color);
}

void DrawLineNV(const std::vector<int> Point, const std::vector<int> color)
{
	DrawLineN(Point[0], Point[1], Point[2], Point[3], color);
}


void DrawLineV(const std::vector<int> Point, const std::vector<int> color)
{
	DrawLine(Point[0], Point[1], Point[2], Point[3], color);
}

PointPosition PointPos(double x1, double y1, double x2, double y2, double x, double y)
{
	double ax = x2 - x1;
	double ay = y2 - y1;
	double bx = x - x1;
	double by = y - y1;
	double s = ax * by - bx * ay;
	if(s > 0.0f) return LEFT;
	if(s < 0.0f) return RIGHT;
	if(ax * bx < 0.0f || ay * by < 0.0f)
		return BEHIND;
	if(ax * ax + ay * ay < bx * bx + by * by)
		return BEYOND;
	if(x2 == x && y2 == y) return DESTINATION;
	if(x1 == x && y1 == y) return ORIGIN;
	return BETWEEN;

}

PointPosition PointPos(const std::vector<int> Points)
{
	if(Points.size() < 6) return NONE;
	return PointPos(Points[0], Points[1], Points[2], Points[3], Points[4], Points[5]);
}

PointPosition PointPos(const std::vector<int> Points, double x, double y)
{
	return PointPos(Points[0], Points[1], Points[2], Points[3], x, y);
}

void DrawPoly(const std::vector<int> Point, const std::vector<int> color  = { 255, 255, 255 }, bool Flag = false)
{
	int PolySize = Point.size();
	if(PolySize % 2 != 0) PolySize = PolySize - 1;
	
	int i = 0, j = 2;
	for(i = 0, j = 2; i < PolySize && j < PolySize; )
	{
		(Flag ? (DrawLineN) : (DrawLine))(Point[i], Point[i + 1], Point[j], Point[j + 1], color);
		i += 2;
	 	j += 2;
	}
	j = j - 2;
	(Flag ? (DrawLineN) : (DrawLine))(Point[0], Point[1], Point[j], Point[j + 1], color);
}

double CalcAngle(std::vector<int> x)
{
	int a = std::max(x[0], x[2]);
	int aa = std::min(x[0], x[2]);
	int b;
	int bb;
	if(a == x[0]) { bb = x[3]; b = x[1]; }
	else { b = x[3]; bb = x[1]; }
	int da = a - aa;
	int db = b - bb;
	return (double)db/(double)da;
}

std::vector<int> DrawSquareN(int x, int y, int d, const std::vector<int> col = { 255, 255, 255 }, double angle = 0)
{
	int a;
	double dist, unit, ddd;
	std::vector<int> result;
	ddd = (double)std::sqrt(d * d + d * d);
	dist = std::max(WINDOW_SIZE_ONE, WINDOW_SIZE_TWO);
	unit = dist/dist;
	
	if(angle != 0)
	{
		double g, h;
		int gg, hh;
		g = (ddd/unit) * std::cos(angle);
		h = (ddd/unit) * std::sin(angle);
		gg = std::floor(std::abs(g));
		hh = std::floor(std::abs(h));		
		result = { x + hh, y - gg, x + gg, y + hh, x - hh, y + gg, x - gg, y - hh};
		DrawPoly({ x + hh, y - gg, x + gg, y + hh, x - hh, y + gg, x - gg, y - hh}, col, true);
		return result;

	}

	a = std::ceil(ddd/unit);
	
	result = { x, y - a, x + a, y, x, y + a, x - a, y};
	DrawPoly({ x, y - a, x + a, y, x, y + a, x - a, y}, col, true);
	return result;
}

void DrawLFunI(ff Fnc, const std::vector<int> p = {WINDOW_SIZE_ONE / 2, -WINDOW_SIZE_ONE / 2 }, const std::vector<int> color = { 255, 255, 255 })
{
	int a = Fnc(p[0]);
	int b = Fnc(p[1]);
	DrawLineN(p[0], a, p[1], b, color);
}

void DrawLFunD(FF Fnc, const std::vector<double> p = {WINDOW_SIZE_ONE / 2.0, -WINDOW_SIZE_ONE / 2.0  }, const std::vector<int> color = { 255, 255, 255 })
{
	int a = std::floor(Fnc(p[0]));
	int b = std::floor(Fnc(p[1]));
	DrawLineN((int)p[0], a, (int)p[1], b, color);
}


void Draw_0(std::vector<int> c, int d = 10, const std::vector<int> a = { 255, 255, 255 }, const std::vector<int> b = { 255, 255, 255 }, unsigned int style = 0, bool coords = false)
{
	
	std::vector<int> H1, H2;
	if(!coords)
	{
		for(int i = 0; i < c.size(); i++)
		{
			if(i % 2 == 0) c[i] =  c[i] - WINDOW_SIZE_ONE / 2;
			if(i % 2 == 1) c[i] =  WINDOW_SIZE_TWO / 2 - c[i];
		}
	}
	
	if(CheckFlag(style, 1)) DrawLineN(c[0], c[1], c[2], c[3], a);
	
	double angle = CalcAngle(c);
	angle = std::atan(angle);

	H1 = DrawSquareN(c[0], c[1], d, b, angle);	
	H2 = DrawSquareN(c[2], c[3], d, b, angle);
	DrawLineN(H1[0], H1[1], H2[0], H2[1], a);
	DrawLineN(H1[4], H1[5], H2[4], H2[5], a);

}

void Draw_1(std::vector<int> c, int d = 10, const std::vector<int> col = { 255, 255, 255 }, const std::vector<int> coll = { 255, 255, 255 }, unsigned int style = 0, bool coords = false)
{
	
	std::vector<int> H1, H2;
	double angle, h, k, b, w, X, Y;
	if(!coords)
	{
		for(int i = 0; i < c.size(); i++)
		{
			if(i % 2 == 0) c[i] =  c[i] - WINDOW_SIZE_ONE / 2;
			if(i % 2 == 1) c[i] =  WINDOW_SIZE_TWO / 2 - c[i];
		}
	}

	if(CheckFlag(style, 1)) DrawLineN(c[0], c[1], c[2], c[3], col);
	angle = CalcAngle(c);
	angle = std::atan(angle);
	h = (double)d / std::cos(angle);
	k = (double)(c[1] - c[3])/(double)(c[0] - c[2]);
	b = (double)c[1] - k * (double)c[0];
	w = 2 * (double)d * std::sin(angle);
	
	X = std::max((double)c[0], (double)c[2]);
	Y = std::min((double)c[0], (double)c[2]);
	X = X - w;
	Y = Y + w;

	DrawLFunD([=](double x)->double{ return k*x + b + h; }, {X, Y}, coll);
	DrawLFunD([=](double x)->double{ return k*x + b - h; }, {X, Y}, coll);
	
}

void Draw_2(std::vector<int> c, int d = 10, const std::vector<int> col = {255, 255, 255}, const std::vector<int> coll = {255, 255, 255}, unsigned int style = 0, bool coords = false)
{
	std::vector<int> H1, H2, H3;
	std::vector<double> HH1, HH2;
	double angle, A, a, X0, Y0, D, dist, unit, cs, sn;
	if(!coords)
	{
		for(int i = 0; i < c.size(); i++)
		{
			if(i % 2 == 0) c[i] =  c[i] - WINDOW_SIZE_ONE / 2;
			if(i % 2 == 1) c[i] =  WINDOW_SIZE_TWO / 2 - c[i];
		}
	}

	if(CheckFlag(style, 1)) DrawLineN(c[0], c[1], c[2], c[3], col);
	angle = CalcAngle(c);
	angle = std::atan(angle);
	
	dist = std::max(WINDOW_SIZE_ONE, WINDOW_SIZE_TWO);
	unit = dist/dist;
	a = std::sqrt((c[0] - c[2])*(c[0] - c[2]) + (c[1] - c[3])*(c[1] - c[3]));
	D = (double)d/unit;
	A = a*0.5;
	X0 = (double)(c[0] + c[2])/2.0;
	Y0 = (double)(c[1] + c[3])/2.0;
	HH1 = {A + D, D, A - D, D, A-D, -D, A+D, -D};
	HH2 = {-A-D, D, -A-D, -D, -A+D, -D, -A+D, D};
	cs = std::cos(angle);
	sn = std::sin(angle);
	
	for(int i = 0; i < HH1.size() - 1; )
	{
		H3.push_back(std::ceil(X0 + cs*HH1[i] - sn*HH1[i + 1]));
		H3.push_back(std::ceil(Y0 + sn*HH1[i] + cs*HH1[i + 1]));
		i = i + 2;
	}
	H1 = H3;
	H3.clear();
	for(int i = 0; i < HH2.size() - 1; )
	{
		H3.push_back(std::ceil(X0 + cs*HH2[i] - sn*HH2[i + 1]));
		H3.push_back(std::ceil(Y0 + sn*HH2[i] + cs*HH2[i + 1]));
		i = i + 2;
	}
	H2 = H3;
	H3.clear();

	DrawPoly(H1, coll, true);
	DrawPoly(H2, coll, true);	
	DrawLineN(H1[2], H1[3], H2[6], H2[7], col);
	DrawLineN(H1[4], H1[5], H2[4], H2[5], col);
}

void ClearScreen(const std::vector<int> color)
{
	glClearColor((float)color[0]/255.0, (float)color[1]/255.0, (float)color[2]/255.0, 0.0f);
	glClearDepth(1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

bool PASSIVE_MOUSE_DISABLED = true;
bool ACTIVE_MOUSE_DISABLED = false;
bool MOUSE_BUTTON_DISABLED = false;

void PrepareMouse(unsigned int Flags)
{
	if(CheckFlag(Flags, ACTIVE_MOUSE_OFF))
		ACTIVE_MOUSE_DISABLED = true;
	else
		ACTIVE_MOUSE_DISABLED = false;
	if(CheckFlag(Flags, PASSIVE_MOUSE_OFF))
		PASSIVE_MOUSE_DISABLED = true;
	else
		PASSIVE_MOUSE_DISABLED = false;
	if(CheckFlag(Flags, MOUSE_BUTTON_OFF))
		MOUSE_BUTTON_DISABLED = true;
	else
		MOUSE_BUTTON_DISABLED = false;
	
}

void DisableActiveMouse() { ACTIVE_MOUSE_DISABLED = true;}
void DisablePassiveMouse() { PASSIVE_MOUSE_DISABLED = true; }
void EnableActiveMouse() { ACTIVE_MOUSE_DISABLED = false; }
void EnablePassiveMouse() { PASSIVE_MOUSE_DISABLED = false; }
void DisableMouseButton() { MOUSE_BUTTON_DISABLED = true; }
void EnableMouseButton() { MOUSE_BUTTON_DISABLED = false; }

void ActiveMouseFunction(int x, int y)
{

	if(ACTIVE_MOUSE_DISABLED)
	{
		return;
	}
	else
	{




	}

}

void PassiveMouseFunction(int x, int y)
{
	if(PASSIVE_MOUSE_DISABLED)
	{
		return;
	}
	else
	{




	}
}

void MouseButton(int button, int state, int x, int y) 
{

	if(MOUSE_BUTTON_DISABLED)
	{
		return;
	}
	else
	{
		if(button == GLUT_LEFT_BUTTON) 
		{

			

			if (state == GLUT_UP) 
			{
				
			}
			else  
			{
				
			}
		}
	}
}

void SetUserSettings()
{
	PrepareMouse(PrepareFlags({PASSIVE_MOUSE_DISABLED}));
}


#define DINFO
#undef DINFO

bool IsConvex(const std::vector<int> Poly)
{
	int size = Poly.size();
	int number_of_edges = size / 2;
	bool Convex = true;
	
	double P;
	double epsilon = 0.000001;
	double Z = VectorProductZ({(double)Poly[0] - (double)Poly[size - 2], (double)Poly[1] - (double)Poly[size - 1], 0.0}, {(double)Poly[2] - (double)Poly[0], (double)Poly[3] - (double)Poly[1], 0.0});

	#ifdef DINFO
		std::cout << "Z: " << Z << "\n";
	#endif

/*	
	if(Z == 0.0f) 
	{
		Z = VectorProductZ({(double)Poly[0] - (double)Poly[size - 2] + epsilon, (double)Poly[1] - (double)Poly[size - 1] - epsilon, 0.0}, {(double)Poly[size - 2] - (double)Poly[0] - epsilon, (double)Poly[size - 1] - (double)Poly[1]  + epsilon, 0.0});
	#ifdef DINFO
		std::cout << "Z(+epsilon): " << Z << "\n";
	#endif
	} 
	
	if(Z != 0.0f) 
		Z = Z/std::abs(Z);

*/
	if(Z == 0.0f) Z = epsilon;
	
	Z = Z/std::abs(Z);
	P = 1.0f;
	

	#ifdef DINFO
		std::cout << "P: " << P << "\n";
	#endif

	double R;
	
	int j = 1;
	for(int i = 0; i < size - 4 && j < size - 4; i+=2, j+=2)
	{
		R = VectorProductZ({(double)Poly[i + 2] - (double)Poly[i], (double)Poly[j + 2] - (double)Poly[j], 0.0}, {-(double)Poly[i + 2] + (double)Poly[i + 4], -(double)Poly[j + 2] + (double)Poly[j + 4], 0.0});

		#ifdef DINFO
		std::cout << "R: " << R << "\n";
		#endif		
		
/*
		if(R == 0.0f)
		{
			R = VectorProductZ({(double)Poly[i] - (double)Poly[i + 2] + epsilon, (double)Poly[j] - (double)Poly[j + 2] - epsilon, 0.0}, {-(double)Poly[i + 2] - (double)Poly[i + 4] - epsilon, -(double)Poly[j + 2] - (double)Poly[j + 4] + epsilon, 0.0});
		}

		if(R != 0.0f)
			R = R/std::abs(R);
*/
		R = R/std::abs(R);
		P = P * Z * R;

		if(P < 0)
		{
			Convex = false;
		}
		
	}
	if(Convex)
	{
		R = VectorProductZ({(double)Poly[size - 2] - (double)Poly[size - 4], (double)Poly[size - 1] - (double)Poly[size - 3], 0.0}, {-(double)Poly[size-2] + (double)Poly[0], -(double)Poly[size - 1] + (double)Poly[1], 0.0});
		R = R/std::abs(R);
		P = P * Z * R;
		if(P < 0)
		{
			Convex = false;
		}
	}
	return Convex;
}

#undef DINFO

bool IsConvex2(const std::vector<int> Poly)
{
	std::vector<std::pair<double, double>> edges;
	int i = 2;
	int j = 3;
	int size = Poly.size();
	
	for(; i < Poly.size(); i+=2, j+=2)
	{
		edges.push_back(std::pair<double, double>((double)(Poly[i] - Poly[i - 2]), (double)(Poly[j] - Poly[j - 2])));
	}
	edges.push_back(std::pair<double, double>((double)(Poly[0] - Poly[size - 2]), (double)(Poly[1] - Poly[size - 1])));
	
	std::vector<double> signs;
	double value;

	for(i = 1; i < edges.size(); i++)
	{
		value = VectorProductZ({ edges[i].first, edges[i].second }, { edges[i - 1].first, edges[i - 1].second });
	//	if(value != 0.0) 
		if(value == 0.0)
		{
			signs.push_back(0.0);
			continue;
		}
			signs.push_back(value/std::abs(value));
	}
	value = VectorProductZ({ edges[0].first, edges[0].second }, { edges[edges.size() - 1].first, edges[edges.size() - 1].second });
	//	if(value != 0.0)
			if(value == 0.0)
			{
				signs.push_back(0.0);
			} 
			else
			{
				signs.push_back(value/std::abs(value));
			}
//	value = ScalarProduct2({ edges[0].first, edges[0].second }, { edges[1].first, edges[1].second });
	//	if(value != 0.0) 
//			signs.push_back(value/std::abs(value));
	
	std::string tmp = "";
	for(i = 0; i < signs.size(); i++)
	{
		if(signs[i] >= 0.0)
			tmp = tmp + "+";
		else
			tmp = tmp + "-";
	}

	std::size_t p = tmp.find("+");
	std::size_t m = tmp.find("-");

	if(!(p == std::string::npos) && !(m == std::string::npos))
	{
		return false;
	}
	else
	{
		return true;
	}

}

bool IsSimple(std::vector<double> Poly)
{
	
	std::vector<std::vector<double>> edges;
	
	int k = 1;
	for(int i = 0; i < Poly.size() - 2 && k < Poly.size() - 3; i+=2, k+=2)
	{
		edges.push_back({Poly[i], Poly[k], Poly[i + 2], Poly[k + 2]});
	
	}
	
	edges.push_back({Poly[Poly.size() - 2], Poly[Poly.size() - 1], Poly[0], Poly[1]});

	IntersectType type;
	double tab, tcd;
	
	for(int i = 0; i < edges.size(); i++)
	{
		for(int j = 0; j < edges.size(); j++)
		{
			if(i == j) continue;
			type = CROSS({(edges[i])[0], (edges[i])[1], (edges[i])[2], (edges[i])[3], (edges[j])[0], (edges[j])[1], (edges[j])[2], (edges[j])[3]}, &tab, &tcd);
			
			if(type == SKEW_CROSS)
			{
				
				
				//std::cout << tab << "\t" << tcd << "\n";
				if((tab != 0.0 && tab != 1.0) && (tcd != 0.0 && tcd != 1.0))
				{
					return false;
				}


			}
		}


	}



	return true;




}

