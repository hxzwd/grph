#include "main.hpp"

void DrawText(std::string text, int x, int y, bool Flag = false, std::vector<int> color = {255, 255, 255}, void *font = GLUT_BITMAP_HELVETICA_12)
{
	glColor3f((float)color[0]/255.0, (float)color[1]/255.0, (float)color[2]/255.0);
	const char *c_style = text.c_str();
	if(!Flag)
	{
		x = -(WINDOW_SIZE_ONE / 2 - x);
		y = -(-WINDOW_SIZE_TWO / 2 + y);
	}
	glRasterPos2f((float)x/(0.5 * WINDOW_SIZE_ONE), (float)y/(0.5 * WINDOW_SIZE_TWO));
	for(int i = 0; i < text.size(); i++)
		glutBitmapCharacter(font, c_style[i]);
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

void DrawCircle(int xc, int yc, int R, std::vector<int> color = {255, 255, 255})
{

	int D, delta;
	int y = R;
	int x = 0;
	int mh, md, mv;
	
	for(; x <= R && y >= 0; )
	{
		mh = std::abs((x + 1)*(x + 1) + y*y - R*R);
		md = std::abs((x + 1)*(x + 1) + (y - 1)*(y - 1) - R*R);
		mv = std::abs(x*x + (y - 1)*(y - 1) - R*R);
		D = (x + 1)*(x + 1) + (y - 1)*(y - 1) - R*R;
		if(D < 0)
		{
			delta = mh - md;
			if(delta <= 0)
			{
				x = x + 1;
				y = y;
				SetPoint(x + xc, y + yc, color);
				SetPoint(-x + xc, -y + yc, color);
				SetPoint(x + xc, -y + yc, color);
				SetPoint(-x + xc, y + yc, color);
			}
			else
			{
				x = x + 1;
				y = y - 1;
				SetPoint(x + xc, y + yc, color);
				SetPoint(-x + xc, -y + yc, color);
				SetPoint(x + xc, -y + yc, color);
				SetPoint(-x + xc, y + yc, color);
			}
		}
		if(D > 0)
		{
			delta = md - mv;
			if(delta <= 0)
			{
				x = x + 1;
				y = y - 1;
				SetPoint(x + xc, y + yc, color);
				SetPoint(-x + xc, -y + yc, color);
				SetPoint(x + xc, -y + yc, color);
				SetPoint(-x + xc, y + yc, color);
			}
			else
			{
				x = x;
				y = y - 1;
				SetPoint(x + xc, y + yc, color);
				SetPoint(-x + xc, -y + yc, color);
				SetPoint(x + xc, -y + yc, color);
				SetPoint(-x + xc, y + yc, color);
			}
		}
		if(D == 0)
		{
			x = x + 1;
			y = y - 1;
			SetPoint(x + xc, y + yc, color);
			SetPoint(-x + xc, -y + yc, color);
			SetPoint(x + xc, -y + yc, color);
			SetPoint(-x + xc, y + yc, color);
		}
	
	}
	
		
}


IntersectType Intersect(double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy, double *t)
{
	double nx = dy - cy;
	double ny = cx - dx;
	PointPosition type;
	double denom = nx*(bx - ax) + ny*(by - ay);
	if(denom == 0.0f)
	{
		type = PointPos(cx, cy, dx, dy, ax, ay);
		if(type == LEFT || type == RIGHT)
			return PARALLEL;
		else
			return COLLINEAR;
	}
	double num = nx*(ax - cx) + ny*(ay - cy);
	*t = -num/denom;
	return SKEW;
}

IntersectType Intersect(const std::vector<double> Points, double *t)
{
	return Intersect(Points[0], Points[1], Points[2], Points[3], Points[4], Points[5], Points[6], Points[7], t);
}

IntersectType CROSS(double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy, double *tab, double *tcd)
{
	IntersectType type = Intersect(ax, ay, bx, by, cx, cy, dx, dy, tab);
	if(type == COLLINEAR || type == PARALLEL) return type;
	if((*tab < 0) || (*tab > 1)) return SKEW_NO_CROSS;
	Intersect(cx, cy, dx, dy, ax, ay, bx, by, tcd);
	if((*tcd < 0) || (*tcd > 1)) return SKEW_NO_CROSS;
	return SKEW_CROSS;
}

IntersectType CROSS(const std::vector<double> Points, double *tab, double *tcd)
{
	return CROSS(Points[0], Points[1], Points[2], Points[3], Points[4], Points[5], Points[6], Points[7], tab, tcd);
}

void FlatCaps(std::vector<int> c, int d, const std::vector<int> col = {255, 255, 255}, unsigned int style = 0, bool coords = false)
{
	std::vector<int> H3;
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
	HH1 = {A, D, A, -D, -A, -D, -A, D};
	
	cs = std::cos(angle);
	sn = std::sin(angle);
	
	for(int i = 0; i < HH1.size() - 1; )
	{
		H3.push_back(std::ceil(X0 + cs*HH1[i] - sn*HH1[i + 1]));
		H3.push_back(std::ceil(Y0 + sn*HH1[i] + cs*HH1[i + 1]));
		i = i + 2;
	}

	DrawPoly(H3, col, true);
}

void QuadCaps(std::vector<int> c, int d = 10, const std::vector<int> col = {255, 255, 255}, const std::vector<int> coll = {255, 255, 255}, unsigned int style = 0, bool coords = false)
{
	Draw_2(c, d, col, coll, style, coords);
}


void EvenOdd(std::vector<int> poly, std::vector<int> color = {255, 255, 255}, bool Flag = false)
{
	std::vector<int> xx;
	std::vector<int> yy;
	int n = poly.size() / 2;
	for(int i = 0; i < n * 2; i+=2)
		xx.push_back(poly[i]);
	for(int i = 1; i < n * 2; i+=2)
		yy.push_back(poly[i]);
	int x_max = *(std::max_element(xx.begin(), xx.end()));
	int x_min = *(std::min_element(xx.begin(), xx.end()));
	int y_max = *(std::max_element(yy.begin(), yy.end()));
	int y_min = *(std::min_element(yy.begin(), yy.end()));
	
	DrawPoly(poly, color, Flag);

	int dc_x, dc_y;
	double tab, tcd;
	IntersectType type;
	int counter;
	for(int x = x_min + 1; x < x_max; x++)
	{

		for(int y = y_min + 1; y < y_max; y++)
		{
			if((y - y_min) < (y_max - y_min)/2) { dc_x = x_max + 5; dc_y = y_max - 5; }
			else { dc_x = x_min - 5; dc_y = y_min + 5; }
			
			counter = 0;
			
			for(int i = 0; i < n - 1; i++)
			{
				type = CROSS(xx[i], yy[i], xx[i + 1], yy[i + 1], x, y, dc_x, dc_y, &tab, &tcd);
				if(type == SKEW_CROSS) counter++;
			}
			type = CROSS(xx[0], yy[0], xx[n - 1], yy[n - 1], x, y, dc_x, dc_y, &tab, &tcd);
				if(type == SKEW_CROSS) counter++;
			if((counter % 2) == 1) (Flag ? (SetPoint) : (SetPixel))(x, y, color);
		}

	}

}

void ShowPolyEdges(std::vector<std::vector<int>> edges)
{
	int counter = 1;
	std::cout << "Num of edges is " << edges.size() << "\n";
	for(auto it = edges.begin(); it != edges.end(); it++, counter++)
	{
		std::cout << counter << ". " << "(" << it->at(0) << ", " << it->at(1) << ")";
		std::cout << "\t--->" << "(" << it->at(2) << ", " << it->at(3) << ")\n";
	}
}

void NZW(std::vector<int> poly, std::vector<int> color = {255, 255, 255}, bool Flag = false)
{
	std::vector<int> xx;
	std::vector<int> yy;
	int n = poly.size() / 2;
	for(int i = 0; i < n * 2; i+=2)
		xx.push_back(poly[i]);
	for(int i = 1; i < n * 2; i+=2)
		yy.push_back(poly[i]);
	int x_max = *(std::max_element(xx.begin(), xx.end()));
	int x_min = *(std::min_element(xx.begin(), xx.end()));
	int y_max = *(std::max_element(yy.begin(), yy.end()));
	int y_min = *(std::min_element(yy.begin(), yy.end()));

	DrawPoly(poly, color, Flag);
	
	int counter;
	int ray_x = x_max + 5;
	int ray_y;
	std::vector<std::vector<int>> edges;
	double tab, tcd;
	IntersectType type;
	PointPosition pos;
	int wn;

	for(int i = 0; i < n; i++)
	{	
		if(i == n - 1)
			edges.push_back({xx[i], yy[i], xx[0], yy[0]});
		else
			edges.push_back({xx[i], yy[i], xx[i + 1], yy[i + 1]});
	}

	ShowPolyEdges(edges);

	for(int x = x_min; x <= x_max; x++)
	{

		for(int y = y_min; y <= y_max; y++)
		{
			wn = 0;
			ray_y = y;
			for(int i = 0; i < n; i++)
			{
				type = CROSS(x, y, ray_x, ray_y, (edges[i])[0], (edges[i])[1], (edges[i])[2], (edges[i])[3], &tab, &tcd);
				if(type == SKEW_CROSS)
				{
					pos = PointPos(edges[i], (double)x, (double)y);
					if(pos == RIGHT) wn++;
					if(pos == LEFT) wn--;
				}
			}
			if(wn != 0) (Flag ? (SetPoint) : (SetPixel))(x, y, color);
		}

	}
}

void NamedPointRN(int x, int y, std::vector<int> color = {255, 255, 255}, std::vector<int> color2 = {255, 255, 255})
{
	SetPoint(x, y, color);
	std::string str = "(";
	str = str + std::to_string(x) + std::string(", ") + std::to_string(y) + std::string(")");
	DrawText(str, x + 1, y + 1, true, color2, GLUT_BITMAP_8_BY_13);
}

void NamedPointLN(int x, int y, std::vector<int> color = {255, 255, 255}, std::vector<int> color2 = {255, 255, 255})
{
	SetPoint(x, y, color);
	std::string str = "(";
	str = str + std::to_string(x) + std::string(", ") + std::to_string(y) + std::string(")");
	int len = str.size();
	DrawText(str, x - 8 * len - 1, y - 9, true, color2, GLUT_BITMAP_8_BY_13);
}

void NamedPoint2N(int x, int y, int a, int b, std::vector<int> color = {255, 255, 255}, std::vector<int> color2 = {255, 255, 255})
{
	int yy, aa, bb;
	int xx = std::max(x, a);
	if(xx == a) { yy = b; aa = x; bb = y; }
	else { yy = y; aa = a; bb = b; }
	if(yy >= bb)
	{
		NamedPointRN(xx, yy, color, color2);
		NamedPointLN(aa, bb, color, color2);
	}
	else
	{
		NamedPointLN(xx, yy, color, color2);
		NamedPointRN(aa, bb, color, color2);
	}
}	

void NamedLineN(int x, int y, int a, int b, std::vector<int> color = {255, 255, 255}, std::vector<int> color2 = {255, 255, 255})
{
	DrawLineN(x, y, a, b, color);
	NamedPoint2N(x, y, a, b, color, color2);
}


void Test_DrawBezier2(std::vector<int> points, std::vector<int> color = {255, 255, 255}, bool Flag = false)
{
	double t = 0.0f;
	double step = 0.01f;
	double x;
	double y;
	int xx, yy;
	int steps = (int)(1.0f/step);
	for(int i = 1; i <= steps; i++)
	{
		x = BernsteinPoly(2, 0, t)*points[0] + BernsteinPoly(2, 1, t)*points[2] + BernsteinPoly(2, 2, t)*points[4];
		y = BernsteinPoly(2, 0, t)*points[1] + BernsteinPoly(2, 1, t)*points[3] + BernsteinPoly(2, 2, t)*points[5];
		xx = std::round(x);
		yy = std::round(y);
		(Flag ? (SetPoint) : (SetPixel))(xx, yy, color);
		t = t + step;
	}
}

void DrawBezier2(std::vector<int> points, std::vector<int> color = {255, 255, 255}, bool Flag = false)
{
	double a, b, c, S, p;
	a = std::sqrt((points[0] - points[2])*(points[0] - points[2]) + (points[1] - points[3])*(points[1] - points[3]));
	b = std::sqrt((points[2] - points[4])*(points[2] - points[4]) + (points[3] - points[5])*(points[3] - points[5]));
	c = std::sqrt((points[0] - points[4])*(points[0] - points[4]) + (points[1] - points[5])*(points[1] - points[5]));
	p = (a + b + c)/2.0;
	S = std::sqrt(p * (p - a) * (p - b) * (p - c));
	if(S <= 0.5)
	{
		DrawLineN(points[0], points[1], points[2], points[3], color);
		DrawLineN(points[2], points[3], points[4], points[5], color);
		std::cout << "Square is " << S << "\n";
		return;
	}
		
	int mx, my, ma, mb;
	mx = (points[0] + points[2])/2;
	my = (points[1] + points[3])/2;
	ma = (points[4] + points[2])/2;
	mb = (points[5] + points[3])/2;
	DrawBezier2({points[0], points[1], mx, my, (mx + ma)/2, (my + mb)/2}, color, Flag);
	DrawBezier2({(mx + ma)/2, (my + mb)/2, ma, mb, points[4], points[5]}, color, Flag);
}


