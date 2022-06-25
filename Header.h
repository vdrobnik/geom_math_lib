#pragma once
#include <iostream>
#include <stdio.h>
#include <conio.h>
#include <windows.h>
#include <gdiplus.h>
#include <math.h>

using namespace std;

const double EPS = 1E-9;

class Intersection
{
public:
	struct pt
	{
		double x, y;

		bool operator< (const pt& p) const
		{
			return x < p.x - EPS || abs(x - p.x) < EPS && y < p.y - EPS;
		}
		bool operator> (const pt& p) const
		{
			return x > p.x - EPS || abs(x - p.x) > EPS && y > p.y - EPS;
		}
	};

	struct line
	{
		double a, b, c;

		line() {}
		line(pt p, pt q)
		{
			a = p.y - q.y;
			b = q.x - p.x;
			c = -a * p.x - b * p.y;
			norm();
		}

		void norm()
		{
			double z = sqrt(a * a + b * b);
			if (abs(z) > EPS)
				a /= z, b /= z, c /= z;
		}

		double dist(pt p) const
		{
			return a * p.x + b * p.y + c;
		}
	};

#define det(a,b,c,d)  (a*d-b*c)

	inline bool betw(double l, double r, double x)
	{
		return min(l, r) <= x + EPS && x <= max(l, r) + EPS;
	}

	inline bool intersect_1d(double a, double b, double c, double d)
	{
		if (a > b)  swap(a, b);
		if (c > d)  swap(c, d);
		return max(a, c) <= min(b, d) + EPS;
	}

	bool intersect_segments(pt a, pt b, pt c, pt d, pt& left, pt& right)
	{
		if (!intersect_1d(a.x, b.x, c.x, d.x) || !intersect_1d(a.y, b.y, c.y, d.y))
			return false;
		line m(a, b);
		line n(c, d);
		double zn = det(m.a, m.b, n.a, n.b);
		if (abs(zn) < EPS) {
			if (abs(m.dist(c)) > EPS || abs(n.dist(a)) > EPS)
				return false;
			if (b < a)
				swap(a, b);
			if (d < c)
				swap(c, d);
			left = max(a, c);
			right = min(b, d);
			return true;
		}
		else {
			left.x = right.x = -det(m.c, m.b, n.c, n.b) / zn;
			left.y = right.y = -det(m.a, m.c, n.a, n.c) / zn;
			return betw(a.x, b.x, left.x)
				&& betw(a.y, b.y, left.y)
				&& betw(c.x, d.x, left.x)
				&& betw(c.y, d.y, left.y);
		}
	}

	bool intersect_with_circle(double r, pt p, pt q)//всегда одна окружность с центром (0, 0)
	{//аргументы функции это корординаты прямой - две точки, и радиус окружности
		line m(p, q);
		double x0 = -m.a * m.c / (m.a * m.a + m.b * m.b), y0 = -m.b * m.c / (m.a * m.a + m.b * m.b);
		if (m.c * m.c > r * r * (m.a * m.a + m.b * m.b) + EPS)
			//puts("no points");
			return 0;
		else if (abs(m.c * m.c - r * r * (m.a * m.a + m.b * m.b)) < EPS) {
			//puts("1 point");
			cout << x0 << ' ' << y0 << '\n';
			return true;
		}
		else {
			double d = r * r - m.c * m.c / (m.a * m.a + m.b * m.b);
			double mult = sqrt(d / (m.a * m.a + m.b * m.b));
			double ax, ay, bx, by;
			ax = x0 + m.b * mult;
			bx = x0 - m.b * mult;
			ay = y0 - m.a * mult;
			by = y0 + m.a * mult;
			//puts("2 points");
			cout << ax << ' ' << ay << '\n' << bx << ' ' << by << '\n';
			return true;
		}
	}
	bool intersect_double_circle(pt second_circle, double r1, double r2)
	{
		//находим персечение окружностей таким путем что 
		// Ax + By + C = 0,
		// A = -2x2,
		// B = -2y2,
		// C = x2^2 + y2^2 + r1^2 - r2^2.
		double a = -2 * second_circle.x;
		double b = -2 * second_circle.y;
		double c = second_circle.x * second_circle.x + second_circle.y * second_circle.y + r1 * r1 - r2 * r2;
		double x0 = -a * c / (a * a + b * b), y0 = -b * c / (a * a + b * b);
		if (c * c > r1 * r1 * (a * a + b * b) + EPS)
			//puts("no points");
			return 0;
		else if (abs(c * c - r1 * r1 * (a * a + b * b)) < EPS) {
			//puts("1 point");
			cout << x0 << ' ' << y0 << '\n';
			return 1;
		}
		else {
			double d = r1 * r1 - c * c / (a * a + b * b);
			double mult = sqrt(d / (a * a + b * b));
			double ax, ay, bx, by;
			ax = x0 + b * mult;
			bx = x0 - b * mult;
			ay = y0 - a * mult;
			by = y0 + a * mult;
			//puts("2 points");
			cout << ax << ' ' << ay << '\n' << bx << ' ' << by << '\n';
			return 1;
		}
	}
private:
};


class Figures : public Intersection
{
public:
	struct Circle
	{
		double r;
		double P, S;
	};
	struct Quadrilateral
	{
		double P, S;
		pt a, b, c, d;
	};
	struct Triangle
	{
		double P, S;
		pt a, b, c;
	};

	double S, P;
	void CreateCircle(double r)
	{
		Circle t_circle;
		t_circle.r = r;
		S = pi * r * r;
		P = 2 * pi * r;
	}
	void CreateRectangle(pt a, pt b, pt c, pt d)
	{
		Quadrilateral t_quadrilateral;
		t_quadrilateral.a = a;
		t_quadrilateral.b = b;
		t_quadrilateral.c = c;
		t_quadrilateral.d = d;
		//нахождение перимертра прямоугольника по координатам
		P = 2 * (abs(t_quadrilateral.a.x - t_quadrilateral.c.x) + abs(t_quadrilateral.a.y - t_quadrilateral.c.y));
		//нахождение площади по координатам
		S = (t_quadrilateral.c.x - t_quadrilateral.a.x) * (t_quadrilateral.c.y - t_quadrilateral.a.y);
	}
	void CreateTriangle(pt a, pt b, pt c)
	{
		Triangle t_triangle;
		t_triangle.a = a;
		t_triangle.b = b;
		t_triangle.c = c;
		double f1, f2, f3;

		f1 = sqrt((t_triangle.a.x - t_triangle.b.x) * (t_triangle.a.x - t_triangle.b.x) + (t_triangle.a.y - t_triangle.b.y) * (t_triangle.a.y - t_triangle.b.y));
		f2 = sqrt((t_triangle.a.x - t_triangle.c.x) * (t_triangle.a.x - t_triangle.c.x) + (t_triangle.a.y - t_triangle.c.y) * (t_triangle.a.y - t_triangle.c.y));
		f3 = sqrt((t_triangle.c.x - t_triangle.b.x) * (t_triangle.c.x - t_triangle.b.x) + (t_triangle.c.y - t_triangle.b.y) * (t_triangle.c.y - t_triangle.b.y));
		//нахождение периметра прямоугольника по координатам
		P = f1 + f2 + f3;
		t_triangle.P = P;
		//нахождение площади формулой герона
		S = sqrt(t_triangle.P / 2 * (t_triangle.P / 2 - f1) * (t_triangle.P / 2 - f2) * (t_triangle.P / 2 - f3));

		//return t_triangle;
	}
	double area_last_figure()
	{
		return S;
	}
	double perimeter_last_figure()
	{
		return P;
	}
private:
	double pi = 22 / 7;
};