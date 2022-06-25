#include <iostream>
#include <math.h>
#include "Header.h"

using namespace std;


int main()
{ 
	//Вывод площади и периметра круга
	Figures Circle;
	Circle.CreateCircle(3);
	cout << Circle.area_last_figure() << " " << Circle.perimeter_last_figure() << "\n";
	//Вывод площади и периметра Прямоугольника
	Figures Rect;//координаты задаются по часовой стрелке либо же против часовой
	Figures::pt a ,b, c, d;
	a.x = 0;
	a.y = 0;

	b.x = 0;
	b.y = 2;

	c.x = 4;
	c.y = 2;

	d.x = 4;
	d.y = 0;
	Rect.CreateRectangle(a,b,c,d);
	cout << Rect.area_last_figure() << " " << Rect.perimeter_last_figure() << "\n";
	//Вывод площади и периметра треугольника
	Figures Triangle;
	Triangle.CreateTriangle(a, b, c);
	cout << Triangle.area_last_figure() << " " << Triangle.perimeter_last_figure() << "\n";
	
	Figures::pt p, q, p1, q1, left, right;
	p.x = 0;
	p.y = 0;

	q.x = 1;
	q.y = 1;

	p1.x = 0;
	p1.y = 1;

	q1.x = 4;
	q1.y = 1;
	//проверка пересечения двух отрезков
	Intersection example;
	cout<<example.intersect_segments(p, q, p1, q1, left, right)<<"\n";

	//точки пересечения прямой и окружности 
	example.intersect_with_circle(2, p, q);

	//точки пересечения окружности (представленной в виде прямой) и окружности 
	example.intersect_double_circle(q, 1, 1);
	return 0;
}