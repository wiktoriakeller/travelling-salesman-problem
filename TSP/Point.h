#pragma once

class Point {
public:
	Point(int x, int y);
	int GetX() const;
	int GetY() const;

private:
	int x;
	int y;
};
