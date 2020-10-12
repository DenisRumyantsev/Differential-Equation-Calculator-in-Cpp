#include <iostream>
#include <vector>
#include <conio.h>
#include <windows.h>

using namespace std;
HWND win = GetConsoleWindow();
HDC hdc = GetDC(win);

HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
COORD coord;

COLORREF color_white = RGB(255, 255, 255);
COLORREF color_red = RGB(255, 0, 0);
COLORREF color_green = RGB(0, 255, 0);
COLORREF color_blue = RGB(0, 0, 255);
COLORREF color_black = RGB(0, 0, 0);
COLORREF color_yellow = RGB(255, 255, 0);
COLORREF color_magenta = RGB(255, 0, 255);
COLORREF color_cyan = RGB(0, 255, 255);
HBRUSH brush_red = CreateSolidBrush(color_red);
HBRUSH brush_black = CreateSolidBrush(color_black);
HPEN pen_red = CreatePen(NULL, 1, color_red);
HPEN pen_green = CreatePen(NULL, 1, color_green);
HPEN pen_blue = CreatePen(NULL, 1, color_blue);
HPEN pen_black = CreatePen(NULL, 1, color_black);

int screen_width = 1366;
int screen_height = 768;

double Fx(double t, double x, double y, double z)
{
	return -10 * (x + y);
}

double Fy(double t, double x, double y, double z)
{
	return -y - 10 * x * z;
}

double Fz(double t, double x, double y, double z)
{
	return 10 * x * y + 4.272;
}

struct TXYZ
{
	vector<double> T, X, Y, Z;
};

TXYZ TXYZ_RK(double t_0, double t_end, double dt, double x_0, double y_0, double z_0)
{
	double
		t = t_0,
		x = x_0,
		y = y_0,
		z = z_0,
		x1, x2, x3,
		y1, y2, y3,
		z1, z2, z3,
		Kx1, Kx2, Kx3, Kx4,
		Ky1, Ky2, Ky3, Ky4,
		Kz1, Kz2, Kz3, Kz4;

	vector<double> T, X, Y, Z;

	T.push_back(t);
	X.push_back(x);
	Y.push_back(y);
	Z.push_back(z);

	while (t < t_end)
	{
		Kx1 = Fx(t, x, y, z);
		Ky1 = Fy(t, x, y, z);
		Kz1 = Fz(t, x, y, z);

		t += dt / 2;
		x1 = x + Kx1 * dt / 2;
		y1 = y + Ky1 * dt / 2;
		z1 = z + Kz1 * dt / 2;
		Kx2 = Fx(t, x1, y1, z1);
		Ky2 = Fy(t, x1, y1, z1);
		Kz2 = Fz(t, x1, y1, z1);

		x2 = x + Kx2 * dt / 2;
		y2 = y + Ky2 * dt / 2;
		z2 = z + Kz2 * dt / 2;
		Kx3 = Fx(t, x2, y2, z2);
		Ky3 = Fy(t, x2, y2, z2);
		Kz3 = Fz(t, x2, y2, z2);

		t += dt / 2;
		x3 = x + Kx3 * dt;
		y3 = y + Ky3 * dt;
		z3 = z + Kz3 * dt;
		Kx4 = Fx(t, x3, y3, z3);
		Ky4 = Fy(t, x3, y3, z3);
		Kz4 = Fz(t, x3, y3, z3);

		x += (Kx1 + 2 * (Kx2 + Kx3) + Kx4) * dt / 6;
		y += (Ky1 + 2 * (Ky2 + Ky3) + Ky4) * dt / 6;
		z += (Kz1 + 2 * (Kz2 + Kz3) + Kz4) * dt / 6;

		T.push_back(t);
		X.push_back(x);
		Y.push_back(y);
		Z.push_back(z);
	}

	TXYZ TXYZ_1;
	TXYZ_1.T = T;
	TXYZ_1.X = X;
	TXYZ_1.Y = Y;
	TXYZ_1.Z = Z;
	return TXYZ_1;
}

TXYZ TXYZ_AB(double t_0, double t_end, double dt, double x_0, double y_0, double z_0)
{
	int i;
	double
		t = t_0 + 3 * dt,
		c_1 = 2.291666666667,
		c_2 = 2.458333333333,
		c_3 = 1.541666666667,
		c_4 = 0.375,
		sum_X, sum_Y, sum_Z,
		x, y, z;
	TXYZ TXYZ_1 = TXYZ_RK(t_0, t, dt, x_0, y_0, z_0);
	vector<double> T, X, Y, Z, F_X, F_Y, F_Z;
	T = TXYZ_1.T;
	X = TXYZ_1.X;
	Y = TXYZ_1.Y;
	Z = TXYZ_1.Z;
	if (T.size() > 4)
	{
		T.erase(T.end() - 1);
		X.erase(X.end() - 1);
		Y.erase(Y.end() - 1);
		Z.erase(Z.end() - 1);
	}
	for (i = 0; i < 4; i++)
	{
		F_X.push_back(Fx(T[i], X[i], Y[i], Z[i]));
		F_Y.push_back(Fy(T[i], X[i], Y[i], Z[i]));
		F_Z.push_back(Fz(T[i], X[i], Y[i], Z[i]));
	}
	x = X[3];
	y = Y[3];
	z = Z[3];
	while (t < t_end)
	{
		sum_X
			= c_1 * F_X[F_X.size() - 1]
			- c_2 * F_X[F_X.size() - 2]
			+ c_3 * F_X[F_X.size() - 3]
			- c_4 * F_X[F_X.size() - 4];

		sum_Y
			= c_1 * F_Y[F_Y.size() - 1]
			- c_2 * F_Y[F_Y.size() - 2]
			+ c_3 * F_Y[F_Y.size() - 3]
			- c_4 * F_Y[F_Y.size() - 4];

		sum_Z
			= c_1 * F_Z[F_Z.size() - 1]
			- c_2 * F_Z[F_Z.size() - 2]
			+ c_3 * F_Z[F_Z.size() - 3]
			- c_4 * F_Z[F_Z.size() - 4];

		t += dt;
		x += sum_X * dt;
		y += sum_Y * dt;
		z += sum_Z * dt;
		T.push_back(t);
		X.push_back(x);
		Y.push_back(y);
		Z.push_back(z);
		F_X.push_back(Fx(t, x, y, z));
		F_Y.push_back(Fy(t, x, y, z));
		F_Z.push_back(Fz(t, x, y, z));
	}
	TXYZ_1.T = T;
	TXYZ_1.X = X;
	TXYZ_1.Y = Y;
	TXYZ_1.Z = Z;
	return TXYZ_1;
}

void grid(int Sx, int Sy, int l)
{
	int i;
	SelectObject(hdc, pen_black);
	SelectObject(hdc, brush_black);
	Rectangle(hdc, 0, 0, screen_width, screen_height);
	SelectObject(hdc, brush_red);
	SelectObject(hdc, pen_blue);
	for (i = -500; i < 500; i++)
	{
		MoveToEx(hdc, Sx + l * i, 0, NULL);
		LineTo(hdc, Sx + l * i, screen_height);
		MoveToEx(hdc, 0, Sy + l * i, NULL);
		LineTo(hdc, screen_width, Sy + l * i);
	}
	SelectObject(hdc, pen_green);
	MoveToEx(hdc, Sx, 0, NULL);
	LineTo(hdc, Sx, screen_height);
	MoveToEx(hdc, 0, Sy, NULL);
	LineTo(hdc, screen_width, Sy);
	SelectObject(hdc, pen_red);
}

void step(double x, double y, int Sx, int Sy, int l, int s)
{
	Rectangle(hdc, l * x + Sx - s, -l * y + Sy - s, l * x + Sx + s, -l * y + Sy + s);
}

TXYZ system_AB(double x_0, double y_0, double z_0,
	double t_0, double t_end, double eps)
{
	int i, i2;
	double
		dt = eps, delta,
		max_delta = eps + 1;

	TXYZ TXYZ_1, TXYZ_2;
	TXYZ_1 = TXYZ_AB(t_0, t_end, dt, x_0, y_0, z_0);

	while (max_delta > eps)
	{
		max_delta = 0;
		dt /= 2;
		TXYZ_2 = TXYZ_AB(t_0, t_end, dt, x_0, y_0, z_0);

		for (i = 0; i < TXYZ_1.T.size() - 1; i++)
		{
			i2 = 2 * i;
			delta = abs(TXYZ_1.X[i] - TXYZ_2.X[i2]);
			if (delta > max_delta)
				max_delta = delta;
			delta = abs(TXYZ_1.Y[i] - TXYZ_2.Y[i2]);
			if (delta > max_delta)
				max_delta = delta;
			delta = abs(TXYZ_1.Z[i] - TXYZ_2.Z[i2]);
			if (delta > max_delta)
				max_delta = delta;
		}
		TXYZ_1 = TXYZ_2;
	}
	return TXYZ_1;
}

void graph_XY(int Sx, int Sy, int l, int s,
	vector<double> T, vector<double> X, vector<double> Y, double st)
{
	int i;
	double t_st = T[0];
	grid(Sx, Sy, l);
	for (i = 0; i < T.size(); i++)
	{
		SetPixel(hdc, l * X[i] + Sx, -l * Y[i] + Sy, color_white);
		if (T[i] >= t_st)
		{
			step(X[i], Y[i], Sx, Sy, l, s);
			t_st += st;
		}
	}
}

void graph_tX(int Sx, int Sy, int l,
	vector<double> T, vector<double> X)
{
	int i;
	grid(Sx, Sy, l);
	for (i = 0; i < T.size(); i++)
		SetPixel(hdc, l * T[i] + Sx, -l * X[i] + Sy, color_white);
}

void graph_t_(int Sx, int Sy, int l,
	vector<double> T, vector<double> X, vector<double> Y, vector<double> Z)
{
	int i;
	grid(Sx, Sy, l);
	for (i = 0; i < T.size(); i++)
	{
		SetPixel(hdc, l * T[i] + Sx, -l * X[i] + Sy, color_yellow);
		SetPixel(hdc, l * T[i] + Sx, -l * Y[i] + Sy, color_magenta);
		SetPixel(hdc, l * T[i] + Sx, -l * Z[i] + Sy, color_cyan);
	}
}

void table_TXYZ(vector<double> T, vector<double> X, vector<double> Y, vector<double> Z, double st)
{
	int i;
	double t_st = T[0];
	cout << endl << "T\tX\tY\tZ" << endl;
	for (i = 0; i < T.size(); i++)
		if (T[i] >= t_st)
		{
			cout << t_st << "\t" << X[i] << "\t" << Y[i] << "\t" << Z[i] << endl;
			t_st += st;
		}
	SetConsoleCursorPosition(handle, coord);
}

void result(TXYZ TXYZ_1, double st)
{
	vector<double> T, X, Y, Z;

	T = TXYZ_1.T;
	X = TXYZ_1.X;
	Y = TXYZ_1.Y;
	Z = TXYZ_1.Z;

	int
		i, Sx, Sy, l, s = 3,
		Sx_XY = 500, Sy_XY = 400, l_XY = 100,
		Sx_XZ = 500, Sy_XZ = 400, l_XZ = 100,
		Sx_YZ = 500, Sy_YZ = 400, l_YZ = 100,
		Sx_tX = 10, Sy_tX = 400, l_tX = 40,
		Sx_tY = 10, Sy_tY = 400, l_tY = 40,
		Sx_tZ = 10, Sy_tZ = 400, l_tZ = 40,
		Sx_t_ = 10, Sy_t_ = 400, l_t_ = 40;

	double t_st, t_0 = T[0];

	coord.X = 0;
	coord.Y = 0;

	cout << fixed;
	cout.precision(1);

	while (true)
		if (_kbhit())
			switch (_getch())
			{
			case (49):
				graph_XY(Sx_XY, Sy_XY, l_XY, s, T, X, Y, st);
				break;

			case (50):
				graph_XY(Sx_XZ, Sy_XZ, l_XZ, s, T, X, Z, st);
				break;

			case (51):
				graph_XY(Sx_YZ, Sy_YZ, l_YZ, s, T, Y, Z, st);
				break;

			case (52):
				graph_tX(Sx_tX, Sy_tX, l_tX, T, X);
				break;

			case (53):
				graph_tX(Sx_tY, Sy_tY, l_tY, T, Y);
				break;

			case (54):
				graph_tX(Sx_tZ, Sy_tZ, l_tZ, T, Z);
				break;

			case(55):
				graph_t_(Sx_t_, Sy_t_, l_t_, T, X, Y, Z);
				break;

			case(56):
				table_TXYZ(T, X, Y, Z, st);
				break;

			default:
				break;
			}
}

int main()
{
	// result(system_AB(x_0, y_0, z_0, t_0, t_end, eps), st);
	result(system_AB(1, 1, 1, 0, 30, 0.1), 5);
	return 0;
}