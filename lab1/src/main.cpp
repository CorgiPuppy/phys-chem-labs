#include <iostream> // для вывода
#include <math.h> // для использования натурального логарифма и возведения числа в степень
#include <fstream> // для вывода в файл

using namespace std;

// defining constant
#define R 8.31 // универсальная газовая постоянная
#define RANGE 11 // количество рассматриваемых температур (задание 1)
#define NUM_COEF1 6 // количество коэффициентов (задание 1)
#define NUM_COEF2 7 // количество коэффициентов (задание 2)
#define P 101325 // нормальное давление (задание 2_2)
#define TEMP 298 // температура (задание 2_2)
#define V_A 2 // объём первого вещества (задание 2_2)
#define V_B 3 // объём второго вещества (задание 2_2)

// defining of temperature limits
enum Temp { T1 = 500, T2 = 1000 };

// defining of coefficients
constexpr double a_C2H4[NUM_COEF1] = {3.952920063E00, -7.57051373E-03,
5.70989993E-05, -6.91588352E-08, 2.69884190E-11, 5.08977598E+03};
constexpr double a_C2H6[NUM_COEF1] = {4.29142572E+00, -5.50154901E-03,
5.99438458E-05, -7.08466469E-08, 2.68685836E-11, -1.15222056E+04};
constexpr double a_C2H6S[NUM_COEF2] = {0.19139966E+01, 0.29420442E-01,
-0.24128528E-04, 0.15495718E-07, -0.50061422E-11, -0.62072425E+04,
0.15648303E+02};

void enthalpy(int [], const double [], double []);
double entropia(const double []);
double entropia_mix();

int main() {
	// объявление массива температур
	int T[RANGE] = {T1};
	for (int i = 0; i < RANGE; i++)
		T[i + 1] = T[i] + (T2 - T1) / (RANGE - 1);
	
	// объявление стандартных энтальпий образований веществ
	double H_C2H4[RANGE] = {0};
	double H_H2;
	double H_C2H6[RANGE] = {0};
	
	enthalpy(T, a_C2H4, H_C2H4);
	H_H2 = 0.0;
	enthalpy(T, a_C2H6, H_C2H6);
	
	// поиск стандартного теплового эффекта
	double H[RANGE] = {0};
	for(int i = 0; i < RANGE; i++)
		H[i] = H_C2H6[i] - H_H2 - H_C2H4[i];
	
	printf("Task_1:\n");
	for(int i = 0; i < RANGE; i++)
		printf("\tdH = %.2lf kJ/m (T = %d K)\n", H[i] / 1000, T[i]);
	
	// вывод в файл для отображения на графике
	ofstream output ("../files/output.txt");
	if (output.is_open()) {
		for (int i = 0; i < RANGE; i++)
			output << T[i] << " " << H[i]/1000 << "\n";
		
		output.close();
	}
	else 
		cout << "Unable to open file";
	
	// поиск энтропии вещества
	double S = entropia(a_C2H6S);
	
	// поиск энтропии при смешивании двух веществ
	double dS_mix = entropia_mix();
	printf("Task_2_1:\n\tS = %.5lf J/m/K\nTask_2_2:\n\tS_mix = %.5lfJ/K\n", S, dS_mix);
	
	return 0;
}

// функция для поиска стандартной энтальпии образования вещества с помощью полинома NASA
void enthalpy(int T[], const double a[], double DH[]) {
	for (int i = 0; i < RANGE; i++) {
		DH[i] = (a[0] +
				 a[1]/2 * T[i] +
				 a[2]/3 * pow(T[i], 2) +
				 a[3]/4 * pow(T[i], 3) +
				 a[4]/5 * pow(T[i], 4) +
				 a[5] / T[i]) * R * T[i];
	}
}

// функция для поиска стандартной энтропии вещества с помощью полинома NASA
double entropia(const double a[]) {
	double Temper = 250;
	return (a[0] * std::log(Temper) +
			a[1] * Temper +
			a[2]/2 * pow(Temper, 2) +
			a[3]/3 * pow(Temper, 3) +
			a[4]/4 * pow(Temper, 4) +
			a[6]) * R;
}

// функция для поиска энтропии при смешивании двух веществ
double entropia_mix() {
	double n_a = P * V_A / (R * TEMP);
	double n_b = P * V_B / (R * TEMP);
	double x_a = n_a / (n_a + n_b);
	double x_b = n_b / (n_a + n_b);

	return -(n_a + n_b) * R * (x_a * std::log(x_a) + x_b * std::log(x_b)) * pow(10, -4);
}
