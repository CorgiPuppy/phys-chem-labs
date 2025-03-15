const std = @import("std");
const print = std.debug.print; // для вывода в консоль
const math = std.math; // для работы с библиотекой math

const R: f32 = 8.31; // универсальная газовая постоянная
const T1_1: f32 = 298.0; // начальная температура (задание 1_1)
const T2_1: f32 = 300.0; // температура повышенная (задание 1_1)
const T_2: f32 = 800.0; // температура для смеси (задание 2)
const p_N_2: f32 = 20.0e4; // парциальное давление азота (задание 2)
const p_H_2: f32 = 20.0e4; // парциальное давление водорода (задание 2)
const p_NH_3: f32 = 10.0e4; // парциальное давление аммиака (задание 2)
const v_COCl_2: f32 = 1.0; // стехиометрический коэффициент фосгена (задание 1)
const v_CO: f32 = 1.0; // стехиометрический коэффициент монооксида (задание 1)
const v_Cl_2: f32 = 1.0; // стехиометрический коэффициент хлора (задание 1)
const v_N_2: f32 = 1.0 / 2.0; // стехиометрический коэффициент азота (задание 2)
const v_H_2: f32 = 3.0 / 2.0; // стехиометрический коэффициент водорода (задание 2)
const v_NH_3: f32 = 1.0; // стехиометрический коэффициент аммиака (задание 2)

const NUM_COEF: i32 = 7; // количество коэффициентов в полиноме NASA

// определение коэффициентов для полиномов NASA (задание 1)
const a_COCl_2 = [NUM_COEF]f32{ 1.70787910E+00, 2.89369464E-02, -4.93289116E-05, 4.16910139E-08, -1.37057391E-11, -2.78350932E+04, 1.76202114E+01 };
const a_CO = [NUM_COEF]f32{ 3.57953347E+00, -6.10353680E-04, 1.01681433E-06, 9.07005884E-10, -9.04424499E-13, -1.43440860E+04, 3.50840928E+00 };
const a_Cl_2 = [NUM_COEF]f32{ 2.73638114E+00, 7.83525699E-03, -1.45104963E-05, 1.25730834E-08, -4.13247143E-12, -1.05880114E+03, 9.44557148E+00 };

// определение коэффициентов для полиномов NASA (задание 2)
const a_N_2 = [NUM_COEF]f32{ 3.53100528E+00, -1.23660988E-04, -5.02999433E-07, 2.43530612E-09, -1.40881235E-12, -1.04697628E+03, 2.96747038E+00 };
const a_H_2 = [NUM_COEF]f32{ 0.23443029E+01, 0.79804248E-02, -0.19477917E-04, 0.20156967E-07, -0.73760289E-11, -0.91792413E+03, 0.68300218E+00 };
const a_NH_3 = [NUM_COEF]f32{ 4.46075151E+00, -5.68781763E-03, 2.11411484E-05, -2.02849980E-08, 6.89500555E-12, -6.70753514E+03, -1.34450793E+00 };

// структура для вывода мольных долей
const MolarFrac = struct {
	title: []const u8,
	meaning: f32,
};

pub fn main() !void {
	print("{s}", .{"Lab: 5. Variant: 7\n"});
	print("{s}", .{"Task_1\n"});
	print("{s}", .{"\tInitial chemical reaction: COCl2 -> CO + Cl2\n"});

	// объявление стандартных энтальпий образований веществ
	var H_COCl_2: f32 = 0.0;
	var H_CO: f32 = 0.0;
	var H_Cl_2: f32 = 0.0;
	
	// поиск стандартных энтальпий образований веществ при T = 298K
	H_COCl_2 = enthalpy(T1_1, &a_COCl_2);
	H_CO = enthalpy(T1_1, &a_CO);
	H_Cl_2 = enthalpy(T1_1, &a_Cl_2);

	var H: f32 = 0.0; // объявление стандартного теплового эффекта
	H = (H_CO + H_Cl_2) - H_COCl_2; // поиск стандартного теплового эффекта

	var K_p: f32 = 0.0; // объявление константы равновесия
	K_p = equilibrium_constant_VantGoff(H, T1_1); // поиск константы равновесия при T = 298K

	print("\tTask_1_1:\n\t\tdH = {any} kJ/mol (T = {any} K)\n\t\tK_p = {any}\n", .{ H / 1000, T1_1, K_p }); // вывод задания 1_1

	// поиск стандартных энтальпий образований веществ при T = 300K
	H_COCl_2 = enthalpy(T2_1, &a_COCl_2);
	H_CO = enthalpy(T2_1, &a_CO);
	H_Cl_2 = enthalpy(T2_1, &a_Cl_2);

	H = (H_CO + H_Cl_2) - H_COCl_2; // поиск стандартного теплового эффекта

	K_p = equilibrium_constant_VantGoff(H, T2_1); // поиск константы равновесия при T = 300K

	print("\t\tdH = {any} kJ/mol (T = {any} K)\n\t\tK_p = {any}\n", .{ H / 1000, T2_1, K_p }); // вывод задания 1_1

	var dv: i32 = 0; // объявление разности между стехиометрическими коэффициентами
	dv = (v_CO + v_Cl_2) - v_COCl_2; // поиск разности между стехиометрическими коэффициентами

	print("\tTask_1_2:\n\t\tdv = {}\n", .{dv}); // вывод задания_1_2
	print("\tTask_1_3:\n\t\tdv = {}\n", .{dv}); // вывод задания_1_3

	print("{s}", .{"Task_2\n"});
	print("{s}", .{"\tInitial chemical reaction: 1/2N2 + 3/2H2 -> NH3\n"});

	// объявление стандартных энтальпий образований веществ
	var H_N_2: f32 = 0.0;
	var H_H_2: f32 = 0.0;
	var H_NH_3: f32 = 0.0;

	// поиск стандартных энтальпий образований веществ при T = 800K
	H_N_2 = enthalpy(T_2, &a_N_2);
	H_H_2 = enthalpy(T_2, &a_H_2);
	H_NH_3 = enthalpy(T_2, &a_NH_3);

	// объявление стандартных энтропий образований веществ
	var S_N_2: f32 = 0.0;
	var S_H_2: f32 = 0.0;
	var S_NH_3: f32 = 0.0;

	// поиск стандартных энтропий образований веществ при T = 800K
	S_N_2 = entropia(T_2, &a_N_2);
	S_H_2 = entropia(T_2, &a_H_2);
	S_NH_3 = entropia(T_2, &a_NH_3);

	var G: f32 = 0.0; // объявление стандартного изменения энергии Гиббса
	G = 1 * H_NH_3 - (1 / 2 * H_N_2 + 3 / 2 * H_H_2) - T_2 * (1 * S_NH_3 - (1 / 2 * S_N_2 + 3 / 2 * S_H_2)); // поиск стандартного изменения энергии Гиббса
	print("\tTask_2_1:\n\t\tdG = {any} kJ/mol\n", .{G / 1000}); // вывод задания_2_1

	// объявление мольных долей веществ
	var x_N_2: f32 = 0.0;
	var x_H_2: f32 = 0.0;
	var x_NH_3: f32 = 0.0;

	// поиск мольных долей веществ
	x_N_2 = molar_fraction(p_N_2, p_N_2 + p_H_2 + p_NH_3);
	x_H_2 = molar_fraction(p_H_2, p_N_2 + p_H_2 + p_NH_3);
	x_NH_3 = molar_fraction(p_NH_3, p_N_2 + p_H_2 + p_NH_3);

	// вывод задания 2_2
	print("\t{s}", .{"Task_2_2:\n"});
	var alloc = std.heap.ArenaAllocator.init(std.heap.page_allocator);
	defer alloc.deinit();
	const allocator = alloc.allocator();
	var list = std.MultiArrayList(MolarFrac){};
	try list.append(allocator, .{ .title = "X_N_2", .meaning = x_N_2 });
	try list.append(allocator, .{ .title = "X_H_2", .meaning = x_H_2 });
	try list.append(allocator, .{ .title = "X_NH_3", .meaning = x_NH_3 });
	for (list.items(.title), list.items(.meaning)) |title, meaning| {
		print("\t\t{s} = {any}\n", .{ title, meaning });
	}
	
	K_p = equilibrium_constant_main(G, T_2); // поиск константы равновесия
	print("\tTask_2_3\n\t\tK_a = {}\n", .{K_p}); // вывод задания_2_3

	// вычисление корней уравнения
	const x_1 = math.pow(f32, 4 / ((p_N_2 + p_H_2 + p_NH_3) * K_p * math.pow(f32, 3.0, 3.0 / 2.0) + 4.0), 1.0 / 2.0);
	const x_2 = -math.pow(f32, 4 / ((p_N_2 + p_H_2 + p_NH_3) * K_p * math.pow(f32, 3.0, 3.0 / 2.0) + 4.0), 1.0 / 2.0);

	// вычисление равновесных концентраций
	const C_N_2 = equilibrium_concentration(x_1, (p_N_2 + p_H_2 + p_NH_3), T_2, v_N_2);
	const C_H_2 = equilibrium_concentration(x_1, (p_N_2 + p_H_2 + p_NH_3), T_2, v_H_2);
	const C_NH_3 = equilibrium_concentration(x_1, (p_N_2 + p_H_2 + p_NH_3), T_2, v_NH_3);

	print("{s}x_1 = {}\n\t\tx_2 = {}\n\t\tC_N_2 = {} mol/l\n\t\tC_H_2 = {} mol/l\n\t\tC_NH_3 = {} mol/l\n", .{ "\tTask_2_4\n\t\t", x_1, x_2, C_N_2, C_H_2, C_NH_3 }); // вывод задания_2_4

	// вычисление равновесной степени превращения
	const alpha = equilibrium_degree_transformation(x_1, v_N_2, v_N_2);
	print("\tTask_2_5\n\t\talpha = {}\n", .{alpha}); // вывод задания_2_5
}

fn enthalpy(T: f32, a: *const [NUM_COEF]f32) f32 {
	return (a[0] + a[1] / 2 * T + a[2] / 3 * math.pow(f32, T, 2) + a[3] / 4 * math.pow(f32, T, 3) + a[4] / 5 * math.pow(f32, T, 4) + a[5] / T) * R * T;
}

fn equilibrium_constant_VantGoff(H: f32, T: f32) f32 {
	return math.exp(H / (R * math.pow(f32, T, 2)));
}

fn molar_fraction(p_i: f32, p_total: f32) f32 {
	return p_i / p_total;
}

fn entropia(T: f32, a: *const [NUM_COEF]f32) f32 {
	return (a[0] * math.log(f32, math.e, T) + a[1] * T + a[2] / 2 * math.pow(f32, T, 2) + a[3] / 3 * math.pow(f32, T, 3) + a[4] / 4 * math.pow(f32, T, 4) + a[6]) * R;
}

fn equilibrium_constant_main(G: f32, T: f32) f32 {
	return math.exp(-(G / (R * T)));
}

fn equilibrium_concentration(x: f32, p: f32, T: f32, v: f32) f32 {
	return v * p * x / ((1 + x) * R * T);
}

fn equilibrium_degree_transformation(x: f32, n: f32, n_0: f32) f32 {
	return n * x / n_0;
}
