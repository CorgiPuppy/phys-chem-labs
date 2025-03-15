import numpy as np
import matplotlib.pyplot as plt

T = 298.15
R_temp = 8.31
X_1 = [0.100, 0.300, 0.500, 0.700, 0.900]
Y_1 = [0.055, 0.275, 0.589, 0.780, 0.951]
P_bar = [0.2420, 0.2216, 0.2287, 0.2546, 0.2887]
P = [i * 1e5 for i in P_bar]
R = 1.987
A_1, A_2 = 16.65, 15.97
B_1, B_2 = 2940.46, 2696.79
C_1, C_2 = -35.93, -46.16
Tc_1, Tc_2 = 508.10, 536.40
Pc_1, Pc_2 = 47.01, 54.72
w_1, w_2 = 0.31, 0.22

def gamma_calc(X_1, X_2, Y_1, Y_2, P_vp_1, P_vp_2):
	gamma_1 = []
	gamma_2 = []
	for i in range(len(X_1)):
		gamma_1.append((Y_1[i] * P[i]) / (X_1[i] * P_vp_1))
		gamma_2.append((Y_2[i] * P[i]) / (X_2[i] * P_vp_2))
	return (gamma_1, gamma_2)

def gibbs_energy_calc(X_1, X_2, gamma_1, gamma_2):
	g_e = []
	for i in range(0, len(X_1)):
		g_e.append(R * T * (X_1[i] * np.log(gamma_1[i]) +
		X_2[i] * np.log(gamma_2[i])))
	return g_e

def gamma_Wilson_calc(X_1, X_2, lambda_12, lambda_21):
	gamma_1 = []
	gamma_2 = []
	for i in range(len(X_1)):
		gamma_1.append(np.exp(-1 * np.log(X_1[i] + lambda_12 * X_2[i]) + X_2[i] * ((lambda_12 / (X_1[i] + lambda_12 * X_2[i])) - (lambda_21 / (X_2[i] + lambda_21 * X_1[i])))))
		gamma_2.append(np.exp(-1 * np.log(X_2[i] + lambda_21 * X_1[i]) - X_1[i] * ((lambda_12 / (X_1[i] + lambda_12 * X_2[i])) - (lambda_21 / (X_2[i] + lambda_21 * X_1[i])))))
	return (gamma_1, gamma_2)

def gibbs_energy_Wilson_calc(X_1, X_2, lambda_12, lambda_21):
	g_e = []
	for i in range(0, len(X_1)):
		g_e.append(R * T * (-X_1[i] * np.log(X_1[i] + lambda_12 * X_2[i]) - X_2[i] * np.log(X_2[i] + lambda_21 * X_1[i])))
	return g_e

def main():
	X_2 = []
	for i in range(0, len(X_1)):
		X_2.append(1 - X_1[i])

	Y_2 = []
	for i in range(0, len(Y_1)):
		Y_2.append(1 - Y_1[i])

	P_vp_1 = np.exp(A_1 - B_1 / (T + C_1))
	P_vp_2 = np.exp(A_2 - B_2 / (T + C_2))

	gamma_1, gamma_2 = gamma_calc(X_1, X_2, Y_1, Y_2, P_vp_1, P_vp_2)

	g_e = gibbs_energy_calc(X_1, X_2, gamma_1, gamma_2)

	X_1_W = np.arange(0, 1.01, 0.01)
	X_2_W = []
	for i in range(0, len(X_1_W)):
		X_2_W.append(1 - X_1_W[i])

	X_1_W = np.arange(0, 1.01, 0.01)
	X_2_W = []
	for i in range(0, len(X_1_W)):
		X_2_W.append(1 - X_1_W[i])

	alpha_12, alpha_21 = 116, -506
	Z_RA_1 = 0.29056 - 0.08775 * w_1
	Z_RA_2 = 0.29056 - 0.08775 * w_2
	Tr_1 = T / Tc_1
	Tr_2 = T / Tc_2
	V_L_1 = (R_temp * Tc_1 / Pc_1) * Z_RA_1 ** (1 + (1 - Tr_1) ** (2 / 7))
	V_L_2 = (R_temp * Tc_2 / Pc_2) * Z_RA_2 ** (1 + (1 - Tr_2) ** (2 / 7))
	lambda_12 = V_L_2 / V_L_1 * np.exp(-1 * alpha_12 / (R * T))
	lambda_21 = V_L_1 / V_L_2 * np.exp(-1 * alpha_21 / (R * T))
	gamma_1_W, gamma_2_W = gamma_Wilson_calc(X_1_W, X_2_W, lambda_12, lambda_21)

	Y_1_W = [0]
	for i in range(1, len(X_1_W)):
		Y_1_W.append(1 / ((gamma_2_W[i] * X_2_W[i] * P_vp_2) / (gamma_1_W[i] * X_1_W[i] * P_vp_1) + 1))
	Y_2_W = []
	for i in range(0, len(Y_1_W)):
		Y_2_W.append(1 - Y_1_W[i])

	plt.figure
	plt.grid(True)
	plt.title("Ацетон (y-x)")
	plt.plot(X_1_W, Y_1_W, color='blue')
	plt.scatter(X_1, Y_1, color='orange')
	plt.xlabel("X")
	plt.ylabel("Y")
	plt.show()

	plt.figure
	plt.grid(True)
	plt.title("Хлороформ (y-x)")
	plt.plot(X_2_W, Y_2_W, color='blue')
	plt.scatter(X_2, Y_2, color='orange')
	plt.xlabel("X")
	plt.ylabel("Y")
	plt.show()

	P_1_W = []
	for i in range(0, len(X_1_W)):
		P_1_W.append(gamma_1_W[i] * X_1_W[i] * P_vp_1)

	P_2_W = []
	for i in range(0, len(X_2_W)):
		P_2_W.append(gamma_2_W[i] * X_2_W[i] * P_vp_2)

	P_W = []
	for i in range(0, len(P_1_W)):
		P_W.append((P_1_W[i] + P_2_W[i]) / 750.1)

	plt.figure
	plt.grid(True)
	plt.title("Ацетон-Хлороформ (P-y-x)")
	plt.plot(X_1_W, P_W, color="blue")
	plt.plot(Y_1_W, P_W, color="orange")
	plt.scatter(X_1, P_bar, color='orange')
	plt.scatter(Y_1, P_bar, color='orange')
	plt.xlabel("X, Y")
	plt.ylabel("Давление (P), бар")
	plt.show()

	g_e_W = gibbs_energy_Wilson_calc(X_1_W, X_2_W, lambda_12, lambda_21)

main()
