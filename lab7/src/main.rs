use std::{fs::File, io::prelude::*};


// объявление структуры исходных данных
#[derive(Debug)]
struct InitialData {
    universal_gas_equation_system: f32, // универсальная газовая постоянная
    temperature_1: f32, // температура газовой смеси (задание 1)
    time: Vec<i32>, // время (задание 1)
    pressure: Vec<f32>, // давление газовой смеси (задание 1)
    time_1: i32, // момент времени после начала реакции (задание 1)
    temperature_2: f32, // температура (задание 2)
    k_1: f32, // константа скорости первой реакции (задание 2)
    k_2: f32, // константа скорости второй реакции (задание 2)
    c_c8h18: f32, // начальная концентрация вещества C8H18 (задание 2)
    c_ic8h18: f32, // начальная концентрация вещества i-C8H18 (задание 2)
    c_c4h10: f32, // начальная концентрация вещества C4H10 (задание 2)
    c_c4h8: f32, // начальная концентрация вещества C4H8 (задание 2)
    dt: f32, // шаг в методе Рунге-Кутты 4-го порядка (задание 2)
    number_of_iterations: Vec<usize>, // количество итераций
}


impl InitialData {
    // вычисление давления вещества A
    fn pressure_calc(&self) -> Vec<f32> {
        let mut temporary: Vec<f32> = vec![];
        for i in 0..self.pressure.len() {
            temporary.push(self.pressure[0] - (2.0 * self.pressure[i] - 2.0 * self.pressure[0]));
        }
        return temporary;
    }


    // вычисление концентрации вещества A
    fn concentration_calc(&self, pressure_a: Vec<f32>) -> Vec<f32> {
        let mut temporary: Vec<f32> = vec![];
        for i in 0..pressure_a.len() {
            temporary.push(pressure_a[i] / (self.temperature_1 * self.universal_gas_equation_system));
        }
        return temporary;
    }


    // вычисление анаморфозы нулевого порядка
    fn anamorphosis0(&self, concentration_a: Vec<f32>) -> Vec<f32> {
        let mut temporary: Vec<f32> = vec![];
        for i in 0..concentration_a.len() {
            temporary.push(concentration_a[0] - concentration_a[i]);
        }
        return temporary;
    }


    // вычисление анаморфозы первого порядка
    fn anamorphosis1(&self, concentration_a: Vec<f32>) -> Vec<f32> {
        let mut temporary: Vec<f32> = vec![];
        for i in 0..concentration_a.len() {
            temporary.push((concentration_a[0] / concentration_a[i]).ln());
        }
        return temporary;
    }


    // вычисление анаморфозы второго порядка
    fn anamorphosis2(&self, concentration_a: Vec<f32>) -> Vec<f32> {
        let mut temporary: Vec<f32> = vec![];
        for i in 0..concentration_a.len() {
            temporary.push(1.0 / concentration_a[i] - 1.0 / concentration_a[0]);
        }
        return temporary;
    }


    // вычисление константы скорости аналитическим методом для 0-го порядка
    fn rate_constant0_calc(&self, concentration_a: Vec<f32>) -> Vec<f32> {
        let mut temporary: Vec<f32> = vec![];
        for i in 0..concentration_a.len() {
            temporary.push((1 as f32 / self.time[i] as f32) * (concentration_a[0] - concentration_a[i]));
        }
        return temporary;
    }


    // вычисление константы скорости аналитическим методом для 1-го порядка
    fn rate_constant1_calc(&self, concentration_a: Vec<f32>) -> Vec<f32> {
        let mut temporary: Vec<f32> = vec![];
        for i in 0..concentration_a.len() {
            temporary.push((1 as f32 / self.time[i] as f32) * (concentration_a[0] / concentration_a[i]).ln());
        }
        return temporary;
    }


    // вычисление константы скорости аналитическим методом для 2-го порядка
    fn rate_constant2_calc(&self, concentration_a: Vec<f32>) -> Vec<f32> {
        let mut temporary: Vec<f32> = vec![];
        for i in 0..concentration_a.len() {
            temporary.push((1 as f32 / self.time[i] as f32) * (1.0 / concentration_a[i] - 1.0 / concentration_a[0]));
        }
        return temporary;
    }


    // вычисление времени полупревращения, концентрации и степени превращения в момент времени t1
    fn time_moment(&self, moment: f32, rate_constant: f32, concentration0: f32) -> (f32, f32, f32) {
        let mut half_turn_time: f32 = 0.0;
        let mut concentration_time_moment: f32 = 0.0;
        let mut tranformation_degree: f32 = 0.0;
        half_turn_time = 1.0 / (rate_constant * concentration0);
        concentration_time_moment = 1.0 / (rate_constant * moment + 1.0 / concentration0);
        tranformation_degree = rate_constant * moment * concentration_time_moment;
        return (half_turn_time, concentration_time_moment, tranformation_degree);
    }


   
    // решение системы дифференциальных уравнений методом Рунге-Кутты 7-го порядка
    fn equation_system(&self, mut c_c8h18_mod: Vec<f32>, mut c_ic8h18_mod: Vec<f32>, mut c_c4h10_mod: Vec<f32>, mut c_c4h8_mod: Vec<f32>) -> (Vec<f32>, Vec<f32>, Vec<f32>, Vec<f32>) {
        for i in 0..self.number_of_iterations.len() {
            let h = self.dt;
            let k1_c8h18 = -self.k_1 * c_c8h18_mod[i];
            let k1_ic8h18 = self.k_1 * c_c8h18_mod[i] - self.k_2 * c_ic8h18_mod[i];
            let k1_c4h10 = self.k_2 * c_ic8h18_mod[i];
            let k1_c4h8 = self.k_2 * c_ic8h18_mod[i];
   
            let k2_c8h18 = -self.k_1 * (c_c8h18_mod[i] + h * (1 as f32 / 12 as f32) * k1_c8h18);
            let k2_ic8h18 = self.k_1 * (c_c8h18_mod[i] + h * (1 as f32 / 12 as f32) * k1_c8h18) - self.k_2 * (c_ic8h18_mod[i] + h * (1 as f32 / 12 as f32) * k1_ic8h18);
            let k2_c4h10 = self.k_2 * (c_ic8h18_mod[i] + h * (1 as f32 / 12 as f32) * k1_ic8h18);
            let k2_c4h8 = self.k_2 * (c_ic8h18_mod[i] + h * (1 as f32 / 12 as f32) * k1_ic8h18);
   
            let k3_c8h18 = -self.k_1 * (c_c8h18_mod[i] + h * (-10.0 * k1_c8h18 + 11.0 * k2_c8h18)/12.0);
            let k3_ic8h18 = self.k_1 * (c_c8h18_mod[i] + h * (-10.0 * k1_c8h18 + 11.0 * k2_c8h18)/12.0) - self.k_2 * (c_ic8h18_mod[i] + h * (-10.0 * k1_ic8h18 + 11.0 * k2_ic8h18)/12.0);
            let k3_c4h10 = self.k_2 * (c_ic8h18_mod[i] + h * (-10.0 * k1_ic8h18 + 11.0 * k2_ic8h18)/12.0);
            let k3_c4h8 = self.k_2 * (c_ic8h18_mod[i] + h * (-10.0 * k1_ic8h18 + 11.0 * k2_ic8h18)/12.0);
   
            let k4_c8h18 = -self.k_1 * (c_c8h18_mod[i] + h * (2 as f32 / 12 as f32) * k3_c8h18);
            let k4_ic8h18 = self.k_1 * (c_c8h18_mod[i] + h * (2 as f32 / 12 as f32) * k3_c8h18) - self.k_2 * (c_ic8h18_mod[i] + h * (2 as f32 / 12 as f32) * k3_ic8h18);
            let k4_c4h10 = self.k_2 * (c_ic8h18_mod[i] + h * (2 as f32 / 12 as f32) * k3_ic8h18);
            let k4_c4h8 = self.k_2 * (c_ic8h18_mod[i] + h * (2 as f32 / 12 as f32) * k3_ic8h18);


            let k5_c8h18 = -self.k_1 * (c_c8h18_mod[i] + h * (157.0 * k1_c8h18 - 318.0 * k1_c8h18 + 4.0 * k1_c8h18 + 160.0 * k1_c8h18)/9.0);
            let k5_ic8h18 = self.k_1 * (c_c8h18_mod[i] + h * (157.0 * k1_c8h18 - 318.0 * k1_c8h18 + 4.0 * k1_c8h18 + 160.0 * k1_c8h18)/9.0) - self.k_2 * (c_ic8h18_mod[i] + h * (157.0 * k1_ic8h18 - 318.0 * k1_ic8h18 + 4.0 * k1_ic8h18 + 160.0 * k1_ic8h18)/9.0);
            let k5_c4h10 = self.k_2 * (c_ic8h18_mod[i] + h * (157.0 * k1_ic8h18 - 318.0 * k1_ic8h18 + 4.0 * k1_ic8h18 + 160.0 * k1_ic8h18)/9.0);
            let k5_c4h8 = self.k_2 * (c_ic8h18_mod[i] + h * (157.0 * k1_ic8h18 - 318.0 * k1_ic8h18 + 4.0 * k1_ic8h18 + 160.0 * k1_ic8h18)/9.0);


            let k6_c8h18 = (-self.k_1 * (c_c8h18_mod[i] + h * (-322.0 * k1_c8h18 + 199.0 * k1_c8h18 + 108.0 * k1_c8h18 - 131.0 * k1_c8h18)/30.0)) as f32;
            let k6_ic8h18 = self.k_1 * (c_c8h18_mod[i] + h * (-322.0 * k1_c8h18 + 199.0 * k1_c8h18 + 108.0 * k1_c8h18 - 131.0 * k1_c8h18)/30.0) - self.k_2 * (c_ic8h18_mod[i] + h * (-322.0 * k1_ic8h18 + 199.0 * k1_ic8h18 + 108.0 * k1_ic8h18 - 131.0 * k1_ic8h18)/30.0);
            let k6_c4h10 = self.k_2 * (c_ic8h18_mod[i] + h * (-322.0 * k1_ic8h18 + 199.0 * k1_ic8h18 + 108.0 * k1_ic8h18 - 131.0 * k1_ic8h18)/30.0);
            let k6_c4h8 = self.k_2 * (c_ic8h18_mod[i] + h * (-322.0 * k1_ic8h18 + 199.0 * k1_ic8h18 + 108.0 * k1_ic8h18 - 131.0 * k1_ic8h18)/30.0);


            let k7_c8h18 = -self.k_1 * (c_c8h18_mod[i] + h * ((3158 as f32 / 45 as f32) * k1_c8h18 - (638 as f32 / 6 as f32) * k2_c8h18 - (23 as f32 / 2 as f32) * k3_c8h18 + (157 as f32 / 3 as f32) * k4_c8h18 + (157 as f32 / 45 as f32) * k6_c8h18));
            let k7_ic8h18 = self.k_1 * (c_c8h18_mod[i] + h * ((3158 as f32 / 45 as f32) * k1_c8h18 - (638 as f32 / 6 as f32) * k2_c8h18 - (23 as f32 / 2 as f32) * k3_c8h18 + (157 as f32 / 3 as f32) * k4_c8h18 + (157 as f32 / 45 as f32) * k6_c8h18)) - self.k_2 * (c_ic8h18_mod[i] + h * ((3158 as f32 / 45 as f32) * k1_ic8h18 - (638 as f32 / 6 as f32) * k2_ic8h18 - (23 as f32 / 2 as f32) * k3_ic8h18 + (157 as f32 / 3 as f32) * k4_ic8h18 + (157 as f32 / 45 as f32) * k6_ic8h18));
            let k7_c4h10 = self.k_2 * (c_ic8h18_mod[i] + h * ((3158 as f32 / 45 as f32) * k1_ic8h18 - (638 as f32 / 6 as f32) * k2_ic8h18 - (23 as f32 / 2 as f32) * k3_ic8h18 + (157 as f32 / 3 as f32) * k4_ic8h18 + (157 as f32 / 45 as f32) * k6_ic8h18));
            let k7_c4h8 = self.k_2 * (c_ic8h18_mod[i] + h * ((3158 as f32 / 45 as f32) * k1_ic8h18 - (638 as f32 / 6 as f32) * k2_ic8h18 - (23 as f32 / 2 as f32) * k3_ic8h18 + (157 as f32 / 3 as f32) * k4_ic8h18 + (157 as f32 / 45 as f32) * k6_ic8h18));
           
            let k8_c8h18 = -self.k_1 * (c_c8h18_mod[i] + h * (- (53 as f32 / 14 as f32) * k1_c8h18 + (38 as f32 / 7 as f32) * k2_c8h18 - (3 as f32 / 14 as f32) * k3_c8h18 - (65 as f32 / 72 as f32) * k5_c8h18 + (29 as f32 / 90 as f32) * k7_c8h18));
            let k8_ic8h18 = self.k_1 * (c_c8h18_mod[i] + h * (- (53 as f32 / 14 as f32) * k1_c8h18 + (38 as f32 / 7 as f32) * k2_c8h18 - (3 as f32 / 14 as f32) * k3_c8h18 - (65 as f32 / 72 as f32) * k5_c8h18 + (29 as f32 / 90 as f32) * k7_c8h18)) 
- self.k_2 * (c_ic8h18_mod[i] + h * (- (53 as f32 / 14 as f32) * k1_ic8h18 + (38 as f32 / 7 as f32) * k2_ic8h18 - (3 as f32 / 14 as f32) * k3_ic8h18 - (65 as f32 / 72 as f32) * k5_ic8h18 + (29 as f32 / 90 as f32) * k7_ic8h18));
            let k8_c4h10 = self.k_2 * (c_ic8h18_mod[i] + h * (- (53 as f32 / 14 as f32) * k1_ic8h18 + (38 as f32 / 7 as f32) * k2_ic8h18 - (3 as f32 / 14 as f32) * k3_ic8h18 - (65 as f32 / 72 as f32) * k5_ic8h18 + (29 as f32 / 90 as f32) * k7_ic8h18));
            let k8_c4h8 = self.k_2 * (c_ic8h18_mod[i] + h * (- (53 as f32 / 14 as f32) * k1_ic8h18 + (38 as f32 / 7 as f32) * k2_ic8h18 - (3 as f32 / 14 as f32) * k3_ic8h18 - (65 as f32 / 72 as f32) * k5_ic8h18 + (29 as f32 / 90 as f32) * k7_ic8h18));

            let k9_c8h18 = -self.k_1 * (c_c8h18_mod[i] + h * ((56 as f32 / 25 as f32) * k1_c8h18 + (283 as f32 / 14 as f32) * k2_c8h18 - (119 as f32 / 6 as f32) * k3_c8h18 - (26 as f32 / 7 as f32) * k4_c8h18 - (13 as f32 / 15 as f32) * k5_c8h18 + (149 as f32 / 32 as f32) * k6_c8h18 - (25 as f32 / 9 as f32) * k7_c8h18 + (27 as f32 / 25 as f32) * k8_c8h18));
            let k9_ic8h18 = self.k_1 * (c_c8h18_mod[i] + h * ((56 as f32 / 25 as f32) * k1_c8h18 + (283 as f32 / 14 as f32) * k2_c8h18 - (119 as f32 / 6 as f32) * k3_c8h18 - (26 as f32 / 7 as f32) * k4_c8h18 - (13 as f32 / 15 as f32) * k5_c8h18 + (149 as f32 / 32 as f32) * k6_c8h18 - (25 as f32 / 9 as f32) * k7_c8h18 + (27 as f32 / 25 as f32) * k8_c8h18)) - self.k_2 * (c_ic8h18_mod[i] + h * ((56 as f32 / 25 as f32) * k1_ic8h18 + (283 as f32 / 14 as f32) * k2_ic8h18 - (119 as f32 / 6 as f32) * k3_ic8h18 - (26 as f32 / 7 as f32) * k4_ic8h18 - (13 as f32 / 15 as f32) * k5_ic8h18 + (149 as f32 / 32 as f32) * k6_ic8h18 - (25 as f32 / 9 as f32) * k7_ic8h18 + (27 as f32 / 25 as f32) * k8_ic8h18));
            let k9_c4h10 = self.k_2 * (c_ic8h18_mod[i] + h * ((56 as f32 / 25 as f32) * k1_ic8h18 + (283 as f32 / 14 as f32) * k2_ic8h18 - (119 as f32 / 6 as f32) * k3_ic8h18 - (26 as f32 / 7 as f32) * k4_ic8h18 - (13 as f32 / 15 as f32) * k5_ic8h18 + (149 as f32 / 32 as f32) * k6_ic8h18 - (25 as f32 / 9 as f32) * k7_ic8h18 + (27 as f32 / 25 as f32) * k8_ic8h18));
            let k9_c4h8 = self.k_2 * (c_ic8h18_mod[i] + h * ((56 as f32 / 25 as f32) * k1_ic8h18 + (283 as f32 / 14 as f32) * k2_ic8h18 - (119 as f32 / 6 as f32) * k3_ic8h18 - (26 as f32 / 7 as f32) * k4_ic8h18 - (13 as f32 / 15 as f32) * k5_ic8h18 + (149 as f32 / 32 as f32) * k6_ic8h18 - (25 as f32 / 9 as f32) * k7_ic8h18 + (27 as f32 / 25 as f32) * k8_ic8h18));

            c_c8h18_mod.push(c_c8h18_mod[i] + h/840.0 * (41.0*k1_c8h18 + 216.0*k4_c8h18 + 27.0*k5_c8h18 + 272.0*k6_c8h18 + 27.0*k7_c8h18 + 216.0*k8_c8h18 + 41.0*k9_c8h18));
            c_ic8h18_mod.push(c_ic8h18_mod[i] + h/840.0 * (41.0*k1_ic8h18 + 216.0*k4_ic8h18 + 27.0*k5_ic8h18 + 272.0*k6_ic8h18 + 27.0*k7_ic8h18 + 216.0*k8_ic8h18 + 41.0*k9_ic8h18));
            c_c4h10_mod.push(c_c4h10_mod[i] + h/840.0 * (41.0*k1_c4h10 + 216.0*k4_c4h10 + 27.0*k5_c4h10 + 272.0*k6_c4h10 + 27.0*k7_c4h10 + 216.0*k8_c4h10 + 41.0*k9_c4h10));
            c_c4h8_mod.push(c_c4h8_mod[i] + h/840.0 * (41.0*k1_c4h8 + 216.0*k4_c4h8 + 27.0*k5_c4h8 + 272.0*k6_c4h8 + 27.0*k7_c4h8 + 216.0*k8_c4h8 + 41.0*k9_c4h8));
        }
        return (c_c8h18_mod, c_ic8h18_mod, c_c4h10_mod, c_c4h8_mod);
    }
}

fn main() {
    let init = InitialData {
        temperature_1: 458.0,
        time: vec![0, 30, 60, 120, 240, 300],
        pressure: vec![4.21e3, 4.51e3, 4.73e3, 5.04e3, 5.40e3, 5.52e3],
        time_1: 250,
        universal_gas_equation_system: 8.31,
        temperature_2: 620.0,
        k_1: 0.12,
        k_2: 0.80,
        c_c8h18: 0.036e3,
        c_ic8h18: 0.0,
        c_c4h10: 0.0,
        c_c4h8: 0.0,
        dt: 0.01,
        number_of_iterations: (0..1000).collect(),
    };

    println!("Task_1:");

    let pressure_a = init.pressure_calc();
    let concentration_a = init.concentration_calc(pressure_a.clone());
    let concentration_anamorphosis_0 = init.anamorphosis0(concentration_a.clone());
    let concentration_anamorphosis_1 = init.anamorphosis1(concentration_a.clone());
    let concentration_anamorphosis_2 = init.anamorphosis2(concentration_a.clone());
    let rate_constant0 = init.rate_constant0_calc(concentration_a.clone());
    let rate_constant1 = init.rate_constant1_calc(concentration_a.clone());
    let rate_constant2 = init.rate_constant2_calc(concentration_a.clone());

    println!("\tTask_1_1:");
    println!("\t\tt, sec\t\tp_A, Pa\t\tC_A, mol/m^3");
    for i in 0..pressure_a.len() {
        println!("\t\tt[{}] = {}\tp_A[{}] = {}\tC_A[{}] = {}", i, init.time[i], i, pressure_a[i], i, concentration_a[i]);
    }
    println!("\n\tTask_1_2:");
    println!("\t\tk0, mol/(m^3 * sec)\tk1, 1/sec\t\tk2, 1/(mol * m^3 * sec))");
    for i in 0..pressure_a.len() {
        println!("\t\tk0[{}] = {}\tk1[{}] = {}\tk2[{}] = {}", i, rate_constant0[i], i, rate_constant1[i], i, rate_constant2[i]);
    }
    println!("\n\tTask_1_3:");
    println!("\t\ta_0, mol/m^3\t\ta_1, mol/m^3\t\ta_2, mol/m^3");
    for i in 0..pressure_a.len() {
        println!("\t\ta_0[{}] = {}\t\ta_1[{}] = {}\ta_2[{}] = {}", i, concentration_anamorphosis_0[i], i, concentration_anamorphosis_1[i], i, concentration_anamorphosis_2[i]);
    }

    write_to_file(concentration_a.clone(), init.time.clone(), "D:/lab_PhysChem/Lab6/C(dt)/C(dt).txt");
    write_to_file(concentration_anamorphosis_0.clone(), init.time.clone(), "D:/lab_PhysChem/Lab6/C_a0(dt)/C_a0(dt).txt");
    write_to_file(concentration_anamorphosis_1.clone(), init.time.clone(), "D:/lab_PhysChem/Lab6/C_a1(dt)/C_a1(dt).txt");
    write_to_file(concentration_anamorphosis_2.clone(), init.time.clone(), "D:/lab_PhysChem/Lab6/C_a2(dt)/C_a2(dt).txt");
	
    println!("\n\tTask_1_4:");
    let (half_turn_time, concentration_time_moment, tranformation_degree) = init.time_moment(init.time[1] as f32, rate_constant2[1], concentration_a[0], concentration_a[1]);
    println!("\n\tt(1/2) = {} sec\n\tC(t[1]) = {} mol/m^3\n\ta(t[1]) = {}", half_turn_time, concentration_time_moment, tranformation_degree);

    let mut c_c8h18_mod: Vec<f32> = vec![];
    c_c8h18_mod.push(init.c_c8h18);
    let mut c_ic8h18_mod: Vec<f32> = vec![];
    c_ic8h18_mod.push(init.c_ic8h18);
    let mut c_c4h10_mod: Vec<f32> = vec![];
    c_c4h10_mod.push(init.c_c4h10);
    let mut c_c4h8_mod: Vec<f32> = vec![];
    c_c4h8_mod.push(init.c_c4h8);
    (c_c8h18_mod, c_ic8h18_mod, c_c4h10_mod, c_c4h8_mod) = init.equation_system(c_c8h18_mod, c_ic8h18_mod, c_c4h10_mod, c_c4h8_mod);
    write_to_file_conc(c_c8h18_mod, init.number_of_iterations.clone(), "D:/lab_PhysChem/Lab6/Concentrations/C_C8H18/c_c8h18(t).txt");
    write_to_file_conc(c_ic8h18_mod, init.number_of_iterations.clone(), "D:/lab_PhysChem/Lab6/Concentrations/C_i-C8H18/c_ic8h18(t).txt");
    write_to_file_conc(c_c4h10_mod, init.number_of_iterations.clone(), "D:/lab_PhysChem/Lab6/Concentrations/C_C4H10/c_c4h10(t).txt");
    write_to_file_conc(c_c4h8_mod, init.number_of_iterations.clone(), "D:/lab_PhysChem/Lab6/Concentrations/C_C4H8/c_c4h8(t).txt");
}

// вывод в файл анаморфоз
fn write_to_file(vector: Vec<f32>, time: Vec<i32>, path: &str) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    for i in 0..time.len() {
        write!(file, "{} {}\n", time[i], vector[i])?;
    }
    Ok(())
}

// вывод в файл системы уравнений
fn write_to_file_conc(vector: Vec<f32>, number_of_iterations: Vec<usize>, path: &str) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    for i in 0..number_of_iterations.len() {
        write!(file, "{} {}\n", number_of_iterations[i], vector[i])?;
    }
    Ok(())
}
