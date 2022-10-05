data {
	// Number of participants
	int<lower=1> N_part;           
	
	// number of tasks completed in the second half
	int<lower=1> N_tasks_H2[N_part];  
	// the maximum number of tasks completed over all participants
	int<lower=1> max_N_tasks_H2;      
	
	// 1 ("Yes, found 3 diffs")
	// 0 ("No, I didn't find the 3 diffs")
	int<lower=0, upper=1> response_option_logit[N_part, max_N_tasks_H2]; 
	
	// participant-specific explanatory variables
	// 1 = stress, 0 = control
	vector[N_part] X_in_stress_group;
	
	// standardized cheating measure in the first half
	vector[N_part] X_cheat_H1;
}

parameters {
	real gamma_00;
	real gamma_std_cheat;
	real gamma_stress;
	real gamma_std_cheat_stress;
}

model {
	gamma_00 ~ normal(0, 10);
	gamma_std_cheat ~ normal(0, 10);
	gamma_stress ~ normal(0, 10);
	gamma_std_cheat_stress ~ normal(0, 10);
	
	for(i in 1:N_part) {
		for(j in 1:N_tasks_H2[i]) {
			response_option_logit[i, j] ~ bernoulli_logit(gamma_00 +
															gamma_std_cheat * X_cheat_H1[i] + 
															gamma_stress * X_in_stress_group[i] +
															gamma_std_cheat_stress * X_cheat_H1[i] * X_in_stress_group[i]);  
		}
	}
}
