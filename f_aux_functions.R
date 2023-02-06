f_calc_std_cheating_H1 <- function(cheatable_tasks) {
	use_cheating_at_T1 <- cheatable_tasks %>%
		filter(study_half_number %in% c("H1")) %>%
		group_by(participant_ID) %>%
		summarise(n_cheating = sum(is_cheating),
			    n_tasks = n()) %>%
		ungroup()
	use_cheating_at_T1 <- use_cheating_at_T1 %>%
		mutate(perc_cheating = n_cheating / n_tasks)
	
	mean_cheating <- mean(use_cheating_at_T1$perc_cheating)
	sd_cheating   <- sd(use_cheating_at_T1$perc_cheating)
	
	use_cheating_at_T1 <- use_cheating_at_T1 %>%
		mutate(std_perc_cheating_H1 = (perc_cheating - mean_cheating) / sd_cheating) %>%
		select(participant_ID,
			 std_perc_cheating_H1)
	
	return(use_cheating_at_T1)
}

f_calc_std_cheat_percentile_group <- function(use_cheating_at_T1) {
	bottom_5perc  <- quantile(use_cheating_at_T1$std_perc_cheating_H1, probs = c(0.05))
	bottom_10perc <- quantile(use_cheating_at_T1$std_perc_cheating_H1, probs = c(0.10))
	bottom_15perc <- quantile(use_cheating_at_T1$std_perc_cheating_H1, probs = c(0.15))
	top_5perc     <- quantile(use_cheating_at_T1$std_perc_cheating_H1, probs = c(0.95))
	top_10perc    <- quantile(use_cheating_at_T1$std_perc_cheating_H1, probs = c(0.90))
	top_15perc    <- quantile(use_cheating_at_T1$std_perc_cheating_H1, probs = c(0.85))
	
	use_cheating_at_T1 <- use_cheating_at_T1 %>%
		mutate(is_between_bottom_top_5perc = (std_perc_cheating_H1 >= bottom_5perc)*
			 	(std_perc_cheating_H1 <= top_5perc),
			 is_between_bottom_top_10perc = (std_perc_cheating_H1 >= bottom_10perc)*
			 	(std_perc_cheating_H1 <= top_10perc),
			 is_between_bottom_top_15perc = (std_perc_cheating_H1 >= bottom_15perc)*
			 	(std_perc_cheating_H1 <= top_15perc))
	return(use_cheating_at_T1)
}

f_filter_percentile_cheat <- function(all_tasks, all_part, exclude_percentile) {
	df_part <- f_calc_std_cheating_H1(cheatable_tasks = all_tasks)
	
	if(exclude_percentile == "no_exclusion") {
		df_part <- df_part %>% 
			select(participant_ID,
				 std_perc_cheating_H1)
	}
	
	if(exclude_percentile == "top_bottom_5") {
		df_part <- f_calc_std_cheat_percentile_group(use_cheating_at_T1 = df_part)
		df_part <- df_part %>% 
			filter(is_between_bottom_top_5perc == 1) %>%
			select(participant_ID,
				 std_perc_cheating_H1)
	}
	
	if(exclude_percentile == "top_bottom_10") {
		df_part <- f_calc_std_cheat_percentile_group(use_cheating_at_T1 = df_part)
		df_part <- df_part %>% 
			filter(is_between_bottom_top_10perc == 1) %>%
			select(participant_ID,
				 std_perc_cheating_H1)
	}
	
	if(exclude_percentile == "top_bottom_15") {
		df_part <- f_calc_std_cheat_percentile_group(use_cheating_at_T1 = df_part)
		df_part <- df_part %>% 
			filter(is_between_bottom_top_15perc == 1) %>%
			select(participant_ID,
				 std_perc_cheating_H1)
	}
	
	return(df_part)
}

f_sum_results <- function(par_names, stan_fit) {
	df_post_draw <- rstan::extract(stan_fit, pars = par_names)
	df_post_draw <- data.frame(df_post_draw)
	df_post_draw <- df_post_draw %>%
		pivot_longer(cols = everything(),
				 names_to = "param_name",
				 values_to = "post_draw")
	df_post_draw <- df_post_draw %>%
		group_by(param_name) %>%
		mutate(post_mean = mean(post_draw),
			 post_sd = sd(post_draw),
			 p_value = f_calc_pvalue(post_draw),
			 lower_CI = quantile(post_draw, probs = c(0.025)),
			 upper_CI = quantile(post_draw, probs = c(0.975))) %>% 
		ungroup()
	df_post_draw <- distinct(df_post_draw %>%
					 	select(-post_draw))
	
	df_post_draw <- df_post_draw %>% 
		select(param_name,
			 post_mean, post_sd, p_value, lower_CI, upper_CI)
	return(df_post_draw)
}

f_calc_pvalue <- function(x) {
	return ((sum(x > 0) * (mean(x) < 0) + sum(x < 0) * (mean(x > 0))) / length(x))
}

f_prep_df_plot_interaction <- function(choice_fit, 
						   grid_std_cheating,
						   mean_cheating, 
						   sd_cheating) {
	
	choice_draws <- rstan::extract(choice_fit, 
						 pars = c("gamma_00", 
						 	   "gamma_std_cheat", 
						 	   "gamma_stress", 
						 	   "gamma_std_cheat_stress"))
	choice_draws <- data.frame(choice_draws)
	
	aux_list <- NULL
	index_list <- 1
	
	for(std_cheat_value in grid_std_cheating) {
		aux_df <- choice_draws %>%
			mutate(u_cheat_control = gamma_00 + gamma_std_cheat * std_cheat_value,
				 u_cheat_stress  = gamma_00 + gamma_std_cheat * std_cheat_value +
				 	gamma_stress + gamma_std_cheat_stress*std_cheat_value)
		aux_df <- aux_df %>%
			mutate(prob_control = exp(u_cheat_control) / (exp(u_cheat_control) + 1),
				 prob_stress  = exp(u_cheat_stress) / (exp(u_cheat_stress) + 1))
		
		aux_list[[index_list]] <- aux_df %>%
			summarise(mean_control = mean(prob_control),
				    LCI_control  = quantile(prob_control, probs = c(0.025)),
				    UCI_control  = quantile(prob_control, probs = c(0.975)),
				    mean_stress = mean(prob_stress),
				    LCI_stress  = quantile(prob_stress, probs = c(0.025)),
				    UCI_stress  = quantile(prob_stress, probs = c(0.975))) %>%
			mutate(std_cheat_value = std_cheat_value,
				 perc_cheat_value = std_cheat_value * sd_cheating + mean_cheating)
		
		index_list <- index_list + 1
	}
	df_plot <- map_dfr(aux_list, rbind)
	
	df_plot <- df_plot %>%
		pivot_longer(cols = c(contains("control"),
					    contains("stress")),
				 values_to = "cheat_prob",
				 names_to = "group_level")
	df_plot <- df_plot %>%
		mutate(experimental_group = if_else(str_detect(group_level, "control"),
								"control", "stress"))
	df_plot <- df_plot %>%
		mutate(aggregation_level = case_when(str_detect(group_level, "mean") ~ "post_mean",
								 str_detect(group_level, "LCI") ~ "post_LCI",
								 str_detect(group_level, "UCI") ~ "post_UCI",
								 TRUE ~ "check_this"))
	
	df_plot <- df_plot %>%
		select(-group_level)
	df_plot <- df_plot %>%
		pivot_wider(values_from = cheat_prob,
				names_from = aggregation_level)
	
	return(df_plot)
}


f_plot_interaction_3in1 <- function(choice_fit, grid_std_cheating,
						mean_cheating, sd_cheating) {
	
	df_plot <- f_prep_df_plot_interaction(choice_fit = choice_fit, 
							  grid_std_cheating = grid_std_cheating,
							  mean_cheating = mean_cheating, 
							  sd_cheating = sd_cheating)
	
	plot1 <- df_plot %>%
		ggplot(aes(x = perc_cheat_value, group = experimental_group)) +
		geom_ribbon(aes(ymin = post_LCI, 
				    ymax = post_UCI, 
				    fill = experimental_group)) +
		geom_line(aes(y = post_mean, color = experimental_group)) + 
		theme_bw() +
		ylab("Probability of cheating at t=2") +
		xlab("Percent cheating at t=1") + 
		scale_fill_manual("Experimental group",
					values = c("control" = "grey75",
						     "stress" = "lightblue")) +
		scale_color_manual("Experimental group",
					 values = c("control" = "grey50",
					 	     "stress" = "blue")) +
		coord_equal(ratio = 1) +
		scale_y_continuous(limits = c(0, 1)) +
		scale_x_continuous(limits = c(0, 1)) +
		theme(legend.position = "bottom",
			plot.title = element_text(size = 12),
			plot.subtitle = element_text(size = 8)) +
		labs(title = "a)")
	plot1
	
	plot2 <- df_plot %>%
		filter(experimental_group == "control") %>%
		ggplot(aes(x = perc_cheat_value, group = experimental_group)) +
		geom_ribbon(aes(ymin = post_LCI, 
				    ymax = post_UCI, 
				    fill = experimental_group)) +
		geom_line(aes(y = post_mean, color = experimental_group)) + 
		theme_bw() +
		ylab("Probability of cheating at t=2") +
		xlab("Percent cheating at t=1") + 
		scale_fill_manual("Experimental group",
					values = c("control" = "grey75",
						     "stress" = "lightblue")) +
		scale_color_manual("Experimental group",
					 values = c("control" = "grey50",
					 	     "stress" = "blue")) +
		coord_equal(ratio = 1) +
		scale_y_continuous(limits = c(0, 1)) +
		scale_x_continuous(limits = c(0, 1)) +
		geom_abline(intercept = 0,
				slope = 1,
				linewidth = 0.5,
				linetype = "dashed") +
		theme(plot.title = element_text(size = 12),
			plot.subtitle = element_text(size = 8),
			legend.position = "none") +
		labs(title = "b)")
	
	plot3 <- df_plot %>%
		filter(experimental_group == "stress") %>%
		ggplot(aes(x = perc_cheat_value, group = experimental_group)) +
		geom_ribbon(aes(ymin = post_LCI, 
				    ymax = post_UCI, 
				    fill = experimental_group)) +
		geom_line(aes(y = post_mean, color = experimental_group)) + 
		theme_bw() +
		ylab("Probability of cheating at t=2") +
		xlab("Percent cheating at t=1") + 
		scale_fill_manual("Experimental group",
					values = c("control" = "grey75",
						     "stress" = "lightblue")) +
		scale_color_manual("Experimental group",
					 values = c("control" = "grey50",
					 	     "stress" = "blue")) +
		coord_equal(ratio = 1) +
		scale_y_continuous(limits = c(0, 1)) +
		scale_x_continuous(limits = c(0, 1)) +
		geom_abline(intercept = 0,
				slope = 1,
				linewidth = 0.5,
				linetype = "dashed") +
		theme(plot.title = element_text(size = 12),
			plot.subtitle = element_text(size = 8),
			legend.position = "none") +
		labs(title = "c)")
	
	lay <- rbind(c(1,1,1,2,2),
			 c(1,1,1,3,3))
	
	gridExtra::grid.arrange(plot1, plot2, plot3,
					layout_matrix = lay)
}

f_plot_interaction_single <- function(choice_fit, 
						  grid_std_cheating,
						  mean_cheating, 
						  sd_cheating) {
	
	df_plot <- f_prep_df_plot_interaction(choice_fit = choice_fit, 
							  grid_std_cheating = grid_std_cheating,
							  mean_cheating = mean_cheating, 
							  sd_cheating = sd_cheating)
	
	df_plot %>%
		ggplot(aes(x = perc_cheat_value, group = experimental_group)) +
		geom_ribbon(aes(ymin = post_LCI, 
				    ymax = post_UCI, 
				    fill = experimental_group)) +
		geom_line(aes(y = post_mean, color = experimental_group)) + 
		geom_abline(intercept = 0,
				slope = 1,
				linewidth = 0.5,
				linetype = "dashed") +
		theme_bw() +
		ylab("Probability of cheating at t=2") +
		xlab("Percent cheating at t=1") + 
		scale_fill_manual("Experimental group",
					values = c("control" = "grey75",
						     "stress" = "lightblue")) +
		scale_color_manual("Experimental group",
					 values = c("control" = "grey50",
					 	     "stress" = "blue")) +
		coord_equal(ratio = 1) +
		scale_y_continuous(limits = c(0, 1)) +
		scale_x_continuous(limits = c(-0.01, 1.01)) +
		theme(legend.position = "bottom",
			plot.title = element_text(size = 12),
			plot.subtitle = element_text(size = 8))
}


f_prep_use_tasks <- function(part_info, 
				     part_and_task_info,
				     exclusion_C1, 
				     exclusion_C2) {
	
	if(!(exclusion_C1 %in% c("no_exclusion", 
					 "Hand_out"))) {
		stop("the exclusion criteria 1 is not specified")
	}
	
	if(!(exclusion_C2 %in% c("no_exclusion", 
					 "top_bottom_5", 
					 "top_bottom_10", 
					 "top_bottom_15"))) {
		stop("the exclusion criteria 2 (based on percentile cheating) is not specified")
	}
	
	if(exclusion_C1 == "no_exclusion") {
		part_info <- part_info %>%
			select(participant_ID,
				 experimental_cond)
	}
	
	if(exclusion_C1 == "Hand_out") {
		cat("you exclude participants who take the hand out")
		part_info <- part_info %>%
			filter(hand_taken_out == 0) %>%
			select(participant_ID,
				 experimental_cond)
	}
	
	part_and_task_info <- part_info %>%
		left_join(part_and_task_info,
			    by = "participant_ID",
			    multiple = "all")
	
	# keep tasks of interest
	part_and_task_info <- part_and_task_info %>%
		filter(type_of_task %in% c("noclick")) %>%
		filter(n_differences %in% c(1, 2)) %>%
		filter(part_response %in% c("No, did not find 3 diffs", "Yes, found 3 diffs")) %>%
		filter(reward %in% c(0.2, 0.4))
	
	# add is_cheating label
	part_and_task_info <- part_and_task_info %>% 
		mutate(is_cheating = if_else(part_response == "Yes, found 3 diffs", 1, 0)) %>%
		mutate(study_half_number = if_else(is_in_second_half == 0, "H1", "H2"))
	
	# calculate cheating behavior at T1
	# this will be used both as X in choice models and to exclude participants
	# now, if I want to filter based on percentile cheating at H1, I do it
	# if I don't exclude anyone, then I'm only adding the std_cheating info
	part_info <- f_filter_percentile_cheat(all_tasks = part_and_task_info,
							   all_part = part_info,
							   exclude_percentile = exclusion_C2)
	
	# keep only observations from the participants I have now
	part_and_task_info <- part_info %>%
		left_join(part_and_task_info,
			    by = "participant_ID",
			    multiple = "all")
	
	return(list(use_tasks = part_and_task_info,
			use_participants = part_info))
}


f_prep_stan_data <- function(use_tasks) {
	
	part_and_task_info <- use_tasks %>%
		filter(is_in_second_half == 1) %>%
		select(participant_ID, 
			 experimental_cond,
			 std_perc_cheating_H1,
			 is_cheating)
	part_info <- part_and_task_info %>%
		group_by(participant_ID) %>%
		mutate(n_tasks_H2 = n()) %>%
		ungroup()
	part_info <- distinct(part_info %>%
				    	select(participant_ID, 
				    		 experimental_cond,
				    		 std_perc_cheating_H1,
				    		 n_tasks_H2))
	part_info <- part_info %>%
		mutate(part_number = row_number())
	part_and_task_info <- part_and_task_info %>%
		left_join(part_info %>%
			    	select(participant_ID,
			    		 part_number),
			    by = "participant_ID",
			    multiple = "all")
	
	N_part <- nrow(part_info)
	N_tasks_H2 <- array(0, dim = c(N_part))
	for(i in 1:N_part) {
		aux_df <- part_info %>%
			filter(part_number == i)
		N_tasks_H2[i] <- aux_df[["n_tasks_H2"]]
	}
	
	response_option_logit <- array(0, dim = c(N_part, max(N_tasks_H2)))
	X_in_stress_group <- array(0, dim = c(N_part))
	X_cheat_H1 <- array(0, dim = c(N_part))
	
	for(i in 1:N_part) {
		aux_ij <- part_and_task_info %>%
			filter(part_number == i)
		
		aux_0j <- part_info %>%
			filter(part_number == i)
		
		response_option_logit[i, 1:nrow(aux_ij)] <- aux_ij[["is_cheating"]]
		X_in_stress_group[i] <- (aux_0j[["experimental_cond"]] == "Stress")*1
		X_cheat_H1[i] <- aux_0j[["std_perc_cheating_H1"]]
	}
	
	stan_data <- list(N_part = N_part,
				N_tasks_H2 = N_tasks_H2,
				max_N_tasks_H2 = max(N_tasks_H2),
				response_option_logit = response_option_logit,
				X_in_stress_group = X_in_stress_group,
				X_cheat_H1 = X_cheat_H1)
	
	return(stan_data)
}

