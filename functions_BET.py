import numpy as np

def long_term_probability(phi_prior, lambda_prior, n, x, sample_length, time_window):
    if time_window > 1:
        raise ValueError('Error: the time window must be equal or lower than 1')
    if x > n:
        raise ValueError('Error: x must be lower or equal than n')
    if lambda_prior < 1:
        raise ValueError('Error: lambda must be equal or higher than 1')
    if type(lambda_prior) != int:
        raise ValueError('Error: lambda must be an integer')
    A = np.array([[phi_prior - 1, phi_prior], [1, 1]]) #equations to find alpha and beta from
    C = np.array([0, lambda_prior + 1]) #the mean and equivalent number of data
    alphb = np.linalg.solve(A, C) #solving the equations to find alpha and beta
    alpha_param = alphb[0] + x #merging alpha and x to find the alpha of the posterior
    beta_param = alphb[1] + n - x #merging beta, n and x to find beta of the posterior
    phi = np.random.beta(alpha_param, beta_param, sample_length) * time_window #posterior long-term
    return phi

def short_term_probability(Monitoring, Parameters, lambda_short_term, a_param, sample_length):
    if lambda_short_term < 1:
        raise ValueError('Error: lambda must be equal or higher than 1')
    if type(lambda_short_term) != int:
        raise ValueError('Error: lambda must be an integer')
    if a_param >= 1:
        raise ValueError("Error: 'a' must be lower than 1")
    P = np.delete(Parameters, 0, axis = 1) #deletes the first element of each row of the parameters input file
    if len(Monitoring[0]) != len(P[0]):
        raise ValueError('Error: the number of parameters does not match the number of weights and thresholds')
    W = P[0] #weights
    T = P[1:3] #thresholds
    deg_z = np.zeros_like(Monitoring) #vectors of zeros of the same shape of the monitoring vector, it's basically a matrix
    #where each row will be updated with the degrees of anomalies at the time of each observation
    for i in range(len(Monitoring)):
        for j in range(len(Monitoring[0])):
            a1 = T[0, j]
            a2 = T[1, j]
            ai = Monitoring[i, j]
            if a1 < a2:    #when the parameter is anomalous for high values        
                if ai <= a1:
                    zi = 0
                elif a1 < ai < a2:
                    zi = 0.5 * (np.sin(np.pi * ((ai - a1) / (a2 - a1)) - np.pi * 0.5) + 1) #equation 1a Marzocchi et al. 2024
                elif ai >= a2:
                    zi = 1
            if a1 > a2:    #when the parameter is anomalous for low values        
                if ai <= a2:
                    zi = 1
                elif a2 < ai < a1:
                    zi = 0.5 * (np.sin(np.pi * ((ai - a1) / (a2 - a1)) + np.pi * 0.5) + 1) #equation 1b Marzocchi et al. 2024
                elif ai >= a1:
                    zi = 0    
            deg_z[i, j] = zi
    Z_score = np.sum(deg_z * W, axis = 1) #finds the anomaly score as the weighted sum of the degrees of anomaly
    #axis = 1 means I'm making the row sum, if I did not specify it would have make the sum over the columns
    Phi_M = 1 - a_param * np.exp(-np.array(Z_score, dtype = float)) #expected value of the short-term 
    #magmatic unrest (equations, 5 and 7 of Marzocchi et al., 2024)
    alphb = np.zeros((2, len(Phi_M))) #vector of zeros to be updated with alphas and betas of the short-term prior over time
    alpha_param_M = np.zeros(len(Monitoring)) #vector of zeros for alphas
    beta_param_M = np.zeros(len(Monitoring)) #vector of zero for betas
    phi = np.zeros((len(Monitoring), sample_length)) #vector of zeros for the short-term sampled distributions of magma
    for i in range(len(Monitoring)):
        Phi_Mi = Phi_M[i]
        A = np.array([[Phi_Mi - 1, Phi_Mi], [1, 1]]) #finding alpha and beta
        C = np.array([0, lambda_short_term + 1])
        alphb[:, i] = np.linalg.solve(A, C) #vector of alphas and betas
        alpha_param_M[i] = alphb[0, i] #alphas
        beta_param_M[i] = alphb[1, i] #betas
        phi[i] = np.random.beta(alpha_param_M[i], beta_param_M[i], sample_length) #sampled short-term distribution 
    return phi

def eruption_probability_rescale_st(phi_3, t0, time_window):
    if time_window > 1:
        raise ValueError('Error: the time window must be equal or lower than 1')
    if t0 + time_window > 1:
        raise ValueError('Error: the sum between t0 and time_window must be equal or lower than 1')
    tau = 1 #original time window (1 month)
    Tc = 1/4 #characteristic time (1 week)
    t = np.linspace(0, tau, 1000) #time vector for the distribution that shows how the probability of eruption
    #evolves inside the time-window (see Selva et al., 2014)
    k = 1 / (Tc * (1 - np.exp(-tau / Tc))) #normalizing factor
    ft = k * np.exp(-t / Tc) #equation 7 of Selva et al., 2014 (pdf)
    Ft = Tc * k - Tc * k * np.array(np.exp(-t / Tc)) #integral of equation 7 (cdf)
    Ft0 = np.interp(t0, t, Ft) #value of Ft at the time of observation
    Ftpdt = np.interp(t0 + time_window, t, Ft) #value of Ft at time of observation + dt
    #recaling the probability of eruption
    for i in range(len(phi_3)):
        phi3i = phi_3[i]
        phi_3[i] = phi3i * (Ftpdt - Ft0) #equation 6 of Selva et al. 2014
    return phi_3

def absolute_probability(phi_2, phi_3, percentiles, phi_1 = [1]):
    if len(phi_1) != 1:
        P_unrest = phi_1
        Mean_unrest = np.mean(P_unrest)
        PRC_unrest = np.percentile(P_unrest, percentiles)
        P_magma = P_unrest * phi_2
        Mean_magma = np.mean(P_magma)
        PRC_magma = np.percentile(P_magma, percentiles)
        P_eruption = P_magma * phi_3
        Mean_eruption = np.mean(P_eruption)
        PRC_eruption = np.percentile(P_eruption, percentiles)
        return Mean_unrest, PRC_unrest, Mean_magma, PRC_magma, Mean_eruption, PRC_eruption
    else: 
        if len(phi_2) != len(phi_3):
            raise ValueError('Error: the Magmatic Unrest and Eruption sheets in the monitoring file must have the same length')
        P_magma = phi_2
        Mean_magma = np.zeros(len(phi_2)) #vector of zeros to be updated with the mean
        for i in range(len(phi_2)):
            Mean_magma[i] = np.mean(P_magma[i]) #expected value of the distributions over time
        PRC_magma = np.zeros((len(phi_2), len(percentiles))) #zeros to be updated with percentiles
        for i in range(len(phi_2)):
            PRC_magma[i] = np.percentile(P_magma[i], percentiles) #variation of the percentiles    
        P_eruption = P_magma * phi_3
        Mean_eruption = np.zeros(len(phi_3)) #vector of zeros to be updated with the mean
        for i in range(len(phi_3)):
            Mean_eruption[i] = np.mean(P_eruption[i]) #expected value of the distributions over time
        PRC_eruption = np.zeros((len(phi_3), len(percentiles))) #zeros to be updated with percentiles
        for i in range(len(phi_3)):
            PRC_eruption[i] = np.percentile(P_eruption[i], percentiles) #variation of the percentiles
        return Mean_magma, PRC_magma, Mean_eruption, PRC_eruption    