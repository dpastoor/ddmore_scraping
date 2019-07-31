# Experimental data from Jansen et al.(2004) J Pharm Biomed Anal 34:585-593 

Time<-c(0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0)
CPL_APAP_mcgL <- c(9.2797, 9.4681, 9.2923, 9.0719, 8.7246, 8.3994, 8.1881, 7.892, 7.3309, 6.9294, 6.5891)
CPL_AG_mcgL <- c(8.3935, 9.3632, 9.8072, 9.982, 10.0493, 9.926, 9.7817, 9.552, 9.1014, 8.6473, 8.2504)
CPL_AS_mcgL <- c(8.317, 8.7968, 8.9727, 8.887, 8.7968, 8.539, 8.4138, 8.1712, 7.6798, 7.2098, 6.8372)
CPL_APAP_mcgL <- exp(CPL_APAP_mcgL); CPL_AG_mcgL <- exp(CPL_AG_mcgL); CPL_AS_mcgL <- exp(CPL_AS_mcgL)
APAP_cal1.1 <- data.frame(Time, CPL_APAP_mcgL, CPL_AG_mcgL, CPL_AS_mcgL)

write.csv(data, file = "Real_APAP_data.csv")