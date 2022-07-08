import numpy as np
import matplotlib.pyplot as plt
import occurrence_rates_poisson as occ
import parameters_for_occurrence_rates as parms


if __name__ == '__main__':
    with open("occurrence_rates_logfile.dat", "a") as logfile:
        for k in range(len(parms.period_range)-1):
            for l in range(len(parms.mass_range)-1):
                plt.figure(l)
                n_pl_detected = occ.number_of_detected_planets(parms.period_range[k],parms.period_range[k+1], parms.mass_range[l], parms.mass_range[l+1], parms.planet_file_single)
                probability_of_frequency = []
                cummulative_probability_of_frequency = []
                cummulative_probability = 0
                for frequency in parms.planet_frequency:
                    number_of_planets = []
                    count = 0
                    while True:
                        if count == parms.number_of_simulations:
                            break
                        number_of_test_planets = np.random.binomial(parms.number_of_stars,frequency) #int(multiplicity_factor*np.random.binomial(number_of_stars, frequency))
                        mass, period = occ.create_planet_sample(number_of_test_planets, parms.period_range[k],parms.period_range[k+1], parms.mass_range[l], parms.mass_range[l+1])
                        number_of_planets.append(occ.calculate_expected_detections(mass, period))
                        count += 1
    
                    unique, counts = np.unique(number_of_planets, return_counts=True)
                    simulation_outcome = dict(zip(unique, counts))
    
                    try:
                        probability_of_frequency.append(simulation_outcome[n_pl_detected])
                        cummulative_probability += simulation_outcome[n_pl_detected]
                        cummulative_probability_of_frequency.append(cummulative_probability)
                    except KeyError:
                        probability_of_frequency.append(0.0)
                        cummulative_probability_of_frequency.append(cummulative_probability)
    
                half_of_simulation =  (np.abs(np.array(cummulative_probability_of_frequency)/cummulative_probability-0.5)).argmin()
                one_percent_of_simulation =  (np.abs(np.array(cummulative_probability_of_frequency)/cummulative_probability-0.018)).argmin()
                ninetynine_percent_of_simulation =  (np.abs(np.array(cummulative_probability_of_frequency)/cummulative_probability-0.982)).argmin()
    
                sixteen_percent_of_simulation =  (np.abs(np.array(cummulative_probability_of_frequency)/cummulative_probability-0.16)).argmin()
                eightyfour_percent_of_simulation =  (np.abs(np.array(cummulative_probability_of_frequency)/cummulative_probability-0.84)).argmin()
       
   
                logfile.write("Period-mass bin with period " + str(parms.period_range[k]) + "d to " + str(parms.period_range[k+1]) + "d and\n")            
                logfile.write("mass " + str(parms.mass_range[l]) + "M_Earth to " + str(parms.mass_range[l+1]) + "M_Earth.\n")
                logfile.write("Number of planets in this bin is " + str(n_pl_detected)+". All planets around " + parms.stellar_mass_bin + "  stars.  \n")
                logfile.write("Most probable frequency is " + str(parms.planet_frequency[half_of_simulation]) + " stars with planets.\n")
                logfile.write("Lower limit is " + str(parms.planet_frequency[sixteen_percent_of_simulation]) + " stars with planets.\n")
                logfile.write("Lowest limit is " + str(parms.planet_frequency[one_percent_of_simulation]) + " stars with planets.\n")
                logfile.write("Upper limit is " + str(parms.planet_frequency[eightyfour_percent_of_simulation]) + " stars with planets.\n")
                logfile.write("Uppest limit is " + str(parms.planet_frequency[ninetynine_percent_of_simulation]) + " stars with planets.\n")
