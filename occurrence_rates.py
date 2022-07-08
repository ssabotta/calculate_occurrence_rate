import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import parameters_for_occurrence_rates as parms

def mass_distribution(mass):
    """Define mass distribution of test planets here.

    Default is two power laws for m sin i > 50 M_Earth and m sin i < 50 M_Earth (see Sabotta+2021)

    For log-uniform change to:     
    weight = mass/mass
    return weight
    """
    
    idx = (np.abs(mass-50)).argmin()
    
    weight_part1 = 319.2509222787752*mass[0:idx+1]**(-1.062040969485353)
    weight_part2 = 7.237015175751411*mass[idx+1:len(mass)]**(-0.170547275)
    #weight = mass/mass
    return np.append(weight_part1,weight_part2) # weight


def period_distribution(period):
    """Define period distribution of test planets here.

    Default is log-uniform.
    """
    weight = period/period
    
    return weight


def create_planet_sample(number_of_planets,period_lower_limit,period_upper_limit, mass_lower_limit, mass_upper_limit):
    """Create a single set of test planets
    
    returns two numpy arrays with mass and period of the planets
    """

    index_mass_low = (np.abs(parms.minimum_mass_in_earth_masses-mass_lower_limit)).argmin()
    index_mass_high = (np.abs(parms.minimum_mass_in_earth_masses-mass_upper_limit)).argmin()
    index_period_low = (np.abs(parms.period_in_days-period_lower_limit)).argmin()
    index_period_high = (np.abs(parms.period_in_days-period_upper_limit)).argmin()
    
    probabilities_mass = mass_distribution(parms.minimum_mass_in_earth_masses[index_mass_low:index_mass_high])
    probabilities_period = period_distribution(parms.period_in_days[index_period_low:index_period_high])
    
    probabilities_mass = probabilities_mass/sum(probabilities_mass)
    probabilities_period = probabilities_period/sum(probabilities_period)
        

    mass_dummy = np.random.choice(parms.minimum_mass_in_earth_masses[index_mass_low:index_mass_high], number_of_planets, p= probabilities_mass)
    period_dummy = np.random.choice(parms.period_in_days[index_period_low:index_period_high], number_of_planets, p= probabilities_period)

    return mass_dummy, period_dummy


def number_of_detected_planets(period_lower_limit,period_upper_limit, mass_lower_limit, mass_upper_limit, planet_list):
    """Find out how many CARMENES planets were detected in your period mass bin 
    
    input csv file with planet list
    returns number of planets as integer
    """
    planets = pd.read_csv(planet_list, delimiter = ",")

    number_of_planets = np.sum((period_lower_limit < planets.P_orb) & (planets.P_orb < period_upper_limit) & (mass_lower_limit < planets.M_sini) & (planets.M_sini < mass_upper_limit))
    return number_of_planets


def calculate_expected_detections(masses,period):
    """Find out how many of the test planet would be detected from the global detection map.
    """
    planet_list = []
    detection_map = pd.read_csv(parms.filename[parms.stellar_mass_bin], delimiter="\t") 

    for i in range(len(masses)):
        test = detection_map.loc[(np.isclose(detection_map["mass"],masses[i], atol=10**-6)) & (np.isclose(detection_map["period"], period[i], atol=10**-6))]
        
        probability = 0.01*float(np.array(test.probability))
        planet_list.append(np.random.choice([0,1],p=[1-probability,probability]))
        
    return np.sum(planet_list)


if __name__ == '__main__':
    with open("occurrence_rates_logfile.dat", "a") as logfile:
        for k in range(len(parms.period_range)-1):
            for l in range(len(parms.mass_range)-1):
                plt.figure(l)
                n_pl_detected = number_of_detected_planets(parms.period_range[k],parms.period_range[k+1], parms.mass_range[l], parms.mass_range[l+1], parms.planet_file)
                probability_of_frequency = []
                cummulative_probability_of_frequency = []
                cummulative_probability = 0
                for frequency in parms.planet_frequency:
                    number_of_planets = []
                    count = 0
                    while True:
                        if count == parms.number_of_simulations:
                            break
                        number_of_test_planets = np.random.poisson(parms.number_of_stars*frequency) #int(multiplicity_factor*np.random.binomial(number_of_stars, frequency))
                        mass, period = create_planet_sample(number_of_test_planets, parms.period_range[k],parms.period_range[k+1], parms.mass_range[l], parms.mass_range[l+1])
                        number_of_planets.append(calculate_expected_detections(mass, period))
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
                logfile.write("Number of planets in this bin is " + str(n_pl_detected)+". All planets around " + parms.stellar_mass_bin + "  stars. \n")
                logfile.write("Most probable frequency is " + str(parms.planet_frequency[half_of_simulation]) + " planets per star.\n")
                logfile.write("Lower limit is " + str(parms.planet_frequency[sixteen_percent_of_simulation]) + " planets per star.\n")
                logfile.write("Lowest limit is " + str(parms.planet_frequency[one_percent_of_simulation]) + " planets per star.\n")
                logfile.write("Upper limit is " + str(parms.planet_frequency[eightyfour_percent_of_simulation]) + " planets per star.\n")
                logfile.write("Uppest limit is " + str(parms.planet_frequency[ninetynine_percent_of_simulation]) + " planets per star.\n")
