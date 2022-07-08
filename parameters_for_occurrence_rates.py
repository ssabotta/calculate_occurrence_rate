import numpy as np

## This is specific for CARMENES detection maps. Do not change.
LOG_MINIMUM_PERIOD = 0 # in days
LOG_MAXIMUM_PERIOD = 4 # in days
LOG_MINIMUM_MASS = -2.502241 # in Jupiter mass
LOG_MAXIMUM_MASS = 0.497759 # in Jupiter mass
GRIDPOINTS_PERIOD = 60
GRIDPOINTS_MASS = 60

period_in_days = np.logspace(LOG_MINIMUM_PERIOD, LOG_MAXIMUM_PERIOD,
                      num=GRIDPOINTS_PERIOD, endpoint=True, base=10.0)
minimum_mass_in_jupiter_masses = np.logspace(LOG_MINIMUM_MASS, LOG_MAXIMUM_MASS,
                      num=GRIDPOINTS_MASS,  endpoint=True, base=10.0)
minimum_mass_in_earth_masses = minimum_mass_in_jupiter_masses*317.828


filename = {"all":"all_stars.dat","low_mass":"all_m_star_lower_0p337msol.dat","high_mass":"all_m_star_higher_0p337msol.dat"}
number_of_stars_dict = {"all":238,"low_mass":90,"high_mass":148}
planet_file_list = {"all":"CARMENES_planets.csv","low_mass":"CARMENES_planets_mstar_lower_0p337.csv", \
           "high_mass":"CARMENES_planets_mstar_higher_0p337.csv"}
   
planet_file_list_single = {"all":"CARMENES_planets_single_planets.csv","low_mass":"CARMENES_planets_single_planets_mstar_lower_0p337.csv", \
                        "high_mass":"CARMENES_planets_single_planets_mstar_higher_0p337.csv"}
    
### Start making changes here

stellar_mass_bin = "all" # or "low_mass" or "high_mass"
    
planet_file = planet_file_list[stellar_mass_bin]
planet_file_single = planet_file_list_single[stellar_mass_bin]

number_of_stars = number_of_stars_dict[stellar_mass_bin]


name_of_writefile = stellar_mass_bin + ".dat"
number_of_simulations = 1000
frequency_resolution = 0.001
frequency_start, frequency_end = 0.001, 0.1# 1.0#2.5

period_range = [1,1000] 
mass_range = [100,1000]

number_of_points = int((frequency_end-frequency_start)/frequency_resolution)
planet_frequency = np.linspace(frequency_start,frequency_end, num=number_of_points, endpoint=True, dtype=None)
