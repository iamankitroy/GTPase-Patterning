#!/Users/roy/anaconda3/bin/python

import numpy as np
import pickle

def initiate_actin(time, pol_time, scale, rates):
	# 110% actin length
	max_length = int((pol_time * rates["k_pol"]) * 1.1)
	# actin filament
	scaled_time = int(time/scale)
	actin = np.zeros(shape = (scaled_time, max_length))

	return actin


def events_per_timescale(rates, scale):
	return {key: value*scale for key, value in zip(rates.keys(), rates.values())}


def polymerize(filament, end_pos):
	if np.random.uniform() <= rates_scaled["k_pol"]:
		filament[end_pos] = 3
		end_pos += 1

	return filament, end_pos


def hydrolyze_atp(filament):
	atp_subunits = (filament == 3).astype(int)

	change_probs = (np.random.uniform(size = filament.size) <= rates_scaled["k_adp_pi"]).astype(int)

	change_state = (atp_subunits * change_probs).astype(bool)

	filament[change_state] = 2

	return filament


def release_internalPi(filament, end_pos):
	adp_pi_subunits = (filament == 2).astype(int)

	change_probs = (np.random.uniform(size = filament.size) <= rates_scaled["k_adp_int"]).astype(int)

	change_state = (adp_pi_subunits * change_probs).astype(bool)
	change_state[end_pos] = False

	filament[change_state] = 1

	return filament


def release_terminalPi(filament, end_pos):
	if (filament[end_pos] == 2) and (np.random.uniform() <= rates_scaled["k_adp_ter"]):
		filament[end_pos] = 1

	return filament


def depolymerize(filament, end_pos):

	prob = np.random.uniform()

	if (end_pos == 3) and (prob <= rates_scaled["k_atp_depol"]):
		filament[end_pos] = 0
		end_pos -= 1

	elif (end_pos == 2) and (prob <= rates_scaled["k_adp_pi_depol"]):
		filament[end_pos] = 0
		end_pos -= 1

	elif prob <= rates_scaled["k_adp_depol"]:
		filament[end_pos] = 0
		end_pos -= 1

	return filament, end_pos


def cycle_polymerization(actin, end_pos, time, pol_time):

	scaled_time = int(time/scale)
	scaled_pol_time = int(pol_time/scale)

	simulation_results = []

	for t in range(1, scaled_time):
		filament = actin[t-1]
		actin[t] = hydrolyze_atp(filament)
		actin[t] = release_internalPi(filament, end_pos)
		actin[t] = release_terminalPi(filament, end_pos)

		if t <= scaled_pol_time:
			actin[t], end_pos = polymerize(filament, end_pos)
		else:
			actin[t], end_pos = depolymerize(filament, end_pos)

		simulation_results.append((t,
		end_pos,
		round(sum(actin[t]==3)/end_pos, 3),
		round(sum(actin[t]==2)/end_pos, 3),
		round(sum(actin[t]==1)/end_pos, 3)))

		print(t,
			end_pos,
			round(sum(actin[t]==3)/end_pos, 3),
			round(sum(actin[t]==2)/end_pos, 3),
			round(sum(actin[t]==1)/end_pos, 3))

	writeOUT(simulation_results, actin)


def writeOUT(simulation_results, actin):

	with open("actin_time-evol.p", "wb") as fh:
		pickle.dump(actin, fh)

	with open("simulation_results.p", "wb") as fh:
		pickle.dump(simulation_results, fh)


def main():

	global actin
	global rates_scaled
	global scale

	rates = {
		"k_pol" : 320,				# monomers per second
		"k_adp_pi" : 0.3,			# per second
		"k_adp_int" : 0.003,		# per second
		"k_adp_ter" : 2,			# per second
		"k_atp_depol" : 1.4,		# per second
		"k_adp_pi_depol" : 0.2,		# per second
		"k_adp_depol" : 5.4			# per second
	}

	# simulation time
	time = 5 * 60					# seconds
	pol_time = 10					# seconds
	scale = 1e-3					# scale to milliseconds

	# initiate end position
	end_pos = 0

	actin = initiate_actin(time, pol_time, scale, rates)
	rates_scaled = events_per_timescale(rates, scale)
	cycle_polymerization(actin, end_pos, time, pol_time)

main()

# Ankit Roy
# 21st February, 2021
