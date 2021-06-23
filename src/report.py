import os
import shutil

def build_subdirs(base_path):
    dirs = ['percent_condensed', 'number_densities', 'mass_balance']
    if os.path.exists(base_path):
        shutil.rmtree(base_path)
    os.mkdir(base_path)
    for i in dirs:
        os.mkdir(base_path + "/{}".format(i))


def write_report(path, iteration, temperature, percent_condensed_dict, number_densities, mass_balance):
    condensed_outfile = open(path + "/percent_condensed/{}.csv".format(iteration), 'w')
    number_outfile = open(path + "/number_densities/{}.csv".format(iteration), 'w')
    mass_balance_outfile = open(path + "/mass_balance/{}.csv".format(iteration), 'w')
    condensed_outfile.write("{}\n".format(temperature))
    number_outfile.write("{}\n".format(temperature))
    mass_balance_outfile.write("{}\n".format(temperature))

    for i in percent_condensed_dict.keys():
        condensed_outfile.write("{},{}\n".format(i, percent_condensed_dict[i][-1]['percent']))
    for i in number_densities.keys():
        number_outfile.write("{},{}\n".format(i, number_densities[i]))
    for i in mass_balance.keys():
        mass_balance_outfile.write("{},{}\n".format(i, mass_balance[i]))
    condensed_outfile.close()
    number_outfile.close()
    mass_balance_outfile.close()
