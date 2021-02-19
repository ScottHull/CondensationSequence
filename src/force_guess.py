def get_guess(molecule, mass_balance, total_n):
    if molecule == 'Fe1_s':
        return 1.e-3 * mass_balance[molecule.strip("1_s")]
    elif molecule == 'Ni1_s':
        return 0.5e-1 * mass_balance[molecule.strip("1_s")]
    elif molecule == 'No1_s':
        return 8.e-14
    elif molecule == 'Ca1Al4O7_s':
        return 3.e-14
    elif molecule == 'Si1O2_s':
        return 1.e-6 * total_n
    elif molecule == 'Ca2Mg1Si2O7_s':
        return 1.e-13
    elif molecule == 'Ca2Al2Si1O7_s':
        return 8.e-14
    else:
        return 5.e-14
