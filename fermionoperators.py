import numpy as np

"""
UCCSD_operator takes a single reference state as a binary string which the
cluster operator will act on. The occupation of the lowest orbitals is represented
by the least significant bits.
"""
def UCCSD_operator(ref_str):

    reversed_ref_str = ref_str[::-1]
    virtual_spin_orbitals = []
    occupied_spin_orbitals = []

    for idx in range(len(ref_str)):
        if reversed_ref_str[idx] == '1':
            occupied_spin_orbitals.append(idx)
        elif reversed_ref_str[idx] == '0':
            virtual_spin_orbitals.append(idx)
        else:
            raise ValueError("String", ref_str, "is not a valid reference state.")

    print("occupied:", occupied_spin_orbitals)
    print("virtual:", virtual_spin_orbitals)

    single_excitations = []

    for occ in occupied_spin_orbitals:
        operators = []
        for virt in virtual_spin_orbitals:
            operators.append( ("a_{}^dagger".format(virt), "a_{}".format(occ)) )
        single_excitations.append(operators)

    print()
    print("Single excitation operators")
    print(single_excitations)
    print()

    double_excitations = []
    double_excitations_flat = []
    # double_ex_arr = np.zeros((4, 4, 4, 4))

    num_occ_spin_orbitals = len(occupied_spin_orbitals)
    num_virt_spin_orbitals = len(virtual_spin_orbitals)

    for occ_idx in range(num_occ_spin_orbitals):
        annihilation_outer = []
        for occ_idx_runner in range(occ_idx+1, num_occ_spin_orbitals):
            annihilation_inner = []
            for virt_idx in range(num_virt_spin_orbitals):
                creation_outer = []
                for virt_idx_runner in range(virt_idx+1, num_virt_spin_orbitals):
                    # print(virtual_spin_orbitals[virt_idx_runner],
                    #       virtual_spin_orbitals[virt_idx],
                    #       occupied_spin_orbitals[occ_idx_runner],
                    #       occupied_spin_orbitals[occ_idx]
                    # )
                    # double_ex_arr[virtual_spin_orbitals[virt_idx_runner],
                    #               virtual_spin_orbitals[virt_idx],
                    #               occupied_spin_orbitals[occ_idx_runner],
                    #               occupied_spin_orbitals[occ_idx]] = 1
                    operator = ("a_{}^dagger".format(virtual_spin_orbitals[virt_idx_runner]),
                                "a_{}^dagger".format(virtual_spin_orbitals[virt_idx]),
                                "a_{}".format(occupied_spin_orbitals[occ_idx_runner]),
                                "a_{}".format(occupied_spin_orbitals[occ_idx])
                    )
                    adjoint_operator = ("a_{}^dagger".format(occupied_spin_orbitals[occ_idx]),
                                        "a_{}^dagger".format(occupied_spin_orbitals[occ_idx_runner]),
                                        "a_{}".format(virtual_spin_orbitals[virt_idx]),
                                        "a_{}".format(virtual_spin_orbitals[virt_idx_runner])
                    )
                    print(operator)
                    print(adjoint_operator)
                    creation_outer.append(operator)
                    double_excitations_flat.append(operator)
                    # print(" ".join(operator))
                if virt_idx+1 != num_virt_spin_orbitals:
                    annihilation_inner.append(creation_outer)
            if virt_idx != num_virt_spin_orbitals:
                annihilation_outer.append(annihilation_inner)
        if occ_idx+1 != num_occ_spin_orbitals:
            double_excitations.append(annihilation_outer)

    # for x in range(len(double_excitations)):
    #     print(double_excitations[x])
    #     print()

    print()
    print("Double excitation operators")
    for op in double_excitations_flat:
        print(" ".join(op))
    print()

    # print(double_ex_arr)


"""

"""
def fermion_hamiltonian_operator(num_occ_spin_orbitals):
    pass

UCCSD_operator("0011")
