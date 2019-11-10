import psi4
import numpy as np

# ==> Set Basic Psi4 Options <==
# Memory specification
psi4.set_memory(int(5e8))
numpy_memory = 2

# ==> Set default program options <==
# Maximum SCF iterations
MAXITER = 100
# Energy convergence criterion
E_conv = 1.0e-8

psi4_options = {"basis": "sto-3g",
# psi4_options = {"basis": "6-31g**",
                "scf_type": "pk",
                "e_convergence": 1e-8}

psi4.set_options(psi4_options)

# Set output file
psi4.core.set_output_file('output.dat', False)

H2_molecule_str = """
    units bohr
    H
    H 1 {bond_separation}
"""

HeH_molecule_str = """
    units bohr
    0 1
    He
    --
    1 1
    H 1 {bond_separation}
"""

Li2_molecule_str = """
    Li
    Li 1 {bond_separation}
"""

# Li_molecule_str = """
#     +1 1
#     Li
# """

molecule_strings = {}
molecule_strings["H2"] = H2_molecule_str
molecule_strings["HeH"] = HeH_molecule_str
molecule_strings["Li2"] = Li2_molecule_str
# molecule_strings["Li"] = Li_molecule_str


def get_mol(molecule_name, bond_separation):
    return psi4.geometry(molecule_strings[molecule_name].format(bond_separation=bond_separation))


mol = get_mol("H2", 1.401)

# ==> Compute static 1e- and 2e- quantities with Psi4 <==
# Class instantiation
wavefunction = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option("basis"))
mints = psi4.core.MintsHelper(wavefunction.basisset())

# Overlap matrix
S = np.asarray(mints.ao_overlap())
print(S)


# Number of basis Functions & doubly occupied orbitals
num_basis_functions = S.shape[0]
num_doubly_occupied_orbitals = wavefunction.nalpha()
print("Number of occupied orbitals: %3d" % (num_doubly_occupied_orbitals))
print("Number of basis functions: %3d" % (num_basis_functions))

# Memory check for ERI tensor
I_size = (num_basis_functions**4) * 8.e-9
print("Size of the ERI tensor will be {:4.2f} GB.".format(I_size))
memory_footprint = I_size * 1.5
if I_size > numpy_memory:
    psi4.core.clean()
    raise Exception("Estimated memory utilization (%4.2f GB) exceeds allotted memory \
                     limit of %4.2f GB." % (memory_footprint, numpy_memory))

# Electron Repulsion Interaction
I = np.asarray(mints.ao_eri())
print("ERI tensor shape:", I.shape)

# 1 electron energies
# kinetic energy
T = np.asarray(mints.ao_kinetic())
# potential energy
V = np.asarray(mints.ao_potential())

print(T)
print()
print(V)
print()
print(T+V)
print()
print(I)

one_body_integrals = T+V
two_body_integrals = I

save_dir = "../molecules/H2/"
one_body_fname = "H2_one_body.npy"
two_body_fname = "H2_two_body.npy"

def write_integrals_to_file():
    np.save(save_dir + one_body_fname, one_body_integrals, allow_pickle=False)
    np.save(save_dir + two_body_fname, two_body_integrals, allow_pickle=False)

    # loadtest_one = np.load(save_dir + one_body_fname)
    # loadtest_two = np.load(save_dir + two_body_fname)

    # print()
    # print(loadtest_one)
    # print()
    # print(loadtest_two)
    #
    # print(np.equal(one_body_integrals, loadtest_one).all())
    # print(np.equal(two_body_integrals, loadtest_two).all())

# for i in range(I.shape[0]):
#     for j in range(I.shape[1]):
#         for m in range(I.shape[2]):
#             for n in range(I.shape[3]):
#                 # print(i+1, j+1, m+1, n+1, I[i,j,m,n])
#                 print(i, j, m, n, '\t', I[i,j,m,n])



# idxs = [(0, 0), (1, 0), (1, 1)]
# idxs = [(0, 1), (1, 0)]
# idxs = [(1, 0), (0, 1)]
# for i in range(len(idxs)):
#     for j in range(i+1):
#         print(i, j)
#         print(idxs[i], idxs[j])
#         print(idxs[i], idxs[j], I[idxs[i][0], idxs[i][1], idxs[j][0], idxs[j][1]])
#         print(idxs[j], idxs[i], I[idxs[j][0], idxs[j][1], idxs[i][0], idxs[i][1]])

# paulis = ['I', 'X', 'Y', 'Z']


# def creation_op(loc, nbf):
#     op_list = [[], []]
#     for idx in range(nbf):
#         if idx < loc:
#             op_list[0].append(paulis[3])
#             op_list[1].append(paulis[3])
#         elif idx == loc:
#             op_list[0].append(paulis[1])
#             op_list[1].append('-'+paulis[2])
#         else:
#             op_list[0].append(paulis[0])
#             op_list[1].append(paulis[0])
#     return op_list
#
#
# def annihilation_op(loc, nbf):
#     op_list = [[], []]
#     for idx in range(nbf):
#         if idx < loc:
#             op_list[0].append(paulis[3])
#             op_list[1].append(paulis[3])
#         elif idx == loc:
#             op_list[0].append(paulis[1])
#             op_list[1].append(paulis[2])
#         else:
#             op_list[0].append(paulis[0])
#             op_list[1].append(paulis[0])
#     return op_list
#
#
# num_fermions = 2
#
# for i in range(I.shape[0]):
#     for j in range(I.shape[1]):
#         for m in range(I.shape[2]):
#             for n in range(I.shape[3]):
#                 print(i, j, m, n, '\t', I[i,j,m,n])
#                 print(annihilation_op(n, nbf))
#                 print(annihilation_op(n, nbf))
#                 print(creation_op(j, nbf))
#                 print(creation_op(i, nbf))
