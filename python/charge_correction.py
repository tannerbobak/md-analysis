"""

~=~=~=~=~=~=~=~=~=~=~ Charge Correction Solver ~=~=~=~=~=~=~=~=~=~=~

Solves for MD charge corrections as described in:
"Calculating the binding free energies of charged species based on 
explicit-solvent simulations employing lattice-sum methods: An 
accurate correction scheme for electrostatic finite-size effects" by
GJ Rocklin et al.

"""
import math


# Parameters
epsilon_S		=				   		   97 # dimensionless, ratio to vacuum permittivity
xi_CB			= 					-2.380077 # dimensionless
xi_LS			= 					-2.837397 # dimensionless
epsilon_0		=      8.854187 * (10 ** -12) # Farads per meter, vacuum permittivity
pi				=                    3.141592 # ...6535897932... duh
gamma_S			= 					   0.0764 # e*nm^2

# Inputs
Q_P 	=      int(input('Q_P: ')) # e
Q_L 	=      int(input('Q_L: ')) # e
L  	    =      float(input('L: ')) # nm
I_P		=    float(input('I_P: ')) # Protein RIP
I_L		=    float(input('I_L: ')) # Ligand RIP
I_LSLV  = float(input('I_LSLV: ')) # Solvent RIP, raw from APBS
N_S		=      int(input('N_S: ')) # Number of solvent molecules
G_MDPBC =  int(input('G_MDPBC: ')) # Binding free energy with periodic boundary conditions. kJ/mol

# Periodicity-induced net-charge interactions
G_NET = ( -(xi_LS)/(8*pi*epsilon_0) ) * ( (Q_P + Q_L)**2 - Q_P**2 ) * ( 1/L )

# Periodicity-induced undersolvation
G_USV = ( (xi_LS)/(8*pi*epsilon_0) ) * ( 1 - (1/epsilon_S) ) * ( (Q_P + Q_L)**2 - Q_P**2 ) * ( 1/L )

# Residual integrated potential
G_RIP = ( (I_P + I_L)*(Q_P + Q_L) - I_P*Q_P ) * ( 1/(L**3) )

# Empirical correction term
I_LSLV = I_L - I_LSLV
R_L = math.sqrt( ( (1/(6*epsilon_0)) * (1-(1/epsilon_S)) * Q_L )**(-1) * I_LSLV )
G_EMP = (-(2*pi)/(45*epsilon_0)) * (1-(1/epsilon_S)) * ( (Q_P + Q_L)**2 - Q_P**2 ) * ((R_L**5)/(L**6))

# Discrete solvent effects
G_DSC = -(gamma_S*Q_L*N_S)/(6*epsilon_0*L**3)

# FINALLY
G_MDNBC = G_MDPBC + G_NET + G_USV + G_RIP + G_EMP + G_DSC
print('(delta)G_{NET}='+ str(G_NET) + ' kJ/mol')
print('(delta)G_{USV}='+ str(G_USV) + ' kJ/mol')
print('(delta)G_{RIP}='+ str(G_RIP) + ' kJ/mol')
print('(delta)G_{EMP}='+ str(G_EMP) + ' kJ/mol')
print('(delta)G_{DSC}='+ str(G_DSC) + ' kJ/mol')
print('=======================')
print('(delta)G_{MD, NBC}='+ str(G_MDNBC) + ' kJ/mol')