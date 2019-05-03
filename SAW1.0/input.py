number_density=0.1
bond_length=1.53
random_seed=9031
mass=14.02
chain_length=[1000]*8
#chain_length=[100,300,500,600,4000]

chain_length.sort(reverse=True)
number_chains=len(chain_length)

input_file=open('input.txt','w')
input_file.write('Title \n')
input_file.write('#Number_Density \n')
input_file.write(str(number_density)+ '\n')
input_file.write('#Bond_Length \n')
input_file.write(str(bond_length)+ '\n')
input_file.write('#Seed \n')
input_file.write(str(random_seed)+ '\n')
input_file.write('#Mass \n')
input_file.write(str(mass)+ '\n')
input_file.write('#number_chains \n')
input_file.write(str(number_chains)+ '\n')
input_file.write('#chain_length \n')
for i in range(number_chains):
    input_file.write(str(chain_length[i])+ '\n')
input_file.write('\n')
input_file.close()
print('input.txt has been generated \n')
print('check the file first \n')
print('press enter to generate lammps data file\n')
import os
os.popen("input.txt")
input()
import subprocess
subprocess.call('SAW.exe')  


