import os
import math
header="""
#------------------------------------------------------------------------
# Perform grid search for hyper parameters
# key kernels are mc_kernel for lattice MC
#                 forster_kernel for forster energy tranfer
# Author:                                Lan Huang headhuanglan@tamu.edu
#                                                               2018,3,11
#------------------------------------------------------------------------

from MC_Forster import lattice as lattice
from MC_Forster import mc_kernel as mc_kernel
from MC_Forster import forster_kernel as forster_kernel
import math

def grid_search(repeat_counter,
                seed,
                atomratio,
                interaction,
                mc_outputfrq_screen,
                lattice_type,
                lattice_length,
                KbT,
                mc_1st_neigh_cutoff,
                atomtypes,
                mc_iteration,
                mc_outputfrq_file,
                lattice_size,
                lattice_real,
                R0,
                rcut,
                cluster_shape
                ):


    # folder setup
    simulation_folder = "simulation_data/"+cluster_shape+"/" + str(atomratio[0]) + "_" + str(atomratio[1]) + "AB_" + str(
        interaction['AB']) + "r_" + str(repeat_counter)


    #----------------------------------------------------------------------------------------------------------------
    #Part1:unitless MC simulation
    #test1 lattice class
    my_lattice= lattice(lattice_length=lattice_length,lattice_size= lattice_size,lattice_type=lattice_type,simulation_folder=simulation_folder)
    print("box size")
    print(my_lattice.getboxsize())
    print("lattice xyz")
    print(my_lattice.getxyz())
    print("lattic  write2file")
    my_lattice.writexyzfile()
    #test2 mc_kernel class    set cutoff of KDTree to find the 1st nerest neighbor
    my_mc_kernel=mc_kernel(iteration=mc_iteration,atomtypes=atomtypes,atomratio=atomratio,interaction=interaction,
                           cutoff= mc_1st_neigh_cutoff,latticeOBJ=my_lattice,KbT=KbT,simulation_folder=simulation_folder,seed=seed)
    my_mc_kernel.simulate(mc_outputfrq_screen,outputfrq_file=mc_outputfrq_file)
    print("final atom types arrary")
    print(my_mc_kernel.getAtomTypesArray())
    atomtypesarrary=my_mc_kernel.getAtomTypesArray()
    #finish unitless simulation

    #------------------------------------------------------------------------------------------------------------------
    #part2:real unit simulation
    #start real unit simulation  unit:A
    my_real_lattice=lattice(lattice_length=lattice_real,lattice_size= lattice_size,lattice_type=lattice_type,simulation_folder=simulation_folder)
    my_real_lattice.setatomtypesarrary(atomtypesarrary)
    my_real_lattice.writexyzfile()
    print("real lattice simulation box size (unit A)")
    print(my_real_lattice.getboxsize())
    my_forster_kernel=forster_kernel(R0=R0,rcut=rcut,latticeOBJ=my_real_lattice,simulation_folder=simulation_folder)
    my_forster_kernel.simulate()



    return my_forster_kernel.getIaverage()

"""

if __name__=='__main__':
    repeat=5
    randomseeds=[345343,123435,43543,23355,546456]
    mc_outputfrq_screen = 10
    lattice_type = 'customized'
    lattice_length = 1
    KbT = 1
    mc_1st_neigh_cutoff = math.sqrt(2.0) / 2.0
    atomtypes = ['A', 'B']
    # mc
    mc_iteration = 10000
    mc_outputfrq_file = 10
    lattice_size = [15, 15, 15]
    ratio_of_A=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1]
    atomratio_list=[[rationA,1.0-rationA] for rationA in ratio_of_A]
    interactionAB={"LargeClusters":{'AA':0.0,'BB':0.0,'AB':2.0},
                   "SmallClusters":{'AA':0.0,'BB':0.0,'AB':0.25},
                   "SmallClusters1":{'AA':0.0,'BB':0.0,'AB':0.5},
                   "SmallClusters2":{'AA':0.0,'BB':0.0,'AB':0.75},
                   "SmallClusters3":{'AA':0.0,'BB':0.0,'AB':1.0},
                   "SmallClusters4":{'AA':0.0,'BB':0.0,'AB':1.25},
                   "SmallClusters5":{'AA':0.0,'BB':0.0,'AB':1.5},
                   "SmallClusters6":{'AA':0.0,'BB':0.0,'AB':1.75},
                   "Random":{'AA':0.0,'BB':0.0,'AB':0.0},
                   "Alternating":{'AA':0.0,'BB':0.0,'AB':-2.0}}
    cluster_shapes=["LargeClusters","SmallClusters","SmallClusters1","SmallClusters2","SmallClusters3","SmallClusters4","SmallClusters5","SmallClusters6","Random","Alternating"]
    #forster
    lattice_real = 27.6517 / 2.0
    R0 = 30.55  # 32.2~28.9
    rcut = 2.5 * R0
    assert(lattice_real*min(lattice_size)/2.0 >rcut)
    #------------------------------------------------------------------
    file_counter=0
    scan_infile_folder="Scan_input"
    for atomratio in atomratio_list:
        for cluster_shape in cluster_shapes:
            for repeat_counter in range(repeat):
                #result_key='_'.join([cluster_shape,str(atomratio[0]),str(atomratio[1]),str(repeat_counter)])
                file_counter+=1
                filename=str(file_counter) + ".py"
                filename=os.path.join(scan_infile_folder,filename)
                with open(filename,"w") as f:
                    f.writelines(header)
                    f.writelines("grid_search(\n"+
                                 "repeat_counter="+str(repeat_counter)+",\n"+
                                 "seed="+str(randomseeds[repeat_counter])+",\n"+
                                 "atomratio="+str(atomratio,)+",\n"+
                                 "interaction="+str(interactionAB[cluster_shape])+",\n"+
                                 "mc_outputfrq_screen="+str(mc_outputfrq_screen)+",\n"+
                                 "lattice_type=\""+str(lattice_type)+"\",\n"+
                                 "lattice_length="+str(lattice_length)+",\n"+
                                 "KbT="+str(KbT)+",\n"+
                                 "mc_1st_neigh_cutoff="+str(mc_1st_neigh_cutoff)+",\n"+
                                 "atomtypes="+str(atomtypes)+",\n"+
                                 "mc_iteration="+str(mc_iteration)+",\n"+
                                 "mc_outputfrq_file="+str(mc_outputfrq_file)+",\n"+
                                 "lattice_size="+str(lattice_size)+",\n"+
                                 "lattice_real="+str(lattice_real)+",\n"+
                                 "R0="+str(R0)+",\n"+
                                 "rcut="+str(rcut)+",\n"+
                                 "cluster_shape=\""+str(cluster_shape)+"\")"

                                 )


                # grid_search(
                # repeat_counter=repeat_counter,
                # seed=randomseeds[repeat_counter],
                # atomratio=atomratio,
                # interaction=interactionAB[cluster_shape],
                # mc_outputfrq_screen=mc_outputfrq_screen,
                # lattice_type=lattice_type,
                # lattice_length=lattice_length,
                # KbT=KbT,
                # mc_1st_neigh_cutoff=mc_1st_neigh_cutoff,
                # atomtypes=atomtypes,
                # mc_iteration=mc_iteration,
                # mc_outputfrq_file=mc_outputfrq_file,
                # lattice_size=lattice_size,
                # lattice_real=lattice_real,
                # R0=R0,
                # rcut=rcut,
                # cluster_shape=cluster_shape)









