#------------------------------------------------------------------------
# MC kernel for A-B  pairs in lattice
# only first nearest neighbor to be considered. Choose cutoff carefully.
# neighbor list output to @debug_neighbor_table.csv. Must check manually.
# Author:                                Lan Huang headhuanglan@tamu.edu
#                                                               2018,3,5
#------------------------------------------------------------------------
import time
import os
import math
import numpy as np
from scipy.spatial.kdtree import KDTree
#------------------------------------------------------------------------
#lattic class
#@param<float>   lattice_length    e.g: 1.0
#@param<list>    lattice_size      e.g: [10,10,10]
#@param<string>  lattice_type      e.g: 'fcc' ,'bcc'
#                IF lattice_type is customized
#                   lattice_type=np.vstack((basis1,basis2,basis...basisN))
#                             where basisX=np.array([x,y,z])
class lattice(object):
    def __init__(self,lattice_length,lattice_size,lattice_type='fcc',simulation_folder="."):
        self.simulation_folder=simulation_folder
        self.lattice_length=lattice_length
        self.lattice_size=lattice_size
        self.lattice_type= lattice_type
        self.atomtypesarray=None
        self.xyz=self.__setxyz()
        self.box=self.__setboxsize()
        self.natoms=len(self.xyz)
        if not os.path.exists(self.simulation_folder):
            os.makedirs(self.simulation_folder)

    def __basis(self):
        if self.lattice_type=='fcc':
            basis1 = np.array([0.0, 0.0, 0.0])
            basis2 = np.array([0.0, 0.5, 0.5])
            basis3 = np.array([0.5, 0.0, 0.5])
            basis4 = np.array([0.5, 0.5, 0.0])
            return np.vstack((basis1,basis2,basis3,basis4))
        elif self.lattice_type=='bcc':
            basis1 = np.array([0.0, 0.0, 0.0])
            basis2 = np.array([0.5, 0.5, 0.5])
            return np.vstack((basis1,basis2))
        elif self.lattice_type=='customized':
            basis1 = np.array([0.5, 0.5, 0.0])
            basis2 = np.array([0.0, 0.5, 0.5])
            basis3 = np.array([0.5, 0.0, 0.5])
            return np.vstack((basis1,basis2,basis3))
        else:
            return self.lattice_type

    def __setxyz(self):
        if not len(self.lattice_size)==3:
            raise Exception("lattice size must be L by M by N")
        basis=self.__basis()*self.lattice_length
        xyz=np.array([])
        first_flag=True
        for i in range(self.lattice_size[0]):
            for j in range(self.lattice_size[1]):
                for k in range(self.lattice_size[2]):
                    if first_flag:
                        xyz=basis+self.lattice_length*np.array([i,j,k])
                        first_flag=False
                    else:
                        xyz=np.vstack([xyz,basis+self.lattice_length*np.array([i,j,k])])
        self.natoms=len(xyz)
        return xyz


    def __setboxsize(self):
        return  self.lattice_length*np.asarray(self.lattice_size)

    def getxyz(self):
        return self.xyz

    def getboxsize(self):
        return self.box

    #foutput xyz to file for debug
    def writexyzfile(self,filename='@debug_lattice.xyz'):
        filename=os.path.join(self.simulation_folder, filename)
        natoms = self.natoms
        xyz = self.xyz
        atomtypesarray=self.atomtypesarray
        if  atomtypesarray is None:
            with open(filename,'w') as f:
                f.write("%d\n" % (natoms,))
                f.write("Lan Huang headhuanglan@tamu.edu\n")
                for atom_xyz in xyz:
                    f.write("C %f %f %f \n" % (atom_xyz[0],atom_xyz[1],atom_xyz[2]))
            print("Write lattice xyz info to file "+ filename)
        else:
            with open(filename, 'w') as f:
                f.write("%d\n" % (natoms,))
                f.write("Lan Huang headhuanglan@tamu.edu\n")
                for atomtype,atom_xyz in zip(atomtypesarray,xyz):
                    f.write("%s %f %f %f \n" % (atomtype, atom_xyz[0], atom_xyz[1], atom_xyz[2]))
            print("Write lattice xyz (with atom type info) to file " + filename)


    #reset latticelength
    def resetlatticelength(self,lattice_length):
        self.lattice_length=lattice_length
        #refresh info
        self.__setboxsize()
        self.__setxyz()

    def setatomtypesarrary(self,atomtypesarrary):
        if not len(atomtypesarrary)==self.natoms:
            raise Exception("len atomtypesarrary != natoms")
        self.atomtypesarray=atomtypesarrary

    def getatomtypesarrary(self):
        return self.atomtypesarray



#------------------------------------------------------------------------
#mc_kernel class
#param<int>      iteration
#param<lattice>  latticeOBJ
#param<list>     atomtypes
#param<list>     atomratio
#param<dict>     interaction
#param<int>      seed      random number genarator seed
class mc_kernel(object):
    def __init__(self,iteration,cutoff,latticeOBJ,atomtypes=['A','B'],atomratio=[0.5,0.5],interaction={'AA':0,'BB':0,'AB':4},KbT=1,simulation_folder=".",seed=12345):
        self.iteration=iteration
        self.cutoff=cutoff
        self.lattice=latticeOBJ
        self.lattice_length=latticeOBJ.lattice_length
        self.box=latticeOBJ.getboxsize()
        self.xyz=latticeOBJ.getxyz()
        self.natoms=len(self.xyz)
        self.neighbor_table=None
        self.atomtypes=atomtypes
        self.atomratio=atomratio
        self.interaction=interaction
        self.KbT=KbT
        self.simulation_folder=simulation_folder
        if not os.path.exists(self.simulation_folder):
            os.makedirs(self.simulation_folder)
        self.seed=seed
        np.random.seed(seed)

    # write xyz trajectroy will not create or close file
    # the only function is to write one frame
    def __writexyztrj_one(self, file_handler, atomtypesarrary):
        xyz=self.xyz
        natoms = len(xyz)
        if not len(atomtypesarrary)==natoms:
            raise Exception("Error| natoms != len atomtypes arrary")
        file_handler.write("%d\n" % (natoms,))
        file_handler.write("Lan Huang headhuanglan@tamu.edu\n")
        for atomtype, atom_xyz in zip(atomtypesarrary, xyz):
            file_handler.write("%s %f %f %f \n" % (atomtype, atom_xyz[0], atom_xyz[1], atom_xyz[2]))
    #pbc query one atom's neigher
    def __pbc_query(self,atomid,xyz_one,my_kdtree):
        cutoff=self.cutoff
        latticelength=self.lattice_length
        box=self.box
        neighborIDs = []
        mask=[[ 0, 0, 0], [ 1, 0, 0], [-1, 0, 0], [ 0, 1, 0], [ 0,-1, 0], [ 0, 0, 1], [ 0, 0,-1],
              [ 1, 1, 0], [ 0, 1, 1], [ 1, 0, 1], [-1,-1, 0], [ 0,-1,-1], [-1, 0,-1], [ 1,-1, 0],
              [-1, 1, 0], [ 0, 1,-1], [ 0,-1, 1], [ 1, 0,-1], [-1, 0, 1], [ 1, 1, 1], [-1, 1, 1],
              [ 1,-1, 1], [ 1, 1,-1], [-1,-1, 1], [ 1,-1,-1], [-1, 1,-1], [-1,-1,-1]]
        mask=np.array(mask)
        #xyz_list contain all posibble image
        xyz_list=xyz_one+mask*box
        #filter image in box, carefully dealing the wraped atoms
        mask_lt_boxhi=(xyz_list < (box+latticelength)).all(axis=1)
        mask_gt_boxlo=(xyz_list > (np.array([0,0,0])-latticelength)).all(axis=1)
        #get query points
        xyz_list=xyz_list[np.logical_and(mask_gt_boxlo,mask_lt_boxhi)]
        for xyz_one_in_list in xyz_list:
             neighborIDs.extend(my_kdtree.query_ball_point(xyz_one_in_list, cutoff))
        #remove ID of myself
        neighborIDs=list(set(neighborIDs)-set([atomid,]))
        return neighborIDs
    #this function query one atom's neigher and get the energy sum(E_ij)
    # i is the current atoms; j is the neigher; sum over all neighbor j
    def __energy_one(self,atomid,atomtypesarray):
        neighbor_table=self.neighbor_table
        interaction=self.interaction
        e_one=0
        atom_me_type=atomtypesarray[atomid]
        me_neighbor_atoms=neighbor_table[atomid]
        for neighbor_atom_id in me_neighbor_atoms:
            atom_nei_type = atomtypesarray[neighbor_atom_id]
            type1type2 = [atom_me_type, atom_nei_type]
            type1type2.sort()
            interactionkey = ''.join(type1type2)
            e_one += interaction[interactionkey]
        return  e_one
    #this function calculates the energy of the entire crystal normalized to per atom
    def __energy(self,atomtypesarray):
        natoms=self.natoms
        interaction=self.interaction
        neighbor_table=self.neighbor_table
        e=0
        for atomid in range(natoms):
            atom_me_type=atomtypesarray[atomid]
            me_neighbor_atoms=neighbor_table[atomid]
            for neighbor_atom_id in me_neighbor_atoms:
                atom_nei_type=atomtypesarray[neighbor_atom_id]
                type1type2=[atom_me_type, atom_nei_type]
                type1type2.sort()
                interactionkey=''.join(type1type2)
                e+=interaction[interactionkey]
        #A-B  pair  when calculate e. A->B   B->A  double counting e, so e=e/2
        e=e/2.0
        #normalize
        e=e/float(self.natoms)
        return e
    #this function build the neighbor table for atoms
    def __build_neighbor_table(self,debug=True):
        xyz=self.xyz
        my_kdtree = KDTree(xyz)
        neighbor_table={}
        if debug:
            FirstFlag=True
            csv_neighbor_table=None
        for atomid, xyz_one in enumerate(xyz):
            neighborIDs = self.__pbc_query(atomid, xyz_one, my_kdtree)
            neighbor_table[atomid]=neighborIDs
            if debug:
                if FirstFlag:
                    FirstFlag=False
                    csv_neighbor_table=np.array(neighborIDs)
                else:
                    csv_neighbor_table=np.vstack((csv_neighbor_table,np.array(neighborIDs)))
        if debug:
            np.savetxt(os.path.join(self.simulation_folder,"@debug_MC_neighbor_table.csv"), csv_neighbor_table, delimiter=",")
        return neighbor_table
    #this function switch the diffrent atom type's postion and retrun atomids
    #for efficency, only diffrent types of atoms are swaped.
    def __swap_atompair(self,atomtypesarray):
        natoms=self.natoms
        while True:
            atomid1,atomid2=np.random.randint(0,natoms,2)
            atom1type=atomtypesarray[atomid1]
            atom2type=atomtypesarray[atomid2]
            if not atom1type==atom2type:
                break
        #cal e_old
        e_old=self.__energy_one(atomid1,atomtypesarray)+self.__energy_one(atomid2,atomtypesarray)
        #do swap
        atomtypesarray[atomid1]=atom2type
        atomtypesarray[atomid2]=atom1type
        #cal e_new
        e_new=self.__energy_one(atomid1,atomtypesarray)+self.__energy_one(atomid2,atomtypesarray)
        #deltaE
        deltaE=e_new-e_old
        #normalize delatE to per atom value
        deltaE=deltaE/float(self.natoms)
        return atomid1,atomid2,atomtypesarray,deltaE
    #this function build inital atomtypesarray
    def __build_atomtypesarrary(self):
        natoms = self.natoms
        natomtypes = len(self.atomtypes)
        atomtypeslist = []
        for i, type in enumerate(self.atomtypes):
            if i == natomtypes - 1:
                n = natoms - len(atomtypeslist)
                atomtypeslist.extend([type] * n)
            else:
                n = int(self.atomratio[i] * natoms)
                atomtypeslist.extend([type] * n)
        atomtypesarray = np.asarray(atomtypeslist)
        np.random.shuffle(atomtypesarray)
        return atomtypesarray


    #this is the main function for mc kernel
    def simulate(self,outputfrq_screen=10,outputfrq_file=10,trjfilename='trj.xyz',logfilename='MC_E.log'):
        natoms = self.natoms
        #initial setup, generate randomly distributed atoms types
        atomtypesarray=self.__build_atomtypesarrary()
        #build neighbor table
        self.neighbor_table=self.__build_neighbor_table()
        #get inital configuration's total_energy per atom
        estart=self.__energy(atomtypesarray)
        print("estart per atom:%f \n" % (estart,))
        elast=estart
        naccept = 0
        # trajectory output
        file_handler = open(os.path.join(self.simulation_folder,trjfilename), 'w')
        log_headler= open(os.path.join(self.simulation_folder,logfilename),'w')
        for i in range(self.iteration):
            #per iteration loop over all atoms
            for _ in range(self.natoms):
                # try to switch atomtypes of pair atoms
                atomid1,atomid2,atomtypesarray,deltaE=self.__swap_atompair(atomtypesarray)
                #if deltaE<0 accectp move
                if deltaE<0:
                    elast+=deltaE
                    naccept += 1
                elif np.random.random() <= math.exp(natoms*(-deltaE)/self.KbT):
                    elast += deltaE
                    naccept += 1
                else:
                    #swap back to original atomtypesarry
                    atomtypesarray[atomid1],atomtypesarray[atomid2]= atomtypesarray[atomid2],atomtypesarray[atomid1]
            if i%outputfrq_screen ==0:
                print("%d %f\n" % (i,elast))
            if i%outputfrq_file == 0:
                self.__writexyztrj_one(file_handler, atomtypesarray)
                log_headler.writelines("%d %f\n" % (i,elast))
        print("Naccept")
        print(naccept)
        print("biased accept percent(only swap diffrent atom type)")
        print(naccept/(self.iteration*natoms))
        print("e last")
        print(elast)
        print("e end")
        print(self.__energy(atomtypesarray))
        file_handler.close()
        log_headler.close()
        #save the final atomtypesarray
        self.atomtypesarray=atomtypesarray

    def getAtomTypesArray(self):
        return self.atomtypesarray

#------------------------------------------------------------------------
#forster_kernel class
#param<lattice>  latticeOBJ
#param<real>     R0
#param<real>     rcut
class forster_kernel(object):
    def __init__(self,R0,rcut,latticeOBJ,simulation_folder="."):
        self.latticeOBJ=latticeOBJ
        self.atomtypesarray=latticeOBJ.getatomtypesarrary()
        self.xyz=latticeOBJ.getxyz()
        self.box=latticeOBJ.box
        self.natoms=latticeOBJ.natoms
        self.lattice_length=latticeOBJ.lattice_length
        self.R0=R0
        self.rcut=rcut
        self.simulation_folder=simulation_folder
        self.Iaverage=None
        if not os.path.exists(self.simulation_folder):
            os.makedirs(self.simulation_folder)

    # pbc query one atom's neigher
    def __pbc_query(self, atomid, xyz_one, my_kdtree):
        cutoff = self.rcut
        latticelength = self.lattice_length
        box = self.box
        neighborIDs = []
        mask = [[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1],
                [1, 1, 0], [0, 1, 1], [1, 0, 1], [-1, -1, 0], [0, -1, -1], [-1, 0, -1], [1, -1, 0],
                [-1, 1, 0], [0, 1, -1], [0, -1, 1], [1, 0, -1], [-1, 0, 1], [1, 1, 1], [-1, 1, 1],
                [1, -1, 1], [1, 1, -1], [-1, -1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, -1]]
        mask = np.array(mask)
        # xyz_list contain all posibble image
        xyz_list = xyz_one + mask * box
        # filter image in box, carefully dealing the wraped atoms
        mask_lt_boxhi = (xyz_list < (box + np.ceil(cutoff/latticelength)*latticelength)).all(axis=1)
        mask_gt_boxlo = (xyz_list > (np.array([0, 0, 0]) - np.ceil(cutoff/latticelength)*latticelength)).all(axis=1)
        # get query points
        xyz_list = xyz_list[np.logical_and(mask_gt_boxlo, mask_lt_boxhi)]
        for xyz_one_in_list in xyz_list:
            neighborIDs.extend(my_kdtree.query_ball_point(xyz_one_in_list, cutoff))
        # remove ID of myself
        neighborIDs = list(set(neighborIDs) - set([atomid, ]))
        return neighborIDs

    def __build_neighbor_table(self, debug=True):
        xyz = self.xyz
        my_kdtree = KDTree(xyz)
        neighbor_table = {}
        if debug:
            FirstFlag = True
            csv_neighbor_table = None
        for atomid, xyz_one in enumerate(xyz):
            neighborIDs = self.__pbc_query(atomid, xyz_one, my_kdtree)
            neighbor_table[atomid] = neighborIDs
            if debug:
                if FirstFlag:
                    FirstFlag = False
                    csv_neighbor_table = np.array(neighborIDs)
                else:
                    csv_neighbor_table = np.vstack((csv_neighbor_table, np.array(neighborIDs)))
        if debug:
            np.savetxt(os.path.join(self.simulation_folder,"@debug_Forster_neighbor_table.csv"), csv_neighbor_table, delimiter=",")
        return neighbor_table


    def __distance(self,atomid1,atomid2):
        xyz=self.xyz
        box=self.box
        halfbox=self.box/2.0
        xyz_one_1=xyz[atomid1]
        xyz_one_2=xyz[atomid2]
        xyz_1_2=xyz_one_1 - xyz_one_2
        distance=np.sqrt(np.sum(xyz_1_2**2))
        if distance>self.rcut:
            abs_xyz_1_2=np.abs(xyz_1_2)
            mask=abs_xyz_1_2 > halfbox
            xyz_1_2= box*mask-abs_xyz_1_2
            distance = np.sqrt(np.sum(xyz_1_2 ** 2))
            if distance>self.rcut:
                raise Exception("Error!!! Somthing wrong!!!")
        return  distance


    #output xyz and I to file
    def __writexyzfile(self,I,filename='result.xyz'):
        filename=os.path.join(self.simulation_folder,filename)
        natoms = self.natoms
        xyz = self.xyz
        atomtypesarray=self.atomtypesarray
        with open(filename, 'w') as f:
            f.write("%d\n" % (natoms,))
            f.write("Lan Huang headhuanglan@tamu.edu\n")
            for atomtype,atom_xyz,ivalue in zip(atomtypesarray,xyz,I):
                f.write("%s %f %f %f %f\n" % (atomtype, atom_xyz[0], atom_xyz[1], atom_xyz[2], ivalue))
        print("Write I value to xyz file: " + filename)

    #Do the forster energy transfer calculation
    #filter neighbor to select A B pairs only
    # Förster resonance energy transfer efficiency
    # E=R0^6/(r^6+R0^6)    I_DA/ID=1-E
    # E decay fast, If rcut=3*R0  E= 0.137%
    # r=0 E=Emax=1    r= (R0^6/E-R0^6)^(1/6)  = R0*(1/E-1)^(1/6)
    def simulate(self):
        self.neighbor_table = self.__build_neighbor_table()
        atomtypesarrary=self.atomtypesarray
        R0=self.R0
        # E=R0^6/(r^6+R0^6)    I_DA/ID=1-E
        E=np.zeros(shape=(self.natoms,1))
        I = np.zeros(shape=(self.natoms, 1))
        distancefun=np.vectorize(self.__distance)
        R0_6=np.power(R0,6)
        for atomid in range(self.natoms):
            neighbor_list=np.array(self.neighbor_table[atomid])
            atomtype = atomtypesarrary[atomid]
            neighbortype=atomtypesarrary[neighbor_list]

            #filter neighbor to select A B pairs only
            filter_neighbor_list=neighbor_list[neighbortype != atomtype]
            if not len(filter_neighbor_list)==0:
                distances=distancefun(atomid,filter_neighbor_list)
                #Förster resonance energy transfer efficiency
                Ei=R0_6 / (np.power(distances, 6) + R0_6)
                I[atomid]=np.multiply.reduce(1-Ei)

        self.I=I
        self.Iaverage=np.sum(self.I)/self.natoms
        self.__writexyzfile(I)

        #output I average
        filename = os.path.join(self.simulation_folder, "I.txt")
        with open(filename, 'w') as f:
            f.write("%f\n" % (self.Iaverage,))

        print("I averaged over all atoms:%f\n" % (self.Iaverage,))
        print("vmd -e vmd.tcl")
        print("topo readvarxyz trj.xyz")

    def getIaverage(self):
        if self.Iaverage is None:
            raise Exception("Error! please run froster kernel simulate() first")
        else:
            return self.Iaverage





def test():


    #setup run--------------------------------------------------------------------------------------------------------
    #lattic MC unitless simulation

    KbT=1
    seed=12345
    lattice_size=[10,10,10]
    lattice_type = 'customized'
    lattice_length = 1
    mc_1st_neigh_cutoff=math.sqrt(2.0)/2.0
    mc_iteration=1#000
    atomtypes = ['A', 'B']
    atomratio = [0.5, 0.5]
    interaction = {'AA': 0, 'BB': 0, 'AB': 2}
    mc_outputfrq_screen = 10
    mc_outputfrq_file = 10

    #forster
    lattice_real = 27.6517/2.0
    R0 = 30.55  # 32.2~28.9
    rcut = 2*R0#3 * R0

    #folder setup
    simulation_folder = "simulation_data/"+str(atomratio[0])+"_"+str(atomratio[1])+"AB_"+str(interaction['AB'])


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


if __name__=='__main__':
    test()

