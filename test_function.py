import glob
species = "./species/"
def Read_Species_List(pr=1):
    global species
    """ If pr is 1, it will print "Species_List". """ 
    read_species =glob.glob(species+"*")   
    selected_species_dic = {}
    backward_selected_species_dic = {}
    number = 0
    for i, species in enumerate(sorted(read_species), start=1):
        selected_species_dic[i] = species.split('/')[-1]
        backward_selected_species_dic[species.split('/')[-1]] = i 
        if pr == 1 :
            print (str(i)+".", species.split('/')[-1])
        number = i
    return selected_species_dic, backward_selected_species_dic, number

print(Read_Species_List())