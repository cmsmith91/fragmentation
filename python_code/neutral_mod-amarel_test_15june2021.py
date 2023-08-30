import numpy as np
import random
# import os
#import multiprocessing as mp
import itertools
#import ray
#import line_profiler
import multiprocessing as mp

# wd='/Users/colleen/Dropbox/Fall_2014/Research/Fragmentation/data/'
# #wd='/Documents/Colleen/'
# os.chdir(wd)

# new_dir='neutralmod_output_25may2021'
# if os.path.exists(new_dir)==False:
#     os.makedirs(new_dir)
    
# os.chdir(new_dir)


def get_neighbors(size_xy,density):

    '''
    this function returns an array with the indices of all individuals' 
    nearest neighbors when the community occupies a landscape of 
    'size_xy'^2 pixels, with the number of individuals per pixel given by 
    'density'.
    
    the nearest neighbors of an individual are defined as all the individuals
    in the 8 pixels surrounding the individual's pixel, plus the indices of 
    all the individuals in that individual's pixel (so 9 pixels total)
    
    in the 'neutral' function below this one, the 3d community 
    (size_xy*sizee_xy*density) is represented by a flattened, 1d array. 
    however, for this get_neighbors function, i use the 
    3d indices bc it's easier
    
    this function starts with the indices in the 1d list, 
    converts them to 3d indices to get the neighbors and then back to 
     1d indices 
    '''
    
    # first, get 1d indices for each individual
    individual_indices=np.arange(size_xy**2*density) 
    
    # from these, get the indices of pixels in 1d list
    pixel_indices=(individual_indices/density).astype(int)
    
    #then get indices of rows and columns for each pixel
    i=pixel_indices%size_xy
    j=(pixel_indices/size_xy).astype(int) 
        
    
    # in the for loop, get these 
    # [i,j]
    # [i+1,j]
    # [i-1,j]
    # [i,j+1]
    # [i,j-1]
    # [i+1,j+1]
    # [i+1,j-1]
    # [i-1,j+1]
    # [i-1,j-1]

    # then if the new i or j is negative or > size_xy, get the neighbor 
    # on other side of matrix instead
    # then  converts the i,j index back to index of pixel in flattened array
   
    nbr_indices_pixels=[]

    for add_me_i in range(-1,2):
        for add_me_j in range(-1,2):
            
            new_coords=np.array([i+add_me_i,j+add_me_j]).transpose()
            
            # if new is or js are negative or greater than size_xy then 
            # [wrap around] to other side of the matrix by taking the remainder
            newer_coords=new_coords%size_xy
    
            #flatten matrix: convert is and js to index in 1d array of pixels
            new_pixel_index=newer_coords[:,0]+size_xy*newer_coords[:,1]
            
            #add to list
            nbr_indices_pixels.append(new_pixel_index)

           
    nbr_indices_t=np.array(nbr_indices_pixels).transpose()
    
    # right now indices are for list of pixels 
    # convert  the pixel indices to indices of individuals in 1d array 
    # assuming there the number of individuals per pixel is 'density'
    nbr_indices_inds=[]
    for nbrs in nbr_indices_t:
        nbr_indices_inds.append([pixel_index*density + k for 
                                 pixel_index in nbrs for k in range(density)])

    return(np.array(nbr_indices_inds))
def get_immediate_neighbors(size_xy,density):

    '''
    this function returns an array with the indices of all individuals' 
    immediate neighbors when the community occupies a landscape of 
    'size_xy'^2 pixels, with the number of individuals per pixel given by 
    'density'.
    
    the immediate neighbors of an individual are defined as all the individuals
    in the same pixel as the individual 
    '''
    
    # first, get 1d indices for each individual
    individual_indices=np.arange(size_xy**2*density) 
    
    # from these, get the indices of pixels in 1d list
    pixel_indices=(individual_indices/density).astype(int)
    

    nbr_indices_inds=[]

    for pixel_index in pixel_indices:
        nbr_indices_inds.append([pixel_index*density + k for k in range(density)])
   

    return(np.array(nbr_indices_inds))

# speciation,dispersal,replicate,n_simulations,print_every,size_xy,density=1e-7,1e-5,10,1000,1000,110,13
def neutral(speciation,dispersal,replicate,n_simulations,print_every,size_xy,density):
    #speciation,dispersal,replicate,n_simulations,print_every,size_xy,density=param_list
    '''
    this function returns a community located in a landscape with 
    size_xy^2 pixels and the number of individuals per pixel given by 
    'density'.
    
    with each of n_simulations, all individuals 
    die and are replaced either by a speciation event with probability 
    speciation, or by dispersal. There's 2 types of dispersal: long-distance 
    and nearest neighbor dispersal. 
    long-distance dispersal occurs with probability 'dispersal' if 
    speciation doesn't occur. if neither speciation or long-distance
    dispersal occur then nearest neighbor dispersal occurs
    
    '''
    
    size=(size_xy**2)*density
    
    #make an array with the indices each individual's neighbors
    nbr_indices=get_immediate_neighbors(size_xy,density)
     
    #start off the simulation with all individuals 
    #being from the same species. 
    #Each species is a different number

    species = np.array([1]*(size))
    species_id = 2
    
    #make the file name for downloading as a csv. This will have
    #all the parameter info for the simulation
    path_name = "spec{}_global{}_sizexy{}_density{}_rep{}".format(speciation,
                       dispersal,size_xy,density,replicate)

    array_falses=np.array([False]*size)
    #loop through all the timesteps
    for j in range(n_simulations):
        
        
        #make copy of species at previous timestep
        species_previous = np.copy(species)
        
        ## first, see where speciation, long-distance dispersal and 
        ## nearest neighbor dispersal happen
    
        # draw randomly from uniform distribution to see where speciation 
        # happens
        where_speciation=np.random.uniform(size=size)<speciation
        numb_new_species=np.count_nonzero(where_speciation)

        
        # draw randomly from uniform distribution to see where long-distance 
        # dispersal happens
        #where_long_distance=np.copy(where_speciation)
        where_long_distance=np.copy(array_falses)
        where_long_distance[where_speciation==False]=np.random.uniform(size=(size-numb_new_species))<dispersal 
        
        #np.random.binomial(
         #   n=1, p=dispersal, size=size) * (where_speciation ==0)
        # nearest-neighbor dispersal happens in all the other locations, 
        everything_else=where_speciation+where_long_distance==0
    
        
        ## next, each individual dies and gets replaced via a speciation event
        # or by dispersal
        
        #1) speciation happens. each new species is a new number
        species[where_speciation==True]=range(
            species_id,species_id+numb_new_species)
        species_id+=numb_new_species
    
        #2) long-distance dispersal happens: random individual chosen from 
        #anywhere
        species[where_long_distance]=random.choice(species_previous)
    
        #3) nearest neighbor dispersal happens
        how_many=np.count_nonzero(everything_else)
        
        # randomly pick which index you'll use to get parent index
        nbr_choices=np.random.choice(range(density), how_many) 
        # then get the actual parent index using array of nbr indices
        # indices_parents1=nbr_indices[everything_else] 
        indices_parents=nbr_indices[range(how_many),nbr_choices] 

        # then replace with individual at the parent index
        species[everything_else]=species_previous[indices_parents] 
        
        #save the new community as a txt file if 
        #it's the first or last round of the simulation
        #or if the timestep is a multiple of print_every       
        if j==0 or ((j+1) % print_every) ==0 or j == (n_simulations-1):            
            
            #add how many timesteps to the file name
            the_path = path_name+"timesteps{}.txt".format((j+1))
            
            #save species list to a txt file
            with open(the_path, 'w') as f:
                for item in species:
                    f.write("%s\n" % item)
                f.close()
    #return(species)
#neutral(1e-7,1e-5,10,1000,1000,110,13)

# import timeit
# %timeit neutral(1e-7,1e-5,10,10,10,110,13)
# %timeit where_long_distance=np.copy(where_speciation)
# %timeit np.array([False]*size)
#make a df of the parameter values i'm checking 
#iteration=0

# all_specs = [1.7e-5,  2.25e-5,2.8e-05,3.35e-05, 3.90e-05, 
#               4.45e-05, 5.00e-05]
# all_specs=[2.5e-6,7.5e-6,1.25e-5]
all_specs=[4.5e-6+n*1.5e-6 for n in range(10)]
# all_specs=[4.5e-6,6e-06,7.5e-6,9e-6,1.05e-5,1.2e-5,1.35e-5,1.5e-5,1.65e-5,1.8e-5]
# all_specs=[9e-6]
# all_specs=1.05e-5
#specs = [all_specs[iteration]]
specs=all_specs[0]

# specs=[3.5e-5]


# all_gds=[10**(-i) for i in range(1,8)]
# all_gds=[1e-4,2e-4,3e-4,4e-4,5e-4,6e-4,7e-4,8e-5,9e-4,1e-3]
all_gds=[10**(-5+.2*n) for n in range(11)]

sims = [1]
print_every=[1]
size=[80]
density=[13]
gds=[1e-4]
gds=all_gds

N=1
replicates = (i+1 for i in range(N))
params=itertools.product(specs, gds,replicates,sims,print_every,size,density)
#param_list=[i for i in params][0]

# # # start 4 tasks in parallel.
# # result_ids = []
# # for i in range(4):
# #     result_ids.append(neutral.remote(i))
# if __name__ == '__main__':
#     ray.init(num_cpus=5)
#     sp_obs=[neutral.remote(i) for i in params]
#     sp=ray.get(sp_obs[0])
#     ray.shutdown()

# #ray.cancel(sp)


#%timeit neutral(1e-6,1e-4,1,10,10,110,13)

#sp=neutral(1e-6,1e-4,1,10000,100,110,13)
#Counter(sp)
#Run through each combo of parameters N times in parallel
cores=4
if __name__ == '__main__':
    #mp.freeze_support()
    pool = mp.get_context('fork').Pool(cores)
    #pool = mp.Pool(cores)

    pool.starmap(neutral, params)
    pool.close()