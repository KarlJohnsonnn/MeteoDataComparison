

def findBasesTops(dbz_m, range_v):
    """
    % FINDBASESSTOPS
    %
    % functionality:
    % find cloud bases and tops from radar reflectivity profiles for up to 10 cloud layers
    % no cloud = NaN
    %
    % input:
    %   dbz       ... reflectivity matrix         [dbz] (range x time)
    %   range     ... radar height vector         [m or km] (range)
    %
    % output:
    %   bases     ... matrix of indices (idx) of cloud bases  (10 x time)
    %   tops      ... matrix of indices (idx) of cloud tops   (10 x time)
    %   base_m    ... matrix of heights of cloud bases  [m or km, same unit as range] (10 x time), 1st base = -1 means no cloud detected
    %   top_m     ... matrix of heights of cloud tops   [m or km, same unit as range] (10 x time), 1st top  = -1 means no cloud detected
    %   thickness ... matrix of cloud thickness         [m or km, same unit as range] (10 x time)
    """


    shape_dbz = dbz_m.shape
    len_time  = shape_dbz[1]
    len_range = len(range_v)

    bases = np_NaN(10,len_time)
    tops = np_NaN(10,len_time)
    thickness = np_NaN(10,len_time)

    top_m  = np_NaN(10,len_time) # max. 10 cloud layers detected
    base_m = np_NaN(10,len_time) # max. 10 cloud layers detected

    if pts:
        print('')
        print(' dBZ(i,:) = ', [dbz_m[i, :] for i in range(shape_dbz[0])])

    for i in range(0,len_time):

        in_cloud  = 0
        layer_idx = 0
        current_base = np.nan

        if pts: print("    Searching for cloud bottom and top ({} of {}) time steps".format(i+1,len_time), end="\r")

        # found the first base in first bin.
        if ( not np.isnan(dbz_m[0,i]) ):
            layer_idx = 1
            current_base = 1
            in_cloud = 1

        for j in range(1,len_range):

            if (in_cloud == 1): # if in cloud

                # cloud top found at (j-1)
                if np.isnan(dbz_m[j,i]):

                    current_top = j-1
                    thickness[layer_idx,i] = range_v[current_top] - range_v[current_base]
                    bases[layer_idx,i]     = current_base  # bases is an idx
                    tops[layer_idx,i]      = current_top    # tops is an idx

                    base_m[layer_idx,i]    = range_v[current_base]  # cloud base in m or km
                    top_m[layer_idx,i]     = range_v[current_top]   # cloud top in m or km

                    print(str(i)+': found '+str(layer_idx)+'. cloud ['+str(bases[layer_idx,i])+', '+\
                        str(tops[layer_idx,i])+'], thickness: '+str(thickness)+'km')

                    in_cloud = 0

            else: # if not in cloud

                # cloud_base found at j
                if ( not np.isnan(dbz_m[j,i]) ):
                    layer_idx += 1
                    current_base = j
                    in_cloud = 1

        # at top height but still in cloud, force top
        if ( in_cloud == 1 ):
            tops[layer_idx,i]  = len(range_v)
            top_m[layer_idx,i] = max(range_v)    #  cloud top in m or km



    ###
    # keep only first 10 cloud layers
    bases     = bases[:10,:]
    tops      = tops [:10,:]
    base_m    = base_m[:10,:]
    top_m     = top_m [:10,:]
    thickness = thickness[:10,:]
    # give clear sky flag when first base_m ==NaN (no dbz detected over all heights),
    # problem: what if radar wasn't working, then dbz would still be NaN!
    loc_nan = np.where(np.isnan(base_m[0,:]))

    if pts:
        print('vor bases = ',[ba for ba in base_m])

    base_m[0,np.where(np.isnan(base_m[0,:]))] = -1
    top_m[0,np.where(np.isnan(top_m[0,:]))]   = -1

    if pts:
        print('')
        print('npisnan = ',[loc_nan[i] for i in range(len(loc_nan))])


        if pts: print('')
        for ba in base_m:
            print(' bases = ',ba)

        for to in top_m:
            print(' tops = ', to)


    return bases, tops, base_m, top_m, thickness