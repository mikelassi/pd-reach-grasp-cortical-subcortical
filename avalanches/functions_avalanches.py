import numpy as np

# avalanches_global_pattern(x, n)
# INPUT: binary matrix x representing active region, structured as regions x %times; n=number_regions
# f is a matrix (containing the consecutive series of activation of one patients (by consecutive I mean that all the clean bits have been attached one after the other)
# of axb, where a is the number of areas and b is the time. The argument n will specify how many regions to take into account. h will provide y vectors, with y the number
# of avalanches, each vector being long n, that is full with logicals where 0 means that that particular regions has never been active in that
# particular avalanches, and 1 means the opposite. If a region has been active only in 1 time bin or in more than 1 time bin, it will be 1 in the vector eitherway.

def avalanches_global_pattern(x, n):
    x = x[:n, :]
    activations_bins = np.where(np.any(x, axis=0))[0]
    activations_bins = np.concatenate((activations_bins, [5, 6]))  # Adding dummy values at the end
    gg = np.diff(activations_bins) != 1
    Index_start = np.where(gg == 1)[0]

    f = []
    h = []
    f.append(x[:, activations_bins[0]:activations_bins[Index_start[0]]+1])

    for ii in range(len(Index_start)-1):
        f.append(x[:, activations_bins[Index_start[ii]+1]:activations_bins[Index_start[ii + 1]]+1]) 

    for kk in range(len(f)):
        a = np.sum(f[kk], axis=1)
        #a[a > 0] = 1 # here to have any activation = 1 (to compute flexibility)
        h.append(a)
    
    flex = len(np.unique(h,axis=0))  
    return f, h, flex

def func_transition_matrix(b):
    """
    Compute the transition matrix from a binarized spatiotemporal activity pattern.

    For each brain region, this function calculates how often it transitions to 
    other regions in the next time step, based on observed activations. The result 
    is a matrix where each entry (i, j) represents the normalized frequency with 
    which region i is followed by region j.

    Parameters
    ----------
    binary_sequence : np.ndarray
        2D binary matrix of shape (n_regions, n_timepoints), where each element is 0 or 1,
        indicating whether a region is inactive or active at a given time.

    Returns
    -------
    transition_matrix : np.ndarray
        2D matrix of shape (n_regions, n_regions), where element (i, j) represents the 
        normalized probability that region i at time t is followed by region j at time t+1.
    """
    num_regions, num_times = b.shape
    out = np.zeros((num_regions, num_regions))
    
    num_activations_per_reg = np.sum(b, axis=1)
    
    for kk1 in range(num_regions):
        g_loop = np.zeros(num_regions)
        for kk2 in range(num_times - 1):
            if b[kk1, kk2] == 1:
                curr_reg = np.repeat(b[kk1, kk2], num_regions)
                next_reg = b[:, kk2 + 1]
                g = curr_reg == next_reg
                g_loop += g
        g_loop[g_loop != 0] /= num_activations_per_reg[kk1]
        out[kk1, :] = g_loop
    
    return out

def compute_ATM(binarized_data):
    """
    Compute the Average Transition Matrix (ATM) and avalanche statistics for a set of binarized brain activity time series.

    Parameters
    ----------
    binarized_data : np.ndarray
        3D binary matrix of shape (n_regions, n_timepoints, n_subjects), representing brain activity
        over time for multiple subjects. Each entry is 0 or 1, indicating inactive or active states.

    Returns
    -------
    mean_avalanche_durations : np.ndarray
        1D array of shape (n_subjects,) containing the mean duration (in timepoints) of avalanches
        for each subject.
    
    subject_ATMs : np.ndarray
        3D array of shape (n_subjects, n_regions, n_regions) containing the Average Transition Matrix
        (ATM) for each subject, computed by averaging the symmetric transition matrices of all
        individual avalanches.

    avalanche_counts : list of int
        Number of avalanches detected for each subject.
    """

    n_subjects = binarized_data.shape[2]
    n_regions = binarized_data.shape[0]

    mean_avalanche_durations = np.zeros(n_subjects)    # Average avalanche duration per subject
    subject_ATMs = []                                   # List of ATMs per subject
    avalanche_counts = []                               # Number of avalanches per subject
    flexibility = []                                    # Flexibility per subject
    
    # Iterate over all subjects
    for subject_idx in range(n_subjects):
        # Detect avalanches for current subject
        avalanche_list, _, flex = avalanches_global_pattern(binarized_data[:, :, subject_idx], n_regions)

        n_avalanches = len(avalanche_list)
        print(f"Subject {subject_idx}: Detected {n_avalanches} avalanches.")
        avalanche_lengths = np.zeros(n_avalanches)  # Store the length of each avalanche
        transition_matrices = np.zeros((n_regions, n_regions, n_avalanches))  # Transition matrices per avalanche

        # Process each avalanche
        for aval_idx, avalanche in enumerate(avalanche_list):
            avalanche_lengths[aval_idx] = avalanche.shape[1]  # Number of timepoints in this avalanche

            transition_matrix = func_transition_matrix(avalanche)  # Compute transition matrix
            symmetric_matrix = (transition_matrix + transition_matrix.T) / 2  # Symmetrize
            transition_matrices[:, :, aval_idx] = symmetric_matrix

        # Store mean duration of avalanches for this subject
        mean_avalanche_durations[subject_idx] = np.mean(avalanche_lengths)
        avalanche_counts.append(n_avalanches)
        flexibility.append(flex)

        atm = np.sum(transition_matrices, axis=2) / np.sum(transition_matrices != 0, axis=2)
        atm[np.isnan(atm)] = 0
        subject_ATMs.append(atm)       
    return mean_avalanche_durations, np.array(subject_ATMs), avalanche_counts,flexibility

def compute_ATM_no_simm(binarized_data):
    """
    Compute the Average Transition Matrix (ATM) and avalanche statistics for a set of binarized brain activity time series.

    Parameters
    ----------
    binarized_data : np.ndarray
        3D binary matrix of shape (n_regions, n_timepoints, n_subjects), representing brain activity
        over time for multiple subjects. Each entry is 0 or 1, indicating inactive or active states.

    Returns
    -------
    mean_avalanche_durations : np.ndarray
        1D array of shape (n_subjects,) containing the mean duration (in timepoints) of avalanches
        for each subject.
    
    subject_ATMs : np.ndarray
        3D array of shape (n_subjects, n_regions, n_regions) containing the Average Transition Matrix
        (ATM) for each subject, computed by averaging the symmetric transition matrices of all
        individual avalanches.

    avalanche_counts : list of int
        Number of avalanches detected for each subject.
    """

    n_subjects = binarized_data.shape[2]
    n_regions = binarized_data.shape[0]

    mean_avalanche_durations = np.zeros(n_subjects)    # Average avalanche duration per subject
    subject_ATMs = []                                   # List of ATMs per subject
    avalanche_counts = []                               # Number of avalanches per subject

    # Iterate over all subjects
    for subject_idx in range(n_subjects):
        # Detect avalanches for current subject
        avalanche_list, _ ,flex = avalanches_global_pattern(binarized_data[:, :, subject_idx], n_regions)

        n_avalanches = len(avalanche_list)
        avalanche_lengths = np.zeros(n_avalanches)  # Store the length of each avalanche
        transition_matrices = np.zeros((n_regions, n_regions, n_avalanches))  # Transition matrices per avalanche

        # Process each avalanche
        for aval_idx, avalanche in enumerate(avalanche_list):
            avalanche_lengths[aval_idx] = avalanche.shape[1]  # Number of timepoints in this avalanche

            transition_matrix = func_transition_matrix(avalanche)  # Compute transition matrix
            transition_matrices[:, :, aval_idx] = transition_matrix

        # Store mean duration of avalanches for this subject
        mean_avalanche_durations[subject_idx] = np.mean(avalanche_lengths)
        avalanche_counts.append(n_avalanches)
        

        atm = np.sum(transition_matrices, axis=2) / np.sum(transition_matrices != 0, axis=2)
        atm[np.isnan(atm)] = 0
        subject_ATMs.append(atm)       
    return mean_avalanche_durations, np.array(subject_ATMs), avalanche_counts

def avalanches_detection(x_binary, x_signal, n=None):
    """
    Detects neural avalanches from binary and real-valued signal matrices.

    Parameters:
    ----------
    x_binary : ndarray (regions x time)
        Binary matrix indicating activation (e.g., x_binary[i, t] = 1 if region i is active at time t).
    x_signal : ndarray (regions x time)
        Original real-valued signal matrix corresponding to x_binary.
    n : int, optional
        Number of regions to consider (if None, all regions are used).

    Returns:
    -------
    avalanches_signal : list of ndarrays
        List where each element is the original signal (regions x duration) during a detected avalanche.
    activation_patterns : list of 1D ndarrays
        List of activation vectors per avalanche (each vector has length = number of regions,
        with values representing how many time bins each region was active during that avalanche).
    avalanche_bounds : list of tuples
        List of (start_index, end_index) tuples representing the time range of each avalanche.
    """
    activations_bins = np.where(np.any(x_binary, axis=0))[0]
    diff = np.diff(activations_bins)
    split_points = np.where(diff != 1)[0]
    start_indices = [activations_bins[0]] + [activations_bins[i+1] for i in split_points]
    end_indices = [activations_bins[i] for i in split_points] + [activations_bins[-1]]

    avalanches_signal = []
    activation_patterns = []
    avalanche_bounds = []

    for start, end in zip(start_indices, end_indices):
        segment_signal = x_signal[:, start:end+1]
        segment_binary = x_binary[:, start:end+1]

        pattern = np.sum(segment_binary, axis=1)
        pattern=np.where(pattern > 0, 1, 0)
        avalanches_signal.append(segment_signal)
        activation_patterns.append(pattern)
        avalanche_bounds.append((start, end))

    return avalanches_signal, activation_patterns, avalanche_bounds

def avalanches_detection_with_dummy(x_binary, x_signal, n=None): 
    if n is not None:
        x_binary = x_binary[:n, :]
        x_signal = x_signal[:n, :]

    activations_bins = np.where(np.any(x_binary, axis=0))[0]
    dummy_tail = np.array([x_binary.shape[1], x_binary.shape[1] + 1])
    activations_bins = np.concatenate((activations_bins, dummy_tail))

    diff = np.diff(activations_bins)
    split_points = np.where(diff != 1)[0]
    start_indices = [activations_bins[0]] + [activations_bins[i+1] for i in split_points]
    end_indices = [activations_bins[i] for i in split_points]

    avalanches_signal = []
    activation_patterns = []
    avalanche_bounds = []

    for start, end in zip(start_indices, end_indices):
        if end >= x_binary.shape[1]:  # skip dummy points
            continue
        segment_signal = x_signal[:, start:end+1]
        segment_binary = x_binary[:, start:end+1]

        pattern = np.sum(segment_binary, axis=1)
        pattern = np.where(pattern > 0, 1, 0)
        avalanches_signal.append(segment_signal)
        activation_patterns.append(pattern)
        avalanche_bounds.append((start, end))

    return avalanches_signal, activation_patterns, avalanche_bounds
