# =============================================================================
# Bram L Gorissen, Massachusetts General Hospital, 2019.
# =============================================================================

import numpy as np
import scipy.io as sio
import scipy.sparse as sps
import h5py

def get_h5py_struct_array_value(f, structname, fieldname, idx):
    ref = f[structname][fieldname][idx,0]
    val = f[ref][()]
    return val

# =============================================================================
# This TROTS-reader can read the 20 proton cases only
# other cases use nonlinear objectives that are not supported by this reader
#
# The output are a matrix A and vectors b and c such that the problem
# becomes min{c'x : Ax <= b, x >= 0}
# =============================================================================

filename = 'Protons_01.mat'
f = h5py.File(filename)

A = {'data':[], 'row_ind':[], 'col_ind':[]}
b = []
c = None

# =============================================================================
# The objectives/constraints are read one-by-one and add elements to A/b/c.
# =============================================================================
row_offset = 0
num_aux_variables = 0
for i in range(len(f['problem']['dataID'])):
    data_id = int(get_h5py_struct_array_value(f, 'problem', 'dataID', i)[0][0])
    is_constraint = bool(get_h5py_struct_array_value(f, 'problem', 'IsConstraint', i)[0][0])
    minimize = bool(get_h5py_struct_array_value(f, 'problem', 'Minimise', i)[0][0])
    weight = float(get_h5py_struct_array_value(f, 'problem', 'Weight', i)[0][0])
    bound = float(get_h5py_struct_array_value(f, 'problem', 'Objective', i)[0][0])
    
    factor = 1 if minimize else -1
    
    # create the (i,j,value) pairs of the nonzeros in the constraint matrix
    ref = f['data']['matrix']['A'][data_id-1,0]
    if f[ref].attrs.get('MATLAB_class') == b'double':
        # dose matrix is in sparse format, so we can simply copy the
        # nonzero elements and their indices
        A['data'].append(factor*f[ref]['data'][()])
        A['row_ind'].append(row_offset + f[ref]['ir'][()])
        jc = f[ref]['jc'][()]
        # decompress compressed column indices
        col_ind = np.zeros(jc[-1], dtype=np.uint64)
        for j in range(len(jc)-1):
            for k in range(jc[j],jc[j+1]):
                col_ind[k] = j
        A['col_ind'].append(col_ind)
        
        num_voxels = int(f[ref].attrs.get('MATLAB_sparse'))
        num_pencil_beams = jc.size-1
        del jc, j, k
        
    else:
        # dose matrix is in dense format and is transposed
        DT = f[ref][()]
        (num_pencil_beams,num_voxels) = DT.shape
        (col_ind,row_ind) = np.nonzero(DT)
        A['data'].append(factor*DT[col_ind, row_ind])
        A['row_ind'].append(row_offset + row_ind)
        A['col_ind'].append(col_ind)
        del DT, col_ind, row_ind
    
    if not(is_constraint):
        # add an extra column to the constraint matrix for the objective variable
        A['data'].append(-factor*np.ones(num_voxels))
        A['row_ind'].append(row_offset + np.arange(num_voxels))
        A['col_ind'].append( (num_pencil_beams+num_aux_variables)*np.ones(num_voxels) )
        num_aux_variables += 1
    
    # create the right hand side
    if is_constraint:
        b.append(factor*bound*np.ones(num_voxels))
    else:
        b.append(np.zeros(num_voxels))
        
    # append the coefficient for the auxiliary variable to the objective vector
    if c is None:
        c = np.zeros(num_pencil_beams)
    if not(is_constraint):
        c = np.append(c, factor*weight)
        
    if False and is_constraint:
        # this OPTIONAL section changes each mean/min/max-dose constraint
        # into a mean under/-overdose constraint
        # add an identity matrix to the coefficient matrix for the the auxiliary variables
        A['data'].append(-np.ones(num_voxels))
        A['row_ind'].append(row_offset + np.arange(num_voxels))
        A['col_ind'].append(num_pencil_beams + num_aux_variables + np.arange(num_voxels))
        
        # add the row sum y_i / n <= 0.1, with y_i the auxiliary variables
        A['data'].append(np.ones(num_voxels) / num_voxels)
        A['row_ind'].append( (row_offset + num_voxels) * np.ones(num_voxels) )
        A['col_ind'].append(num_pencil_beams + num_aux_variables + np.arange(num_voxels))
        b.append([0.1])
        c = np.append(c, np.zeros(num_voxels))
        num_aux_variables += num_voxels
        row_offset += 1
        
       
    row_offset += num_voxels

del f, i, data_id, is_constraint, minimize, weight, bound, factor, num_voxels
     
# set up the primal problem min{c'x : Ax <= b, x >= 0}
num_variables = num_pencil_beams + num_aux_variables
A['data']     = np.concatenate(A['data'])
A['row_ind']  = np.concatenate(A['row_ind'])
A['col_ind']  = np.concatenate(A['col_ind'])
A = sps.csr_matrix((A['data'], (A['row_ind'], A['col_ind'])), shape=(row_offset, num_variables))
b = np.concatenate(b)

# Here, use A, b and c
#solve(A, b, c)
