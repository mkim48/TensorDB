#!/usr/bin/python
#
# comp_cov_sparse.py computes the covariance matrices for each mode of an input tensor in sparse representation for Dynamic Tensor Decomposition (DTA: Dynamic Tensor Analysis), which is used to incrementally update the tensor decomposition.
#
# Parameters
#  - tensor_name: input tensor
#  - tensor_size: input tensor size, e.g., 1000,1000,1000 for a tensor of size 1000x1000x1000
#  - chunk_size: chunk size for the input tensor, e.g., 100,100,100 for a chunk size 100x100x100
#  - old_tensor (optional): an old tensor on which the new tensor is incrementally updated.
#  - debug (optional): 1 for debugging mode (not running the operations but showing the commands), default value: 0
#
# Output
#  - The outputs include covariance matrices for each mode of the input tensor. The covariance matrices are <tensor_name>_cov0 for mode-1, <tensor_name>_cov1 for mode-2, etc.
#
# Examples
#  comp_cov_sparse.py tensor1 1000,1000,1000 100,100,100
#
#   - The example takes a 3-mode tensor of size 1000x1000x1000 with the chunk size 100x100x100 and computes the covariance matrices of the tensor. The output of this example includes covariance matrices (tensor1_cov0, tensor1_cov1, and tensor1_cov2).
#
#  comp_cov_sparse.py tensor2 1000,1000,1000 100,100,100 tensor1
#
#   - The example takes a 3-mode tensor of size 1000x1000x1000 with the chunk size 100x100x100 and incrementally update the covariance matrices of the old tensor tensor1 with the new tensor tensor2. The output of this example includes covariance matrices (tensor2_cov0, tensor2_cov1, and tensor2_cov2).
#
# BEGIN_COPYRIGHT
#
# This is the TensorDB, 2014.
# Reference:
#   Mijung Kim (2014). TensorDB and Tensor-based Relational Model (TRM) for Efficient Tensor-Relational Operations. Ph.D. Thesis. Arizona State University.
#
# This file is part of SciDB.
# Copyright (C) 2008-2011 SciDB, Inc.
#
# SciDB is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SciDB is distributed "AS-IS" AND WITHOUT ANY WARRANTY OF ANY KIND,
# INCLUDING ANY IMPLIED WARRANTY OF MERCHANTABILITY,
# NON-INFRINGEMENT, OR FITNESS FOR A PARTICULAR PURPOSE. See
# the GNU General Public License for the complete license terms.
#
# You should have received a copy of the GNU General Public License
# along with SciDB.  If not, see <http://www.gnu.org/licenses/>.
# END_COPYRIGHT
import csv
import os
import sys
import traceback
import datetime
sys.path.append(os.getcwd()) # NOCHECKIN 
sys.path.append('/opt/scidb/12.12' + '/lib')
import scidbapi as scidb


try:
    db = scidb.connect("localhost", 1239)
except Exception, inst:
    handleException(inst, True, op="connecting")
def handleException(inst, exitWhenDone, op=None):
    traceback.print_exc()
    if op:
        print >> sys.stderr, "Exception while ", op
    print >> sys.stderr, "     Exception Type: %s" % type(inst)     # the exception instance
    print >> sys.stderr, "     Exception Value: %r" % inst 
    print >> sys.stderr, ""
    if(exitWhenDone):
        exit(2)


def main():
	argcnt = len(sys.argv)
	if argcnt < 4:
		print "comp_cov_sparse.py <tensor_name> <tensor_size> <chunk_size> <old_tensor_name> <debug>"
		exit(2)

	start = datetime.datetime.now();

	tensz = sys.argv[2].split(',')
	tenchunksz = sys.argv[3].split(',')
	nmode=len(tensz)

	new_tensor=sys.argv[1]
	new_cov = new_tensor + '_cov'
	mat = new_tensor + '_mat'

	if argcnt >= 5:
		old_tensor=sys.argv[4]
		old_cov = old_tensor + '_cov'

	debug = 0
	if argcnt == 6:
		debug = 1

	print 'debug='+str(debug)

	if debug == 0:
		resfile = 'res_' + new_tensor
		f = open(resfile,'a')
		msg = 'comp_cov start'
		f.write(msg+'\n')

	for i in range(nmode):	
		prod_sz = 1
		prod_chunksz = 1
		for j in range(nmode):
			if i==j:
				continue
			prod_sz = prod_sz*int(tensz[j])
			prod_chunksz = prod_chunksz*int(tenchunksz[j])

		query="store(filter(matricize("+new_tensor+","+str(i+1)+"),matricize <> 0),"+mat+str(i)+")"
		if debug == 0:
			print query
			f.write(query+'\n')
			result=db.executeQuery(query,"afl")
			db.completeQuery(result.queryID)
		else:
			print query
		query="multiply_row("+mat+str(i)+","+mat+str(i)+")"

		query="store("+query+","+new_cov+str(i)+")"
		if debug == 0:
			print query
			result=db.executeQuery(query,"afl")
			db.completeQuery(result.queryID)
			f.write(query+'\n')
		else:
			print query

    	#if old_cov, update old_cov
		if argcnt >= 5:
			query="select a.multiply+b.multiply into "+new_cov+str(i)+" from "+new_cov+str(i)+" a, "+old_cov+str(i)+" b"
			if debug == 0:
				print query
				f.write(query+'\n')
				result=db.executeQuery(query,"aql")
				db.completeQuery(result.queryID)
			else:
				print query

	end = datetime.datetime.now();
	minutes, seconds = divmod((end-start).seconds, 60)
	microsec = (end-start).microseconds
	
	if debug == 0:
		msg = 'comp_cov Time elapsed: %d min %d.%d sec.' % (minutes, seconds, microsec)
		f.write(msg+'\n\n')
		f.close()

    #Disconnect from the SciDB server.
	db.disconnect() 
	sys.exit(0) #success

if __name__ == "__main__":
    main()
