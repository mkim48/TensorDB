#!/usr/bin/python
#
<<<<<<< HEAD
# dta.py runs Dynamic Tensor Decomposition (DTA: Dynamic Tensor Analysis). Given the input tensor, incrementally update the tensor decomposition. This includes operations for computing the covariance matrices for an input tensor in dense representation.  
=======
# dta.py runs Dynamic Tensor Decomposition (DTA: Dynamic Tensor Analysis). Given the input tensor, incrementally update the tensor decomposition. This operation requires the covariance matrices of the input tensor, which are computed by comp_cov_dense.py, comp_cov_sparse.py, or comp_cov_comp.py. 
>>>>>>> 179692c168a8120cac788efb45c7d537185433e6
#
# Parameters
#  - tensor_name: input tensor
#  - tensor_size: input tensor size, e.g., 1000,1000,1000 for a tensor of size 1000x1000x1000
#  - chunk_size: chunk size for the input tensor, e.g., 100,100,100 for a chunk size 100x100x100
#  - rank: target ranks for each mode, e.g., 5,10,15
<<<<<<< HEAD
#  - old tensor (optional): an old tensor on which the new tensor is incrementally updated.
=======
>>>>>>> 179692c168a8120cac788efb45c7d537185433e6
#  - debug (optional): 1 for debugging mode (not running the operations but showing the commands), default value: 0
#
# Output
#  - The outputs include factor matrices and a core. The factor matrices are <tensor_name>_fac_0, <tensor_name>_fac_1, ..., and the core is <tensor_name>_core.
#
# Examples
#  dta.py tensor1 1000,1000,1000 100,100,100 2,3,4
#
<<<<<<< HEAD
#   - The example takes a 3-mode tensor of size 1000x1000x1000 with the chunk size 100x100x100 and decompose the tensor. The output of this example includes factor matrices (tensor2_fac_0, tensor2_fac_1, and tensor2_fac_2) and a core (tensor2_core). 
#
#  dta.py tensor2 1000,1000,1000 100,100,100 2,3,4 tensor1
#
#   - The example takes a 3-mode tensor of size 1000x1000x1000 with the chunk size 100x100x100 and incrementally update the tensor decomposition of the old tensor tensor1 with the new tensor tensor2. The output of this example includes factor matrices (tensor2_fac_0, tensor2_fac_1, and tensor2_fac_2) and a core (tensor2_core).
#
# Note: the initial DTA (without the old tensor to be incrementally updated) is used as the input of the Tucker decomposition.
=======
#   - The example takes a 3-mode tensor of size 1000x1000x1000 with the chunk size 100x100x100 and decompose the tensor. The output of this example includes fac
tor matrices (tensor1_fac_0, tensor1_fac_1, and tensor1_fac_2) and a core (tensor1_core).
>>>>>>> 179692c168a8120cac788efb45c7d537185433e6
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
#
import csv
import os
import sys
import traceback
import datetime
import math
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
	if argcnt < 5:
<<<<<<< HEAD
		print "dta.py <tensor_name> <tensor_size> <chunk_size> <rank> <old_tensor> <debug>"
=======
		print "dta.py <tensor_name> <tensor_size> <chunk_size> <rank> <debug>"
>>>>>>> 179692c168a8120cac788efb45c7d537185433e6
		exit(2)

	start = datetime.datetime.now();

	tensor=sys.argv[1]
	facmat = tensor + '_fac_'
	core = tensor + '_core'
	cov = tensor + '_cov'
	tmp = tensor + '_tmp'

<<<<<<< HEAD
	if argcnt >= 6:
		old_tensor=sys.argv[5]
		old_cov = old_tensor + '_cov'

	debug = 0
	if argcnt == 7:
		debug = 1

	#print 'debug='+str(debug)
=======
	debug = 0
	if argcnt == 6:
		debug = 1

	print 'debug='+str(debug)
>>>>>>> 179692c168a8120cac788efb45c7d537185433e6

	if debug == 0:
		resfile = 'res_' + tensor
		f = open(resfile,'a')
<<<<<<< HEAD
		msg = 'dta start'
=======
		msg = 'dtd start'
>>>>>>> 179692c168a8120cac788efb45c7d537185433e6
		f.write(msg+'\n')

	facsz=sys.argv[2].split(',')
	facchunksz=sys.argv[3].split(',')
	rank=sys.argv[4].split(',')
	nmode=len(facsz)

<<<<<<< HEAD
	new_tensor=sys.argv[1]
	new_cov = new_tensor + '_cov'
	mat = new_tensor + '_mat'

	for i in range(nmode):
		prod_sz = 1
		prod_chunksz = 1
		for j in range(nmode):
			if i==j:
				continue
			prod_sz = prod_sz*int(facsz[j])
			prod_chunksz = prod_chunksz*int(facchunksz[j])

		query="store(matricize("+new_tensor+","+str(i+1)+"),"+mat+str(i)+")"
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
		if argcnt >= 6:
			query="select a.multiply+b.multiply into "+new_cov+str(i)+" from "+new_cov+str(i)+" a, "+old_cov+str(i)+" b"
			if debug == 0:
				print query
				f.write(query+'\n')
				result=db.executeQuery(query,"aql")
				db.completeQuery(result.queryID)
			else:
				print query

=======
>>>>>>> 179692c168a8120cac788efb45c7d537185433e6
	#create factor matrices (the chunk size is same as the tensor size)
	for i in range(nmode):
		query="create array "+facmat+str(i)+"<val:double>[i=1:"+str(rank[i])+","+str(rank[i])+",0,j=1:"+str(facsz[i])+","+str(facchunksz[i])+",0]"
		if debug == 0:
			print query
			f.write(query+'\n')
			try:
				result=db.executeQuery(query,"aql")
				db.completeQuery(result.queryID)
			except Exception, inst:
<<<<<<< HEAD
				if debug == 1:
					print >> sys.stderr, "     Exception Type: %s" % type(inst)     # the exception instance
					print >> sys.stderr, "     Exception Value: %r" % inst
=======
				print >> sys.stderr, "     Exception Type: %s" % type(inst)     # the exception instance
				print >> sys.stderr, "     Exception Value: %r" % inst
>>>>>>> 179692c168a8120cac788efb45c7d537185433e6
		else:
			print query

	#create core
	str_idx = ""
	for i in range(nmode):
		str_idx = str_idx + "i"+str(i)+"=1:"+str(rank[i])+","+str(rank[i])+",0"
		if i < len(facsz)-1:
			str_idx = str_idx + ","	

	query="create array "+core+" <val:double>["+str_idx+"]"
	if debug == 0:
		print query
		f.write(query+'\n')
		try:
			result=db.executeQuery(query,"aql")
			db.completeQuery(result.queryID)
		except Exception, inst:
<<<<<<< HEAD
			if debug == 1:
				print >> sys.stderr, "     Exception Type: %s" % type(inst)     # the exception instance
				print >> sys.stderr, "     Exception Value: %r" % inst

=======
			print >> sys.stderr, "     Exception Type: %s" % type(inst)     # the exception instance
			print >> sys.stderr, "     Exception Value: %r" % inst
>>>>>>> 179692c168a8120cac788efb45c7d537185433e6
	else:
		print query
	
	#compute factor matrices
	for i in range(nmode):	
		if facsz[i]==facchunksz[i]:
			query="copyArray(eigen("+cov+str(i)+","+rank[i]+"),"+facmat+str(i)+")"
		else:
			query = "repart("+cov+str(i)+",<val:double>[i=1:"+str(facsz[i])+","+str(facsz[i])+",0,j=1:"+str(facsz[i])+","+str(facsz[i])+",0])"
			query="repart(eigen("+query+","+rank[i]+"),<val:double>[i=1:"+str(facsz[i])+","+str(facsz[i])+",0,j=1:"+str(facsz[i])+","+str(facchunksz[i])+",0])"
			query="copyArray("+query+","+facmat+str(i)+")"
		if debug == 0:
			print query
			result=db.executeQuery(query,"afl")
			db.completeQuery(result.queryID)
			f.write(query+'\n')
		else:
			print query
				
	#compute core=ttm(tensor,facmat)
	prod_sz = 1
	prod_chunksz = 1
	for i in range(1,nmode):
		prod_sz = prod_sz * int(facsz[i])
		prod_chunksz = prod_chunksz * int(facchunksz[i])
	query = "store(reshape("+tensor+",<val:double>[i=1:"+str(facsz[0])+","+str(facchunksz[0])+",0,j=1:"+str(prod_sz)+","+str(prod_chunksz)+",0]),"+tmp+str(0)+")"
	if debug == 0:
		print query
		result=db.executeQuery(query,"afl")
		db.completeQuery(result.queryID)
		f.write(query+'\n');
	else:
		print query
	query = tmp+str(0)

	for i in range(1,nmode):	
<<<<<<< HEAD
=======
		print i
>>>>>>> 179692c168a8120cac788efb45c7d537185433e6
		str_idx1 = "i=1:"+str(facsz[i])+","+str(facchunksz[i])+",0"
		prod_sz = 1
		prod_chunksz = 1
		for j in range(i+1,nmode):
			prod_sz = prod_sz * int(facsz[j])
			prod_chunksz = prod_chunksz * int(facchunksz[j])
		for j in range(0,i):
			prod_sz = prod_sz * int(rank[j])
			prod_chunksz = prod_chunksz * int(rank[j])
		str_idx2 = "j=1:"+str(prod_sz)+","+str(prod_chunksz)+",0"
		query = "reshape(multiply_row(transpose("+query+")," + facmat+str(i-1)+"),<val:double>["+str_idx1+","+str_idx2+"])"	
		query="store("+query+","+tmp+str(i)+")";
		if debug == 0:
			print query
			result=db.executeQuery(query,"afl")
			db.completeQuery(result.queryID)
			f.write(query+'\n');
		else:
			print query
		query = tmp+str(i)

	#core
	str_idx = "i0=1:"+str(rank[0])+","+str(rank[0])+",0"
	for j in range(1,nmode):
		str_idx = str_idx + ",i"+str(j)+"=1:"+str(rank[j])+","+str(rank[j])+",0"
	query = "reshape(multiply_row(transpose("+query+"),"+facmat+str(nmode-1)+"),<val:double>["+str_idx+"])"
	query="store("+query+","+core+")";
	if debug == 0:
		print query
		result=db.executeQuery(query,"afl")
		db.completeQuery(result.queryID)
		f.write(query+'\n');
	else:
		print query

	#compute norm^2 for input tensor
	query="select sum(pow(val,2)) from " + tensor
	if debug == 0:
		print query
		f.write(query+'\n')
		result=db.executeQuery(query,"aql")
		dataitem = result.array.getConstIterator(0).getChunk().getConstIterator().getItem()
		attrs = result.array.getArrayDesc().getAttributes()
		normx =  scidb.getTypedValue(dataitem, attrs[0].getType())
		db.completeQuery(result.queryID)
	else:
		print query

	#compute norm^2 for core 
	query="select sum(pow(val,2)) from " + core
	if debug == 0:
		print query
		f.write(query+'\n')
		result=db.executeQuery(query,"aql")
		dataitem = result.array.getConstIterator(0).getChunk().getConstIterator().getItem()
		attrs = result.array.getArrayDesc().getAttributes()
		normcore =  scidb.getTypedValue(dataitem, attrs[0].getType())
		db.completeQuery(result.queryID)
	else:
		print query

	# frobenius norm of x - x^
	if debug == 0: 
		f.write("normx="+str(normx)+"\n")
		f.write("normcore="+str(normcore)+"\n")
		norm_residual = math.sqrt(normx - normcore)
		f.write("norm_residual="+str(norm_residual)+"\n")
		fit = 1-norm_residual/math.sqrt(normx)
		f.write("fit="+str(fit)+'\n')
<<<<<<< HEAD
		print "\nfit="+str(fit)+'\n'
=======
>>>>>>> 179692c168a8120cac788efb45c7d537185433e6

	end = datetime.datetime.now();
	minutes, seconds = divmod((end-start).seconds, 60)
	microsec = (end-start).microseconds
	
	if debug == 0:
<<<<<<< HEAD
		msg = 'dta Time elapsed: %d min %d.%d sec.' % (minutes, seconds, microsec)
		f.write(msg+'\n\n')
		print msg+'\n\n'
=======
		msg = 'dtd Time elapsed: %d min %d.%d sec.' % (minutes, seconds, microsec)
		f.write(msg+'\n\n')
>>>>>>> 179692c168a8120cac788efb45c7d537185433e6
		f.close()

    #Disconnect from the SciDB server.
	db.disconnect() 
	sys.exit(0) #success

if __name__ == "__main__":
    main()
