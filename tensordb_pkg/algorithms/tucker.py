#!/usr/bin/python
#
# tucker.py runs Tucker Decomposition. Given the input tensor, perform the tucker decomposition. The Tucker decomposition takes HOSVD of the input tensor as the initial factor matrices and core. HOSVD of the input tensor is obtained by running dta_full.py without an old tensor.
#
# Parameters
#  - tensor_name: input tensor
#  - tensor_size: input tensor size, e.g., 1000,1000,1000 for a tensor of size 1000x1000x1000
#  - chunk_size: chunk size for the input tensor, e.g., 100,100,100 for a chunk size 100x100x100
#  - rank: target ranks for each mode, e.g., 2,3,4
#  - max_iter (optional): the maximum iteration count, default value: 50
#  - debug (optional): 1 for debugging mode (not running the operations but showing the commands), default value: 0
#
# Output 
#  - The outputs include factor matrices and a core. The factor matrices are <tensor_name>_fac_0, <tensor_name>_fac_1, ..., and the core is <tensor_name>_core.
#
# Example 
#  tucker.py tensor1 1000,1000,1000 100,100,100 2,3,4
#
#   - The example takes a 3-mode tensor of size 1000x1000x1000 with the chunk size 100x100x100 and decompose the tensor. The output of the example includes 3 factor matrices (tensor1_fac_0, tensor1_fac_1, and tensor1_fac_2) and a core (tensor1_core).
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
		print "tucker.py <tensor_name> <tensor_size> <chunk_size> <rank> <max_iter> <debug>"
		exit(2)

	start = datetime.datetime.now();

	tensor=sys.argv[1]
	facmat = tensor + '_fac_'
	core = tensor + '_core'
	ztensor = tensor + '_z'
	tmp = tensor + '_tmpp'
	ttm = tensor + '_ttm'
	
	if argcnt >= 6:
		max_iter=int(sys.argv[5])
	else:
		max_iter=50

	debug = 0
	if argcnt == 7:
		debug = 1

	print 'debug='+str(debug)

	if debug == 0:
		resfile = 'res_' + tensor
		f = open(resfile,'a')
		msg = 'tucker start'
		f.write(msg+'\n')

	facsz=sys.argv[2].split(',')
	facchunksz=sys.argv[3].split(',')
	rank=sys.argv[4].split(',')
	nmode=len(facsz)

	fit = 0
	fitchangetol = 1e-4

	for r in range(max_iter):
		fitold = fit
		#compute z=ttm(tensor,facmat,-mode)
		for m in range(nmode):
			if m == 0:
				prod_sz = 1
				prod_chunksz = 1
				for i in range(0,nmode):
					if i == 1:
						continue
					prod_sz = prod_sz * int(facsz[i])
					prod_chunksz = prod_chunksz * int(facchunksz[i])
				query = "store(reshape("+tensor+",<val:double>[i=1:"+str(facsz[1])+","+str(facchunksz[1])+",0,j=1:"+str(prod_sz)+","+str(prod_chunksz)+",0]),"+ttm+str(m)+")"
				if debug == 0:
					print query
					result=db.executeQuery(query,"afl")
					db.completeQuery(result.queryID)
					f.write(query+'\n');
				else:
					print query
			else:
				prod_sz = 1
				prod_chunksz = 1
				for i in range(1,nmode):
					prod_sz = prod_sz * int(facsz[i])
					prod_chunksz = prod_chunksz * int(facchunksz[i])
				query = "store(reshape("+tensor+",<val:double>[i=1:"+str(facsz[0])+","+str(facchunksz[0])+",0,j=1:"+str(prod_sz)+","+str(prod_chunksz)+",0]),"+ttm+str(m)+")"
				print '#1'
				print prod_sz
				if debug == 0:
					print query
					result=db.executeQuery(query,"afl")
					db.completeQuery(result.queryID)
					f.write(query+'\n');
				else:
					print query
			query = ttm+str(m)

			for i in range(nmode):	
				if m == i:
					continue

				if m != nmode and i == nmode:
					query = "multiply_row(transpose("+query+"),"+facmat+str(i)+")"	
				else:
					#next = i+1
					#if m == i+1:
					#		next = i+2
					#if next != nmode:
					prod_sz = prod_sz / int(facsz[i])
					prod_chunksz = prod_chunksz / int(facchunksz[i])
					str_idx1 = "i=1:"+str(facsz[i])+","+str(facchunksz[i])+",0"
						#prod_sz = 1
						#prod_chunksz = 1
						#for j in range(next+1,nmode):
						#for j in range(0,next):
					prod_sz = prod_sz * int(rank[i])
					prod_chunksz = prod_chunksz * int(rank[i])
					print prod_sz
					str_idx2 = "j=1:"+str(prod_sz)+","+str(prod_chunksz)+",0"
					#else:
					#	str_idx1 = "i=1:"+str(facsz[i])+","+str(facchunksz[i])+",0"
						
					query = "reshape(multiply_row(transpose("+query+"),"+facmat+str(i)+"),<val:double>["+str_idx1+","+str_idx2+"])"	
				query="store("+query+","+ttm+str(m)+str(i)+")";
				if debug == 0:
					print query
					result=db.executeQuery(query,"afl")
					db.completeQuery(result.queryID)
					f.write(query+'\n');
				else:
					print query
				query = ttm+str(m)+str(i)

			# compute ztensor 
			query="store(multiply_row("+query+","+query+"),"+ztensor+str(m)+")";
			if debug == 0:
				print query
				result=db.executeQuery(query,"afl")
				db.completeQuery(result.queryID)
				f.write(query+'\n');
			else:
				print query

		#remove factor matrices (the chunk size is same as the tensor size)
		for i in range(nmode):
			query="remove("+facmat+str(i)+")"
			if debug == 0:
				print query
				f.write(query+'\n')
				try:
					result=db.executeQuery(query,"afl")
					db.completeQuery(result.queryID)
				except Exception, inst:
					print >> sys.stderr, "     Exception Type: %s" % type(inst)     # the exception instance
					print >> sys.stderr, "     Exception Value: %r" % inst
			else:
				print query
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
					print >> sys.stderr, "     Exception Type: %s" % type(inst)     # the exception instance
					print >> sys.stderr, "     Exception Value: %r" % inst
			else:
				print query

		#compute factor matrices
		for i in range(nmode):	
			if facsz[i]==facchunksz[i]:
				query="copyArray(eigen("+ztensor+str(i)+","+rank[i]+"),"+facmat+str(i)+")"
			else:
				query = "repart("+ztensor+str(i)+",<val:double>[i=1:"+str(facsz[i])+","+str(facsz[i])+",0,j=1:"+str(facsz[i])+","+str(facsz[i])+",0])"
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
			print i
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
			query = "reshape(multiply_row(transpose("+query+"),"+facmat+str(i-1)+"),<val:double>["+str_idx1+","+str_idx2+"])"	

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

	end = datetime.datetime.now();
	minutes, seconds = divmod((end-start).seconds, 60)
	microsec = (end-start).microseconds
	
	if debug == 0:
		msg = 'tucker Time elapsed: %d min %d.%d sec.' % (minutes, seconds, microsec)
		f.write(msg+'\n\n')
		f.close()

    #Disconnect from the SciDB server.
	db.disconnect() 
	sys.exit(0) #success

if __name__ == "__main__":
    main()
