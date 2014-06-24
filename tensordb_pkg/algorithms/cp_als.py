#!/usr/bin/python
#
# cp_als.py runs the CP ALS (Alternating Least Squares) algorithm
#
# Parameters
#  - tensor name: input tensor
#  - tensor size: input tensor size, e.g., 1000,1000,1000 for a tensor of size 1000x1000x1000
#  - chunk size: chunk size for the input tensor, e.g., 100,100,100 for a chunk size 100x100x100
#  - target rank: target rank, e.g.,: 10
#  - modes of materialized matricization (optional): modes which use materialized matricization, e.g., 1,2 - mode-1 and 2 use the materialized matricization
#  - max_iter (optional): the maximum iteration count, default value: 50
#  - debug (optional): 1 for debugging mode (not running the operations but showing the commands), default value: 0
#
# Output 
#  - The outputs include factor matrices and a core. The factor matrices are <tensor_name>_fac_0, <tensor_name>_fac_1, ..., and the core is <tensor_name>_core. 
#
# Example 
#  cp als.py tensor1000 1000,1000,1000 100,100,100 10 1,2
#
#  - The example takes a 3-mode tensor of size 1000x1000x100 with the chunk size 100x100x100 and runs the rank-10 CP ALS algorithm with materialized matricization on two modes (mode-1,2). The materialized matricization is created by cp_mat.py. The output of the example includes factor matrices (tensor1000_fac_0, tensor1000_fac_1, and tensor1000_fac_2) and a core (tensor1000_core). 
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
	if argcnt < 5:
		print "cp_als.py <tensor_name> <tensor_size> <chunk_size> <target_rank> <modes_of_materialized_matricization> <max_iter> <debug>"
		exit(2)

	start = datetime.datetime.now();

	tensor=sys.argv[1]
	matt = tensor+'_'
	facmat = tensor+'_fac_'
	core = tensor + '_core'
	rank=int(sys.argv[4])

	debug = 0
	if argcnt == 8:
		debug = int(sys.argv[7]) 

	#print 'debug='+str(debug)

	if debug == 0:
		resfile = 'res_' + tensor
		f = open(resfile,'a')
		msg = 'cp_als start'
		f.write(msg+'\n')

	if argcnt >= 7:
		max_iter=int(sys.argv[6])
	else:
		max_iter=50

	#print 'max_iter='+str(max_iter)

	# initialize factor matrices
	facsz=sys.argv[2].split(',')
	facchunksz=sys.argv[3].split(',')
	modes = []
	if argcnt >= 6:
		modes=sys.argv[5].split(',')
	#print modes
	nmode=len(facsz)

	#facsz=[8,184,184,5632]
	#facchunksz=[8,184,184,200]
	#facchunksz=[5,3952,3020]

	# initialize factor matrices
	for i in range(len(facsz)):
		query="create array "+facmat+str(i)+"<val:double>[i=1:"+str(rank)+","+str(rank)+",0,j=1:"+str(facsz[i])+","+str(facchunksz[i])+",0]"
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
		query="copyArray(build(<val:double>[i=1:"+str(rank)+","+str(rank)+",0,j=1:"+str(facsz[i])+","+str(facchunksz[i])+",0],random()%9/10.0),"+facmat+str(i) + ")"
		if debug == 0:
			print query
			f.write(query+'\n')
			result=db.executeQuery(query,"afl")
			db.completeQuery(result.queryID)
		else:
			print query
	query="create array "+core+" <val:double>[i=1:1,1,0,j=1:"+str(rank)+","+str(rank)+",0]"
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

	#compute norm for input
	query="select sqrt(sum(pow(val,2))) from " + tensor
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
	
	fit = 0;
	fitchangetol = 1e-4;
	
	# iterate until the fit converges or for max_iter 
	for r in range(max_iter):
		fitold = fit;
		for i in range(nmode):	
			# matricized Khatri-Rao product
			query = ""
			for k in range(nmode):
				if k == i:
					continue
				if query == "":
					query=facmat+str(k)
				else:
					query="khatrirao("+query+","+facmat+str(k)+")"

			if str(i) in modes:
				query="multiply_row("+matt+str(i)+","+query+")"
			else:
				query="multiply_row(matricize("+tensor+","+str(i+1)+"),"+query+")"

			query2="multiply_row(pinverse("
			# multiplication of each factor matrix and their transpose matrix
			query1 = ""
			for k in range(nmode):
				if i == k:
					continue
				if query1 == "":
					query1="multiply_row("+facmat+str(k)+","+facmat+str(k)+")"
				else:
					query1="hadamard("+query1+",multiply_row("+facmat+str(k)+","+facmat+str(k)+"))"
			query=query2+query1+"),"+query+")"
			query="copyArray("+query+","+facmat+str(i)+")"
			if debug == 0:
				print query
				result=db.executeQuery(query,"afl")
				db.completeQuery(result.queryID)
				f.write(query+'\n')
			else:
				print query

			#normalize
			if r == 0:
				if i == nmode-1:
					query = "twoNormCoreArray("+facmat + str(i) + ","+core+")"
				else:
					query = "twoNormArray("+facmat + str(i) + ")"
			else:
				if i== nmode-1:
					query = "maxNormCoreArray("+facmat + str(i) + ","+core+")"
				else:
					query = "maxNormArray("+facmat + str(i) + ")"
			if debug == 0: 
				print query
				result=db.executeQuery(query,"afl")
				db.completeQuery(result.queryID)
				f.write(query+'\n')
			else:
				print query
				
		# frobenius norm of x - x^
		query="khatrirao("+facmat+"0,"+facmat+"1)"
		for k in range(2,nmode):
			query="khatrirao("+query+","+facmat+str(k)+")"
		query="multiply_row("+core+",transpose("+query+"))"

		if str(0) in modes:
			query="norm_residual("+matt+"0,"+query+")";
		else:
			query="norm_residual(matricize("+tensor+",1),"+query+")";
		if debug == 0:
			print query
			f.write(query+'\n');
		else:
			print query
			
		if debug == 0: 
			print query
			result=db.executeQuery(query,"afl")
			dataitem = result.array.getConstIterator(0).getChunk().getConstIterator().getItem()
			attrs = result.array.getArrayDesc().getAttributes()
			norm_residual =  scidb.getTypedValue(dataitem, attrs[0].getType())
			f.write("norm_residual="+str(norm_residual)+"\n")
			f.write("normx="+str(normx)+'\n')
			fit = 1-norm_residual/normx
			print "\n##iteration="+str(r+1)+",fit="+str(fit)+'\n'
			f.write("\n##iteration="+str(r+1)+",fit="+str(fit)+'\n')
			fitchange = abs(fitold - fit);
			db.completeQuery(result.queryID)

			if r > 1 and fitchange < fitchangetol:
				break

	end = datetime.datetime.now();
	minutes, seconds = divmod((end-start).seconds, 60)
	microsec = (end-start).microseconds
	
	if debug == 0:
		msg = 'cp_als Time elapsed: %d min %d.%d sec.' % (minutes, seconds, microsec)
		print msg+'\n\n'
		f.write(msg+'\n\n')
		f.close()

    #Disconnect from the SciDB server.
	db.disconnect() 
	sys.exit(0) #success

if __name__ == "__main__":
    main()
