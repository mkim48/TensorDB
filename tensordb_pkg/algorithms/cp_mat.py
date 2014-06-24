#!/usr/bin/python
#
# cp_mat.py creates materialized matricizations on the modes of the tensor. The materialized matricizations are used to speed up the cp als algorithm (cp_als.py).
#
# Parameters:
#  - tensor name: input tensor
#  - tensor size: input tensor size, e.g., 1000,1000,1000 for a tensor of size 1000x1000x1000
#  - chunk size: chunk size for the input tensor, e.g., 100,100,100 for a chunk size 100x100x100
#  - matricized modes (optional): modes which matricize, e.g., 1,2 - mode-1 and 2 to matricize
#  - debug (optional): 1 for debugging mode (not running the operations but showing the commands), default value: 0
#
# Output: 
#  - The outputs include materialized matricization, <input tensor name>_0, <input tensor name>_1, etc.
#
# Example: 
#	cp mat.py tensor1000 1000,1000,1000 100,100,100 10 1,2
#
#   - The example takes a 3-mode tensor of size 1000x1000x1000 with the chunk size 100x100x100 and materialize the matricization of the tensor on two modes (mode-1,2). The output of the example includes materialized matricizations (tensor1000_0, tensor1000_1).
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
		print "cp_mat.py <tensor_name> <tensor_size> <chunk_size> <matricized_modes> <debug>"
		exit(2)

	start = datetime.datetime.now();
	tensor=sys.argv[1]
	matt = tensor + '_'
	sz=sys.argv[2].split(',')
	print sz
	chunksz=sys.argv[3].split(',')
	modes =sys.argv[4].split(',')
	nmode = len(sz)

	debug = 0
	if argcnt == 6:
		debug = 1

	print 'debug='+str(debug)

	if debug == 0:
		resfile = 'res_' + tensor
		f = open(resfile,'a')
		msg = 'matricize start'
		f.write(msg+'\n')

	tot_sz = 1
	tot_chunksz= 1
	for i in range(nmode):
		tot_sz = tot_sz*int(sz[i])
		tot_chunksz = tot_chunksz*int(chunksz[i])

	mode_sz = [0]*nmode
	mode_chunksz = [0]*nmode
	for i in range(nmode):
		mode_sz[i] = tot_sz/int(sz[i])
		mode_chunksz[i] = tot_chunksz/int(chunksz[i])

	for i in range(nmode):
		if str(i) in modes:
			if i == nmode-1:
				l = nmode-2 #second last dimension
			else:
				l = nmode-1  #last dimension

			# perform matricize operation and store the result array in the created array for the i-mode matricized tensor
			query="store(matricise("+tensor+","+str(i)+","+str(l)+"),"+matt+str(i) + ")"

			if debug == 0:
				f.write(query+'\n')
				result=db.executeQuery(query,"afl")
				db.completeQuery(result.queryID)
			else:
				print query

	end = datetime.datetime.now();
	minutes, seconds = divmod((end-start).seconds, 60)
	microsec = (end-start).microseconds

	if debug == 0:
		msg = 'matricize Time elapsed: %d min %d.%d sec.' % (minutes, seconds, microsec)
		f.write(msg+'\n\n')
		f.close()

    #Disconnect from the SciDB server.
	db.disconnect() 
	sys.exit(0) #success

if __name__ == "__main__":
    main()
