<How to run TensorDB>

- Prerequisites
1. SciDB 12.12 (www.scidb.org)  
2. Armadillo C++ linear algebra library (http://arma.sourceforge.net/)
3. FFTW library (http://www.fftw.org/): a shared library (.so) file should be installed

- Files 
1. ErrorCodes.h

:Append the error codes in the file to <scidb directory>/include/system/ErrorCodes.h (e.g., scidb-12.12.0/include/system/ErrorCodes.h) 

2. LongErrorsList.h

:Append the error messsages in the file to <scidb directory>/src/system/LongErrorsList.h 

3. ops/

  - BuildInOps.inc  
  : Append all to <scidb director>/src/query/ops/BuildInOps.inc 

  - CMakeLists.txt
  : Append all within "set()" of <scidb director>/src/query/ops/CMakeLists.txt
  : Put libraries in target_link_libraries
   E.g., target_link_libraries(ops_lib compression_lib network_lib MurmurHash_lib /usr/lib/libarmadillo.so /usr/lib/libfftw3.so) 

  - Copy <operator directories> to <scidb directory>/src/query/ops    


- 1,2, and 3 steps are done, compile and restart SciDB.
> <scidb directory>/sudo make install
> scidb.py stopall <dbname>
> scidb.py initall <dbname>
> scidb.py startall <dbname>

- Create tensor and import data using SciDB ('create array' annd 'load' commands). Note: when creating a tensor, the attribute name for each element value of the tensor must be 'val'. 

// example for creating an random tensor of size 100x100x100
iquery -naq "create array rand<val:double>[i=1:100,100,0,j=1:100,100,0,k=1:100,100,0]"
iquery -naq "store(build(<val:double>[i=1:100,100,0,j=1:100,100,0,k=1:100,100,0],random()%2/1),rand)"
////////////////////////////////////////////////////////////

- Run algorithms/*.py 

