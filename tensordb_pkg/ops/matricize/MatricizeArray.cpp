/*
**
* BEGIN_COPYRIGHT
*
* This file is part of TensorDB.
* Copyright (C) 2014 Mijung Kim, K. Selcuk Candan All right Reserved.
*
* Reference:
*   Mijung Kim (2014). TensorDB and Tensor-based Relational Model (TRM) for Efficient Tensor-Relational Operations. Ph.D. Thesis. Arizona State University.
*
*
* This file is part of SciDB.
* Copyright (C) 2008-2012 SciDB, Inc.
*
* SciDB is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation version 3 of the License.
*
* SciDB is distributed "AS-IS" AND WITHOUT ANY WARRANTY OF ANY KIND,
* INCLUDING ANY IMPLIED WARRANTY OF MERCHANTABILITY,
* NON-INFRINGEMENT, OR FITNESS FOR A PARTICULAR PURPOSE. See
* the GNU General Public License for the complete license terms.
*
* You should have received a copy of the GNU General Public License
* along with SciDB.  If not, see <http://www.gnu.org/licenses/>.
*
* END_COPYRIGHT
*/

/*
 * MatricizeArray.cpp
 *
 * @author: Mijung Kim (mijung.kim.1@asu.edu)
 */

#include <log4cxx/logger.h>
#include "array/DelegateArray.h"
#include "query/ops/matricize/MatricizeArray.h"

namespace scidb
{
    // Logger for operator. static to prevent visibility of variable outside of file
    static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("scidb.query.ops.matricize"));

    //
    // MatricizeArrayIterator methods
    //
    void MatricizeArrayIterator::operator ++()
    {
        if (!hasCurrent)
            throw USER_EXCEPTION(SCIDB_SE_EXECUTION, SCIDB_LE_NO_CURRENT_ELEMENT);
        Dimensions const& dims = array.dims;
        currPos[1] += dims[1].getChunkInterval();
		if(currPos[1] > dims[1].getEndMax()) {
			currPos[0] += dims[0].getChunkInterval();
			if(currPos[0] > dims[0].getEndMax()) {
				hasCurrent = false;
				return;
			}
			currPos[1] = dims[1].getStart();
		}
		return;
    }

    bool MatricizeArrayIterator::end()
    {
        return !hasCurrent;
    }

    Coordinates const& MatricizeArrayIterator::getPosition()
    {
        return currPos;
    }

    bool MatricizeArrayIterator::setPosition(Coordinates const& pos)
    {
        Dimensions const& dims = array.dims;
        for (size_t i = 0, n = currPos.size(); i < n; i++) {
            if (pos[i] < dims[i].getStart() || pos[i] > dims[i].getEndMax()) {
                return hasCurrent = false;
            }
        }
        currPos = pos;
        array.desc.getChunkPositionFor(currPos);
        return hasCurrent = true;
    }

    void MatricizeArrayIterator::reset()
    {
        Dimensions const& dims = array.dims;
        hasCurrent = true;
        for (size_t i = 0, n = currPos.size(); i < n; i++) {
            currPos[i] = dims[i].getStart();
            hasCurrent &= dims[i].getLength() != 0;
        }
    }

    /***
     * This is the main method in MatricizeArrayIterator that computes a whole output chunk
     * of the matricization and sends it up the pipeline
     */
    ConstChunk const& MatricizeArrayIterator::getChunk()
    {
        if (!hasCurrent)
            throw USER_EXCEPTION(SCIDB_SE_EXECUTION, SCIDB_LE_NO_CURRENT_ELEMENT);
        // Initialize the chunk
        if (chunk.isInitialized() && currPos == chunk.getFirstPosition(false)) {
           return chunk;
        }
        Address addr(0, currPos);
        chunk.initialize(&array, &array.desc, addr, 0);
        chunk.setRLE(false);

        boost::shared_ptr<ConstArrayIterator> inputIter = array.inputArray->getConstIterator(attr);
        boost::shared_ptr<ChunkIterator> outIter;

    	Dimensions const& inputDims = array.inputArray->getArrayDesc().getDimensions();
    	Coordinate ndims = inputDims.size();
    	Coordinates inputPos(ndims);
    	inputPos[array.rowmode] = currPos[0];

		Coordinate cn = (currPos[1]-array.jStart)/array.jChunkLen;
		Coordinate ci = cn/array.nChunk[array.colmode];
		inputPos[array.colmode] = cn%array.nChunk[array.colmode]*inputDims[array.colmode].getChunkInterval() + 1; 
    	Coordinate cj = array.restNumChunk; 
    	for(Coordinate i=0;i<ndims;i++) {
            if(i == array.rowmode || i == array.colmode) continue;

        	cj = cj/array.nChunk[i];
        	inputPos[i] = (ci/cj)*inputDims[i].getChunkInterval() + 1;

        	ci = ci%cj;
    	}
        boost::shared_ptr<Query> query(_query.lock());
        if (inputIter->setPosition(inputPos)) {
        	if (!outIter) {
            	outIter = chunk.getIterator(query, ChunkIterator::IGNORE_OVERLAPS);
        	}
        	array.matricizeChunks(inputIter, outIter, currPos, inputPos);
		}
 		if (!outIter) {
           outIter = chunk.getIterator(query, ChunkIterator::SPARSE_CHUNK);
        }
       	outIter->flush();

        return chunk;
    }


    MatricizeArrayIterator::MatricizeArrayIterator(MatricizeArray const& arr, AttributeID attrID)
    : array(arr),
      attr(attrID),
      currPos(arr.dims.size()),
      _query(arr._query)
    {
        reset();
    }

    //
    // Matricize array methods
    //

    ArrayDesc const& MatricizeArray::getArrayDesc() const
    {
        return desc;
    }

    boost::shared_ptr<ConstArrayIterator> MatricizeArray::getConstIterator(AttributeID attr) const
    {
        return result ? result->getConstIterator(attr)
            : boost::shared_ptr<ConstArrayIterator>((ConstArrayIterator*)new MatricizeArrayIterator(*this, attr));

    }

    void MatricizeArray::matricizeChunks(boost::shared_ptr<ConstArrayIterator> const& inputIter, boost::shared_ptr<ChunkIterator> const& outIter, Coordinates currPos, Coordinates inputPos) const
    {
        ConstChunk const* inputChunk = &inputIter->getChunk();
        MemChunk inputMatChunk;
        if (!inputChunk->isMaterialized()) {
            MaterializedArray::materialize(inputMatChunk, *inputChunk, MaterializedArray::DenseFormat);
            inputChunk = &inputMatChunk;
        }
        boost::shared_ptr<Query> query(_query.lock());

        void* resultData = outIter->getChunk().getData();

        PinBuffer s1(*inputChunk);
        void* inputData = inputChunk->getData();
        if (inputChunk->isRLE()) { 
           ConstRLEPayload inputPayload((char const*)inputData);
           if (inputPayload.nSegments() != 1 || inputPayload.getSegment(0).same || inputPayload.getSegment(0).length() != totalInputChunkLen) { 
              if (inputPayload.nSegments() == 1 && inputPayload.getSegment(0).same) { 
                   Value val;
                   inputPayload.getValueByIndex(val, inputPayload.getSegment(0).valueIndex);
                   if (val.isZero()) { // result is 0
                        return;
                    }
                }
                MaterializedArray::materialize(inputMatChunk, *inputChunk, MaterializedArray::DenseFormat);
                inputData = inputMatChunk.getData();
            } else { 
                inputData = inputPayload.getRawValue(inputPayload.getSegment(0).valueIndex);
            }
        }
		Coordinate ci = currPos[0];
		Coordinate cj = currPos[1];

		Coordinate iUpper = Coordinate(ci+iChunkLen) < Coordinate(iStart + iLength) ? ci+iChunkLen : iStart + iLength;
    	Coordinate jUpper = Coordinate(cj+jChunkLen) < Coordinate(jStart + jLength) ? cj+jChunkLen : jStart + jLength;

		double* result = (double*)resultData;
		double* input = (double*)inputData;
		if(inputIncChunkLen[rowmode] == jChunkLen && inputIncChunkLen[colmode] == 1) {
			memcpy(result,input,inputIncChunkLen[rowmode]*(iUpper-ci)*sizeof(double*));
				
		} else {

			for(Coordinate i=ci;i<iUpper;i++) {
       			for(Coordinate j=cj;j<jUpper;j+=inputChunkLen[colmode]) {
				
					Coordinate b = (j-cj)/inputChunkLen[colmode];
					Coordinate rest = 0;

					Coordinate a = restInputChunkLen;
        			for(Coordinate k=0;k<inputPos.size();k++) {

               			if(k == rowmode || k == colmode) continue;

        				a = a/inputChunkLen[k];

						
						Coordinate pos = b/a;
						rest += pos*inputIncChunkLen[k];

        				b = b%a;
    				}

					
					Coordinate p = (i-ci)*inputIncChunkLen[rowmode] + rest; 
					if(inputIncChunkLen[colmode] != 1) {
        				for(Coordinate k=0;k<inputChunkLen[colmode];k++) {
							Coordinate q = p+k*inputIncChunkLen[colmode];
							Coordinate y = (i-ci)*jChunkLen + (j-cj) + k;
       						result[y] = input[q];
						}
					} else {
						Coordinate y = (i-ci)*jChunkLen + (j-cj);
						memcpy((double*)(result+y),(double*)(input+p),inputChunkLen[colmode]*sizeof(double*));
					}
    			}
			}
		}
	}

    MatricizeArray::MatricizeArray(ArrayDesc aDesc, boost::shared_ptr<Array> const& in, boost::shared_ptr<Query> const& query, Coordinate rmode, Coordinate cmode)
        : desc(aDesc),
          dims(desc.getDimensions()),
	      inputArray(in),
          _query(query),
	      rowmode(rmode),
	 	  colmode(cmode)
    {
		ArrayDesc const& inputDesc = inputArray->getArrayDesc();
        Dimensions const& inputDims = inputDesc.getDimensions();
		totalInputChunkLen = 1;
		for(Coordinate i=0;i<inputDims.size();i++) {
			totalInputChunkLen *= inputDims[i].getChunkInterval();
		}
		Coordinates d(inputDims.size());
		Coordinates c(inputDims.size());
		Coordinates a(inputDims.size());
		d[inputDims.size()-1] = 1;
		size_t len = inputDims[inputDims.size()-1].getChunkInterval();
		for(Coordinate i=inputDims.size()-2;i>=0;i--) {
			d[i] = len;
			len *= inputDims[i].getChunkInterval();
		}
		for(Coordinate i=inputDims.size()-1;i>=0;i--) {
			c[i] = inputDims[i].getChunkInterval();
			a[i] = inputDims[i].getLength()/inputDims[i].getChunkInterval();
			Coordinate res = inputDims[i].getLength()%inputDims[i].getChunkInterval();
			if(res > 0) a[i] += 1;
		}
		inputIncChunkLen = d;
		inputChunkLen = c;
		nChunk = a;
       	restInputChunkLen = 1; 
       	restNumChunk = 1; 
        for(Coordinate i=0;i<inputDims.size();i++) {
            if(i == rowmode || i == colmode) continue;
            restInputChunkLen *= inputDims[i].getChunkInterval();
            restNumChunk *= nChunk[i]; 
        }

        iChunkLen = dims[0].getChunkInterval();
        jChunkLen = dims[1].getChunkInterval();
        iLength = dims[0].getLength();
        jLength = dims[1].getLength();
        iStart = dims[0].getStart();
        jStart = dims[1].getStart();

    }

} //namespace scidb
