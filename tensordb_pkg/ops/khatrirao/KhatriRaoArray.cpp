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
 * KhatriRaoArray.cpp 
 *
 * @author: Mijung Kim (mijung.kim.1@asu.edu)
 */

#include <log4cxx/logger.h>
#include "array/DelegateArray.h"
#include "query/ops/khatrirao/KhatriRaoArray.h"

namespace scidb
{
    const int ROW = 0;
    const int COL = 1;

    // Logger for operator. static to prevent visibility of variable outside of file
    static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("scidb.query.ops.khatrirao"));

    //
    // KhatriRaoArrayIterator methods
    //
    void KhatriRaoArrayIterator::operator ++()
    {
        if (!hasCurrent)
            throw USER_EXCEPTION(SCIDB_SE_EXECUTION, SCIDB_LE_NO_CURRENT_ELEMENT);
        Dimensions const& dims = array.dims;
        currPos[COL] += dims[COL].getChunkInterval();
		if(currPos[COL] > dims[COL].getEndMax()) {
			currPos[ROW] += dims[ROW].getChunkInterval();
			if(currPos[ROW] > dims[ROW].getEndMax()) {
				hasCurrent = false;
				return;
			}
			currPos[COL] = dims[COL].getStart();
		}
		return;
    }

    bool KhatriRaoArrayIterator::end()
    {
        return !hasCurrent;
    }

    Coordinates const& KhatriRaoArrayIterator::getPosition()
    {
        return currPos;
    }

    bool KhatriRaoArrayIterator::setPosition(Coordinates const& pos)
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

    void KhatriRaoArrayIterator::reset()
    {
        Dimensions const& dims = array.dims;
        hasCurrent = true;
        for (size_t i = 0, n = currPos.size(); i < n; i++) {
            currPos[i] = dims[i].getStart();
            hasCurrent &= dims[i].getLength() != 0;
        }
    }

    /***
     * This is the main method in KhatriRaoArrayIterator that computes a whole output chunk
     * of the khatri-rao product and sends it up the pipeline
     */
    ConstChunk const& KhatriRaoArrayIterator::getChunk()
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

        boost::shared_ptr<ConstArrayIterator> leftIter = array.leftArray->getConstIterator(attr);
        boost::shared_ptr<ConstArrayIterator> rightIter = array.rightArray->getConstIterator(attr);
        boost::shared_ptr<ChunkIterator> outIter;

        Coordinate cn = (currPos[COL]-array.jStart)/array.jChunkLen;
        Coordinate ci = cn/array.jRightNumChunk*array.jLeftChunkLen;
        Coordinate ck = cn%array.jRightNumChunk*array.jRightChunkLen;
        Coordinate cj = currPos[ROW]-array.iStart;
        Coordinates leftPos(2);
        Coordinates rightPos(2);
        leftPos[COL] = array.jLeftStart + ci; // /array.iRightLength;
        rightPos[COL] = array.jRightStart + ck; // %array.iRightLength;
        leftPos[ROW] = array.iLeftStart+cj;
		rightPos[ROW] = array.iRightStart+cj;
        //
        // Iterate over aligned blocks from the input, compute khatri-rao product for each i & j block
        //
        boost::shared_ptr<Query> query(_query.lock());
        if (leftIter->setPosition(leftPos) && rightIter->setPosition(rightPos)) {
            if (!outIter) {
                outIter = chunk.getIterator(query, ChunkIterator::IGNORE_OVERLAPS);
            }
            array.khatriraoChunks(leftIter, rightIter, outIter, ci, cj, ck);
        }
        if (!outIter) {
           outIter = chunk.getIterator(query, ChunkIterator::SPARSE_CHUNK);
        }
        outIter->flush();

        return chunk;
    }


    KhatriRaoArrayIterator::KhatriRaoArrayIterator(KhatriRaoArray const& arr, AttributeID attrID)
    : array(arr),
      attr(attrID),
      currPos(arr.dims.size()),
      _query(arr._query)
    {
        reset();
    }

    //
    // KhatriRao array methods
    //

    ArrayDesc const& KhatriRaoArray::getArrayDesc() const
    {
        return desc;
    }

    boost::shared_ptr<ConstArrayIterator> KhatriRaoArray::getConstIterator(AttributeID attr) const
    {
        return result ? result->getConstIterator(attr)
            : boost::shared_ptr<ConstArrayIterator>((ConstArrayIterator*)new KhatriRaoArrayIterator(*this, attr));

    }

    void KhatriRaoArray::khatriraoChunks(boost::shared_ptr<ConstArrayIterator> const& rowArrayIter,
                                       boost::shared_ptr<ConstArrayIterator> const& colArrayIter,
                                       boost::shared_ptr<ChunkIterator> const& outIter, Coordinate ci, Coordinate cj, Coordinate ck) const
    {
        ConstChunk const* leftChunk = &rowArrayIter->getChunk();
        ConstChunk const* rightChunk = &colArrayIter->getChunk();
        MemChunk leftMatChunk, rightMatChunk;
        if (!leftChunk->isMaterialized()) {
            MaterializedArray::materialize(leftMatChunk, *leftChunk, MaterializedArray::DenseFormat);
            leftChunk = &leftMatChunk;
        }
        if (!rightChunk->isMaterialized()) {
            MaterializedArray::materialize(rightMatChunk, *rightChunk, MaterializedArray::DenseFormat);
            rightChunk = &rightMatChunk;
        }
        boost::shared_ptr<Query> query(_query.lock());

        vector< boost::shared_ptr<KhatriRaoJob> > jobs(nCPUs);

        void* resultData = outIter->getChunk().getData();

        PinBuffer s1(*leftChunk);
        void* leftData = leftChunk->getData();
        if (leftChunk->isRLE()) { 
            ConstRLEPayload leftPayload((char const*)leftData);
            if (leftPayload.nSegments() != 1 || leftPayload.getSegment(0).same || leftPayload.getSegment(0).length() != iLeftChunkLen*jChunkLen) { 
                if (leftPayload.nSegments() == 1 && leftPayload.getSegment(0).same) { 
                    Value val;
                    leftPayload.getValueByIndex(val, leftPayload.getSegment(0).valueIndex);
                    if (val.isZero()) { // result is 0
                        return;
                    }
                }
                MaterializedArray::materialize(leftMatChunk, *leftChunk, MaterializedArray::DenseFormat);
                leftData = leftMatChunk.getData();
            } else { 
                leftData = leftPayload.getRawValue(leftPayload.getSegment(0).valueIndex);
            }
        }

        PinBuffer s2(*rightChunk);
        void* rightData = rightChunk->getData();
        if (rightChunk->isRLE()) { 
            ConstRLEPayload rightPayload((char const*)rightData);
            if (rightPayload.nSegments() != 1 || rightPayload.getSegment(0).same || rightPayload.getSegment(0).length() != iRightChunkLen*jChunkLen) {
                if (rightPayload.nSegments() == 1 && rightPayload.getSegment(0).same) { 
                    Value val;
                    rightPayload.getValueByIndex(val, rightPayload.getSegment(0).valueIndex);
                    if (val.isZero()) { // result is 0
                        return;
                    }
                }
                MaterializedArray::materialize(rightMatChunk, *rightChunk, MaterializedArray::DenseFormat);
                rightData = rightMatChunk.getData();
            } else { 
                rightData = rightPayload.getRawValue(rightPayload.getSegment(0).valueIndex);
            }
        }
 
        for (int i = 0; i < nCPUs; i++) {
           jobs[i] = boost::shared_ptr<KhatriRaoJob>(new KhatriRaoJob(*this, leftData, rightData, resultData, ci, cj, ck, i, query));
           queue->pushJob(jobs[i]);
        }
        int errorJob = -1;
        for (int i = 0; i < nCPUs; i++) {
            if (!jobs[i]->wait()) {
                errorJob = i;
            }
        }
        if (errorJob >= 0) {
            jobs[errorJob]->rethrow();
        }
    }

    void KhatriRaoArray::khatriraoChunks(void const* leftData,
                                       void const* rightData,
                                       void* resultData,
                                       Coordinate ci, Coordinate cj, Coordinate ck, Coordinate offset) const
    {
        Coordinate iUpper = Coordinate(ci+jLeftChunkLen) < Coordinate(jLeftStart + jLeftLength) ? ci+jLeftChunkLen : jLeftStart + jLeftLength;
        Coordinate jUpper = Coordinate(cj+iChunkLen) < Coordinate(iStart + iLength) ? cj+iChunkLen : iStart + iLength;
        Coordinate kUpper = Coordinate(ck+jRightChunkLen) < Coordinate(jRightStart + jRightLength) ? ck+jRightChunkLen : jRightStart + jRightLength;
        Coordinate step = nCPUs;
        StatisticsScope sScope(_statistics);

        //
        //  PGB: Naive sum-of-products like this has numerical
        //       stability issues. Won't fix for now.
        if (TID_DOUBLE == attrType )
        {
            double* left = (double*)leftData;
            double* right = (double*)rightData;
            double* result = (double*)resultData;
           	for (Coordinate j=cj+offset; j<jUpper; j += step)
           	{
            	for (Coordinate i=ci; i<iUpper; i += 1)
            	{
            		for (Coordinate k=ck; k<kUpper; k += 1)
            		{
                		result[(j-cj)*jChunkLen + (i-ci)*jRightChunkLen + k-ck] = 
								left[(j-cj)*jLeftChunkLen+(i-ci)] * right[(j-cj)*jRightChunkLen + (k-ck)];
                	} 	
            	}
			}
        } else if (TID_FLOAT == attrType ) {
            float* left = (float*)leftData;
            float* right = (float*)rightData;
            float* result = (float*)resultData;
            for (Coordinate j=cj+offset; j<jUpper; j += step)
            {
                for (Coordinate i=ci; i<iUpper; i += 1)
                {
                    for (Coordinate k=ck; k<kUpper; k += 1)
                    {
                        result[(j-cj)*jChunkLen + (i-ci)*jRightChunkLen + k-ck] = 
                                left[(j-cj)*jLeftChunkLen+(i-ci)] * right[(j-cj)*jRightChunkLen + (k-ck)];
                    }
                }
            }
        } else if (TID_INT32 == attrType ) {
            int32_t* left = (int32_t*)leftData;
            int32_t* right = (int32_t*)rightData;
            int32_t* result = (int32_t*)resultData;
            for (Coordinate j=cj+offset; j<jUpper; j += step)
            {
                for (Coordinate i=ci; i<iUpper; i += 1)
                {
                    for (Coordinate k=ck; k<kUpper; k += 1)
                    {
                        result[(j-cj)*jChunkLen + (i-ci)*jRightChunkLen + k-ck] = 
                                left[(j-cj)*jLeftChunkLen+(i-ci)] * right[(j-cj)*jRightChunkLen + (k-ck)];
                    }
                }
            }
        } else if (TID_INT16 == attrType ) {
            int16_t* left = (int16_t*)leftData;
            int16_t* right = (int16_t*)rightData;
            int16_t* result = (int16_t*)resultData;
            for (Coordinate j=cj+offset; j<jUpper; j += step)
            {
                for (Coordinate i=ci; i<iUpper; i += 1)
                {
                    for (Coordinate k=ck; k<kUpper; k += 1)
                    {
                        result[(j-cj)*jChunkLen + (i-ci)*jRightChunkLen + k-ck] = 
                                left[(j-cj)*jLeftChunkLen+(i-ci)] * right[(j-cj)*jRightChunkLen + (k-ck)];
                    }
                }
            }
        } else if (TID_INT64 == attrType ) {
            int64_t* left = (int64_t*)leftData;
            int64_t* right = (int64_t*)rightData;
            int64_t* result = (int64_t*)resultData;
            for (Coordinate j=cj+offset; j<jUpper; j += step)
            {
                for (Coordinate i=ci; i<iUpper; i += 1)
                {
                    for (Coordinate k=ck; k<kUpper; k += 1)
                    {
                        result[(j-cj)*jChunkLen + (i-ci)*jRightChunkLen + k-ck] = 
                                left[(j-cj)*jLeftChunkLen+(i-ci)] * right[(j-cj)*jRightChunkLen + (k-ck)];
                    }
                }
            }
        } else if (TID_UINT32 == attrType ) {
            uint32_t* left = (uint32_t*)leftData;
            uint32_t* right = (uint32_t*)rightData;
            uint32_t* result = (uint32_t*)resultData;
            for (Coordinate j=cj+offset; j<jUpper; j += step)
            {
                for (Coordinate i=ci; i<iUpper; i += 1)
                {
                    for (Coordinate k=ck; k<kUpper; k += 1)
                    {
                        result[(j-cj)*jChunkLen + (i-ci)*jRightChunkLen + k-ck] = 
                                left[(j-cj)*jLeftChunkLen+(i-ci)] * right[(j-cj)*jRightChunkLen + (k-ck)];
                    }
                }
            }
        } else if (TID_UINT16 == attrType ) {
            uint16_t* left = (uint16_t*)leftData;
            uint16_t* right = (uint16_t*)rightData;
            uint16_t* result = (uint16_t*)resultData;
            for (Coordinate j=cj+offset; j<jUpper; j += step)
            {
                for (Coordinate i=ci; i<iUpper; i += 1)
                {
                    for (Coordinate k=ck; k<kUpper; k += 1)
                    {
                        result[(j-cj)*jChunkLen + (i-ci)*jRightChunkLen + k-ck] = 
                                left[(j-cj)*jLeftChunkLen+(i-ci)] * right[(j-cj)*jRightChunkLen + (k-ck)];
                    }
                }
            }
        } else if (TID_UINT64 == attrType ) {
            uint64_t* left = (uint64_t*)leftData;
            uint64_t* right = (uint64_t*)rightData;
            uint64_t* result = (uint64_t*)resultData;
            for (Coordinate j=cj+offset; j<jUpper; j += step)
            {
                for (Coordinate i=ci; i<iUpper; i += 1)
                {
                    for (Coordinate k=ck; k<kUpper; k += 1)
                    {
                        result[(j-cj)*jChunkLen + (i-ci)*jRightChunkLen + k-ck] = 
                                left[(j-cj)*jLeftChunkLen+(i-ci)] * right[(j-cj)*jRightChunkLen + (k-ck)];
                    }
                }
            }
        } else {
            throw SYSTEM_EXCEPTION(SCIDB_SE_INFER_SCHEMA, SCIDB_LE_OP_MULTIPLY_ERROR1);
        }
    }


    KhatriRaoArray::KhatriRaoJob::KhatriRaoJob
        (KhatriRaoArray const& arr,
         void const* left,
         void const* right,
         void* result,
         Coordinate ci, Coordinate cj, Coordinate ck, Coordinate offset,
         const boost::shared_ptr<Query>& query)
        : Job(query),
          array(arr),
          leftData(left),
          rightData(right),
          resultData(result)
    {
        this->ci = ci;
        this->cj = cj;
        this->ck = ck;
        this->offset = offset;
    }

    void KhatriRaoArray::KhatriRaoJob::run()
    {
        array.khatriraoChunks(leftData, rightData, resultData, ci, cj, ck, offset);
    }

    KhatriRaoArray::KhatriRaoArray(ArrayDesc aDesc, boost::shared_ptr<Array> const& leftInputArray, boost::shared_ptr<Array> const& rightInputArray, const boost::shared_ptr<Query>& query)
        : desc(aDesc),
          dims(desc.getDimensions()),
          _query(query)
    {
        nCPUs = Sysinfo::getNumberOfCPUs();
        queue = boost::shared_ptr<JobQueue>(new JobQueue());
        threadPool = boost::shared_ptr<ThreadPool>(new ThreadPool(nCPUs, queue));
        threadPool->start();

        ArrayDesc const& leftDesc = leftInputArray->getArrayDesc();
        ArrayDesc const& rightDesc = rightInputArray->getArrayDesc();
        Dimensions const& leftDims = leftDesc.getDimensions();
        Dimensions const& rightDims = rightDesc.getDimensions();
        iChunkLen = dims[ROW].getChunkInterval();
        jChunkLen = dims[COL].getChunkInterval();
        iLength = dims[ROW].getLength();
        jLength = dims[COL].getLength();
        iStart = dims[ROW].getStart();
        jStart = dims[COL].getStart();
        iLeftLength = leftDims[ROW].getLength();
        jLeftLength = leftDims[COL].getLength();
        iRightLength = rightDims[ROW].getLength();
        jRightLength = rightDims[COL].getLength();
        iLeftChunkLen = leftDims[ROW].getChunkInterval();
        jLeftChunkLen = leftDims[COL].getChunkInterval();
        iRightChunkLen = rightDims[ROW].getChunkInterval();
        jRightChunkLen = rightDims[COL].getChunkInterval();
		jRightNumChunk = rightDims[COL].getLength()/rightDims[COL].getChunkInterval();
		Coordinate res = rightDims[COL].getLength()%rightDims[COL].getChunkInterval();
		if(res > 0) jRightNumChunk += 1;
        iLeftStart = leftDims[ROW].getStart();
        iRightStart = rightDims[ROW].getStart();
        jLeftStart = leftDims[COL].getStart();
        jRightStart = rightDims[COL].getStart();
        attrType = desc.getAttributes()[0].getType();

        leftArray = leftInputArray;
        rightArray = rightInputArray;
    }

} //namespace scidb
