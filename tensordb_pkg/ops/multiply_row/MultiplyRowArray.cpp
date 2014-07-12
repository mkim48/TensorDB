/*
**
* BEGIN_COPYRIGHT
*
* This file is part of SciDB.
* Copyright (C) 2008-2012 SciDB, Inc.
*
* Reference:
*   Mijung Kim (2014). TensorDB and Tensor-based Relational Model (TRM) for Efficient Tensor-Relational Operations. Ph.D. Thesis. Arizona State University.
*
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
 * MultiplyRowArray.cpp
 *
 * @author: Mijung Kim (mijung.kim.1@asu.edu)
 */
#include <unistd.h>
#include <sys/time.h>
#include <log4cxx/logger.h>
#include <boost/scope_exit.hpp>
#include "array/DelegateArray.h"
#include "query/ops/multiply_row/MultiplyRowArray.h"
#include "network/NetworkManager.h"
#include "query/FunctionDescription.h"
#include "query/Statistics.h"

namespace scidb
{
    const int ROW = 0;
    const int COL = 1;

    const size_t SPARSE_MULTIPLICATION_THRESHOLD = 10;

    // Logger for operator. static to prevent visibility of variable outside of file
    static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("scidb.query.ops.multiply_row"));

    const size_t incrementalMultiplyRowThreshold = 1024*1024;

    //
    // MultiplyRowArrayIterator methods
    //
    void MultiplyRowArrayIterator::operator ++()
    {
        if (!hasCurrent)
            throw USER_EXCEPTION(SCIDB_SE_EXECUTION, SCIDB_LE_NO_CURRENT_ELEMENT);
        Dimensions const& dims = array.dims;
        int i = currPos.size();
        while (true) {
            i -= 1;
            currPos[i] += dims[i].getChunkInterval();
            if (currPos[i] > dims[i].getEndMax()) {
                if (i == 0) {
                    hasCurrent = false;
                    return;
                }
                currPos[i] = dims[i].getStart();
            } else {
                if (array.isSelfChunk(currPos)) {
                    hasCurrent = true;
                    return;
                }
                i = currPos.size();
            }
        }
    }

    bool MultiplyRowArrayIterator::end()
    {
        return !hasCurrent;
    }

    Coordinates const& MultiplyRowArrayIterator::getPosition()
    {
        return currPos;
    }

    bool MultiplyRowArrayIterator::setPosition(Coordinates const& pos)
    {
        Dimensions const& dims = array.dims;
        for (size_t i = 0, n = currPos.size(); i < n; i++) {
            if (pos[i] < dims[i].getStart() || pos[i] > dims[i].getEndMax()) {
                return hasCurrent = false;
            }
        }
        currPos = pos;
        array.desc.getChunkPositionFor(currPos);
        return hasCurrent = array.isSelfChunk(currPos);
    }

    void MultiplyRowArrayIterator::reset()
    {
        Dimensions const& dims = array.dims;
        hasCurrent = true;
        for (size_t i = 0, n = currPos.size(); i < n; i++) {
            currPos[i] = dims[i].getStart();
            hasCurrent &= dims[i].getLength() != 0;
        }
        if (hasCurrent && !array.isSelfChunk(currPos)) {
            ++(*this);
        }
    }

    /***
     * This is the main method in MultiplyRowArrayIterator that computes a whole output chunk
     * of the multiplication and sends it up the pipeline
     */
    ConstChunk const& MultiplyRowArrayIterator::getChunk()
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

        boost::shared_ptr<ConstArrayIterator> rowIter = array.leftArray->getConstIterator(attr);
        boost::shared_ptr<ConstArrayIterator> colIter = array.rightArray->getConstIterator(attr);
        boost::shared_ptr<ChunkIterator> outIter;

        Coordinate ci = currPos[ROW];
        Coordinate cj = currPos[COL];
        Coordinates rowPos(2);
        Coordinates colPos(2);
        rowPos[ROW] = ci;
        colPos[ROW] = cj;

        //
        // Iterate over aligned blocks from the input, compute partial multiply for each i & j block
        //
        boost::shared_ptr<Query> query(_query.lock());
        for (Coordinate ck = 0; ck < Coordinate(array.kLength); ck += array.kChunkLen)
        {
            rowPos[COL] = array.kStartLeft + ck;
            colPos[COL] = array.kStartRight + ck;
            if (rowIter->setPosition(rowPos) && colIter->setPosition(colPos)) {
                if (!outIter) {
                   outIter = chunk.getIterator(query, ChunkIterator::IGNORE_OVERLAPS);
                }
                array.multiplyChunks(rowIter, colIter, outIter, ci, cj, ck);
            }
        }
        if (!outIter) {
           outIter = chunk.getIterator(query, ChunkIterator::SPARSE_CHUNK);
        }
        outIter->flush();

        return chunk;
    }


    MultiplyRowArrayIterator::MultiplyRowArrayIterator(MultiplyRowArray const& arr, AttributeID attrID)
    : array(arr),
      attr(attrID),
      currPos(arr.dims.size()),
      _query(arr._query)
    {
        reset();
    }

    //
    // Mutliply array methods
    //

    ArrayDesc const& MultiplyRowArray::getArrayDesc() const
    {
        return desc;
    }

    boost::shared_ptr<Array> MultiplyRowArray::doIterativeMultiplyRow(AttributeID attr)
    {
        boost::shared_ptr<Array> result = boost::shared_ptr<Array>(new MemArray(desc));
        boost::shared_ptr<ArrayIterator> outArrayIter = result->getIterator(0);
        {
            map<Coordinates, boost::shared_ptr<ChunkIterator>, CoordinatesLess> chunkIterators;
            Coordinates pos(2);
            size_t iterNo = 0;
            InstanceID neighborInstance = (instanceId + 1) % nInstances;
            timeval begin, end;
            time_t delta;
            boost::shared_ptr<Query> query(_query.lock());

            BOOST_SCOPE_EXIT((&chunkIterators)) 
            {
                // This is finally block: always executed in case of normal or abnormal method termination
                // Flush iterators to unpin chunks
                for (map<Coordinates, boost::shared_ptr<ChunkIterator>, CoordinatesLess>::iterator ci = chunkIterators.begin();
                     ci != chunkIterators.end();
                     ++ci)
                {
                    ci->second->flush();
                }
            }
            BOOST_SCOPE_EXIT_END

            if (iLength > jLength) { // shift left matrix
                while (true) {
                    uint64_t numMultipliedChunks = 0;
                    gettimeofday(&begin, NULL);
                    if (iterNo != 0) {
                        delta = (begin.tv_sec - end.tv_sec)*1000000 + (begin.tv_usec - end.tv_usec);
                        LOG4CXX_DEBUG(logger, "Shift left matrix: " << delta << " microseconds");
                    }

                    boost::shared_ptr<ConstArrayIterator> leftIterator = leftArray->getConstIterator(attr);
                    boost::shared_ptr<ConstArrayIterator> rightIterator = rightArray->getConstIterator(attr);
                    //printf("Shift left matrix iteration %d:\n", (int)iterNo);
                    while (!rightIterator->end()) {
                        Coordinates const& rightPos = rightIterator->getPosition();
                        Coordinate cj = rightPos[ROW];
                        Coordinate ck = rightPos[COL];
                        for (uint64_t i = 0; i < iLength; i += iChunkLen)
                        {
                            Coordinate ci = iStart + i;
                            pos[ROW] = ci;
                            pos[COL] = ck;
                            if (leftIterator->setPosition(pos)) {
                                pos[ROW] = cj;
                                boost::shared_ptr<ChunkIterator>& outIter = chunkIterators[pos];
                                if (!outIter) {
                                    Chunk& outChunk = outArrayIter->newChunk(pos);
                                    outChunk.setRLE(false);
                                    outIter = outChunk.getIterator(query, ChunkIterator::IGNORE_OVERLAPS);
                                }
                                //printf("[%lld,%lld] * [%lld,%lld]\n", ci, ck, ck, cj);
                                multiplyChunks(leftIterator, rightIterator, outIter, ci, cj, ck);
                                numMultipliedChunks ++;
                            }
                        }
                        ++(*rightIterator);
                    }
                    gettimeofday(&end, NULL);
                    delta = (end.tv_sec - begin.tv_sec)*1000000 + (end.tv_usec - begin.tv_usec);
                    LOG4CXX_DEBUG(logger, "Iteration " << (iterNo+1) << ": " << delta << " microseconds; multiplied "<<numMultipliedChunks<<" chunk pairs");
                    if (++iterNo == nInstances) {
                        break;
                    }
                    leftArray = redistribute(leftArray, query, psLocalInstance, "", neighborInstance);
                }
            } else { // Shift right matrix
                while (true) {
                    uint64_t numMultipliedChunks = 0;
                    gettimeofday(&begin, NULL);
                    if (iterNo != 0) {
                        delta = (begin.tv_sec - end.tv_sec)*1000000 + (begin.tv_usec - end.tv_usec);
                        LOG4CXX_DEBUG(logger, "Shift right matrix: " << delta << " microseconds");
                    }
                    boost::shared_ptr<ConstArrayIterator> leftIterator = leftArray->getConstIterator(attr);
                    boost::shared_ptr<ConstArrayIterator> rightIterator = rightArray->getConstIterator(attr);
                    //printf("Shift right matrix iteration %d:\n", (int)iterNo);
                    while (!leftIterator->end()) {
                        Coordinates const& leftPos = leftIterator->getPosition();
                        Coordinate ci = leftPos[ROW];
                        Coordinate ck = leftPos[COL];
                        for (uint64_t j = 0; j < jLength; j += jChunkLen)
                        {
                            Coordinate cj = jStart + j;
                            pos[COL] = ck;
                            pos[ROW] = cj;
                            if (rightIterator->setPosition(pos)) {
                                pos[ROW] = ci;
                                boost::shared_ptr<ChunkIterator>& outIter = chunkIterators[pos];
                                if (!outIter) {
                                    Chunk& outChunk = outArrayIter->newChunk(pos);
                                    outChunk.setRLE(false);
                                    outIter = outChunk.getIterator(query, ChunkIterator::IGNORE_OVERLAPS);
                                }
                                //printf("[%lld,%lld] * [%lld,%lld]\n", ci, ck, ck, cj);
                                multiplyChunks(leftIterator, rightIterator, outIter, ci, cj, ck);
                                numMultipliedChunks ++;
                            }
                        }
                        ++(*leftIterator);
                    }
                    gettimeofday(&end, NULL);
                    delta = (end.tv_sec - begin.tv_sec)*1000000 + (end.tv_usec - begin.tv_usec);
                    LOG4CXX_DEBUG(logger, "Iteration " << (iterNo+1) << ": " << delta << " microseconds; multiplied "<<numMultipliedChunks<<" chunk pairs");
                    if (++iterNo == nInstances) {
                        break;
                    }
                    rightArray = redistribute(rightArray, query, psLocalInstance, "", neighborInstance);
                }
            }
            leftArray.reset();
            rightArray.reset();
        }
        return result;
    }

    boost::shared_ptr<ConstArrayIterator> MultiplyRowArray::getConstIterator(AttributeID attr) const
    {
        return result ? result->getConstIterator(attr)
            : boost::shared_ptr<ConstArrayIterator>(isSparse
                                             ? (ConstArrayIterator*)new MultiplyRowSparseArrayIterator(*this, attr)
                                             : (ConstArrayIterator*)new MultiplyRowArrayIterator(*this, attr));

    }

    void MultiplyRowArray::multiplyChunks(boost::shared_ptr<ConstArrayIterator> const& rowArrayIter,
                                       boost::shared_ptr<ConstArrayIterator> const& colArrayIter,
                                       boost::shared_ptr<ChunkIterator> const& outIter, Coordinate ci, Coordinate cj, Coordinate ck) const
    {
        ConstChunk const* leftChunk = &rowArrayIter->getChunk();
        ConstChunk const* rightChunk = &colArrayIter->getChunk();
        if (isBuiltinType(attrType) && leftChunk->getSize()*SPARSE_MULTIPLICATION_THRESHOLD < leftChunk->getNumberOfElements(false)*sizeof(double)
            && rightChunk->getSize()*SPARSE_MULTIPLICATION_THRESHOLD < rightChunk->getNumberOfElements(false)*sizeof(double))
        {
            SparseRowMatrix left(ci, iChunkLen);
            { 
                boost::shared_ptr<ConstChunkIterator> chunkIterator = leftChunk->getConstIterator(ChunkIterator::IGNORE_EMPTY_CELLS|ChunkIterator::IGNORE_OVERLAPS|ChunkIterator::IGNORE_DEFAULT_VALUES); 
                while (!chunkIterator->end()) { 
                    Value const& v = chunkIterator->getItem();
                    if (!v.isZero()) {                         
                        Coordinates const& pos = chunkIterator->getPosition();
                        left.add(pos[ROW], pos[COL], ValueToDouble(attrType, v));
                    }
                    ++(*chunkIterator);
                }
            }
            SparseColumnMatrix right(cj, jChunkLen);
            { 
                boost::shared_ptr<ConstChunkIterator> chunkIterator = rightChunk->getConstIterator(ChunkIterator::IGNORE_EMPTY_CELLS|ChunkIterator::IGNORE_OVERLAPS|ChunkIterator::IGNORE_DEFAULT_VALUES); 
                while (!chunkIterator->end()) { 
                    Value const& v = chunkIterator->getItem();
                    if (!v.isZero()) {                         
                        Coordinates const& pos = chunkIterator->getPosition();
                        right.add(pos[COL], pos[ROW], ValueToDouble(attrType, v));
                    }
                    ++(*chunkIterator);
                }
            }
            Coordinate iUpper = Coordinate(ci+iChunkLen) < Coordinate(iStart + iLength) ? ci+iChunkLen : iStart + iLength;
            Coordinate jUpper = Coordinate(cj+jChunkLen) < Coordinate(jStart + jLength) ? cj+jChunkLen : jStart + jLength;
            Coordinates outCoords(2);
            Value value(TypeLibrary::getType(attrType));
            
            for (Coordinate i=ci; i<iUpper; i++)
            {
                for (Coordinate j=cj; j<jUpper; j++)
                {
                    SparseMatrixCell* row = left[i];
                    SparseMatrixCell* col = right[j];
                    double partialProd = 0.0;
                    if (row != NULL && col != NULL) {
                        while (true) {
                            if (col->coord > row->coord) {
                                if ((col = col->next) == NULL) {
                                    break;
                                }
                            } else if (col->coord < row->coord) {
                                if ((row = row->next) == NULL) {
                                    break;
                                }
                            } else {
                                partialProd += row->value * col->value;
                                if ((col = col->next) == NULL || (row = row->next) == NULL) {
                                    break;
                                }
                            }
                        }
                    }
                    if (partialProd != 0) {
                        outCoords[ROW] = i;
                        outCoords[COL] = j;
                        bool rc = outIter->setPosition(outCoords);
                        if (!rc) assert(false);
                        if (!outIter->isEmpty()) { 
                            Value& v = outIter->getItem();
                            partialProd += ValueToDouble(attrType, v);
                        }
                        DoubleToValue(attrType, partialProd, value);
                        outIter->writeItem(value);
                    }
                }
            }
            return;
        }    
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
//        if (leftChunk.isPlain() && rightChunk.isPlain()) 
        if (true)
        {
            vector< boost::shared_ptr<MultiplyRowJob> > jobs(nCPUs);

            void* resultData = outIter->getChunk().getData();

            PinBuffer s1(*leftChunk);
            void* leftData = leftChunk->getData();
            if (leftChunk->isRLE()) { 
                ConstRLEPayload leftPayload((char const*)leftData);
                if (leftPayload.nSegments() != 1 || leftPayload.getSegment(0).same || leftPayload.getSegment(0).length() != iChunkLen*kChunkLen) { 
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
                if (rightPayload.nSegments() != 1 || rightPayload.getSegment(0).same || rightPayload.getSegment(0).length() != jChunkLen*kChunkLen) {
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
               jobs[i] = boost::shared_ptr<MultiplyRowJob>(new MultiplyRowJob(*this, leftData, rightData, resultData, ci, cj, ck, i, query));
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
        } else {
            Coordinate iUpper = Coordinate(ci+iChunkLen) < Coordinate(iStart + iLength) ? ci+iChunkLen : iStart + iLength;
            Coordinate jUpper = Coordinate(cj+jChunkLen) < Coordinate(jStart + jLength) ? cj+jChunkLen : jStart + jLength;
            Coordinate kUpper = Coordinate(ck+kChunkLen) < Coordinate(kStartLeft + kLength) ? ck+kChunkLen : kStartLeft + kLength;

            boost::shared_ptr<ConstChunkIterator> rowIter = leftChunk->getConstIterator(ChunkIterator::IGNORE_EMPTY_CELLS|ChunkIterator::IGNORE_NULL_VALUES|ChunkIterator::IGNORE_OVERLAPS);
            boost::shared_ptr<ConstChunkIterator> colIter = rightChunk->getConstIterator(ChunkIterator::IGNORE_EMPTY_CELLS|ChunkIterator::IGNORE_NULL_VALUES|ChunkIterator::IGNORE_OVERLAPS);

            Value multiplyResultValue;
            Value* multiplyResult = &multiplyResultValue;
            vector<Value> operandValues(5, Value(TypeLibrary::getType(attrType)));
            Value* operands[5] = {&operandValues[0], &operandValues[1], &operandValues[2],
                                  &operandValues[3], &operandValues[4]};

            Coordinates outCoords(2);
            Coordinates rowCoords(2);
            Coordinates colCoords(2);

            for (Coordinate i=ci; i<iUpper; i += 1)
            {
                for (Coordinate j=cj; j<jUpper; j += 1)
                {
                    outCoords[ROW] = i;
                    outCoords[COL] = j;

                    rowCoords[ROW] = i;
                    colCoords[ROW] = j;

                    bool empty = true;

                    for (Coordinate k=ck; k<kUpper; k++)
                    {
                        rowCoords[COL] = kStartLeft + k;
                        colCoords[COL] = kStartRight + k;

                        // get items from the chunks
                        if (rowIter->setPosition(rowCoords) && colIter->setPosition(colCoords))
                        {
                            *operands[0] = rowIter->getItem();
                            *operands[1] = colIter->getItem();
                            if (attrType != multiplyResultType) {
                                multiply((const Value**)&operands[0], multiplyResult, NULL);
                                if (empty) {
                                    convert((const Value**)&multiplyResult, operands[2], NULL);
                                    empty = false;
                                } else {
                                    convert((const Value**)&multiplyResult, operands[3], NULL);
                                    add((const Value**)&operands[2], operands[4], NULL);
                                    *operands[2] = *operands[4];
                                }
                            } else {
                                if (empty) {
                                    multiply((const Value**)&operands[0], operands[2], NULL);
                                    empty = false;
                                } else {
                                    multiply((const Value**)&operands[0], operands[3], NULL);
                                    add((const Value**)&operands[2], operands[4], NULL);
                                    *operands[2] = *operands[4];
                                }
                            }
                        }
                    }
                    // store the partial product in the output matrix
                    // Write the value to the output chunk
                    if (!empty) {
                        const bool rc = outIter->setPosition(outCoords);
                        if (!rc) assert(false);
                        if (!outIter->isEmpty()) {
                            *operands[3] = outIter->getItem();
                            add((const Value**)&operands[2], operands[4], NULL);
                            *operands[2] = *operands[4];
                        }
                        outIter->writeItem(*operands[2]);
                    }
                }
            }
        }
    }

    void MultiplyRowArray::multiplyChunks(void const* leftData,
                                       void const* rightData,
                                       void* resultData,
                                       Coordinate ci, Coordinate cj, Coordinate ck, Coordinate offset) const
    {
        Coordinate iUpper = Coordinate(ci+iChunkLen) < Coordinate(iStart + iLength) ? ci+iChunkLen : iStart + iLength;
        Coordinate jUpper = Coordinate(cj+jChunkLen) < Coordinate(jStart + jLength) ? cj+jChunkLen : jStart + jLength;
        Coordinate kUpper = Coordinate(ck+kChunkLen) < Coordinate(kStartLeft + kLength) ? ck+kChunkLen : kStartLeft + kLength;
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
            for (Coordinate i=ci + offset; i<iUpper; i += step)
            {
                for (Coordinate j=cj; j<jUpper; j += 1)
                {
                    double partialProd = 0.0;
                    for (Coordinate k=ck; k<kUpper; k++)
                    {
                        double v1 = left[(i-ci)*kChunkLen + (k-ck)];
                        double v2 = right[(j-cj)*kChunkLen + (k-ck)];
                        // add the product to the partial product
                        partialProd += v1 * v2;
                    }
                    result[(i-ci)*jChunkLen + (j-cj)] += partialProd;
                }
            }
        } else if (TID_FLOAT == attrType ) {
            float* left = (float*)leftData;
            float* right = (float*)rightData;
            float* result = (float*)resultData;
            for (Coordinate i=ci + offset; i<iUpper; i += step)
            {
                for (Coordinate j=cj; j<jUpper; j += 1)
                {
                    float partialProd = 0.0;
                    for (Coordinate k=ck; k<kUpper; k++)
                    {
                        float v1 = left[(i-ci)*kChunkLen + (k-ck)];
                        float v2 = right[(j-cj)*kChunkLen + (k-ck)];
                        // add the product to the partial product
                        partialProd += v1 * v2;
                    }
                    result[(i-ci)*jChunkLen + (j-cj)] += partialProd;
                }
            }
        } else if (TID_INT32 == attrType ) {
            int32_t* left = (int32_t*)leftData;
            int32_t* right = (int32_t*)rightData;
            int32_t* result = (int32_t*)resultData;
            for (Coordinate i=ci + offset; i<iUpper; i += step)
            {
                for (Coordinate j=cj; j<jUpper; j += 1)
                {
                    int32_t partialProd = 0.0;
                    for (Coordinate k=ck; k<kUpper; k++)
                    {
                        int32_t v1 = left[(i-ci)*kChunkLen + (k-ck)];
                        int32_t v2 = right[(j-cj)*kChunkLen + (k-ck)];
                        // add the product to the partial product
                        partialProd += v1 * v2;
                    }
                    result[(i-ci)*jChunkLen + (j-cj)] += partialProd;
                }
            }
        } else if (TID_INT16 == attrType ) {
            int16_t* left = (int16_t*)leftData;
            int16_t* right = (int16_t*)rightData;
            int16_t* result = (int16_t*)resultData;
            for (Coordinate i=ci + offset; i<iUpper; i += step)
            {
                for (Coordinate j=cj; j<jUpper; j += 1)
                {
                    int16_t partialProd = 0.0;
                    for (Coordinate k=ck; k<kUpper; k++)
                    {
                        int16_t v1 = left[(i-ci)*kChunkLen + (k-ck)];
                        int16_t v2 = right[(j-cj)*kChunkLen + (k-ck)];
                        // add the product to the partial product
                        partialProd += v1 * v2;
                    }
                    result[(i-ci)*jChunkLen + (j-cj)] += partialProd;
                }
            }
        } else if (TID_INT64 == attrType ) {
            int64_t* left = (int64_t*)leftData;
            int64_t* right = (int64_t*)rightData;
            int64_t* result = (int64_t*)resultData;
            for (Coordinate i=ci + offset; i<iUpper; i += step)
            {
                for (Coordinate j=cj; j<jUpper; j += 1)
                {
                    int64_t partialProd = 0.0;
                    for (Coordinate k=ck; k<kUpper; k++)
                    {
                        int64_t v1 = left[(i-ci)*kChunkLen + (k-ck)];
                        int64_t v2 = right[(j-cj)*kChunkLen + (k-ck)];
                        // add the product to the partial product
                        partialProd += v1 * v2;
                    }
                    result[(i-ci)*jChunkLen + (j-cj)] += partialProd;
                }
            }
        } else if (TID_UINT32 == attrType ) {
            uint32_t* left = (uint32_t*)leftData;
            uint32_t* right = (uint32_t*)rightData;
            uint32_t* result = (uint32_t*)resultData;
            for (Coordinate i=ci + offset; i<iUpper; i += step)
            {
                for (Coordinate j=cj; j<jUpper; j += 1)
                {
                    uint32_t partialProd = 0.0;
                    for (Coordinate k=ck; k<kUpper; k++)
                    {
                        uint32_t v1 = left[(i-ci)*kChunkLen + (k-ck)];
                        uint32_t v2 = right[(j-cj)*kChunkLen + (k-ck)];
                        // add the product to the partial product
                        partialProd += v1 * v2;
                    }
                    result[(i-ci)*jChunkLen + (j-cj)] += partialProd;
                }
            }
        } else if (TID_UINT16 == attrType ) {
            uint16_t* left = (uint16_t*)leftData;
            uint16_t* right = (uint16_t*)rightData;
            uint16_t* result = (uint16_t*)resultData;
            for (Coordinate i=ci + offset; i<iUpper; i += step)
            {
                for (Coordinate j=cj; j<jUpper; j += 1)
                {
                    uint16_t partialProd = 0.0;
                    for (Coordinate k=ck; k<kUpper; k++)
                    {
                        uint16_t v1 = left[(i-ci)*kChunkLen + (k-ck)];
                        uint16_t v2 = right[(j-cj)*kChunkLen + (k-ck)];
                        // add the product to the partial product
                        partialProd += v1 * v2;
                    }
                    result[(i-ci)*jChunkLen + (j-cj)] += partialProd;
                }
            }
        } else if (TID_UINT64 == attrType ) {
            uint64_t* left = (uint64_t*)leftData;
            uint64_t* right = (uint64_t*)rightData;
            uint64_t* result = (uint64_t*)resultData;
            for (Coordinate i=ci + offset; i<iUpper; i += step)
            {
                for (Coordinate j=cj; j<jUpper; j += 1)
                {
                    uint64_t partialProd = 0.0;
                    for (Coordinate k=ck; k<kUpper; k++)
                    {
                        uint64_t v1 = left[(i-ci)*kChunkLen + (k-ck)];
                        uint64_t v2 = right[(j-cj)*kChunkLen + (k-ck)];
                        // add the product to the partial product
                        partialProd += v1 * v2;
                    }
                    result[(i-ci)*jChunkLen + (j-cj)] += partialProd;
                }
            }
        } else {
            throw SYSTEM_EXCEPTION(SCIDB_SE_INFER_SCHEMA, SCIDB_LE_OP_MULTIPLY_ERROR1);
        }
    }


    MultiplyRowArray::MultiplyRowJob::MultiplyRowJob
        (MultiplyRowArray const& arr,
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

    void MultiplyRowArray::MultiplyRowJob::run()
    {
        array.multiplyChunks(leftData, rightData, resultData, ci, cj, ck, offset);
    }

    MultiplyRowArray::MultiplyRowArray(ArrayDesc aDesc, boost::shared_ptr<Array> const& leftInputArray,
                                 boost::shared_ptr<Array> const& rightInputArray, const boost::shared_ptr<Query>& query, Algorithm algorithm)
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
        kChunkLen = leftDims[COL].getChunkInterval();
        kLength = leftDims[COL].getLength();
        kStartLeft = leftDims[COL].getStart();
        kStartRight = rightDims[COL].getStart();
        iStart = leftDims[ROW].getStart();
        jStart = rightDims[ROW].getStart();
        NetworkManager& netMgr = *NetworkManager::getInstance();
        instanceId = query->getInstanceID();
        nInstances = query->getInstancesCount();
        boost::shared_ptr<ConstArrayIterator> leftIterator = leftInputArray->getConstIterator(0);
        attrType = desc.getAttributes()[0].getType();

        vector<TypeId> inputTypes(2);
        inputTypes[0] = attrType;
        inputTypes[1] = attrType;
        FunctionDescription functionDesc;
        vector<FunctionPointer> converters;
        if (!FunctionLibrary::getInstance()->findFunction("*", inputTypes, functionDesc, converters, false)
            || converters.size() != 0)
            throw USER_EXCEPTION(SCIDB_SE_EXECUTION, SCIDB_LE_OPERATION_NOT_FOUND) << "*" << attrType;
        multiply = functionDesc.getFuncPtr();
        multiplyResultType = functionDesc.getOutputArg();
        convert = (multiplyResultType != attrType)
            ? FunctionLibrary::getInstance()->findConverter(multiplyResultType, attrType, false)
            : NULL;
        if (!FunctionLibrary::getInstance()->findFunction("+", inputTypes, functionDesc, converters, false)
                || converters.size() != 0 || functionDesc.getOutputArg() != attrType)
            throw USER_EXCEPTION(SCIDB_SE_EXECUTION, SCIDB_LE_OPERATION_NOT_FOUND) << "+" << attrType;
        add = functionDesc.getFuncPtr();
        isSparse = isBuiltinType(attrType)
            && (algorithm == Sparse 
                || (algorithm == Auto
                    && !leftIterator->end() && leftIterator->getChunk().isSparse()
                    && leftInputArray->getArrayDesc().getAttributes()[0].getDefaultValue().isZero()
                    && rightInputArray->getArrayDesc().getAttributes()[0].getDefaultValue().isZero()));
        if (nInstances > 1) {
            if (instanceId == 0) {
                for (InstanceID instance = 0; instance < nInstances; instance++) {
                    if (instance != instanceId) {
                        netMgr.send(instance, shared_ptr<SharedBuffer>(new MemoryBuffer(&isSparse, sizeof isSparse)), query);
                    }
                }
            } else {
                shared_ptr<SharedBuffer> sb = netMgr.receive(0, query);
                isSparse = *(bool*)sb->getData();
            }
        }
        if (!isSparse && nInstances > 1 && algorithm != Dense) {
            size_t minOperandSize = (iLength < jLength ? iLength : jLength)*kLength;
            size_t resultSize = iLength*jLength / nInstances;
            if (algorithm == Iterative || (algorithm == Auto && minOperandSize >= resultSize*2 && minOperandSize > incrementalMultiplyRowThreshold))
            {
                leftArray = redistribute(leftInputArray, query, psByRow);
                rightArray = redistribute(rightInputArray, query, psByCol);
                result = doIterativeMultiplyRow(0);
                return;
            }
        }
        rangePartitionBy = iLength <= jLength ? COL : ROW; // replicate smaller matrix and partition larger one
        if (iLength > iChunkLen || jLength > jChunkLen || kLength > kChunkLen) {
            if (rangePartitionBy == ROW) {
                leftArray = redistribute(leftInputArray, query, psByRow);
                rightArray = redistribute(rightInputArray, query, psReplication);
            } else {
                leftArray = redistribute(leftInputArray, query, psReplication);
                rightArray = redistribute(rightInputArray, query, psByCol);
            }
        } else {
            leftArray = leftInputArray;
            rightArray = rightInputArray;
        }
    }

    //
    // MultiplyRowSparseArrayIterator
    //
    ConstChunk const& MultiplyRowSparseArrayIterator::getChunk()
    {
        if (!hasCurrent)
            throw USER_EXCEPTION(SCIDB_SE_EXECUTION, SCIDB_LE_NO_CURRENT_ELEMENT);
        // Initialize the chunk
        if (chunk.isInitialized() && currPos == chunk.getFirstPosition(false)) {
            return chunk;
        }
        Address addr(0, currPos);
        chunk.initialize(&array, &array.desc, addr, 0);

        boost::shared_ptr<Query> query(_query.lock());
        boost::shared_ptr<ChunkIterator> outIter = chunk.getIterator(query, ChunkIterator::SPARSE_CHUNK);
        Coordinate ci = currPos[ROW];
        Coordinate cj = currPos[COL];
        Coordinate iUpper = Coordinate(ci+array.iChunkLen) < Coordinate(array.iStart + array.iLength) ? ci+array.iChunkLen : array.iStart + array.iLength;
        Coordinate jUpper = Coordinate(cj+array.jChunkLen) < Coordinate(array.jStart + array.jLength) ? cj+array.jChunkLen : array.jStart + array.jLength;
        Coordinates outCoords(2);
        TypeId attrType = array.attrType;
        Value value(TypeLibrary::getType(attrType));

        for (Coordinate i=ci; i<iUpper; i++)
        {
            for (Coordinate j=cj; j<jUpper; j++)
            {
                SparseMatrixCell* row = left[i];
                SparseMatrixCell* col = right[j];
                double partialProd = 0.0;
                if (row != NULL && col != NULL) {
                    while (true) {
                        if (col->coord > row->coord) {
                            if ((col = col->next) == NULL) {
                                break;
                            }
                        } else if (col->coord < row->coord) {
                            if ((row = row->next) == NULL) {
                                break;
                            }
                        } else {
                            partialProd += row->value * col->value;
                            if ((col = col->next) == NULL || (row = row->next) == NULL) {
                                break;
                            }
                        }
                    }
                }
                if (partialProd != 0) {
                    outCoords[ROW] = i;
                    outCoords[COL] = j;
                    bool rc = outIter->setPosition(outCoords);
                    if (!rc) assert(false);
                    DoubleToValue(attrType, partialProd, value);
                    outIter->writeItem(value);
                }
            }
        }
        outIter->flush();

        return chunk;
    }

    MultiplyRowSparseArrayIterator::MultiplyRowSparseArrayIterator(MultiplyRowArray const& array, AttributeID id)
    : MultiplyRowArrayIterator(array, id),
      left(array.dims[ROW].getStart(), array.dims[ROW].getLength()),
      right(array.dims[ROW].getStart(), array.dims[ROW].getLength())
    {
        TypeId attrType = array.attrType;
        {
            boost::shared_ptr<ConstArrayIterator> arrayIterator = array.leftArray->getConstIterator(id);
            while (!arrayIterator->end()) {
                {
                    boost::shared_ptr<ConstChunkIterator> chunkIterator = arrayIterator->getChunk().getConstIterator(ChunkIterator::IGNORE_EMPTY_CELLS|ChunkIterator::IGNORE_OVERLAPS|ChunkIterator::IGNORE_DEFAULT_VALUES);
                    while (!chunkIterator->end()) {
                        Value const& v = chunkIterator->getItem();
                        if (!v.isZero()) {                         
                            Coordinates const& pos = chunkIterator->getPosition();
                            left.add(pos[ROW], pos[COL], ValueToDouble(attrType, v));
                        }
                        ++(*chunkIterator);
                    }
                }
                ++(*arrayIterator);
            }
        }
        {
            boost::shared_ptr<ConstArrayIterator> arrayIterator = array.rightArray->getConstIterator(id);
            while (!arrayIterator->end()) {
                {
                    boost::shared_ptr<ConstChunkIterator> chunkIterator = arrayIterator->getChunk().getConstIterator(ChunkIterator::IGNORE_EMPTY_CELLS|ChunkIterator::IGNORE_OVERLAPS|ChunkIterator::IGNORE_DEFAULT_VALUES);
                    while (!chunkIterator->end()) {
                        Value const& v = chunkIterator->getItem();
                        if (!v.isZero()) {                         
                            Coordinates const& pos = chunkIterator->getPosition();
                            right.add(pos[COL], pos[ROW], ValueToDouble(attrType, v));
                        }
                        ++(*chunkIterator);
                    }
                }
                ++(*arrayIterator);
            }
        }
    }

    bool MultiplyRowArray::isSelfChunk(Coordinates const& pos) const
    {
        return InstanceID((pos[rangePartitionBy] - dims[rangePartitionBy].getStart()) / dims[rangePartitionBy].getChunkInterval()
                      / ((( dims[rangePartitionBy].getLength() + dims[rangePartitionBy].getChunkInterval() - 1) / dims[rangePartitionBy].getChunkInterval() + nInstances - 1) / nInstances)) == instanceId;
    }


} //namespace scidb
