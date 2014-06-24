/*
**
* BEGIN_COPYRIGHT
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
 * MultiplyArrayComp.h
 *
 *  Created on: Mar 9, 2010
 */
#include "query/Operator.h"
#include "array/Metadata.h"
#include "array/MemArray.h"
#include "util/JobQueue.h"
#include "util/ThreadPool.h"
#include "SparseMatrix.h"

#ifndef MULTIPLY_COMP_ARRAY_H
#define MULTIPLY_COMP_ARRAY_H

namespace scidb
{

class MultiplyArrayComp;

class MultiplyArrayCompIterator : public ConstArrayIterator
{
  public:
        virtual ConstChunk const& getChunk();
        virtual bool end();
        virtual void operator ++();
        virtual Coordinates const& getPosition();
        virtual bool setPosition(Coordinates const& pos);
        virtual void reset();

        MultiplyArrayCompIterator(MultiplyArrayComp const& array, AttributeID id);

  protected:
    bool isSelfChunk() const;

    MultiplyArrayComp const& array;
    AttributeID attr;
    Coordinates currPos;
    bool hasCurrent;
    MemChunk chunk;
    boost::weak_ptr<Query> _query;
};

class MultiplySparseArray1CompIterator : public MultiplyArrayCompIterator
{
  public:
        virtual ConstChunk const& getChunk();
        MultiplySparseArray1CompIterator(MultiplyArrayComp const& array, AttributeID id);

  private:
    SparseRowMatrix left;
    SparseColumnMatrix right;
};

class MultiplyArrayComp : public Array
{
    friend class MultiplyArrayCompIterator;
    friend class MultiplySparseArray1CompIterator;

    class MultiplyJob : public Job
    {
        MultiplyArrayComp const& array;
        void const* leftData;
        void const* rightData;
        void* resultData;
        Coordinate ci, cj, ck, offset;
      public:
        MultiplyJob(MultiplyArrayComp const& array,
                    void const* leftData,
                    void const* rightData,
                    void* resultData,
                    Coordinate ci, Coordinate cj, Coordinate ck, Coordinate offset,
                    const boost::shared_ptr<Query>& query);
        virtual void run();
    };

  public:
    enum Algorithm { 
        Auto,
        Sparse,                                 
        Dense,                                 
        Iterative
    };

    virtual ArrayDesc const& getArrayDesc() const;
    virtual boost::shared_ptr<ConstArrayIterator> getConstIterator(AttributeID attr) const;
    MultiplyArrayComp(ArrayDesc desc,
                  boost::shared_ptr<Array> const& leftArray,
                  boost::shared_ptr<Array> const& rightArray,
                  const boost::shared_ptr<Query>& query,
                  Algorithm algorithm
        );
    
  protected:
    boost::shared_ptr<Array> doIterativeMultiply(AttributeID attr);
    void multiplyChunks(boost::shared_ptr<ConstArrayIterator> const& rowArrayIter,
                        boost::shared_ptr<ConstArrayIterator> const& colArrayIter,
                        boost::shared_ptr<ChunkIterator> const& outIter, Coordinate ci, Coordinate cj, Coordinate ck) const;
    
    void multiplyChunks(void const* leftData,
                        void const* rightData,
                        void* resultData,
                        Coordinate ci, Coordinate cj, Coordinate ck, Coordinate offset) const;
    void compressedMultiplyChunks(void const* leftData,
                        void const* rightData,
                        void* resultData,
                        Coordinate ci, Coordinate cj, Coordinate ck, Coordinate offset) const;
    
    bool isSelfChunk(Coordinates const& pos) const;

    ArrayDesc desc;
    Dimensions const& dims;
    boost::shared_ptr<Array> leftArray;
    boost::shared_ptr<Array> rightArray;
    size_t iChunkLen, jChunkLen, kChunkLen;
    uint64_t iLength, jLength, kLength;
    Coordinate iStart, jStart, kStartLeft, kStartRight;
    bool isSparse;
    InstanceID instanceId;
    size_t nInstances;
    boost::weak_ptr<Query> _query;
    boost::shared_ptr<Array> result;
    int rangePartitionBy;
    TypeId attrType;
    TypeId multiplyResultType;
    FunctionPointer multiply;
    FunctionPointer add;
    FunctionPointer convert;
    boost::shared_ptr<JobQueue> queue;
    boost::shared_ptr<ThreadPool> threadPool;
    long int nCPUs;
};


} //namespace scidb

#endif /* MUTLIPLY_ARRAY_H */
