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
 * MultiplyRowArray.h
 *
 * @author: Mijung Kim (mijung.kim.1@asu.edu) 
 */
#include "query/Operator.h"
#include "array/Metadata.h"
#include "array/MemArray.h"
#include "util/JobQueue.h"
#include "util/ThreadPool.h"
#include "SparseMatrix.h"

#ifndef MULTIPLY_ROW_ARRAY_H
#define MULTIPLY_ROW_ARRAY_H

namespace scidb
{

class MultiplyRowArray;

class MultiplyRowArrayIterator : public ConstArrayIterator
{
  public:
        virtual ConstChunk const& getChunk();
        virtual bool end();
        virtual void operator ++();
        virtual Coordinates const& getPosition();
        virtual bool setPosition(Coordinates const& pos);
        virtual void reset();

        MultiplyRowArrayIterator(MultiplyRowArray const& array, AttributeID id);

  protected:
    bool isSelfChunk() const;

    MultiplyRowArray const& array;
    AttributeID attr;
    Coordinates currPos;
    bool hasCurrent;
    MemChunk chunk;
    boost::weak_ptr<Query> _query;
};

class MultiplyRowSparseArrayIterator : public MultiplyRowArrayIterator
{
  public:
        virtual ConstChunk const& getChunk();
        MultiplyRowSparseArrayIterator(MultiplyRowArray const& array, AttributeID id);

  private:
    SparseRowMatrix left;
    SparseColumnMatrix right;
};

class MultiplyRowArray : public Array
{
    friend class MultiplyRowArrayIterator;
    friend class MultiplyRowSparseArrayIterator;

    class MultiplyRowJob : public Job
    {
        MultiplyRowArray const& array;
        void const* leftData;
        void const* rightData;
        void* resultData;
        Coordinate ci, cj, ck, offset;
      public:
        MultiplyRowJob(MultiplyRowArray const& array,
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
    MultiplyRowArray(ArrayDesc desc,
                  boost::shared_ptr<Array> const& leftArray,
                  boost::shared_ptr<Array> const& rightArray,
                  const boost::shared_ptr<Query>& query,
                  Algorithm algorithm
        );
    
  protected:
    boost::shared_ptr<Array> doIterativeMultiplyRow(AttributeID attr);
    void multiplyChunks(boost::shared_ptr<ConstArrayIterator> const& rowArrayIter,
                        boost::shared_ptr<ConstArrayIterator> const& colArrayIter,
                        boost::shared_ptr<ChunkIterator> const& outIter, Coordinate ci, Coordinate cj, Coordinate ck) const;
    
    void multiplyChunks(void const* leftData,
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

#endif /* MUTLIPLY_ROW_ARRAY_H */
