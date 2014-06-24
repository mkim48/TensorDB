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
 * HadamardArray.h
 *
 * @author: Mijung Kim (mijung.kim.1@asu.edu) 
 */
#include "query/Operator.h"
#include "array/Metadata.h"
#include "array/MemArray.h"
#include "util/JobQueue.h"
#include "util/ThreadPool.h"

#ifndef ELEMPROD_ARRAY_H
#define ELEMPROD_ARRAY_H

namespace scidb
{

class HadamardArray;

class HadamardArrayIterator : public ConstArrayIterator
{
  public:
        virtual ConstChunk const& getChunk();
        virtual bool end();
        virtual void operator ++();
        virtual Coordinates const& getPosition();
        virtual bool setPosition(Coordinates const& pos);
        virtual void reset();

        HadamardArrayIterator(HadamardArray const& array, AttributeID id);

  protected:

    HadamardArray const& array;
    AttributeID attr;
    Coordinates currPos;
    bool hasCurrent;
    MemChunk chunk;
    boost::weak_ptr<Query> _query;
};

class HadamardArray : public Array
{
    friend class HadamardArrayIterator;

    class HadamardJob : public Job
    {
        HadamardArray const& array;
        void const* leftData;
        void const* rightData;
        void* resultData;
        Coordinate ci, cj, offset;
      public:
        HadamardJob(HadamardArray const& array,
                    void const* leftData,
                    void const* rightData,
                    void* resultData,
                    Coordinate ci, Coordinate cj, Coordinate offset,
                    const boost::shared_ptr<Query>& query);
        virtual void run();
    };

  public:

    virtual ArrayDesc const& getArrayDesc() const;
    virtual boost::shared_ptr<ConstArrayIterator> getConstIterator(AttributeID attr) const;
    HadamardArray(ArrayDesc desc,
                  boost::shared_ptr<Array> const& leftArray,
                  boost::shared_ptr<Array> const& rightArray,
                  const boost::shared_ptr<Query>& query
        );
    
  protected:
    boost::shared_ptr<Array> doIterativeKrProd(AttributeID attr);
    void hadamardChunks(boost::shared_ptr<ConstArrayIterator> const& rowArrayIter, boost::shared_ptr<ConstArrayIterator> const& colArrayIter, boost::shared_ptr<ChunkIterator> const& outIter, Coordinate ci, Coordinate cj) const;
    
    void hadamardChunks(void const* leftData,
                        void const* rightData,
                        void* resultData,
                        Coordinate ci, Coordinate cj, Coordinate offset) const;
    

    ArrayDesc desc;
    Dimensions const& dims;
    boost::shared_ptr<Array> leftArray;
    boost::shared_ptr<Array> rightArray;
    size_t iChunkLen, jChunkLen, iLeftChunkLen, iRightChunkLen;
    uint64_t iLength, jLength, iLeftLength, iRightLength;
    Coordinate iStart, jStart, iLeftStart, iRightStart, jLeftStart, jRightStart;
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

#endif /* ELEMPROD_ARRAY_H */
