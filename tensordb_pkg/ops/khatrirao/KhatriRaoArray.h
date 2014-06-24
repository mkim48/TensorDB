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
 * KhatriRaoArray.h
 *
 * @author: Mijung Kim (mijung.kim.1@asu.edu) 
 */
#include "query/Operator.h"

#ifndef KHATRIRAO_ARRAY_H
#define KHATRIRAO_ARRAY_H

namespace scidb
{

class KhatriRaoArray;

class KhatriRaoArrayIterator : public ConstArrayIterator
{
  public:
        virtual ConstChunk const& getChunk();
        virtual bool end();
        virtual void operator ++();
        virtual Coordinates const& getPosition();
        virtual bool setPosition(Coordinates const& pos);
        virtual void reset();

        KhatriRaoArrayIterator(KhatriRaoArray const& array, AttributeID id);

  protected:

    KhatriRaoArray const& array;
    AttributeID attr;
    Coordinates currPos;
    bool hasCurrent;
    MemChunk chunk;
    boost::weak_ptr<Query> _query;
};

class KhatriRaoArray : public Array
{
    friend class KhatriRaoArrayIterator;

    class KhatriRaoJob : public Job
    {
        KhatriRaoArray const& array;
        void const* leftData;
        void const* rightData;
        void* resultData;
        Coordinate ci, cj, ck, offset;
      public:
        KhatriRaoJob(KhatriRaoArray const& array,
                    void const* leftData,
                    void const* rightData,
                    void* resultData,
                    Coordinate ci, Coordinate cj, Coordinate ck, Coordinate offset,
                    const boost::shared_ptr<Query>& query);
        virtual void run();
    };

  public:
    virtual ArrayDesc const& getArrayDesc() const;
    virtual boost::shared_ptr<ConstArrayIterator> getConstIterator(AttributeID attr) const;
    KhatriRaoArray(ArrayDesc desc,
                  boost::shared_ptr<Array> const& leftArray,
                  boost::shared_ptr<Array> const& rightArray,
                  const boost::shared_ptr<Query>& query
        );
    
  protected:
    boost::shared_ptr<Array> doIterativeKhatriRao(AttributeID attr);
    void khatriraoChunks(boost::shared_ptr<ConstArrayIterator> const& rowArrayIter,
                        boost::shared_ptr<ConstArrayIterator> const& colArrayIter,
                        boost::shared_ptr<ChunkIterator> const& outIter, Coordinate ci, Coordinate cj, Coordinate ck) const;
    
    void khatriraoChunks(void const* leftData,
                        void const* rightData,
                        void* resultData,
                        Coordinate ci, Coordinate cj, Coordinate ck, Coordinate offset) const;
    

    ArrayDesc desc;
    Dimensions const& dims;
    boost::shared_ptr<Array> leftArray;
    boost::shared_ptr<Array> rightArray;
    size_t iChunkLen, jChunkLen, iLeftChunkLen, jLeftChunkLen,iRightChunkLen,jRightChunkLen, jRightNumChunk;
    uint64_t iLength, jLength, iLeftLength, jLeftLength,iRightLength, jRightLength;
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

#endif /* KHATRIRAO_ARRAY_H */
