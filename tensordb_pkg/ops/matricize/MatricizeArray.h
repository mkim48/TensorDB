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

#include "query/Operator.h"

#ifndef MATRICIZE_ARRAY_H
#define MATRICIZE_ARRAY_H

namespace scidb
{

class MatricizeArray;

class MatricizeArrayIterator : public ConstArrayIterator
{
  public:
        virtual ConstChunk const& getChunk();
        virtual bool end();
        virtual void operator ++();
        virtual Coordinates const& getPosition();
        virtual bool setPosition(Coordinates const& pos);
        virtual void reset();

        MatricizeArrayIterator(MatricizeArray const& array, AttributeID id);

  protected:
    MatricizeArray const& array;
    AttributeID attr;
    Coordinates currPos;
    bool hasCurrent;
    MemChunk chunk;
    boost::weak_ptr<Query> _query;
};

class MatricizeArray : public Array
{
    friend class MatricizeArrayIterator;

public:
    virtual ArrayDesc const& getArrayDesc() const;
    virtual boost::shared_ptr<ConstArrayIterator> getConstIterator(AttributeID attr) const;
    MatricizeArray(ArrayDesc desc,
                  boost::shared_ptr<Array> const& inputArray,
                  boost::shared_ptr<Query> const& query,
                  Coordinate rmode, 
                  Coordinate cmode 
        );
    
  protected:
    void matricizeChunks(boost::shared_ptr<ConstArrayIterator> const& inputArrayIter, boost::shared_ptr<ChunkIterator> const& outIter, Coordinates currPos, Coordinates inputPos) const;
    
    ArrayDesc desc;
    Dimensions const& dims;
    boost::shared_ptr<Array> inputArray;
    size_t iChunkLen, jChunkLen, totalInputChunkLen, restInputChunkLen, restInputIncChunkLen;
    uint64_t iLength, jLength;
    Coordinate iStart, jStart;
    Coordinates inputIncChunkLen;
    Coordinates inputChunkLen;
    boost::weak_ptr<Query> _query;
    boost::shared_ptr<Array> result;
    Coordinate rowmode;
    Coordinate colmode;
	Coordinates nChunk;
	size_t restNumChunk;
};


} //namespace scidb

#endif /* MATRICIZE_ARRAY_H */
