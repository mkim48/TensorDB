/*
**
* BEGIN_COPYRIGHT
*
* This file is part of TensorDB.
* Copyright (C) 2014 Mijung Kim and K. Selcuk Candan. All right Reserved.
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
 * PhysicalTwoNormCoreArray.cpp
 *
 * @author: Mijung Kim (mijung.kim.1@asu.edu) 
 */
#include <float.h>
#include <string.h>
#include "query/Operator.h"
#include "array/Metadata.h"
#include "array/MemArray.h"
#include "array/DelegateArray.h"

namespace scidb
{

class PhysicalTwoNormCoreArray : public  PhysicalOperator
{
public:
    PhysicalTwoNormCoreArray(std::string const& logicalName,
                     std::string const& physicalName,
                     Parameters const& parameters,
                     ArrayDesc const& schema):
	    PhysicalOperator(logicalName, physicalName, parameters, schema)
	{
	}

    virtual bool changesDistribution(std::vector<ArrayDesc> const&) const
    {
        return true;
    }

    virtual ArrayDistribution getOutputDistribution(
            std::vector<ArrayDistribution> const&,
            std::vector<ArrayDesc> const&) const
    {
        return ArrayDistribution(psLocalInstance);
    }


	/***
	 */
	boost::shared_ptr< Array> execute(std::vector< boost::shared_ptr< Array> >& inputArrays,
            boost::shared_ptr<Query> query)
	{
		boost::shared_ptr<Array> array = inputArrays[0];
		boost::shared_ptr<Array> coreArray = inputArrays[1];
        if ( query->getInstancesCount() > 1) {
           uint64_t coordinatorID = (int64_t)query->getCoordinatorID() == -1 ? query->getInstanceID() : query->getCoordinatorID();
            array = redistribute(array, query, psLocalInstance, "", coordinatorID);
            if ( query->getInstanceID() != coordinatorID) {
                return boost::shared_ptr<Array>(new MemArray(_schema));
            }
        }

		AttributeID attr = 0;
		boost::shared_ptr<ConstArrayIterator> iter = array->getConstIterator(attr);

        Dimensions const& dims = array->getArrayDesc().getDimensions();
        Coordinate iLen = dims[0].getLength();
        Coordinate iChunkLen = dims[0].getChunkInterval();
        Coordinate jChunkLen = dims[1].getChunkInterval();
        Coordinate iEnd = dims[0].getEndMax();
        Coordinate jEnd = dims[1].getEndMax();

		std::vector<double> aggvals(iLen);
		for(Coordinate i=0;i<iLen;i++) {
			aggvals[i] = 0.0;
		}
		while (!iter->end()) {
			Coordinates const& start = iter->getPosition();
			ConstChunk const* chunk = &iter->getChunk();
		 	PinBuffer s1(*chunk);
			void* dat = (void*)chunk->getData();

			if (chunk->isRLE()) {

                MemChunk leftMatChunk;
                ConstRLEPayload leftPayload((char const*)dat);
                if (leftPayload.nSegments() != 1 || leftPayload.getSegment(0).same || leftPayload.getSegment(0).length() != iChunkLen*jChunkLen) {
                    if (leftPayload.nSegments() == 1 && leftPayload.getSegment(0).same) {
                        Value val;
                        leftPayload.getValueByIndex(val, leftPayload.getSegment(0).valueIndex);
                    }
                    MaterializedArray::materialize(leftMatChunk, *chunk, MaterializedArray::DenseFormat);
                    dat = leftMatChunk.getData();
                } else {
                    dat = leftPayload.getRawValue(leftPayload.getSegment(0).valueIndex);
                }
            }
			double* val = (double*) dat;
			Coordinates last(2);
            last[0] = Coordinate(start[0]+iChunkLen-1) < Coordinate(iEnd) ? start[0]+iChunkLen-1 : iEnd;
            last[1] = Coordinate(start[1]+jChunkLen-1) < Coordinate(jEnd) ? start[1]+jChunkLen-1 : jEnd;
			
			for(Coordinate i=start[0];i<=last[0];i++) {
				for(Coordinate j=start[1];j<=last[1];j++) {
					aggvals[i-start[0]] += pow(val[(i-start[0])*jChunkLen+(j-start[1])],2); //2-norm
				}
			}	
			++(*iter);
        }
		for(Coordinate i=0;i<iLen;i++) {
			aggvals[i] = sqrt(aggvals[i]);
		}
        boost::shared_ptr<ConstArrayIterator> coreIter = coreArray->getConstIterator(0);
        Coordinates pos(2);
        pos[0]=1;
        pos[1]=1;
        Dimensions const& coredims = coreArray->getArrayDesc().getDimensions();
        Coordinate jCoreChunkLen = coredims[1].getChunkInterval();
        double* coreData;
        if(!coreIter->setPosition(pos)){
            boost::shared_ptr<ArrayIterator> coreArrayIter = coreArray->getIterator(attr);
            Chunk* coreChunk = &coreArrayIter->newChunk(pos);
            PinBuffer s2(*coreChunk);
            coreChunk->setRLE(false);
            coreChunk->allocate(jCoreChunkLen*sizeof(double));
            coreData = (double*)coreChunk->getData();
        } else {
            ConstChunk const* coreChunk = &coreIter->getChunk();
            PinBuffer s2(*coreChunk);
            coreData = (double*)coreChunk->getData();
        }
        for(Coordinate i=0;i<iLen;i++) {
            coreData[i] = aggvals[i];
        }

		iter->reset();
        while (!iter->end()) {
            Coordinates const& start = iter->getPosition();
            ConstChunk const* chunk = &iter->getChunk();
            PinBuffer s1(*chunk);
            void* dat = (void*)chunk->getData();

            if (chunk->isRLE()) {
                MemChunk leftMatChunk;
                ConstRLEPayload leftPayload((char const*)dat);
                if (leftPayload.nSegments() != 1 || leftPayload.getSegment(0).same || leftPayload.getSegment(0).length() != iChunkLen*jChunkLen) {
                    if (leftPayload.nSegments() == 1 && leftPayload.getSegment(0).same) {
                        Value val;
                        leftPayload.getValueByIndex(val, leftPayload.getSegment(0).valueIndex);
                    }
                    MaterializedArray::materialize(leftMatChunk, *chunk, MaterializedArray::DenseFormat);
                    dat = leftMatChunk.getData();
                } else {
                    dat = leftPayload.getRawValue(leftPayload.getSegment(0).valueIndex);
                }
            }
            double* val = (double*) dat;
            Coordinates last(2);
            last[0] = Coordinate(start[0]+iChunkLen-1) < Coordinate(iEnd) ? start[0]+iChunkLen-1 : iEnd;
            last[1] = Coordinate(start[1]+jChunkLen-1) < Coordinate(jEnd) ? start[1]+jChunkLen-1 : jEnd;

            for(Coordinate i=start[0];i<=last[0];i++) {
                for(Coordinate j=start[1];j<=last[1];j++) {
					val[(i-start[0])*jChunkLen+(j-start[1])] =  val[(i-start[0])*jChunkLen+(j-start[1])]/aggvals[i-start[0]];
                }
            }      
			++(*iter);
        }
		return array;
	}
};

DECLARE_PHYSICAL_OPERATOR_FACTORY(PhysicalTwoNormCoreArray, "twoNormCoreArray", "PhysicalTwoNormCoreArray")

} //namespace scidb
