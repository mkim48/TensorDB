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
 * PhysicalCopyArray.cpp
 *
 * @author: Mijung Kim (mijung.kim.1@asu.edu) 
 */
#include <string.h>
#include "query/Operator.h"
#include "array/Metadata.h"
#include "array/MemArray.h"
#include "array/DelegateArray.h"

namespace scidb
{

class PhysicalCopyArray : public  PhysicalOperator
{
public:
    PhysicalCopyArray(std::string const& logicalName,
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
            std::vector<ArrayDistribution> const& inputDistributions,
            std::vector<ArrayDesc> const& inputSchemas) const
    {
        DimensionDesc const& d1 = inputSchemas[0].getDimensions()[0];
        DimensionDesc const& d2 = inputSchemas[1].getDimensions()[1];
        PartitioningSchema ps = d1.getLength() <= d2.getLength() ? psByCol : psByRow;
        return ArrayDistribution(ps);
    }

	/***
	 */
	boost::shared_ptr< Array> execute(std::vector< boost::shared_ptr< Array> >& inputArrays,
            boost::shared_ptr<Query> query)
	{

		boost::shared_ptr<Array> srcArray = inputArrays[0];
		boost::shared_ptr<Array> dstArray = inputArrays[1];

		ArrayDesc arrayDesc = srcArray->getArrayDesc();
		for (AttributeID attr=0; attr<arrayDesc.getAttributes().size(); attr++)
    	{
			boost::shared_ptr<ConstArrayIterator> srcIter = srcArray->getConstIterator(attr);
			boost::shared_ptr<ConstArrayIterator> dstIter = dstArray->getConstIterator(attr);

			while (!srcIter->end()) {

				Coordinates const& srcPos = srcIter->getPosition();
            	Dimensions const& dims = srcArray->getArrayDesc().getDimensions();
            	size_t chunkLen = 1;
            	Coordinates last(dims.size());
            	for(uint32_t i=0;i<dims.size();i++) {
               		Coordinate ci = srcPos[i];
               		Coordinate iChunkLen = dims[i].getChunkInterval();
					chunkLen *= iChunkLen;
            	}

				ConstChunk const* leftChunk = &srcIter->getChunk();
		 		PinBuffer s1(*leftChunk);
        		void* leftData = leftChunk->getData();
				MemChunk inputMatChunk;
        		if (leftChunk->isRLE()) {
           			ConstRLEPayload inputPayload((char const*)leftData);
           			if (inputPayload.nSegments() != 1 || inputPayload.getSegment(0).same || inputPayload.getSegment(0).length() != chunkLen) {
              			if (inputPayload.nSegments() == 1 && inputPayload.getSegment(0).same) {
                   			Value val;
                   			inputPayload.getValueByIndex(val, inputPayload.getSegment(0).valueIndex);
                   			if (val.isZero()) { // result is 0
                        		return dstArray;
                    		}
                		}
                		MaterializedArray::materialize(inputMatChunk, *leftChunk, MaterializedArray::DenseFormat);
                		leftData = inputMatChunk.getData();
            		} else {
                		leftData = inputPayload.getRawValue(inputPayload.getSegment(0).valueIndex);
            		}
        		}

            	size_t chunkSize = chunkLen*sizeof(double*);
				double* rightData;

				if(!dstIter->setPosition(srcPos)){ 
					boost::shared_ptr<ArrayIterator> dstArrayIter = dstArray->getIterator(attr);
					Chunk* rightChunk = &dstArrayIter->newChunk(srcPos);
					PinBuffer s2(*rightChunk);
					rightChunk->setRLE(false);
					rightChunk->allocate(chunkSize);
					rightData = (double*)rightChunk->getData();
				} else {
					ConstChunk const* rightChunk = &dstIter->getChunk();
        			if (rightChunk->isRLE()) {
						throw USER_EXCEPTION(SCIDB_SE_EXECUTION, SCIDB_LE_OP_COPYARRAY_ERROR1); 
					}
					PinBuffer s2(*rightChunk);
					rightData = (double*)rightChunk->getData();
				}
			
				double* left = (double*)leftData;
				memcpy(rightData, left, chunkSize); 

				++(*srcIter);
			}
		}
		return dstArray;
	}
};

DECLARE_PHYSICAL_OPERATOR_FACTORY(PhysicalCopyArray, "copyArray", "PhysicalCopyArray")

} //namespace scidb
