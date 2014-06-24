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
 * PhysicalNormResidual.cpp
 *
 * @author: Mijung Kim (mijung.kim.1@asu.edu) 
 */
#include <math.h>
#include "query/Operator.h"
#include "array/TupleArray.h"
#include "array/Metadata.h"
#include "array/MemArray.h"
#include "array/DelegateArray.h"

namespace scidb
{

class PhysicalNormResidual : public  PhysicalOperator
{
public:
    PhysicalNormResidual(std::string const& logicalName,
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

		boost::shared_ptr<Array> leftArray = inputArrays[0];
        boost::shared_ptr<Array> rightArray = inputArrays[1];	

        AttributeID attr = 0;
        boost::shared_ptr<ConstArrayIterator> leftIter = leftArray->getConstIterator(attr);
        boost::shared_ptr<ConstArrayIterator> rightIter = rightArray->getConstIterator(attr);
	
		ArrayDesc const& rightDesc = rightArray->getArrayDesc();
        Dimensions const& rightDims = rightDesc.getDimensions();
		size_t iRightChunkLen = rightDims[0].getChunkInterval();
		size_t jRightChunkLen = rightDims[1].getChunkInterval();
		ArrayDesc const& leftDesc = leftArray->getArrayDesc();
        Dimensions const& leftDims = leftDesc.getDimensions();
		size_t iLeftChunkLen = leftDims[0].getChunkInterval();
		size_t jLeftChunkLen = leftDims[1].getChunkInterval();
		double res = 0;
		MemChunk leftMatChunk;
		while (!leftIter->end()) {
			ConstChunk const* leftChunk = &leftIter->getChunk();
        	PinBuffer s1(*leftChunk);
			double* leftData = (double*)leftChunk->getData();

			if (leftChunk->isRLE()) {
                ConstRLEPayload leftPayload((char const*)leftData);
                if (leftPayload.nSegments() != 1 || leftPayload.getSegment(0).same || leftPayload.getSegment(0).length() != iLeftChunkLen*jLeftChunkLen) {
                    if (leftPayload.nSegments() == 1 && leftPayload.getSegment(0).same) {
                        Value val;
                        leftPayload.getValueByIndex(val, leftPayload.getSegment(0).valueIndex);
                    }
                    MaterializedArray::materialize(leftMatChunk, *leftChunk, MaterializedArray::DenseFormat);
                    leftData = (double*)leftMatChunk.getData();
                } else {
                    leftData = (double*)leftPayload.getRawValue(leftPayload.getSegment(0).valueIndex);
                }
        	}
			ConstChunk const* rightChunk = &rightIter->getChunk();
            PinBuffer s2(*rightChunk);

			double* rightData = (double*)rightChunk->getData();
			for(uint64_t i=0;i<iLeftChunkLen*jLeftChunkLen;i++) {
				res += pow((leftData[i] - rightData[i]),2); 
			}
			++(*leftIter);
			++(*rightIter);
		}
		vector< boost::shared_ptr<Tuple> > tuples(1);
        Tuple& tuple = *new Tuple(1);
        tuple[0].setDouble(sqrt(res));
		tuples[0] = boost::shared_ptr<Tuple>(&tuple);
        return boost::shared_ptr<Array>(new TupleArray(_schema, tuples, Coordinate(query->getInstanceID())));
	}
};

DECLARE_PHYSICAL_OPERATOR_FACTORY(PhysicalNormResidual, "norm_residual", "PhysicalNormResidual")

} //namespace scidb
