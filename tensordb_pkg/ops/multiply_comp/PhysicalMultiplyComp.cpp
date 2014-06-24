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
 * PhysicalMultiplyComp.cpp
 *
 *  Created on: Mar 9, 2010
 */
#include "query/Operator.h"
#include "array/Metadata.h"
#include "array/MemArray.h"
#include "query/ops/multiply_comp/MultiplyArrayComp.h"

namespace scidb
{

class PhysicalMultiplyComp : public  PhysicalOperator
{
public:
    PhysicalMultiplyComp(std::string const& logicalName,
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
        DimensionDesc const& d2 = inputSchemas[1].getDimensions()[0];
        PartitioningSchema ps = d1.getLength() <= d2.getLength() ? psByCol : psByRow;
        return ArrayDistribution(ps);
    }

	/***
	 */
	boost::shared_ptr< Array> execute(std::vector< boost::shared_ptr< Array> >& inputArrays,
            boost::shared_ptr<Query> query)
	{
        MultiplyArrayComp::Algorithm algorithm = MultiplyArrayComp::Auto;
        if (_parameters.size() > 0) { 
            string a = ((boost::shared_ptr<OperatorParamPhysicalExpression>&)_parameters[0])->getExpression()->evaluate().getString();
            if (a == "iterative") { 
                algorithm = MultiplyArrayComp::Iterative;
            } else if (a == "dense") { 
                algorithm = MultiplyArrayComp::Dense;
            } else if (a == "sparse") { 
                algorithm = MultiplyArrayComp::Sparse;
            } else if (a == "auto") { 
                algorithm = MultiplyArrayComp::Auto;
            }
        }
        return boost::shared_ptr<Array>(new MultiplyArrayComp(_schema, inputArrays[0], inputArrays[1], query, algorithm));
	}
};

DECLARE_PHYSICAL_OPERATOR_FACTORY(PhysicalMultiplyComp, "multiply_comp", "PhysicalMultiplyComp")

} //namespace scidb
