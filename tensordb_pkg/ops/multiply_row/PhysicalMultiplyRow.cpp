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
 * PhysicalMultiplyRow.cpp
 *
 * @author: Mijung Kim (mijung.kim.1@asu.edu) 
 */
#include "query/Operator.h"
#include "array/Metadata.h"
#include "array/MemArray.h"
#include "query/ops/multiply_row/MultiplyRowArray.h"

namespace scidb
{

class PhysicalMultiplyRow : public  PhysicalOperator
{
public:
    PhysicalMultiplyRow(std::string const& logicalName,
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
        MultiplyRowArray::Algorithm algorithm = MultiplyRowArray::Auto;
        if (_parameters.size() > 0) { 
            string a = ((boost::shared_ptr<OperatorParamPhysicalExpression>&)_parameters[0])->getExpression()->evaluate().getString();
            if (a == "iterative") { 
                algorithm = MultiplyRowArray::Iterative;
            } else if (a == "dense") { 
                algorithm = MultiplyRowArray::Dense;
            } else if (a == "sparse") { 
                algorithm = MultiplyRowArray::Sparse;
            } else if (a == "auto") { 
                algorithm = MultiplyRowArray::Auto;
            }
        }
        return boost::shared_ptr<Array>(new MultiplyRowArray(_schema, inputArrays[0], inputArrays[1], query, algorithm));
	}
};

DECLARE_PHYSICAL_OPERATOR_FACTORY(PhysicalMultiplyRow, "multiply_row", "PhysicalMultiplyRow")

} //namespace scidb
