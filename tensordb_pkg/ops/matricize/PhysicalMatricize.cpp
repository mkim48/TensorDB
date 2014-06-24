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
 * PhysicalMatricize.cpp
 *
 * @author: Mijung Kim (mijung.kim.1@asu.edu)
 */

#include "query/Operator.h"
#include "array/Metadata.h"
#include "array/MemArray.h"
#include "query/ops/matricize/MatricizeArray.h"

namespace scidb
{

class PhysicalMatricize : public  PhysicalOperator
{
public:
    PhysicalMatricize(std::string const& logicalName,
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
        return ArrayDistribution(psUndefined);
    }

	/***
	 */
	boost::shared_ptr< Array> execute(std::vector< boost::shared_ptr< Array> >& inputArrays, boost::shared_ptr<Query> query)
	{
		boost::shared_ptr<Array> inputArray = inputArrays[0];
		Dimensions const& inputDims = inputArray->getArrayDesc().getDimensions();
        Coordinate ndims = inputDims.size();
		Coordinate rowmode = 0;
		Coordinate colmode = 1;
        if (_parameters.size() == 1) { 
            rowmode = ((boost::shared_ptr<OperatorParamPhysicalExpression>&)_parameters[0])->getExpression()->evaluate().getInt64();
			if(rowmode < ndims)
                colmode = ndims;
            else
                colmode = ndims-1;

            rowmode -= 1;
            colmode -= 1;
        }
        return boost::shared_ptr<Array>(new MatricizeArray(_schema, inputArray, query,  rowmode, colmode));
	}
};

DECLARE_PHYSICAL_OPERATOR_FACTORY(PhysicalMatricize, "matricize", "PhysicalMatricize")

} //namespace scidb
