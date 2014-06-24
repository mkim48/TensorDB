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
 * LogicalTwoNormCoreArray.cpp
 *
 * @author: Mijung Kim (mijung.kim.1@asu.edu) 
 */

#include "query/Operator.h"


namespace scidb
{

inline bool hasSingleAttribute(ArrayDesc const& desc)
{
    return desc.getAttributes(true).size() == 1;
}

/**
 * @brief The operator: TwoNormCoreArray().
 *
 * @par Synopsis:
 *   TwoNormCoreArray( array, core )
 *
 * @par Summary:
 *   Normalize array by 2-norm and the weight is copied into core
 *
 */
class LogicalTwoNormCoreArray : public  LogicalOperator
{
public:
	LogicalTwoNormCoreArray(const std::string& logicalName, const std::string& alias):
	    LogicalOperator(logicalName, alias)
	{
		ADD_PARAM_INPUT()
		ADD_PARAM_INPUT()
	}

    ArrayDesc inferSchema(std::vector< ArrayDesc> schemas, boost::shared_ptr< Query> query)
    {
        assert(schemas.size() == 2);

        return schemas[0]; 
	}

};

DECLARE_LOGICAL_OPERATOR_FACTORY(LogicalTwoNormCoreArray, "twoNormCoreArray")

} //namespace scidb
