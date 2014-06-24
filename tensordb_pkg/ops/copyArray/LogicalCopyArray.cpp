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

#include <boost/foreach.hpp>
#include <map>

#include "query/Operator.h"
#include "system/SystemCatalog.h"
#include "system/Exceptions.h"
#include <smgr/io/Storage.h>

/*
 * LogicalCopyArray.cpp
 *
 * @author: Mijung Kim (mijung.kim.1@asu.edu)
 */


namespace scidb
{

inline bool hasSingleAttribute(ArrayDesc const& desc)
{
    return desc.getAttributes(true).size() == 1;
}

/**
 * @brief The operator: copyArray().
 *
 * @par Synopsis:
 *		copyArray(operator(args), array)
 *
 * @par Summary:
 *   This operator copies the result of 'operator(args)' to a temporary array 'array' and is used for updating the intermediate results (e.g., in the in-database CP decomposition, for updating the factor matrices, which get updated in each iteration).  
 *
 */
class LogicalCopyArray : public  LogicalOperator
{
public:
	LogicalCopyArray(const std::string& logicalName, const std::string& alias):
	    LogicalOperator(logicalName, alias)
	{
		ADD_PARAM_INPUT()
		ADD_PARAM_INPUT()
        ADD_PARAM_VARIES()
	}

    std::vector<boost::shared_ptr<OperatorParamPlaceholder> > nextVaryParamPlaceholder(const std::vector<ArrayDesc> &schemas)
    {
        std::vector<boost::shared_ptr<OperatorParamPlaceholder> > res;
        res.push_back(PARAM_CONSTANT(TID_STRING));
        res.push_back(END_OF_VARIES_PARAMS());
        return res;
    }

    ArrayDesc inferSchema(std::vector< ArrayDesc> schemas, boost::shared_ptr< Query> query)
    {
		assert(schemas.size() == 2);

		return schemas[1];
	}

};

DECLARE_LOGICAL_OPERATOR_FACTORY(LogicalCopyArray, "copyArray")

} //namespace scidb
