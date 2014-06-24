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
 * LogicalEigen.cpp
 *
 * @author: Mijung Kim (mijung.kim.1@asu.edu)
 */

#include "query/Operator.h"


namespace scidb
{

/**
 * @brief The operator: eigen().
 *
 * @par Synopsis:
 *   eigen( srcArray )
 *
 * @par Summary:
 *   Produces the matrix eigen decomposition of a square matrix.
 *
 */
class LogicalEigen : public  LogicalOperator
{
public:
	LogicalEigen(const std::string& logicalName, const std::string& alias):
	    LogicalOperator(logicalName, alias)
	{
		ADD_PARAM_INPUT()
		ADD_PARAM_VARIES()
	}

    std::vector<boost::shared_ptr<OperatorParamPlaceholder> > nextVaryParamPlaceholder(const std::vector<ArrayDesc> &schemas)
    {
        std::vector<boost::shared_ptr<OperatorParamPlaceholder> > res;
        res.push_back(END_OF_VARIES_PARAMS());
        if(_parameters.size() == 0 || _parameters.size() == 1) {
            res.push_back(PARAM_CONSTANT("int64"));
        }
        return res;
    }

    ArrayDesc inferSchema(std::vector< ArrayDesc> schemas, boost::shared_ptr< Query> query)
    {
        assert(schemas.size() == 1);
        assert(_parameters.size() == 0);

        if (schemas[0].getAttributes(true).size() != 1)
            throw USER_EXCEPTION(SCIDB_SE_INFER_SCHEMA, SCIDB_LE_OP_EIGEN_ERROR1);
        if (schemas[0].getDimensions().size() != 2)
            throw USER_EXCEPTION(SCIDB_SE_INFER_SCHEMA, SCIDB_LE_OP_EIGEN_ERROR2);
        if ( schemas[0].getDimensions()[0].getLength() != schemas[0].getDimensions()[1].getLength())
            throw USER_EXCEPTION(SCIDB_SE_INFER_SCHEMA, SCIDB_LE_OP_EIGEN_ERROR3);

		Attributes atts;
		AttributeDesc multAttr((AttributeID)0, "v",  TID_DOUBLE, 0, 0);
		atts.push_back(multAttr);

		int rank;
		if(_parameters.size() == 1) {
			rank = evaluate(((boost::shared_ptr<OperatorParamLogicalExpression>&)_parameters[0])->getExpression(), query, TID_INT64).getInt64();
		} else {
			throw USER_EXCEPTION(SCIDB_SE_INFER_SCHEMA, SCIDB_LE_OP_EIGEN_ERROR5);
		}

		Dimensions dims;
		DimensionDesc d1 = schemas[0].getDimensions()[0];
		DimensionDesc d2 = schemas[0].getDimensions()[1];
		d1.setEndMax(rank);

        if (d2.getLength() == INFINITE_LENGTH || d2.getLength() == INFINITE_LENGTH)
            throw USER_EXCEPTION(SCIDB_SE_INFER_SCHEMA, SCIDB_LE_OP_EIGEN_ERROR4);

        dims.push_back(DimensionDesc(d1.getBaseName(), d1.getNamesAndAliases(), d1.getStartMin(), d1.getCurrStart(), d1.getCurrEnd(), rank, rank, 0));
        dims.push_back(DimensionDesc(d2.getBaseName(), d2.getNamesAndAliases(), d1.getStartMin(), d1.getCurrStart(), d2.getCurrEnd(), d2.getEndMax(), d2.getChunkInterval(), 0));

        return ArrayDesc("eigen", atts, dims);
	}

};

DECLARE_LOGICAL_OPERATOR_FACTORY(LogicalEigen, "eigen")

} //namespace scidb
