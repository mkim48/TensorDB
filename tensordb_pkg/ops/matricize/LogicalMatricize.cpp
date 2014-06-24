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
 * LogicalMatricize.cpp
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
 * @brief The operator: matricize().
 *
 * @par Synopsis:
 *    matricize(tensor, m)
 *      
 * @par Summary:
 *    Matricization transforms a tensor into a matrix along the given mode, m (m=1..N), N is the number of modes
 *
 */
class LogicalMatricize : public  LogicalOperator
{
public:
	LogicalMatricize(const std::string& logicalName, const std::string& alias):
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

		Attributes atts(1);
        TypeId type = schemas[0].getAttributes()[0].getType();
		AttributeDesc multAttr((AttributeID)0, "matricize", type, 0, 0);
		atts[0] = multAttr;
		Dimensions const& inputDims = schemas[0].getDimensions();
		Coordinate ndims = inputDims.size();

		Coordinate rowmode = 0;
		Coordinate colmode = 1;
        if (_parameters.size() == 1) {
            rowmode = evaluate(((boost::shared_ptr<OperatorParamLogicalExpression>&)_parameters[0])->getExpression(), query, TID_INT64).getInt64();
			if(rowmode < 1 || rowmode > ndims) {
				throw USER_EXCEPTION(SCIDB_SE_INFER_SCHEMA, SCIDB_LE_OP_MATRICIZE_ERROR1);
			}
 
			if(rowmode < ndims)
				colmode = ndims;
			else
				colmode = ndims-1;
			
			rowmode -= 1;
			colmode -= 1;
        } else {
			throw USER_EXCEPTION(SCIDB_SE_INFER_SCHEMA, SCIDB_LE_OP_MATRICIZE_ERROR2);
		}
        Coordinate len = 1;
        Coordinate chunklen = 1;
        for(Coordinate i=ndims-1;i>=0;i--) {
                if(i == rowmode) continue;
                len *= inputDims[i].getLength();
                chunklen *= inputDims[i].getChunkInterval();
        }
		Dimensions dims(2);
		dims[0] = inputDims[rowmode];
        dims[1] = DimensionDesc(inputDims[rowmode].getBaseName(),
                                1, 
                                len,
                                chunklen, 
                                0,
                                inputDims[rowmode].getType(),
                                inputDims[rowmode].getFlags(),
                                inputDims[rowmode].getMappingArrayName(),
                                inputDims[rowmode].getComment());

        return ArrayDesc("Matricize",atts,dims);
	}

};

DECLARE_LOGICAL_OPERATOR_FACTORY(LogicalMatricize, "matricize")

} //namespace scidb
