/*
**
* BEGIN_COPYRIGHT
*
* This file is part of TensorDB.
* Copyright (C) 2014 Mijung Kim, K. Selcuk Candan All right Reserved.
*
* Reference:
*  [1] Mijung Kim (2014). TensorDB and Tensor-based Relational Model (TRM) for Efficient Tensor-Relational Operations. Ph.D. Thesis. Arizona State University.
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
 * LogicalMultiplyComp.cpp
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
 * @brief The operator: multiply_comp().
 *
 * @par Synopsis:
 *   multiply_comp( leftArray, rightArray )
 *
 * @par Summary:
 *   Approximate matrix multiplication using compressed matrix multiplication.
 *   Note that unlike 'multiply_row', 'multiply_comp' performs matrix multiplication of entries of 'columns' of each matrix 
 *
 */
class LogicalMultiplyComp : public  LogicalOperator
{
public:
	LogicalMultiplyComp(const std::string& logicalName, const std::string& alias):
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

        if (!hasSingleAttribute(schemas[0]) || !hasSingleAttribute(schemas[1]))
            throw USER_EXCEPTION(SCIDB_SE_INFER_SCHEMA, SCIDB_LE_OP_MULTIPLY_ERROR2);
        if (schemas[0].getDimensions().size() != 2 || schemas[1].getDimensions().size() != 2)
            throw USER_EXCEPTION(SCIDB_SE_INFER_SCHEMA, SCIDB_LE_OP_MULTIPLY_ERROR3);

        if (schemas[0].getDimensions()[0].getLength() == INFINITE_LENGTH
                || schemas[0].getDimensions()[1].getLength() == INFINITE_LENGTH
                || schemas[1].getDimensions()[0].getLength() == INFINITE_LENGTH
                || schemas[1].getDimensions()[1].getLength() == INFINITE_LENGTH)
            throw USER_EXCEPTION(SCIDB_SE_INFER_SCHEMA, SCIDB_LE_OP_MULTIPLY_ERROR4);

        if (schemas[0].getDimensions()[1].getLength() != schemas[1].getDimensions()[0].getLength()
                || schemas[0].getDimensions()[1].getStart() != schemas[1].getDimensions()[0].getStart())
            throw USER_EXCEPTION(SCIDB_SE_INFER_SCHEMA, SCIDB_LE_OP_MULTIPLY_ERROR5);

        // FIXME: This condition needs to go away later
        if (schemas[0].getDimensions()[1].getChunkInterval() != schemas[1].getDimensions()[0].getChunkInterval())
            throw USER_EXCEPTION(SCIDB_SE_INFER_SCHEMA, SCIDB_LE_OP_MULTIPLY_ERROR6);
        if (schemas[0].getAttributes()[0].getType() != schemas[1].getAttributes()[0].getType())
            throw USER_EXCEPTION(SCIDB_SE_INFER_SCHEMA, SCIDB_LE_OP_MULTIPLY_ERROR7);
        if (schemas[0].getAttributes()[0].isNullable() || schemas[1].getAttributes()[0].isNullable())
            throw USER_EXCEPTION(SCIDB_SE_INFER_SCHEMA, SCIDB_LE_OP_MULTIPLY_ERROR8);

		Attributes atts(1);
        TypeId type = schemas[0].getAttributes()[0].getType();
		AttributeDesc multAttr((AttributeID)0, "multiply_comp", type, 0, 0);
		atts[0] = multAttr;

		Dimensions dims(2);
		DimensionDesc const& d1 = schemas[0].getDimensions()[0];
		dims[0] = DimensionDesc(d1.getBaseName(), 
                                d1.getNamesAndAliases(), 
                                d1.getStartMin(), 
                                d1.getCurrStart(), 
                                d1.getCurrEnd(), 
                                d1.getEndMax(), 
                                d1.getChunkInterval(), 
                                0, 
                                d1.getType(), 
                                d1.getFlags(), 
                                d1.getMappingArrayName(), 
                                d1.getComment(),
                                d1.getFuncMapOffset(),
                                d1.getFuncMapScale());

		DimensionDesc const& d2 = schemas[1].getDimensions()[1];
		dims[1] = DimensionDesc(d1.getBaseName() == d2.getBaseName() ? d1.getBaseName() + "2" : d2.getBaseName(), 
                                d2.getNamesAndAliases(), 
                                d2.getStartMin(), 
                                d2.getCurrStart(), 
                                d2.getCurrEnd(), 
                                d2.getEndMax(), 
                                d2.getChunkInterval(), 
                                0, 
                                d2.getType(), 
                                d2.getFlags(), 
                                d2.getMappingArrayName(), 
                                d2.getComment(),
                                d2.getFuncMapOffset(),
                                d2.getFuncMapScale());

        return ArrayDesc("MultiplyComp",atts,dims);
	}

};

DECLARE_LOGICAL_OPERATOR_FACTORY(LogicalMultiplyComp, "multiply_comp")

} //namespace scidb
