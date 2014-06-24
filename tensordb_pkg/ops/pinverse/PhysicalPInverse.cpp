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
 * PhysicalPInverse.cpp
 *
 * @author: Mijung Kim (mijung.kim.1@asu.edu) 
 */
#include "query/Operator.h"
#include "array/Metadata.h"
#include "array/DelegateArray.h"
#include "network/NetworkManager.h"
#include <armadillo>

namespace scidb
{

class PhysicalPInverse: public PhysicalOperator
{
public:
	PhysicalPInverse(const std::string& logicalName, const std::string& physicalName, const Parameters& parameters, const ArrayDesc& schema):
	    PhysicalOperator(logicalName, physicalName, parameters, schema)
	{
		_schema = schema;
	}

	//This will be needed when we support inverse on multiple instances
    virtual bool changesDistribution(std::vector<ArrayDesc> const&) const
    {
        return true;
    }

    virtual ArrayDistribution getOutputDistribution(
            std::vector<ArrayDistribution> const&,
            std::vector<ArrayDesc> const&) const
    {
        return ArrayDistribution(psLocalInstance);
    }

	/***
	 */
	boost::shared_ptr< Array> execute(std::vector< boost::shared_ptr< Array> >& inputArrays,
            boost::shared_ptr<Query> query)
	{
        boost::shared_ptr<Array> input = inputArrays[0];
        if ( query->getInstancesCount() > 1) {
           uint64_t coordinatorID = (int64_t)query->getCoordinatorID() == -1 ? query->getInstanceID() : query->getCoordinatorID();
            input = redistribute(input, query, psLocalInstance, "", coordinatorID);
            if ( query->getInstanceID() != coordinatorID) { 
                return boost::shared_ptr<Array>(new MemArray(_schema));
            }
        }
        Dimensions const& dims = _schema.getDimensions();
        size_t length = dims[0].getLength();
        Coordinates first(2);
        Coordinates last(2);
        first[0] = dims[0].getStart();
        first[1] = dims[1].getStart();
        last[0] = dims[0].getEndMax();
        last[1] = dims[1].getEndMax();
        double* matrix = new double[length*length];
        input->extractData(0, matrix, first, last);
	
		arma::mat A(matrix,length,length,false);
		arma::mat B = arma::pinv(A);
		double* p = (double*)B.memptr();
		for(int i=0;i<length*length;i++) {
			matrix[i] = p[i];
		}
        shared_array<char> buf(reinterpret_cast<char*>(matrix));
		return boost::shared_ptr<Array>(new SplitArray(_schema, buf, first, last));
	}

private:
	ArrayDesc _schema;
};

DECLARE_PHYSICAL_OPERATOR_FACTORY(PhysicalPInverse, "pinverse", "PhysicalPInverse")

} //namespace scidb
