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
 * MultiplyArray.h
 *
 *  Created on: Sep 12, 2010
 */

#ifndef __SPARSE_MATRIX_H__
#define __SPARSE_MATRIX_H__

#include "array/Array.h"
#include "util/FixAlloc.h"

namespace scidb {

struct SparseMatrixCell 
{
    double value;
    Coordinate coord;
    SparseMatrixCell* next;
};


class SparseRowMatrix 
{
    SparseMatrixCell** rows;
    Coordinate base;
    Allocator<SparseMatrixCell> allocator;
    
  public:
    SparseMatrixCell* operator[](Coordinate coord)
    { 
        return rows[coord - base];
    }

    void add(Coordinate row, Coordinate col, double value) 
    { 
        SparseMatrixCell* cell = allocator.allocate();
        cell->next = rows[row - base];
        rows[row - base] = cell;
        cell->value = value;
        cell->coord = col;
    }

    SparseRowMatrix(Coordinate first, uint64_t length) 
    {
        base = first;
        rows = (SparseMatrixCell**)calloc((size_t)length, sizeof(SparseMatrixCell*));
    }

    ~SparseRowMatrix() 
    {
        free(rows);
    }
};

class SparseColumnMatrix 
{
    SparseMatrixCell** cols;
    Coordinate base;
    Allocator<SparseMatrixCell> allocator;
    
  public:
    SparseMatrixCell* operator[](Coordinate coord)
    { 
        return cols[coord - base];
    }

    void add(Coordinate row, Coordinate col, double value) 
    { 
        SparseMatrixCell* cell = allocator.allocate();
        cell->next = cols[col - base];
        cols[col - base] = cell;
        cell->value = value;
        cell->coord = row;
    }

    SparseColumnMatrix(Coordinate first, uint64_t length) 
    { 
        base = first;
        cols = (SparseMatrixCell**)calloc((size_t)length, sizeof(SparseMatrixCell*));
    }

    ~SparseColumnMatrix() 
    {
        free(cols);
    }
};

}

#endif
        
