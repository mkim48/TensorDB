ERRMSG(SCIDB_LE_OP_COPYARRAY_ERROR1,                "This target matrix cannot be RLE");
ERRMSG(SCIDB_LE_OP_COPYARRAY_ERROR2,                "Matrix must have 1 attrbiute");
ERRMSG(SCIDB_LE_OP_COPYARRAY_ERROR3,                "Cannot copy unbounded matrices");

ERRMSG(SCIDB_LE_OP_MATRICIZE_ERROR1,                "The input mode (m) for mode-m matricization should be >= 1 and <= N, the number of the modes");
ERRMSG(SCIDB_LE_OP_MATRICIZE_ERROR2,                "One input mode (m) for mode-m matricization should be entered");

ERRMSG(SCIDB_LE_OP_EIGEN_ERROR1,                  "Matrix must contain 1 attribute");
ERRMSG(SCIDB_LE_OP_EIGEN_ERROR2,                  "Only 2-D matrices can be input");
ERRMSG(SCIDB_LE_OP_EIGEN_ERROR3,                  "Eigen decomposition only works on square matrices");
ERRMSG(SCIDB_LE_OP_EIGEN_ERROR4,                  "Cannot perform on an unbounded matrix");
ERRMSG(SCIDB_LE_OP_EIGEN_ERROR5,                  "Usage: eigen(<array>, <target rank)");
