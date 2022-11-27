# need some internal greta functions accessible

dag_class <- .internals$inference$dag_class

check_dims <- .internals$checks$check_dims
check_in_family <- .internals$checks$check_in_family
check_positive <- .internals$checks$check_positive
check_unit <- .internals$checks$check_unit
check_tf_version <- .internals$checks$check_tf_version

has_representation <- .internals$greta_arrays$has_representation
representation <- .internals$greta_arrays$representation

distrib <- .internals$nodes$constructors$distrib
distribution_node <- .internals$nodes$node_classes$distribution_node

as.greta_array <- .internals$greta_arrays$as.greta_array

as_tf_function <- .internals$utils$greta_array_operations$as_tf_function

tf_lchoose <- .internals$tensors$tf_lchoose
tf_lbeta <- .internals$tensors$tf_lbeta
tf_flatten <- .internals$tensors$tf_flatten
tf_iprobit <- .internals$tensors$tf_iprobit
tf_as_float <- .internals$tensors$tf_as_float
tf_rowsums <- .internals$tensors$tf_rowsums
tf_sweep <- .internals$tensors$tf_sweep
tf_colsums <- .internals$tensors$tf_colsums
tf_float <- .internals$utils$misc$tf_float
tf_as_integer <- .internals$tensors$tf_as_integer

fl <- .internals$utils$misc$fl
op <- .internals$nodes$constructors$op
to_shape <- .internals$utils$misc$to_shape
expand_to_batch <- .internals$utils$misc$expand_to_batch
dummy_greta_array <- .internals$utils$dummy_arrays$dummy_greta_array
