trace('calculateLW', edit = T, where = asNamespace("monocle3"))
trace("jaccard_coeff", edit = T, where = asNamespace("monocle3"))
change Matrix::rBind to rbind


```R

my_calculateLW <- function (cds, k, neighbor_graph, reduction_method, verbose = FALSE)
{
    if (verbose) {
        message("retrieve the matrices for Moran's I test...")
    }
    knn_res <- NULL
    principal_g <- NULL
    cell_coords <- reducedDims(cds)[[reduction_method]]
    if (neighbor_graph == "knn") {
        knn_res <- RANN::nn2(cell_coords, cell_coords, min(k +
            1, nrow(cell_coords)), searchtype = "standard")[[1]]
    }
    else if (neighbor_graph == "principal_graph") {
        pr_graph_node_coords <- cds@principal_graph_aux[[reduction_method]]$dp_mst
        principal_g <- igraph::get.adjacency(cds@principal_graph[[reduction_method]])[colnames(pr_graph_node_coords),
            colnames(pr_graph_node_coords)]
    }
    exprs_mat <- exprs(cds)
    if (neighbor_graph == "knn") {
        if (is.null(knn_res)) {
            knn_res <- RANN::nn2(cell_coords, cell_coords, min(k +
                1, nrow(cell_coords)), searchtype = "standard")[[1]]
        }
        links <- jaccard_coeff(knn_res[, -1], F)
        links <- links[links[, 1] > 0, ]
        relations <- as.data.frame(links)
        colnames(relations) <- c("from", "to", "weight")
        knn_res_graph <- igraph::graph.data.frame(relations,
            directed = T)
        knn_list <- lapply(1:nrow(knn_res), function(x) knn_res[x,
            -1])
        region_id_names <- colnames(cds)
        id_map <- 1:ncol(cds)
        names(id_map) <- id_map
        points_selected <- 1:nrow(knn_res)
        knn_list <- lapply(points_selected, function(x) id_map[as.character(knn_res[x,
            -1])])
    }
    else if (neighbor_graph == "principal_graph") {
        cell2pp_map <- cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_closest_vertex
        if (is.null(cell2pp_map)) {
            stop(paste("Error: projection matrix for each cell to principal",
                "points doesn't exist, you may need to rerun learn_graph"))
        }
        cell2pp_map <- cell2pp_map[row.names(cell2pp_map) %in%
            row.names(colData(cds)), , drop = FALSE]
        cell2pp_map <- cell2pp_map[colnames(cds), ]
        if (verbose) {
            message("Identify connecting principal point pairs ...")
        }
        knn_res <- RANN::nn2(cell_coords, cell_coords, min(k +
            1, nrow(cell_coords)), searchtype = "standard")[[1]]
        principal_g_tmp <- principal_g
        diag(principal_g_tmp) <- 1
        cell_membership <- as.factor(cell2pp_map)
        uniq_member <- sort(unique(cell_membership))
        membership_matrix <- Matrix::sparse.model.matrix(~cell_membership +
            0)
        colnames(membership_matrix) <- levels(uniq_member)
        feasible_space <- membership_matrix %*% Matrix::tcrossprod(principal_g_tmp[as.numeric(levels(uniq_member)),
            as.numeric(levels(uniq_member))], membership_matrix)
        links <- jaccard_coeff(knn_res[, -1], F)
        links <- links[links[, 1] > 0, ]
        relations <- as.data.frame(links)
        colnames(relations) <- c("from", "to", "weight")
        knn_res_graph <- igraph::graph.data.frame(relations,
            directed = T)
        tmp_a <- igraph::get.adjacency(knn_res_graph)
        block_size <- 10000
        num_blocks = ceiling(nrow(tmp_a)/block_size)
        if (verbose) {
            message("start calculating valid kNN graph ...")
        }
        tmp <- NULL
        for (j in 1:num_blocks) {
            if (j < num_blocks) {
                block_a <- tmp_a[((((j - 1) * block_size) + 1):(j *
                  block_size)), ]
                block_b <- feasible_space[((((j - 1) * block_size) +
                  1):(j * block_size)), ]
            }
            else {
                block_a <- tmp_a[((((j - 1) * block_size) + 1):(nrow(tmp_a))),
                  ]
                block_b <- feasible_space[((((j - 1) * block_size) +
                  1):(nrow(tmp_a))), ]
            }
            cur_tmp <- block_a * block_b
            if (is.null(tmp)) {
                tmp <- cur_tmp
            }
            else {
                tmp <- rbind(tmp, cur_tmp)
            }
        }
        if (verbose) {
            message("Calculating valid kNN graph, done ...")
        }
        region_id_names <- colnames(cds)
        id_map <- 1:ncol(cds)
        names(id_map) <- id_map
        knn_list <- slam::rowapply_simple_triplet_matrix(slam::as.simple_triplet_matrix(tmp),
            function(x) {
                res <- which(as.numeric(x) > 0)
                if (length(res) == 0)
                  res <- 0L
                res
            })
    }
    else {
        stop("Error: unrecognized neighbor_graph option")
    }
    names(knn_list) <- id_map[names(knn_list)]
    class(knn_list) <- "nb"
    attr(knn_list, "region.id") <- region_id_names
    attr(knn_list, "call") <- match.call()
    lw <- spdep::nb2listw(knn_list, zero.policy = TRUE)
    lw
}
```
getAnywhere("calculateLW")
unlockBinding("calculateLW", getNamespace("monocle3"))
assign("calculateLW", my_calculateLW, envir = asNamespace("monocle3"))
lockBinding("calculateLW", getNamespace("monocle3"))
