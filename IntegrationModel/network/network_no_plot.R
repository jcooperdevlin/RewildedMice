network_no_plot=function (mat, comp = NULL, blocks = c(1, 2), cutoff = NULL, 
    row.names = TRUE, col.names = TRUE, block.var.names = TRUE, 
    color.node = NULL, shape.node = NULL, cex.node.name = 1, 
    color.edge = color.GreenRed(100), lty.edge = "solid", lwd.edge = 1, 
    show.edge.labels = FALSE, cex.edge.label = 1, show.color.key = TRUE, 
    symkey = TRUE, keysize = c(1, 1), keysize.label = 1, breaks, 
    interactive = FALSE, layout.fun = NULL, save = NULL, name.save = NULL) 
{
    arg.call = match.call()
    user.arg = names(arg.call)[-1]
    err = tryCatch(mget(names(formals()), sys.frame(sys.nframe())), 
        error = function(e) e)
    if ("simpleError" %in% class(err)) 
        stop(err[[1]], ".", call. = FALSE)
    function.arg = names(mget(names(formals()), sys.frame(sys.nframe())))
    not.arg = !(user.arg %in% function.arg)
    if (any(not.arg)) {
        unused.arg = user.arg[not.arg]
        not.arg = which(not.arg) + 1
        output = rep("", length(not.arg))
        for (i in 1:length(not.arg)) {
            output[i] = paste0(unused.arg[i], " = ", arg.call[[not.arg[i]]])
        }
        output = paste0("(", paste(output, collapse = ", "), 
            ").")
        msg = "unused argument "
        if (length(not.arg) > 1) 
            msg = "unused arguments "
        stop(msg, output, call. = FALSE)
    }
    if (!is.null(save)) {
        if (!save %in% c("jpeg", "tiff", "png", "pdf")) 
            stop("'save' must be one of 'jpeg', 'png', 'tiff' or 'pdf'.", 
                call. = FALSE)
    }
    if (!is.null(name.save)) {
        if (!is.character(name.save) || length(name.save) > 1) 
            stop("'name.save' must be a character.", call. = FALSE)
    }
    else {
        if (!is.null(save)) 
            name.save = paste0("network_", gsub(".", "_", deparse(substitute(mat)), 
                fixed = TRUE))
    }
    if (!is.null(save)) {
        while (dev.cur() > 1) dev.off()
        if (save == "jpeg") 
            jpeg(filename = paste0(name.save, ".jpeg"), res = 600, 
                width = 4000, height = 4000)
        if (save == "png") 
            jpeg(filename = paste0(name.save, ".png"), res = 600, 
                width = 4000, height = 4000)
        if (save == "tiff") 
            tiff(filename = paste0(name.save, ".tiff"), res = 600, 
                width = 4000, height = 4000)
        if (save == "pdf") 
            pdf(file = paste0(name.save, ".pdf"))
    }
    class.object = class(mat)
    object.pls = c("mixo_pls", "mixo_spls", "mixo_mlspls")
    object.rcc = "rcc"
    object.blocks = c("sgcca", "rgcca")
    if (!any(class.object %in% c(object.pls, object.rcc, object.blocks, 
        "matrix"))) 
        stop(" 'network' is only implemented for the following objects: matrix, pls, plsda, spls, splsda, rcc, sgcca, rgcca, sgccda", 
            call. = FALSE)
    if (any(class.object %in% c(object.rcc, object.pls))) {
        p = ncol(mat$X)
        if (any(class.object == "DA")) 
            mat$Y = mat$ind.mat
        q = ncol(mat$Y)
        n = nrow(mat$X)
        ncomp = mat$ncomp
        if (is.null(comp)) 
            comp = 1:mat$ncomp
        if (length(comp) == 1) {
            if (comp > ncomp) {
                stop("the elements of 'comp' must be smaller than or equal to ", 
                  ncomp, ".", call. = FALSE)
            }
            else if (!is.numeric(comp) || comp <= 0) {
                stop("invalid value for 'comp'.", call. = FALSE)
            }
        }
        if (length(comp) > 1) {
            if (length(comp) > ncomp) 
                stop("the length of 'comp' must be smaller than or equal to ", 
                  ncomp, ".", call. = FALSE)
            if (!is.numeric(comp) || any(comp < 1)) 
                stop("invalid vector for 'comp'.", call. = FALSE)
            if (any(comp > ncomp)) 
                stop("the elements of 'comp' must be smaller or equal than ", 
                  ncomp, ".", call. = FALSE)
        }
        comp = round(comp)
        row.names.plot = TRUE
        if (is.logical(row.names)) {
            if (!isTRUE(row.names)) {
                row.names.plot = FALSE
            }
            row.names = mat$names$colnames$X
        }
        else {
            row.names = as.vector(row.names)
            if (length(unique(row.names)) != p) 
                stop("'row.names' must be a character vector of ", 
                  p, " unique entries.", call. = FALSE)
        }
        if (row.names.plot == TRUE) {
            row.names.plot = row.names
        }
        else {
            row.names.plot = rep("", p)
        }
        col.names.plot = TRUE
        if (is.logical(col.names)) {
            if (!isTRUE(col.names)) {
                col.names.plot = FALSE
            }
            col.names = mat$names$colnames$Y
        }
        else {
            col.names = as.vector(col.names)
            if (length(col.names) != q) 
                stop("'col.names' must be a character vector of ", 
                  q, " unique entries.", call. = FALSE)
        }
        if (col.names.plot == TRUE) {
            col.names.plot = col.names
        }
        else {
            col.names.plot = rep("", q)
        }
        if (any(class.object %in% object.rcc)) {
            bisect = mat$variates$X[, comp] + mat$variates$Y[, 
                comp]
            cord.X = cor(mat$X, bisect, use = "pairwise")
            cord.Y = cor(mat$Y, bisect, use = "pairwise")
            mat = cord.X %*% t(cord.Y)
        }
        else if (any(class.object %in% object.pls)) {
            if (all(class(mat) %in% "mixo_pls")) {
                keep.X = rep(TRUE, p)
                keep.Y = rep(TRUE, q)
            }
            else {
                keep.X = apply(abs(mat$loadings$X[, comp, drop = FALSE]), 
                  1, sum) > 0
                keep.Y = apply(abs(mat$loadings$Y[, comp, drop = FALSE]), 
                  1, sum) > 0
                row.names = row.names[keep.X]
                col.names = col.names[keep.Y]
            }
            if (mat$mode == "canonical") {
                cord.X = cor(mat$X[, keep.X], mat$variates$X[, 
                  comp], use = "pairwise")
                cord.Y = cor(mat$Y[, keep.Y], mat$variates$Y[, 
                  comp], use = "pairwise")
            }
            else {
                cord.X = cor(mat$X[, keep.X], mat$variates$X[, 
                  comp], use = "pairwise")
                cord.Y = cor(mat$Y[, keep.Y], mat$variates$X[, 
                  comp], use = "pairwise")
            }
            mat = cord.X %*% t(cord.Y)
        }
    }
    else if (any(class.object %in% object.blocks)) {
        if (any(class.object == "DA")) {
            mat$names$blocks = mat$names$blocks[-mat$indY]
            mat$names$colnames = mat$names$colnames[-mat$indY]
            mat$ncomp = mat$ncomp[-mat$indY]
        }
        if (is.null(blocks)) {
            if (any(mat$ncomp > 1)) {
                blocks = mat$names$blocks[which(mat$ncomp > 1)]
            }
            else {
                stop(("The number of components for each block is 1. The number of components must be superior or equal to 2."), 
                  call. = FALSE)
            }
        }
        else if (is.numeric(blocks) & min(blocks) > 0 & max(blocks) <= 
            length(mat$names$blocks)) {
            blocks = mat$names$blocks[blocks]
        }
        else if (is.character(blocks)) {
            if (!all(blocks %in% mat$names$blocks)) 
                stop("One element of 'blocks' does not match with the names of the blocks")
        }
        else {
            stop("Incorrect value for 'blocks", call. = FALSE)
        }
        if (is.null(comp)) {
            comp = vector("list", length(blocks))
            names(comp) = blocks
            for (i in blocks) comp[[i]] = 1:mat$ncomp[i]
        }
        if (is.list(comp)) {
            if (length(comp) != length(blocks)) 
                stop("'comp' must be either NULL a list of length ", 
                  length(blocks), ".", call. = FALSE)
            if (!all(blocks %in% names(comp))) 
                stop("names of 'comp' must be from {", paste(blocks, 
                  collapse = ", "), "}.", call. = FALSE)
            for (i in blocks) {
                if (any(!is.finite(comp[[i]]))) 
                  stop("invalid value for 'comp' of the block '", 
                    i, "'.", call. = FALSE)
                if (any(comp[[i]] > mat$ncomp[i])) 
                  stop("the elements of 'comp' for block '", 
                    i, "' must be smaller or equal than ", mat$ncomp[i], 
                    ".", call. = FALSE)
                if (any(comp[[i]] < 1)) 
                  stop("invalid value for 'comp' of the block '", 
                    i, "'.", call. = FALSE)
            }
        }
        else {
            stop("'comp' must be either NULL or a list of length ", 
                length(blocks), ".", call. = FALSE)
        }
        num.var = unlist(lapply(mat$X[blocks], ncol))
        if (is.logical(block.var.names)) {
            if (length(block.var.names) == 1) 
                block.var.names = rep(block.var.names, length(blocks))
            if (length(block.var.names) != length(blocks)) 
                stop("'block.var.names' must be a logical vector of length 1 or ", 
                  length(blocks), ", or a list of length ", length(blocks), 
                  ".", call. = FALSE)
            vec = (which(block.var.names == FALSE))
            block.var.names = mat$names$colnames
            for (i in 1:length(blocks)) {
                if (i %in% vec) 
                  block.var.names[[blocks[i]]] = rep(" ", length(mat$names$colnames[[blocks[i]]]))
            }
        }
        else {
            if (is.list(block.var.names)) {
                if (length(block.var.names) != length(blocks)) {
                  stop("'block.var.names' must be a logical vector or a list of length ", 
                    length(blocks), ".", call. = FALSE)
                }
                else {
                  if (!all(unlist(lapply(block.var.names, is.vector)))) 
                    stop("each component of 'block.var.names' must be a vector.", 
                      call. = FALSE)
                  block.var.names.length = unlist(lapply(block.var.names, 
                    length))
                  if (any(block.var.names.length != num.var)) 
                    stop("components of 'block.var.names' must be vectors of length ", 
                      paste(num.var, collapse = ", "), ".", call. = FALSE)
                }
            }
            else {
                stop("'block.var.names' must be either a logical value or a list of length ", 
                  length(blocks), ".", call. = FALSE)
            }
        }
        coord = M_block = list()
        j = 1
        if (is(mat, "sgcca")) {
            for (k in blocks) {
                if (length(comp[[k]]) > 1) {
                  keep = (apply(abs(mat$loadings[[k]][, comp[[k]]]), 
                    1, sum) > 0)
                }
                else {
                  keep = abs(mat$loadings[[k]][, comp[[k]]]) > 
                    0
                }
                coord[[j]] = cor(mat$X[[k]][, keep], mat$variates[[k]][, 
                  comp[[k]]], use = "pairwise")
                j = j + 1
            }
        }
        else {
            for (k in blocks) {
                coord[[j]] = cor(mat$X[[k]], mat$variates[[k]][, 
                  comp[[k]]], use = "pairwise")
                j = j + 1
            }
        }
        node.X = node.Y = w = NULL
        l = 1
        for (j in 1:(length(blocks) - 1)) {
            for (k in (j + 1):length(blocks)) {
                if (!any(comp[[blocks[j]]] %in% comp[[blocks[k]]])) 
                  stop("comp of block ", blocks[j], " is ", comp[[blocks[j]]], 
                    " but comp of block ", blocks[k], " is ", 
                    comp[[blocks[k]]], call. = FALSE)
                int.comp = intersect(comp[[blocks[j]]], comp[[blocks[k]]])
                object = coord[[j]][, comp[[blocks[j]]] %in% 
                  int.comp] %*% t(coord[[k]][, comp[[blocks[k]]] %in% 
                  int.comp])
                M_block[[l]] = object
                l = l + 1
                X = rownames(coord[[j]])
                Y = rownames(coord[[k]])
                rep.X = rep(X, each = length(Y))
                rep.Y = rep(Y, length(X))
                node.X = c(node.X, rep.X)
                node.Y = c(node.Y, rep.Y)
                w = c(w, as.vector(t(object)))
            }
        }
    }
    return(mat)
}

