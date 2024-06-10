
###################################
# Auxiliary functions for cellQDA #
###################################

# (copied from VCR_auxiliary in classmap package)
# Hidden from the user:
checkLabels <- function(y, n, training = TRUE, levels = NULL) {
  if (training == TRUE) { # for training data
    if (!is.null(levels)) {
      stop("For training data, the levels argument is not used.")
    }
    whichy <- "y"
  } else {# for new data
    if (is.null(levels)) {
      stop("For new data, the argument levels is required")
    }
    whichy <- "ynew"
  }
  if (!(is.factor(y))) {
    stop(paste0("\n The response ", whichy, " is not a factor."))
  }
  if (is.null(attr(y, "levels"))) {
    stop(paste0("\n The response factor ", whichy, " has no levels.",
                "\n Please use as.factor() or factor() first."))
  }
  levelsattr <- attr(y, "levels")
  if (sum(is.na(levelsattr)) > 0) {
    stop(paste0("\n The levels attribute of ", whichy,
                " has at least one NA."))
  }
  dupls <- duplicated(levelsattr)
  if (sum(dupls) > 0) stop(paste0(
    "\n The label ", levelsattr[dupls == TRUE], " occurs more than ",
    " once in the levels attribute of ", whichy))
  if (length(y) != n) {
    stop(paste0("\n The response y should have length ", n,
                ", the number of cases."))
  }
  yv <- as.character(y) # this yv is not a factor any more.
  # When the factor y is malformed, the function stops here.
  # This happens when not all entries of y belong to levelsattr.
  indsv <- which(!is.na(yv)) # INDiceS of y that can be Visualized
  if (training == TRUE) { # for training data:
    if (length(indsv) == 0) {
      stop(paste0("The response factor y only contains NA's,"
                  , "\n so training is not possible."))
    }
    yv <- yv[indsv] # thus has at least 1 entry
    uniqy <- sort(unique(yv))
    for (g in seq_len(length(uniqy))) { # g = 1
      if (!(uniqy[g] %in% levelsattr)) {
        stop(paste0("\n The response label ", uniqy[g], " does not",
                    " appear in the levels attribute of y."))
      }
    }
    # From here on we know that all the yv are in levelsattr.
    if (length(uniqy) < 2) {
      stop(paste0("\n Only a single level occurs ( ", uniqy[1],
                  " ) so training is not possible."))
    }
    for (g in seq_len(length(levelsattr))) { # g = 1
      if (!(levelsattr[g] %in% uniqy)) {
        wnq(paste0("\nThe level ", levelsattr[g], " does not occur",
                   " among the training cases. A model trained",
                   "\non these data will not be able to",
                   " assign any objects to this class.\n"))
      }
    }
    # Create array "levels" with the labels that actually occur:
    levels <- levelsattr[which(levelsattr %in% uniqy)]; levels
    #
  } else {# for new data:
    if (length(indsv) > 0) {
      yv <- yv[indsv] # thus has at least 1 entry
      uniqy <- sort(unique(yv))
      for (g in seq_len(length(uniqy))) { # g = 1
        if (!(uniqy[g] %in% levelsattr)) {
          stop(paste0("\n The response label ", uniqy[g], " does not",
                      " appear in the levels attribute of ynew."))
        }
      }
      # From here on we know that all the yv are in levelsattr.
      badlevels <- uniqy[which(!(uniqy %in% levels))]; badlevels
      if (length(badlevels) > 0) { wnq(paste0(
        "\n The level ", badlevels, " occurs in ynew but",
        " was not in the training data.",
        "\n Such cases will be treated as if their ynew is NA",
        " so their response", "\n can be predicted, but",
        " these cases will not be in the class map.\n"))
        indsv <- indsv[which(!(yv %in% badlevels))]
        # the resulting indsv may again be empty
      }
    }
  }
  xlevelz <- levels # for use in the following function:
  lab2int <- function(labs){
    # labs is a vector of labels, that may contain NA's
    ints <- rep(NA, length(labs))
    indv <- which(!is.na(labs))
    labv <- labs[indv] # uses the final indsv above
    for (g in seq_len(length(xlevelz))) { # g = 1
      clinds <- indv[which(labv == xlevelz[g])] # in all of labs
      ints[clinds] <- g
    }
    ints
  }
  return(list(lab2int = lab2int, levels = levels, indsv = indsv))
}


# (copied from VCR_auxiliary in classmap package)
# Hidden from the user:
numericalCond <- function(S) {
  # Computes the inverse condition number of a matrix, which
  # is the ratio: smallest eigenvalue/largest eigenvalue .
  # Unlike the CN, this avoids division by (near) zero.
  # I call this the "numerical condition" (NC).
  #
  S <- as.matrix(S)
  if (length(which(is.na(S))) > 0) stop(" S contains NA's.")
  d <- nrow(S)
  if (ncol(S) != d) stop(" S is not square.")
  maxS <- max(abs(as.vector(S)))
  if (!isSymmetric(S, tol = 1e-10 * maxS)) stop(" S is not symmetric.")
  eigvals <- eigen(S, only.values = TRUE)$values
  eigvals[d] / eigvals[1]
}


# Hidden from the user: 
cellMCD <- function(X) {
  # Wrapper around cellWise::cellMCD
  cellmcd.out <- cellWise::cellMCD(X, noCits = 100, alpha = 0.75) # default: alpha = 0.75
  return(list(m = cellmcd.out$mu, S = cellmcd.out$S))
}

