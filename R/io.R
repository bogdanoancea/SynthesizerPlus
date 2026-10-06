# Reading and writing data in many formats -----------------------------------

.format_table <- function() {
  data.frame(
    format = c("csv", "tsv", "rds", "rdata", "parquet", "feather", "hdf5",
               "json", "ndjson", "excel", "spss", "stata", "sas", "xpt", "fst"),
    extensions = c("csv", "tsv, tab, txt", "rds", "rdata, rda", "parquet, pq",
                   "feather, arrow, ipc", "h5, hdf5", "json", "ndjson, jsonl",
                   "xlsx, xls", "sav, zsav", "dta", "sas7bdat", "xpt", "fst"),
    read = c("base (data.table if installed)", "base (data.table if installed)",
             "base", "base", "arrow or nanoparquet", "arrow", "hdf5r",
             "jsonlite", "jsonlite", "readxl", "haven", "haven", "haven",
             "haven", "fst"),
    write = c("base (data.table if installed)", "base (data.table if installed)",
              "base", "base", "arrow or nanoparquet", "arrow", "hdf5r",
              "jsonlite", "jsonlite", "writexl or openxlsx", "haven", "haven",
              "-", "haven", "fst"),
    stringsAsFactors = FALSE
  )
}

.ext_map <- c(
  csv = "csv", tsv = "tsv", tab = "tsv", txt = "tsv", rds = "rds",
  rdata = "rdata", rda = "rdata", parquet = "parquet", pq = "parquet",
  feather = "feather", arrow = "feather", ipc = "feather", h5 = "hdf5",
  hdf5 = "hdf5", he5 = "hdf5", json = "json", ndjson = "ndjson",
  jsonl = "ndjson", xlsx = "excel", xls = "excel", sav = "spss",
  zsav = "spss", dta = "stata", sas7bdat = "sas", xpt = "xpt", fst = "fst"
)

#' Supported file formats
#'
#' Lists the file formats understood by [read_data()], [write_data()] and
#' [generate_to_file()], with their file extensions and the package used as
#' backend. Backends other than base R are optional (`Suggests`) and only
#' needed for the corresponding format. Text formats (`csv`, `tsv`, `json`,
#' `ndjson`) can additionally be compressed by appending `.gz`, `.bz2` or
#' `.xz` to the file name.
#'
#' @param available_only If `TRUE`, only formats whose backend is installed
#'   are returned.
#' @return A data frame with columns `format`, `extensions`, `read`,
#'   `write` and `available`.
#' @examples
#' supported_formats()
#' @export
supported_formats <- function(available_only = FALSE) {
  tab <- .format_table()
  tab$available <- vapply(tab$format, .backend_available, logical(1))
  if (available_only) tab <- tab[tab$available, , drop = FALSE]
  rownames(tab) <- NULL
  tab
}

.backend_available <- function(format) {
  has <- function(p) requireNamespace(p, quietly = TRUE)
  switch(format,
    csv = , tsv = , rds = , rdata = TRUE,
    parquet = has("arrow") || has("nanoparquet"),
    feather = has("arrow"),
    hdf5 = has("hdf5r"),
    json = , ndjson = has("jsonlite"),
    excel = has("readxl") || has("writexl") || has("openxlsx"),
    spss = , stata = , sas = , xpt = has("haven"),
    fst = has("fst"),
    FALSE
  )
}

# Split "file.csv.gz" into format "csv" and compression "gz"
.detect_format <- function(path, format = NULL) {
  comp <- NULL
  base <- basename(path)
  m <- regmatches(base, regexec("\\.(gz|bz2|xz)$", base, ignore.case = TRUE))[[1L]]
  if (length(m)) {
    comp <- tolower(m[2L])
    base <- sub("\\.(gz|bz2|xz)$", "", base, ignore.case = TRUE)
  }
  if (is.null(format)) {
    ext <- tolower(tools::file_ext(base))
    format <- unname(.ext_map[ext])
    if (is.na(format)) {
      stop(sprintf("Cannot infer the file format of '%s'. Use the 'format' argument; see supported_formats().",
                   path), call. = FALSE)
    }
  } else {
    format <- tolower(format)
    if (format %in% names(.ext_map)) format <- unname(.ext_map[format])
    if (!format %in% .format_table()$format) {
      stop("Unknown format '", format, "'. See supported_formats().", call. = FALSE)
    }
  }
  if (!is.null(comp) && !format %in% c("csv", "tsv", "json", "ndjson")) {
    stop("Compression suffixes are only supported for csv, tsv, json and ndjson.", call. = FALSE)
  }
  list(format = format, compression = comp)
}

.open_con <- function(path, compression, mode) {
  switch(compression %||% "none",
    gz = gzfile(path, mode),
    bz2 = bzfile(path, mode),
    xz = xzfile(path, mode),
    file(path, mode)
  )
}

#' Read a data set from a file
#'
#' @description
#' `read_data()` reads tabular data from many formats commonly used in data
#' science (see [supported_formats()]) and returns a plain `data.frame`. The
#' format is inferred from the file extension unless `format` is given.
#'
#' If `path` is a directory (for example as written by [generate_to_file()]),
#' all files in it with a recognised extension are read and row-bound.
#'
#' Because text formats (CSV, JSON, Excel) do not store column types, a
#' `template` data frame can be supplied: columns are then coerced to the
#' classes of the template (factor levels, `Date`, `POSIXct`, integer, ...)
#' with [match_types()].
#'
#' @param path File or directory path.
#' @param format Optional format name or extension, e.g. `"parquet"`.
#' @param template Optional data frame whose column classes are imposed on
#'   the result.
#' @param object For `rdata` files: name of the object to load (default: the
#'   first data frame in the file).
#' @param dataset For `hdf5` files: name of the group (written by
#'   [write_data()]) or of a 1- or 2-dimensional dataset to read.
#'   Default `"data"`.
#' @param sheet For Excel files: sheet name or number.
#' @param ... Further arguments passed to the backend reader.
#' @return A `data.frame`.
#' @seealso [write_data()], [supported_formats()], [match_types()]
#' @examples
#' f <- tempfile(fileext = ".csv")
#' write_data(iris, f)
#' x <- read_data(f, template = iris)
#' str(x)
#'
#' f2 <- tempfile(fileext = ".rds")
#' write_data(iris, f2)
#' identical(read_data(f2), iris)
#' @export
read_data <- function(path, format = NULL, template = NULL, object = NULL,
                      dataset = "data", sheet = 1L, ...) {
  if (!is.character(path) || length(path) != 1L) {
    stop("'path' must be a single file path.", call. = FALSE)
  }
  if (dir.exists(path)) {
    files <- sort(list.files(path, full.names = TRUE))
    files <- files[vapply(files, function(f) {
      !inherits(try(.detect_format(f, format), silent = TRUE), "try-error")
    }, logical(1))]
    if (!length(files)) stop("No readable data files found in directory '", path, "'.", call. = FALSE)
    parts <- lapply(files, read_data, format = format, template = template,
                    object = object, dataset = dataset, sheet = sheet, ...)
    out <- do.call(rbind, parts)
    attr(out, "row.names") <- .set_row_names(nrow(out))
    return(out)
  }
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  fmt <- .detect_format(path, format)
  out <- switch(fmt$format,
    csv = .read_delim(path, ",", fmt$compression, ...),
    tsv = .read_delim(path, "\t", fmt$compression, ...),
    rds = readRDS(path),
    rdata = .read_rdata(path, object),
    parquet = .read_parquet(path, ...),
    feather = {
      .require("arrow", "to read Feather/Arrow IPC files")
      arrow::read_feather(path, ...)
    },
    hdf5 = .read_hdf5(path, dataset),
    json = {
      .require("jsonlite", "to read JSON files")
      con <- .open_con(path, fmt$compression, "r")
      on.exit(close(con), add = TRUE)
      jsonlite::fromJSON(paste(readLines(con, warn = FALSE), collapse = "\n"),
                         simplifyVector = TRUE, ...)
    },
    ndjson = {
      .require("jsonlite", "to read NDJSON files")
      con <- .open_con(path, fmt$compression, "r")
      on.exit(close(con), add = TRUE)
      jsonlite::stream_in(con, verbose = FALSE, ...)
    },
    excel = {
      .require("readxl", "to read Excel files")
      readxl::read_excel(path, sheet = sheet, ...)
    },
    spss = .read_haven(path, "spss", ...),
    stata = .read_haven(path, "stata", ...),
    sas = .read_haven(path, "sas", ...),
    xpt = .read_haven(path, "xpt", ...),
    fst = {
      .require("fst", "to read fst files")
      fst::read_fst(path, ...)
    }
  )
  out <- .finalize_df(out)
  if (!is.null(template)) out <- match_types(out, template)
  out
}

.finalize_df <- function(x) {
  if (is.matrix(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!is.data.frame(x)) {
    x <- tryCatch(as.data.frame(x, stringsAsFactors = FALSE),
                  error = function(e) stop("The file does not contain tabular data.", call. = FALSE))
  }
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  for (nm in names(x)) {
    col <- x[[nm]]
    if (inherits(col, "IDate")) {
      col <- structure(as.numeric(unclass(col)), class = "Date")
      x[[nm]] <- col
    }
    # strip label/format attributes left by readers (haven, arrow)
    if (!is.null(attr(col, "label")) || !is.null(attr(col, "format.spss")) ||
        !is.null(attr(col, "format.stata")) || !is.null(attr(col, "format.sas"))) {
      attr(col, "label") <- NULL
      attr(col, "format.spss") <- NULL
      attr(col, "format.stata") <- NULL
      attr(col, "format.sas") <- NULL
      attr(col, "display_width") <- NULL
      x[[nm]] <- col
    }
  }
  x
}

.read_delim <- function(path, sep, compression, ...) {
  if (is.null(compression) && requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(data.table::fread(path, sep = sep, data.table = FALSE,
                                           showProgress = FALSE,
                                           na.strings = c("NA", ""), ...)))
  }
  con <- .open_con(path, compression, "r")
  on.exit(close(con), add = TRUE)
  utils::read.table(con, header = TRUE, sep = sep, quote = "\"",
                    stringsAsFactors = FALSE, na.strings = c("NA", ""),
                    check.names = FALSE, comment.char = "", ...)
}

.read_rdata <- function(path, object) {
  env <- new.env(parent = emptyenv())
  objs <- load(path, envir = env)
  if (!is.null(object)) {
    if (!object %in% objs) stop("Object '", object, "' not found in ", path, call. = FALSE)
    return(get(object, envir = env))
  }
  is_df <- vapply(objs, function(o) is.data.frame(get(o, envir = env)), logical(1))
  if (!any(is_df)) stop("No data frame found in ", path, call. = FALSE)
  get(objs[is_df][1L], envir = env)
}

.read_parquet <- function(path, ...) {
  if (requireNamespace("arrow", quietly = TRUE)) {
    return(arrow::read_parquet(path, ...))
  }
  if (requireNamespace("nanoparquet", quietly = TRUE)) {
    return(nanoparquet::read_parquet(path, ...))
  }
  stop("Reading Parquet files requires the 'arrow' or 'nanoparquet' package.", call. = FALSE)
}

.read_haven <- function(path, kind, ...) {
  .require("haven", "to read SPSS/Stata/SAS files")
  x <- switch(kind,
    spss = haven::read_sav(path, ...),
    stata = haven::read_dta(path, ...),
    sas = haven::read_sas(path, ...),
    xpt = haven::read_xpt(path, ...)
  )
  x <- as.data.frame(x)
  for (nm in names(x)) {
    if (inherits(x[[nm]], "haven_labelled")) {
      x[[nm]] <- haven::as_factor(x[[nm]], levels = "labels")
    }
    x[[nm]] <- haven::zap_formats(haven::zap_label(x[[nm]]))
    if (is.character(x[[nm]])) x[[nm]][x[[nm]] == ""] <- NA
  }
  x
}

#' Write a data set to a file
#'
#' @description
#' `write_data()` writes a data frame to any format listed by
#' [supported_formats()]. The format is inferred from the file extension
#' unless `format` is given; text formats are compressed when the file name
#' ends in `.gz`, `.bz2` or `.xz`.
#'
#' HDF5 files are written with one dataset per column inside a group
#' (`dataset`, default `"data"`), with attributes that record column order,
#' factor levels, `Date`/`POSIXct` classes, time zones and missing strings,
#' so that [read_data()] restores the original column types.
#'
#' Date-times (`POSIXct`) keep their instant in every format. Text formats
#' (CSV, TSV, JSON, NDJSON) store them as ISO 8601 strings in UTC; Excel,
#' SPSS, Stata and SAS files, which store clock time without a time zone,
#' receive the UTC clock time; Parquet, Feather, fst and HDF5 store the
#' instant together with the time zone. Pass the original data as `template`
#' to [read_data()] to restore the original time zone.
#'
#' @param x A data frame (matrices are converted).
#' @param path Output file path.
#' @param format Optional format name or extension.
#' @param overwrite Overwrite an existing file?
#' @param object For `rdata`: name under which the object is saved.
#' @param dataset For `hdf5`: name of the group to write.
#' @param ... Further arguments passed to the backend writer.
#' @return The path, invisibly.
#' @seealso [read_data()], [generate_to_file()]
#' @examples
#' f <- tempfile(fileext = ".tsv.gz")
#' write_data(mtcars, f)
#' head(read_data(f))
#' @export
write_data <- function(x, path, format = NULL, overwrite = TRUE,
                       object = "data", dataset = "data", ...) {
  if (is.matrix(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!is.data.frame(x)) stop("'x' must be a data frame or matrix.", call. = FALSE)
  if (!is.character(path) || length(path) != 1L) stop("'path' must be a single file path.", call. = FALSE)
  if (file.exists(path) && !.check_flag(overwrite, "overwrite")) {
    stop("File exists: ", path, " (use overwrite = TRUE).", call. = FALSE)
  }
  fmt <- .detect_format(path, format)
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  switch(fmt$format,
    csv = .write_delim(x, path, ",", fmt$compression, append = FALSE, ...),
    tsv = .write_delim(x, path, "\t", fmt$compression, append = FALSE, ...),
    rds = saveRDS(x, path, ...),
    rdata = {
      env <- new.env(parent = emptyenv())
      assign(object, x, envir = env)
      save(list = object, envir = env, file = path, ...)
    },
    parquet = .write_parquet(x, path, ...),
    feather = {
      .require("arrow", "to write Feather/Arrow IPC files")
      arrow::write_feather(x, path, ...)
    },
    hdf5 = .write_hdf5(x, path, dataset),
    json = {
      .require("jsonlite", "to write JSON files")
      con <- .open_con(path, fmt$compression, "w")
      on.exit(close(con), add = TRUE)
      writeLines(jsonlite::toJSON(.posix_to_iso(x), dataframe = "rows", digits = NA, na = "null",
                                  POSIXt = "ISO8601", Date = "ISO8601", ...), con)
    },
    ndjson = .write_ndjson(x, path, fmt$compression, append = FALSE, ...),
    excel = {
      if (tolower(tools::file_ext(path)) == "xls") {
        stop("Writing legacy .xls files is not supported; use .xlsx.", call. = FALSE)
      }
      # Excel stores clock time without a zone and readxl reads it as UTC;
      # writexl >= 2.0 writes local wall-clock time, so write UTC explicitly
      x <- .posix_to_utc(x)
      if (requireNamespace("writexl", quietly = TRUE)) {
        writexl::write_xlsx(x, path, ...)
      } else if (requireNamespace("openxlsx", quietly = TRUE)) {
        openxlsx::write.xlsx(x, path, overwrite = TRUE, ...)
      } else {
        stop("Writing Excel files requires the 'writexl' or 'openxlsx' package.", call. = FALSE)
      }
    },
    spss = .write_haven(x, path, "spss", ...),
    stata = .write_haven(x, path, "stata", ...),
    xpt = .write_haven(x, path, "xpt", ...),
    sas = stop("Writing sas7bdat files is not supported; use format 'xpt' instead.", call. = FALSE),
    fst = {
      .require("fst", "to write fst files")
      fst::write_fst(x, path, ...)
    }
  )
  invisible(path)
}

.write_delim <- function(x, path, sep, compression, append = FALSE, ...) {
  x <- .posix_to_iso(x)
  if (is.null(compression) && requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(x, path, sep = sep, append = append,
                       col.names = !append, ...)
    return(invisible(path))
  }
  con <- .open_con(path, compression, if (append) "a" else "w")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, con, sep = sep, row.names = FALSE, col.names = !append,
                     quote = TRUE, na = "", qmethod = "double", ...)
  invisible(path)
}

# date-times are written as unambiguous ISO 8601 UTC strings
.posix_to_iso <- function(x) {
  for (nm in names(x)) {
    if (inherits(x[[nm]], "POSIXct")) {
      x[[nm]] <- format(x[[nm]], "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    }
  }
  x
}

.write_ndjson <- function(x, path, compression, append = FALSE, ...) {
  .require("jsonlite", "to write NDJSON files")
  con <- .open_con(path, compression, if (append) "a" else "w")
  on.exit(close(con), add = TRUE)
  jsonlite::stream_out(.posix_to_iso(x), con, verbose = FALSE, digits = NA, na = "null",
                       POSIXt = "ISO8601", Date = "ISO8601", ...)
  invisible(path)
}

.write_parquet <- function(x, path, ...) {
  if (requireNamespace("arrow", quietly = TRUE)) {
    return(arrow::write_parquet(x, path, ...))
  }
  if (requireNamespace("nanoparquet", quietly = TRUE)) {
    return(nanoparquet::write_parquet(x, path, ...))
  }
  stop("Writing Parquet files requires the 'arrow' or 'nanoparquet' package.", call. = FALSE)
}

# Stata variable names: letters, digits and '_', not starting with a digit,
# at most 32 characters and no reserved words.
.stata_names <- function(nm) {
  reserved <- c("_all", "_b", "byte", "_coef", "_cons", "double", "float", "if",
                "in", "int", "long", "_n", "_N", "_pi", "_pred", "_rc", "_skip",
                "strL", "using", "with")
  new <- gsub("[^A-Za-z0-9_]", "_", nm)
  new <- ifelse(grepl("^[A-Za-z_]", new), new, paste0("v", new))
  new <- substr(new, 1L, 32L)
  bad <- new %in% reserved | grepl("^str[0-9]+$", new)
  new[bad] <- substr(paste0(new[bad], "_"), 1L, 32L)
  new <- make.unique(new, sep = "_")
  changed <- new != nm
  if (any(changed)) {
    warning("Column names changed to valid Stata names: ",
            paste(sprintf("%s -> %s", nm[changed], new[changed]), collapse = ", "),
            call. = FALSE)
  }
  new
}

# Formats that store clock time without a time zone (Excel, SPSS, Stata,
# SAS) are written in UTC; the instant is unchanged, only the zone label.
.posix_to_utc <- function(x) {
  for (nm in names(x)) {
    if (inherits(x[[nm]], "POSIXct")) attr(x[[nm]], "tzone") <- "UTC"
  }
  x
}

.write_haven <- function(x, path, kind, ...) {
  .require("haven", "to write SPSS/Stata/SAS transport files")
  x <- .posix_to_utc(x)
  for (nm in names(x)) {
    if (is.logical(x[[nm]])) x[[nm]] <- as.integer(x[[nm]])
  }
  if (kind == "stata") names(x) <- .stata_names(names(x))
  switch(kind,
    spss = haven::write_sav(x, path, ...),
    stata = haven::write_dta(x, path, ...),
    xpt = haven::write_xpt(x, path, ...)
  )
}

# HDF5 ------------------------------------------------------------------------

.write_hdf5 <- function(x, path, dataset = "data") {
  .require("hdf5r", "to write HDF5 files")
  if (file.exists(path)) file.remove(path)
  h5 <- hdf5r::H5File$new(path, mode = "w")
  on.exit(h5$close_all(), add = TRUE)
  grp <- h5$create_group(dataset)
  hdf5r::h5attr(grp, "column_names") <- names(x)
  hdf5r::h5attr(grp, "nrow") <- nrow(x)
  for (j in seq_along(x)) {
    col <- x[[j]]
    nm <- sprintf("col%05d", j)
    cls <- .value_class(col)
    if (is.na(cls)) cls <- "character"
    na_idx <- which(is.na(col))
    vals <- switch(cls,
      factor = , ordered = as.integer(col),
      logical = as.integer(col),
      Date = , POSIXct = , numeric_classed = as.numeric(col),
      character = {
        v <- as.character(col)
        v[is.na(v)] <- ""
        v
      },
      integer = col,
      numeric = col,
      as.character(col)
    )
    if (cls %in% c("factor", "ordered", "logical")) vals[is.na(vals)] <- -1L
    if (cls == "numeric_classed") cls <- "numeric"
    ds <- grp$create_dataset(nm, robj = vals)
    hdf5r::h5attr(ds, "name") <- names(x)[j]
    hdf5r::h5attr(ds, "r_class") <- cls
    if (cls %in% c("factor", "ordered")) hdf5r::h5attr(ds, "levels") <- levels(col)
    if (cls == "POSIXct") hdf5r::h5attr(ds, "tz") <- attr(col, "tzone") %||% ""
    if (cls == "character" && length(na_idx)) hdf5r::h5attr(ds, "na_index") <- na_idx
    ds$close()
  }
  grp$close()
  invisible(path)
}

.read_hdf5 <- function(path, dataset = "data") {
  .require("hdf5r", "to read HDF5 files")
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  if (!h5$exists(dataset)) {
    stop(sprintf("Object '%s' not found in HDF5 file. Available: %s", dataset,
                 paste(names(h5), collapse = ", ")), call. = FALSE)
  }
  obj <- h5[[dataset]]
  if (inherits(obj, "H5D")) {
    v <- obj$read()
    return(if (is.null(dim(v))) data.frame(x = v) else as.data.frame(v))
  }
  members <- sort(names(obj))
  attrs <- hdf5r::h5attr_names(obj)
  cols <- lapply(members, function(m) {
    ds <- obj[[m]]
    on.exit(ds$close(), add = TRUE)
    v <- ds$read()
    an <- hdf5r::h5attr_names(ds)
    get_attr <- function(a) if (a %in% an) hdf5r::h5attr(ds, a) else NULL
    cls <- get_attr("r_class")
    nm <- get_attr("name") %||% m
    if (!is.null(cls)) {
      v <- switch(cls,
        factor = , ordered = {
          v[v < 0] <- NA
          factor(get_attr("levels")[v], levels = get_attr("levels"),
                 ordered = cls == "ordered")
        },
        logical = {
          v[v < 0] <- NA
          as.logical(v)
        },
        Date = as.Date(v, origin = "1970-01-01"),
        POSIXct = as.POSIXct(v, origin = "1970-01-01", tz = get_attr("tz") %||% ""),
        character = {
          idx <- get_attr("na_index")
          if (length(idx)) v[idx] <- NA
          v
        },
        integer = as.integer(v),
        v
      )
    }
    list(name = nm, value = v)
  })
  out <- lapply(cols, `[[`, "value")
  names(out) <- vapply(cols, `[[`, character(1), "name")
  if ("column_names" %in% attrs) {
    ord <- hdf5r::h5attr(obj, "column_names")
    if (all(ord %in% names(out))) out <- out[ord]
  }
  obj$close()
  lens <- unique(lengths(out))
  if (length(lens) != 1L) stop("Datasets in the HDF5 group have different lengths.", call. = FALSE)
  attr(out, "row.names") <- .set_row_names(lens)
  class(out) <- "data.frame"
  out
}

# Type alignment --------------------------------------------------------------

#' Coerce the columns of a data frame to the classes of a template
#'
#' Useful after reading formats that do not store R classes (CSV, JSON,
#' Excel): factors regain their levels (and ordering), dates and date-times
#' are parsed, and integers/logicals are restored. Columns not present in
#' the template are left unchanged.
#'
#' Date-time strings ending in `Z` are read as UTC, strings with an explicit
#' offset (e.g. `+02:00`) are converted from that offset, and strings without
#' zone information are read as local times in the time zone of the template
#' column.
#'
#' @param x A data frame.
#' @param template A data frame with the desired column classes.
#' @return `x` with coerced columns.
#' @examples
#' f <- tempfile(fileext = ".csv")
#' write.csv(iris, f, row.names = FALSE)
#' x <- match_types(read.csv(f), iris)
#' identical(levels(x$Species), levels(iris$Species))
#' @export
match_types <- function(x, template) {
  if (!is.data.frame(x) || !is.data.frame(template)) {
    stop("'x' and 'template' must be data frames.", call. = FALSE)
  }
  for (nm in intersect(names(x), names(template))) {
    tc <- template[[nm]]
    v <- x[[nm]]
    cls <- .value_class(tc)
    if (is.na(cls)) next
    x[[nm]] <- switch(cls,
      factor = , ordered = .as_factor_like(v, tc, cls == "ordered"),
      character = as.character(v),
      logical = if (is.logical(v)) {
        v
      } else if (is.numeric(v)) {
        as.logical(v)
      } else {
        u <- toupper(trimws(as.character(v)))
        ifelse(u %in% c("TRUE", "T", "1", "YES"), TRUE,
               ifelse(u %in% c("FALSE", "F", "0", "NO"), FALSE, NA))
      },
      integer = as.integer(v),
      numeric = as.numeric(v),
      numeric_classed = .restore_numeric_classed(v, tc),
      Date = if (inherits(v, "Date")) v else if (is.numeric(v)) {
        as.Date(v, origin = "1970-01-01")
      } else {
        as.Date(substr(as.character(v), 1L, 10L))
      },
      POSIXct = {
        tz <- attr(tc, "tzone") %||% ""
        if (inherits(v, "POSIXct")) {
          attr(v, "tzone") <- tz
          v
        } else if (is.numeric(v)) {
          as.POSIXct(v, origin = "1970-01-01", tz = tz)
        } else {
          .parse_datetime(as.character(v), tz)
        }
      },
      v
    )
  }
  x
}

.as_factor_like <- function(v, tc, ordered) {
  lev <- levels(tc)
  vc <- as.character(v)
  ok <- !is.na(vc)
  # integer codes (e.g. from formats that drop value labels)
  if (is.numeric(v) && any(ok) && !any(vc[ok] %in% lev) &&
      all(v[ok] == round(v[ok])) && all(v[ok] >= 1 & v[ok] <= length(lev))) {
    vc <- lev[v]
  }
  factor(vc, levels = lev, ordered = ordered)
}

# Parse date-time strings. Strings ending in "Z" are UTC and strings with an
# explicit offset ("+02:00", "-0500") are converted from that offset; strings
# without zone information are local times in the template's time zone.
.parse_datetime <- function(v, tz) {
  v <- trimws(sub("T", " ", as.character(v), fixed = TRUE))
  out <- rep(NA_real_, length(v))
  utc <- grepl("Z$", v)
  off_re <- "([+-])([0-9]{2}):?([0-9]{2})$"
  has_off <- !utc & grepl(paste0("[0-9]", off_re), v)
  local <- !utc & !has_off
  if (any(utc)) out[utc] <- .parse_dt_formats(sub("Z$", "", v[utc]), "UTC")
  if (any(has_off)) {
    m <- regmatches(v[has_off], regexec(off_re, v[has_off]))
    secs <- vapply(m, function(z) {
      (if (z[2L] == "-") -1 else 1) * (as.numeric(z[3L]) * 3600 + as.numeric(z[4L]) * 60)
    }, numeric(1))
    out[has_off] <- .parse_dt_formats(sub(off_re, "", v[has_off]), "UTC") - secs
  }
  if (any(local)) out[local] <- .parse_dt_formats(v[local], tz)
  structure(out, class = c("POSIXct", "POSIXt"), tzone = tz)
}

.parse_dt_formats <- function(x, tz) {
  formats <- c("%Y-%m-%d %H:%M:%OS", "%Y-%m-%d %H:%M", "%Y-%m-%d",
               "%Y/%m/%d %H:%M:%OS", "%Y/%m/%d %H:%M", "%Y/%m/%d")
  out <- rep(NA_real_, length(x))
  for (f in formats) {
    todo <- is.na(out) & !is.na(x) & nzchar(x)
    if (!any(todo)) break
    out[todo] <- as.numeric(as.POSIXct(x[todo], tz = tz, format = f))
  }
  out
}

# Chunked generation to files ----------------------------------------------------

#' Generate a large synthetic data set directly to disk
#'
#' @description
#' Generates `n` synthetic rows in chunks of `chunk_size`, so that data sets
#' much larger than memory can be produced. For appendable formats (`csv`,
#' `tsv`, `ndjson`, optionally compressed) all chunks go to one file. For
#' other formats (Parquet, Feather, HDF5, fst, RDS, ...) a directory `path`
#' is created holding one part file per chunk (`part-00001.parquet`, ...),
#' which [read_data()] (or e.g. `arrow::open_dataset()`) reads as a whole.
#'
#' @param object A fitted `sp_synthesizer`.
#' @param n Total number of rows.
#' @param path Output file (appendable formats) or directory (other formats).
#' @param format Format name; inferred from `path` for appendable formats.
#'   Required when `path` is a directory name without extension.
#' @param chunk_size Rows per chunk.
#' @param seed Optional random seed (the whole output is reproducible).
#' @param progress Print progress messages?
#' @param ... Passed to [generate()] (e.g. `dependence`).
#' @return Invisibly, a list with the `path`, `format`, number of rows and
#'   the files written.
#' @examples
#' fit <- fit_synthesizer(iris)
#' out <- tempfile(fileext = ".csv")
#' generate_to_file(fit, n = 1000, path = out, chunk_size = 300, seed = 1)
#' nrow(read_data(out))
#'
#' dir <- tempfile()
#' generate_to_file(fit, n = 1000, path = dir, format = "rds", chunk_size = 400)
#' list.files(dir)
#' @export
generate_to_file <- function(object, n, path, format = NULL, chunk_size = 1e5,
                             seed = NULL, progress = interactive(), ...) {
  if (!inherits(object, "sp_synthesizer")) {
    stop("'object' must be an sp_synthesizer.", call. = FALSE)
  }
  n <- .check_count(n)
  chunk_size <- .check_count(chunk_size, "chunk_size", allow_zero = FALSE)
  fmt <- .detect_format(path, format)
  appendable <- fmt$format %in% c("csv", "tsv", "ndjson")
  nchunks <- max(1L, ceiling(n / chunk_size))
  sizes <- rep(chunk_size, nchunks)
  sizes[nchunks] <- n - chunk_size * (nchunks - 1L)
  ext <- .canonical_ext(fmt$format)

  files <- character(0)
  if (!appendable) {
    if (file.exists(path) && !dir.exists(path)) stop("'path' exists and is not a directory.", call. = FALSE)
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    old <- list.files(path, pattern = "^part-[0-9]+\\.", full.names = TRUE)
    if (length(old)) file.remove(old)
  }

  .with_seed(seed, {
    for (i in seq_len(nchunks)) {
      chunk <- generate(object, n = sizes[i], ...)
      if (object$input == "vector") chunk <- data.frame(x = chunk)
      if (appendable) {
        if (fmt$format == "ndjson") {
          .write_ndjson(chunk, path, fmt$compression, append = i > 1L)
        } else {
          .write_delim(chunk, path, if (fmt$format == "csv") "," else "\t",
                       fmt$compression, append = i > 1L)
        }
      } else {
        f <- file.path(path, sprintf("part-%05d.%s", i, ext))
        write_data(chunk, f, format = fmt$format)
        files <- c(files, f)
      }
      if (progress) message(sprintf("chunk %d/%d written", i, nchunks))
    }
  })
  invisible(list(path = path, format = fmt$format, n = n,
                 files = if (appendable) path else files))
}

.canonical_ext <- function(format) {
  switch(format, rdata = "rda", feather = "feather", hdf5 = "h5", excel = "xlsx",
         spss = "sav", stata = "dta", format)
}
