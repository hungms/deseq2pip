#' Log renv Snapshot   
#'
#' This function logs the renv snapshot to the specified directory.
#'
#' @param save_dir Directory to save the renv snapshot. Default is the current working directory
#' @return Nothing, but logs the renv snapshot to the specified directory
#' @examples
#' \dontrun{
#' # Log renv snapshot
#' log_renv(save_dir = "results")
#' }
#' @export
log_renv <- function(save_dir) {
  # create logs directory
  validate_paths(save_dir)
  dir.create(paste0(save_dir, "/logs"), showWarnings = FALSE, recursive = TRUE)

  # log sessionInfo
  if(!dir.exists(paste0(save_dir, "/logs/renv"))){
    message("renv not initialized, initializing renv")
    writeLines(capture.output(sessionInfo()), paste0(save_dir, "/logs/sessionInfo.Rmd"))
    renv::init(project = paste0(save_dir, "/logs"), bare = TRUE, force = TRUE)
    renv::snapshot(project = paste0(save_dir, "/logs"), prompt = FALSE)}
  else{
    message("renv already initialized, restoring renv")
    renv::restore(project = paste0(save_dir, "/logs"), prompt = FALSE)
    }
}

#' Log output
#' 
#' This function logs the script to the specified directory.
#'
#' @param save_dir Directory to save the script. Default is the current working directory
#' @param expr Expression to evaluate and log. If NULL, only logs the function call and arguments
#' @return Nothing, but logs the script to the specified directory
#' @examples
#' \dontrun{
#' # Log script
#' log_script(save_dir = "results")
#' }
#' @export
log_output <- function(call, save_dir, expr = NULL) {

  # create logs directory
  validate_paths(save_dir)
  dir.create(paste0(save_dir, "/logs"), showWarnings = FALSE, recursive = TRUE)
  
  # log renv snapshot
  log_renv(save_dir)

  # time
  time <- c(
    paste0("Date : ", Sys.time()), 
    " ")
  # function arguments
  args_list <- as.list(call)[-1]
  args <- sapply(names(args_list), function(n) {
    expr <- call[[n]]
    # Only evaluate if it's a symbol or a constant
    if (is.symbol(expr) || is.atomic(expr)) {
      val <- eval(expr, envir = parent.frame())
      if (is.atomic(val) && (is.character(val) || is.numeric(val) || is.logical(val))) {
        paste(n, "=", paste(deparse(val), collapse = " "))
      } else {
        paste(n, "=", paste(deparse(expr), collapse = " "))
      }
    } else {
      # For calls or complex objects, just deparse the expression
      paste(n, "=", paste(deparse(expr), collapse = " "))
    }
  })
  args <- c(
    "Arguments : ",
    args,
    " ")
  
  # message
  message <- c(
    "Log Messages : ")

  # Combine all logs
  log_vec <- c(time, args, message)

  # Get the next log file number
  log_files <- list.files(paste0(save_dir, "/logs/"), pattern = "^logfile\\.[0-9]+$")
  log_nums <- as.numeric(sub("logfile\\.", "", log_files))
  log_num <- if (length(log_nums) > 0) max(log_nums, na.rm = TRUE) + 1 else 1
  log_file <- paste0(save_dir, "/logs/logfile.", log_num)
  writeLines(log_vec, log_file)

  # If expr is provided, capture all output/messages/warnings/errors
  if (!is.null(expr)) {
    con <- file(log_file, open = "a")
    sink(con, append = TRUE, split = TRUE)           # split=TRUE for output
    sink(con, append = TRUE, type = "message")      # split not allowed for message
    on.exit({
      sink(type = "message")
      sink()
      close(con)
    }, add = TRUE)
    withCallingHandlers(
      eval(expr, envir = parent.frame()),
      warning = function(w) {
        message("WARNING: ", conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    # If an error occurs, it will stop here and propagate up
  }
}


#' Suppress messages
#' 
#' This function suppresses messages.
#' 
#' @param message Message to suppress.
#' @return Nothing, but suppresses the message.
#' @keywords internal
#' @export
quiet <- function(expr) {
  invisible(
    suppressWarnings(
      suppressMessages(
        capture.output({
          result <- eval.parent(substitute(expr))
        })
      )
    )
  )
  result
}
