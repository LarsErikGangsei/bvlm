#' Hello from BMLM
#' @param name Name of the person you want to Halloo
#' @return Print hello message to the screen
#' @examples hello("user")
#' @import simrel
#' @importFrom reshape2 melt
#' @export
hello <- function(name) {
  print(paste0("Hello, ", name, "!"))
}
