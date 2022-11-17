.fixHOME <- function(x) {
    gsub('~', Sys.getenv('HOME'), x)
}
