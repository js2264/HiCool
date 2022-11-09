.fix_HOME <- function(x) {
    gsub('~', Sys.getenv('HOME'), x)
}
