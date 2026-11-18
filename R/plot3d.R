################################################################################
# plot3D is optional (Suggests).                                               #
#                                                                              #
# plot3D loads tcltk (through misc3d), and tcltk needs XQuartz on macOS, so     #
# importing plot3D made every macOS user install XQuartz just to load Rsafd.   #
# The package code calls persp3D() and mesh() by name; these internal          #
# wrappers forward to plot3D when it is installed and otherwise explain how    #
# to install it, so only the 3D plots depend on plot3D.                        #
################################################################################

.plot3d_available <- function() {
  requireNamespace("plot3D", quietly = TRUE)
}

.need_plot3d <- function() {
  if (!.plot3d_available()) {
    stop("3D plots need the plot3D package. Install it with ",
         "install.packages(\"plot3D\") (on macOS, XQuartz is also required: ",
         "https://www.xquartz.org).", call. = FALSE)
  }
}

persp3D <- function(...) {
  .need_plot3d()
  plot3D::persp3D(...)
}

mesh <- function(...) {
  .need_plot3d()
  plot3D::mesh(...)
}
