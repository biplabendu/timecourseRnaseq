#' Generate a Zero-Centered Color Palette
#'
#' Creates a custom color palette and breaks for heatmap plotting. The scale
#' is centered at zero (white), with negative values transitioning to the
#' first color and positive values to the third color.
#'
#' @param mat A numeric matrix used to determine the range (min/max) of the
#'   data for calculating breaks.
#' @param color_set A character vector of length 3 defining the low, mid
#'   (zero), and high colors. Defaults to \code{c("maroon", "white",
#'   "orange")}.
#' @param n_colors Integer. The desired total number of colors in the palette.
#'
#' @return A list containing two elements:
#'   \item{p_colors}{A vector of color codes.}
#'   \item{p_breaks}{A numeric vector of break points for the color scale.}
#'
#' @export
#'
generate_color_pallete <- function(mat,
                                   color_set = c("maroon", "white", "orange"),
                                   n_colors) {
  # Define the color palette
  # For values > 0: yellow to white
  # For values < 0: white to purple
  colors <- colorRampPalette(color_set)(n_colors)

  # Determine the breaks for the color scale to ensure 0 is exactly white
  # Find the range of the data
  min_val <- min(mat)
  max_val <- max(mat)

  # Calculate the proportion of negative and positive values relative to the total range
  neg_proportion <- abs(min_val) / (max_val - min_val)
  pos_proportion <- max_val / (max_val - min_val)

  # Adjust the number of colors for each segment to ensure white is at 0
  # This is a common approach to center the color at a specific value like 0
  num_colors_neg <- round(neg_proportion * 100)
  num_colors_pos <- round(pos_proportion * 100)

  # Create the specific color palette for the breaks
  min_to_zero <- colorRampPalette(color_set[1:2])(num_colors_neg + 1)
  zero_to_max <- colorRampPalette(color_set[2:3])(num_colors_pos + 1)

  # Combine the palettes, removing the duplicate 'white' in the middle
  p_colors <- c(min_to_zero[1:num_colors_neg], zero_to_max)
  p_breaks <- seq(min_val, max_val, length.out = length(p_colors))

  list(p_colors, p_breaks) |>
    setNames(
      c("p_colors", "p_breaks")
    )
}
