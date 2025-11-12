# Process points layer for Seurat spatial plots

Process points layer for Seurat spatial plots

## Usage

``` r
.seurat_points_layer(
  object,
  fov = NULL,
  boundaries = NULL,
  x = "x",
  y = "y",
  swap_xy = TRUE,
  image,
  args,
  crop,
  points_data,
  ext_unscaled,
  scale_factor,
  group_by,
  shape,
  features,
  layer,
  legend.position,
  legend.direction,
  flip_y,
  ext
)
```

## Arguments

- object:

  A Seurat object or a Giotto object.

- fov:

  The name of the field of view (FOV) to plot, only works for Seurat
  objects.

- boundaries:

  The name of the boundaries to plot, only works for Seurat objects.

- x:

  Internal use only, the name of the x coordinate column in the data.
  Used to adopt different data types.

- y:

  Internal use only, the name of the y coordinate column in the data.
  Used to adopt different data types.

- swap_xy:

  Logical, whether to swap x and y coordinates in the points data.
  Seurat objects, loaded with Visium, Xenium or SlideSeq data, have x
  and y coordinates in the opposite order. But when loaded with FOV
  data, the coordinates are in the correct order.

- image:

  The name of the image to plot. Possible values are:

  - NULL: For Seurat objects with Visium data, the first image will be
    used. For Giotto objects and Seurat objects with other spatial data,
    no image will be plotted.

  - image name(s): the name of the image(s) to plot.

  - color name: a color to use as a background for the plot.

  - TRUE: For Seurat objects with Visium data, the first image will be
    used. For Seurat objects with other spatial data, an error will be
    raised. For Giotto objects with FOV, all (non-overlapping) images
    will be plotted; otherwise first image will be used.

  - FALSE: disable image plotting.

- crop:

  Whether to crop the plot to the extent of the available data. Similar
  to `crop` argument in
  [`Seurat::SpatialDimPlot()`](https://satijalab.org/seurat/reference/SpatialPlot.html).
  Defaults to TRUE.

- scale_factor:

  Internal use only. The scale factor to use for the image, which will
  be extracted from the object.

- group_by:

  The name of the metadata column to group the points by. Should be a
  character or factor column. A special value "molecules" can be used to
  plot molecules in the FOV.

- shape:

  The shape of the points, alias of `points_shape`. See
  <https://ggplot2.tidyverse.org/reference/aes_linetype_size_shape.html>
  for more details.

- features:

  A character vector of feature names to plot. If provided, the points
  will be colored by the features. For
  [`SpatDimPlot()`](https://pwwang.github.io/scplotter/reference/SpatDimPlot.md),
  this will be used to plot the molecules in the FOV. For
  [`SpatFeaturePlot()`](https://pwwang.github.io/scplotter/reference/SpatFeaturePlot.md),
  the plots will be faceted by the features.

- layer:

  The layer to use for the feature expression data. Applicable for both
  Seurat and Giotto objects. Defaults to "data" for Seurat objects, and
  "normalized" for Giotto objects. For Giotto objects, it can also be
  "scaled", "raw", "counts", or "custom". For Seurat objects, it can be
  "data", "scale.data", or "counts".

- legend.position:

  The position of the legend. Defaults to "right".

- legend.direction:

  The direction of the legend. Defaults to "vertical".

- flip_y:

  Internal use mostly, unless you want to flip the y-axis of the plot.

- ext:

  The extent of the plot. If NULL, the extent will be calculated from
  the data. If a numeric vector of length 4, it should be in the format
  c(xmin, xmax, ymin, ymax). It can also be an object created by
  [`terra::ext()`](https://rspatial.github.io/terra/reference/ext.html).

## Value

A list containing the ggplot2 layer object and the facet_by variable if
applicable.
