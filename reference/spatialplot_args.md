# Arguments for spatial plot functions

Arguments for spatial plot functions

## Arguments

- object:

  A Seurat object or a Giotto object.

- fov:

  The name of the field of view (FOV) to plot, only works for Seurat
  objects.

- boundaries:

  The name of the boundaries to plot, only works for Seurat objects.

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

- masks:

  Logical, whether to plot masks. Not supported for Seurat or Giotto
  objects for now.

- shapes:

  Plot shapes on the spatial plot. Not supported for Seurat objects with
  Visium or SlideSeq data. For Seurat objects with FOV, when TRUE,
  `boundaries` will be used as boundaries for the shapes. Otherwise
  itself will be used as boundaries. For Giotto objects or Seurat
  objects with FOV, this defaults to TRUE when `shapes_fill_by` is
  provided. Set to FALSE to disable shapes plotting.

- points:

  Logical, whether to plot points. If TRUE, the points will be plotted
  using the coordinates from the object. Defaults to TRUE.

- ext:

  The extent of the plot. If NULL, the extent will be calculated from
  the data. If a numeric vector of length 4, it should be in the format
  c(xmin, xmax, ymin, ymax). It can also be an object created by
  [`terra::ext()`](https://rspatial.github.io/terra/reference/ext.html).

- crop:

  Whether to crop the plot to the extent of the available data. Similar
  to `crop` argument in
  [`Seurat::SpatialDimPlot()`](https://satijalab.org/seurat/reference/SpatialPlot.html).
  Defaults to TRUE.

- group_by:

  The name of the metadata column to group the points by. Should be a
  character or factor column. A special value "molecules" can be used to
  plot molecules in the FOV.

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

- scale_factor:

  Internal use only. The scale factor to use for the image, which will
  be extracted from the object.

- layers:

  A character vector of layers to plot. Possible values are:

  - "image": plot the image as a background, which should be the first
    layer if provided.

  - "masks": plot the masks

  - "shapes": plot the shapes

  - "points": plot the points The order of the layers matters, as the
    first layer will be plotted first. You can also use it to disable
    some layers by excluding them from the vector.

- flip_y:

  Internal use mostly, unless you want to flip the y-axis of the plot.

- padding:

  The padding to add to the extent of the plot, only available when
  `crop = TRUE` and `ext = NULL`. For Seurat objects with FOV, this
  defaults to 0. In other cases, this defaults to 0 when image is
  plotted, and 0.05 otherwise.

- image_scale:

  Choose the scale factor ("lowres"/"hires") to apply in order to
  matchthe plot with the specified `image`. Similar to `image.scale`
  argument in
  [`Seurat::SpatialDimPlot()`](https://satijalab.org/seurat/reference/SpatialPlot.html).

- x:

  Internal use only, the name of the x coordinate column in the data.
  Used to adopt different data types.

- y:

  Internal use only, the name of the y coordinate column in the data.
  Used to adopt different data types.

- nmols:

  Max number of each molecule specified in `features` for dim plot
  Similar to `nmols` argument in
  [`Seurat::ImageDimPlot()`](https://satijalab.org/seurat/reference/ImageDimPlot.html).
  It also applied to Giotto objects.

- shapes_fill_by:

  The name of the variable to fill the shapes by. It can also be a color
  name, in which case the shapes will be filled with that color. When
  this is provided, the `shapes` argument will be set to TRUE by
  default.

- graph:

  The name of the graph to use for the spatial plot. Currently only
  supported for Giotto objects. The graph data is obtained using
  [`GiottoClass::getSpatialNetwork()`](https://drieslab.github.io/GiottoClass/reference/getSpatialNetwork.html).
  When TRUE, the default graph will be used. When given as a character,
  it should be the name of the graph to use. If there is ":" in the
  name, the first part will be used as spat_unit, and the second part as
  the graph name.

- shape:

  The shape of the points, alias of `points_shape`. See
  <https://ggplot2.tidyverse.org/reference/aes_linetype_size_shape.html>
  for more details.

- legend.position:

  The position of the legend. Defaults to "right".

- legend.direction:

  The direction of the legend. Defaults to "vertical".

- theme:

  The theme to use for the plot. Defaults to `"theme_box"`. It can be
  the name of the theme (e.g. "ggplot2::theme_bw") or the function
  itself. There are three themes that can be passed without namespace:
  "theme_box", "theme_this" and "theme_blank", which are actually
  aliases of
  [`plotthis::theme_box()`](https://pwwang.github.io/plotthis/reference/theme_box.html),
  [`plotthis::theme_this()`](https://pwwang.github.io/plotthis/reference/theme_this.html)
  and
  [`ggplot2::theme_void()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)
  (without braces).

- theme_args:

  A list of arguments to pass to the theme function.

- title:

  The title of the plot. If NULL, no title will be added.

- subtitle:

  The subtitle of the plot. If NULL, no subtitle will be added.

- xlab:

  The label for the x-axis. If NULL, no label will be added.

- ylab:

  The label for the y-axis. If NULL, no label will be added.

- facet_scales:

  The scales to use for the facets. Defaults to "free". Can be "free",
  "fixed", "free_x", "free_y".

- facet_nrow:

  The number of rows to use for the facets. Defaults to NULL, which
  means the number of rows will be calculated automatically.

- facet_ncol:

  The number of columns to use for the facets. Defaults to NULL, which
  means the number of columns will be calculated automatically.

- facet_byrow:

  Logical, whether to facet by row. Defaults to FALSE.

- feat_type:

  feature type of the features (e.g. "rna", "dna", "protein"), only
  applied to Giotto objects.

- use_overlap:

  use polygon and feature coordinates overlap results, only applied to
  Giotto objects.

- shapes_feat_type:

  feature type of the features to use for shapes (e.g. "rna", "dna",
  "protein"), only applied to Giotto objects.

- shapes_alpha:

  The alpha value to use for the shapes. When "points" are plotted, this
  defaults to 0.5; otherwise it defaults to 1.

- spat_unit:

  The spatial unit to use for the plot. Only applied to Giotto objects.

- spat_loc_name:

  The name of the spatial location to use for the plot. Only applied to
  Giotto objects.

- spat_enr_names:

  The names of the spatial enrichment results to use for the plot. Only
  applied to Giotto objects.

- ...:

  Additional arguments that will be passed to the spatial plot function.
  With the `image_` prefix, these arguments will be used to plot the
  image
  ([`plotthis::SpatImagePlot()`](https://pwwang.github.io/plotthis/reference/spatialplots.html)).
  With the `masks_` prefix, these arguments will be used to plot the
  masks
  ([`plotthis::SpatMasksPlot()`](https://pwwang.github.io/plotthis/reference/spatialplots.html)).
  With the `shapes_` prefix, these arguments will be used to plot the
  shapes
  ([`plotthis::SpatShapesPlot()`](https://pwwang.github.io/plotthis/reference/spatialplots.html)).
  With the `points_` prefix, these arguments will be used to plot the
  points
  ([`plotthis::SpatPointsPlot()`](https://pwwang.github.io/plotthis/reference/spatialplots.html)).
  If the prefix is not provided, the arguments will be used as points
  arguments, but with lower priority than the `points_` prefixed
  arguments.
