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

  A Seurat object (with spatial data loaded via SeuratObject) or a
  Giotto object (created with Giotto). The spatial technology is
  auto-detected from the object's image class.

- fov:

  The name of the field of view (FOV) to plot. For Seurat FOV-based
  objects (Xenium, CosMx, etc.), defaults to
  [`SeuratObject::DefaultFOV()`](https://satijalab.github.io/seurat-object/reference/DefaultFOV.html).
  Not applicable to Visium or Slide-seq objects.

- boundaries:

  The name of the segmentation boundaries within the FOV to use for cell
  outlines. For Seurat FOV-based objects, defaults to
  [`SeuratObject::DefaultBoundary()`](https://satijalab.github.io/seurat-object/reference/Boundaries.html).
  Not applicable to Visium or Slide-seq objects.

- x:

  Internal use. The name of the x-coordinate column in the spatial data.
  Auto-detected based on the spatial technology (`"imagerow"` for Visium
  V1, varies for other types).

- y:

  Internal use. The name of the y-coordinate column in the spatial data.
  Auto-detected based on the spatial technology (`"imagecol"` for Visium
  V1, varies for other types).

- swap_xy:

  Logical, whether to swap x and y coordinates in the points data.
  Seurat objects, loaded with Visium, Xenium or SlideSeq data, have x
  and y coordinates in the opposite order. But when loaded with FOV
  data, the coordinates are in the correct order.

- image:

  Controls the image/background layer. Possible values:

  - `NULL` — Default behavior: for Visium, the first image is used; for
    Giotto and FOV objects, no image is plotted.

  - A character string naming an image in the object — that specific
    image is plotted.

  - A color name (e.g., `"white"`, `"lightgray"`) — fills the background
    with a solid color rectangle.

  - `TRUE` — For Visium: uses the first image. For Giotto FOV: plots all
    non-overlapping images. For Seurat FOV: raises an error (no single
    default image).

  - `FALSE` — Disables the image layer entirely.

- crop:

  Logical. Whether to crop the plot to the extent of the tissue/spots.
  When `TRUE` (default), the plot is automatically zoomed to the data
  extent with optional `padding`. Analogous to the `crop` argument in
  [`Seurat::SpatialDimPlot()`](https://satijalab.org/seurat/reference/SpatialPlot.html).

- scale_factor:

  Internal use. The image scale factor extracted from the object, used
  to map between pixel and tissue coordinate spaces. Automatically
  determined from the object's image data.

- group_by:

  A metadata column name used to color the points. Must be a character
  or factor column in the object's metadata. For
  [`SpatDimPlot()`](https://pwwang.github.io/scplotter/reference/SpatDimPlot.md),
  if `NULL` and the object has FOV data with features, defaults to
  `"molecules"`; otherwise defaults to `"Identity"` (the active cluster
  identities). For
  [`SpatFeaturePlot()`](https://pwwang.github.io/scplotter/reference/SpatFeaturePlot.md),
  `group_by` is ignored — use
  [`SpatDimPlot`](https://pwwang.github.io/scplotter/reference/SpatDimPlot.md)
  for categorical grouping. The special value `"molecules"` enables
  molecule-level visualization in FOV-based data.

- shape:

  Numeric. The point shape (ggplot2 shape aesthetic). Default: `16`
  (filled circle). See
  <https://ggplot2.tidyverse.org/reference/aes_linetype_size_shape.html>
  for the full shape palette.

- features:

  A character vector of feature names to visualize. For
  [`SpatFeaturePlot()`](https://pwwang.github.io/scplotter/reference/SpatFeaturePlot.md),
  each feature is plotted as a separate facet (or combined theme when a
  single feature is given), with expression values coloring the points.
  For
  [`SpatDimPlot()`](https://pwwang.github.io/scplotter/reference/SpatDimPlot.md),
  features are treated as molecule names to plot at single-molecule
  resolution (FOV data only). Can include gene names, metadata column
  names, or dimension reduction components.

- layer:

  The assay layer from which to extract feature expression values. For
  Seurat objects, one of `"data"` (default), `"scale.data"`, or
  `"counts"`. For Giotto objects, one of `"normalized"` (default),
  `"scaled"`, `"raw"`, `"counts"`, or `"custom"`.

- legend.position:

  Character. Legend position. One of `"right"` (default), `"left"`,
  `"top"`, `"bottom"`, or `"none"`.

- legend.direction:

  Character. Legend direction. One of `"vertical"` (default) or
  `"horizontal"`.

- flip_y:

  Logical. Whether to flip the y-axis. This is primarily for internal
  coordinate system alignment — Visium/Slide-seq objects default to
  `TRUE`, FOV-based objects to `FALSE`. In most cases you do not need to
  set this manually.

- ext:

  The spatial extent (bounding box) of the plot. If `NULL`, the extent
  is calculated automatically from the data or, when `crop = TRUE`, from
  the tissue coordinates. Can be a numeric vector in the format
  `c(xmin, xmax, ymin, ymax)` or a
  [`terra::SpatExtent`](https://rspatial.github.io/terra/reference/ext.html)
  object.

## Value

A list containing the ggplot2 layer object and the facet_by variable if
applicable.
